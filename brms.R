# Polyommatus Animal Model - GxE with brms

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(tidyverse,brms,nadiv,emmeans,tidybayes,broom,broom.mixed,
               patchwork,ggh4x,ggtext,MetBrewer)
#ERROR loading broom.mixed 

# Backend: cmdstanr is recommended (faster, more stable) but rstan works too.
# To install cmdstanr (not on CRAN):
  # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
  # cmdstanr::check_cmdstan_toolchain(fix = TRUE)
  # cmdstanr::install_cmdstan()
# Then uncomment the next line:
options(mc.cores = 4,
        brms.backend = "cmdstanr")
# If cmdstanr is not installed, brms falls back to rstan automatically.

# Set some global Stan options
CHAINS <- 4
ITER <- 4000
WARMUP <- 2000
BAYES_SEED <- 1234

# Use the Johnson color palette
clrs <- MetBrewer::met.brewer("Johnson")

# Tell bayesplot to use the Johnson palette (for things like pp_check())
bayesplot::color_scheme_set(c("grey30", clrs[2], clrs[1], clrs[3], clrs[5], clrs[4]))

# Custom ggplot theme to make pretty plots
# Get the font at https://fonts.google.com/specimen/Jost
theme_nice <- function() {
  theme_minimal(base_family = "Jost") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(family = "Jost", face = "bold"),
          axis.title = element_text(family = "Jost Medium"),
          strip.text = element_text(family = "Jost", face = "bold",
                                    size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA))
}

######################################
# 1. PREPROCESSING
######################################

#load data
bluesum <- read.csv("bluesum.csv")

str(bluesum)
# convert categorical variables
bluesum$motherID = as.character(bluesum$motherID)
bluesum$region = as.factor(bluesum$region)
bluesum$temp = as.factor(bluesum$temp)
bluesum$sex = as.factor(bluesum$sex)

#subset data
blueFdata <- subset(bluesum, sex=="F")
blueMdata <- subset(bluesum, sex=="M")


n_families <- length(unique(blueFdata$motherID)) # 42 families
n_offspring <- length(blueFdata$offspringID) # 327 daughters

# Create data frame
data <- blueFdata %>%
  mutate(
    animal = paste0("ind", offspringID),
    MotherID = as.factor(motherID),
    Temperature = as.factor(temp_label),
    Region = as.factor(region_label)
  ) %>%
  select(animal, MotherID, Temperature, Region, TotalArea = avg_total_mm, 
         Blueness = avg_blue_mm)

# Add mothers to pedigree
mothers <- blueFdata %>%
  group_by(motherID) %>%
  summarise(
    Temperature = first(temp_label),
    Region = first(region_label),
    TotalArea = mean(avg_total_mm, na.rm = TRUE),
    Blueness = mean(avg_blue_mm, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    animal = paste0("mom", motherID),
    MotherID = NA,
    Temperature = NA,
    Region = NA,
    TotalArea = NA,
    Blueness = NA
  ) %>%
  select(animal, MotherID, Temperature, Region, TotalArea, Blueness)


##############################################
# 2. Create Pedigree
##############################################
pedigree <- data.frame(
  animal = c(mothers$animal, data$animal),
  dam = c(rep(NA, n_families), paste0("mom", data$MotherID)),
  sire = NA
)

# Create relatedness matrix for brms
Amat <- as.matrix(nadiv::makeA(pedigree))


##############################################
# 4. Fit Animal Model
##############################################

# HURDLE LOGNORMAL
# Intercept only model
model <- brm(
  bf(Blueness ~ TotalArea + Temperature + Region + TotalArea:Temperature + Temperature:Region + 
    (1 | gr(animal, cov = Amat)), # random effect structure
    hu ~ 1), # intercept-only hurdle component
  data = data, # data set (without mothers)
  data2 = list(Amat = Amat),
  family = hurdle_lognormal(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  silent = 2
)
  # Does not converge

# Add mcmc random effect structure
model <- brm(
  bf(
    Blueness ~ TotalArea + Temperature + Region +
      TotalArea:Temperature + Temperature:Region +
      (1 | gr(animal, cov = Amat)) +
      (1 | MotherID) +
      (0 + Temperature || gr(animal, cov = Amat)) +
      (0 + Region || gr(animal, cov = Amat)),
    hu ~ 1  # intercept-only hurdle component
  ),
  data = data,
  data2 = list(Amat = Amat),
  family = hurdle_lognormal(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  silent = 2
)

summary(model)
plot(model)
pp_check(model, ndraws = 100)

# Add hurdle structure: account for zero-inflation in Blueness
model <- brm(
  bf(Blueness ~ TotalArea + Temperature + Region +
       TotalArea:Temperature + Temperature:Region +
       (1 | gr(animal, cov = Amat)) +
       (1 | MotherID) +
       (0 + Temperature || gr(animal, cov = Amat)) +
       (0 + Region || gr(animal, cov = Amat)), 
     hu ~ Temperature + Region),
  data = data, # data set (without mothers)
  data2 = list(Amat = Amat),
  family = hurdle_lognormal(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  silent = 2
)


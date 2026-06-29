# Polyommatus Animal Model - GxE with brms

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(tidyverse,brms,nadiv,emmeans,tidybayes,broom,broom.mixed,
               patchwork,ggh4x,ggtext,MetBrewer,dplyr,grid,gridExtra)

PLOT_DIR <- file.path("plots")

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
ITER <- 6000
WARMUP <- 3000
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

PAL_REGION  <- c("Öland" = "#2CA02C", "Skåne" = "#9467BD")

######################################
# 1. PREPROCESSING
######################################

#load data
bluesum <- read.csv("bluesum.csv")

# convert categorical variables
bluesum$motherID = as.character(bluesum$motherID)
bluesum$region = as.factor(bluesum$region)
bluesum$temp = as.factor(bluesum$temp)
bluesum$sex = as.factor(bluesum$sex)

#subset data
blueFdata <- subset(bluesum, sex=="F")

# z-score blue area
blueFdata <- blueFdata %>%
  mutate(
    z_blue_mm = scale(avg_blue_mm)
  )

# oblueFdata <- subset(blueFdata, region_label=="Öland")
# sblueFdata <- subset(blueFdata, region_label=="Skåne")

# # mutate to set Blueness values below 0.5% to zero
# bluesum <- bluesum %>%
#   mutate(
#     avg_blue_mm   = if_else(avg_prop_blue < 0.005, 0, avg_blue_mm),
#     avg_prop_blue = if_else(avg_prop_blue < 0.005, 0, avg_prop_blue)
#   )

# # add logged blueness and add small constant to zeros
# bluesum <- bluesum %>%
#   mutate(
#     avg_blue_mm = if_else(avg_blue_mm == 0, 0.001, avg_blue_mm),
#     logblue = log(avg_blue_mm)
#   )


n_families <- length(unique(blueFdata$motherID)) # 42 families
n_offspring <- length(blueFdata$offspringID) # 327 daughters

# region subsets
# on_families <- length(unique(oblueFdata$motherID)) # 25 families
# on_offspring <- length(oblueFdata$offspringID) # 175 daughters
# sn_families <- length(unique(sblueFdata$motherID)) # 17 families
# sn_offspring <- length(sblueFdata$offspringID) # 152 daughters

# Create data frame
data <- blueFdata %>%
  mutate(
    animal = paste0("ind", offspringID),
    MotherID = (paste0("mom", motherID)),
    FatherID = (gsub("mom", "dad", MotherID)), # Assuming fatherID is "dad" + motherID)
    Temperature = as.factor(temp_label),
    Region = as.factor(region_label)
  ) %>%
  dplyr::select(animal, MotherID, FatherID, Temperature, Region, TotalArea = avg_total_mm, 
         Blueness = avg_blue_mm, ZBlueness = z_blue_mm)

#region subsets
# data <- sblueFdata %>%
#   mutate(
#     animal = paste0("ind", offspringID),
#     MotherID = (paste0("mom", motherID)),
#     FatherID = (gsub("mom", "dad", MotherID)), # Assuming fatherID is "dad" + motherID)
#     Temperature = as.factor(temp_label)
#   ) %>%
#   dplyr::select(animal, MotherID, FatherID, Temperature, TotalArea = avg_total_mm, 
#                 Blueness = avg_blue_mm, ZBlueness = z_blue_mm)

##############################################
# 2. Create Pedigree
##############################################
pedigree <- data.frame(
  animal = c(unique(data$MotherID), unique(data$FatherID), data$animal),
  dam = c(rep(NA, n_families*2), data$MotherID),
  sire = c(rep(NA, n_families*2), data$FatherID)
)

# Create relatedness matrix for brms
Amat <- as.matrix(nadiv::makeA(pedigree)) 
# Forces full-sibs- still overestimate but more accurate than assuming half-sibs
# for this experimental design where we don't have paternal information and 
# can't separate maternal from additive genetic effects.


##############################################
# 4. Fit Animal Model
##############################################

# HURDLE LOGNORMAL - BAD

# Intercept only model
# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region + 
#        TotalArea:Temperature + Temperature:Region + 
#     (1 | gr(animal, cov = Amat)), # random effect structure
#     hu ~ 1), # intercept-only hurdle component
#   data = data, # data set (without mothers)
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
#   # Does not converge - BAD 23% divergence
# 
# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region + 
#        TotalArea:Temperature + Temperature:Region + 
#        (1 | gr(animal, cov = Amat)),
#      hu ~ (1 | gr(animal, cov = Amat))),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
# # BEST SO FAR - Rhats all </= 1.01 except sigma
# 
# 
# # Add mcmc random effect structure
# model <- brm(
#   bf(
#     Blueness ~ TotalArea + Temperature + Region +
#       TotalArea:Temperature + Temperature:Region +
#       (1 | gr(animal, cov = Amat)) +
#       (1 | MotherID) +
#       (0 + Temperature || gr(animal, cov = Amat)) +
#       (0 + Region || gr(animal, cov = Amat)),
#     hu ~ 1  # intercept-only hurdle component
#   ),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
# 
# # Add hurdle structure: account for zero-inflation in Blueness
# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region +
#        TotalArea:Temperature + Temperature:Region +
#        (1 | gr(animal, cov = Amat)) +
#        (1 | MotherID) +
#        (0 + Temperature || gr(animal, cov = Amat)) +
#        (0 + Region || gr(animal, cov = Amat)), 
#      hu ~ Temperature + Region),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
# 
# # Add full model into hurdle
# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region +
#        TotalArea:Temperature + Temperature:Region +
#        (1 | gr(animal, cov = Amat)) +
#        (1 | MotherID) +
#        (0 + Temperature || gr(animal, cov = Amat)) +
#        (0 + Region || gr(animal, cov = Amat)), 
#      hu ~ TotalArea + Temperature + Region + 
#        TotalArea:Temperature + Temperature:Region +
#        (1 | gr(animal, cov = Amat)) +
#        (1 | MotherID) +
#        (0 + Temperature || gr(animal, cov = Amat)) +
#        (0 + Region || gr(animal, cov = Amat))),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
# # OOF bad ppcheck, bad convergence
# 
# 
# # Mother as fixed effect in hurdle, only random individual? - NO
# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region +
#        TotalArea:Temperature + Temperature:Region +
#        (1 | gr(animal, cov = Amat)) +
#        (1 | MotherID) +
#        (0 + Temperature || gr(animal, cov = Amat)) +
#        (0 + Region || gr(animal, cov = Amat)), 
#      hu ~ TotalArea + Temperature + Region + #MotherID +
#        TotalArea:Temperature + Temperature:Region +
#        (1 | gr(animal, cov = Amat)) +
#        (1 | MotherID)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_lognormal(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )

#####################################
# HURDLE GAMMA
#####################################

# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region + 
#        TotalArea:Temperature + Temperature:Region + 
#        (1 | gr(animal, cov = Amat)),
#      hu ~ (1 | gr(animal, cov = Amat))),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )

# model <- brm(
#   bf(Blueness ~ TotalArea + Temperature + Region +
#        TotalArea:Temperature + Temperature:Region +
#        (0 + Temperature | gr(animal, cov = Amat)) + (1 | MotherID),
#      hu ~ Temperature * Region + (1 | MotherID)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12),
#   silent = 2
# )
# BEST SO FAR - Rhats all </= 1.01, pp_check looks good 
# Remove MotherID- included in Amat, forced full-sibs
# Include priors and shape 
# Add random slopes for region too?

# model <- brm(
#   bf(
#     Blueness ~ TotalArea + Temperature + Region + 
#       TotalArea:Temperature + Temperature:Region + 
#       (0 + Temperature + Region | gr(animal, cov = Amat)),
#     hu ~ Temperature * Region + (0 + Temperature + Region | gr(animal, cov = Amat)),
#     shape ~ Temperature + Region
#      ),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   prior = c(
#     prior(normal(0, 2), class = "Intercept"),
#     prior(normal(0, 1), class = "b"),
#     prior(normal(0, 2), class = "Intercept", dpar = "hu"),
#     prior(normal(0, 1.5), class = "b", dpar = "hu"),
#     prior(normal(0, 2), class = "Intercept", dpar = "shape"),
#     prior(normal(0, 1), class = "b", dpar = "shape"),
#     prior(lkj(2), class = "cor")   # weakly informative on genetic correlation
#   ),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.99, max_treedepth = 15)
# )
# 0 divergences, shape Region Rhat is 1.02, pp_check looks good, model is very complex and may be overfitting.

# Change Temperature + Region to Temperature:Region random slope to capture GxE
# VERY complex model
# model <- brm(
#   bf(
#     Blueness ~ TotalArea + Temperature + Region + 
#       TotalArea:Temperature + Temperature:Region + 
#       (0 + Temperature:Region | gr(animal, cov = Amat)),
#     hu ~ Temperature * Region + (0 + Temperature:Region | gr(animal, cov = Amat)),
#     shape ~ Temperature + Region
#   ),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   prior = c(
#     prior(normal(0, 2), class = "Intercept"),
#     prior(normal(0, 1), class = "b"),
#     prior(normal(0, 2), class = "Intercept", dpar = "hu"),
#     prior(normal(0, 1.5), class = "b", dpar = "hu"),
#     prior(normal(0, 2), class = "Intercept", dpar = "shape"),
#     prior(normal(0, 1), class = "b", dpar = "shape"),
#     prior(lkj(4), class = "cor") # stronger prior on genetic correlation to help with convergence
#   ),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.99, max_treedepth = 15)
# )
# 
# # Remove random region slope from hurdle to help convergence
# model <- brm(
#   bf(
#     Blueness ~ TotalArea + Temperature + Region + 
#       TotalArea:Temperature + Temperature:Region + 
#       (0 + Temperature:Region | gr(animal, cov = Amat)),
#     hu ~ Temperature * Region + (0 + Temperature | gr(animal, cov = Amat)),
#     shape ~ Temperature + Region
#   ),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   prior = c(
#     prior(normal(0, 2), class = "Intercept"),
#     prior(normal(0, 1), class = "b"),
#     prior(normal(0, 2), class = "Intercept", dpar = "hu"),
#     prior(normal(0, 1.5), class = "b", dpar = "hu"),
#     prior(normal(0, 2), class = "Intercept", dpar = "shape"),
#     prior(normal(0, 1), class = "b", dpar = "shape"),
#     prior(lkj(2), class = "cor") # weaker prior on genetic correlation
#   ),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.99, max_treedepth = 15)
# )
# 
# # Remove region random slope and shape to simplify model and improve convergence
# # CONSERVATIVE MODEL - ONLY GXE FOR TEMP
# model <- brm(
#   bf(
#     Blueness ~ TotalArea + Temperature + Region + 
#       TotalArea:Temperature + Temperature:Region + 
#       (0 + Temperature | gr(animal, cov = Amat)),
#     hu ~ Temperature * Region + (0 + Temperature | gr(animal, cov = Amat)),
#   ),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = hurdle_gamma(),
#   prior = c(
#     prior(normal(0, 2), class = "Intercept"),
#     prior(normal(0, 1), class = "b"),
#     prior(normal(0, 2), class = "Intercept", dpar = "hu"),
#     prior(normal(0, 1.5), class = "b", dpar = "hu")
#   ),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )


#########################################
# BACK TO GAUSSIAN - SIMPLE MODEL
#===============
# REGION EFFECTS
#===============

# Raw blueness
model <- brm(
  Blueness ~ TotalArea + Temperature + Region + 
    TotalArea:Temperature + Temperature:Region + 
    (0 + Temperature:Region | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 9/12000 (0%) divergence, 3/4 chains had E-BFMI < 0.3
# pp check looks ok, 5 Rhats > 1.00
# some estimate errors quite high
# VA OCold: 70 [49, 96], VA OWarm: 21 [10,35]
# VA SCold: 54 [35, 77], VA SWarm: 8 [2, 16]

# z-score blueness
# model <- brm(
#   ZBlueness ~ TotalArea + Temperature + Region + 
#     TotalArea:Temperature + Temperature:Region + 
#     (0 + Temperature:Region | gr(animal, cov = Amat)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )

#logged outcomes - modified data so 7 zeros were replaced with small constant
# model <- brm(
#   LogBlueness ~ TotalArea + Temperature + Region + 
#     TotalArea:Temperature + Temperature:Region + 
#     (0 + Temperature:Region | gr(animal, cov = Amat)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
# # 22/12000 (0%) divergences, 4/4 chains had E-BFMI < 0.3
# # pp check looks good, 7 Rhats > 1.00, 3 = 1.02
# # estimate errors also high

# Total Area
model <- brm(
  TotalArea ~ Temperature + Region + Temperature:Region +
    (0 + Temperature:Region | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 12/12000 (0%) divergences, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, 8 Rhats (& sigma) > 1.00
# VA OCold: 128 [72, 197], VA OWarm: 133 [50,241]
# VA SCold: 123 [58, 208], VA SWarm: 163 [71, 278]

# Global model, no random slopes
model <- brm(
  Blueness ~ TotalArea + Temperature + Region + 
    TotalArea:Temperature + Temperature:Region + 
    (1 | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 29/12000 (0%) divergence, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, but Intercept and sigma Rhat very high 1.05, 1.06
# partially not converged
# VA global: 44.024 [25.464, 63.389]

# area global model
model <- brm(
  TotalArea ~ Temperature + Region + Temperature:Region + 
    (1 | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 1/12000 (0%) divergence, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, Intercept and sigma Rhat 1.01
# VA global: 126.291 [70.123, 200.810]

# region subsets blue
model <- brm(
  Blueness ~ TotalArea + Temperature + TotalArea:Temperature + 
    (0 + Temperature | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# Oland
# 314/12000 (3%) divergence, 4/4 chains had E-BFMI < 0.3
# pp check looks good, 4 Rhats > 1.00, temp random slopes and sigma >1.01
# temp estimate error quite high
# VA OCold: 67 [40, 98], VA OWarm: 21 [7, 41]
# Skane
# 170/12000 (1%) divergence, 4/4 chains had E-BFMI < 0.3
# pp check looks wonky, sigma Rhat 1.14 temp random slopes 1.05-1.09
# VA SCold: 58 [36, 84], VA SWarm: 12 [3, 23]


# region subsets area
model <- brm(
  TotalArea ~ Temperature + (0 + Temperature | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# Oland
# 42/12000 (0%) divergences, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, All 4 Rhats (& sigma) >> 1.00
# VA OCold: 149 [78, 227], VA OWarm: 160 [54, 273]
# Skane
# 21/12000 (0%) divergences, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, 3 Rhats (& sigma) >> 1.00
# VA SCold: 113 [43, 211], VA SWarm: 146 [45, 282]

#=========================
# FAMILY EFFECTS
#=========================

model <- brm(
  Blueness ~ TotalArea + Temperature + Region + 
    TotalArea:Temperature + Temperature:Region + 
    (0 + Temperature | MotherID) +
    (0 + Temperature:Region | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
saveRDS(model, "full_model.rds")
# 3/12000 (0%) divergence, 4/4 chains had E-BFMI < 0.3
# pp check looks good, no > 1.01, 5 Rhats > 1
# some estimate errors quite high - temp and intercept
# BEST SO FAR

# LOO model comparison
model <- brm(
  Blueness ~ TotalArea + Temperature + Region + 
    TotalArea:Temperature + Temperature:Region + 
    (1 | MotherID) + # families constrained to have same slope
    (0 + Temperature:Region | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 0/12000 (0%) divergence
# pp check looks ok, 6 Rhats 1.01
# comparable fit

loo_full <- loo(model_full)
loo_reduced <- loo(model_reduced)
loo_compare(loo_full, loo_reduced)


model <- brm(
  Blueness ~ TotalArea + Temperature + Region + 
    TotalArea:Temperature + Temperature:Region + 
    (0 + Temperature | MotherID) +
    (1 | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 701/12000 (6%) divergence, 2/4 chains had E-BFMI < 0.3
# pp check bad, 2 > 1.01


# mother as a fixed effect - NO REGION
model <- brm(
  Blueness ~ TotalArea + Temperature +
    TotalArea:Temperature + Temperature:MotherID +
    (0 + Temperature | MotherID) +
    (1 | gr(animal, cov = Amat)),
  data = data,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
# 258/12000 (2%) divergence, 167/12000 (1%) hit max treedpeth
# 4/4 chains had E-BFMI < 0.3
# pp check looks ok, 7 Rhats >> 1.01


# model <- brm(
#   Blueness ~ TotalArea + Temperature + Region + 
#     TotalArea:Temperature + Temperature:MotherID + 
#     (0 + Temperature | MotherID) +
#     (0 + Temperature:Region | gr(animal, cov = Amat)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
# 105/12000 (1%) divergence, 11895/12000 (99%) hit max treedepth
# 4/4 chains had E-BFMI < 0.3
# pp check looks good, Rhats horrible
# BAD


# model <- brm(
#   Blueness ~ TotalArea + Temperature +
#     TotalArea:Temperature + Temperature:MotherID +
#     (0 + Temperature | MotherID),
#   data = data,
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
# 2624/12000 (22%) divergence BAD, 50/12000 hit max treedepth
# pp check BAD, Rhats horrible
# DO NOT REMOVE ANIMAL MODEL

# model <- brm(
#   Blueness ~ TotalArea + Temperature +
#     TotalArea:Temperature + Temperature:MotherID +
#     (1 | gr(animal, cov = Amat)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )
# 152/12000 (1%) divergence BAD, 572/12000 (5%) hit max treedepth
# 4/4 chains had E-BFMI < 0.3
# pp check looks ok, 7 Rhats >> 1.01
# not great either



#########################################
# CHECK CONVERGENCE AND FIT
#########################################
model <- readRDS("full_model.rds")

n_div <- sum(subset(nuts_params(model), Parameter == "divergent__")$Value)
cat("Divergences:", n_div, "/",
    4 * 3000, "\n\n")

plot(model)
pp_check(model, ndraws = 100)

summary(model)

# Helper function for formatting posterior summaries
fmt <- function(x) sprintf("%.3f [%.3f, %.3f]",
                           mean(x), quantile(x, 0.025), quantile(x, 0.975))

#################################################

# Extract posterior draws for variance components
draws <- as_draws_df(model)

# fixed effects probability of direction
bayestestR::p_direction(model)

# DO FAMILY SLOPES DIFFER
# Extract posterior draws of the variance components
slope_var <- as_draws_df(model) |>
  mutate(
    slope_sd = sqrt(
      `sd_MotherID__TemperatureWarm26°C`^2 +
        `sd_MotherID__TemperatureCold18°C`^2 -
        2 * `cor_MotherID__TemperatureCold18°C__TemperatureWarm26°C` *
        `sd_MotherID__TemperatureCold18°C` *
        `sd_MotherID__TemperatureWarm26°C`
    )
  )

fmt(slope_var$slope_sd)
# 1.925 [0.338,4.202] - excludes 0, suggests family-level GxE for temperature

# # Always extract animal SD column names for fallback use
# animal_sd_cols <- grep("^sd_animal__", names(draws), value = TRUE)
# 
# # If the column names don't match, try to find them
# if (is.null(va_ocold_sd)) {
#   # brms may name factor levels differently — search for the right columns
#   cat("\nAnimal SD columns found:", paste(animal_sd_cols, collapse = ", "), "\n")
#   
#   if (length(animal_sd_cols) >= 2) {
#     va_cold_sd <- draws[[animal_sd_cols[1]]]
#     va_warm_sd <- draws[[animal_sd_cols[2]]]
#     cat("Using:", animal_sd_cols[1], "as VA_cold_sd\n")
#     cat("Using:", animal_sd_cols[2], "as VA_warm_sd\n")
#   }
# }

#########################################
# RAW GENETIC CORRELATIONS (within region)
# rG(Cold, Warm) from family mean correlations
#########################################

# Family means per temperature × region
family_means <- data %>%
  group_by(MotherID, Region, Temperature) %>%
  summarise(
    mean_blue = mean(Blueness, na.rm = TRUE),
    n_obs     = n(),
    .groups   = "drop"
  )

# Pivot so each temperature is a column
family_wide <- family_means %>%
  pivot_wider(
    id_cols    = c(MotherID, Region),
    names_from = Temperature,
    values_from = mean_blue
  )

# Pearson correlation of family means between temperatures, by region
rG_raw <- family_wide %>%
  group_by(Region) %>%
  summarise(
    n_families = n(),
    rG         = cor(`Cold (18°C)`, `Warm (26°C)`, use = "pairwise.complete.obs"),
    .groups    = "drop"
  )
print(rG_raw)
# Öland: rG = 0.653 (n=25), Skåne: rG = 0.825 (n=17)

# Bootstrap 95% CIs (families are the sampling unit)
set.seed(BAYES_SEED)
n_boot <- 5000

boot_rG <- function(df, region_name) {
  d <- filter(df, Region == region_name)
  replicate(n_boot, {
    s <- d[sample(nrow(d), replace = TRUE), ]
    cor(s[[3]], s[[4]], use = "pairwise.complete.obs")  # cols 3,4 are the two temps
  })
}

boot_oland <- boot_rG(family_wide, "Öland")
boot_skane <- boot_rG(family_wide, "Skåne")

cat("rG Öland (Cold↔Warm):", fmt(boot_oland), "\n")
cat("rG Skåne (Cold↔Warm):", fmt(boot_skane), "\n")
# rG Öland (Cold↔Warm): 0.659 [0.413, 0.852]
# rG Skåne (Cold↔Warm): 0.810 [0.448, 0.952]

# Genetic correlations between family effects across temperatures
mom_cors <- draws$`cor_MotherID__TemperatureCold18°C__TemperatureWarm26°C`
fmt(mom_cors)
# 0.528 [-0.813,0.993] - CI includes 0, but mean suggests some positive 
# correlation between family effects across temperatures, so not strong GxE 
# at family level

# Extract all 6 genetic correlations from the Temperature:Region animal random effect
rG_draws <- as_draws_df(model) |>
  select(starts_with("cor_animal__"))

# Rename to readable labels (check column names first with: names(rG_draws))
rG_summary <- rG_draws |>
  rename(
      rG_OCold_OWarm  =
  `cor_animal__TemperatureCold18°C:RegionÖland__TemperatureWarm26°C:RegionÖland`,
      rG_SCold_SWarm  =
  `cor_animal__TemperatureCold18°C:RegionSkåne__TemperatureWarm26°C:RegionSkåne`,
      
      rG_OCold_SCold  =
  `cor_animal__TemperatureCold18°C:RegionÖland__TemperatureCold18°C:RegionSkåne`,
      rG_OWarm_SWarm  =
  `cor_animal__TemperatureWarm26°C:RegionÖland__TemperatureWarm26°C:RegionSkåne`,
      
      rG_OCold_SWarm  =
  `cor_animal__TemperatureCold18°C:RegionÖland__TemperatureWarm26°C:RegionSkåne`,
      rG_OWarm_SCold  =
  `cor_animal__TemperatureWarm26°C:RegionÖland__TemperatureCold18°C:RegionSkåne`
    ) |>
    pivot_longer(everything(), names_to = "Correlation", values_to = "rG") |>
    group_by(Correlation) |>
    summarise(
      mean  = mean(rG),
      lower = quantile(rG, 0.025),
      upper = quantile(rG, 0.975),
      p_lt1 = mean(rG < 1),    # P(GxE exists)
      p_lt0.8 = mean(rG < 0.8) # P(biologically meaningful GxE)
    )
print(rG_summary)


#########################################
# EXTRACT GENETIC PARAMETERS: VA
#########################################

# Additive genetic SDs by temperature and region (random slope SDs) 
va_ocold_sd <- draws$`sd_animal__TemperatureCold18°C:RegionÖland`
va_scold_sd <- draws$`sd_animal__TemperatureCold18°C:RegionSkåne`
va_owarm_sd <- draws$`sd_animal__TemperatureWarm26°C:RegionÖland`
va_swarm_sd <- draws$`sd_animal__TemperatureWarm26°C:RegionSkåne`

# va_sd <- draws$'sd_animal__Intercept'
# va_global <- va_sd^2
# fmt(va_global)
# 
# # region subsets
# va_cold_sd <- draws$`sd_animal__TemperatureCold18°C`
# va_warm_sd <- draws$`sd_animal__TemperatureWarm26°C`
# va_cold <- va_cold_sd^2
# va_warm <- va_warm_sd^2
# fmt(va_cold)
# fmt(va_warm)

# mother SDs
va_cold_mother_sd <- draws$`sd_MotherID__TemperatureCold18°C`
va_warm_mother_sd <- draws$`sd_MotherID__TemperatureWarm26°C`
va_cold_mother <- va_cold_mother_sd^2
va_warm_mother <- va_warm_mother_sd^2
fmt(va_cold_mother)
fmt(va_warm_mother)


if (!is.null(va_ocold_sd)) {
  va_ocold <- va_ocold_sd^2
  va_scold <- va_scold_sd^2
  va_owarm <- va_owarm_sd^2
  va_swarm <- va_swarm_sd^2
  
  cat("  VA at Cold (18°C):Region Öland:      ", fmt(va_ocold), "\n")
  cat("  VA at Warm (26°C):Region Öland:      ", fmt(va_owarm), "\n")
  cat("  VA at Cold (18°C):Region Skåne:      ", fmt(va_scold), "\n")
  cat("  VA at Warm (26°C):Region Skåne:      ", fmt(va_swarm), "\n")
  cat("  sqrt(VA) ÖlandCold:          ", fmt(va_ocold_sd), "\n")
  cat("  sqrt(VA) ÖlandWarm:          ", fmt(va_owarm_sd), "\n")
  cat("  sqrt(VA) SkåneCold:          ", fmt(va_scold_sd), "\n")
  cat("  sqrt(VA) SkåneWarm:          ", fmt(va_swarm_sd), "\n")
  
  if (exists("rG")) {
    cat("  Genetic correlation rG: ", fmt(rG), "\n")
    cat("  P(rG < 1):              ", sprintf("%.4f", mean(rG < 1)), "\n")
    cat("  P(rG < 0.8):            ", sprintf("%.4f", mean(rG < 0.8)), "\n")
    cat("\n  rG < 1 indicates genotype × environment interaction.\n")
    cat("  rG < 0.8 is often considered biologically meaningful G×E.\n")
  }
}

va_raw <- tryCatch({
  tibble(
    `Öland Cold` = draws$`sd_animal__TemperatureCold18°C:RegionÖland`^2,
    `Öland Warm` = draws$`sd_animal__TemperatureWarm26°C:RegionÖland`^2,
    `Skåne Cold` = draws$`sd_animal__TemperatureCold18°C:RegionSkåne`^2,
    `Skåne Warm` = draws$`sd_animal__TemperatureWarm26°C:RegionSkåne`^2
  )
}, error = function(e) NULL)

if (is.null(va_raw) || any(sapply(va_raw, is.null))) {
  cat("Direct column names didn't match; using positional extraction.\n")
  cat("Columns found (check order matches OCold/OWarm/SCold/SWarm):\n")
  print(animal_sd_cols)
  va_raw <- tibble(
    `Öland Cold` = draws[[animal_sd_cols[1]]]^2,
    `Öland Warm` = draws[[animal_sd_cols[2]]]^2,
    `Skåne Cold` = draws[[animal_sd_cols[3]]]^2,
    `Skåne Warm` = draws[[animal_sd_cols[4]]]^2
  )
}

va_draws <- va_raw %>%
  pivot_longer(everything(), names_to = "Group", values_to = "VA") %>%
  mutate(Group = factor(Group, levels = c("Öland Cold", "Öland Warm",
                                          "Skåne Cold", "Skåne Warm")))

# Posterior mean and 95% CrI for each group
va_summary <- va_draws %>%
  group_by(Group) %>%
  summarise(
    post_mean = mean(VA),
    lo  = quantile(VA, 0.025),
    hi  = quantile(VA, 0.975),
    .groups = "drop"
  )

print(va_summary)

# Plot
va_violin <- ggplot(va_draws, aes(x = Group, y = VA, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
  geom_linerange(
    data = va_summary,
    aes(x = Group, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    linetype = "dashed", colour = "grey30"
  ) +
  geom_errorbar(
    data = va_summary,
    aes(x = Group, ymin = post_mean, ymax = post_mean),
    inherit.aes = FALSE,
    width = 0.25, linewidth = 0.8
  ) +
  geom_vline(xintercept = 2.5, colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = c("Öland Cold" = "#3B7DD8", 
                               "Öland Warm" = "#E8712A", 
                               "Skåne Cold" = "#3B7DD8", 
                               "Skåne Warm" = "#E8712A")) +
  labs(x = NULL,
       y = expression("Additive genetic variance (" * V[A] * ")"),
       title = NULL) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

va_violin

ggsave(file.path(PLOT_DIR,"va_violin.png"), va_violin, width = 7, height = 5, dpi = 300)


#########################################
# PHENOTYPIC VARIANCE (VP) AND HERITABILITY (h²)
#########################################

# Observed VP per environment
vp_summary <- blueFdata %>%
  group_by(region_label, temp_label) %>%
  summarise(
    n       = n(),
    VP      = var(avg_blue_mm, na.rm = TRUE),
    SD      = sd(avg_blue_mm, na.rm = TRUE),
    .groups = "drop"
  )

print(vp_summary)

# Save VP as individual objects
# Adjust the string values below if they don't match your temp_label/region_label levels
vp_ocold <- vp_summary$VP[vp_summary$region_label == "Öland" & grepl("Cold", vp_summary$temp_label)]
vp_owarm <- vp_summary$VP[vp_summary$region_label == "Öland" & grepl("Warm", vp_summary$temp_label)]
vp_scold <- vp_summary$VP[vp_summary$region_label == "Skåne" & grepl("Cold", vp_summary$temp_label)]
vp_swarm <- vp_summary$VP[vp_summary$region_label == "Skåne" & grepl("Warm", vp_summary$temp_label)]

sd_ocold <- sqrt(vp_ocold)
sd_owarm <- sqrt(vp_owarm)
sd_scold <- sqrt(vp_scold)
sd_swarm <- sqrt(vp_swarm)

cat("VP Öland Cold:", vp_ocold, sd_ocold, "\n")
cat("VP Öland Warm:", vp_owarm, sd_owarm, "\n")
cat("VP Skåne Cold:", vp_scold, sd_scold, "\n")
cat("VP Skåne Warm:", vp_swarm, sd_swarm, "\n")

# # mother:temp VPs
# vp_cold <- vp_summary$VP[grepl("Cold", vp_summary$temp_label)]
# vp_warm <- vp_summary$VP[grepl("Warm", vp_summary$temp_label)]
# 
# sd_cold <- sqrt(vp_cold)
# sd_warm <- sqrt(vp_warm)
# 
# cat("VP Cold:", vp_cold, sd_cold, "\n")
# cat("VP Warm:", vp_warm, sd_warm, "\n")

# h² = VA/VP — divide each posterior draw of VA by the scalar VP
h2_ocold <- va_ocold / vp_ocold
h2_owarm <- va_owarm / vp_owarm
h2_scold <- va_scold / vp_scold
h2_swarm <- va_swarm / vp_swarm

cat("\nHeritability (h² = VA/VP):\n")
cat("  h² Öland Cold:", fmt(h2_ocold), "\n")
cat("  h² Öland Warm:", fmt(h2_owarm), "\n")
cat("  h² Skåne Cold:", fmt(h2_scold), "\n")
cat("  h² Skåne Warm:", fmt(h2_swarm), "\n")

# # mother:temp h²
# h2_cold <- va_cold_mother / vp_cold
# h2_warm <- va_warm_mother / vp_warm
# 
# cat("  h² Cold:", fmt(h2_cold), "\n")
# cat("  h² Warm:", fmt(h2_warm), "\n")

#============================================

sigma_draws <- draws$sigma

vp_ocold_draws <- va_ocold + sigma_draws^2
vp_owarm_draws <- va_owarm + sigma_draws^2
vp_scold_draws <- va_scold + sigma_draws^2
vp_swarm_draws <- va_swarm + sigma_draws^2

h2_ocold <- va_ocold / vp_ocold_draws
h2_owarm <- va_owarm / vp_owarm_draws
h2_scold <- va_scold / vp_scold_draws
h2_swarm <- va_swarm / vp_swarm_draws


# Posterior distribution of h² for plotting
h2_draws <- tibble(
  `Öland Cold` = h2_ocold,
  `Öland Warm` = h2_owarm,
  `Skåne Cold` = h2_scold,
  `Skåne Warm` = h2_swarm
) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "h2") %>%
  mutate(Group = factor(Group, levels = c("Öland Cold", "Öland Warm",
                                          "Skåne Cold", "Skåne Warm")))

h2_summary <- h2_draws %>%
  group_by(Group) %>%
  summarise(
    post_mean = mean(h2),
    lo  = quantile(h2, 0.025),
    hi  = quantile(h2, 0.975),
    .groups = "drop"
  )

print(h2_summary)

h2_violin <- ggplot(h2_draws, aes(x = Group, y = h2, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
  geom_linerange(
    data = h2_summary,
    aes(x = Group, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    linetype = "dashed", colour = "grey30"
  ) +
  geom_errorbar(
    data = h2_summary,
    aes(x = Group, ymin = post_mean, ymax = post_mean),
    inherit.aes = FALSE,
    width = 0.25, linewidth = 0.8
  ) +
  geom_vline(xintercept = 2.5, colour = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = c("Öland Cold" = "#3B7DD8", 
                               "Öland Warm" = "#E8712A", 
                               "Skåne Cold" = "#3B7DD8", 
                               "Skåne Warm" = "#E8712A")) +
  labs(x = NULL,
       y = expression("Heritability (" * h^2 * ")"),
       title = NULL) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

h2_violin

ggsave(file.path(PLOT_DIR,"h2_violin.png"), h2_violin, width = 7, height = 5, dpi = 300)


################################
# EVOLVABILITY
################################

# VA/mean blueness per environment
env_means <- data %>%
  group_by(Region, Temperature) %>%
  summarise(Mean = mean(Blueness),.groups = "drop") %>%
  arrange(Region, Temperature) %>%
  as.data.frame()

mean_ocold <- env_means[1,3]
mean_owarm <- env_means[2,3]
mean_scold <- env_means[3,3]
mean_swarm <- env_means[4,3]

# Evolvability = VA / mean^2
e_ocold <- va_ocold / (mean_ocold^2)
e_scold <- va_scold / (mean_scold^2)
e_owarm <- va_owarm / (mean_owarm^2)
e_swarm <- va_swarm / (mean_swarm^2)

cat("\nEvolvability (e = VA/mean^2):\n")
cat("  Evolvability Öland Cold:", fmt(e_ocold), "\n")
cat("  Evolvability Öland Warm:", fmt(e_owarm), "\n")
cat("  Evolvability Skåne Cold:", fmt(e_scold), "\n")
cat("  Evolvability Skåne Warm:", fmt(e_swarm), "\n")

#============================================

e_draws <- tibble(
  `Öland Cold` = e_ocold,
  `Öland Warm` = e_owarm,
  `Skåne Cold` = e_scold,
  `Skåne Warm` = e_swarm
) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "e") %>%
  mutate(Group = factor(Group, levels = c("Öland Cold", "Öland Warm",
                                          "Skåne Cold", "Skåne Warm")))

e_summary <- e_draws %>%
  group_by(Group) %>%
  summarise(
    post_mean = mean(e),
    lo  = quantile(e, 0.025),
    hi  = quantile(e, 0.975),
    .groups = "drop"
  )

print(e_summary)

e_violin <- ggplot(e_draws, aes(x = Group, y = e, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
  geom_linerange(
    data = e_summary,
    aes(x = Group, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    linetype = "dashed", colour = "grey30"
  ) +
  geom_errorbar(
    data = e_summary,
    aes(x = Group, ymin = post_mean, ymax = post_mean),
    inherit.aes = FALSE,
    width = 0.25, linewidth = 0.8
  ) +
  geom_vline(xintercept = 2.5, colour = "grey40", linewidth = 0.5) +
  scale_y_continuous(limits = c(0,9)) +
  scale_fill_manual(values = c("Öland Cold" = "#3B7DD8", 
                               "Öland Warm" = "#E8712A", 
                               "Skåne Cold" = "#3B7DD8", 
                               "Skåne Warm" = "#E8712A")) +
  labs(x = NULL,
       y = expression("Evolvability (e)"),
       title = NULL) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

e_violin

ggsave(file.path(PLOT_DIR,"e_violin.png"), e_violin, width = 7, height = 5, dpi = 300)

combined <- grid.arrange(va_violin, h2_violin, e_violin, ncol = 3)
ggsave(file.path(PLOT_DIR,"va_h_e_violin.png"), combined, width = 16, height = 5, dpi = 300)


####################################
# MOTHER GXE PREDICTIVE RXN NORMS
####################################

slopes <- spread_draws(model, r_MotherID[MotherID, Temperature]) |>
  pivot_wider(names_from = Temperature, values_from = r_MotherID) |>
  mutate(slope = `TemperatureWarm26°C` - `TemperatureCold18°C`) |>
  group_by(MotherID) |>
  summarise(
    slope_mean  = mean(slope),
    slope_lower = quantile(slope, 0.025),
    slope_upper = quantile(slope, 0.975)
  )
print(slopes, n=Inf) # ALL FAMILY SLOPES

mother_region <- data |>
  distinct(MotherID, Region)

temp_labels <- c(
  "TemperatureCold18°C" = "Cold (18°C)",
  "TemperatureWarm26°C" = "Warm (26°C)"
)

# predictive rxn norm plot - deviations from average expected blueness
# accounting for all fixed effect and random effects at zero
pred_dev_norms <- spread_draws(model, r_MotherID[MotherID, Temperature]) |>
  group_by(MotherID, Temperature) |>
  summarise(mean_re = mean(r_MotherID), .groups = "drop") |>
  left_join(mother_region, by = "MotherID") |>
  ggplot(aes(x = Temperature, y = mean_re, group = MotherID, colour = Region)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_REGION) +
  scale_x_discrete(labels = temp_labels, expand = expansion(mult = 0.1)) +
  labs(x = "Temperature treatment",
       y = "Family-level deviation",
       colour = "Region") +
  theme_classic() +
  theme(text = element_text(size = 13),
        legend.position = "right",
        legend.background = element_rect(colour = "grey80"))

ggsave(file.path(PLOT_DIR, "pred_dev_norms.png"), pred_dev_norms, width = 6, height = 5, dpi = 200)


# Build prediction grid - one row per MotherID × Temperature combination
# Hold TotalArea at its mean, use modal Region per mother
pred_grid <- data |>
  group_by(MotherID) |>
  mutate(Region = names(sort(table(Region), decreasing = TRUE))[1]) |>
  ungroup() |>
  distinct(MotherID, Temperature, Region) |>
  mutate(
    TotalArea = mean(data$TotalArea),
    animal    = NA  # exclude animal/genetic effect
  )

# Draw posterior expected values per MotherID × Temperature
pred_rxn_norms <- pred_grid |>
  add_epred_draws(
    model,
    re_formula = ~ (0 + Temperature | MotherID),  # include only MotherID random effects
    allow_new_levels = FALSE
  ) |>
  group_by(MotherID, Temperature) |>
  summarise(mean_pred = mean(.epred), .groups = "drop") |>
  left_join(mother_region, by = "MotherID") |>
  ggplot(aes(x = Temperature, y = mean_pred, group = MotherID, colour = Region)) +
  geom_line(alpha = 0.45, linewidth = 0.5) +
  geom_point(alpha = 0.45, size = 1.5) +
  scale_colour_manual(values = PAL_REGION) +
  scale_x_discrete(labels = temp_labels, expand = expansion(mult = 0.1)) +
  labs(x = "Temperature treatment",
       y = "Predicted blue wing area (mm²)",
       colour = "Region") +
  theme_classic() +
  theme(text = element_text(size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"))

ggsave(file.path(PLOT_DIR, "pred_rxn_norms.png"), pred_rxn_norms, width = 6, height = 5, dpi = 200)


####################################
# REGIONAL PREDICTIVE RXN NORMS
####################################

region_grid <- expand.grid(
  Temperature = levels(data$Temperature),
  Region      = levels(data$Region)
) |>
  mutate(
    TotalArea = mean(data$TotalArea),
    animal    = NA,
    MotherID  = NA
  )

region_epred <- posterior_epred(model, newdata = region_grid,
                                re_formula = NA, allow_new_levels = TRUE)

region_pred_df <- region_grid |>
  mutate(
    post_mean = apply(region_epred, 2, mean),
    lo        = apply(region_epred, 2, quantile, 0.025),
    hi        = apply(region_epred, 2, quantile, 0.975)
  )

pred_region_norms <- ggplot(region_pred_df,
                            aes(x = Temperature, y = post_mean,
                                colour = Region, group = Region)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, linewidth = 0.8) +
  scale_colour_manual(values = PAL_REGION) +
  scale_x_discrete(labels = temp_labels, expand = expansion(mult = 0.15)) +
  labs(x = "Temperature treatment",
       y = expression("Predicted blue wing area (mm"^2*")"),
       colour = "Region") +
  theme_classic() +
  theme(text = element_text(size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"))

pred_region_norms

ggsave(file.path(PLOT_DIR, "pred_region_norms.png"), pred_region_norms,
       width = 6, height = 5, dpi = 300)


#########################################
# VARIANCE PARTITIONING
#########################################

# --- R² via bayes_R2() ---

# Marginal R²: fixed effects only (TotalArea, Temperature, Region, interactions)
r2_marginal <- bayes_R2(model, re_formula = NA, summary = FALSE)

# Conditional R²: fixed + all random effects
r2_conditional <- bayes_R2(model, summary = FALSE)

# Random effects R² (combined MotherID + animal)
r2_random <- r2_conditional[, "R2"] - r2_marginal[, "R2"]

cat("Marginal R² (fixed effects):    ", fmt(r2_marginal[, "R2"]), "\n")
cat("Conditional R² (fixed + random):", fmt(r2_conditional[, "R2"]), "\n")
cat("Random effects R² (combined):   ", fmt(r2_random), "\n")

# --- Manual variance decomposition ---
# Partition random effects into MotherID vs animal (VA) vs residual

# Residual variance
v_resid <- draws$sigma^2

# Average VA across 4 environments (or report per-environment)
va_mean <- (draws$`sd_animal__TemperatureCold18°C:RegionÖland`^2 +
              draws$`sd_animal__TemperatureWarm26°C:RegionÖland`^2 +
              draws$`sd_animal__TemperatureCold18°C:RegionSkåne`^2 +
              draws$`sd_animal__TemperatureWarm26°C:RegionSkåne`^2) / 4

# Average MotherID variance across 2 temperatures
v_mother <- (draws$`sd_MotherID__TemperatureCold18°C`^2 +
               draws$`sd_MotherID__TemperatureWarm26°C`^2) / 2

# Fixed effect variance: var of fitted values with no random effects
fe_pred  <- posterior_epred(model, re_formula = NA)        # draws × obs matrix
v_fixed  <- apply(fe_pred, 1, var)                         # one value per draw

# Total VP (sum of all components)
v_total <- v_fixed + v_mother + va_mean + v_resid

# Proportions of variance
vp_df <- tibble(
  Component = c("Fixed effects", "MotherID (maternal)", "Animal (VA)", "Residual"),
  Proportion = list(v_fixed / v_total, v_mother / v_total,
                    va_mean / v_total,  v_resid / v_total)
) |>
  mutate(
    mean  = sapply(Proportion, mean),
    lower = sapply(Proportion, quantile, 0.025),
    upper = sapply(Proportion, quantile, 0.975)
  ) |>
  select(-Proportion)

print(vp_df)

#########################################
# VARIANCE DECOMPOSITION BY TEMPERATURE
# Matching analysis.R fig5 style
#########################################

# Indices by temperature (check levels(data$Temperature) if these fail)
cold_idx <- which(grepl("Cold", data$Temperature))
warm_idx <- which(grepl("Warm", data$Temperature))

# Fixed effect variance per temperature (subset fe_pred columns by temperature)
# fe_pred computed above: posterior_epred(model, re_formula = NA)
v_fixed_cold <- apply(fe_pred[, cold_idx], 1, var)
v_fixed_warm <- apply(fe_pred[, warm_idx], 1, var)

# Average VA across regions per temperature
va_cold_draws <- (va_ocold + va_scold) / 2
va_warm_draws <- (va_owarm + va_swarm) / 2

# Posterior means of each component per temperature
comp_cold <- c(
  Residual     = mean(v_resid),
  Maternal     = mean(va_cold_mother),
  `Animal (VA)`= mean(va_cold_draws),
  Fixed        = mean(v_fixed_cold)
)
comp_warm <- c(
  Residual     = mean(v_resid),
  Maternal     = mean(va_warm_mother),
  `Animal (VA)`= mean(va_warm_draws),
  Fixed        = mean(v_fixed_warm)
)

comp_names <- c("Residual", "Maternal", "Animal (VA)", "Fixed")

# Build bar data frame (actual variance values, matching analysis.R fig5)
bar_df <- data.frame(
  temp      = rep(c("Cold (18°C)", "Warm (26°C)"), each = length(comp_names)),
  component = rep(comp_names, times = 2),
  value     = c(comp_cold, comp_warm),
  stringsAsFactors = FALSE
) |>
  mutate(
    temp      = factor(temp, levels = c("Cold (18°C)", "Warm (26°C)")),
    component = factor(component, levels = comp_names)
  )

label_df <- bar_df |>
  group_by(temp) |>
  mutate(pct = sprintf("%.1f%%", value / sum(value) * 100)) |>
  ungroup()

vp_temp_bar <- ggplot(bar_df, aes(x = temp, y = value, fill = component)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(
    data = label_df,
    aes(label = pct),
    position = position_stack(vjust = 0.5),
    colour = "white", fontface = "bold", size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Residual"     = "#AAAAAA",
      "Maternal"     = "#E8A838",
      "Animal (VA)"  = "#56A156",
      "Fixed"        = "#5A6EA1"
    )
  ) +
  labs(
    x     = NULL,
    y     = expression("Variance in blue wing area"),
    fill  = ""
  ) +
  theme_classic() +
  theme(
    legend.position   = "inside",
    legend.position.inside = c(0.78, 0.82),
    legend.title = element_blank(),
    legend.background = element_rect(colour = "grey80"),
    legend.text       = element_text(size = 9),
    axis.title        = element_text(size = 13),
    axis.text         = element_text(size = 12)
  )

vp_temp_bar

ggsave(file.path(PLOT_DIR, "var_decomp.png"), vp_temp_bar,
       width = 5, height = 5, dpi = 300)

#########################################
# VARIANCE DECOMPOSITION DOT-AND-RANGE
#########################################

vp_draws_long <- bind_rows(
  tibble(temp = "Cold (18°C)", component = "Fixed",       value = v_fixed_cold),
  tibble(temp = "Cold (18°C)", component = "Maternal",    value = va_cold_mother),
  tibble(temp = "Cold (18°C)", component = "Animal (VA)", value = va_cold_draws),
  tibble(temp = "Cold (18°C)", component = "Residual",    value = v_resid),
  tibble(temp = "Warm (26°C)", component = "Fixed",       value = v_fixed_warm),
  tibble(temp = "Warm (26°C)", component = "Maternal",    value = va_warm_mother),
  tibble(temp = "Warm (26°C)", component = "Animal (VA)", value = va_warm_draws),
  tibble(temp = "Warm (26°C)", component = "Residual",    value = v_resid)
) |>
  mutate(
    temp      = factor(temp, levels = c("Cold (18°C)", "Warm (26°C)")),
    component = factor(component, levels = rev(comp_names))
  )

vp_range_summary <- vp_draws_long |>
  group_by(temp, component) |>
  summarise(
    post_mean = mean(value),
    lo        = quantile(value, 0.025),
    hi        = quantile(value, 0.975),
    .groups   = "drop"
  )

vp_range <- ggplot(vp_range_summary,
                   aes(x = post_mean, y = component, colour = temp,
                       xmin = lo, xmax = hi)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.5) +
  scale_colour_manual(values = c("Cold (18°C)" = "#3B7DD8", "Warm (26°C)" = "#E8712A")) +
  labs(
    x      = expression("Variance in blue wing area"),
    y      = NULL,
    colour = NULL
  ) +
  theme_classic() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.80, 0.85),
    legend.background      = element_rect(colour = "grey80"),
    axis.title             = element_text(size = 13),
    axis.text              = element_text(size = 12)
  )

vp_range

ggsave(file.path(PLOT_DIR, "var_ranges.png"), vp_range, width = 6, height = 4, dpi = 300)

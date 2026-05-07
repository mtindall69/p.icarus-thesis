# Polyommatus Animal Model - GxE with brms

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(tidyverse,brms,nadiv,emmeans,tidybayes,broom,broom.mixed,
               patchwork,ggh4x,ggtext,MetBrewer,dplyr)

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

#subset data
blueFdata <- subset(bluesum, sex=="F")
blueMdata <- subset(bluesum, sex=="M")


n_families <- length(unique(blueFdata$motherID)) # 42 families
n_offspring <- length(blueFdata$offspringID) # 327 daughters

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
         Blueness = avg_blue_mm)#, LogBlueness = logblue)

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
# Forces full-sibs- gives conservative estimate of heritability


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
#########################################

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
# 
# # Now for total area to get correlation of blueness va to area va
# model <- brm(
#   TotalArea ~ Temperature + Region + Temperature:Region + 
#     (0 + Temperature:Region | gr(animal, cov = Amat)),
#   data = data,
#   data2 = list(Amat = Amat),
#   family = gaussian(),
#   chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
#)
# 12/12000 (0%) divergences, 4/4 chains had E-BFMI < 0.3
# pp check looks ok, 8 Rhats (& sigma) > 1.00
# VA OCold: 128 [72, 197], VA OWarm: 133 [50,241]
# VA SCold: 123 [58, 208], VA SWarm: 163 [71, 278]


#########################################
# CHECK CONVERGENCE AND FIT
#########################################

n_div <- sum(subset(nuts_params(model), Parameter == "divergent__")$Value)
cat("Divergences:", n_div, "/",
    4 * 3000, "\n\n")

plot(model)
pp_check(model, ndraws = 100)

#########################################
# EXTRACT GENETIC PARAMETERS
#########################################

summary(model)
vc <- VarCorr(model)
print(vc)

# Helper function for formatting posterior summaries
fmt <- function(x) sprintf("%.3f [%.3f, %.3f]",
                           mean(x), quantile(x, 0.03), quantile(x, 0.97))

# Extract posterior draws for variance components
draws <- as_draws_df(model)

# Additive genetic SDs by temperature and region (random slope SDs) 
va_ocold_sd <- draws$`sd_animal__TemperatureCold18°C:RegionÖland`
va_scold_sd <- draws$`sd_animal__TemperatureCold18°C:RegionSkåne`
va_owarm_sd <- draws$`sd_animal__TemperatureWarm26°C:RegionÖland`
va_swarm_sd <- draws$`sd_animal__TemperatureWarm26°C:RegionSkåne`

# Always extract animal SD column names for fallback use
animal_sd_cols <- grep("^sd_animal__", names(draws), value = TRUE)

# If the column names don't match, try to find them
if (is.null(va_ocold_sd)) {
  # brms may name factor levels differently — search for the right columns
  cat("\nAnimal SD columns found:", paste(animal_sd_cols, collapse = ", "), "\n")
  
  if (length(animal_sd_cols) >= 2) {
    va_cold_sd <- draws[[animal_sd_cols[1]]]
    va_warm_sd <- draws[[animal_sd_cols[2]]]
    cat("Using:", animal_sd_cols[1], "as VA_cold_sd\n")
    cat("Using:", animal_sd_cols[2], "as VA_warm_sd\n")
  }
}

# Genetic correlation
cor_cols <- grep("^cor_animal__", names(draws), value = TRUE)
if (length(cor_cols) > 0) {
  rG <- draws[[cor_cols[1]]]
  cat("\nUsing:", cor_cols[1], "as genetic correlation\n")
}

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
    .groups = "drop"
  )

print(vp_summary)

# Save VP as individual objects
# Adjust the string values below if they don't match your temp_label/region_label levels
vp_ocold <- vp_summary$VP[vp_summary$region_label == "Öland" & grepl("Cold", vp_summary$temp_label)]
vp_owarm <- vp_summary$VP[vp_summary$region_label == "Öland" & grepl("Warm", vp_summary$temp_label)]
vp_scold <- vp_summary$VP[vp_summary$region_label == "Skåne" & grepl("Cold", vp_summary$temp_label)]
vp_swarm <- vp_summary$VP[vp_summary$region_label == "Skåne" & grepl("Warm", vp_summary$temp_label)]

cat("VP Öland Cold:", vp_ocold, "\n")
cat("VP Öland Warm:", vp_owarm, "\n")
cat("VP Skåne Cold:", vp_scold, "\n")
cat("VP Skåne Warm:", vp_swarm, "\n")

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

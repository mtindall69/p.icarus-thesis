# Polyommatus Animal Model - GxE with MCMCglmm

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm)

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


#==============================================================
# REAL DATA
#==============================================================

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

##############################################
# 3. Define Priors
##############################################
prior <- list(
  G = list(
    G1 = list(V = 1, nu = 0.002),  # Animal (VA)
    G2 = list(V = 1, nu = 0.002),  # MotherID (maternal)
    G3 = list(V = diag(2), nu = 0.002), # GxE for Temperature
    G4 = list(V = diag(2), nu = 0.002)
  ),
  R = list(V = 1, nu = 0.002)      # Residual
)

##############################################
# 4. Fit Animal Model
##############################################
model <- MCMCglmm(
  Blueness ~ TotalArea + Temperature + Region + TotalArea:Temperature + Temperature:Region,
  random = ~ animal + MotherID + idh(Temperature):animal + idh(Region):animal,
  pedigree = pedigree,
  data = data,
  family = "gaussian",
  prior = prior,
  nitt = 130000, burnin = 30000, thin = 100
)

summary(model) # intercept and temp returning nonsignificant

# CHECK FOR COVERGENCE
plot(model$Sol)
plot(model$VCV) # does not converge

##############################################
# 5. Extract Genetic Parameters
##############################################

VA <- model$VCV[,"animal"]
Vres <- model$VCV[,"units"]
Vm <- model$VCV[,"MotherID"]

# Heritability
h2 <- VA / (VA + Vres + Vm)

# Evolvability
meanBlueness <- mean(data$Blueness)
evolvability <- VA / (meanBlueness^2)

cat("Posterior mean VA:", mean(VA), "\n")
cat("Posterior mean h2:", mean(h2), "\n")
cat("Posterior mean evolvability:", mean(evolvability), "\n")

# --- Helper: full summary for any posterior chain ---
posterior_summary <- function(chain, param_name) {
  hpd <- HPDinterval(chain, prob = 0.95)
  # pMCMC: proportion of posterior on the side of zero, x2 (two-tailed)
  # For variance components (always > 0), this will always be ~1; 
  p_lower <- mean(chain <= 0)
  p_upper <- mean(chain >= 0)
  
  data.frame(
    Parameter  = param_name,
    Post.Mean  = mean(chain),
    Post.Mode  = as.numeric(posterior.mode(chain)),
    SD         = sd(chain),           # posterior SD = Bayesian analog of SE
    CI.lower   = hpd[1, "lower"],
    CI.upper   = hpd[1, "upper"]
  )
}

# --- Build results table ---
results <- rbind(
  posterior_summary(VA,          "VA"),
  posterior_summary(h2,          "h2"),
  posterior_summary(evolvability,"Evolvability (e)")
)

print(results, digits = 4)

##############################################
# 6. GxE Interaction
##############################################
# idh(Temperature):animal gives separate VA for cold and warm
VA_cold <- model$VCV[,"TemperatureCold (18°C).animal"]
VA_warm <- model$VCV[,"TemperatureWarm (26°C).animal"]

VA_oland <- model$VCV[,"RegionSkåne.animal"]
VA_skane <- model$VCV[,"RegionÖland.animal"]

cat("VA cold:", mean(VA_cold), "\n")
cat("VA warm:", mean(VA_warm), "\n")

results_GxE <- rbind(
  posterior_summary(VA_cold, "VA_cold"),
  posterior_summary(VA_warm, "VA_warm"),
  posterior_summary(VA_oland, "VA_Oland"),
  posterior_summary(VA_skane, "VA_Skane")
)

print(results_GxE, digits = 4)


##############################################
# 6b. OPTIONAL: pMCMC for fixed effects (already in summary)
##############################################
# summary(model)$solutions gives: post.mean, l-95% CI, u-95% CI, eff.samp, pMCMC
# for all fixed effects automatically.
print(summary(model)$solutions)


#EXTENDED PLOT FOR VISUALIZATIONS

##############################################
# 7. Posterior Density Plots
##############################################

# Function to plot posterior distributions
plot_posterior <- function(samples, title, xlab){
  df <- data.frame(value = samples)
  ggplot(df, aes(x = samples)) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    theme_bw() +
    labs(title = title, x = xlab, y = "Density") +
    geom_vline(xintercept = mean(samples), color = "red", linetype = "dashed") +
    annotate("text", x = mean(samples), y = 0, 
             label = paste0("Mean = ", round(mean(samples), 3)), 
             hjust = -0.1, color = "red")
}


# Plot VA
p <- plot_posterior(VA, "Posterior of Additive Genetic Variance (VA)", "VA")
p
#ggsave(file.path("plots", "posteriorVA.png"), p)

# Plot h²
p <- plot_posterior(h2, "Posterior of Heritability (h²)", "h²")
p
#ggsave(file.path("plots", "posteriorh2.png"), p)

# Plot Evolvability
p <- plot_posterior(evolvability, "Posterior of Evolvability (e)", "e")
p
#ggsave(file.path("plots", "posteriorE.png"), p)

# Plot VA for cold and warm
p <- plot_posterior(VA_cold, "Posterior of VA in Cold Environment", "VA_cold")
p
#ggsave(file.path("plots", "posteriorVAc.png"), p)

p <- plot_posterior(VA_warm, "Posterior of VA in Warm Environment", "VA_warm")
p
#ggsave(file.path("plots", "posteriorVAw.png"), p)


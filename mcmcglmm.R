# Polyommatus Animal Model - GxE 

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


#############################
## SIMULATED DATA
#############################

# set.seed(123)
# 
# n_families <- 40
# offspring_per_family <- 10
# n_offspring <- n_families * offspring_per_family
# 
# # Factors
# MotherID <- factor(rep(1:n_families, each = offspring_per_family))
# Temperature <- factor(rep(rep(c("cold","warm"), each = offspring_per_family/2), times = n_families))
# Region <- factor(rep(sample(c("Skåne","Öland"), n_families, replace = TRUE), each = offspring_per_family))
# 
# # Simulate genetic and environmental effects
# VA_true <- 0.8
# Vm_true <- 0.3
# Ve_true <- 1.0
# 
# # Mother genetic effect
# mother_effect <- rnorm(n_families, 0, sqrt(VA_true))
# 
# # Maternal effect
# maternal_effect <- rnorm(n_families, 0, sqrt(Vm_true))
# 
# # Temperature effect
# temp_effect <- ifelse(Temperature == "warm", 0.5, -0.5)
# 
# # Region effect
# region_effect <- ifelse(Region == "Skåne", 0.3, -0.3)
# 
# # Individual residual
# residual <- rnorm(n_offspring, 0, sqrt(Ve_true))
# 
# # Phenotype
# Blueness <- 3 + mother_effect[MotherID] + maternal_effect[MotherID] + temp_effect + region_effect + residual
# 
# # Create data frame
# data <- data.frame(animal = paste0("ind", 1:n_offspring),
#                    MotherID = MotherID,
#                    Temperature = Temperature,
#                    Region = Region,
#                    Blueness = Blueness)
# 
# # Add mothers to pedigree
# mothers <- data.frame(animal = paste0("mom", 1:n_families),
#                       MotherID = NA,
#                       Temperature = NA,
#                       Region = NA,
#                       Blueness = NA)
# 
# data_full <- rbind(data, mothers)
# 
# ##############################################
# # 2. Create Pedigree
# ##############################################
# pedigree <- data.frame(
#   animal = c(mothers$animal, data$animal),
#   dam = c(rep(NA, n_families), data$MotherID),
#   sire = NA
# )
# 
# # Replace MotherID with mother names
# pedigree$dam <- ifelse(is.na(pedigree$dam), NA, paste0("mom", pedigree$dam))
# 
# ##############################################
# # 3. Define Priors
# ##############################################
# prior <- list(
#   G = list(
#     G1 = list(V = 1, nu = 0.002),  # Animal (VA)
#     G2 = list(V = 1, nu = 0.002),  # MotherID (maternal)
#     G3 = list(V = diag(2), nu = 0.002) # GxE for Temperature
#   ),
#   R = list(V = 1, nu = 0.002)      # Residual
# )
# 
# ##############################################
# # 4. Fit Animal Model
# ##############################################
# model <- MCMCglmm(
#   Blueness ~ Temperature * Region,
#   random = ~ animal + MotherID + idh(Temperature):animal,
#   pedigree = pedigree,
#   data = data,
#   family = "gaussian",
#   prior = prior,
#   nitt = 130000, burnin = 30000, thin = 100
# )
# 
# summary(model)
# 
# ##############################################
# # 5. Extract Genetic Parameters
# ##############################################
# VA <- model$VCV[,"animal"]
# Vres <- model$VCV[,"units"]
# Vm <- model$VCV[,"MotherID"]
# 
# # Heritability
# h2 <- VA / (VA + Vres + Vm)
# 
# # Evolvability
# meanBlueness <- mean(data$Blueness)
# evolvability <- VA / (meanBlueness^2)
# 
# cat("Posterior mean VA:", mean(VA), "\n")
# cat("Posterior mean h2:", mean(h2), "\n")
# cat("Posterior mean evolvability:", mean(evolvability), "\n")
# 
# ##############################################
# # 6. GxE Interaction
# ##############################################
# # idh(Temperature):animal gives separate VA for cold and warm
# VA_cold <- model$VCV[,"Temperaturecold.animal"]
# VA_warm <- model$VCV[,"Temperaturewarm.animal"]
# 
# cat("VA cold:", mean(VA_cold), "\n")
# cat("VA warm:", mean(VA_warm), "\n")
# 
# 
# #EXTENDED PLOT FOR VISUALIZATIONS
# 
# ##############################################
# # Extended Script: Posterior Plots
# ##############################################
# 
# ##############################################
# # 5. Extract Genetic Parameters (same as before)
# ##############################################
# VA <- model$VCV[,"animal"]
# Vres <- model$VCV[,"units"]
# Vm <- model$VCV[,"MotherID"]
# 
# h2 <- VA / (VA + Vres + Vm)
# meanBlueness <- mean(data$Blueness)
# evolvability <- VA / (meanBlueness^2)
# 
# VA_cold <- model$VCV[,"Temperaturecold.animal"]
# VA_warm <- model$VCV[,"Temperaturewarm.animal"]
# 
# ##############################################
# # 7. Posterior Density Plots
# ##############################################
# 
# # Function to plot posterior distributions
# plot_posterior <- function(samples, title, xlab){
#   df <- data.frame(value = samples)
#   ggplot(df, aes(x = samples)) +
#     geom_density(fill = "steelblue", alpha = 0.6) +
#     theme_minimal() +
#     labs(title = title, x = xlab, y = "Density") +
#     geom_vline(xintercept = mean(samples), color = "red", linetype = "dashed") +
#     annotate("text", x = mean(samples), y = 0, label = paste0("Mean = ", round(mean(samples), 3)), hjust = -0.1, color = "red")
# }
# 
# # Plot VA
# plot_posterior(VA, "Posterior of Additive Genetic Variance (VA)", "VA")
# 
# # Plot h²
# plot_posterior(h2, "Posterior of Heritability (h²)", "h²")
# 
# # Plot Evolvability
# plot_posterior(evolvability, "Posterior of Evolvability (e)", "e")
# 
# # Plot VA for cold and warm
# plot_posterior(VA_cold, "Posterior of VA in Cold Environment", "VA_cold")
# plot_posterior(VA_warm, "Posterior of VA in Warm Environment", "VA_warm")


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

data_full <- rbind(data, mothers)

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
    G3 = list(V = diag(2), nu = 0.002) # GxE for Temperature
  ),
  R = list(V = 1, nu = 0.002)      # Residual
)

##############################################
# 4. Fit Animal Model
##############################################
model <- MCMCglmm(
  Blueness ~ TotalArea + Temperature + Region + TotalArea:Temperature + Temperature:Region,
  random = ~ animal + MotherID + idh(Temperature):animal,
  pedigree = pedigree,
  data = data,
  family = "gaussian",
  prior = prior,
  nitt = 130000, burnin = 30000, thin = 100
)

summary(model)

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

##############################################
# 6. GxE Interaction
##############################################
# idh(Temperature):animal gives separate VA for cold and warm
VA_cold <- model$VCV[,"TemperatureCold (18°C).animal"]
VA_warm <- model$VCV[,"TemperatureWarm (26°C).animal"]

cat("VA cold:", mean(VA_cold), "\n")
cat("VA warm:", mean(VA_warm), "\n")


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

# CHECK FOR COVERGENCE
plot(model$Sol)
plot(model$VCV)


# Plot VA
p <- plot_posterior(VA, "Posterior of Additive Genetic Variance (VA)", "VA")
p
ggsave(file.path("plots", "posteriorVA.png"), p)

# Plot h²
p <- plot_posterior(h2, "Posterior of Heritability (h²)", "h²")
p
ggsave(file.path("plots", "posteriorh2.png"), p)

# Plot Evolvability
p <- plot_posterior(evolvability, "Posterior of Evolvability (e)", "e")
p
ggsave(file.path("plots", "posteriorE.png"), p)

# Plot VA for cold and warm
p <- plot_posterior(VA_cold, "Posterior of VA in Cold Environment", "VA_cold")
p
ggsave(file.path("plots", "posteriorVAc.png"), p)

p <- plot_posterior(VA_warm, "Posterior of VA in Warm Environment", "VA_warm")
p
ggsave(file.path("plots", "posteriorVAw.png"), p)


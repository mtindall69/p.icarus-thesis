# Polyommatus icarus project: Blue area analyses

#clear environment
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm,psych)

#load data
bluedata <- read.csv("blueness.csv")

glimpse(bluedata)
# convert motherID and population to categorical variables
bluedata$motherID = as.factor(bluedata$motherID)
bluedata$pop = as.factor(bluedata$pop)

#=======================================================================
# sanity check, is blue score a good predictor of blue area?
m <- lm(blue_mm ~ total_mm + daughterscore, data=bluedata)
summary(m) #YES

#mean scale blue area
plot(bluedata$daughterscore, bluedata$blue_mm)

#=======================================================================
# MODELS
#=======================================================================

# Is blueness explained by temp (plasticity), population and family (genetics)?
m1 <- lmer(blue_mm ~ total_mm + temp + pop + (1|motherID), data=bluedata)
summary(m1)
Anova(m1)
AIC(m1) #8523.431

# Adding interactions
# Is the interaction between temp and pop significant? 
m2 <- lmer(blue_mm ~ total_mm + temp * pop + (1|motherID), data=bluedata)
summary(m2)
Anova(m2) #YES VERY SIGNIFICANT
AIC(m2) #8408.129

m3 <- lmer(blue_mm ~ total_mm * temp * pop + (1|motherID), data=bluedata)
summary(m3)
Anova(m3) #temp*pop YES, area*temp YES, area*pop YES, but area*temp*pop NO
AIC(m3) #8321.624

m4 <- lmer(blue_mm ~ (total_mm + temp + pop)^2 + (1 | motherID), data = bluedata)
summary(m4)
Anova(m4)
AIC(m4) #8311.56 best so far, all interactions highly significant

#subtract area*pop interaction
m5 <- lmer(blue_mm ~ total_mm + temp + pop + (1 | motherID), data = bluedata)
summary(m5)
Anova(m5) 
AIC(m5)

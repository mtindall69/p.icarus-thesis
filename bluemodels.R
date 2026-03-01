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

# Adding interactions
# Is the interaction between temp and pop significant? 
m2 <- lmer(blue_mm ~ total_mm + temp * pop + (1|motherID), data=bluedata)
summary(m2)
Anova(m1) #YES VERY SIGNIFICANT



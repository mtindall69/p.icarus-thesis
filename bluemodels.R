# Polyommatus icarus project: Blue area analyses

#clear environment
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm,psych,effects)

#==================================================================
# LOAD AND PREPROCESS DATA
#==================================================================

#load data
bluedata <- read.csv("blueness.csv")

str(bluedata)
glimpse(bluedata)
# convert categorical variables
bluedata$motherID = as.factor(bluedata$motherID)
bluedata$region = as.factor(bluedata$region)
bluedata$temp = as.factor(bluedata$temp)
bluedata$sex = as.factor(bluedata$sex)
bluedata$temp_label <- ifelse(bluedata$temp == 18, "Cold (18°C)", "Warm (26°C)")
bluedata$region_label  <- ifelse(bluedata$region == "O", "Öland", "Skåne")
bluedata$sex_label <- ifelse(bluedata$sex == "F", "Female", "Male")


#create aggregate dataframe of averages per individual
# avg wing area, avg blue area, avgs for forewing and hindwing

aggregate_to_individual <- function(bluedata) {
  group_cols <- c("offspringID", "motherID", "region", "region_label",
                  "temp", "temp_label", "sex", "sex_label", "motherscore", 
                  "daughterscore", "pupation_weight_g", "adult_weight_g", 
                  "start_day", "pupa_day", "adult_day", "pupation_length")
  measure_cols <- c("total_mm", "blue_mm", "prop_blue")
  
  # Overall averages (all 4 wings)
  agg_all <- bluedata %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(across(all_of(measure_cols), mean), .groups = "drop")
  
  # Forewing averages (FL and FR only)
  agg_fw <- bluedata %>%
    filter(wing %in% c("FL", "FR")) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      avg_fw_total_mm  = mean(total_mm,  na.rm = TRUE),
      avg_fw_blue_mm   = mean(blue_mm,   na.rm = TRUE),
      avg_fw_prop_blue = mean(prop_blue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Hindwing averages (HL and HR only)
  agg_hw <- bluedata %>%
    filter(wing %in% c("HL", "HR")) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      avg_hw_total_mm  = mean(total_mm,  na.rm = TRUE),
      avg_hw_blue_mm   = mean(blue_mm,   na.rm = TRUE),
      avg_hw_prop_blue = mean(prop_blue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join all three together
  agg <- agg_all %>%
    left_join(agg_fw, by = group_cols) %>%
    left_join(agg_hw, by = group_cols) %>%
    rename(
      avg_total_mm  = total_mm,
      avg_blue_mm   = blue_mm,
      avg_prop_blue = prop_blue
    ) %>%
    as.data.frame()
  
  agg
}

bluesum <- aggregate_to_individual(bluedata)
write.csv(bluesum, "bluesum.csv", row.names = FALSE)


#subset data
blueFdata <- subset(bluesum, sex=="F")
blueMdata <- subset(bluesum, sex=="M")

#=======================================================================
# sanity check, is blue score a good predictor of blue area?
m <- lm(blue_mm ~ total_mm + daughterscore, data=blueFdata)
summary(m) #YES

# adjust blue area for total area, not just prop
plot(blueFdata$daughterscore, blueFdata$prop_blue)

#=======================================================================
# MODELS
#=======================================================================

# Is blueness explained by temp (plasticity), population and family (genetics)?
m1 <- lmer(blue_mm ~ total_mm + temp + region + (1|motherID), data=blueFdata)
summary(m1)
Anova(m1)
AIC(m1) #8523.431



# Adding interactions
# Is the interaction between temp and pop significant? 
m2 <- lmer(blue_mm ~ total_mm + temp * region + (1|motherID), data=blueFdata)
summary(m2)
Anova(m2) #YES VERY SIGNIFICANT
AIC(m2) #8408.129
plot(allEffects(m2))
ef<-effect("temp:region", m2)
summary(ef)

# blue males
m2b <- lmer(blue_mm ~ total_mm + temp + region + (1|motherID), data=blueMdata)
summary(m2b)
Anova(m2b) # temp and pop neither significant for male blueness
plot(allEffects(m2b))
ef<-effect("temp:region", m2b)

#area
ma <- lmer(total_mm ~ temp * region + (1|motherID)+ (1|start_day), data=blueFdata)
summary(m2b)
Anova(ma) # temp significant to size
plot(allEffects(ma))
ef<-effect("temp:region", ma)
summary(ef)

# area males
m2a <- lmer(total_mm ~ temp * region + (1|motherID)+ (1|start_day), data=blueMdata)
summary(m2a)
Anova(m2a) # temp significant to size
plot(allEffects(m2a))
ef<-effect("temp:region", m2a)
summary(ef)

#mother as fixed effect
mb <- lm(blue_mm ~ total_mm + temp * motherID, data=blueFdata)
summary(mb)
Anova(mb) 
AIC(mb) #8224.183
#plot(allEffects(mb))
ef<-effect("temp:motherID", mb)
summary(ef)
        
#both sexes
m3 <- lmer(blue_mm ~ total_mm + sex + temp + region + (1|motherID), data=bluedata)
summary(m3)
Anova(m3)
AIC(m3)

m4 <- lmer(blue_mm ~ total_mm + sex + temp * region + (1|motherID), data=bluedata)
summary(m4)
Anova(m4)
AIC(m4)

m5 <- lmer(blue_mm ~ total_mm + sex * temp * region + (1|motherID), data=bluedata)
summary(m5)
Anova(m5) #all interactions significant for blue area
AIC(m5)

m6 <- lmer(total_mm ~ sex + temp * region + (1|motherID), data=bluedata)
summary(m6)
Anova(m6)
AIC(m6)

m7 <- lmer(total_mm ~ sex * temp * region + (1|motherID), data=bluedata)
summary(m7)
Anova(m7) #temp*region significant for total area but other interactions NOT
AIC(m7)

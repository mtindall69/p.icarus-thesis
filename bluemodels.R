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

#load aggregated data
bluesum <- read.csv("bluesum.csv")
glimpse(bluesum)

bluesum$motherID = as.factor(bluesum$motherID)
bluesum$region = as.factor(bluesum$region)
bluesum$temp = as.factor(bluesum$temp)
bluesum$sex = as.factor(bluesum$sex)
bluesum$temp_label <- ifelse(bluesum$temp == 18, "Cold (18°C)", "Warm (26°C)")
bluesum$region_label  <- ifelse(bluesum$region == "O", "Öland", "Skåne")
bluesum$sex_label <- ifelse(bluesum$sex == "F", "Female", "Male")

#subset data
blueFdata <- subset(bluesum, sex=="F")
blueMdata <- subset(bluesum, sex=="M")

#=======================================================================
# sanity check, is blue score a good predictor of blue area?
n <- lm(avg_blue_mm ~ avg_total_mm + daughterscore, data=blueFdata)
summary(n) #YES

# adjust blue area for total area, not just prop
plot(blueFdata$daughterscore, blueFdata$avg_prop_blue)


#=======================================================================
# MODELS
#=======================================================================

# Is blueness explained by temp (plasticity), region and family (genetics)?
m <- lmer(avg_blue_mm ~ avg_total_mm + temp + region + (1|motherID), data=blueFdata)
summary(m)
Anova(m)
AIC(m) #2185.787


# Adding interactions
# Is the interaction between temp and region significant? 
m1 <- lmer(avg_blue_mm ~ avg_total_mm + temp * region + (1|motherID), data=blueFdata)
summary(m1)
Anova(m1) #YES temp:region VERY SIGNIFICANT
AIC(m1) #2157.852
plot(allEffects(m1))
ef<-effect("temp:region", m1)
summary(ef)


#mother as fixed effect
m2 <- lm(avg_blue_mm ~ avg_total_mm + temp * motherID, data=blueFdata)
summary(m2)
Anova(m2) # motherID and temp:motherID significant
AIC(m2) #2153.591
ef<-effect("temp:motherID", m2)
summary(ef) # output for reaction norms


mb <- lm(avg_blue_mm ~ avg_total_mm + temp * region + motherID, data=blueFdata)
summary(mb)
Anova(mb) 
AIC(mb) #2125.188

mb <- lmer(avg_blue_mm ~ avg_total_mm + temp * region + motherID + (1|motherID), data=blueFdata)
summary(mb)
Anova(mb) # motherID as fixed effect not significant
AIC(mb) #1987.349

mb <- lmer(avg_blue_mm ~ avg_total_mm + temp * region * motherID + (1|motherID), data=blueFdata)
summary(mb)
Anova(mb) # motherID nor interactions significant
AIC(mb) #1845.377

mb <- lm(avg_blue_mm ~ avg_total_mm + temp + region + temp:region + motherID 
         + temp:motherID, data=blueFdata)
summary(mb)
Anova(mb) # region and temp:region not significant
AIC(mb) #2153.591


mb <- lmer(avg_blue_mm ~ avg_total_mm + temp + region + motherID 
         + temp:region + temp:motherID + (1|motherID), data=blueFdata)
summary(mb)
Anova(mb) 
AIC(mb) #1845.377



# MALE BLUENESS

m2b <- lmer(avg_blue_mm ~ avg_total_mm + temp * region + (1|motherID), data=blueMdata)
summary(m2b)
Anova(m2b) # temp significant, region not, interaction not
AIC(m2b) #638.6096
plot(allEffects(m2b))
ef<-effect("temp:region", m2b)
summary(ef)

m2b <- lmer(avg_blue_mm ~ avg_total_mm + temp + region + (1|motherID), data=blueMdata)
summary(m2b)
Anova(m2b) # temp significant, region not
AIC(m2b) #638.1577
plot(allEffects(m2b))

m2b <- lmer(avg_blue_mm ~ avg_total_mm + temp + (1|motherID), data=blueMdata)
summary(m2b)
Anova(m2b) # temp significant, region not
AIC(m2b) #637.2989
plot(allEffects(m2b))




# AREA

ma <- lmer(avg_total_mm ~ temp * region + (1|motherID)+ (1|start_day), data=blueFdata)
summary(m2b)
Anova(ma) # temp significant to size
plot(allEffects(ma))
ef<-effect("temp:region", ma)
summary(ef)

# area males
m2a <- lmer(avg_total_mm ~ temp * region + (1|motherID)+ (1|start_day), data=blueMdata)
summary(m2a)
Anova(m2a) # temp significant to size
plot(allEffects(m2a))
ef<-effect("temp:region", m2a)
summary(ef)


        

#BOTH SEXES

#Blueness
m3 <- lmer(avg_blue_mm ~ avg_total_mm + sex + temp + region + (1|motherID), data=bluesum)
summary(m3)
Anova(m3)
AIC(m3) #3244.193

m4 <- lmer(avg_blue_mm ~ avg_total_mm + sex + temp * region + (1|motherID), data=bluesum)
summary(m4)
Anova(m4) # temp:region significant
AIC(m4) #3229.128

m5 <- lmer(avg_blue_mm ~ avg_total_mm + sex * temp * region + (1|motherID), data=bluesum)
summary(m5)
Anova(m5) # sex:temp and sex:temp:region not significant
AIC(m5) #3138.222

m5 <- lmer(avg_blue_mm ~ avg_total_mm * sex * temp * region + (1|motherID), data=bluesum)
summary(m5)
Anova(m5) # avg_total_mm:temp:region and avg_total_mm:sex:temp:region not significant
AIC(m5) #2995.347

m5 <- lmer(avg_blue_mm ~ avg_total_mm + sex + temp + region + 
             avg_total_mm:sex + avg_total_mm:temp + avg_total_mm:region +
             sex:temp + sex:region + temp:region +
             sex:temp:region + avg_total_mm:sex:region +
             (1|motherID), data=bluesum)
summary(m5)
Anova(m5) #all interactions significant for blue area
AIC(m5) #2985.876

m5 <- lmer(avg_blue_mm ~ avg_total_mm + sex + temp + region + 
             avg_total_mm:sex + avg_total_mm:temp + avg_total_mm:region +
             sex:temp + sex:region + temp:region +
             sex:temp:region + (1|motherID), data=bluesum)
summary(m5)
Anova(m5) #all interactions significant for blue area
AIC(m5) #2986.796

m5 <- lmer(avg_blue_mm ~ avg_total_mm + sex + temp + region + 
             avg_total_mm:sex + avg_total_mm:temp + avg_total_mm:region +
             sex:temp + sex:region + temp:region + (1|motherID), data=bluesum)
summary(m5)
Anova(m5) #all interactions significant for blue area
AIC(m5) #2997.341


#Area
m6 <- lmer(avg_total_mm ~ sex * temp * region + (1|motherID), data=bluesum)
summary(m6)
Anova(m6)#temp:region significant for total area but other interactions NOT
AIC(m6) #3679.995
plot(allEffects(m6))

m7 <- lmer(avg_total_mm ~ sex + temp * region + (1|motherID), data=bluesum)
summary(m7)
Anova(m7) 
AIC(m7) #3686.801
plot(allEffects(m7))



#================================================================
# VISUALIZATION
#================================================================

plot(blueFdata$motherscore, blueFdata$avg_prop_blue)

plot(blueFdata$avg_total_mm, blueFdata$avg_blue_mm)
#total 
t <- lm(avg_blue_mm ~ avg_total_mm, data=blueFdata)
abline(t)
#temps
WFdata <- subset(blueFdata, temp =="26") 
CFdata <- subset(blueFdata, temp =="18")
w <- lm(avg_blue_mm ~ avg_total_mm, data=WFdata)
c <- lm(avg_blue_mm ~ avg_total_mm, data=CFdata)
plot(blueFdata$avg_total_mm, blueFdata$avg_blue_mm, col=ifelse(blueFdata$temp == 26, "red", "blue"))
abline(w, col="red")
abline(c, col="blue")

#regions
OFdata <- subset(blueFdata, region =="O")
SFdata <- subset(blueFdata, region =="S")
o <- lm(avg_blue_mm ~ avg_total_mm, data=OFdata)
s <- lm(avg_blue_mm ~ avg_total_mm, data=SFdata)
plot(blueFdata$avg_total_mm, blueFdata$avg_blue_mm, col=ifelse(blueFdata$region == "O", "purple", "green"))
abline(o, col="purple")
abline(s, col="green")

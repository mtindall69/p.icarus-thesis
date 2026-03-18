# Polyommatus Weight Loss and Metabolic Cost of Color

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm)

#load data
bluefemales = read_excel("Blue pupae.xlsx")
blueoffspring = read_excel("Blue offspring.xlsx")

#subset data into regions and temperatures, SEXES
#FEMALES
ofemales = subset(bluefemales, region=='O')
sfemales = subset(bluefemales, region=='S')
coldfemales = subset(bluefemales, Temp=='18')
warmfemales = subset(bluefemales, Temp=='26')
ocold = subset(coldfemales, region=='O')
owarm = subset(warmfemales, region=='O')
scold = subset(coldfemales, region=='S')
swarm = subset(warmfemales, region=='S')

#ALL
ooff = subset(blueoffspring, region=='O')
soff = subset(blueoffspring, region=='S')
coldoff = subset(blueoffspring, Temp=='18')
warmoff = subset(blueoffspring, Temp=='26')
allocold = subset(coldoff, region=='O')
allowarm = subset(warmoff, region=='O')
allscold = subset(coldoff, region=='S')
allswarm = subset(warmoff, region=='S')

#MALES
bluemales <- subset(blueoffspring, Sex=='M')
omales = subset(bluemales, region=='O')
smales = subset(bluemales, region=='S')
coldmales = subset(bluemales, Temp=='18')
warmmales = subset(bluemales, Temp=='26')
ocoldm = subset(coldmales, region=='O')
owarmm = subset(warmmales, region=='O')
scoldm = subset(coldmales, region=='S')
swarmm = subset(warmmales, region=='S')

#================================================================
#################################
# MEAN SIZES
#################################

#ALL OFFSPRING MEAN SIZES
mean(blueoffspring$Adult_weight_g) # avg for all offspring is 0.02290g
mean(ooff$Adult_weight_g) # 0.02243g
mean(soff$Adult_weight_g) # 0.02357g
#Skane slightly larger
mean(coldoff$Adult_weight_g) # 0.02595g
mean(warmoff$Adult_weight_g) # 0.01875g
#Size decreases with temp, 26 much smaller than 18
mean(allocold$Adult_weight_g) # 0.02566g
mean(allscold$Adult_weight_g) # 0.02640g
#Skane slightly larger in 18
mean(allowarm$Adult_weight_g) # 0.01754g
mean(allswarm$Adult_weight_g) # 0.02022g
#Skane larger in 26


#FEMALE Mean sizes
mean(bluefemales$Adult_weight_g) # avg for all females is 0.02315g
mean(ofemales$Adult_weight_g) # 0.02316g
mean(sfemales$Adult_weight_g) # 0.02314g
#regions do not differ much in weight, Oland marginally larger
mean(coldfemales$Adult_weight_g) # 0.02569g
mean(warmfemales$Adult_weight_g) # 0.01912g
#Size decreases with temp, 26 much smaller than 18
mean(ocold$Adult_weight_g) # 0.02579g
mean(scold$Adult_weight_g) # 0.02554g 
#Oland are slightly larger in 18
mean(owarm$Adult_weight_g) # 0.01824g
mean(swarm$Adult_weight_g) # 0.01995g
#Skane larger in 26, Oland size more sensitive to warm
# significantly different? chi square? sig figs and error?

#MALE Mean sizes
mean(bluemales$Adult_weight_g) # avg for all males is 0.02259g 
        # Average male is smaller/lighter than female
mean(omales$Adult_weight_g) # 0.02168g
mean(smales$Adult_weight_g) # 0.02424g
# Skane males larger than Oland
mean(coldmales$Adult_weight_g) # 0.02632g
mean(warmmales$Adult_weight_g) # 0.01837g
#Size decreases with temp, 26 much smaller than 18
    # Males larger than females in cold but smaller in warm
mean(ocoldm$Adult_weight_g) # 0.02551g
mean(scoldm$Adult_weight_g) # 0.02798g
#Skane larger in 18
mean(owarmm$Adult_weight_g) # 0.01699g
mean(swarmm$Adult_weight_g) # 0.02059g
#Skane much larger in 26

#Size decreases with temp BUT
#Oland females are larger than Skane in cold, smaller in warm
#Oland males are smaller than Skane in cold and much smaller in warm
      # Oland size more sensitive to warm
#Males generally smaller, slightly larger in cold, slightly smaller in warm BUT
#Skane males are bigger than Skane females, larger in cold, about same in warm
#Oland males are smaller than Oland females, about same in cold, smaller in warm

#===============================================================

# WEIGHT LOSS AND METABOLIC COST OF COLOR

par(mfrow=c(2,1))

#Females
plot(bluefemales$Pupation_weight_g, bluefemales$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Female Pupa Weight to Adult Weight")
m <- lm(bluefemales$Adult_weight_g ~ bluefemales$Pupation_weight_g)
abline(m, col="red")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.458g adult/g pupa, adult about 46% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.46")),"p < 0.01"))

#Males
plot(bluemales$Pupation_weight_g, bluemales$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Male Pupa Weight to Adult Weight")
m <- lm(bluemales$Adult_weight_g ~ bluemales$Pupation_weight_g)
abline(m, col="blue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.429g adult/g pupa, adult about 43% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.43")),"p < 0.01"))
  #Males lose more weight, blue more expensive? or male weight not as important

par(mfrow=c(1,1))

#Offspring
colors = c("red", "blue")
plot(blueoffspring$Pupation_weight_g, blueoffspring$Adult_weight_g, 
     col=colors[factor(blueoffspring$Sex)], pch=20,
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Pupa Weight to Adult Weight")
m <- lm(blueoffspring$Adult_weight_g ~ blueoffspring$Pupation_weight_g)
abline(m)
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.439g adult/g pupa, adult about 44% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.44")),"p < 0.01"))
legend("topleft", legend=c("Female", "Male"), pch=20, 
       col=c("red","blue"))



#FOR REGIONS AND TEMPS

par(mfrow=c(2,1))

#Cold
plot(coldoff$Pupation_weight_g, coldoff$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Cold Pupa Weight to Adult Weight")
m <- lm(coldoff$Adult_weight_g ~ coldoff$Pupation_weight_g)
abline(m, col="blue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.418g adult/g pupa, adult about 42% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.42")),"p < 0.01"))

#Warm
plot(warmoff$Pupation_weight_g, warmoff$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Warm Pupa Weight to Adult Weight")
m <- lm(warmoff$Adult_weight_g ~ warmoff$Pupation_weight_g)
abline(m, col="red")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.405g adult/g pupa, adult about 40% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.40")),"p < 0.01"))

#Oland
plot(ooff$Pupation_weight_g, ooff$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Oland Pupa Weight to Adult Weight")
m <- lm(ooff$Adult_weight_g ~ ooff$Pupation_weight_g)
abline(m, col="steelblue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.452g adult/g pupa, adult about 46% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.45")),"p < 0.01"))

#Skane
plot(soff$Pupation_weight_g, soff$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Skane Pupa Weight to Adult Weight")
m <- lm(soff$Adult_weight_g ~ soff$Pupation_weight_g)
abline(m, col="orange")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.429g adult/g pupa, adult about 43% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.43")),"p < 0.01"))

par(mfrow=c(1,1))


#=========================================================
# Mixed Models - 
# IS THERE A METABOLIC COST TO PRODUCE DIFFERENT COLORS?
#=========================================================

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
#======================================================================

# Do blue or brown lose more weight?

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp + daughterscore + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #daughter score NOT significant
plot(allEffects(mw))
summary(effect("daughterscore", mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp + avg_prop_blue + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #prop blueness NOT significant
plot(allEffects(mw))
summary(effect("avg_prop_blue", mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * temp * daughterscore + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #daughter score NOT significant, pupation_weight:daughterscore IS
plot(allEffects(mw))
summary(effect("daughterscore", mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * temp * avg_prop_blue + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #prop blueness NOT significant nor any interactions
plot(allEffects(mw))
summary(effect("avg_prop_blue", mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp * daughterscore + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw)
plot(allEffects(mw))
summary(effect("daughterscore", mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp * avg_prop_blue + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw)
plot(allEffects(mw))
summary(effect("avg_prop_blue", mw))


mw2 <- lmer(adult_weight_g ~ pupation_weight_g + temp + daughterscore + 
              pupation_weight_g:temp + pupation_weight_g:daughterscore + 
              (1|motherID), data=blueFdata)
summary(mw2)
Anova(mw2) # daughterscore NOT significant, interactions are

mw2 <- lmer(adult_weight_g ~ pupation_weight_g + temp + avg_prop_blue + 
              pupation_weight_g:temp + pupation_weight_g:avg_prop_blue + 
              (1|motherID), data=blueFdata)
summary(mw2)
Anova(mw2) # avg_prop_blue NOT significant nor any interactions


# Do regions differ in weight loss? 
mw <- lmer(adult_weight_g ~ pupation_weight_g * temp * region + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #temp:region significant but region is not
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp * region + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #temp:region significant but region is not
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g + temp + region + (1|motherID), 
           data=blueFdata)
summary(mw)
Anova(mw) #region NOT significant
plot(allEffects(mw))


#Do sexes differ in weight loss?
mw <- lmer(adult_weight_g ~ pupation_weight_g * sex + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw) # sex significant, pupation_weight:sex weakly significant
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex + temp + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw) # YES
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * temp + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw) # no interactions significant
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * region + temp + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw) 
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * temp * region + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw)
plot(allEffects(mw))


mw <- lmer(adult_weight_g ~ pupation_weight_g * sex + temp * region + (1|motherID), 
           data=bluesum)
summary(mw)
Anova(mw) # All effects significant
plot(allEffects(mw)) # pupation weight:sex and temp:region interactions significant

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

#subset data into populations and temperatures, SEXES
#FEMALES
ofemales = subset(bluefemales, MotherPop=='O')
sfemales = subset(bluefemales, MotherPop=='S')
coldfemales = subset(bluefemales, Temp=='18')
warmfemales = subset(bluefemales, Temp=='26')
ocold = subset(coldfemales, MotherPop=='O')
owarm = subset(warmfemales, MotherPop=='O')
scold = subset(coldfemales, MotherPop=='S')
swarm = subset(warmfemales, MotherPop=='S')

#ALL
ooff = subset(blueoffspring, MotherPop=='O')
soff = subset(blueoffspring, MotherPop=='S')
coldoff = subset(blueoffspring, Temp=='18')
warmoff = subset(blueoffspring, Temp=='26')
allocold = subset(coldoff, MotherPop=='O')
allowarm = subset(warmoff, MotherPop=='O')
allscold = subset(coldoff, MotherPop=='S')
allswarm = subset(warmoff, MotherPop=='S')

#MALES
bluemales <- subset(blueoffspring, Sex=='M')
omales = subset(bluemales, MotherPop=='O')
smales = subset(bluemales, MotherPop=='S')
coldmales = subset(bluemales, Temp=='18')
warmmales = subset(bluemales, Temp=='26')
ocoldm = subset(coldmales, MotherPop=='O')
owarmm = subset(warmmales, MotherPop=='O')
scoldm = subset(coldmales, MotherPop=='S')
swarmm = subset(warmmales, MotherPop=='S')

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
#Populations do not differ much in weight, Oland marginally larger
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
#Skane males are bigger than Skane females, larger in coold, about same in warm
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
legend(0.08,0.02, legend=c(expression(paste(beta, " = 0.46")),"p < 0.01"), bty="n")

#Males
plot(bluemales$Pupation_weight_g, bluemales$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Male Pupa Weight to Adult Weight")
m <- lm(bluemales$Adult_weight_g ~ bluemales$Pupation_weight_g)
abline(m, col="blue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.429g adult/g pupa, adult about 43% pupal weight
legend(0.08,0.02, legend=c(expression(paste(beta, " = 0.43")),"p < 0.01"), bty="n")
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
legend(0.08,0.02, legend=c(expression(paste(beta, " = 0.44")),"p < 0.01"), bty="n")
legend("bottomright", legend=c("Female", "Male"), pch=20, 
       col=c("red","blue"))


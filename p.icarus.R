#Polyommatus icarus project: Female polymorphism, thermal plasticity, and GxE

#clear environments
rm(list=ls()) 

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm)
               #picante,geiger,phytools,ggtree,treeio,openxlsx,janitor,phylolm)

#load data
bluefemales = read_excel("Blue pupae.xlsx")
glimpse(bluefemales)
#how to treat larval jar? not continuous, factor? ASK
bluefemales$Larval_jar = as.factor(bluefemales$Larval_jar)

#subset data into populations and temperatures
ofemales = subset(bluefemales, MotherPop=='O')
sfemales = subset(bluefemales, MotherPop=='S')
coldfemales = subset(bluefemales, Temp=='18')
warmfemales = subset(bluefemales, Temp=='26')

ocold = subset(coldfemales, MotherPop=='O')
owarm = subset(warmfemales, MotherPop=='O')
scold = subset(coldfemales, MotherPop=='S')
swarm = subset(warmfemales, MotherPop=='S')

#============================================================

#Distribution of Blue scores
bluecount <- table(bluefemales$Blue_score)
print(bluecount)  # n=337 1: 106, 2: 40, 3: 53, 4: 54, 5: 84

ocount <- table(ofemales$Blue_score) # n=181
scount <- table(sfemales$Blue_score) # n=156
coldcount <- table(coldfemales$Blue_score) # n=207
warmcount <- table(warmfemales$Blue_score) # n=130

ocoldcount <- table(ocold$Blue_score) # 1: 3, 2: 5, 3: 20, 4: 28, 5: 62 (118)
owarmcount <- table(owarm$Blue_score) # 1: 19, 2: 9, 3: 10, 4: 15, 5: 10 (63)
scoldcount <- table(scold$Blue_score) # 1: 34, 2: 19, 3: 18, 4: 7, 5: 11 (89)
swarmcount <- table(swarm$Blue_score) # 1: 50, 2: 7, 3: 5, 4: 4, 5: 1 (67)

#Make table of temp counts by pop
#coldpops <- print(ocold$Blue_score, scold$Blue_score)
#warmpops <- print(owarm$Blue_score, swarm$Blue_score)

#================================================================

#Mean sizes
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
#Skane larger in warm, Oland size more sensitive to warm
# significantly different? chi square? sig figs and error?

plot(bluefemales$Pupation_weight_g, bluefemales$Adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="Adult Weight (g)",
     main="Pupa Weight to Adult Weight")
m <- lm(bluefemales$Adult_weight_g ~ bluefemales$Pupation_weight_g)
abline(m)
#add y=x line
summary(m)$coef # 0.458g adult/g pupa, adult about 46% pupal weight
#INCLUDE MALES

# IS THERE A METABOLIC COST TO PRODUCE DIFFERENT COLORS?
mw <- lmer(Adult_weight_g ~ Pupation_weight_g * Temp * Blue_score + (1|MotherID), 
           data=bluefemales)
summary(mw)
Anova(mw)
plot(allEffects(mw)) # scores do not seem to lose different weights BUT
      # Blue score is significant p=0.0115, interaction is not?
#summary(effect(Blue_score, mw)) #doesn't work?

mw2 <- lmer(Adult_weight_g ~ Pupation_weight_g + Temp + Blue_score + 
              Pupation_weight_g:Temp + Pupation_weight_g:Blue_score + 
              (1|MotherID), data=bluefemales)
summary(mw2)
Anova(mw2)

#================================================================

## Preliminary data visualization
hist(bluefemales$Blue_score, breaks=6, xlab="Blue Score", 
     main="Distribution of Blue Scores") # 1 and 5 modes

hist(ofemales$Blue_score, xlab="Blue Score",
     main="Blue Scores of Oland Females") # Oland females are bluer

hist(sfemales$Blue_score, xlab="Blue Score",
     main="Blue Scores of Skane Females") # Skane females are browner

hist(coldfemales$Blue_score, xlab="Blue Score",
     main="Blue Scores of Cold Treatment Females") # Females are bluer in cold

hist(warmfemales$Blue_score, xlab="Blue Score",
     main="Blue Scores of Warm Treatment Females") # Females are browner in warm


#Distribution colored by population- compare total distribution and temp treatments
#cols=c("steelblue","orangered") FIX COLORS
ggplot(bluefemales, aes(x=Blue_score, color=MotherPop)) +
  geom_histogram(binwidth=1) +
  theme_bw()

ggplot(coldfemales, aes(x=Blue_score, color=MotherPop)) +
  geom_histogram(binwidth=1) +
  ggtitle("Blue Scores in Cold Treatment by Population") +
  theme_bw()

ggplot(warmfemales, aes(x=Blue_score, color=MotherPop)) +
  geom_histogram(binwidth=1) +
  ggtitle("Blue Scores in Warm Treatment by Population") +
  theme_bw()

#Side by side histograms of blue scores by pop
par(mfrow=c(2,2))
hist(ocold$Blue_score, xlab="Blue Score",
     main="Oland Cold Treatment Females")
hist(owarm$Blue_score, xlab="Blue Score",
     main="Oland Warm Treatment Females")
hist(scold$Blue_score, xlab="Blue Score",
     main="Skane Cold Treatment Females")
hist(swarm$Blue_score, xlab="Blue Score",
     main="Skane Warm Treatment Females")
par(mfrow=c(1,1))

#temperature effect on blue score by pop
# ggplot(data=bluefemales, aes(x=Temp, y=Blue_score, color=MotherPop)) +
#   geom_point() +
#   geom_line(aes(x=Temp, y=Blue_score, group=MotherID)) +
#   labs(x="Temperature", y="Blue Score") +
#   ggtitle("Temperature effect on Blueness by Population") +
#   theme_bw()
#doesn't really work


#===================================================================
## MOTHER FITNESS * BLUENESS

bluemothers = read_excel("Blue butterflies.xlsx")

# Fecundity vs blueness
plot(bluemothers$Blue_score, bluemothers$Total_Eggs,
     xlab="Mother Blue Score", ylab="Fecundity (eggs laid)",
     main="Blueness vs Fecundity")
o <- lm(Total_Eggs ~ Blue_score, data=bluemothers)
summary(o)
abline(o)
# 46 +/- 18 eggs per blue score, p=0.0132

# Survival vs blueness
plot(bluemothers$Blue_score, bluemothers$Offspring_matured,
     xlab="Mother Blue Score", ylab="Mature Offspring",
     main="Mother Blueness vs Offspring Survival")
v <- lm(Offspring_matured ~ Blue_score, data=bluemothers)
summary(v)
abline(v)
# 2.18 more offspring per blue score, p=0.00499

#===================================================================

# Linear Mixed Models
null <- lm(Blue_score ~ 1, data=bluefemales)
summary(null)

mod <- lmer(Blue_score ~ Temp + MotherPop + (1|MotherID), data=bluefemales)
Anova(mod)
summary(mod)

# POPULATION INTERACTION EFFECT
mod1 <- lmer(Blue_score ~ Temp * MotherPop + (1|MotherID), data=bluefemales)
Anova(mod1)
summary(mod1) # interaction of temp and motherpop is significant
plot(allEffects(mod1))


## MOTHER-OFFSPRING REGRESSION
ggplot(data=bluefemales, aes(x=MotherScore, y=Blue_score, color=MotherPop)) +
    geom_point() +
    geom_jitter(width=0.3, height=0.3) +
    labs(x="Mother Score", y="Daughter Score") +
    ggtitle("Mother-Offspring Blueness") +
    theme_bw()

#mean offspring vs mother score- slope*2= heritability
meanscores <- bluefemales %>%
  group_by(MotherID) %>%
  summarise(Blue_score = mean(Blue_score), MotherScore = mean(MotherScore))
print(n=42,meanscores)

plot(meanscores$MotherScore, meanscores$Blue_score, xlim=c(1,5), ylim=c(1,5),
     xlab="Mother Score", ylab="Mean Offspring Score",
     main="Mother Blue Score to Mean Offspring Blue Score")
b <- lm(Blue_score ~ MotherScore, data=meanscores)
abline(b)
#abline(a=0,b=1,lty=2) #extra line 
legend(0.75,3.1, legend=c(expression(paste(beta, " = 0.32")),
                        "p = 0.088",expression(paste("R"^2," = 0.0479"))), bty="n")
summary(b)

# Heritability h2
h2 <- summary(b)$coef[2]*2 #slope is 0.32 so 
h2 # heritability is 0.65
# Phenotypic variance Vp
vp <- var(bluefemales$Blue_score)
vp # 2.5337
# Additive genetic variance Va
va <- h2*vp 
va # 1.6437
# Evolvability (mean scale) e
e <- va/(mean(bluefemales$Blue_score))
e

meanscoresbytemp <- bluefemales %>%
  group_by(MotherID,Temp) %>%
  summarise(Blue_score = mean(Blue_score), MotherScore = mean(MotherScore))
print(n=42,meanscoresbytemp)
meanscoresbytemp$bluenesschange <- meanscoresbytemp$Blue_score - meanscoresbytemp$MotherScore
coldmeans <- subset(meanscoresbytemp, Temp==18)
warmmeans <- subset(meanscoresbytemp, Temp==26)

## HERITABILITY OF BLUENESS BY TEMP
par(mfrow=c(1,2))
plot(coldmeans$MotherScore, coldmeans$Blue_score, 
     xlim=c(1,5), ylim=c(1,5),
     xlab="Mother Score", ylab="Mean Offspring Score",
     main="Cold Mother-Offspring Regression")
k <- lm(Blue_score ~ MotherScore, data=coldmeans)
abline(k)
summary(k)
legend(0.75,4.1, legend=c(expression(paste(beta, " = 0.40")),
                          "p = 0.043",expression(paste("R"^2," = 0.0775"))), bty="n")


plot(warmmeans$MotherScore, warmmeans$Blue_score, 
     xlim=c(1,5), ylim=c(1,5),
     xlab="Mother Score", ylab="Mean Offspring Score",
     main="Warm Mother-Offspring Regression")
v <- lm(Blue_score ~ MotherScore, data=warmmeans)
abline(v)
summary(v)
legend(0.75,3.1, legend=c(expression(paste(beta, " = 0.39")),
                          "p = 0.044",expression(paste("R"^2," = 0.0833"))), bty="n")

par(mfrow=c(1,1))

## CHANGE IN BLUENESS BY MOTHER
ggplot(data=meanscoresbytemp, aes(x=Temp, y=Blue_score, color=MotherID)) +
  geom_point() +
  geom_jitter(width=0.5) +
  geom_line(aes(x=Temp, y=Blue_score, group=MotherID)) +
  ggtitle("Change in Mean Blueness by Mother") +
  theme_bw()

# Change in blueness from mother
ggplot(data=meanscoresbytemp, aes(x=MotherID, y=bluenesschange, fill=Temp)) +
  geom_bar(position="dodge", stat="identity")
barplot(height=meanscoresbytemp$bluenesschange, beside=TRUE, col=c("blue","red"))

#mcmcglmm - gxe interaction
#wombat - g-matrix?

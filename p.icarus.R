#Polyommatus icarus project: data exploration

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
# region, motherID, temp, larval jar are categorical variables
bluefemales$region = as.factor(bluefemales$region)
bluefemales$motherID = as.factor(bluefemales$motherID)
bluefemales$temp = as.factor(bluefemales$temp)
bluefemales$larval_jar = as.factor(bluefemales$larval_jar)

#subset data into regions and temperatures
ofemales = subset(bluefemales, region=='O')
sfemales = subset(bluefemales, region=='S')
coldfemales = subset(bluefemales, temp=='18')
warmfemales = subset(bluefemales, temp=='26')

ocold = subset(coldfemales, region=='O')
owarm = subset(warmfemales, region=='O')
scold = subset(coldfemales, region=='S')
swarm = subset(warmfemales, region=='S')

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


#Distribution colored by region- compare total distribution and temp treatments
#cols=c("steelblue","orangered") FIX COLORS
ggplot(bluefemales, aes(x=Blue_score, color=region)) +
  geom_histogram(binwidth=1) +
  theme_bw()

ggplot(coldfemales, aes(x=Blue_score, color=region)) +
  geom_histogram(binwidth=1) +
  ggtitle("Blue Scores in Cold Treatment by region") +
  theme_bw()

ggplot(warmfemales, aes(x=Blue_score, color=region)) +
  geom_histogram(binwidth=1) +
  ggtitle("Blue Scores in Warm Treatment by region") +
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
# ggplot(data=bluefemales, aes(x=temp, y=Blue_score, color=region)) +
#   geom_point() +
#   geom_line(aes(x=temp, y=Blue_score, group=MotherID)) +
#   labs(x="temperature", y="Blue Score") +
#   ggtitle("temperature effect on Blueness by region") +
#   theme_bw()
#doesn't really work


#===========================================================
# FITNESS * BLUENESS
#===========================================================

bluemothers = read_excel("Blue butterflies.xlsx")
omothers = subset(bluemothers, region=='O')
smothers = subset(bluemothers, region=='S')

# FECUNDITY VS BLUENESS
plot(bluemothers$Blue_score, bluemothers$Total_Eggs,
     xlab="Mother Blue Score", ylab="Fecundity (eggs laid)",
     main="Blueness vs Fecundity")
o <- lm(Total_Eggs ~ Blue_score, data=bluemothers)
summary(o)
abline(o)
legend("topleft", legend=c(expression(paste(beta, " = 46" %+-% "18")), 
                           "p = 0.0132"))
# 46 +/- 18 eggs per blue score, p=0.0132

#by region
par(mfrow=c(1,2))

plot(omothers$Blue_score, omothers$Total_Eggs,
     xlab="Mother Blue Score", ylab="Fecundity (eggs laid)",
     main="Oland Blueness vs Fecundity")
o <- lm(Total_Eggs ~ Blue_score, data=omothers)
summary(o)
abline(o, lty=2)
legend("topleft", legend=c(expression(paste(beta, " = 10" %+-% "38")), 
                           "p = 0.796"))

plot(smothers$Blue_score, smothers$Total_Eggs,
     xlab="Mother Blue Score", ylab="Fecundity (eggs laid)",
     main="Skane Blueness vs Fecundity")
o <- lm(Total_Eggs ~ Blue_score, data=smothers)
summary(o)
abline(o, lty=2)
legend("topleft", legend=c(expression(paste(beta, " = 20" %+-% "17")), 
                           "p = 0.229"))

par(mfrow=c(1,1))

# SURVIVAL VS BLUENESS
plot(bluemothers$Blue_score, bluemothers$Offspring_matured,
     xlab="Mother Blue Score", ylab="Mature Offspring",
     main="Mother Blueness vs Offspring Survival")
v <- lm(Offspring_matured ~ Blue_score, data=bluemothers)
summary(v)
abline(v)
legend("topleft", legend=c(expression(paste(beta, " = 2.18" %+-% "0.75")), 
                           "p = 0.005"))
# 2.18 more offspring per blue score, p=0.00499

#by region
par(mfrow=c(1,2))

plot(omothers$Blue_score, omothers$Offspring_matured,
     xlab="Mother Blue Score", ylab="Mature Offspring",
     main="Oland Mother Blueness vs Survivorship")
v <- lm(Offspring_matured ~ Blue_score, data=omothers)
summary(v)
abline(v)
legend("topleft", legend=c(expression(paste(beta, " = 2.4" %+-% "1.3")), 
                           "p = 0.068"))

plot(bluemothers$Blue_score, bluemothers$Offspring_matured,
     xlab="Mother Blue Score", ylab="Mature Offspring",
     main="Skane Mother Blueness vs Survivorship")
v <- lm(Offspring_matured ~ Blue_score, data=smothers)
summary(v)
abline(v)
legend("topleft", legend=c(expression(paste(beta, " = 1.8" %+-% "1.2")), 
                           "p = 0.14"))

par(mfrow=c(1,1))

# BLUE SCORE VS ADULT WEIGHT
plot(bluefemales$Blue_score, bluefemales$Adult_weight_g,
     xlab="Blue Score", ylab="Adult Weight (g)",
     main="Blueness to Adult Weight")
m <- lm(bluefemales$Adult_weight_g ~ bluefemales$Blue_score)
abline(m)
summary(m) # Weakly positive but significant- 0.001g/bluescore

#===================================================================
# LINEAR MIXED MODELS
#===================================================================
null <- lm(Blue_score ~ 1, data=bluefemales)
summary(null)

mod <- lmer(Blue_score ~ temp + region + (1|MotherID), data=bluefemales)
Anova(mod)
summary(mod)

# region INTERACTION EFFECT
mod1 <- lmer(Blue_score ~ temp * region + (1|MotherID), data=bluefemales)
Anova(mod1)
summary(mod1) # interaction of temp and region is significant
plot(allEffects(mod1))

#=====================================================================
## MOTHER-OFFSPRING REGRESSION
#=====================================================================

# All daughters visualized, colored by pop
ggplot(data=bluefemales, aes(x=MotherScore, y=Blue_score, color=region)) +
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

# MOTHER-MEAN-OFFSPRING REGRESSION
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

## HERITABILITY OF BLUENESS BY temp
meanscoresbytemp <- bluefemales %>%
  group_by(MotherID,temp) %>%
  summarise(Blue_score = mean(Blue_score), MotherScore = mean(MotherScore))
print(n=42,meanscoresbytemp)
meanscoresbytemp$bluenesschange <- meanscoresbytemp$Blue_score - meanscoresbytemp$MotherScore
coldmeans <- subset(meanscoresbytemp, temp==18)
warmmeans <- subset(meanscoresbytemp, temp==26)

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


## HERITABILITY OF BLUENESS BY POP
meanscoresbypop <- bluefemales %>%
  group_by(MotherID,region) %>%
  summarise(Blue_score = mean(Blue_score), MotherScore = mean(MotherScore))
print(n=42,meanscoresbypop)
meanscoresbypop$bluenesschange <- meanscoresbypop$Blue_score - meanscoresbypop$MotherScore
omeans <- subset(meanscoresbypop, region=='O')
smeans <- subset(meanscoresbypop, region=='S')

plot(omeans$MotherScore, omeans$Blue_score, 
     xlim=c(1,5), ylim=c(1,5),
     xlab="Mother Score", ylab="Mean Offspring Score",
     main="Oland Mother-Offspring Regression")
k <- lm(Blue_score ~ MotherScore, data=omeans)
abline(k, lty=2)
summary(k)
legend("bottomleft", legend=c(expression(paste(beta, " = -0.28")),
                          "p = 0.284"), bty="n")
#not significant, sample size?

plot(smeans$MotherScore, smeans$Blue_score, 
     xlim=c(1,5), ylim=c(1,5),
     xlab="Mother Score", ylab="Mean Offspring Score",
     main="Skane Mother-Offspring Regression")
v <- lm(Blue_score ~ MotherScore, data=smeans)
abline(v, lty=2)
summary(v)
legend("topleft", legend=c(expression(paste(beta, " = -0.003")),
                          "p = 0.99"), bty="n")
#not significant

par(mfrow=c(1,1))

#=====================================================================

## CHANGE IN BLUENESS BY MOTHER
ggplot(data=meanscoresbytemp, aes(x=temp, y=Blue_score, color=MotherID)) +
  geom_point() +
  geom_jitter(width=0.5) +
  geom_line(aes(x=temp, y=Blue_score, group=MotherID)) +
  ggtitle("Change in Mean Blueness by Mother") +
  theme_bw()

# Change in blueness from mother
ggplot(data=meanscoresbytemp, aes(x=MotherID, y=bluenesschange, fill=temp)) +
  geom_bar(position="dodge", stat="identity")
barplot(height=meanscoresbytemp$bluenesschange, beside=TRUE, col=c("blue","red"))

#wombat - g-matrix?


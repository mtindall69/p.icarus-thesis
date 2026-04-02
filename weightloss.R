# Polyommatus Weight Loss and Metabolic Cost of Color

#clear environments
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm,grid,gridExtra)

#load data
blueoffspring = read_excel("Blue offspring.xlsx") # all 610 offspring
bluewings = read_excel("Blue pupae.xlsx") # wing areas n=470

blueoffspring$motherID = as.character(blueoffspring$motherID)
blueoffspring$region = as.factor(blueoffspring$region)
blueoffspring$temp = as.factor(blueoffspring$temp)
blueoffspring$sex = as.factor(blueoffspring$sex)
blueoffspring$temp_label <- ifelse(blueoffspring$temp == 18, "Cold (18°C)", "Warm (26°C)")
blueoffspring$region_label  <- ifelse(blueoffspring$region == "O", "Öland", "Skåne")
blueoffspring$sex_label <- ifelse(blueoffspring$sex == "F", "Female", "Male")

#subset data
blueFdata <- subset(blueoffspring, sex=="F")
blueMdata <- subset(blueoffspring, sex=="M")

bluewings$motherID = as.character(bluewings$motherID)
bluewings$region = as.factor(bluewings$region)
bluewings$temp = as.factor(bluewings$temp)
bluewings$sex = as.factor(bluewings$sex)
bluewings$temp_label <- ifelse(bluewings$temp == 18, "Cold (18°C)", "Warm (26°C)")
bluewings$region_label  <- ifelse(bluewings$region == "O", "Öland", "Skåne")
bluewings$sex_label <- ifelse(bluewings$sex == "F", "Female", "Male")

#subset data
blueFwings <- subset(bluewings, sex=="F")
blueMwings <- subset(bluewings, sex=="M")
#======================================================================

#subset data into regions and temperatures, SEXES
#FEMALES
ofemales = subset(blueFdata, region=='O')
sfemales = subset(blueFdata, region=='S')
coldfemales = subset(blueFdata, temp=='18')
warmfemales = subset(blueFdata, temp=='26')
ocold = subset(coldfemales, region=='O')
owarm = subset(warmfemales, region=='O')
scold = subset(coldfemales, region=='S')
swarm = subset(warmfemales, region=='S')

#ALL
ooff = subset(blueoffspring, region=='O')
soff = subset(blueoffspring, region=='S')
coldoff = subset(blueoffspring, temp=='18')
warmoff = subset(blueoffspring, temp=='26')
allocold = subset(coldoff, region=='O')
allowarm = subset(warmoff, region=='O')
allscold = subset(coldoff, region=='S')
allswarm = subset(warmoff, region=='S')

#MALES
omales = subset(blueMdata, region=='O')
smales = subset(blueMdata, region=='S')
coldmales = subset(blueMdata, temp=='18')
warmmales = subset(blueMdata, temp=='26')
ocoldm = subset(coldmales, region=='O')
owarmm = subset(warmmales, region=='O')
scoldm = subset(coldmales, region=='S')
swarmm = subset(warmmales, region=='S')


#======================================================================

#subset data into regions and temperatures, SEXES
#FEMALES
wofemales = subset(blueFwings, region=='O')
wsfemales = subset(blueFwings, region=='S')
wcoldfemales = subset(blueFwings, temp=='18')
wwarmfemales = subset(blueFwings, temp=='26')
wocold = subset(wcoldfemales, region=='O')
wowarm = subset(wwarmfemales, region=='O')
wscold = subset(wcoldfemales, region=='S')
wswarm = subset(wwarmfemales, region=='S')

#ALL
wooff = subset(bluewings, region=='O')
wsoff = subset(bluewings, region=='S')
wcoldoff = subset(bluewings, temp=='18')
wwarmoff = subset(bluewings, temp=='26')
wallocold = subset(wcoldoff, region=='O')
wallowarm = subset(wwarmoff, region=='O')
wallscold = subset(wcoldoff, region=='S')
wallswarm = subset(wwarmoff, region=='S')

#MALES
womales = subset(blueMwings, region=='O')
wsmales = subset(blueMwings, region=='S')
wcoldmales = subset(blueMwings, temp=='18')
wwarmmales = subset(blueMwings, temp=='26')
wocoldm = subset(wcoldmales, region=='O')
wowarmm = subset(wwarmmales, region=='O')
wscoldm = subset(wcoldmales, region=='S')
wswarmm = subset(wwarmmales, region=='S')

#================================================================
#################################
# MEAN WEIGHTS
#################################

#ALL OFFSPRING MEAN SIZES
mean(blueoffspring$adult_weight_g) # avg for all offspring is 0.02290g
mean(ooff$adult_weight_g) # 0.02243g
mean(soff$adult_weight_g) # 0.02357g
  #Skane slightly heavier
mean(coldoff$adult_weight_g) # 0.02595g
mean(warmoff$adult_weight_g) # 0.01875g
  #Warmer is smaller, 26 << 18
mean(allocold$adult_weight_g) # 0.02566g
mean(allscold$adult_weight_g) # 0.02640g
  #Skane slightly heavier in 18
mean(allowarm$adult_weight_g) # 0.01754g
mean(allswarm$adult_weight_g) # 0.02022g
  #Skane larger in 26


#FEMALE Mean sizes
mean(blueFdata$adult_weight_g) # avg for all females is 0.02315g
mean(ofemales$adult_weight_g) # 0.02316g
mean(sfemales$adult_weight_g) # 0.02314g
  #females about same between regions, Oland marginally heavier
mean(coldfemales$adult_weight_g) # 0.02569g
mean(warmfemales$adult_weight_g) # 0.01912g
  #Warmer is lighter, 26 << 18
mean(ocold$adult_weight_g) # 0.02579g
mean(scold$adult_weight_g) # 0.02554g 
  #Females about same in cold, Oland marginally heavier in 18
mean(owarm$adult_weight_g) # 0.01824g
mean(swarm$adult_weight_g) # 0.01995g
  #Skane heavier in 26, Oland size more sensitive to warm
# significantly different? SD, t-test, chi square?

#MALE Mean sizes
mean(blueMdata$adult_weight_g) # avg for all males is 0.02259g 
  # Average male is lighter than female
mean(omales$adult_weight_g) # 0.02168g
mean(smales$adult_weight_g) # 0.02424g
  # Skane males larger than Oland
mean(coldmales$adult_weight_g) # 0.02632g
mean(warmmales$adult_weight_g) # 0.01837g
  #Warmer is lighter, 26 << 18
    # Males larger than females in cold but smaller in warm
    # Male size less constrained?
mean(ocoldm$adult_weight_g) # 0.02551g
mean(scoldm$adult_weight_g) # 0.02798g
  #Skane heavier in 18
mean(owarmm$adult_weight_g) # 0.01699g
mean(swarmm$adult_weight_g) # 0.02059g
  #Skane much heavier in 26

#Warmer always lighter, Skane generally heavier BUT
#Females about same in cold, Skane is heavier in warm
#Skane males are always heavier but much heavier in warm
      # Oland size more sensitive to warm
#Males generally smaller BUT
#slightly larger than females in cold, slightly smaller in warm
#Skane males are bigger than Skane females, larger in cold, about same in warm
#Oland males are smaller than Oland females, about same in cold, smaller in warm

#################################
# MEAN WING AREAS
#################################

#ALL OFFSPRING MEAN WING AREAS
mean(bluewings$aTWA) # avg for all offspring is 61.8664mm^2
mean(wooff$aTWA) # 58.78866mm^2
mean(wsoff$aTWA) # 65.74318mm^2
  #Skane larger
mean(wcoldoff$aTWA) # 64.75711mm^2
mean(wwarmoff$aTWA) # 57.75386mm^2
  #Warmer is smaller, 26 << 18
mean(wallocold$aTWA) # 62.63478mm^2
mean(wallscold$aTWA) # 67.51613mm^2
  #Skane larger in 18
mean(wallowarm$aTWA) # 53.12832mm^2
mean(wallswarm$aTWA) # 63.32552mm^2
  #Skane larger in 26


#FEMALE Mean sizes
mean(blueFwings$aTWA) # avg for all females is 61.15329mm^2
mean(wofemales$aTWA) # 58.43012mm^2
mean(wsfemales$aTWA) # 64.28853mm^2
  #Skane larger
mean(wcoldfemales$aTWA) # 63.43108mm^2
mean(wwarmfemales$aTWA) # 57.56623mm^2
  #Warmer is smaller, 26 << 18
mean(wocold$aTWA) # 61.58016mm^2
mean(wscold$aTWA) # 65.83514mm^2
  #Skane larger in 18
mean(wowarm$aTWA) # 52.6889mm^2
mean(wswarm$aTWA) # 62.21845mm^2
  #Skane much larger in 26, Oland size more sensitive to warm


#MALE Mean sizes
mean(blueMwings$aTWA) # avg for all males is 63.49708mm^2
  # Average male is larger than female
mean(womales$aTWA) # 59.50986mm^2
mean(wsmales$aTWA) # 69.69151mm^2
  # Skane males larger than Oland
mean(wcoldmales$aTWA) # 68.24665mm^2
mean(wwarmmales$aTWA) # 58.1095mm^2
  #Warmer is smaller, 26 << 18
# Males much larger than females in cold but about same in warm
# Male size less constrained?
mean(wocoldm$aTWA) # 65.40623mm^2
mean(wscoldm$aTWA) # 71.94781mm^2
  #Skane larger in 18
mean(wowarmm$aTWA) # 53.74751mm^2
mean(wswarmm$aTWA) # 66.45419mm^2
  #Skane much larger in 26
#Oland less constrained

#Warmer always smaller, Skane always larger BUT
#Skane females slightly larger in cold, Skane males larger in cold
#Males always larger BUT about same in warm
#Skane males are always larger but Oland sexes about same in warm
  #Oland size more sensitive to warm


#===============================================================
# VISUALIZATION
#===============================================================

PAL_TEMP <- c("Cold (18°C)" = "#3B7DD8", "Warm (26°C)" = "#E8712A")
PAL_REGION  <- c("Öland" = "#2CA02C", "Skåne" = "#9467BD")
TEMP_LABELS <- c("Cold (18°C)", "Warm (26°C)")
PAL_SEX <- c("Female" = "#D81B60", "Male" = "#1C05B3")

p1 <- ggplot(blueoffspring, aes(x = region_label, y = adult_weight_g, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP, guide = "none") +
  labs(x = "Region", y = "Adult weight (g)",
       title = "A) Adult Weight by Region") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- ggplot(blueoffspring, aes(x = sex_label, y = adult_weight_g, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP) +
  labs(x = "Sex", y = "Adult weight (g)",
       title = "B) Adult Weight by Sex", fill = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p3 <- ggplot(blueoffspring, aes(x = region_label, y = adult_weight_g, fill = sex_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_SEX) +
  labs(x = "Region", y = "Adult weight (g)",
       title = "C) Adult Weight- Sex x Region", fill = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, p3, ncol = 3)

p1 <- ggplot(blueoffspring, aes(x = region_label, y = per_weight_lost, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP, guide = "none") +
  labs(x = "Region", y = "Weight lost (%)",
       title = "A) Weight loss by Region") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- ggplot(blueoffspring, aes(x = sex_label, y = per_weight_lost, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP) +
  labs(x = "Sex", y = "Weight lost (%)",
       title = "B) Weight loss by Sex", fill = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p3 <- ggplot(blueoffspring, aes(x = region_label, y = per_weight_lost, fill = sex_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_SEX) +
  labs(x = "Region", y = "Weight lost (%)",
       title = "C) Weight loss- Sex x Region", fill = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, p3, ncol = 3)

p1 <- ggplot(bluewings, aes(x = region_label, y = aTWA, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP, guide = "none") +
  labs(x = "Region", y = "Average Wing Area (mm^2)",
       title = "A) Wing Area by Region") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p2 <- ggplot(bluewings, aes(x = sex_label, y = aTWA, fill = temp_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_TEMP) +
  labs(x = "Sex", y = "Average Wing Area (mm^2)",
       title = "B) Wing Area by Sex", fill = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p3 <- ggplot(bluewings, aes(x = region_label, y = aTWA, fill = sex_label)) +
  geom_boxplot(outlier.size = 1, linewidth = 0.5) +
  scale_fill_manual(values = PAL_SEX) +
  labs(x = "Region", y = "Average Wing Area (mm^2)",
       title = "C) Wing Area- Sex x Region", fill = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, p3, ncol = 3)

p1 <- ggplot(blueFwings, aes(x = aTWA, y = aTBWA, colour = temp_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_TEMP) +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       title = "Female wing size vs. blueness", colour = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.13, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p2 <- ggplot(blueFwings, aes(x = aTWA, y = aTBWA, colour = region_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_REGION) +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       title = "B) Female wing size vs. blueness", colour = "Region") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p3 <- ggplot(blueMwings, aes(x = aTWA, y = aTBWA, colour = temp_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_TEMP) +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       title = "C) Male wing size vs. blueness", colour = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, p3, ncol = 3)

p2 <- ggplot(wofemales, aes(x = aTWA, y = aTBWA, colour = temp_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_TEMP) +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       title = "B) Oland wing size vs. blueness", colour = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p3 <- ggplot(wsfemales, aes(x = aTWA, y = aTBWA, colour = temp_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = PAL_TEMP) +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       title = "C) Skane wing size vs. blueness", colour = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, p3, ncol = 3)


# Does weight predict wing size?
p1 <- ggplot(bluewings, aes(x = adult_weight_g, y = aTWA, colour = temp_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  scale_colour_manual(values = PAL_TEMP) +
  labs(x = "Adult weight (g)", y = "Average wing area (mm²)",
       title = "Adult weight vs. wing size", colour = "Temperature") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

p2 <- ggplot(bluewings, aes(x = adult_weight_g, y = aTWA, colour = sex_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
  scale_colour_manual(values = PAL_SEX) +
  labs(x = "Adult weight (g)", y = "Average wing area (mm²)",
       title = "Adult weight vs. wing size", colour = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.85),
        legend.background = element_rect(colour = "grey80"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

combined <- grid.arrange(p1, p2, ncol = 2)

m <- lm(aTWA ~ adult_weight_g + sex + temp, data=bluewings)
summary(m) # yes


# WEIGHT LOSS AND METABOLIC COST OF COLOR

par(mfrow=c(2,1))

#Females
plot(blueFdata$pupation_weight_g, blueFdata$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Female Pupa Weight to adult Weight")
m <- lm(blueFdata$adult_weight_g ~ blueFdata$pupation_weight_g)
abline(m, col="red")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.458g adult/g pupa, adult about 46% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.46")),"p < 0.01"))

#Males
plot(blueMdata$pupation_weight_g, blueMdata$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Male Pupa Weight to adult Weight")
m <- lm(blueMdata$adult_weight_g ~ blueMdata$pupation_weight_g)
abline(m, col="blue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.429g adult/g pupa, adult about 43% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.43")),"p < 0.01"))
  #Males lose marginally more weight, male weight not as important? significant?

par(mfrow=c(1,1))

#Offspring
colors = c("red", "blue")
plot(blueoffspring$pupation_weight_g, blueoffspring$adult_weight_g, 
     col=colors[factor(blueoffspring$sex)], pch=20,
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Pupa Weight to adult Weight")
m <- lm(blueoffspring$adult_weight_g ~ blueoffspring$pupation_weight_g)
abline(m)
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.439g adult/g pupa, adult about 44% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.44")),"p < 0.01"))
legend("topleft", legend=c("Female", "Male"), pch=20, 
       col=c("red","blue"))

m <- lm(adult_weight_g ~ pupation_weight_g + sex, data=blueoffspring)
summary(m) #sex is significant but effect size is small

#FOR REGIONS AND tempS

par(mfrow=c(2,1))

#Cold
plot(coldoff$pupation_weight_g, coldoff$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Cold Pupa Weight to adult Weight")
m <- lm(coldoff$adult_weight_g ~ coldoff$pupation_weight_g)
abline(m, col="blue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.418g adult/g pupa, adult about 42% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.42")),"p < 0.01"))

#Warm
plot(warmoff$pupation_weight_g, warmoff$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Warm Pupa Weight to adult Weight")
m <- lm(warmoff$adult_weight_g ~ warmoff$pupation_weight_g)
abline(m, col="red")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.405g adult/g pupa, adult about 40% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.40")),"p < 0.01"))

#Oland
plot(ooff$pupation_weight_g, ooff$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Oland Pupa Weight to adult Weight")
m <- lm(ooff$adult_weight_g ~ ooff$pupation_weight_g)
abline(m, col="steelblue")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.452g adult/g pupa, adult about 46% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.45")),"p < 0.01"))

#Skane
plot(soff$pupation_weight_g, soff$adult_weight_g, 
     xlim=c(0,0.11), ylim=(c(0,0.05)),
     xlab="Pupa Weight (g)", ylab="adult Weight (g)",
     main="Skane Pupa Weight to adult Weight")
m <- lm(soff$adult_weight_g ~ soff$pupation_weight_g)
abline(m, col="orange")
abline(a=0, b=1, lty=2)
summary(m)$coef # 0.429g adult/g pupa, adult about 43% pupal weight
legend("bottomright", legend=c(expression(paste(beta, " = 0.43")),"p < 0.01"))

par(mfrow=c(1,1))


#=========================================================
# Mixed Models - 
# IS THERE A METABOLIC COST TO PRODUCE DIFFERENT COLORS?
#=========================================================

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
plot(allEffects(mw2))

mw2 <- lmer(adult_weight_g ~ pupation_weight_g + temp + avg_prop_blue + 
              pupation_weight_g:temp + pupation_weight_g:avg_prop_blue + 
              (1|motherID), data=blueFdata)
summary(mw2)
Anova(mw2) # avg_prop_blue NOT significant nor any interactions
plot(allEffects(mw2))

mw2 <- lmer(adult_weight_g ~ pupation_weight_g + temp + avg_prop_blue + 
              pupation_weight_g:temp + pupation_weight_g:avg_prop_blue + 
              temp:avg_prop_blue + (1|motherID), data=blueFdata)
summary(mw2)
Anova(mw2) # avg_prop_blue NOT significant nor any interactions
plot(allEffects(mw2))


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
           data=blueoffspring)
summary(mw)
Anova(mw) # sex significant, pupation_weight:sex weakly significant
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex + temp + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw) # YES
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * temp + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw) # no interactions significant
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * region + temp + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw) 
plot(allEffects(mw))

mw <- lmer(adult_weight_g ~ pupation_weight_g * sex * temp * region + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw)
plot(allEffects(mw))


mw <- lmer(adult_weight_g ~ pupation_weight_g * sex + temp * region + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw) # All effects significant
plot(allEffects(mw)) # pupation weight:sex and temp:region interactions significant


mw <- lmer(pupation_weight_g ~  sex * temp * region + (1|motherID), 
           data=blueoffspring)
summary(mw)
Anova(mw) 
plot(allEffects(mw))


# Weight loss

m <- lmer(weight_loss_g ~ pupation_weight_g + sex + dev_length + (1|motherID), data=blueoffspring)
summary(m)
Anova(m) 


# Polyommatus icarus project: Blue area analyses

#clear environment
rm(list=ls())

#set directory
setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis/")

#load packages
pacman::p_load(readxl,dplyr,tidyverse,ggplot2,ggimage,lme4,car,effects,
               RColorBrewer,ape,MCMCglmm,psych,effects,MASS,MuMIn,glmmTMB,rlang,
               interactions,ggfortify,sjPlot,DHARMa,performance)

PAL_TEMP <- c("Cold (18°C)" = "#3B7DD8", "Warm (26°C)" = "#E8712A")
PAL_REGION  <- c("Öland" = "#2CA02C", "Skåne" = "#9467BD")
TEMP_LABELS <- c("Cold (18°C)", "Warm (26°C)")
PAL_SEX <- c("Female" = "#D81B60", "Male" = "#1C05B3")

#==================================================================
# LOAD AND PREPROCESS DATA
#==================================================================

#load data
bluedata <- read.csv("blueness.csv")

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
                  "start_day", "pupa_day", "adult_day", "pupation_length", "dev_length")
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
# FAMILY MEANS PER TEMP
#=======================================================================
family_means <- blueFdata %>%
  group_by(motherID, temp) %>%
  summarise(N = n(), Mean = mean(avg_prop_blue), SD = sd(avg_prop_blue),
            .groups = "drop") %>%
  arrange(motherID, temp) %>%
  as.data.frame()
write.csv(family_means, "propmeans.csv", row.names = FALSE)


#=======================================================================
# MODELS
#=======================================================================

# FEMALE BLUENESS
####################

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

#Add areaxtemp interaction - Friberg et al (2025)
m1 <- lmer(avg_blue_mm ~ avg_total_mm + temp + region + avg_total_mm:temp + 
             temp:region + (1|motherID), data=blueFdata)
summary(m1)
Anova(m1) # ALL VERY SIGNIFICANT
AIC(m1) #2147.6
plot(allEffects(m1))


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
ef<-effect("temp:motherID", m2)
summary(ef) # output for reaction norms


####################
# MALE BLUENESS
####################

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
Anova(m2b) 
AIC(m2b) #637.2989
plot(allEffects(m2b))

        
###################
#BOTH SEXES
###################

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


#################
# TOTAL AREA
#################

# Females
m <- lm(avg_total_mm ~ temp + region + motherID + temp:region + temp:motherID, 
        data=blueFdata)
summary(m)

m <- lmer(avg_total_mm ~ temp + region + motherID + temp:region + temp:motherID 
          + (1|motherID), data=blueFdata)
summary(m)

ma <- lmer(avg_total_mm ~ temp * region + (1|motherID), data=blueFdata)
summary(ma)
Anova(ma) # temp significant to size
plot(allEffects(ma))
ef<-effect("temp:region", ma)
summary(ef)



# Both
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


# Males
m2a <- lmer(avg_total_mm ~ temp * region + (1|motherID)+ (1|start_day), data=blueMdata)
summary(m2a)
Anova(m2a) # temp significant to size
plot(allEffects(m2a))
ef<-effect("temp:region", m2a)
summary(ef)


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

# area vs blue
summary(lm(avg_blue_mm ~ avg_total_mm + temp, data = blueFdata))
# adjusted R^2 = 0.213

summary(lm(avg_blue_mm ~ avg_total_mm, 
           data = subset(blueFdata, temp_label == "Cold (18°C)")))
summary(lm(avg_blue_mm ~ avg_total_mm,
           data = subset(blueFdata, temp_label == "Warm (26°C)")))
  # no relationship in warm

r2_cold <- summary(lm(avg_blue_mm ~ avg_total_mm, 
                      data = subset(blueFdata, temp_label == "Cold (18°C)")))$adj.r.squared
r2_warm <- summary(lm(avg_blue_mm ~ avg_total_mm,
                      data = subset(blueFdata, temp_label == "Warm (26°C)")))$adj.r.squared

cold_lbl <- paste0("Cold (18°C)  R² = ", round(r2_cold, 3))
warm_lbl <- paste0("Warm (26°C)  R² = ", round(r2_warm, 3))

blueFdata$temp_r2_label <- factor(
    ifelse(blueFdata$temp_label == "Cold (18°C)", cold_lbl, warm_lbl),
    levels = c(cold_lbl, warm_lbl)
    )

pal_temp_r2 <- setNames(c("#3B7DD8", "#E8712A"), c(cold_lbl, warm_lbl))

pa <- ggplot(blueFdata, aes(x = avg_total_mm, y = avg_blue_mm, colour = temp_r2_label,
                      fill = temp_r2_label, linetype = temp_r2_label)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
  scale_colour_manual(values = pal_temp_r2) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(cold_lbl, warm_lbl)),
                        guide = "none") +
  guides(
    colour = guide_legend(override.aes = list(fill = PAL_TEMP, alpha = 0.15,
                                              linetype = c("solid", "dashed"))),
    fill   = "none"
  ) +
  scale_fill_manual(values = pal_temp_r2, guide = "none") +
  labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
       colour = "Temperature") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.20, 0.90),
        legend.background = element_rect(colour = "grey80"),
        axis.title = element_text(size=13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave(file.path(PLOT_DIR,"wingvblue.png"), pa, width = 6, height = 5, dpi = 200)

#=====================================================================
# FITNESS
#=====================================================================

#Fecundity
momfit <- read_excel("momfit.xlsx")
momfit$motherID = as.factor(momfit$motherID)
momfit$region = as.factor(momfit$region)
momfit$motherscore = as.factor(momfit$motherscore)
momfit$eggs_per_day <- momfit$total_eggs / momfit$days_alive
momfit$z_eggs <- scale(momfit$total_eggs)
momfit$region_label  <- ifelse(momfit$region == "O", "Öland", "Skåne")
momfit$region_label = as.factor(momfit$region_label)

m <- lm(total_eggs ~ motherscore*mom_aTWA*days_alive*region, data=momfit)
summary(m) # AIC 710, motherscore significant 
library(ggfortify)
autoplot(m)

m <- lm(eggs_per_day ~ motherscore, data=momfit)
summary(m) # AIC 461, near significant

# GLM for non-normal count data
m <- glm(total_eggs ~ motherscore, data=momfit, family=poisson)
summary(m) # AIC 7066, motherscore very significant

m <- glm(total_eggs ~ motherscore + mom_aTWA, data=momfit, family=poisson)
summary(m) # AIC 7057
# motherscore and area significant

m <- glm(total_eggs ~ motherscore + mom_aTWA + days_alive, data=momfit, family=poisson)
summary(m) # AIC 5951
# all effects very significant

m <- glm(total_eggs ~ motherscore + mom_aTWA + days_alive + region, data=momfit, family=poisson)
summary(m) # AIC 4770
#motherscore not significant, region is

m <- glm(total_eggs ~ mom_aTWA + days_alive + region, data=momfit, family=poisson)
summary(m) # AIC 4768 not much better without mother score

m <- glm(total_eggs ~ motherscore + mom_aTWA + region, data=momfit, family=poisson)
summary(m) # AIC 5199
# all effects significant but model worse without days alive

m <- glm(total_eggs ~ region, data=momfit, family=poisson)
summary(m) # AIC 5635, better than just motherscore, less residual deviance

m <- glm(total_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region, data=momfit, family=poisson)
summary(m) # AIC 4751, high R^2 values, all effects significant

m <- glm(total_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region, data=momfit, family=poisson)
summary(m) # AIC 4741, region no longer significant

m <- glm(total_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region + motherscore:mom_aTWA, 
         data=momfit, family=poisson)
summary(m) # AIC 4713, region and score:region no longer significant

m <- glm(total_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region + motherscore:mom_aTWA +
           days_alive:region, 
         data=momfit, family=poisson)
summary(m) # AIC 4310, days_alive near significant, all others significant

m <- glm(total_eggs ~ motherscore * mom_aTWA * days_alive * region, data=momfit, family=poisson)
summary(m) #3041.5

m2 <- glmer(total_eggs ~ motherscore * mom_aTWA * days_alive * region + (1|motherID),
            data=momfit, family=poisson, control_params)
control_params <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))

m3 <- glmmTMB(total_eggs ~ motherscore * region_label + days_alive + (1|motherID),
              data=momfit, ziformula = ~1, family=nbinom1())

m3a <- update(m3,.~.-motherscore:region_label)
anova(m3a,m3)
  
m3b <- update(m3a,.~.+motherscore:days_alive)
anova(m3a,m3b)

library(DHARMa)
simout <- simulateResiduals(fittedModel=m)
plot(simout)
testOutliers(simout)
testOverdispersion(simout)
testZeroInflation(simout)

Anova(m3b)
summary(m3b) #motherscore near significant

library(sjPlot)
library(interactions)

cat_plot(m3b, data=momfit, modx=region_label, pred=days_alive, mod2=motherscore, interval.geom=c("linerange"))

p <- plot_model(m3b, terms=c("motherscore", "region_label"), show.data=TRUE, type="pred") +
  geom_line() +
  labs(x = "Mother Score", y = "Total Eggs", title=NULL, colour = "Region") +
  scale_color_manual(values = PAL_REGION) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        legend.text = element_text(size = 9))
ggsave(file.path(PLOT_DIR, "fecundityxblue.png"), p,
       width = 7, height = 5, dpi = 200)

#z-scale eggs?
m <- glm(z_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region, data=momfit, family=gaussian)
summary(m) # AIC 139, motherscore, region not significant, R^2 0.42

m <- glm(z_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region + motherscore:mom_aTWA, 
         data=momfit, family=gaussian)
summary(m) # AIC 140, R^2 0.42

m <- glm(z_eggs ~ motherscore + days_alive + mom_aTWA + region +
           motherscore:region + mom_aTWA:region + motherscore:mom_aTWA +
           days_alive:region, 
         data=momfit, family=gaussian)
summary(m) # AIC 140, no significance, R^2 0.45 ---------- BEST?

m <- glm(z_eggs ~ days_alive + mom_aTWA + region, data=momfit, family=gaussian)
summary(m) # AIC 135, R^2 0.40

m <- glm(z_eggs ~ motherscore + days_alive + mom_aTWA + region, data=momfit, family=gaussian)
summary(m) # AIC 137, R^2 0.40

r.squaredGLMM(m)
1- (m$deviance/m$null.deviance)

dev   <- deviance(m)
df    <- df.residual(m)
p_value <- 1-pchisq(dev,df)
cbind.data.frame(c("Deviance GOF", "D", "df", "p-value"),c(" ", round(dev,4), df, p_value))

pr <- sum(residuals(m, type="pearson")^2) # get Pearson Chi2
pchisq(pr, m$df.residual, lower=F) # calc p-value
dv <- m$deviance
pchisq(dv, m$df.residual, lower= F) # calc p-vl

P_disp <- function(x) {
  pr <- sum(residuals(x, type="pearson")^2)
  dispersion <- pr/x$df.residual
  cat("\n Pearson Chi2 = ", pr , "\n Dispersion = ", dispersion, "\n")
}
P_disp(m)

#not enough data to say if blueness is related to fecundity

plot(momfit$motherscore, momfit$total_eggs, 
     xlab="Mother Score", ylab="Total Eggs", main="Fecundity vs Blueness")
abline(lm(total_eggs ~ motherscore, data=momfit), lty=2)

# Simple bivariate model for plot annotation
  m_plot <- lm(total_eggs ~ motherscore, data = momfit)
  sm      <- summary(m_plot)
  slope   <- coef(sm)["motherscore", "Estimate"]
  slope_se <- coef(sm)["motherscore", "Std. Error"]
  r2      <- sm$r.squared
  p_val   <- coef(sm)["motherscore", "Pr(>|t|)"]

  label_txt <- paste0(
    "slope = ", round(slope, 2), " ± ", round(slope_se, 2), "\n",
    "R² = ",    round(r2, 3),    "\n",
    "p = ",     format.pval(p_val, digits = 3)
  )

p <- ggplot(momfit, aes(x = motherscore, y = total_eggs)) +
    geom_smooth(method = "lm", colour = "black", fill = "grey70", alpha = 0.4, lty=2) +
    geom_point(size = 2.5, alpha = 0.7) +
    annotate("label",
             x = 2.0, y = 800,
             label     = label_txt,
             hjust = 1.05, vjust = 1.3,
             size  = 3.5, label.padding = unit(0.4, "lines"), label.size = 0.3) +
    scale_x_continuous(breaks = 1:5) +
    labs(
      x = "Mother Blue Score",
      y = "Total Eggs",
      title = "Fecundity vs Blueness",
    ) +
    theme_bw(base_size = 13)+
    theme(plot.background = element_blank(),
          panel.background = element_blank())
  
  ggsave("fecundity.png", plot = p, bg = "transparent")
  
#============================================
#Surivival
survival <- read_excel("survival.xlsx")

survival$motherID = as.factor(survival$motherID)
survival$pupaeID = as.factor(survival$pupaeID)
survival$temp = as.factor(survival$temp)
survival$region = as.factor(survival$region)
survival$sex = as.factor(survival$sex)
#survival$motherscore = as.factor(survival$motherscore)
survival$region_label  <- ifelse(survival$region == "O", "Öland", "Skåne")
survival$temp_label <- ifelse(survival$temp == 18, "Cold (18°C)", "Warm (26°C)")
survival$temp_label = as.factor(survival$temp_label)
survival$region_label = as.factor(survival$region_label)


m <- glm(survived ~ motherscore, data = survival, family = binomial)
summary(m) # AIC:1820, motherscore is significant

m <- glm(survived ~ motherscore + temp + region, data = survival, family = binomial)
summary(m) # AIC:1792, motherscore near significant, temp strongly
#R^2= 0.01996

m <- glm(survived ~ motherscore * temp * region, data = survival, family = binomial)
summary(m) # AIC:1785


m <- glmmTMB(survived ~ motherscore + temp + region + (1|motherID), data = survival, family = binomial)
summary(m) # AIC:1694, motherscore NOT significant, temp strongly
r2(m)

m <- glmmTMB(survived ~ motherscore * temp * region + (1|motherID), data = survival, family = binomial)
summary(m) #AIC: 1680
r2(m)

m <- glmmTMB(survived ~ motherscore + temp_label + region_label + 
               motherscore:temp_label + temp_label:region_label + (1|motherID), 
             data = survival, family = binomial)
summary(m) #AIC: 1677
r2(m)
anova(m1,m2)
plot(allEffects(m))

m_cold_surv <- glm(survived ~ motherscore,
                   data = subset(survival, temp_label == "Cold (18°C)"),
                   family = binomial)
m_warm_surv <- glm(survived ~ motherscore,
                   data = subset(survival, temp_label == "Warm (26°C)"),
                   family = binomial)

r2_tjur(m_cold_surv)
r2_tjur(m_warm_surv)
r2_mcfadden(m_cold_surv)
r2_mcfadden(m_warm_surv)

cat_plot(m, data=survival, modx=temp_label, pred=motherscore, interval.geom=c("linerange"))

library(ggeffects)                                                               
pred_surv <- ggpredict(m, terms = c("motherscore [all]", "temp_label"))

p <- ggplot(pred_surv, aes(x = x, y = predicted, colour = group, fill = group,
                           linetype = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_jitter(data = survival,
              aes(x = motherscore, y = survived, colour = temp_label),
              inherit.aes = FALSE,
              width = 0.35, height = 0.03,
              size = 2.0, alpha = 0.4) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_fill_manual(values = PAL_TEMP, guide = "none") +
  scale_linetype_manual(values = c("Cold (18°C)" = "dashed", "Warm (26°C)" = "solid"),
                        guide = "none") +
  labs(x = "Mother Score", y = "Probability of Survival", title = NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.position = "top")

ggsave(file.path(PLOT_DIR, "survivalxblue.png"), p,
       width = 7, height = 5, dpi = 200)


m <- glm(survived ~ motherscore, data = survival, family = binomial)

coefs    <- summary(m)$coef
invlogit <- function(x) 1 / (1 + exp(-x))
threshold <- -coefs[1, 1] / coefs[2, 1]
cat("Motherscore at 50% survival probability:", round(threshold, 3), "\n")

# Prediction ribbon
x_pred  <- seq(min(survival$motherscore), max(survival$motherscore), by = 0.01)
pred    <- predict(m,
                   newdata = data.frame(motherscore = x_pred),
                   type    = "link",
                   se.fit  = TRUE)

ribbon  <- data.frame(
  motherscore = x_pred,
  fit         = invlogit(pred$fit),
  lwr         = invlogit(pred$fit - 1.96 * pred$se.fit),
  upr         = invlogit(pred$fit + 1.96 * pred$se.fit)
)

# ── Threshold at mean survival probability ─────────────────────────────────────
mean_survival <- mean(survival$survived)
mean_logit    <- log(mean_survival / (1 - mean_survival))   # logit of mean survival
threshold     <- (mean_logit - coefs[1, 1]) / coefs[2, 1]   # solve for x

cat("Mean survival probability:", round(mean_survival, 3), "\n")
cat("Motherscore at mean survival probability:", round(threshold, 3), "\n")

# Plotting
ggplot(survival, aes(x = motherscore, y = survived)) +
  
  # 95% CI ribbon
  geom_ribbon(data = ribbon,
              aes(x=motherscore, y = fit, ymin = lwr, ymax = upr),
              fill = "grey70", alpha = 0.4, inherit.aes = FALSE) +
  
  # Logistic curve
  geom_line(data = ribbon,
            aes(x=motherscore, y = fit),
            colour = "black", linewidth = 0.9, inherit.aes = FALSE) +
  
  # Reference lines at threshold
  geom_hline(yintercept = mean(survival$survived), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey40") +
  
  # Jittered raw data coloured by temperature
  geom_jitter(aes(colour = temp),
              width = 0.35, height = 0.03,
              size = 2.0, alpha = 0.4) +
  
  scale_colour_manual(values = c("#3B7DD8","#E8712A"),
                      labels = c("Cold", "Warm"),
                      name   = "Temperature") +
  
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  
  labs(
    subtitle = paste0("Logistic regression  |  Motherscore at mean survival (",
                      round(mean_survival * 100, 1), "%) = ",
                      round(threshold, 2)),
    x        = "Mother Blue Score",
    y        = "Probability of Survival"
  ) +
  
  theme_classic(base_size = 13) +
  theme(legend.position = "top")


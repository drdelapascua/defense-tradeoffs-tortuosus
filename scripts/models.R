### Danielle De La Pascua
### Ch2 Testing Intraspecifc Defense Trade-offs
### 1-28-24

### Linear Mixed Models

oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
install.packages("emmeans")
install.packages('effects')

#libraries
library(lme4)
library(nlme)
library(Matrix)
library(emmeans)
library(dplyr)
library(tidyr)
library(lattice)
install.packages("mgcv")
library(mgcv)
library(ggplot2)
library(effects)

#pull data
#paired_means <- read.csv("./data/paired_means.csv")
dw <-  read.csv("./data/dw.csv")
head(dw)

#question 1 - i want to know whether populations that have high underlying defenses invest more or less in induced defenses
# Q2 - growth-defense trade-offs & association with environment

#first build model

# random effects model
# fixed effect: compound 
dw$Population = as.factor(dw$Population)
b1 <- lmer(Butenyl_12.1 ~ Elevation*treatment*biomass + (1|Population), data = dw)
summary(b1)

# individual crossed effects

#running population as fixed effect crossed with treatment 

### 3MSO
x3mso_m1 <- lm(X3MSO_5.2 ~ Population*treatment, data = dw) 
summary(x3mso_m1) # a few significant interactions by population

# worth looking into iowa hill and a few other pops

#with elevation*treatment
x3mso_m2 <- lm(X3MSO_5.2 ~ Elevation*treatment, data = dw)
summary(x3mso_m2) # overall slight negative effect of elevation on 3MSO

### 4MSO
# pop*trt
x4mso_m1 <- lm(X4MSO_7.1 ~ Population*treatment, data = dw) 
summary(x4mso_m1) # no significant interactions by population

# worth looking into iowa hill and a few other pops

#with elevation*treatment
x4mso_m2 <- lm(X4MSO_7.1 ~ Elevation*treatment, data = dw)
summary(x4mso_m2) # overall slight negative effect of elevation on 3MSO

### Allyl
# pop*trt
allyl_m1 <- lm(Allyl_7.4 ~ Population*treatment, data = dw) 
summary(allyl_m1) # no significant interactions by population, SQ2 has more allyls than other pops

#with elevation*treatment
butenyl_m2 <- lm(Butenyl_12.1 ~ Elevation*treatment, data = dw)
summary(butenyl_m2) # overall positive effect of elevation on butenyls


### Butenyl
butenyl_m1 <- lm(Butenyl_12.1 ~ Population*treatment, data = dw) 
summary(butenyl_m1) # lots of significant interactions by population

#with elevation*treatment
butenyl_m2 <- lm(Butenyl_12.1 ~ Elevation*treatment, data = dw)
summary(butenyl_m2) # overall positive effect of elevation on butenyls

### I3M

#population*treatment
I3M_m1 <- lm(I3M_16.7 ~ Population*treatment, data = dw)
summary(I3M_m1) 

#with elevation*treatment
I3M_m2 <- lm(I3M_16.7 ~ Elevation*treatment, data = dw)
summary(I3M_m2)

### Indole

#population*treatment
Indole_m1 <- lm(Indole_18.8 ~ Population*treatment, data = dw)
summary(Indole_m1)  #here there are some significant interactions among some populations & treatments

# use EMMEANS to look at slope of global model

# use EMMEANS to look at values for the following populations:

# WL1
# YOSE8

#with elevation*treatment
Indole_m2 <- lm(Indole_18.8~ Elevation*treatment, data = dw)
summary(Indole_m2) # no significant interaction with elevation

# OLD CODE

#question 1 - i want to know whether populations that have high underlying defenses invest more or less in induced defenses

#Response variable: each glucosinolate compounds (indoles most inducible, I3m & Indole)
#Fixed effects:  treatment (categorical), elevation (continuous)
#Random effects: population, rack
#ANOVA: differences between intercepts of treatment lines & slope of the line
#If hypothesis is true, both slopes will be positive, the induced line will have a steeper slope than the control line across elevations, and the induced line will be above the control line. At high elevations, they should be significantly different, but may not be significantly different at low elevations. Populations at high elevations should have higher overall levels of the compound than lower elevations

#example model

#I3M
colnames(paired_means)

I3M_q1 <- lm(CW_I3M_16.7 ~ Control_I3M_16.7*Elevation.x, data = paired_means, na.action = na.exclude)
summary(I3M_q1) #interaction is not significant, also not significant if by population
I3M_q2 <- lm(CW_I3M_16.7 ~ Control_I3M_16.7 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(I3M_q2)
I3M_q3 <- lm(CW_I3M_16.7 ~  Control_I3M_16.7, data = paired_means, na.action = na.exclude)
summary(I3M_q3)

#plot to see the relationship
ggplot(data = paired_means, aes(x = Control_I3M_16.7, y = CW_I3M_16.7)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Compute EMMs
emm <- emmeans(I3M_q1, c(""))
# Perform post-hoc tests for interactions
emm_pairs_specific <- pairs(emm)
print(emm)
emm_pairs_specific

### Testing for Interactions between constitutive defense level & elevation

# 3MSO
x3mso_q1 <- lm(CW_3MSO_5.2 ~ Control_3MSO_5.2*Elevation.x, data = paired_means, na.action = na.exclude)
summary(x3mso_q1)

# OH-Alkenyl
ohalkenyl_q1 <- lm(CW_OH.Alkenyl_6 ~ Control_OH.Alkenyl_6*Elevation.x, data = paired_means, na.action = na.exclude)
summary(ohalkenyl_q1)

# 4MSO
x4mso_q1 <- lm(CW_4MSO_7.1 ~ Control_4MSO_7.1*Elevation.x, data = paired_means, na.action = na.exclude)
summary(x4mso_q1)

# Allyl
allyl_q1 <- lm(CW_Allyl_7.4 ~ Control_Allyl_7.4*Elevation.x, data = paired_means, na.action = na.exclude)
summary(allyl_q1) # significant

plot(allEffects(allyl_q1))

#use EMMEANS to extract slopes at each elevation - something different happening across elevation

# 5MSO
x5mso_q1 <- lm(CW_5MSO_10.2 ~ Control_5MSO_10.2*Elevation.x, data = paired_means, na.action = na.exclude)
summary(x5mso_q1)

# Butenyl
butenyl_q1 <- lm(CW_Butenyl_12.1 ~ Control_Butenyl_12.1*Elevation.x, data = paired_means, na.action = na.exclude)
summary(butenyl_q1)

# 3MT
x3mt_q1 <- lm(CW_3MT_13.6 ~ Control_3MT_13.6*Elevation.x, data = paired_means, na.action = na.exclude)
summary(x3mt_q1)

# MSOO
msoo_q1 <- lm(CW_MSOO_13.8 ~ Control_MSOO_13.8*Elevation.x, data = paired_means, na.action = na.exclude)
summary(msoo_q1)

# OH-I3M
ohi3m_q1 <- lm(CW_OH.I3M_15.1 ~ Control_OH.I3M_15.1*Elevation.x, data = paired_means, na.action = na.exclude)
summary(ohi3m_q1)

# 4MT
x4mt_q1 <- lm(CW_4MT._15.5 ~ Control_4MT._15.5*Elevation.x, data = paired_means, na.action = na.exclude)
summary(x4mt_q1)

# Flavonol-16
flavonol16_q1 <- lm(CW_Flavonol_16.1 ~ Control_Flavonol_16.1*Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol16_q1)

# I3M
I3M_q1 <- lm(CW_I3M_16.7 ~ Control_I3M_16.7*Elevation.x, data = paired_means, na.action = na.exclude)
summary(I3M_q1)

# Flavonol-17
flavonol17_q1 <- lm(CW_Flavonol_17.5 ~ Control_Flavonol_17.5*Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol17_q1)

# Flavonol-18
flavonol18_q1 <- lm(CW_Flavonol_18.5 ~ Control_Flavonol_18.5*Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol18_q1)

# Indole 
indole_q1 <- lm(CW_Indole_18.8 ~ Control_Indole_18.8*Elevation.x, data = paired_means, na.action = na.exclude)
summary(indole_q1)

### Models without interactions 

# 3MSO
x3mso_q1.2 <- lm(CW_3MSO_5.2 ~ Control_3MSO_5.2 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(x3mso_q1.2)

# OH-Alkenyl
ohalkenyl_q1.2 <- lm(CW_OH.Alkenyl_6 ~ Control_OH.Alkenyl_6 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(ohalkenyl_q1.2)

# 4MSO
x4mso_q1.2 <- lm(CW_4MSO_7.1 ~ Control_4MSO_7.1 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(x4mso_q1.2)

# 5MSO
x5mso_q1.2 <- lm(CW_5MSO_10.2 ~ Control_5MSO_10.2 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(x5mso_q1.2)

# Butenyl
butenyl_q1.2 <- lm(CW_Butenyl_12.1 ~ Control_Butenyl_12.1 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(butenyl_q1.2)

# 3MT
x3mt_q1.2 <- lm(CW_3MT_13.6 ~ Control_3MT_13.6  + Elevation.x, data = paired_means, na.action = na.exclude)
summary(x3mt_q1.2)

# MSOO
msoo_q1.2 <- lm(CW_MSOO_13.8 ~ Control_MSOO_13.8 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(msoo_q1.2)

# OH-I3M
ohi3m_q1.2 <- lm(CW_OH.I3M_15.1 ~ Control_OH.I3M_15.1 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(ohi3m_q1.2)

# 4MT
x4mt_q1.2 <- lm(CW_4MT._15.5 ~ Control_4MT._15.5 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(x4mt_q1.2)

# Flavonol-16
flavonol16_q1.2 <- lm(CW_Flavonol_16.1 ~ Control_Flavonol_16.1 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol16_q1.2)

# I3M
I3M_q1.2 <- lm(CW_I3M_16.7 ~ Control_I3M_16.7 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(I3M_q1.2)

# Flavonol-17
flavonol17_q1.2 <- lm(CW_Flavonol_17.5 ~ Control_Flavonol_17.5 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol17_q1.2)

# Flavonol-18
flavonol18_q1.2 <- lm(CW_Flavonol_18.5 ~ Control_Flavonol_18.5 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(flavonol18_q1.2)

# Indole 
indole_q1.2 <- lm(CW_Indole_18.8 ~ Control_Indole_18.8 + Elevation.x, data = paired_means, na.action = na.exclude)
summary(indole_q1.2)

### plot significant effects

# Butenyl
ggplot(data = paired_means, aes(x = Control_Butenyl_12.1, y = CW_Butenyl_12.1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = paired_means, aes(x = Elevation.x, y = CW_Butenyl_12.1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = paired_means, aes(x = Elevation.x, y = Control_Butenyl_12.1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = paired_means, aes(x = Control_Butenyl_12.1, y = CW_Butenyl_12.1, color = Elevation.x, label = Elevation.x)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "purple")

ggplot(data = dl, aes(x = biomass, y = Flavonol_16.1, color = Elevation, label = Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Add a single trend line
  theme_minimal() + 
  labs(x = "End Biomass (g)", y = "Flavonoid 16", color = "Elevation") + 
  scale_color_gradient(low = "orange", high = "purple")


# MSOO
ggplot(data = paired_means, aes(x = Control_MSOO_13.8, y = CW_MSOO_13.8)) + 
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(data = paired_means, aes(x = Elevation.x, y = CW_MSOO_13.8)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# Flaonol 17
ggplot(data = paired_means, aes(x = Control_Flavonol_17.5, y = CW_Flavonol_17.5)) + 
  geom_point() + 
  geom_smooth(method = "lm")

# interpreting the coeffs
# B0 = intrcept
# B1 = slope (effect) for x1 (Control I3M) 
# B2 = difference in intercept between X1 and x2
# B3 = difference in slope between x1 and x2 

anova_result1 <- anova(m1)
print(anova_result1)

# Compute EMMs
emm <- emmeans(m1, c("treatment", "Elevation"))
# Perform post-hoc tests for interactions
emm_pairs_specific <- pairs(emm, by = c("treatment", "Elevation"))
print(emm)
print(emm_pairs_specific)


# Q2 - growth-defense trade-offs & association with environment

#Fixed effects: growth x elevation
#Random effects: population, rack
#Statistical test - report slope of the regression line - r, r2, and p value
#If hypothesis is supported, positive slope across populations, populations that have higher biomass will have higher defense investment. The effect size of the slop should be bigger than 0.

colnames(dl)

# 3MSO
x3mso1 <- lm(x3MSO_5.2 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(x3mso1)

# OH-Alkenyl
ohalkenyl1 <- lm(OH.Alkenyl_6 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(ohalkenyl1)

# 4MSO
x4mso1 <- lm(x4MSO_7.1 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(x4mso1)

# Allyl
allyl1 <- lm(Allyl_7.4 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(allyl1)

# 5MSO
x5mso1 <- lm(x5MSO_10.2 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(x5mso1)

# Butenyl
butenyl1 <- lm(Butenyl_12.1 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(butenyl1)

# 3MT
x3mt1 <- lm(x3MT_13.6 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(x3mt1)

# MSOO
msoo1 <- lm(MSOO_13.8 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(msoo1)

# OH-I3M
ohi3m1 <- lm(OH.I3M_15.1 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(ohi3m1)

# 4MT
x4mt1 <- lm(x4MT._15.5 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(x4mt1)

# Flavonol 16
flavonol16 <- lm(Flavonol_16.1 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(flavonol16)

# I3M
i3m1 <- lm(I3M_16.7 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(i3m1)

# Flavonol 17
flavonol17 <- lm(Flavonol_17.5 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(flavonol17)

# Flavonol 18
flavonol18 <- lm(Flavonol_18.5 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(flavonol18)

# Indole 
indole1 <- lm(Indole_18.8 ~ biomass*Elevation, data = dl, na.action = na.exclude)
summary(indole1)

### no interactions q2

# 3MSO
x3mso2 <- lm(x3MSO_5.2 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(x3mso2)

# OH-Alkenyl
ohalkenyl2 <- lm(OH.Alkenyl_6 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(ohalkenyl2)

# 4MSO
x4mso2 <- lm(x4MSO_7.1 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(x4mso2)

# Allyl
allyl2 <- lm(Allyl_7.4 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(allyl2)

# 5MSO
x5mso2 <- lm(x5MSO_10.2 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(x5mso2)

# Butenyl
butenyl2 <- lm(Butenyl_12.1 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(butenyl2)

# 3MT
x3mt2 <- lm(x3MT_13.6 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(x3mt2)

# MSOO
msoo2 <- lm(MSOO_13.8 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(msoo2)

# OH-I3M
ohi3m2 <- lm(OH.I3M_15.1 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(ohi3m2)

# 4MT
x4mt2 <- lm(x4MT._15.5 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(x4mt2)

# Flavonol 16
flavonol16.1 <- lm(Flavonol_16.1 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(flavonol16.1)

# I3M
i3m2 <- lm(I3M_16.7 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(i3m2)

# Flavonol 17
flavonol17.2 <- lm(Flavonol_17.5 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(flavonol17.2)

# Flavonol 18
flavonol18.2 <- lm(Flavonol_18.5 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(flavonol18.2)

# Indole 
indole2 <- lm(Indole_18.8 ~ biomass + Elevation, data = dl, na.action = na.exclude)
summary(indole2)

### Visualize relationship

# Flavonol 16 by elevation
ggplot(data = dl, aes(x = biomass, y = Flavonol_16.1, color = Elevation, label = Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Add a single trend line
  theme_minimal() + 
  labs(x = "End Biomass (g)", y = "Flavonoid 16", color = "Elevation") + 
  scale_color_gradient(low = "orange", high = "purple") # Custom color palette from orange to purple

# Flavonol 16 by population
ggplot(data = dl, aes(x = biomass, y = Flavonol_16.1, color = Population, label = Population)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Add a single trend line
  theme_minimal() + 
  labs(x = "End Biomass (g)", y = "Flavonoid 16", color = "Elevation") 


ggplot(data = dl, aes(x = Elevation, y =Flavonol_16.1)) + 
  geom_point() + 
  geom_smooth(method = "lm")

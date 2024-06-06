### Danielle De La Pascua
### Ch2 Testing Intraspecifc Defense Trade-offs
### 1-28-24

### Linear Mixed Models

oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
options(oo)
install.packages("emmeans")

#libraries
library(lme4)
library(nlme)
library(Matrix)
library(emmeans)
library(dplyr)
library(tidyr)

#pull data
paired_means <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/paired_means.csv")

#question 1 - i want to know whether populations that have high underlying defenses invest more or less in induced defenses

#Response variable: each glucosinolate compounds (indoles most inducible, I3m & Indole)
#Fixed effects:  treatment (categorical), elevation (continuous)
#Random effects: population, rack
#ANOVA: differences between intercepts of treatment lines & slope of the line
#If hypothesis is true, both slopes will be positive, the induced line will have a steeper slope than the control line across elevations, and the induced line will be above the control line. At high elevations, they should be significantly different, but may not be significantly different at low elevations. Populations at high elevations should have higher overall levels of the compound than lower elevations

#Indoles

#OHI3M
colnames(paired_means)

OHI3M_q1 <- lmer(Control_I3M_16.7 ~ CW_I3M_16.7 + (1|Population), data = paired_means, na.action = na.exclude)
summary(OHI3M_q1)

anova_result1 <- anova(m1)
print(anova_result1)

# Compute EMMs
emm <- emmeans(m1, c("treatment", "Elevation"))
# Perform post-hoc tests for interactions
emm_pairs_specific <- pairs(emm, by = c("treatment", "Elevation"))
print(emm)
print(emm_pairs_specific)


#I3M

m2 <- lmer(I3M_16.7 ~ treatment*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(m2)

anova_result2 <- anova(m2)
print(anova_result2)

#why are these tables the same?

#Indole

m3 <- lmer(Indole_18.8 ~ treatment*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(m3)

#error message: boundary (singular) fit: see help('isSingular')

# Compute EMMs
emm <- emmeans(m1, c("treatment", "Elevation"))
# Perform post-hoc tests for interactions
emm_pairs_specific <- pairs(emm, by = c("treatment", "Elevation"))


print(emm_pairs_specific)
#nothing is showing up...


# Q2 - growth-defense trade-offs & association with environment

#Fixed effects: growth x elevation
#Random effects: population, rack
#Statistical test - report slope of the regression line - r, r2, and p value
#If hypothesis is supported, positive slope across populations, populations that have higher biomass will have higher defense investment. The effect size of the slop should be bigger than 0.

#OH-I3M
x1 <- lmer( OH.I3M_15.1 ~ biomass*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(x1)

anova_result1 <- anova(x1)
print(anova_result1)

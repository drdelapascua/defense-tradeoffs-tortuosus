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
data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")

#question 1 - i want to know whether populations that have high underlying defenses invest more or less in induced defenses

#Response variable: each glucosinolate compounds (indoles most inducible, I3m & Indole)
#Fixed effects:  treatment (categorical), elevation (continuous)
#Random effects: population, rack
#ANOVA: differences between intercepts of treatment lines & slope of the line
#If hypothesis is true, both slopes will be positive, the induced line will have a steeper slope than the control line across elevations, and the induced line will be above the control line. At high elevations, they should be significantly different, but may not be significantly different at low elevations. Populations at high elevations should have higher overall levels of the compound than lower elevations


OHI3Mm1 <- lmer( OH.I3M_15.1 ~ treatment*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(OHI3Mm1)

m2 <- lm(OH.I3M_15.1 ~ treatment*Elevation, data = data, na.action = na.exclude)
summary(m2)

anova_result <- anova(OHI3Mm1)
print(anova_result)

# Compute EMMs
emm <- emmeans(m2, c("treatment", "Elevation"))
# Perform post-hoc tests for interactions
emm_pairs_specific <- pairs(emm, by = c("treatment", "Elevation"))
#nothing is showing up...


print(emm_pairs_specific)

str(emm)

# Trying with other compounds

#X3MSOm1 <- lmer(X3MSO_5.2 ~ treatment*Elevation + (1|Population), data = data, na.action = na.exclude)
#summary(X3MSOm1)

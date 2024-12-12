### Ch2 Testing Intraspecifc Defense Trade-offs
### Danielle De La Pascua
### 1-28-24

### Linear Mixed Models ----

#libraries
library(lme4)
library(nlme)
library(Matrix)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)

#pull data
data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dw.csv")

head(data)

growth_data <- read.csv("./data/mf_means.csv") %>%
  select(Population, mf, treatment, biomass, logGSL, logindoles, logaliphatics) %>%
  filter(treatment == "C") %>%
  filter(!logindoles == "-Inf")

# Growth ~ defense models - q3 ----

#### total GSLs ----
hist(growth_data$logGSL)
growth_totalGSL_m1 <- lme(logGSL ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_totalGSL_m1)

# model diagnostics
plot(growth_totalGSL_m1) # scatering around 0-ish
qqnorm(residuals(growth_totalGSL_m1))
qqline(residuals(growth_totalGSL_m1))

#### total indoles ----
hist(growth_data$logindoles)
growth_indoles_m1 <- lme(logindoles ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_indoles_m1)

# model diagnostics
plot(growth_indoles_m1) # scatering around 0-ish
qqnorm(residuals(growth_indoles_m1))
qqline(residuals(growth_indoles_m1))

#### total aliphatics ----
hist(growth_data$logaliphatics)
growth_aliphatics_m1 <- lme(logaliphatics ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_aliphatics_m1)

# model diagnostics
plot(growth_indoles_m1) # scatering around 0-ish
qqnorm(residuals(growth_indoles_m1))
qqline(residuals(growth_indoles_m1))


# Climate & totals ----

mf_means_with_climate <- read.csv("./data/mf_means_with_clim.csv") %>%
  select(-X) %>%
  filter(treatment == "C") %>%
  filter(Population != "MtSH") %>%
  filter(Population != "YO10")

# simple model with pop as random effect

# Total GSLs
test_m1 <- lme(logGSL ~ PC1, random = ~1 | Population, data = mf_means_with_climate)
summary(test_m1)

ggplot(data = mf_means_with_climate, aes(x = PC1, y = logGSL, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Total GSLs & PC2
test_m6 <- lmer(logGSL ~ PC2 + (1|Population), data = mf_means_with_climate)
test_m6_nore <- lm(logGSL ~ PC2, data = mf_means_with_climate)
summary(test_m6)
summary(test_m6_nore)

ggplot(data = mf_means_with_climate, aes(x = PC2, y = logGSL, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Total indoles & PC 1
# drop inf
mf_means_with_climate_infdrop <- mf_means_with_climate %>%
  filter_all(all_vars(. != -Inf))
test_m2 <- lme(logindoles ~ PC1, random = ~1 | Population, data = mf_means_with_climate_infdrop)
summary(test_m2)


ggplot(data = mf_means_with_climate_infdrop, aes(x = PC1, y = logindoles, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Total indoles & PC 2
test_m3 <- lmer(logindoles ~ PC2 + (1|Population), data = mf_means_with_climate_infdrop)
test_m3_nore <- lm(logindoles ~ PC2, data = mf_means_with_climate_infdrop)
summary(test_m3)
summary(test_m3_nore)

ggplot(data = mf_means_with_climate_infdrop, aes(x = PC2, y = logindoles, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Total log aliphatics
test_m4 <- lme(logaliphatics ~ PC1, random = ~1 | Population, data = mf_means_with_climate)
summary(test_m4)

ggplot(data = mf_means_with_climate, aes(x = PC1, y = logaliphatics, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Total log aliphatics & PC2
test_m5 <- lmer(logaliphatics ~ PC2 + (1|Population), data = mf_means_with_climate)
test_m5_nore <- lm(logaliphatics ~ PC2, data = mf_means_with_climate)
summary(test_m5)
summary(test_m5_nore)

ggplot(data = mf_means_with_climate, aes(x = PC2, y = logaliphatics, color = elevation, label + elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# Question 1 ----

#i want to know whether populations that have high underlying defenses invest more or less in induced defenses

#Response variable: each glucosinolate compounds (indoles most inducible, I3m & Indole)
#Fixed effects:  treatment (categorical), elevation (continuous)
#Random effects: population, rack
#ANOVA: differences between intercepts of treatment lines & slope of the line
#If hypothesis is true, both slopes will be positive, the induced line will have a steeper slope than the control line across elevations, and the induced line will be above the control line. At high elevations, they should be significantly different, but may not be significantly different at low elevations. Populations at high elevations should have higher overall levels of the compound than lower elevations

# > Total GSLs ----


# > Indoles ----

#total idoles



#OHI3M
m1 <- lmer( OH.I3M_15.1 ~ treatment*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(m1)

anova_result1 <- anova(m1)
print(anova_result1)


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

emtrends <- emtrends(TotalGSL_q1, ~1, var = "totalGSL_C")
emtrends # slope below 1, not different than 0

#plot to see the relationship
ggplot(data = paired_means, aes(x = totalGSL_C, y = totalGSL_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# > Total Aliphatic ----
# assess trade-off
TotalAliphatic_q1 <- lmer(totalaliphatic_CW ~ totalaliphatic_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(TotalAliphatic_q1) 
anova(TotalAliphatic_q1) 

# Compute EMMs
emtrends <- emtrends(TotalAliphatic_q1, ~1, var = "totalaliphatic_C")
emtrends # slope below 1, not different than 0

#plot to see the relationship
ggplot(data = paired_means, aes(x = totalaliphatic_C, y = totalaliphatic_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# > Total Indoles ----
# assess trade-off
TotalIndole_q1 <- lmer(totalindole_CW ~ totalindole_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(TotalIndole_q1) 
anova(TotalIndole_q1) 

# Compute EMMs
emtrends <- emtrends(TotalIndole_q1, ~1, var = "totalindole_C")
emtrends # indoles have a slope below 1, not different than 0

#plot to see the relationship
ggplot(data = paired_means, aes(x = totalindole_C, y = totalindole_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# > individual compounds ----

# Indole
# assess trade-off
Indole_q1 <- lmer(Indole_CW ~ Indole_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(Indole_q1) #interaction is not significant, also not significant if by population
anova(Indole_q1) # not sig different across pops

# Compute EMMs
emtrends <- emtrends(Indole_q1, ~1, var = "Indole_C")
emtrends # lower than 1 overlaps 0

#plot to see the relationship
ggplot(data = paired_means, aes(x = Indole_C, y = Indole_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

### OH-I3M
#build model

OHI3m_q1 <- lmer(OHI3M_CW ~ OHI3M_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(OHI3m_q1) 

# Compute EMMs - 
emtrends <- emtrends(OHI3m_q1, ~1, var = "OHI3M_C")
emtrends # overlaps 1 and 0, large variation

#plot to see the relationship
ggplot(data = paired_means, aes(x = OHI3M_C, y = OHI3M_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")


### I3M
#build model
colnames(paired_means)
paired_means <- as.data.frame(paired_means)
I3M_q1 <- lmer(I3M_CW ~ I3M_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(I3M_q1)

# Compute EMMs
emtrends <- emtrends(I3M_q1, ~1, var = "I3M_C")
emtrends # lower than 1, overlaps 0

#plot to see the relationship
ggplot(data = paired_means, aes(x = I3M_C, y = I3M_CW, color = Elevation, label + Elevation)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")


# Question 2 ----

# growth-defense trade-offs & association with environment

#Fixed effects: growth x elevation
#Random effects: population, rack
#Statistical test - report slope of the regression line - r, r2, and p value
#If hypothesis is supported, positive slope across populations, populations that have higher biomass will have higher defense investment. The effect size of the slop should be bigger than 0.

#OH-I3M

summary(means)

# total GSLs
ggplot(data = means, aes(x = biomass, y = totalGSL, color = Elevation, label + Elevation)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# indole GSLs
ggplot(data = means, aes(x = biomass, y = totalindole, color = Elevation, label + Elevation)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

# aliphatic GSLs
ggplot(data = means, aes(x = biomass, y = totalaliphatic, color = Elevation, label + Elevation)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")

x1 <- lm(totalGSL ~ biomass , data = means, na.action = na.exclude)
summary(x1)

x2  <- lm(totalindole ~ biomass , data = means, na.action = na.exclude)
summary(x2)

x3 <- lm(totalaliphatic ~ biomass , data = means, na.action = na.exclude)
summary(x3)

anova_result1 <- anova(x1)
print(anova_result1)

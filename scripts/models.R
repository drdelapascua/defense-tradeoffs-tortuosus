### Danielle De La Pascua
### Ch2 Testing Intraspecifc Defense Trade-offs
### 1-28-24

### Linear Mixed Models ----

#libraries
library(lme4)
library(nlme)
library(Matrix)
library(emmeans)
library(dplyr)
library(tidyr)

# > pull data ----
paired_means <- read.csv("./data/paired_means.csv")
controls <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C")
head(controls)
dw <- read.csv("./data/dw.csv")

colnames(paired_means)
paired_means <- as.data.frame(paired_means)
# Question 1 ----

# i want to know whether populations that have high underlying defenses invest more or less in induced defenses

#Response variable: each glucosinolate compounds (indoles most inducible, I3m & Indole)
#Fixed effects:  treatment (categorical), elevation (continuous)
#Random effects: population, rack
#ANOVA: differences between intercepts of treatment lines & slope of the line
#If hypothesis is true, both slopes will be positive, the induced line will have a steeper slope than the control line across elevations, and the induced line will be above the control line. At high elevations, they should be significantly different, but may not be significantly different at low elevations. Populations at high elevations should have higher overall levels of the compound than lower elevations

# > Total GSLs ----
# assess trade-off
TotalGSL_q1 <- lmer(totalGSL_CW ~ totalGSL_C + (1|Population), data = paired_means, na.action = na.exclude)
summary(TotalGSL_q1) 
anova(TotalGSL_q1) 

# Compute EMMs
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

# OH-I3M
x1 <- lmer( OH.I3M_15.1 ~ biomass*Elevation + (1|Population), data = data, na.action = na.exclude)
summary(x1)

anova_result1 <- anova(x1)
print(anova_result1)

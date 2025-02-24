### Ch2 Testing Intraspecifc Defense Trade-offs
### Danielle De La Pascua
### 1-28-24

### Linear Mixed Models ----

#libraries
library(lme4)
library(nlme)
library(lmerTest)
library(Matrix)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

#pull data
data <- read.csv("./data/dw.csv")

head(data)

growth_data <- read.csv("./data/mf_means.csv") %>%
  select(Population, mf, treatment, biomass, logGSL, logindoles, logaliphatics, logflavonoids, shannon_diversity) %>%
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
qqline(residuals(growth_totalGSL_m1))  #JG: check patterning around the qqline: it deviates towards high values, so it looks like you have extreme high values

#### total indoles ----
hist(growth_data$logindoles)
growth_indoles_m1 <- lme(logindoles ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_indoles_m1)

# model diagnostics
plot(growth_indoles_m1) # scatering around 0-ish
qqnorm(residuals(growth_indoles_m1))
qqline(residuals(growth_indoles_m1))  #JG: looks like you may have some extreme low values, and one high

#### total aliphatics ----
hist(growth_data$logaliphatics)
growth_aliphatics_m1 <- lme(logaliphatics ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_aliphatics_m1)

# model diagnostics
plot(growth_indoles_m1) # scatering around 0-ish
qqnorm(residuals(growth_indoles_m1))
qqline(residuals(growth_indoles_m1)) #JG: looks like you may have some extreme low values, and one high

#### total flavonoids ----
hist(growth_data$logflavonoids)
growth_flavonoid_m1 <- lme(logflavonoids ~ biomass, random = ~1 | Population, data = growth_data)
summary(growth_flavonoid_m1) # positive scaling, but insig

# model diagnostics
plot(growth_flavonoid_m1) # scatering around 0-ish
qqnorm(residuals(growth_flavonoid_m1))
qqline(residuals(growth_flavonoid_m1)) 

# Climate & totals ----

# laod data

# pop means
pop_means_with_climate <- read.csv("./data/pop_means_with_clim.csv") %>%
  select(-X) %>%
  filter(treatment == "C") %>% 
  filter(is.finite(pop_means_with_climate$GSL_logindoles))
  

str(pop_means_with_climate)

pop_means_induced <- read.csv("./data/pop_means_with_clim.csv") %>%
  select(-X) %>%
  filter(treatment == "CW")

# mf means 
mf_means_with_climate <- read.csv("./data/mf_means_with_clim.csv") %>%
  select(-X) %>%
  filter(treatment == "C") 

# > Models ----

# total gsl and pc1
hist(pop_means_with_climate$GSL_logGSL)
GSL_contemporary <- lm(GSL_logGSL ~ contemporary_PC1 + contemporary_PC2, data = pop_means_with_climate)
GSL_historic <- lm(GSL_logGSL ~ historic_PC1 + historic_PC2, data = pop_means_with_climate)

summary(GSL_historic) # PC1 & PC2 p > 0.05
summary(GSL_contemporary) # PC1 & PC2 p > 0.05

# log indoles and pcs
hist(pop_means_with_climate$GSL_logindoles)

pop_means_with_climate <- pop_means_with_climate 
  

str(pop_means_with_climate)

indole_contemporary <- lm(GSL_logindoles ~ contemporary_PC1 + contemporary_PC2, data = pop_means_with_climate)
indole_historic <- glm(GSL_logindoles ~ historic_PC1 + historic_PC2, data = pop_means_with_climate)

summary(indole_historic) # PC1 & PC2 p > 0.05
summary(indole_contemporary) # PC1 & PC2 p > 0.05

# aliphatic and pc1
str(pop_means_with_climate)

aliphatic_contemporary <- lm(GSL_logaliphatics ~ contemporary_PC1 + contemporary_PC2, data = pop_means_with_climate)
aliphatic_historic <- glm(GSL_logaliphatics ~ historic_PC1 + historic_PC2, data = pop_means_with_climate)

summary(aliphatic_historic) # PC1 & PC2 p > 0.05
summary(aliphatic_contemporary) # PC1 & PC2 p > 0.05



# flavonoid and pc1

# aliphatic and pc1
str(pop_means_with_climate)

flavonoid_contemporary <- lm(GSL_logflavonoids ~ contemporary_PC1 + contemporary_PC2, data = pop_means_with_climate)
flavonoid_historic <- lm(GSL_logflavonoids ~ historic_PC1 + historic_PC2, data = pop_means_with_climate)

summary(flavonoid_historic) # PC1 & PC2 p > 0.05
summary(flavonoid_contemporary) # PC1 & PC2 p > 0.05


# shannon and PCs
str(pop_means_with_climate)

shannon_contemporary <- lm(GSL_shannon_diversity ~ contemporary_PC1 + contemporary_PC2, data = pop_means_with_climate)
shannon_historic <- lm(GSL_shannon_diversity ~ historic_PC1 + historic_PC2, data = pop_means_with_climate)

summary(shannon_historic) # PC1 & PC2 p > 0.05
summary(shannon_contemporary) # PC1 & PC2 p > 0.05

#theyre different - summaries are just similar
pop_means_with_climate$GSL_logflavonoids
pop_means_with_climate$GSL_shannon_diversity

# > Visualize data ----

# total GSL pop means
str(pop_means_with_climate)

# Define a function to generate scatterplots with regression lines
plot_pca_vs_gsl <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "GSL_logGSL")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "Log GSL",
      title = paste("logGSL vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1 <- plot_pca_vs_gsl(pop_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2 <- plot_pca_vs_gsl(pop_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3 <- plot_pca_vs_gsl(pop_means_with_climate, "historic_PC1", "Historic PC1")
plot4 <- plot_pca_vs_gsl(pop_means_with_climate, "historic_PC2", "Historic PC2")

# total GSL mf means

# Define a function to generate scatterplots with regression lines
plot_pca_vs_gsl <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "logGSL")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "mf means Log GSL",
      title = paste("mf means logGSL vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1mf <- plot_pca_vs_gsl(mf_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2mf <- plot_pca_vs_gsl(mf_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3mf <- plot_pca_vs_gsl(mf_means_with_climate, "historic_PC1", "Historic PC1")
plot4mf <- plot_pca_vs_gsl(mf_means_with_climate, "historic_PC2", "Historic PC2")

# Display the plots
print(plot1mf)
print(plot3mf)
print(plot2mf)
print(plot4mf)

### aliphatics

# aliphatic pop means
str(pop_means_with_climate)

# Define a function to generate scatterplots with regression lines
plot_pca_vs_aliphatic <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "GSL_logaliphatics")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "Log aliphatic",
      title = paste("log aliphatic vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1_aliphatic <- plot_pca_vs_aliphatic(pop_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2_aliphatic <- plot_pca_vs_aliphatic(pop_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3_aliphatic <- plot_pca_vs_aliphatic(pop_means_with_climate, "historic_PC1", "Historic PC1")
plot4_aliphatic <- plot_pca_vs_aliphatic(pop_means_with_climate, "historic_PC2", "Historic PC2")

# print plots
plot1_alipihatic
plot3_aliphatic
plot2_aliphatic
plot4_aliphatic

# total GSL mf means

# MFs - Define a function to generate scatterplots with regression lines
plot_pca_vs_aliphatics_mf <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "logaliphatics")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "mf means Log GSL",
      title = paste("mf means log aliphatics vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1mf_aliphatics <- plot_pca_vs_aliphatics_mf(mf_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2mf_aliphatics <- plot_pca_vs_aliphatics_mf(mf_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3mf_aliphatics <- plot_pca_vs_aliphatics_mf(mf_means_with_climate, "historic_PC1", "Historic PC1")
plot4mf_aliphatics <- plot_pca_vs_aliphatics_mf(mf_means_with_climate, "historic_PC2", "Historic PC2")

# Display the plots
print(plot1mf_aliphatics)
print(plot3mf_aliphatics)
print(plot2mf_aliphatics)
print(plot4mf_aliphatics)


### indoles

# indoles pop means
str(pop_means_with_climate)

# Define a function to generate scatterplots with regression lines
plot_pca_vs_indoles <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "GSL_logindoles")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "Log indoles",
      title = paste("log indoles vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1_indoles <- plot_pca_vs_indoles(pop_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2_indoles <- plot_pca_vs_indoles(pop_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3_indoles <- plot_pca_vs_indoles(pop_means_with_climate, "historic_PC1", "Historic PC1")
plot4_indoles <- plot_pca_vs_indoles(pop_means_with_climate, "historic_PC2", "Historic PC2")

# print plots
plot1_indoles
plot3_indoles
plot2_indoles
plot4_indoles

# total GSL mf means

# Define a function to generate scatterplots with regression lines
plot_pca_vs_indoles_mf <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "logindoles")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "mf means log indoles",
      title = paste("mf means log indoles vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1mf_indoles <- plot_pca_vs_indoles_mf(mf_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2mf_indoles <- plot_pca_vs_indoles_mf(mf_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3mf_indoles <- plot_pca_vs_indoles_mf(mf_means_with_climate, "historic_PC1", "Historic PC1")
plot4mf_indoles <- plot_pca_vs_indoles_mf(mf_means_with_climate, "historic_PC2", "Historic PC2")

# Display the plots
print(plot1mf_indoles)
print(plot3mf_indoles)
print(plot2mf_indoles)
print(plot4mf_indoles)


### flavonoids 

# flavonoids pop means
str(pop_means_with_climate)

# Define a function to generate scatterplots with regression lines
plot_pca_vs_flavonoids <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "GSL_logflavonoids")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "Log flavonoids",
      title = paste("log flavonoids vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1_flavonoids <- plot_pca_vs_flavonoids(pop_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2_flavonoids <- plot_pca_vs_flavonoids(pop_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3_flavonoids <- plot_pca_vs_flavonoids(pop_means_with_climate, "historic_PC1", "Historic PC1")
plot4_flavonoids <- plot_pca_vs_flavonoids(pop_means_with_climate, "historic_PC2", "Historic PC2")

# print plots
plot1_flavonoids
plot3_flavonoids
plot2_flavonoids
plot4_flavonoids

# total GSL mf means

# Define a function to generate scatterplots with regression lines
plot_pca_vs_flavonoids_mf <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "logflavonoids")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "mf means log flavonoids",
      title = paste("mf means log flavonoids vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1mf_flavonoids <- plot_pca_vs_flavonoids_mf(mf_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2mf_flavonoids <- plot_pca_vs_flavonoids_mf(mf_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3mf_flavonoids <- plot_pca_vs_flavonoids_mf(mf_means_with_climate, "historic_PC1", "Historic PC1")
plot4mf_flavonoids <- plot_pca_vs_flavonoids_mf(mf_means_with_climate, "historic_PC2", "Historic PC2")

# Display the plots
print(plot1mf_flavonoids)
print(plot3mf_flavonoids)
print(plot2mf_flavonoids)
print(plot4mf_flavonoids)



### shannon diversity

# shannon pop means
str(pop_means_with_climate)

# Define a function to generate scatterplots with regression lines
plot_pca_vs_shannon_diversity <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "GSL_shannon_diversity")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "shannon_diversity",
      title = paste("shannon_diversity vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1_shannon_diversity <- plot_pca_vs_shannon_diversity(pop_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2_shannon_diversity <- plot_pca_vs_shannon_diversity(pop_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3_shannon_diversity <- plot_pca_vs_shannon_diversity(pop_means_with_climate, "historic_PC1", "Historic PC1")
plot4_shannon_diversity <- plot_pca_vs_shannon_diversity(pop_means_with_climate, "historic_PC2", "Historic PC2")

# print plots
plot1_shannon_diversity
plot3_shannon_diversity
plot2_shannon_diversity
plot4_shannon_diversity

# total GSL mf means

# Define a function to generate scatterplots with regression lines
plot_pca_vs_shannon_diversity_mf <- function(data, x_var, x_label) {
  ggplot(data, aes_string(x = x_var, y = "shannon_diversity")) +
    geom_point(size = 3, color = "blue", alpha = 0.7) +  # Scatter points
    labs(
      x = x_label,
      y = "mf mean shannon_diversity",
      title = paste("mf mean shannon_diversity vs", x_label)
    ) +
    theme_minimal(base_size = 14)
}

# Generate plots for each climate PC
plot1mf_shannon_diversity <- plot_pca_vs_shannon_diversity_mf(mf_means_with_climate, "contemporary_PC1", "Contemporary PC1")
plot2mf_shannon_diversity <- plot_pca_vs_shannon_diversity_mf(mf_means_with_climate, "contemporary_PC2", "Contemporary PC2")
plot3mf_shannon_diversity <- plot_pca_vs_shannon_diversity_mf(mf_means_with_climate, "historic_PC1", "Historic PC1")
plot4mf_shannon_diversity <- plot_pca_vs_shannon_diversity_mf(mf_means_with_climate, "historic_PC2", "Historic PC2")

# Display the plots
print(plot1mf_shannon_diversity)
print(plot3mf_shannon_diversity)
print(plot2mf_shannon_diversity)
print(plot4mf_shannon_diversity)

### Visualize pairwise Pst


# > Total GSLs ----

# see how climate PC scores relate
plot(pop_means_with_climate$contemporary_PC1 ~ pop_means_with_climate$historic_PC1) # extremely correlatred
plot(pop_means_with_climate$contemporary_PC2 ~ pop_means_with_climate$historic_PC2) # extremely correlatred


# plot it

# Generate predictions from the lme model

# contemporary
pop_means_with_climate$predicted_logGSL_contemporary <- predict(GSL_contemporary)  # Fixed effects only

pop_means_with_climate <- pop_means_with_climate %>%
  arrange(contemporary_PC1)

str(pop_means_with_climate)

# Generate new data for smooth predictions
new_PC1_contemporary <- expand.grid(
  contemporary_PC1 = seq(min(pop_means_with_climate$contemporary_PC1), max(mf_means_with_climate$contemporary_PC1), length.out = 100),
  Population = unique(pop_means_with_climate$Population) # Keep all Population levels!
)

str(new_PC1_contemporary)

# Predict using both fixed + random effects
new_PC1_contemporary$predicted_logGSL <- predict(GSL_contemporary, newdata = new_PC1_contemporary, re.form = NULL)  


# Merge elevation data so each Population has its elevation value
elevation_data <- mf_means_with_climate %>% 
  select(Population, elevation.x) %>% 
  distinct()  # Keep unique Population-elevation pairs

new_PC1 <- new_PC1 %>% left_join(elevation_data, by = "Population")

# Plot again with group-specific fitted lines
ggplot(data = mf_means_with_climate, aes(x = PC1, y = logGSL, color = elevation.x)) + 
  geom_point(size = 3) + 
  geom_line(data = new_PC1, aes(x = PC1, y = predicted_logGSL, group = Population, color = elevation.x), 
            linewidth = 1) +  
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")  # Adjust color scale


# > Total indoles ----

# Total indoles & PC 1
# drop inf
mf_means_with_climate_infdrop <- mf_means_with_climate %>%
  filter_all(all_vars(. != -Inf))

indole_PC1 <- lmer(logindoles ~ PC1 + (1|Population), data = mf_means_with_climate_infdrop)
summary(indole_PC1)

indole_PC1_2 <- lmer(logindoles ~ + (1|Population), data = mf_means_with_climate_infdrop)

anova(indole_PC1_2, indole_PC1, test = "Chi") #p=0.57

summary(indole_PC1_2)

# Total indoles & PC 2
indole_PC2 <- lmer(logindoles ~ PC2 + (1|Population), data = mf_means_with_climate_infdrop)
summary(indole_PC2)

indole_PC2_2 <- lmer(logindoles ~ (1|Population), data = mf_means_with_climate_infdrop)

anova(indole_PC2, indole_PC2_2, test = "Chi") #p=0.68

# > Total aliphatics ----
aliphatics_PC1 <- lmer(logaliphatics ~ PC1 + (1|Population), data = mf_means_with_climate)
summary(aliphatics_PC1)

aliphatics_PC1_2 <- lmer(logaliphatics ~ + (1|Population), data = mf_means_with_climate)

anova(aliphatics_PC1_2, aliphatics_PC1, test = "LRT") #p=0.05 - significant

# simple beter
summary(aliphatics_PC1_2)

# Generate new data for smooth predictions
new_PC1_aliphatics_PC1 <- expand.grid(
  PC1 = seq(min(mf_means_with_climate$PC1), max(mf_means_with_climate$PC1), length.out = 100),
  Population = unique(mf_means_with_climate$Population) # Keep all Population levels!
)


# Predict using both fixed + random effects
new_PC1_aliphatics_PC1$predicted_aliphatics_PC1 <- predict(aliphatics_PC1 , newdata = new_PC1_aliphatics_PC1, re.form = NULL)  

# Merge elevation data so each Population has its elevation value
elevation_data <- mf_means_with_climate %>% 
  select(Population, elevation.x) %>% 
  distinct()  # Keep unique Population-elevation pairs

new_PC1_aliphatics_PC1  <- new_PC1_aliphatics_PC1 %>% left_join(elevation_data, by = "Population")

# Plot again with group-specific fitted lines
ggplot(data = mf_means_with_climate, aes(x = PC1, y = logaliphatics, color = elevation.x)) + 
  geom_point(size = 3) + 
  geom_line(data = new_PC1, aes(x = PC1, y = predicted_aliphatics_PC1, group = Population, color = elevation.x), 
            linewidth = 1) +  
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")  # Adjust color scale


# Total log aliphatics & PC2
aliphatics_PC2 <- lmer(logaliphatics ~ PC2 + (1|Population), data = mf_means_with_climate)
summary(aliphatics_PC2)

aliphatics_PC2_2 <- lmer(logaliphatics ~ (1|Population), data = mf_means_with_climate)

anova(aliphatics_PC2, aliphatics_PC2_2, test = "Chi") #p=0.66


# > Total flavonoids ----

hist(mf_means_with_climate$logGSL)
hist(mf_means_with_climate$logindoles)
hist(mf_means_with_climate$logaliphatics)
hist(mf_means_with_climate$logflavonoids)
hist(mf_means_with_climate$shannon_diversity)

### PC1
flav_PC1 <- lmer(logflavonoids ~ PC1 + (1|Population), data = mf_means_with_climate)
summary(flav_PC1)

#JG need likelihood ratio test for significance
flav_PC1_2 <- lmer(logflavonoids ~ (1|Population), data = mf_means_with_climate)
anova(flav_PC1, flav_PC1_2, test = "Chi") #p=0.27


# Generate new data for smooth predictions
new_PC1_flavonoid <- expand.grid(
  PC1 = seq(min(mf_means_with_climate$PC1), max(mf_means_with_climate$PC1), length.out = 100),
  Population = unique(mf_means_with_climate$Population) # Keep all Population levels!
)


# Predict using both fixed + random effects
new_PC1_flavonoid$predicted_flavonoid_PC1 <- predict(flav_PC1 , newdata = new_PC1_flavonoid, re.form = NULL)  


# Merge elevation data so each Population has its elevation value
elevation_data <- mf_means_with_climate %>% 
  select(Population, elevation) %>% 
  distinct()  # Keep unique Population-elevation pairs



new_PC1_flavonoid  <- new_PC1_flavonoid %>% left_join(elevation_data, by = "Population")

# Plot again with group-specific fitted lines
ggplot(data = mf_means_with_climate, aes(x = PC1, y = logflavonoids, color = elevation)) + 
  geom_point(size = 3) + 
  geom_line(data = new_PC1_flavonoid, aes(x = PC1, y = predicted_flavonoid_PC1, group = Population, color = elevation), 
            linewidth = 1) +  
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")  # Adjust color scale


### PC2
flav_PC2 <- lmer(logflavonoids ~ PC2 + (1|Population), data = mf_means_with_climate)
summary(flav_PC2)

flav_PC2_2 <- lmer(logflavonoids ~  (1|Population), data = mf_means_with_climate)

anova(flav_PC2, flav_PC2_2, test = "Chi") #p =0.03

# Generate new data for smooth predictions
new_PC2_flavonoid <- expand.grid(
  PC2 = seq(min(mf_means_with_climate$PC2), max(mf_means_with_climate$PC2), length.out = 100),
  Population = unique(mf_means_with_climate$Population) # Keep all Population levels!
)


# Predict using both fixed + random effects
new_PC2_flavonoid$predicted_flavonoid_PC2 <- predict(flav_PC2 , newdata = new_PC2_flavonoid, re.form = NULL)  


# Merge elevation data so each Population has its elevation value
elevation_data <- mf_means_with_climate %>% 
  select(Population, elevation) %>% 
  distinct()  # Keep unique Population-elevation pairs



new_PC2_flavonoid  <- new_PC2_flavonoid %>% left_join(elevation_data, by = "Population")

# Plot again with group-specific fitted lines
ggplot(data = mf_means_with_climate, aes(x = PC2, y = logflavonoids, color = elevation)) + 
  geom_point(size = 3) + 
  geom_line(data = new_PC2_flavonoid, aes(x = PC2, y = predicted_flavonoid_PC2, group = Population, color = elevation), 
            linewidth = 1) +  
  theme_minimal() + 
  scale_color_gradient(low = "orange", high = "blue")  # Adjust color scale


# > total shannon ----

### PC1
shan_PC1 <- lmer(shannon_diversity ~ PC1 + (1|Population), data = mf_means_with_climate)
summary(shan_PC1)

shan_PC1_2 <- lmer(shannon_diversity ~ + (1|Population), data = mf_means_with_climate)

anova(shan_PC1, shan_PC1_2, test = "Chi") #p =0.42

### PC2
shan_PC2 <- lmer(shannon_diversity ~ PC2 + (1|Population), data = mf_means_with_climate)
summary(shan_PC2)

shan_PC2_2 <- lmer(shannon_diversity ~ + (1|Population), data = mf_means_with_climate)

anova(shan_PC2, shan_PC2_2, test = "Chi") #p =0.47

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

### Shannon Diversity ----

# shannon & climate - controls
str(mf_means_with_climate)

shannon_climate_controls_1 <- lme(shannon_diversity ~ PC1, random = ~1 | Population, data = mf_means_with_climate)
summary(shannon_climate_controls_1)

shannon_climate_controls_2 <- lme(shannon_diversity ~ PC2, random = ~1 | Population, data = mf_means_with_climate)
summary(shannon_climate_controls_2)

# shannon & climate - induced
shannon_climate_cw_1 <- lme(shannon_diversity ~ PC1, random = ~1 | Population, data = mf_means_induced)
summary(shannon_climate_cw_1)

shannon_climate_cw_2 <- lme(shannon_diversity ~ PC2, random = ~1 | Population, data = mf_means_induced)
summary(shannon_climate_cw_2)

# shannon & elevation - controls

shannon_elevation_controls <- lme(shannon_diversity ~ elevation.x, random = ~1 | Population, data = mf_means_with_climate)
summary(shannon_elevation_controls)

# shannon & elevation - induced
shannon_elevation_cw <- lme(shannon_diversity ~ elevation.x, random = ~1 | Population, data = mf_means_induced)
summary(shannon_elevation_cw)



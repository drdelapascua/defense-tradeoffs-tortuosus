### Ch2 Testing Intraspecifc Defense Trade-offs
### Danielle De La Pascua
### 1-28-24

### Data Visualization ----

### > libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)
library("FactoMineR")
library(ggcorrplot)
library('corrr')
library(factoextra)
library(vegan)
library(missForest)
library(viridis)

### > load data ----
pop_means_long <- read.csv("./data/pop_means_long.csv") %>%
  select(-compound.type) %>%
  select(-X)
head(pop_means_long)

### > Stacked bar plot ----

# build big barplot
ggplot(pop_means_long, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal()

#build barplot with just butenyl & allyl
ab_means <- pop_means_long %>%
  filter(compound %in% c("Allyl", "Butenyl"))

ggplot(ab_means, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal()

# build barplot with all other cmpds
small_cmpd_means <- pop_means_long %>%
  filter(!compound %in% c("Allyl", "Butenyl"))

ggplot(small_cmpd_means, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal()

# build a barplot with allyl, butenyl, and then all other compounds
head(pop_means_long)

other_means <- pop_means_long %>%
  filter(!compound %in% c("Allyl", "Butenyl")) %>% # filter out allyl and butenyl
  group_by(Population, treatment, Elevation) %>%
  summarise(value = sum(value), .groups = 'drop')

head(other_means)
other_means <- as.data.frame(other_means)
head(other_means)

#make df
head(ab_means)
head(other_means)

#merge ab and other means
abother_means <- bind_rows(other_means, ab_means)%>%
  mutate(compound = if_else(is.na(compound), "All other compounds", compound))

head(abother_means)

abother_means2 <- abother_means %>%
  mutate(Population = factor(Population, levels = unique(Population[order(Elevation)])))

ggplot(abother_means2, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

ggplot(abother_means2, aes(x = Population, y = value, fill = compound)) + 
  geom_bar(stat = "identity") +
 # facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# Make a plot with only indole compounds

# make df
indole_cmpd_means <- pop_means_long %>%
  filter(compound %in% c("OHI3M", "I3M", "Indole"))
indole_cmpd_means
#order by elevation
indole_cmpd_means <- indole_cmpd_means %>%
  mutate(Population = factor(Population, levels = unique(Population[order(Elevation)])))

ggplot(indole_cmpd_means, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

#take out I3M
indole_cmpd_means_noI3M <- pop_means_long %>%
  filter(compound %in% c("OHI3M", "Indole"))
indole_cmpd_means_noI3M

#order by elevation
indole_cmpd_means_noI3M <- indole_cmpd_means_noI3M %>%
  mutate(Population = factor(Population, levels = unique(Population[order(Elevation)])))

ggplot(indole_cmpd_means_noI3M, aes(x = treatment, y = value, fill = compound)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Population) +  # Create separate plots for each treatment
  labs(x = "Population", y = "Mean Value", fill = "Compound") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")


### > PCA for induced and non-induced

# create filter two variables so we have one induced and one control

inducedPCA  = data %>%
  filter(data$treatment == "CW")

controlPCA = data %>%
  filter(data$treatment == "C")

#define vars of interest
variables_of_interest <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "X5MSO_10.2", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")

# subset data
inducedPCA <- inducedPCA[, c("Population", variables_of_interest)]
controlPCA <- controlPCA[, c("Population", variables_of_interest)]

#replace NAs with 0s
inducedPCA[is.na(inducedPCA)] <- 0
controlPCA[is.na(controlPCA)] <- 0

# make PCA
induced_pca_result <- prcomp(inducedPCA[, variables_of_interest])
control_pca_result <- prcomp(controlPCA[, variables_of_interest])

summary(induced_pca_result)
fviz_pca_var(induced_pca_result, col.var = "black")

summary(control_pca_result)
fviz_pca_var(control_pca_result, col.var = "black")

#shape be treatment, color elevation - by next meeting

### > NMDS ----

# creat a dataframe with only control plants
controls = data %>% 
  filter(data$treatment == "C")

#impute missing data
imputed_data <- missForest(controls[, c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1","X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")])

# Get data is in a dissimilarity matrix format
dist_matrix <- vegdist(imputed_data[, c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1","X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")], na.rm = TRUE)
dist_matrix <- vegdist(imputed_data$ximp)

# run NMDS with grouping
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100, trace = FALSE,  distance = "bray", subset = controls$Population, na.rm = TRUE)

env_data <- controls$Elevation
envfit_result <- envfit(nmds_result, env_data, na.rm = TRUE)
print(env_data)

# plot results
data$Population <- factor(data$Population)
levels(data$Population)
plot(nmds_result, type = "n")  # creates an empty plot
points(nmds_result, col = data$Population)  # adds points with colors based on the grouping variable
legend("topright", legend = levels(data$Population), col = data$Population, pch = 1, title = "Population", cex = 0.34)


### > GSL compounds by population & treatment ----

# 3MSO

X3MSO_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = X3MSO_5.2, fill = treatment)) + 
  geom_boxplot()

X3MSO_by_pop

# OH-Alkenyl

OH.Alkenyl_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = OH.Alkenyl_6, fill = treatment)) + 
  geom_boxplot()

OH.Alkenyl_by_pop

# 4MSO

X4MSO_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = X4MSO_7.1, fill = treatment)) + 
  geom_boxplot()

X4MSO_by_pop

# Allyl

Allyl_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Allyl_7.4, fill = treatment)) + 
  geom_boxplot()

Allyl_by_pop

# 5MSO

X5MSO_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = X5MSO_10.2, fill = treatment)) + 
  geom_boxplot()

X5MSO_by_pop

# Butenyl

Butenyl_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Butenyl_12.1, fill = treatment)) + 
  geom_boxplot()

Butenyl_by_pop 

# 3MT

X3MT_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = X3MT_13.6, fill = treatment)) + 
  geom_boxplot()

X3MT_by_pop

# MSOO

MSOO_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = MSOO_13.8, fill = treatment)) + 
  geom_boxplot()

MSOO_by_pop

# OH-I3M

OHI3M_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = OH.I3M_15.1, fill = treatment)) + 
  geom_boxplot()

OHI3M_by_pop

# 4MT

X4MT_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = X4MT._15.5, fill = treatment)) + 
  geom_boxplot()

X4MT_by_pop

# Flavonol 16

Flavonol16_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Flavonol_16.1, fill = treatment)) + 
  geom_boxplot()

Flavonol16_by_pop

# I3M

I3M_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = I3M_16.7, fill = treatment)) + 
  geom_boxplot()

I3M_by_pop

# Flavonol 17

Flavonol17_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Flavonol_17.5, fill = treatment)) + 
  geom_boxplot()

Flavonol17_by_pop

# Flavonol 18

Flavonol18_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Flavonol_18.5, fill = treatment)) + 
  geom_boxplot()

Flavonol18_by_pop

# indole 

Indole_by_pop <- ggplot(data = data, aes(x = reorder(Population, Elevation), y = Indole_18.8, fill = treatment)) + 
  geom_boxplot()

Indole_by_pop

# > each compound by biomass & population ----

# q - use control or cw or both?

# separate by population

# 3MSO

growth_defense_3MSO <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = X3MSO_5.2)) + 
  geom_point(aes(colour=Elevation)) + 
  stat_smooth(method=lm)

growth_defense_3MSO

# OH-Alkenyl

growth_defense_OHAlkenyl <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = OH.Alkenyl_6)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_OHAlkenyl

# 4MSO

growth_defense_4MSO <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = X4MSO_7.1)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_4MSO 

# Allyl

growth_defense_Allyl <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Allyl_7.4)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Allyl 

# 5MSO

growth_defense_5MSO <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = X5MSO_10.2)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_5MSO

# Butenyl

growth_defense_Butenyl <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Butenyl_12.1)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Butenyl

# 3MT

growth_defense_3MT <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = X3MT_13.6)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_3MT


# MSOO

growth_defense_MSOO <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = MSOO_13.8)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_MSOO

# OH-I3M

growth_defense_OH_I3M <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = OH.I3M_15.1)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_OH_I3M

# 4MT

growth_defense_4MT <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = X4MT._15.5)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_4MT

# Flavonol 16

growth_defense_Flavonol_16 <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Flavonol_16.1)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Flavonol_16

# I3M

growth_defense_Flavonol_I3M <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = I3M_16.7)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Flavonol_I3M

# Flavonol 17

growth_defense_Flavonol_17 <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Flavonol_17.5)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Flavonol_17

# Flavonol 18

growth_defense_Flavonol_18 <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Flavonol_18.5)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_Flavonol_18

# indole 

growth_defense_indole <- ggplot(subset(data, treatment %in% "C"), aes(x = biomass, y = Indole_18.8)) + 
  geom_point(aes(colour=Elevation)) +
  stat_smooth(method=lm)

growth_defense_indole

# > Multivariate analysis - PCA of compounds ----



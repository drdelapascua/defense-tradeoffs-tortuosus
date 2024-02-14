### Ch2 Testing Intraspecifc Defense Trade-offs
### Danielle De La Pascua
### 1-28-24

### Data Visualization ----

### > libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)

### > load data ----
data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")

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



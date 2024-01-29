### Ch2 Testing Intraspecifc Defense Trade-offs
### Danielle De La Pascua
### 1-28-24

### Data Visualization ----

### > libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)

### > load data ----
data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")

# filter so only the racks with D1-D6 are included
to_keep <- c("D1", "D2", "D3", "D4", "D5", "D6")
d = data %>%
  filter(data$Rack.. %in% to_keep)

### > GSL compounds by population & treatment

# 3MSO

X3MSO_by_pop <- ggplot(data = d, aes(x = Population, y = X3MSO_5.2, fill = treatment)) + 
  geom_boxplot()

X3MSO_by_pop

# OH-Alkenyl

OH.Alkenyl_by_pop <- ggplot(data = d, aes(x = Population, y = OH.Alkenyl_6, fill = treatment)) + 
  geom_boxplot()

OH.Alkenyl_by_pop

# 4MSO

X4MSO_by_pop <- ggplot(data = d, aes(x = Population, y = X4MSO_7.1, fill = treatment)) + 
  geom_boxplot()

X4MSO_by_pop

# Allyl

Allyl_by_pop <- ggplot(data = d, aes(x = Population, y = Allyl_7.4, fill = treatment)) + 
  geom_boxplot()

Allyl_by_pop

# 5MSO

X5MSO_by_pop <- ggplot(data = d, aes(x = Population, y = X5MSO_10.2, fill = treatment)) + 
  geom_boxplot()

X5MSO_by_pop

# Butenyl

Butenyl_by_pop <- ggplot(data = d, aes(x = Population, y = Butenyl_12.1, fill = treatment)) + 
  geom_boxplot()

Butenyl_by_pop

# 3MT

X3MT_by_pop <- ggplot(data = d, aes(x = Population, y = X3MT_13.6, fill = treatment)) + 
  geom_boxplot()

X3MT_by_pop

# MSOO

MSOO_by_pop <- ggplot(data = d, aes(x = Population, y = MSOO_13.8, fill = treatment)) + 
  geom_boxplot()

MSOO_by_pop

# OH-I3M

OHI3M_by_pop <- ggplot(data = d, aes(x = Population, y = OH.I3M_15.1, fill = treatment)) + 
  geom_boxplot()

OHI3M_by_pop

# 4MT

X4MT_by_pop <- ggplot(data = d, aes(x = Population, y = X4MT._15.5, fill = treatment)) + 
  geom_boxplot()

X4MT_by_pop

# Flavonol 16

Flavonol16_by_pop <- ggplot(data = d, aes(x = Population, y = Flavonol_16.1, fill = treatment)) + 
  geom_boxplot()

Flavonol16_by_pop

# I3M

I3M_by_pop <- ggplot(data = d, aes(x = Population, y = I3M_16.7, fill = treatment)) + 
  geom_boxplot()

I3M_by_pop

# Flavonol 17

Flavonol17_by_pop <- ggplot(data = d, aes(x = Population, y = Flavonol_17.5, fill = treatment)) + 
  geom_boxplot()

Flavonol17_by_pop

# Flavonol 18

Flavonol18_by_pop <- ggplot(data = d, aes(x = Population, y = Flavonol_18.5, fill = treatment)) + 
  geom_boxplot()

Flavonol18_by_pop

# indole 

Indole_by_pop <- ggplot(data = d, aes(x = Population, y =Indole_18.8, fill = treatment)) + 
  geom_boxplot()

Indole_by_pop

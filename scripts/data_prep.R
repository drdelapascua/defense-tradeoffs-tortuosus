### Ch2 Testing Intraspecifc Defense Trade-offs 
### Danielle De La Pascua
### 1-28-24

### Data Preparation ----

### > libraries ----
library(tidyr)
library(dplyr)
library(readxl)

### > load data ----

# GSL data

GSL <- read_xlsx("~/GitHub/defense-tradeoffs-tortuosus/data/Tortuosus.xlsx")
GSL$Location <- as.integer(GSL$Location)

# biomass & location data

biomass <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/mass.csv")
colnames(biomass) <- c("Rack #", "Location", "Population", "mf", "treatment", "biomass")

# site elevation, locality, and seed year data 

pop_data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/elevation.csv")

# join mass data & location data with GSL data

data <- left_join(GSL, biomass)

data <- left_join(data, pop_data)

### > prep data ----

# filter so only the racks with D1-D6 are included
to_keep <- c("D1", "D2", "D3", "D4", "D5", "D6")
data = data %>%
  filter(data$`Rack #` %in% to_keep)

# filter dataset so only induced leaf is included
d_induced = data %>%
  filter(data$leaf_type == "induced")

# create new variable with all compounds added together
data$X3MSO_5.2 + data$OH.Alkenyl_6 + data$X4MSO_7.1 + data$Allyl_7.4 + data$X5MSO_10.2 + data$Butenyl_12.1 + data$X3MT_13.6 + data$MSOO_13.8 + data$OH.I3M_15.1 + data$X4MT._15.5 + data$Flavonol_16.1 + data$I3M_16.7 + data$Flavonol_17.5 + data$Flavonol_18.5 + data$Indole_18.8
head(data)

### > save big data table ----

write.csv(data, "~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")
write.csv(d_induced, "~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")

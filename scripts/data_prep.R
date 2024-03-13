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

# create new variable with all compounds added together
head(data)

columns_to_sum <- c("3MSO_5.2", "OH-Alkenyl_6", "4MSO_7.1", "Allyl_7.4", "5MSO_10.2", "Butenyl_12.1", "3MT_13.6", "MSOO_13.8", "OH I3M_15.1", "4MT _15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")
sums <- rowSums(data[, columns_to_sum], na.rm = TRUE)

data$totalGSL <- sums

head(data)

# filter dataset so only induced leaf is included
d_induced = data %>%
  filter(data$leaf_type == "induced")


### > save big data table ----

write.csv(data, "~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")
write.csv(d_induced, "~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")

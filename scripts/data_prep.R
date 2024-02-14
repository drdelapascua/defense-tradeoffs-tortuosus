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

### > save big data table ----

write.csv(data, "~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")
write.csv(d_induced, "~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")

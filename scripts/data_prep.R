### Ch2 Testing Intraspecifc Defense Trade-offs 
### Danielle De La Pascua
### 1-28-24

### Data Preparation ----

### > libraries ----
library(tidyr)
library(dplyr)

### > load data ----

# GSL data

GSL <- read_xlsx("~/GitHub/defense-tradeoffs-tortuosus/data/Tortuosus.xlsx")
GSL$Location <- as.integer(GSL$Location)

# biomass & location data

biomass <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/mass.csv")
colnames(biomass) <- c("Rack #", "Location", "Population", "mf", "treatment", "biomass")

# join mass data & location data with GSL data

data <- left_join(GSL, biomass)

### > save data table ----

write.csv(data, "~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")



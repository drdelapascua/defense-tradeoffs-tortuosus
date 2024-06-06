### Ch2 Testing Intraspecifc Defense Trade-offs 
### Danielle De La Pascua
### 1-28-24

### Data Preparation ----

### > libraries ----
library(tidyr)
library(dplyr)
library(readxl)
library(tidyverse)

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

### manipulate data so we have a df with mf means for CW and C treated plants
dl <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")
head(dl)
dl <- dl[11:35]

# aggregate across mfs and by population
d_induced <- as.data.frame(d_induced)
colnames(dl)

x3MSO_mf <- aggregate(x3MSO_5.2 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
OHAlkenyl_mf <- aggregate(OH.Alkenyl_6 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x4MSO_mf <- aggregate(x4MSO_7.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Allyl_mf <- aggregate(Allyl_7.4 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x5MSO_mf <- aggregate(x5MSO_10.2 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Butenyl_mf <- aggregate(Butenyl_12.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x3MT_mf <- aggregate(x3MT_13.6 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
MSOO_mf <- aggregate(MSOO_13.8 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
OHI3M_mf <- aggregate(OH.I3M_15.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol16_mf <- aggregate(Flavonol_16.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x4MT_mf <- aggregate(x4MT._15.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
I3M_mf <- aggregate(I3M_16.7 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol17_mf <- aggregate(Flavonol_17.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol18_mf <- aggregate(Flavonol_18.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Indole_mf <- aggregate(Indole_18.8 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))

# join aggregated mf dfs together 
mf_means <- left_join(x = x3MSO_mf, y = OHAlkenyl_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = x4MSO_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Allyl_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = x5MSO_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Butenyl_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = x3MT_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = MSOO_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = OHI3M_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Flavonol16_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = x4MT_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = I3M_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Flavonol17_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Flavonol18_mf, by = c("mf", "Population", "treatment"))
mf_means <- left_join(x = mf_means, y = Indole_mf, by = c("mf", "Population", "treatment"))

# separate df into two dfs by trt
control_df <- subset(d_induced, treatment == 'C')
cw_df <- subset(d_induced, treatment == 'CW')

# manually clean col names
write.csv(control_df, "~/GitHub/defense-tradeoffs-tortuosus/data/c_df.csv")
write.csv(cw_df, "~/GitHub/defense-tradeoffs-tortuosus/data/cw_df.csv")

c_df <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/c_df.csv")
cw_df <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/cw_df.csv")

head(c_df)


# merge dfs
paired_means <- left_join(x = c_df, y = cw_df, by = c("mf", "Population"))

# save df
write.csv(paired_means, "~/GitHub/defense-tradeoffs-tortuosus/data/paired_means.csv")

### Defense Trade-offs: Qst-Fst Analysis
### Danielle De La Pascua
### 9-3-2024

### libraries ----
library(tidyr)
library(dplyr)
library(tidyverse)

### pull data ----

# genetic data
fst <- read.csv("./data/populations.fst_summary_full_dist_matrix.csv") %>%
  select(-c("YO11", "WV", "LVTR", "LV3")) #sxclude pops not in my study (columns) %>%
  filter(!X %in% c("YO11", "WV", "LVTR", "LV3")) #exclude pops not in my study (rows)

# change row X to rownames
rownames(fst) <- fst$X

# get rid of extra col with IDs
fst <- fst %>%
  select(-X) #wasnt working as just one pipe for some reason? 

# make table
fst <- as.matrix(fst)
fst <- as.data.frame(as.table(fst))

head(fst)

# get rid of duplicates
fst <- fst %>%
  filter(`Population 1` != `Population 2`) %>%
  arrange(`Population 1`, `Population 2`)

#  rename columns
colnames(fst) <- c("Population 1", "Population 2", "Fst")

# trait data
dw <-  read.csv("./data/dw.csv") %>%
  select("Population", "mf", "rep", "Elevation", "treatment", "X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")

# filter so only the fst populations are here 
dw <- dw %>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SHA", "SQ1", "SQ3", "WL1", "WL2","WL3", "YO1", "YO10"))


### calculate Qst ----

# Calculate variance components



# Extract between-population and total variances

# Calculate Qst for each trait

# print Qst values
### Ch2 Testing Intraspecifc Defense Trade-offs 
### Danielle De La Pascua
### 7-12-24

### Data Preparation ----

### > libraries ----
library(tidyr)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggfortify)
library(ggrepel)
library(viridis)
library(ggbiplot)
library(export)

### > load data ----

# load GSL data
#GSL <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/gsl_data.csv") %>%
GSL <- read.csv("./data/gsl_data.csv") %>% #JRG: this makes your path in relation to the github project, which should allow it to transfer across computers better - go up one folder indicated by ".", two levels "..", so on 
  select(-(Run.Index:Plate.Position)) %>% # removes first few columns from HPLC output
  select(-starts_with("Junk")) %>% #remove any columns that start with Junk
  filter(leaf_type == "induced") %>% # only shows focal leaf in df - Danielle still needs to change this to focal & clamped in df. This also removes NAs because NAs were not assigned "clamped" or "induced" because the label was left blank by accident
  mutate(leaf_type = recode(leaf_type, induced = "focal")) %>% #specifies that the "induced" to "focal"
  mutate(across(c(4:18), ~replace_na(.x, 0))) # replaces NAs with 0s for GSL compound cols - compound was not present in the sample

summary(GSL)
head(GSL)
unique(GSL$leaf_type)

#Make GSL Location col an integer for merging purposes
GSL$Location <- as.integer(GSL$Location)
summary(GSL)

# load biomass & experiment location & ID data
#biomass <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/biomass.csv")
biomass <- read.csv("./data/biomass.csv") %>% #with relational path now
           select(Population = "Pop", treatment = "trt", biomass = "mass..g.", everything()) #this renames your columns and tells it to include the rest

summary(biomass)

# load df with site elevation, locality, and seed year data 
#pop_data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/elevation.csv")
pop_data <- read.csv("./data/elevation.csv") %>% #now with relational pathway
  mutate(Population = recode(Population, YOSE10 = "YO10"))
  

unique(pop_data$Population)

summary(pop_data)

### > join data together ----

# join mass data & location data with GSL data

#join gsl and mass/rack location data by rack & location
data <- left_join(biomass, GSL, by = c("Rack", "Location"))
summary(data)
dim(biomass)
dim(GSL)
dim(data)

### > clean up data ----

# Create combined keys
GSL$key <- paste(GSL$Rack, GSL$Location)
biomass$key <- paste(biomass$Rack, biomass$Location)

# Find non-matching combinations
non_matching <- GSL$key[!GSL$key %in% biomass$key] #samples in GSL data not in biomass data
non_matching2 <- biomass$key[!biomass$key %in% GSL$key] #samples in biomass data not in GSL data

non_matching #samples in GSL not biomass: D4-25 - REWEIGH, i still have these samples

non_matching2 #samples in biomass not GSL: D4-19, D4-26 - note these were samples that did not have labels on microcentrifuge tube when running HPLC


#drop non-matching combinations
data$key <- paste(data$Rack, data$Location)
data <- data %>%
  filter(key != 'D4 25') %>%
  filter(key != 'D4 26') %>%
  filter(key != 'D4 19')

summary(data) 

#join df with elevation & pop loc data
data <- left_join(data, pop_data, by = "Population") 
dim(data)

table(is.na(data)) # means that there is still an NA in the dataframe

# the NA is for biomass of population TFC, treatment C, rack D1, location 29, mf 1, rep 1
# Danielle checked original datasheet, cell was left blank 7/23/24 - kept tissue, re-weigh?

# for now filter data point out

data <- data %>%
  filter(key != 'D1 29') 

table(is.na(data)) # FALSE, no NAs

### > calculating GSL sums ----

# total GSLs
head(data)
#add all to a col
GSL_to_sum <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "I3M_16.7", "Indole_18.8")
GSL_sums <- rowSums(data[, GSL_to_sum], na.rm = TRUE)
data$totalGSL <- GSL_sums 
head(data)

hist(x = data$totalGSL) # pretty normal, a little skewed left

# total indolic
indole_to_sum <- c("OH.I3M_15.1", "I3M_16.7", "Indole_18.8")
indole_sums <- rowSums(data[, indole_to_sum], na.rm = TRUE)
data$totalindole <- indole_sums 
head(data)

hist(x = data$totalindole) # most have low amounts of indoles, few have high amount, 

# total aliphatic
ali_to_sum <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "X4MT._15.5")
ali_sums <- rowSums(data[, ali_to_sum], na.rm = TRUE)
data$totalaliphatic <- ali_sums # pretty normal, a little skewed left
head(data)

hist(x = data$totalaliphatic) # more normal


# total flavonoids
flav_to_sum <- c("Flavonol_16.1", "Flavonol_17.5", "Flavonol_18.5")
flav_sums <- rowSums(data[, flav_to_sum], na.rm = TRUE)
data$totalflavonoid <- flav_sums # pretty normal, a little skewed left

hist(x = data$totalflavonoid) # skewed left

# Assess distributions
hist(x = data$totalGSL) # pretty normal, a little skewed left
hist(x = data$totalindole) # most have low amounts of indoles, few have high amount, 
hist(x = data$totalaliphatic) # more normal

# log transform all

data$logGSL <- log(data$totalGSL)

data$logindoles <- log(data$totalindole)

data$logaliphatics <- log(data$totalaliphatic)

data$logflavonoids <- log(data$totalflavonoid)

# all way more normally distributed
hist(x = data$logGSL)
hist(x = data$logindoles)
hist(x = data$logaliphatics)
hist(x = data$logflavonoids)

### > Calculate Shannon diversity across compounds ----

# make a df for calculating shannon diversity

shannon_df <- data %>%
  select("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8") 

# Function to calculate Shannon Diversity Index
shannon_diversity <- function(profile) {
  # Calculate the proportion of each compound
  proportions <- profile / sum(profile)
  
  # Shannon diversity formula: -sum(p_i * log(p_i))
  H <- -sum(proportions * log(proportions), na.rm = TRUE)
  
  return(H)
}

# Apply the function to each plant (row)
diversity_scores <- apply(shannon_df, 1, shannon_diversity)

# Print the diversity scores for each plant
print(diversity_scores)

# add to the big table
data$shannon_diversity <- diversity_scores

hist(data$shannon_diversity) # pretty normal!

### > save big data table ----

write.csv(data, "./data/dw.csv")

### > create a df without extreme values of total gsls, aliphatics, and indoles ----

plot(data$totalGSL) # no extreme values
plot(data$totalindole) # one extreme value
plot(data$totalaliphatic) 

data_ex_rem <- data %>%
  filter(totalindole < 700)

plot(data_ex_rem$totalindole)


### > aggregating means ----

head(data)

mf_means <- data %>% 
  # Summarize by maternal family
  group_by(Population, treatment, mf) %>% 
  summarise(
    biomass = mean(biomass),
    X3MSO = mean(X3MSO_5.2),
    OHAlkenyl = mean(OH.Alkenyl_6),
    X4MSO = mean(X4MSO_7.1),
    Allyl = mean(Allyl_7.4),
    X5MSO = mean(X5MSO_10.2),
    Butenyl = mean(Butenyl_12.1),
    X3MT = mean(X3MT_13.6),
    MSOO = mean(MSOO_13.8),
    OHI3M = mean(OH.I3M_15.1),
    X4MT = mean(X4MT._15.5),
    Flavonol16 = mean(Flavonol_16.1),
    I3M = mean(I3M_16.7),
    Flavonol17 = mean(Flavonol_17.5),
    Flavonol18 = mean(Flavonol_18.5),
    Indole = mean(Indole_18.8),
    totalaliphatic = mean(totalaliphatic),
    totalindole = mean(totalindole),
    totalGSL = mean(totalGSL),
    totalflavonoid = mean(totalflavonoid),
    logGSL = mean(logGSL),
    logindoles = mean(logindoles),
    logaliphatics = mean(logaliphatics),
    logflavonoids = mean(logflavonoids),
    shannon_diversity = mean(shannon_diversity)
  )

# Use distinct to see if additional grouping variables are present
distinct_data <- data %>% select(Population, treatment, mf) %>% distinct()
print(nrow(distinct_data)) 
print(distinct_data)

# make pop means
pop_means <- mf_means %>% 
  # Summarize by population
  group_by(Population, treatment) %>% 
  summarise(
    GSL_biomass = mean(biomass),
    GSL_X3MSO = mean(X3MSO),
    GSL_OHAlkenyl = mean(OHAlkenyl),
    GSL_X4MSO = mean(X4MSO),
    GSL_Allyl = mean(Allyl),
    GSL_X5MSO = mean(X5MSO),
    GSL_Butenyl = mean(Butenyl),
    GSL_X3MT = mean(X3MT),
    GSL_MSOO = mean(MSOO),
    GSL_OHI3M = mean(OHI3M),
    GSL_X4MT = mean(X4MT),
    GSL_Flavonol16 = mean(Flavonol16),
    GSL_I3M = mean(I3M),
    GSL_Flavonol17 = mean(Flavonol17),
    GSL_Flavonol18 = mean(Flavonol18),
    GSL_Indole = mean(Indole),
    GSL_totalaliphatic = mean(totalaliphatic),
    GSL_totalindole = mean(totalindole),
    GSL_totalGSL = mean(totalGSL),
    GSL_totalflavonoid = mean(totalflavonoid),
    GSL_logGSL = mean(logGSL),
    GSL_logindoles = mean(logindoles),
    GSL_logaliphatics = mean(logaliphatics),
    GSL_logflavonoids = mean(logflavonoids),
    GSL_shannon_diversity = mean(shannon_diversity)
    )

head(pop_means)
dim(pop_means)

# Use distinct to see if additional grouping variables are present
distinct_data <- data %>% select(Population, treatment, mf) %>% distinct()
print(nrow(distinct_data)) 
dim(pop_means)

pop_means <- as.data.frame(pop_means) 

# save dfs
write.csv(mf_means, "./data/mf_means.csv")
write.csv(pop_means, "./data/pop_means.csv")

# make long version of means df
pop_means_long <- pop_means %>%
  pivot_longer(
    cols = starts_with("GSL_"),  # Select columns that start with "chemical_"
    names_to = "compound",
    values_to = "value") %>%
  separate(compound, into = c("compound type", "compound"), sep = "_")

head(pop_means_long)

pop_means_long <- as.data.frame(pop_means_long)

# save long version of data
write.csv(pop_means_long, "./data/pop_means_long.csv")

### make agreggates for mf means without extreme values ----

### > aggregating means ----

head(data_ex_rem)

mf_means_ex_rem <- data_ex_rem %>% 
  # Summarize by maternal family
  group_by(Population, treatment, mf) %>% 
  summarise(
    X3MSO = mean(X3MSO_5.2),
    OHAlkenyl = mean(OH.Alkenyl_6),
    X4MSO = mean(X4MSO_7.1),
    Allyl = mean(Allyl_7.4),
    X5MSO = mean(X5MSO_10.2),
    Butenyl = mean(Butenyl_12.1),
    X3MT = mean(X3MT_13.6),
    MSOO = mean(MSOO_13.8),
    OHI3M = mean(OH.I3M_15.1),
    X4MT = mean(X4MT._15.5),
    Flavonol16 = mean(Flavonol_16.1),
    I3M = mean(I3M_16.7),
    Flavonol17 = mean(Flavonol_17.5),
    Flavonol18 = mean(Flavonol_18.5),
    Indole = mean(Indole_18.8),
    totalaliphatic = mean(totalaliphatic),
    totalindole = mean(totalindole),
    totalGSL = mean(totalGSL),
    logGSL = mean(logGSL),
    logindoles = mean(logindoles),
    logaliphatics = mean(logaliphatics),
    biomass = mean(biomass),
    totalflavonoid = mean(totalflavonoid),
    logflavonoids = mean(logflavonoids),
    shannon_diversity = mean(shannon_diversity)
  )

# Use distinct to see if additional grouping variables are present
distinct_data <- data_ex_rem %>% select(Population, treatment, mf) %>% distinct()
print(nrow(distinct_data)) 
print(distinct_data)

# make pop means
pop_means_ex_rem <- mf_means_ex_rem %>% 
  # Summarize by population
  group_by(Population, treatment) %>% 
  summarise(
    GSL_X3MSO = mean(X3MSO),
    GSL_OHAlkenyl = mean(OHAlkenyl),
    GSL_X4MSO = mean(X4MSO),
    GSL_Allyl = mean(Allyl),
    GSL_X5MSO = mean(X5MSO),
    GSL_Butenyl = mean(Butenyl),
    GSL_X3MT = mean(X3MT),
    GSL_MSOO = mean(MSOO),
    GSL_OHI3M = mean(OHI3M),
    GSL_X4MT = mean(X4MT),
    GSL_Flavonol16 = mean(Flavonol16),
    GSL_I3M = mean(I3M),
    GSL_Flavonol17 = mean(Flavonol17),
    GSL_Flavonol18 = mean(Flavonol18),
    GSL_Indole = mean(Indole),
    GSL_totalaliphatic = mean(totalaliphatic),
    GSL_totalindole = mean(totalindole),
    GSL_totalGSL = mean(totalGSL),
    GSL_logindoles = mean(logindoles),
    GSL_logaliphatics = mean(logaliphatics),
    GSL_biomass = mean(biomass),
    GSL_totalflavonoid = mean(totalflavonoid),
    GSL_logflavonoids = mean(logflavonoids),
    GSL_shannon_diversity = mean(shannon_diversity)
  )

head(pop_means_ex_rem)
dim(pop_means_ex_rem)

# Use distinct to see if additional grouping variables are present
distinct_data <- data %>% select(Population, treatment, mf) %>% distinct()
print(nrow(distinct_data)) 
dim(pop_means)

pop_means <- as.data.frame(pop_means) 

# save dfs
write.csv(mf_means_ex_rem, "./data/mf_means_ex_rem.csv")
write.csv(pop_means_ex_rem, "./data/pop_means_ex_rem.csv")

### make df with no OH-I3M ----

data_rem_ohi3m <- data_ex_rem %>%
  select(!c("OH.I3M_15.1","totalGSL","totalaliphatic","logGSL", "logindoles", "totalindole", "logaliphatics"))
  
# total GSLs
head(data_rem_ohi3m)
#add all to a col
GSL_to_sum_rem_ohi3m <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "X4MT._15.5", "I3M_16.7", "Indole_18.8")
GSL_sums_rem_ohi3m <- rowSums(data_rem_ohi3m[, GSL_to_sum_rem_ohi3m ], na.rm = TRUE)
data_rem_ohi3m$totalGSL <- GSL_sums_rem_ohi3m 
head(data)

hist(x = data$totalGSL) # pretty normal, a little skewed left

# total indolic
indole_to_sum_rem_ohi3m <- c("I3M_16.7", "Indole_18.8")
indole_sums_rem_ohi3m <- rowSums(data_rem_ohi3m[, indole_to_sum_rem_ohi3m], na.rm = TRUE)
data_rem_ohi3m$totalindole <- indole_sums_rem_ohi3m 
head(data_rem_ohi3m)

hist(x = data$totalindole) # most have low ampounts of indoles, few have high amount, 

# total aliphatic
ali_to_sum_rem_ohi3m <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "X4MT._15.5")
ali_sums_rem_ohi3m <- rowSums(data_rem_ohi3m[, ali_to_sum_rem_ohi3m], na.rm = TRUE)
data_rem_ohi3m$totalaliphatic <- ali_sums_rem_ohi3m # pretty normal, a little skewed left
head(data_rem_ohi3m)

# total flavonoid
flav_to_sum_rem_ohi3m <- c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "X4MT._15.5")
flav_sums_rem_ohi3m <- rowSums(data_rem_ohi3m[, flav_to_sum_rem_ohi3m], na.rm = TRUE)
data_rem_ohi3m$totalflavonoid <- flav_sums_rem_ohi3m # pretty normal, a little skewed left


  
# log transform all

data_rem_ohi3m$logGSL <- log(data_rem_ohi3m$totalGSL)

data_rem_ohi3m$logindoles <- log(data_rem_ohi3m$totalindole)

data_rem_ohi3m$logaliphatics <- log(data_rem_ohi3m$totalaliphatic)

data_rem_ohi3m$logflavonoid <- log(data_rem_ohi3m$totalflavonoid)


hist(x = data_rem_ohi3m$logGSL)
hist(x = data_rem_ohi3m$logindoles)
hist(x = data_rem_ohi3m$logaliphatics)

### > aggregating means ----

mf_means_rem_ohi3m <- data_rem_ohi3m %>% 
  # Summarize by maternal family
  group_by(Population, treatment, mf) %>% 
  summarise(
    X3MSO = mean(X3MSO_5.2),
    OHAlkenyl = mean(OH.Alkenyl_6),
    X4MSO = mean(X4MSO_7.1),
    Allyl = mean(Allyl_7.4),
    X5MSO = mean(X5MSO_10.2),
    Butenyl = mean(Butenyl_12.1),
    X3MT = mean(X3MT_13.6),
    MSOO = mean(MSOO_13.8),
    X4MT = mean(X4MT._15.5),
    Flavonol16 = mean(Flavonol_16.1),
    I3M = mean(I3M_16.7),
    Flavonol17 = mean(Flavonol_17.5),
    Flavonol18 = mean(Flavonol_18.5),
    Indole = mean(Indole_18.8),
    totalaliphatic = mean(totalaliphatic),
    totalindole = mean(totalindole),
    totalGSL = mean(totalGSL),
    logGSL = mean(logGSL),
    logindoles = mean(logindoles),
    logaliphatics = mean(logaliphatics),
    biomass = mean(biomass), 
    totalflavonoid = mean(totalflavonoid),
    logflavonoids = mean(logflavonoid), 
    shannon_diversity = mean(shannon_diversity)
  )

# make pop means
pop_means_rem_ohi3m <- mf_means_rem_ohi3m %>% 
  # Summarize by population
  group_by(Population, treatment) %>% 
  summarise(
    GSL_X3MSO = mean(X3MSO),
    GSL_OHAlkenyl = mean(OHAlkenyl),
    GSL_X4MSO = mean(X4MSO),
    GSL_Allyl = mean(Allyl),
    GSL_X5MSO = mean(X5MSO),
    GSL_Butenyl = mean(Butenyl),
    GSL_X3MT = mean(X3MT),
    GSL_MSOO = mean(MSOO),
    GSL_X4MT = mean(X4MT),
    GSL_Flavonol16 = mean(Flavonol16),
    GSL_I3M = mean(I3M),
    GSL_Flavonol17 = mean(Flavonol17),
    GSL_Flavonol18 = mean(Flavonol18),
    GSL_Indole = mean(Indole),
    GSL_totalaliphatic = mean(totalaliphatic),
    GSL_totalindole = mean(totalindole),
    GSL_totalGSL = mean(totalGSL),
    GSL_logindoles = mean(logindoles),
    GSL_logaliphatics = mean(logaliphatics),
    GSL_biomass = mean(biomass)
  )

head(pop_means_rem_ohi3m)
dim(pop_means_rem_ohi3m)

# save dfs
write.csv(mf_means_rem_ohi3m, "./data/mf_means_rem_ohi3m.csv")
write.csv(pop_means_rem_ohi3m, "./data/pop_means_rem_ohi3m.csv")


### Climate data ----

all_locs = read.csv("./data/localities.csv") %>%
  mutate(id = if_else(id == "Ben Hur", "BH", id)) %>%
  mutate(id = if_else(id == "TM2 (was TM P)", "TM2", id))
#  select(-(collection_date:locality)) %>%
#  select(-pop_gen)

head(all_locs)

my_locs = read.csv("./data/elevation.csv") %>%
  select(id = "Population", latitude = "Lat", longitude = "Long", elevation = "Elevation")

my_locs$taxon_name <- "Streptanthus tortuosus"

head(my_locs)
summary(my_locs)

climate = read_csv("./data/flintbcm_climate_tall_herbarium.csv") %>% 
  filter(clim_year > 1950, clim_year < 2000) %>% 
  mutate(pck = abs(pck)) %>% 
  group_by(id) %>% 
  dplyr::summarize(cwd = sum(cwd), ppt_mm = sum(ppt_mm), pck = sum(pck), snw = sum(snw), tmin = mean(tmin), tmax = mean(tmax)) %>% 
  mutate(id = if_else(id == "Ben Hur", "BH", id)) %>%
  mutate(id = if_else(id == "TM2 (was TM P)", "TM2", id)) %>%
  left_join(., all_locs) %>% 
  filter(!is.na(cwd), taxon_name %in% c("Streptanthus tortuosus", "Streptanthus tortuosus var. tortuosus"))

my_climate = read_csv("./data/flintbcm_climate_tall_herbarium.csv") %>% 
  filter(clim_year > 1950, clim_year < 2000) %>% 
  mutate(pck = abs(pck)) %>% 
  group_by(id) %>% 
  dplyr::summarize(cwd = sum(cwd), ppt_mm = sum(ppt_mm), pck = sum(pck), snw = sum(snw), tmin = mean(tmin), tmax = mean(tmax)) %>% 
  mutate(id = if_else(id == "Ben Hur", "BH", id)) %>%
  mutate(id = if_else(id == "TM2 (was TM P)", "TM2", id)) %>%
  mutate(id = if_else(id == "UCD161863", "MtSH", id)) %>%
  mutate(id = if_else(id == "UCD164328", "CALO", id)) %>%
  mutate(id = if_else(id == "RSA795268", "LC", id)) %>%
  mutate(id = if_else(id == "RSA0085001", "TFC", id)) %>%
  mutate(id = if_else(id == "CHSC98313", "RB", id)) %>%
  full_join(., my_locs) %>% 
  filter(!is.na(cwd), taxon_name %in% c("Streptanthus tortuosus", "Streptanthus tortuosus var. tortuosus"))

# get rid of the 2 pops that were excluded form experimment (MtSH and ) 

climate_for_pc = my_climate %>% 
  select(cwd, pck, ppt_mm, snw, tmin, tmax)

summary(climate_for_pc)

pc = prcomp(climate_for_pc, scale = TRUE)

pc_data = data.frame(pc$x)

locs_pc = cbind(my_climate, pc_data)

# Add identifying information to PCA scores
pc_data_with_id <- cbind(ID = locs_pc$id, pc_data[, c("PC1", "PC2")])

loadings = data.frame(varnames=rownames(pc$rotation), pc$rotation)

# Extract PCA scores for the first two principal components
pca_scores <- as.data.frame(pc$x[, 1:2])

# Add sample names (row names) to the PCA scores
pca_scores$names <- my_climate$id

# Extract PCA loadings (variable contributions)
pca_loadings <- as.data.frame(pc$rotation[, 1:2])
colnames(pca_loadings) <- c("PC1", "PC2")
pca_loadings$variable <- rownames(pca_loadings)  # Add variable names

scaling_factor <- 3  # Adjust this value as needed
pca_loadings <- pca_loadings %>%
  mutate(PC1 = PC1 * scaling_factor, PC2 = PC2 * scaling_factor)

# Ensure `elevation` is part of the PCA scores data
pca_scores <- cbind(pca_scores, elevation = my_climate$elevation)

# Plot PCA
ggplot() +
  # Plot sample points with elevation as color
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = elevation), size = 3) +
  # Add labels for sample points
  geom_text(data = pca_scores, aes(x = PC1, y = PC2, label = names), vjust = -0.5, size = 3) +
  # Add loading vectors as arrows
  geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  # Add labels for loadings
  geom_text(data = pca_loadings, aes(x = PC1, y = PC2, label = variable), 
            color = "red", vjust = -0.5, size = 3) +
  # Add axis labels and title
  labs(title = "PCA Plot with Loadings and Elevation", x = "PC1 (79.46%)", y = "PC2 (14.91%)", color = "Elevation") +
  scale_color_gradient(low="orange", high="blue") +
  theme_minimal()

# Add PC1 and PC2 values to the big dataframe
pca_scores_df <- pca_scores %>%
  rename(Population = names)

my_climate_df <- my_climate %>%
  rename(Population = id)

my_climate_df2 <- my_climate_df %>%
  left_join(pca_scores_df, by = "Population") %>%
  select(-taxon_name)

# save big climate df 
write.csv(my_climate_df2, file = "./data/climate_PC_data.csv")

# merge climate data into the big data frame
data_with_clim <- data %>% 
  left_join(my_climate_df2, by = "Population") 

# merge climate data with the mf means
mf_means_with_clim <- mf_means %>% 
  left_join(my_climate_df2, by = "Population") 

# merge climate data with the pop means
pop_means_with_clim <- pop_means %>%
  left_join(my_climate_df2, by = "Population") 
  
# save big df
write.csv(data_with_clim, file = "./data/dw_with_clim.csv")

# save mf means
write.csv(mf_means_with_clim, file = "./data/mf_means_with_clim.csv")

# save pops means
write.csv(pop_means_with_clim, file = "./data/pop_means_with_clim.csv")

controls <- data %>%
  filter(treatment == "C")

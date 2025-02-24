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
library("factoextra")


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

#all_locs = read.csv("./data/localities.csv") %>%
#  mutate(id = if_else(id == "Ben Hur", "BH", id)) %>%
#  mutate(id = if_else(id == "TM2 (was TM P)", "TM2", id))
#  select(-(collection_date:locality)) %>%
#  select(-pop_gen)

#head(all_locs)

my_locs = read.csv("./data/elevation.csv") %>%
  select(pop = "Population", latitude = "Lat", longitude = "Long", elevation = "Elevation") %>%
  mutate(pop = gsub("MtSH", "MSH", pop))

str(my_locs)
my_locs$taxon_name <- "Streptanthus tortuosus"

# make filter for my pops
pops_to_keep <- my_locs$pop

# load data & filter to only my pops 

# contemporary
my_climate_contemporary = read_csv("./data/climate/Dimensions_All_1895-2022.csv") %>% 
  select(pop, year, month, cwd, ppt, str, pck, tmn, tmx) %>%
  filter(year > 1989, year < 2021) %>% 
  mutate(pck = abs(pck)) %>% 
  filter(pop %in% pops_to_keep) %>%
  group_by(pop) %>%
  left_join(my_locs, by = "pop")

# historic
my_climate_historic = read_csv("./data/climate/Dimensions_All_1895-2022.csv") %>% 
  select(pop, year, month, cwd, ppt, str, pck, tmn, tmx) %>%
  filter(year < 1990) %>% 
  mutate(pck = abs(pck)) %>% 
  filter(pop %in% pops_to_keep) %>%
  group_by(pop) %>%
  left_join(my_locs, by = "pop")

str(my_climate_historic)

### load in langes crossing data

# Read each file into a separate dataframe
cwd_df <- read_csv("./data/climate/Danielle_cwd_1895-2022.csv")
pck_df <- read_csv("./data/climate/Danielle_pck_1895-2022.csv")
ppt_df <- read_csv("./data/climate/Danielle_ppt_1895-2022.csv")
str_df <- read_csv("./data/climate/Danielle_str_1895-2022.csv")
tmx_df <- read_csv("./data/climate/Danielle_tmx_1895-2022.csv")
tmn_df <- read_csv("./data/climate/Danielle_tmn_1895-2022.csv")

# Ensure column names are consistent
dfs <- list(cwd_df, pck_df, ppt_df, str_df, tmx_df, tmn_df)
dfs <- lapply(dfs, function(df) {
  colnames(df) <- tolower(colnames(df))  # Standardize column names
  return(df)
})

# Merge datasets using left joins by 'year' and 'month'
LC_climate_contemporary <- reduce(dfs, left_join, by = c("year", "month")) %>%
  filter(year > 1989, year < 2021)
  
# View the merged dataframe
print(head(LC_climate_contemporary))

# Add a new column "pop" with value "LC"
LC_climate_contemporary_2 <- LC_climate_contemporary %>%
  mutate(pop = "LC") %>%
  left_join(my_locs, by = "pop")

# View the updated dataframe
print(head(LC_climate_contemporary_2))

# Save merged climate data
write_csv(LC_climate_contemporary_2, "./data/climate/LC_merged_climate_data.csv")

# add to the big climate df

all_climate_contemporary <- bind_rows(my_climate_contemporary, LC_climate_data_2)

### historics LC data
LC_climate_historic <- reduce(dfs, left_join, by = c("year", "month")) %>%
  filter(year < 1990)

# View the merged dataframe
print(head(LC_climate_historic))

# Add a new column "pop" with value "LC"
LC_climate_historic_2 <- LC_climate_historic %>%
  mutate(pop = "LC") %>%
  left_join(my_locs, by = "pop")

# View the updated dataframe
print(head(LC_climate_historic_2))

# Save merged climate data
write_csv(LC_climate_historic_2, "./data/climate/LC_merged_climate_data_historic.csv")

# add to the big climate df

all_climate_historic <- bind_rows(my_climate_historic, LC_climate_historic_2)

print(all_climate_historic)
### > growing season ----

### Contemporary

# Define site groups
high_elev_sites <- c("LV3", "SQ3", "WL2", "YO1", "SQ1", "WL3", "SQ2", "YOSE8", "CP2", "YOSE10", "LV1", "LV2")
low_elev_sites <- c("WL1", "BH", "TM2", "KC2", "IH", "SC", "RB", "MSH", "CALO", "TFC", "DPR", "LC")

# Function to calculate growing season length for high-elevation sites
calc_growing_season_high <- function(df) {
  df %>%
    arrange(year, month) %>%
    group_by(pop, year) %>%
    summarize(
      start_month = first(month[pck == 0], default = NA),  # First snow-free month
      end_month = first(month[pck > 0 & !is.na(pck)], default = NA),  # First snow-covered month
      growing_season = case_when(
        is.na(start_month) | is.na(end_month) ~ NA_real_,  # Exclude incomplete cases
        end_month < start_month ~ (12 - start_month) + end_month,  # Carry over to next year
        TRUE ~ end_month - start_month  # Normal case
      )
    ) %>%
    ungroup()
}
# Function to calculate growing season length for low-elevation sites
calc_growing_season_low <- function(df) {
  df %>%
    arrange(year, month) %>%
    group_by(pop, year) %>%
    summarize(
      start_month = first(month[ppt > 25], default = NA),  # First month with ppt > 25 mm
      end_month = last(month[ppt > 0], default = NA),  # Last month with any precipitation
      growing_season = ifelse(!is.na(start_month) & !is.na(end_month), end_month - start_month, NA)
    ) %>%
    ungroup()
}


# Compute growing season length for high-elevation sites
high_elev_contemp_gs <- all_climate_contemporary %>%
  filter(pop %in% high_elev_sites) %>%
  calc_growing_season_high()

# Compute growing season length for low-elevation sites
low_elev_contemp_gs <- all_climate_contemporary %>%
  filter(pop %in% low_elev_sites) %>%
  calc_growing_season_low()

# Combine and calculate mean growing season length (1990-2020)
growing_season_summary_contemp <- bind_rows(high_elev_contemp_gs, low_elev_contemp_gs) %>%
  group_by(pop) %>%
  summarize(avg_growing_season = mean(growing_season, na.rm = TRUE))

# Print summary
print(growing_season_summary_contemp)

### Historic

# Compute growing season length for high-elevation sites
high_elev_hist_gs <- all_climate_historic %>%
  filter(pop %in% high_elev_sites) %>%
  calc_growing_season_high()

print(all_climate_historic)

# Compute growing season length for low-elevation sites
low_elev_hist_gs <- all_climate_historic %>%
  filter(pop %in% low_elev_sites) %>%
  calc_growing_season_low()

# Combine and calculate mean growing season length (1990-2020)
growing_season_summary_historic <- bind_rows(high_elev_hist_gs, low_elev_hist_gs) %>%
  group_by(pop) %>%
  summarize(avg_growing_season = mean(growing_season, na.rm = TRUE))

# Print summary
print(growing_season_summary_historic)

### make df for pc data

# contemp
climate_for_pc_contemporary = all_climate_contemporary %>% 
  select(pop, year, month, cwd, pck, ppt, str, tmn, tmx, elevation) %>%
  mutate(water_year = ifelse(month >= 9, year + 1, year)) %>% # Shift year starting from September
  dplyr::summarize(cwd = sum(cwd), ppt = sum(ppt), pck = sum(pck), tmin = mean(tmn), tmax = mean(tmx), elevation = mean(elevation), str = sum(str)) %>% 
  left_join(growing_season_summary_contemp, by = "pop")

climate_for_pc_contemporary_values <- climate_for_pc_contemporary %>%
  select(-"pop")

str(climate_for_pc_contemporary)

# historic
climate_for_pc_historic = all_climate_historic %>% 
  select(pop, year, month, cwd, pck, ppt, str, tmn, tmx, elevation) %>%
  mutate(water_year = ifelse(month >= 9, year + 1, year)) %>% # Shift year starting from September
  dplyr::summarize(cwd = sum(cwd), ppt = sum(ppt), pck = sum(pck), tmin = mean(tmn), tmax = mean(tmx), elevation = mean(elevation), str = sum(str)) %>% 
  left_join(growing_season_summary_historic, by = "pop")

climate_for_pc_historic_values <- climate_for_pc_historic %>%
  select(-"pop")

str(climate_for_pc_historic)
str(climate_for_pc_historic_values)

# > make PCA ----

# contemporary

print(climate_for_pc_contemporary)

pc_contemporary = prcomp(climate_for_pc_contemporary_values, scale = TRUE)

summary(pc_contemporary)

pc_contemporary_scores = data.frame(pc_contemporary$x)

locs_pc_contemporary = cbind(climate_for_pc_contemporary, pc_contemporary_scores)

#scree plot
fviz_eig(pc_contemporary)

contemporary_loadings <- pc_contemporary[["rotation"]]

print(contemporary_loadings)

#historic

print(climate_for_pc_historic_values)

pc_historic = prcomp(climate_for_pc_historic_values, scale = TRUE)

summary(pc_historic)

pc_historic_scores = data.frame(pc_historic$x)

locs_pc_historic = cbind(climate_for_pc_historic, pc_historic_scores)

fviz_eig(pc_historic)
loadings_historic <- pc_historic[["rotation"]]

# graph of individuals - similar profiles grouped together

### contemporary

# PCA biplot with loading arrows and labeled points
fviz_pca_biplot(pc_contemporary,
                geom.ind = "point",        # Use points for individuals
                label = "var",             # Label variables (loadings)
                col.var = "black",           # Loadings (arrows) in red
                col.ind = locs_pc_contemporary$elevation, # Color points by elevation
                pointsize = 4,             # Make points larger
                repel = TRUE               # Avoid overlapping labels
) +
  # Custom color gradient: lower elevation (orange) → higher elevation (blue)
  scale_color_gradient(low = "orange", high = "blue", name = "Elevation (m)") +
    
  # Add text labels for locations
  geom_text(data = locs_pc_contemporary, aes(x = PC1, y = PC2, label = pop), 
            vjust = -0.5, size = 3, color = "black") +
  
# Set custom axis labels
labs(x = "Conemporary Climate PC1", 
     y = "Conemporary Climate PC2")

# graph of vars 

fviz_pca_var(pc_contemporary,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

### historic

# PCA biplot with loading arrows and labeled points
fviz_pca_biplot(pc_historic,
                geom.ind = "point",        # Use points for individuals
                label = "var",             # Label variables (loadings)
                col.var = "black",           # Loadings (arrows) in red
                col.ind = locs_pc_historic$elevation, # Color points by elevation
                pointsize = 4,             # Make points larger
                repel = TRUE               # Avoid overlapping labels
) +
  # Custom color gradient: lower elevation (orange) → higher elevation (blue)
  scale_color_gradient(low = "orange", high = "blue", name = "Elevation (m)") +
  
  # Add text labels for locations
  geom_text(data = locs_pc_historic, aes(x = PC1, y = PC2, label = pop), 
            vjust = -0.5, size = 3, color = "black")+
  
  # Set custom axis labels
  labs(x = "Historic Climate PC1", 
       y = "Historic Climate PC2")

# graph of vars

fviz_pca_var(pc_historic,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# Add identifying information to PCA scores

### contemporary
pc_data_with_id_contemporary <- cbind(ID = locs_pc_contemporary$pop, pc_contemporary_scores[, c("PC1", "PC2")])

loadings_contemporary = data.frame(varnames=rownames(pc_contemporary$rotation), pc_contemporary$rotation)

# Extract PCA scores for the first two principal components
pc1_2_scores_contemporary <- as.data.frame(pc_contemporary$x[, 1:2])

str(pc1_2_scores_contemporary)

# Add sample names (row names) to the PCA scores
pc1_2_scores_contemporary$names <- climate_for_pc_contemporary$pop

str(pc1_2_scores_contemporary)

# Add PC1 and PC2 values to the big dataframe
pca_scores_df_contemporary <- pc1_2_scores_contemporary %>%
  rename(Population = names)

str(pca_scores_df_contemporary)

all_climate_df_contemporary <- climate_for_pc_contemporary %>%
  rename(Population = pop) %>%
  left_join(pca_scores_df_contemporary, by = "Population") %>%
  mutate(contemporary_PC1 = PC1, contemporary_PC2 = PC2) %>%
  select(-c("PC1", "PC2"))

str(all_climate_df_contemporary)

### historic

pc_data_with_id_historic <- cbind(ID = locs_pc_historic$pop, pc_historic_scores[, c("PC1", "PC2")])

str(pc_data_with_id_historic)

loadings_historic = data.frame(varnames=rownames(pc_historic$rotation), pc_historic$rotation)

str(loadings_historic)

# Extract PCA scores for the first two principal components
pc1_2_scores_historic <- as.data.frame(pc_historic$x[, 1:2])

str(pc1_2_scores_historic)

# Add sample names (row names) to the PCA scores
pc1_2_scores_historic$names <- climate_for_pc_historic$pop

str(pc1_2_scores_historic)

# Add PC1 and PC2 values to the big dataframe
pca_scores_df_historic <- pc1_2_scores_historic %>%
  rename(Population = names)

str(pca_scores_df_historic)

all_climate_df_historic <- climate_for_pc_historic %>%
  rename(Population = pop) %>%
  left_join(pca_scores_df_historic, by = "Population") %>%
  mutate(historic_PC1 = PC1, historic_PC2 = PC2) %>%
  select(-c("PC1", "PC2"))

str(all_climate_df_historic)

# make a df with just the 2 PCa and pop

all_climate_df_historic_for_merge <- all_climate_df_historic %>%
  dplyr::select(Population, historic_PC1, historic_PC2)

str(all_climate_df_historic_for_merge)

### Merge dfs 

# merge contemporary and historic climate pcs

all_climate_and_PCs <- all_climate_df_contemporary %>%
  left_join(all_climate_df_historic_for_merge, by = "Population")

str(all_climate_and_PCs)

# merge climate data into the big data frame

data_with_clim <- data %>% 
  mutate(Population = replace(Population, Population == "YO10", "YOSE10")) %>%
  mutate(Population = replace(Population, Population == "MtSH", "MSH")) %>%
  left_join(all_climate_and_PCs, by = "Population") %>%
  drop_na()

print(data_with_clim)

# merge climate data with the mf means
mf_means_with_clim <- mf_means %>% 
  mutate(Population = replace(Population, Population == "YO10", "YOSE10")) %>%
  mutate(Population = replace(Population, Population == "MtSH", "MSH")) %>%
  left_join(all_climate_and_PCs, by = "Population") %>%
  drop_na()

# merge climate data with the pop means
pop_means_with_clim <- pop_means %>%
  mutate(Population = replace(Population, Population == "YO10", "YOSE10")) %>%
  mutate(Population = replace(Population, Population == "MtSH", "MSH")) %>%
  left_join(all_climate_and_PCs, by = "Population") %>%
  drop_na()

### SAVE DFs

# save mf means
write.csv(mf_means_with_clim, file = "./data/mf_means_with_clim.csv")

# save pops means
write.csv(pop_means_with_clim, file = "./data/pop_means_with_clim.csv")

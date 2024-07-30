### Ch2 Testing Intraspecifc Defense Trade-offs 
### Danielle De La Pascua
### 7-12-24

### Data Preparation ----

### > libraries ----
library(tidyr)
library(tidyverse)
library(dplyr)

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

# join mass data & location data with GSL data

#join gsl and mass/rack location data by rack & location
data <- left_join(biomass, GSL, by = c("Rack", "Location"))
summary(data)
dim(biomass)
dim(GSL)
dim(data)

#why are these not the same
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


### > save big data table ----

#write.csv(data, "~/GitHub/defense-tradeoffs-tortuosus/data/dl.csv")
write.csv(data, "./data/dw.csv")

# aggregating means

head(data)

mf_means <-data %>% 
  # Summarize by maternal family
  group_by(Population, treatment, mf) %>% 
  summarise(
    X5MSO = mean(X5MSO_10.2),
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
    Indole = mean(Indole_18.8)
  )

# make pop means
pop_means <- mf_means %>% 
  # Summarize by population
  group_by(Population, treatment) %>% 
  summarise(
    GSL_X5MSO = mean(X5MSO),
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
    GSL_Indole = mean(Indole)
  )
  
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

# save long version of data
write.csv(pop_means_long, "./data/pop_means_long.csv")


########## OLD CODE - Aggregating means (do in tidy way if we need means in the future)

### manipulate data so we have a df with mf means for CW and C treated plants
#JRG: why load this data here?  Is this left over from a previous version?  Why not keep working with the data you have been using above?
#dl <- read.csv("./data/dl-induced.csv")
#head(dl)
#dl <- dl[11:35]

# aggregate across mfs and by population
# across mfs
x3MSO_mf <- aggregate(X3MSO_5.2 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
OHAlkenyl_mf <- aggregate(OH.Alkenyl_6 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x4MSO_mf <- aggregate(X4MSO_7.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Allyl_mf <- aggregate(Allyl_7.4 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x5MSO_mf <- aggregate(X5MSO_10.2 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Butenyl_mf <- aggregate(Butenyl_12.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x3MT_mf <- aggregate(X3MT_13.6 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
MSOO_mf <- aggregate(MSOO_13.8 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
OHI3M_mf <- aggregate(OH.I3M_15.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol16_mf <- aggregate(Flavonol_16.1 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
x4MT_mf <- aggregate(X4MT._15.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
I3M_mf <- aggregate(I3M_16.7 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol17_mf <- aggregate(Flavonol_17.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol18_mf <- aggregate(Flavonol_18.5 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))
Indole_mf <- aggregate(Indole_18.8 ~ mf + Population + treatment, data = dl, FUN = function(x) mean(x, na.rm = TRUE))

# join aggregated mf dfs together Allyl_mf Butenyl_mf
mf_means <- full_join(x = Butenyl_mf, y = Allyl_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = x4MSO_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = OHAlkenyl_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = x5MSO_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = x3MSO_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = x3MT_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = MSOO_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = OHI3M_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = Flavonol16_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = x4MT_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = I3M_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = Flavonol17_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = Flavonol18_mf, by = c("mf", "Population", "treatment"))
mf_means <- full_join(x = mf_means, y = Indole_mf, by = c("mf", "Population", "treatment"))

# across pops
x3MSO_pop <- aggregate(X3MSO_5.2 ~ Population + treatment, data = x3MSO_mf, FUN = function(x) mean(x, na.rm = TRUE))
OHAlkenyl_pop <- aggregate(OH.Alkenyl_6 ~ Population + treatment, data = OHAlkenyl_mf, FUN = function(x) mean(x, na.rm = TRUE))
x4MSO_pop <- aggregate(X4MSO_7.1 ~ Population + treatment, data = x4MSO_mf, FUN = function(x) mean(x, na.rm = TRUE))
Allyl_pop <- aggregate(Allyl_7.4 ~ Population + treatment, data = Allyl_mf, FUN = function(x) mean(x, na.rm = TRUE))
x5MSO_pop <- aggregate(X5MSO_10.2 ~ Population + treatment, data = x5MSO_mf, FUN = function(x) mean(x, na.rm = TRUE))
Butenyl_pop<- aggregate(Butenyl_12.1 ~ Population + treatment, data = Butenyl_mf, FUN = function(x) mean(x, na.rm = TRUE))
x3MT_pop <- aggregate(X3MT_13.6 ~ Population + treatment, data = x3MT_mf, FUN = function(x) mean(x, na.rm = TRUE))
MSOO_pop <- aggregate(MSOO_13.8 ~ Population + treatment, data = MSOO_mf, FUN = function(x) mean(x, na.rm = TRUE))
OHI3M_pop<- aggregate(OH.I3M_15.1 ~ Population + treatment, data = OHI3M_mf, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol16_pop <- aggregate(Flavonol_16.1 ~ Population + treatment, data = Flavonol16_mf, FUN = function(x) mean(x, na.rm = TRUE))
x4MT_pop <- aggregate(X4MT._15.5 ~ Population + treatment, data = x4MT_mf, FUN = function(x) mean(x, na.rm = TRUE))
I3M_pop <- aggregate(I3M_16.7 ~ Population + treatment, data = I3M_mf, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol17_pop<- aggregate(Flavonol_17.5 ~ Population + treatment, data = Flavonol17_mf, FUN = function(x) mean(x, na.rm = TRUE))
Flavonol18_pop <- aggregate(Flavonol_18.5 ~ Population + treatment, data = Flavonol18_mf, FUN = function(x) mean(x, na.rm = TRUE))
Indole_pop <- aggregate(Indole_18.8 ~ Population + treatment, data = Indole_mf, FUN = function(x) mean(x, na.rm = TRUE))

# merge population means

# join aggregated mf dfs together x3MSO_pop
pop_means <- full_join(x = Allyl_pop, y = OHAlkenyl_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = x4MSO_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = x3MSO_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = x5MSO_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = Butenyl_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = x3MT_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = MSOO_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = OHI3M_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = Flavonol16_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = x4MT_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = I3M_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = Flavonol17_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = Flavonol18_pop, by = c("Population", "treatment"))
pop_means <- full_join(x = pop_means, y = Indole_pop, by = c("Population", "treatment"))

# separate df into two dfs by trt
control_df <- subset(d_induced, treatment == 'C')
cw_df <- subset(d_induced, treatment == 'CW')

# manually clean col names
write.csv(c_df, "~/GitHub/defense-tradeoffs-tortuosus/data/c_df.csv")
write.csv(cw_df, "~/GitHub/defense-tradeoffs-tortuosus/data/cw_df.csv")
write.csv(pop_means, "~/GitHub/defense-tradeoffs-tortuosus/data/pop_means.csv")
write.csv(mf_means, "~/GitHub/defense-tradeoffs-tortuosus/data/mf_means.csv")

# 
c_df <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/c_df.csv")
cw_df <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/cw_df.csv")

head(c_df)


# merge dfs
paired_means <- left_join(x = c_df, y = cw_df, by = c("mf", "Population"))

# save df
write.csv(paired_means, "~/GitHub/defense-tradeoffs-tortuosus/data/paired_means.csv")


############################### OLD CODE ################################


# create new variable with all compounds added together
#head(data)

#columns_to_sum <- c("3MSO_5.2", "OH-Alkenyl_6", "4MSO_7.1", "Allyl_7.4", "5MSO_10.2", "Butenyl_12.1", "3MT_13.6", "MSOO_13.8", "OH I3M_15.1", "4MT _15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")
#sums <- rowSums(data[, columns_to_sum], na.rm = TRUE)

#data$totalGSL <- sums

#head(data)

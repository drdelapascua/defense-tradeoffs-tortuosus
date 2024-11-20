### Defense Trade-offs: Qst-Fst Analysis
### Danielle De La Pascua
### 9-3-2024

### libraries ----
library(tidyr)
library(dplyr)
library(tidyverse)
library(QstFstComp)
library(lme4)

### pull & prepare data ----


### Genetic data (Fst)

# Load the CSV file
csv_data <- read.csv("data/ind_ids_reassigned.csv") %>%
  # select only whats needed
  select("ID", "tube_label") %>%
  # separate out the ID column
  separate(ID, into = c("pop", "mf", "rep"), sep = "-") %>%
  # filter out the ones we dont need
  select("pop", "tube_label", "mf")

# Load the RDS file
rds_data <- readRDS("data/genotype_data.rds")
rds_data <- as.data.frame(rds_data)

# Convert row names to a column named "ID"
#dummy_df <- rownames_to_column(dummy_df, var = "tube_label")
rds_data <- rownames_to_column(rds_data, var = "tube_label")

# Merge the data frames based on the "ID" column
merged_data <- merge(rds_data, csv_data, by = "tube_label") %>%
  select(-c("tube_label"))
#merged_data <- merge(dummy_df, csv_data, by = "tube_label") %>%
#  select(-c("tube_label"))

# filter out pops that are not in my experiment
overlap <- read.csv("data/overlap.csv")
pops_to_keep <- overlap$pop
filtered_merged_data <- merged_data %>%
  filter(pop %in% pops_to_keep)

# move pop and mf to frint of df
filtered_merged_data <- filtered_merged_data %>%
  select(mf, everything()) %>%
  select(pop, everything())

# replace pops with numbers
key <- overlap
colnames(key) <- c("pop_number", "pop_name") 
#key$pop <- key$pop_name #merging purposes
colnames(key)
colnames(filtered_merged_data)

#replace pop names with numbers using key
fst.dat <- right_join(filtered_merged_data, key, by = c("pop" = "pop_name")) %>%
  mutate(pop = ifelse(is.na(pop_number), pop, pop_number)) %>%
  select(-pop_number)
  
# combine pop and mf to make each mf a unique ID
fst.dat <- fst.dat %>%
  mutate(merged_col = paste(pop, mf, sep = "")) %>%
  select(-mf) %>%
  rename(mf = merged_col) %>%
  select(mf, everything()) %>%
  select(pop, everything())

### Trait Data (Qst)

#load trait data

#total gsl - danielle 11-19-24 test
GSL_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

#total aliphatic
aliphatic_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalaliphatic")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

#total indole
indole_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logindoles")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

# replace pops with numbers using key

# fix rownames

# total GSLs
colnames(GSL_totals) = c("pop", "mf", "totalGSL")

# Aliphatics
colnames(aliphatic_totals) = c("pop", "mf", "totalaliphatics")

# Indoles 
colnames(indole_totals) = c("pop", "mf", "totalindoles")

### merge with key to replace pops with #s

# total GSL
gsl.dat <- right_join(GSL_totals, key, by = c("pop" = "pop_name")) %>%
  mutate(pop = ifelse(is.na(pop_number), pop, pop_number)) %>%
  select(-pop_number)

# total aliphatics
aliphatic.dat <- right_join(aliphatic_totals, key, by = c("pop" = "pop_name")) %>%
  mutate(pop = ifelse(is.na(pop_number), pop, pop_number)) %>%
  select(-pop_number)

# total indoles
indole.dat <- right_join(indole_totals, key, by = c("pop" = "pop_name")) %>%
  mutate(pop = ifelse(is.na(pop_number), pop, pop_number)) %>%
  select(-pop_number)

# combine pop and mf to make each mf a unique ID

# total GSLs
gsl.dat <- gsl.dat %>%
  mutate(merged_col = paste(pop, mf, sep = "")) %>%
  select(pop, merged_col, totalGSL)
colnames(gsl.dat) <- c("pop", "mf", "totalGSL")

# total indoles
indole.dat <- indole.dat %>%
  mutate(merged_col = paste(pop, mf, sep = "")) %>%
  select(pop, merged_col, totalindoles)
colnames(indole.dat) <- c("pop", "mf", "totalindoles")

# total aliphatics
aliphatic.dat <- aliphatic.dat %>%
  mutate(merged_col = paste(pop, mf, sep = "")) %>%
  select(pop, merged_col, totalaliphatics)
colnames(aliphatic.dat) <- c("pop", "mf", "totalaliphatics")

### Run Qst-Fst Comp ----

# make sure everything is in numerical form
fst.dat[] <- lapply(fst.dat, as.numeric)
gsl.dat[] <- lapply(gsl.dat, as.numeric)
indole.dat[] <- lapply(indole.dat, as.numeric)
aliphatic.dat[] <- lapply(aliphatic.dat, as.numeric)

# run qst fst
fst.dat <- na.omit(fst.dat)
gsl.dat <- na.omit(gsl.dat)
any(is.na(fst.dat))
any(is.na(indole.dat))
head(fst.dat)

# only keep a subset of pops that are well sampled
pops_to_exclude <- c("13")
pops_to_include <- c("2", "6", "7", "8", "9")

gsl.dat.f <- gsl.dat %>% 
  filter(pop %in% pops_to_include)

fst.dat.f <- fst.dat %>% 
  filter(pop %in% pops_to_include)

numpops <- length(unique(gsl.dat.f$pop)) 
numpops <- length(unique(fst.dat.f$pop)) 


# troubleshooting

result.gsl <- QstFstComp::QstFstComp(fst.dat = fst.dat, qst.dat = gsl.dat, numpops = 13, breeding.design = "half.sib.dam")

result.indole <- QstFstComp::QstFstComp(fst.dat = fst.dat, qst.dat = indole.dat, numpops = 13, breeding.design = "half.sib.dam")

result.aliphatic <- QstFstComp::QstFstComp(fst.dat = fst.dat, qst.dat = aliphatic.dat, numpops = 13, breeding.design = "half.sib.dam")

# add more sims

result.gsl.10000sims <- QstFstComp::QstFstComp(fst.dat = fst.dat, qst.dat = gsl.dat, numpops = 13, nsim = 10000, breeding.design = "half.sib.dam")

result.gsl.50000sims <- QstFstComp::QstFstComp(fst.dat = fst.dat, qst.dat = gsl.dat, numpops = 13, nsim = 50000, breeding.design = "half.sib.dam")

# other method to calculate Qst ----

# total GSL
gsl.dat$pop <- factor(gsl.dat$pop)
fit = lmer(totalGSL ~ (1|pop), data = gsl.dat)
summary(fit)

# Get variance components
var_components = as.data.frame(VarCorr(fit))
sigma_B = var_components$vcov[1]  # Between-population variance
sigma_W = attr(VarCorr(fit), "sc")^2  # Within-population variance

# Calculate Qst
Qst_total = sigma_B / (sigma_B + 2 * sigma_W)
print(Qst_total)

# indoles
pops_to_exclude <- c("5", "13")
indole.dat.filtered <- indole.dat %>% 
  filter(!pop %in% pops_to_exclude)
indole.dat$pop <- factor(indole.dat$pop)
fit = lmer(totalindoles ~ (1|pop), data = indole.dat)
help('isSingular')
isSingular(fit) #what does this mean
summary(fit)

# Get variance components
var_components = as.data.frame(VarCorr(fit))
sigma_B = var_components$vcov[1]  # Between-population variance
sigma_W = attr(VarCorr(fit), "sc")^2  # Within-population variance

# Calculate Qst
Qst_indole = sigma_B / (sigma_B + 2 * sigma_W)
print(Qst_indole)

# aliphatics

aliphatic.dat$pop <- factor(aliphatic.dat$pop)
fit = lmer(totalaliphatics ~ (1|pop), data = aliphatic.dat)
summary(fit)

# Get variance components
var_components = as.data.frame(VarCorr(fit))
sigma_B = var_components$vcov[1]  # Between-population variance
sigma_W = attr(VarCorr(fit), "sc")^2  # Within-population variance

# Calculate Qst
Qst_aliphatic = sigma_B / (sigma_B + 2 * sigma_W)
print(Qst_aliphatic)


# Calculate confidence interval ? bootstrapping
# shuffle pop levels, do it 100 or 1000 times, 2-tailed, test where does this Qst fall
# Calculate pairwise estimates for all pairwise Qst, confidence interval for distribution 
# Calculate pairwise Fst and subtract them all - would be good to recalculate with my 13 pops



### calculate Qst ----

# trying package
library(devtools)
install_github("kjgilbert/QstFstComp")
library(QstFstComp)



###############################################################
# PAIRWISE COMPARISONS ----
###############################################################

# Calculate mean and variance within each population
indole_population_stats <- mf_means %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(Indole),
    PopVar = var(Indole),
    .groups = 'drop'
  )

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(indole_population_stats$Population),
                                Pop2 = unique(indole_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise Qst results
pairwise_qst <- data.frame(Pop1 = character(), Pop2 = character(), Qst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate Qst
for (i in 1:nrow(population_pairs)) {
  pair <- population_pairs[i, ]
  pop1_stats <- indole_population_stats %>% filter(Population == pair$Pop1)
  pop2_stats <- indole_population_stats %>% filter(Population == pair$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means <- mean(c(pop1_stats$PopMean, pop2_stats$PopMean))
  between_pop_var <- sum((c(pop1_stats$PopMean, pop2_stats$PopMean) - mean_of_pop_means)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var <- mean(c(pop1_stats$PopVar, pop2_stats$PopVar))
  
  # Calculate Qst
  qst <- between_pop_var / (between_pop_var + within_pop_var)
  
  # Store the result
  indole_pairwise_qst <- rbind(pairwise_qst, data.frame(Pop1 = pair$Pop1, Pop2 = pair$Pop2, Qst = qst))
}

left_join(indole_pairwise_qst, fst, by = c("Pop1", "Pop2"))

# do for all other compounds?

# outstanding questions & next steps
# > this assumes balanced half sib design - how do i modify the loop to account for unbalanced design?
# > how is this method different from pop-level fst? Literature says pairwise is ok, but having trouble finding examples of code
# >> this package may help? https://github.com/kjgilbert/QstFstComp 
# > driftsel - this package can do qpc type analysis on unbalanced half-sib design on different sets of individuals (genotype vs phenotype)


### other method? ----
install.packages("devtools")
library(devtools)
install_github("kjgilbert/QstFstComp")
library(QstFstComp)

?QstFstComp

# Calculate Qst for each trait
qst_result <- QstFstComp(data = mf_means, 
                         trait = "Indole", 
                         pop = "Population", 
                         model = "half.sib.dam", 
                         qst = TRUE) 

# print Qst values
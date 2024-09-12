### Defense Trade-offs: Qst-Fst Analysis
### Danielle De La Pascua
### 9-3-2024

### libraries ----
library(tidyr)
library(dplyr)
library(tidyverse)

### pull & prepare data ----

# move to prep script??

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
colnames(fst) <- c("Pop1", "Pop2", "Fst")

# trait data
mf_means <-  read.csv("./data/mf_means.csv") %>%
  select("Population", "mf", "treatment", "X3MSO", "OHAlkenyl", "X4MSO", "Allyl", "X5MSO", "Butenyl", "X3MT", "MSOO", "OHI3M", "X4MT", "Flavonol16", "I3M", "Flavonol17", "Flavonol18", "Indole")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SHA", "SQ1", "SQ3", "WL1", "WL2","WL3", "YO1", "YO10")) %>%
  filter(treatment == "C") %>% # filter for the controls
  select(-"treatment")

# 3 pops missing - yo1, sq3, sha  (low germination)

### calculate Qst ----

# trying package
library(devtools)
install_github("kjgilbert/QstFstComp")
library(QstFstComp)



# go through with one compound 

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
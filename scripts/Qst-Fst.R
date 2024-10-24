### Defense Trade-offs: Qst-Fst Analysis
### Danielle De La Pascua
### 9-3-2024

### libraries ----
library(tidyr)
library(dplyr)
library(tidyverse)

### pull & prepare data ----


### Genetic data (Fst)

genotype_data <- readRDS("~/GitHub/defense-tradeoffs-tortuosus/data/genotype_data.RDS")

genotype_data$label <- rownames(genotype_data)

ind_ids_reassigned <- read_csv("data/ind_ids_reassigned.csv")


### Trait Data (Qst)

#load trait data

#total gsl
GSL_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SHA", "SQ1", "SQ3", "WL1", "WL2","WL3", "YO1", "YO10"))

#total aliphatic
aliphatic_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalaliphatic")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SHA", "SQ1", "SQ3", "WL1", "WL2","WL3", "YO1", "YO10"))

#total indole
indole_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalindole")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SHA", "SQ1", "SQ3", "WL1", "WL2","WL3", "YO1", "YO10"))

# 3 pops missing - yo1, sq3, sha  (low germination)

### calculate Qst ----

# trying package
library(devtools)
install_github("kjgilbert/QstFstComp")
library(QstFstComp)

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
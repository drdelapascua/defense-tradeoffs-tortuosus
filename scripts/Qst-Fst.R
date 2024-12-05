### Defense Trade-offs: Qst-Fst Analysis
### Danielle De La Pascua
### 9-3-2024

### libraries ----
library(tidyr)
library(dplyr)
library(tidyverse)
library(QstFstComp)
library(lme4)
library(boot)
library(nlme)
library(ggplot2)

# Pst code ----

#### Total GSLs ----

# load data
GSL_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# Fit the random effects model

# log total GSLs
hist(GSL_totals$logGSL)

model_totalGSL <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals)

# PST for dist.r = 1
varpop.obs <- as.numeric(VarCorr(model_totalGSL)[1]) # Between-population variance
varres.obs <- as.numeric(VarCorr(model_totalGSL)[2]) # Within-population variance
pst <- varpop.obs / (varpop.obs + (2 * varres.obs))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals))
MS.bootpop<-matrix(NA,nrow=nperm,ncol=1)
MS.bootres<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs <- do.call(rbind, lapply(split(GSL_totals,GSL_totals$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs)
  
  # Store variances and calculate PST
  variance_components_totalGSLb <- as.data.frame(VarCorr(model_totalGSLboot))
  MS.bootpop[i] <- as.numeric(VarCorr(model_totalGSLboot)[1]) 
  MS.bootres[i] <- as.numeric(VarCorr(model_totalGSLboot)[2])
  tablePSTboot[i] <- MS.bootpop[i] / (MS.bootpop[i] + 2 * MS.bootres[i])
}


# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST <- mean(tablePSTboot, na.rm = TRUE)
se_PST <- sd(tablePSTboot, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot)))
pst_ci <- quantile(tablePSTboot, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST, "\n")
cat("Standard Error:", se_PST, "\n")
cat("95% Confidence Interval:", pst_ci, "\n")

# Calculate confidence intervals
quantpop <- quantile(MS.bootpop, c(0.025, 0.975))
MS.bootpop.min <- as.numeric(quantpop[1])
MS.bootpop.max <- as.numeric(quantpop[2])

quantres <- quantile(MS.bootres, c(0.025, 0.975))
MS.bootres.min <- as.numeric(quantres[1])
MS.bootres.max <- as.numeric(quantres[2])

# Calculate min and max Pst
PST.boot.min <- MS.bootpop.min / (MS.bootpop.min + 2 * MS.bootres.min)


# Plot PST as a function of r
PSTobs <- (dist.r * varpop.obs) / ((dist.r * varpop.obs) + (2 * varres.obs))
PSTobs.ord <- sort(PSTobs)
plot(
  dist.r, PSTobs.ord, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)

# PST min and max as functions of r
PSTmin <- (dist.r * min(MS.bootpop)) / ((dist.r * min(MS.bootpop)) + 2 * min(MS.bootres))
PSTmin.ord <- sort(PSTmin)
lines(dist.r, PSTmin.ord, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax <- (dist.r * max(MS.bootpop)) / ((dist.r * max(MS.bootpop)) + 2 * max(MS.bootres))
PSTmax.ord <- sort(PSTmax)
lines(dist.r, PSTmax.ord, type = "l", lty = "dotted", lwd = 2, col = "grey")

# Add horizontal lines
abline(h = 0.033, lty = "dashed", lwd = 2, col = "black")

# Plot histogram of bootstrapped PST distribution
hist(tablePSTboot, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = mean_PST, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Mean PST"), col = c("red", "blue"), lty = c("dashed", "solid"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### total indoles #### ----

# Read and preprocess data
indole_totals <- read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)

# Plot histogram of logindoles
hist(indole_totals$logindoles)

# Boxplot to visualize logindoles by population
ggplot(aes(x = Population, y = logindoles), data = indole_totals) +
  geom_boxplot()

# Fit mixed model for total indoles
model_totalindoles <- lme(logindoles ~ 1, random = ~1 | Population, data = indole_totals)

# Calculate observed variances
varpop.obs.indole <- as.numeric(VarCorr(model_totalindoles)[1])  # Between-population variance
varres.obs.indole <- as.numeric(VarCorr(model_totalindoles)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed <- varpop.obs.indole / (varpop.obs.indole + 2 * varres.obs.indole)

# Print PST for verification
print(paste("Observed PST:", pst_observed))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
indole_bootpop <- matrix(NA, nrow = nperm, ncol = 1)
indole_bootres <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_indole <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_indole <- do.call(rbind, lapply(split(indole_totals, indole_totals$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalindoles_boot <- lme(logindoles ~ 1, random = ~1 | Population, data = resampled_data_indole)
  
  # Store variances and calculate PST for the bootstrap sample
  indole_bootpop[i] <- as.numeric(VarCorr(model_totalindoles_boot)[1])
  indole_bootres[i] <- as.numeric(VarCorr(model_totalindoles_boot)[2])
  table_pst_boot_indole[i] <- indole_bootpop[i] / (indole_bootpop[i] + 2 * indole_bootres[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_indole_pst <- mean(table_pst_boot_indole, na.rm = TRUE)
se_indole_pst <- sd(table_pst_boot_indole, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_indole)))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_indole_pst))
print(paste("Standard Error Bootstrapped PST:", se_indole_pst))

# Calculate confidence intervals for bootstrapped variances
quantpop_indole <- quantile(indole_bootpop, c(0.025, 0.975))
quantres_indole <- quantile(indole_bootres, c(0.025, 0.975))

indole_bootpop_min <- as.numeric(quantpop_indole[1])
indole_bootpop_max <- as.numeric(quantpop_indole[2])

indole_bootres_min <- as.numeric(quantres_indole[1])
indole_bootres_max <- as.numeric(quantres_indole[2])

# Calculate min and max PST for confidence intervals
pst_boot_min_indole <- indole_bootpop_min / (indole_bootpop_min + 2 * indole_bootres_max)
pst_boot_max_indole <- indole_bootpop_max / (indole_bootpop_max + 2 * indole_bootres_min)

# Plot PST as a function of r
pst_obs_indole <- (dist_r * varpop.obs.indole) / ((dist_r * varpop.obs.indole) + 2 * varres.obs.indole)
pst_obs_ordered_indole <- sort(pst_obs_indole)

plot(
  dist_r, pst_obs_ordered_indole, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 Populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed"), lwd = c(3, 2, 2),
  legend = c("PST", "95% CI (Bootstrap)", "Upper Bound of 95% CI (FST)"),
  col = c("black", "grey", "black")
)

# Calculate min and max PST for confidence intervals
pst_boot_min_indole <- indole_bootpop_min / (indole_bootpop_min + 2 * indole_bootres_max)
pst_boot_max_indole <- indole_bootpop_max / (indole_bootpop_max + 2 * indole_bootres_min)

plot_data <- data.frame(
  PST = table_pst_boot_indole,
  LineType = rep("Bootstrapped PST", length(table_pst_boot_indole))
)

# Histogram plot of bootstrapped PST values with additional lines and labels
ggplot(data = plot_data, aes(x = PST)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = pst_observed, color = "Observed PST", linetype = "Observed PST"), size = 1) +
  geom_vline(aes(xintercept = mean_indole_pst, color = "Mean Bootstrapped PST", linetype = "Mean Bootstrapped PST"), size = 1) +
  geom_vline(aes(xintercept = pst_boot_min_indole, color = "95% CI Lower Bound", linetype = "95% CI Lower Bound"), size = 1) +
  geom_vline(aes(xintercept = pst_boot_max_indole, color = "95% CI Upper Bound", linetype = "95% CI Upper Bound"), size = 1) +
  scale_color_manual(values = c("Observed PST" = "red", 
                                "Mean Bootstrapped PST" = "blue", 
                                "95% CI Lower Bound" = "orange", 
                                "95% CI Upper Bound" = "orange")) +
  scale_linetype_manual(values = c("Observed PST" = "dashed", 
                                   "Mean Bootstrapped PST" = "solid", 
                                   "95% CI Lower Bound" = "dotted", 
                                   "95% CI Upper Bound" = "dotted")) +
  guides(
    color = guide_legend(title = "Legend", override.aes = list(size = 1, linetype = c("dashed", "solid", "dotted", "dotted"), 
                                                               color = c("red", "blue", "orange", "orange"))),
    linetype = guide_legend(override.aes = list(color = c("red", "blue", "orange", "orange")))
  ) +
  xlab("PST Value") + ylab("Frequency") +
  ggtitle("Histogram of Bootstrapped PST Values with Observed PST and Confidence Intervals") +
  theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10))

#### aliphatics ####

#total aliphatic
aliphatic_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logaliphatics")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# Fit the random effects model

# log total GSLs
hist(aliphatic_totals$logaliphatics)

model_totalaliphatics <- lme(logaliphatics ~ 1, random = ~1|Population, data = aliphatic_totals)

# PST for dist.r = 1
varpop.obs.aliphatic <- as.numeric(VarCorr(model_totalaliphatics)[1]) # Between-population variance
varres.obs.aliphatic <- as.numeric(VarCorr(model_totalaliphatics)[2]) # Within-population variance
aliphatic.pst <- varpop.obs.aliphatic / (varpop.obs.aliphatic + (2 * varres.obs.aliphatic))

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
aliphaticboot<- matrix(NA, nrow = nperm, ncol = nrow(aliphatic_totals))
aliphatic.bootpop<-matrix(NA,nrow=nperm,ncol=1)
aliphatic.bootres<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.aliphatic<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalaliphatic <- do.call(rbind, lapply(split(aliphatic_totals,aliphatic_totals$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalaliphatic.boot <- lme(logaliphatics ~ 1, random = ~1|Population, data = resampled_data_totalaliphatic)
  
  # Store variances and calculate PST
  aliphatic.bootpop[i] <- as.numeric(VarCorr(model_totalaliphatic.boot)[1]) 
  aliphatic.bootres[i] <- as.numeric(VarCorr(model_totalaliphatic.boot)[2])
  tablePSTboot.aliphatic[i] <- aliphatic.bootpop[i] / (aliphatic.bootpop[i] + 2 * aliphatic.bootres[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_aliphatic_PST <- mean(tablePSTboot.aliphatic, na.rm = TRUE)
se_aliphatic_PST <- sd(tablePSTboot.aliphatic, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.aliphatic)))

# Calculate confidence intervals
quantpop.aliphatic <- quantile(aliphatic.bootpop, c(0.025, 0.975))
aliphatic.bootpop.min <- as.numeric(quantpop.aliphatic[1])
aliphatic.bootpop.max <- as.numeric(quantpop.aliphatic[2])

quantres.aliphatic <- quantile(aliphatic.bootres, c(0.025, 0.975))
aliphatic.bootres.min <- as.numeric(quantres.aliphatic[1])
aliphatic.bootres.max <- as.numeric(quantres.aliphatic[2])

# Calculate min and max Pst
PST.boot.min.aliphatic <- aliphatic.bootpop.min / (aliphatic.bootpop.min + 2 * aliphatic.bootres.min)
PST.boot.max.aliphatic <- aliphatic.bootpop.max / (aliphatic.bootpop.max + 2 * aliphatic.bootres.max)


# Plot PST as a function of r
PSTobs.aliphatic <- (dist.r * varpop.obs.aliphatic) / ((dist.r * varpop.obs.aliphatic) + (2 * varres.obs.aliphatic))
PSTobs.ord.aliphatic <- sort(PSTobs.aliphatic)
plot(
  dist.r, PSTobs.ord.aliphatic, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.aliphatic <- (dist.r * aliphatic.bootpop.min) / ((dist.r * aliphatic.bootpop.min) + 2 * aliphatic.bootres.min)
PSTmin.ord.aliphatic <- sort(PSTmin.aliphatic)
lines(dist.r, PSTmin.ord.aliphatic, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.aliphatic <- (dist.r * aliphatic.bootpop.max) / ((dist.r * aliphatic.bootpop.max) + 2 * aliphatic.bootres.max)
PSTmax.ord.aliphatic <- sort(PSTmax.aliphatic)
lines(dist.r, PSTmax.ord.aliphatic, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.02, lty = "dotted", lwd = 2, col = "black")


### Qst Fst Comp Package ----

#### pull & prepare data ----

### Load Genetic data (Fst)

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

#### Trait Data (Qst) ----

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

#### Run Qst-Fst Comp ----

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

result.gsl$

qst.MS <- MeanSq.unbalanced.dam(qst.dat)




# Rishav and my work below

# total GSL

gsl.dat$pop <- factor(gsl.dat$pop)
fit = lmer(totalGSL ~ 1 + (1|pop) , data = gsl.dat)
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



### PAIRWISE COMPARISONS ----

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

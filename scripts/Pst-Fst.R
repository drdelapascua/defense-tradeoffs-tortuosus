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
library(adegenet)
library(hierfstat)
library(tidyverse)
library(ggthemes)
library(ggtext)
library(lme4)
library(ggrepel)
library(vegan)

### Constitutive defenses ----

#### Total GSLs ----

#### > load data ----

GSL_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# no extreme values 
GSL_totals_ex_rem <-  read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# no OH-I3M
GSL_totals_rem_ohi3m <-  read.csv("./data/mf_means_rem_ohi3m.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))


# Fit the random effects model

#### > with extreme values ----

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
  variance_components_totalGSL <- as.data.frame(VarCorr(model_totalGSLboot))
  MS.bootpop[i] <- as.numeric(VarCorr(model_totalGSLboot)[1]) 
  MS.bootres[i] <- as.numeric(VarCorr(model_totalGSLboot)[2])
  tablePSTboot[i] <- MS.bootpop[i] / (MS.bootpop[i] + 2 * MS.bootres[i])
}

# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST <- mean(tablePSTboot, na.rm = TRUE)
se_PST <- sd(tablePSTboot, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot)))
pst_ci <- quantile(tablePSTboot, c(0.025, 0.5, 0.975), na.rm = TRUE)

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST, "\n")
cat("Standard Error:", se_PST, "\n")
cat("95% Confidence Interval:", pst_ci, "\n")

# Calculate confidence intervals
quantpop <- quantile(tablePSTboot, c(0.025, 0.975), na.rm = TRUE)
MS.bootpop.min <- as.numeric(quantpop[1])
MS.bootpop.max <- as.numeric(quantpop[2])

quantres <- quantile(tablePSTboot, c(0.025, 0.975), na.rm = TRUE)
MS.bootres.min <- as.numeric(quantres[1])
MS.bootres.max <- as.numeric(quantres[2])

# order PST observations
PSTobs.ord.df <- as.data.frame(PSTobs.ord)

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
PSTmin <- dist.r * MS.bootpop.min
PSTmin.ord <- sort(PSTmin)
lines(dist.r, PSTmin.ord, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax <- dist.r * MS.bootpop.max
PSTmax.ord <- sort(PSTmax)
lines(dist.r, PSTmax.ord, type = "l", lty = "dotted", lwd = 2, col = "grey")

# Add horizontal lines
abline(h = 0.03, lty = "dashed", lwd = 2, col = "black")

# Plot histogram of bootstrapped PST distribution
hist(tablePSTboot, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = pst, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > Without extreme values ----

# log total GSLs 
hist(GSL_totals_ex_rem$logGSL)

model_totalGSL_ex_rem <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals_ex_rem)

# PST for dist.r = 1
varpop.obs.ex.rem <- as.numeric(VarCorr(model_totalGSL_ex_rem)[1]) # Between-population variance
varres.obs.ex.rem <- as.numeric(VarCorr(model_totalGSL_ex_rem)[2]) # Within-population variance
pst_ex_rem <- varpop.obs.ex.rem / (varpop.obs.ex.rem + (2 * varres.obs.ex.rem))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot.ex.rem<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals_ex_rem))
MS.bootpop.ex.rem<-matrix(NA,nrow=nperm,ncol=1)
MS.bootres.ex.rem<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.ex.rem<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs_ex_rem <- do.call(rbind, lapply(split(GSL_totals_ex_rem,GSL_totals_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot_ex_rem <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs_ex_rem)
  
  # Store variances and calculate PST
  variance_components_totalGSLb_ex_rem <- as.data.frame(VarCorr(model_totalGSLboot_ex_rem))
  MS.bootpop.ex.rem[i] <- as.numeric(VarCorr(model_totalGSLboot_ex_rem)[1]) 
  MS.bootres.ex.rem[i] <- as.numeric(VarCorr(model_totalGSLboot_ex_rem)[2])
  tablePSTboot.ex.rem[i] <- MS.bootpop.ex.rem[i] / (MS.bootpop.ex.rem[i] + 2 * MS.bootres.ex.rem[i])
}

# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST_ex_rem <- mean(tablePSTboot.ex.rem, na.rm = TRUE)
se_PST_ex_rem <- sd(tablePSTboot.ex.rem, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.ex.rem)))
pst_ci_ex_rem <- quantile(tablePSTboot.ex.rem, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST_ex_rem, "\n")
cat("Standard Error:", se_PST_ex_rem, "\n")
cat("95% Confidence Interval:", pst_ci_ex_rem, "\n")

#### > Without OH-I3M ----

# log total GSLs 
hist(GSL_totals_rem_ohi3m$logGSL)

model_totalGSL_rem_ohi3m <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals_rem_ohi3m)

# PST for dist.r = 1
varpop.obs.rem.ohi3m <- as.numeric(VarCorr(model_totalGSL_rem_ohi3m)[1]) # Between-population variance
varres.obs.rem.ohi3m <- as.numeric(VarCorr(model_totalGSL_rem_ohi3m)[2]) # Within-population variance
pst_rem_ohi3m <- varpop.obs.rem.ohi3m / (varpop.obs.rem.ohi3m + (2 * varres.obs.ex.rem))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot.ex.rem.ohi3m<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals_rem_ohi3m))
MS.bootpop.rem.ohi3m <-matrix(NA,nrow=nperm,ncol=1)
MS.bootres.rem.ohi3m <-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.rem.ohi3m<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs_rem_ohi3m <- do.call(rbind, lapply(split(GSL_totals_rem_ohi3m,GSL_totals_rem_ohi3m$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot_rem_ohi3m <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs_rem_ohi3m)
  
  # Store variances and calculate PST
  variance_components_totalGSLb_rem_ohi3m <- as.data.frame(VarCorr(model_totalGSLboot_rem_ohi3m))
  MS.bootpop.rem.ohi3m[i] <- as.numeric(VarCorr(model_totalGSLboot_rem_ohi3m)[1]) 
  MS.bootres.rem.ohi3m[i] <- as.numeric(VarCorr(model_totalGSLboot_rem_ohi3m)[2])
  tablePSTboot.rem.ohi3m[i] <- MS.bootpop.rem.ohi3m[i] / (MS.bootpop.rem.ohi3m[i] + 2 * MS.bootres.rem.ohi3m[i])
}

# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST_rem_ohi3m <- mean(tablePSTboot.rem.ohi3m, na.rm = TRUE)
se_PST_rem_ohi3m <- sd(tablePSTboot.rem.ohi3m, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.rem.ohi3m)))
pst_ci_rem_ohi3m <- quantile(tablePSTboot.rem.ohi3m, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST_rem_ohi3m, "\n")
cat("Standard Error:", se_PST_rem_ohi3m, "\n")
cat("95% Confidence Interval:", pst_ci_rem_ohi3m, "\n")


#### total indoles ----

#### > load data ----

# Read and preprocess data
indole_totals <- read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)

indole_totals_ex_rem <- read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)
 
indole_totals_rem_ohi3m <- read.csv("./data/mf_means_rem_ohi3m.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)

#### > with extreme values -----

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
pst_observed_indole <- varpop.obs.indole / (varpop.obs.indole + 2 * varres.obs.indole)

# Print PST for verification
print(paste("Observed PST:", pst_observed_indole))

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

# PST min and max as functions of r
PSTmin_indole <- dist.r*pst_boot_min_indole
PSTmin.ord_indole <- sort(PSTmin_indole)
lines(dist.r, PSTmin.ord_indole, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax_indole <- dist.r * pst_boot_max_indole
PSTmax.ord_indole <- sort(PSTmax_indole)
lines(dist.r, PSTmax.ord_indole, type = "l", lty = "dotted", lwd = 2, col = "grey")

# Add horizontal lines
abline(h = 0.03, lty = "dashed", lwd = 2, col = "black")

# Histogram plot of bootstrapped PST values with additional lines and labels
pst_ci_indole <- quantile(table_pst_boot_indole, c(0.025, 0.5, 0.975))
hist(table_pst_boot_indole, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci_indole, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = pst_observed_indole, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > without extreme values ----

# Plot histogram of logindoles
hist(indole_totals_ex_rem$logindoles)

# Boxplot to visualize logindoles by population
ggplot(aes(x = Population, y = logindoles), data = indole_totals_ex_rem) +
  geom_boxplot()

# Fit mixed model for total indoles
model_totalindoles_ex_rem <- lme(logindoles ~ 1, random = ~1 | Population, data = indole_totals_ex_rem)

# Calculate observed variances
varpop.obs.indole.ex.rem <- as.numeric(VarCorr(model_totalindoles_ex_rem)[1])  # Between-population variance
varres.obs.indole.ex.rem <- as.numeric(VarCorr(model_totalindoles_ex_rem)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_indole_ex_rem <- varpop.obs.indole.ex.rem / (varpop.obs.indole.ex.rem + 2 * varres.obs.indole.ex.rem)

# Print PST for verification
print(paste("Observed PST:", pst_observed_indole_ex_rem))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
indole_bootpop_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
indole_bootres_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_indole_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_indole_ex_rem <- do.call(rbind, lapply(split(indole_totals_ex_rem, indole_totals_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalindoles_boot_ex_rem <- lme(logindoles ~ 1, random = ~1 | Population, data = resampled_data_indole_ex_rem)
  
  # Store variances and calculate PST for the bootstrap sample
  indole_bootpop_ex_rem[i] <- as.numeric(VarCorr(model_totalindoles_boot_ex_rem)[1])
  indole_bootres_ex_rem[i] <- as.numeric(VarCorr(model_totalindoles_boot_ex_rem)[2])
  table_pst_boot_indole_ex_rem[i] <- indole_bootpop_ex_rem[i] / (indole_bootpop_ex_rem[i] + 2 * indole_bootres_ex_rem[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_indole_pst_ex_rem <- mean(table_pst_boot_indole_ex_rem, na.rm = TRUE)
se_indole_pst_ex_rem <- sd(table_pst_boot_indole_ex_rem, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_indole_ex_rem)))
pst_ci_indole_ex_rem <- quantile(table_pst_boot_indole_ex_rem, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_indole_pst_ex_rem))
print(paste("Standard Error Bootstrapped PST:", se_indole_pst_ex_rem))
cat("95% Confidence Interval:", pst_ci_indole_ex_rem, "\n")


#### total aliphatics ----

#### > load data ----

#total aliphatic
aliphatic_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logaliphatics")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

aliphatic_totals_ex_rem <- read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalaliphatic", "logaliphatics") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

aliphatic_totals_rem_ohi3m <- read.csv("./data/mf_means_rem_ohi3m.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalaliphatic", "logaliphatics") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

#### > with extreme values ----

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
quantpop.aliphatic <- quantile(tablePSTboot.aliphatic, c(0.025, 0.975))
aliphatic.bootpop.min <- as.numeric(quantpop.aliphatic[1])
aliphatic.bootpop.max <- as.numeric(quantpop.aliphatic[2])

quantres.aliphatic <- quantile(tablePSTboot.aliphatic, c(0.025, 0.975))
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
PSTmin.aliphatic <- dist.r*aliphatic.bootpop.min
PSTmin.ord.aliphatic <- sort(PSTmin.aliphatic)
lines(dist.r, PSTmin.ord.aliphatic, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.aliphatic <- dist.r * aliphatic.bootpop.max
PSTmax.ord.aliphatic <- sort(PSTmax.aliphatic)
lines(dist.r, PSTmax.ord.aliphatic, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")

# Histogram plot of bootstrapped PST values with additional lines and labels
pst_ci_aliphatic <- quantile(tablePSTboot.aliphatic, c(0.025, 0.5, 0.975))
hist(tablePSTboot.aliphatic, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci_aliphatic, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = aliphatic.pst, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > without extreme values ----

# Plot histogram of logaliphatics
hist(aliphatic_totals_ex_rem$logaliphatics)
hist(aliphatic_totals$logaliphatics)

# Boxplot to visualize logaliphatics by population
ggplot(aes(x = Population, y = logaliphatics), data = aliphatic_totals_ex_rem) +
  geom_boxplot()

# normal data for comparison
ggplot(aes(x = Population, y = logaliphatics), data = aliphatic_totals_ex_rem) +
  geom_boxplot()

# Fit mixed model for total aliphatics
model_totalaliphatics_ex_rem <- lme(logaliphatics ~ 1, random = ~1 | Population, data = aliphatic_totals_ex_rem)

# Calculate observed variances
varpop.obs.aliphatic.ex.rem <- as.numeric(VarCorr(model_totalaliphatics_ex_rem)[1])  # Between-population variance
varres.obs.aliphatic.ex.rem <- as.numeric(VarCorr(model_totalaliphatics_ex_rem)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_aliphatic_ex_rem <- varpop.obs.aliphatic.ex.rem / (varpop.obs.aliphatic.ex.rem + 2 * varres.obs.aliphatic.ex.rem)

# Print PST for verification
print(paste("Observed PST:", pst_observed_aliphatic_ex_rem))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
aliphatic_bootpop_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
aliphatic_bootres_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_aliphatic_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_aliphatic_ex_rem <- do.call(rbind, lapply(split(aliphatic_totals_ex_rem, aliphatic_totals_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalaliphatics_boot_ex_rem <- lme(logaliphatics ~ 1, random = ~1 | Population, data = resampled_data_aliphatic_ex_rem)
  
  # Store variances and calculate PST for the bootstrap sample
  aliphatic_bootpop_ex_rem[i] <- as.numeric(VarCorr(model_totalaliphatics_boot_ex_rem)[1])
  aliphatic_bootres_ex_rem[i] <- as.numeric(VarCorr(model_totalaliphatics_boot_ex_rem)[2])
  table_pst_boot_aliphatic_ex_rem[i] <- aliphatic_bootpop_ex_rem[i] / (aliphatic_bootpop_ex_rem[i] + 2 * aliphatic_bootres_ex_rem[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_aliphatic_pst_ex_rem <- mean(table_pst_boot_aliphatic_ex_rem, na.rm = TRUE)
se_aliphatic_pst_ex_rem <- sd(table_pst_boot_aliphatic_ex_rem, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_aliphatic_ex_rem)))
pst_ci_aliphatic_ex_rem <- quantile(table_pst_boot_aliphatic_ex_rem, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_aliphatic_pst_ex_rem))
print(paste("Standard Error Bootstrapped PST:", se_aliphatic_pst_ex_rem))
cat("95% Confidence Interval:", pst_ci_aliphatic_ex_rem, "\n")

#### > without OH-I3M ----

# Plot histogram of logaliphatics
hist(aliphatic_totals_rem_ohi3m$logaliphatics)
# compare to regular plot
hist(aliphatic_totals$logaliphatics)

# Boxplot to visualize logaliphatics by population
ggplot(aes(x = Population, y = logaliphatics), data = aliphatic_totals_rem_ohi3m) +
  geom_boxplot()
# compare to regular plot
ggplot(aes(x = Population, y = logaliphatics), data = aliphatic_totals) +
  geom_boxplot()

# Fit mixed model for total aliphatics
model_totalaliphatics_rem_ohi3m <- lme(logaliphatics ~ 1, random = ~1 | Population, data = aliphatic_totals_rem_ohi3m)

# Calculate observed variances
varpop.obs.aliphatic.rem.ohi3m <- as.numeric(VarCorr(model_totalaliphatics_rem_ohi3m)[1])  # Between-population variance
varres.obs.aliphatic.rem.ohi3m <- as.numeric(VarCorr(model_totalaliphatics_rem_ohi3m)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_aliphatic_rem_ohi3m <- varpop.obs.aliphatic.rem.ohi3m / (varpop.obs.aliphatic.rem.ohi3m + 2 * varres.obs.aliphatic.rem.ohi3m)

# Print PST for verification
print(paste("Observed PST:", pst_observed_aliphatic_rem_ohi3m))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
aliphatic_bootpop_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)
aliphatic_bootres_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_aliphatic_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_aliphatic_rem_ohi3m <- do.call(rbind, lapply(split(aliphatic_totals_rem_ohi3m, aliphatic_totals_rem_ohi3m$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalaliphatics_boot_rem_ohi3m <- lme(logaliphatics ~ 1, random = ~1 | Population, data = resampled_data_aliphatic_rem_ohi3m)
  
  # Store variances and calculate PST for the bootstrap sample
  aliphatic_bootpop_rem_ohi3m[i] <- as.numeric(VarCorr(model_totalaliphatics_boot_rem_ohi3m)[1])
  aliphatic_bootres_rem_ohi3m[i] <- as.numeric(VarCorr(model_totalaliphatics_boot_rem_ohi3m)[2])
  table_pst_boot_aliphatic_rem_ohi3m[i] <- aliphatic_bootpop_rem_ohi3m[i] / (aliphatic_bootpop_rem_ohi3m[i] + 2 * aliphatic_bootres_rem_ohi3m[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_aliphatic_pst_rem_ohi3m <- mean(table_pst_boot_aliphatic_rem_ohi3m, na.rm = TRUE)
se_aliphatic_pst_rem_ohi3m <- sd(table_pst_boot_aliphatic_rem_ohi3m, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_aliphatic_rem_ohi3m)))
pst_ci_aliphatic_rem_ohi3m <- quantile(table_pst_boot_aliphatic_rem_ohi3m, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_aliphatic_pst_rem_ohi3m))
print(paste("Standard Error Bootstrapped PST:", se_aliphatic_pst_rem_ohi3m))
cat("95% Confidence Interval:", pst_ci_aliphatic_rem_ohi3m, "\n")

#### total flavonoid ----


#### > load data ----

flavonoid_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logflavonoids")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

flavonoid_totals_ex_rem <- read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "C") %>%  # Filter for control treatment
  select("Population", "mf", "totalflavonoid", "logflavonoids") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))


#### > with extreme values ----

# Fit the random effects model

hist(flavonoid_totals$logflavonoids)

model_totalflavonoids <- lme(logflavonoids ~ 1, random = ~1|Population, data = flavonoid_totals)

# PST for dist.r = 1
varpop.obs.flavonoid <- as.numeric(VarCorr(model_totalflavonoids)[1]) # Between-population variance
varres.obs.flavonoid <- as.numeric(VarCorr(model_totalflavonoids)[2]) # Within-population variance
flavonoid.pst <- varpop.obs.flavonoid / (varpop.obs.flavonoid + (2 * varres.obs.flavonoid))

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
flavonoidboot<- matrix(NA, nrow = nperm, ncol = nrow(flavonoid_totals))
flavonoid.bootpop<-matrix(NA,nrow=nperm,ncol=1)
flavonoid.bootres<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.flavonoid<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalflavonoid <- do.call(rbind, lapply(split(flavonoid_totals,flavonoid_totals$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalflavonoid.boot <- lme(logflavonoids ~ 1, random = ~1|Population, data = resampled_data_totalflavonoid)
  
  # Store variances and calculate PST
  flavonoid.bootpop[i] <- as.numeric(VarCorr(model_totalflavonoid.boot)[1]) 
  flavonoid.bootres[i] <- as.numeric(VarCorr(model_totalflavonoid.boot)[2])
  tablePSTboot.flavonoid[i] <- flavonoid.bootpop[i] / (flavonoid.bootpop[i] + 2 * flavonoid.bootres[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_flavonoid_PST <- mean(tablePSTboot.flavonoid, na.rm = TRUE)
se_flavonoid_PST <- sd(tablePSTboot.flavonoid, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.flavonoid)))

# Calculate confidence intervals
boot.quant.flavonoid <- quantile(tablePSTboot.flavonoid, c(0.025, 0.975))
flavonoid.bootpop.min <- as.numeric(boot.quant.flavonoid[1])
flavonoid.bootpop.max <- as.numeric(boot.quant.flavonoid[2])

# plot it
PSTobs.flavonoid <- (dist.r * varpop.obs) / ((dist.r * varpop.obs) + (2 * varres.obs))
PSTobs.ord <- sort(PSTobs.flavonoid)
plot(
  dist.r, PSTobs.ord.flavonoid, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.flavonoid <- dist.r*flavonoid.bootpop.min
PSTmin.ord.flavonoid <- sort(PSTmin.flavonoid)
lines(dist.r, PSTmin.ord.flavonoid, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.flavonoid <- dist.r * flavonoid.bootpop.max
PSTmax.ord.flavonoid <- sort(PSTmax.flavonoid)
lines(dist.r, PSTmax.ord.flavonoid, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")

# Histogram plot of bootstrapped PST values with additional lines and labels
#pst_ci_flavonoid <- quantile(tablePSTboot.flavonoid, c(0.025, 0.5, 0.975))
#hist(tablePSTboot.flavonoid, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
#     xlab = "PST", border = "black")
#abline(v = pst_ci_flavonoid, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
#abline(v = flavonoid.pst, col = "blue", lwd = 2, lty = "solid")
#legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
#abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > without extreme values ----

# Plot histogram of logflavonoids
hist(flavonoid_totals_ex_rem$logflavonoids)
hist(flavonoid_totals$logflavonoids)

# Boxplot to visualize logflavonoids by population
ggplot(aes(x = Population, y = logflavonoids), data = flavonoid_totals_ex_rem) +
  geom_boxplot()

# normal data for comparison
ggplot(aes(x = Population, y = logflavonoids), data = flavonoid_totals_ex_rem) +
  geom_boxplot()

# Fit mixed model for total flavonoids
model_totalflavonoids_ex_rem <- lme(logflavonoids ~ 1, random = ~1 | Population, data = flavonoid_totals_ex_rem)

# Calculate observed variances
varpop.obs.flavonoid.ex.rem <- as.numeric(VarCorr(model_totalflavonoids_ex_rem)[1])  # Between-population variance
varres.obs.flavonoid.ex.rem <- as.numeric(VarCorr(model_totalflavonoids_ex_rem)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_flavonoid_ex_rem <- varpop.obs.flavonoid.ex.rem / (varpop.obs.flavonoid.ex.rem + 2 * varres.obs.flavonoid.ex.rem)

# Print PST for verification
print(paste("Observed PST:", pst_observed_flavonoid_ex_rem))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
flavonoid_bootpop_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
flavonoid_bootres_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_flavonoid_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_flavonoid_ex_rem <- do.call(rbind, lapply(split(flavonoid_totals_ex_rem, flavonoid_totals_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalflavonoids_boot_ex_rem <- lme(logflavonoids ~ 1, random = ~1 | Population, data = resampled_data_flavonoid_ex_rem)
  
  # Store variances and calculate PST for the bootstrap sample
  flavonoid_bootpop_ex_rem[i] <- as.numeric(VarCorr(model_totalflavonoids_boot_ex_rem)[1])
  flavonoid_bootres_ex_rem[i] <- as.numeric(VarCorr(model_totalflavonoids_boot_ex_rem)[2])
  table_pst_boot_flavonoid_ex_rem[i] <- flavonoid_bootpop_ex_rem[i] / (flavonoid_bootpop_ex_rem[i] + 2 * flavonoid_bootres_ex_rem[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_flavonoid_pst_ex_rem <- mean(table_pst_boot_flavonoid_ex_rem, na.rm = TRUE)
se_flavonoid_pst_ex_rem <- sd(table_pst_boot_flavonoid_ex_rem, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_flavonoid_ex_rem)))
pst_ci_flavonoid_ex_rem <- quantile(table_pst_boot_flavonoid_ex_rem, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_flavonoid_pst_ex_rem))
print(paste("Standard Error Bootstrapped PST:", se_flavonoid_pst_ex_rem))
cat("95% Confidence Interval:", pst_ci_flavonoid_ex_rem, "\n")


#### Shannon diversity ----


#### > load data ----

shannon_totals <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "shannon_diversity")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

str(shannon_totals)

#### > with extreme values ----

# Fit the random effects model

model_shannon <- lme(shannon_diversity ~ 1, random = ~1|Population, data = shannon_totals)


as.numeric(VarCorr(model_shannon)[1]) 

# PST for dist.r = 1
varpop.obs.shannon <- as.numeric(VarCorr(model_shannon)[1]) # Between-population variance
                                
varres.obs.shannon <- as.numeric(VarCorr(model_shannon)[2]) # Within-population variance
shannon.pst <- varpop.obs.shannon / (varpop.obs.shannon + (2 * varres.obs.shannon))

print(shannon.pst)

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
shannonboot<- matrix(NA, nrow = nperm, ncol = nrow(shannon_totals))
shannon.bootpop<-matrix(NA,nrow=nperm,ncol=1)
shannon.bootres<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.shannon<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_shannon <- do.call(rbind, lapply(split(shannon_totals, shannon_totals$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_shannon.boot <- lme(shannon_diversity ~ 1, random = ~1|Population, data = resampled_data_shannon)
  
  # Store variances and calculate PST
  shannon.bootpop[i] <- as.numeric(VarCorr(model_shannon.boot)[1]) 
  shannon.bootres[i] <- as.numeric(VarCorr(model_shannon.boot)[2])
  tablePSTboot.shannon[i] <- shannon.bootpop[i] / (shannon.bootpop[i] + 2 * shannon.bootres[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_shannon_PST <- mean(tablePSTboot.shannon, na.rm = TRUE)
se_shannon_PST <- sd(tablePSTboot.shannon, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.flavonoid)))

# Calculate confidence intervals
boot.quant.shannon <- quantile(tablePSTboot.shannon, c(0.025, 0.975))
shannon.bootpop.min <- as.numeric(boot.quant.shannon[1])
shannon.bootpop.max <- as.numeric(boot.quant.shannon[2])
print(shannon.bootpop.max)
print(shannon.bootpop.min)

# plot it
varpop.obs.shannon / (varpop.obs.shannon + (2 * varres.obs.shannon))
PSTobs.shannon <- (dist.r * varpop.obs.shannon) / ((dist.r * varpop.obs.shannon) + (2 * varres.obs.shannon))
PSTobs.ord <- sort(PSTobs.shannon)
plot(
  dist.r, PSTobs.shannon, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.shannon <- dist.r*shannon.bootpop.min
PSTmin.ord.shannon <- sort(PSTmin.shannon)
lines(dist.r, PSTmin.ord.shannon, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.shannon <- dist.r * shannon.bootpop.max
PSTmax.ord.shannon <- sort(PSTmax.shannon)
lines(dist.r, PSTmax.ord.shannon, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")


# Induced defenses ----

#### Total GSLs ----

#### > load data ----
GSL_totals_ind <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# no extreme values 
GSL_totals_ind_ex_rem <-  read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# no OH-I3M
GSL_totals_ind_rem_ohi3m <-  read.csv("./data/mf_means_rem_ohi3m.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))


#### > with extreme values ----

# log total GSLs
hist(GSL_totals_ind$logGSL)

model_totalGSL_ind <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals_ind)

# PST for dist.r = 1
varpop.obs.ind <- as.numeric(VarCorr(model_totalGSL_ind)[1]) # Between-population variance
varres.obs.ind <- as.numeric(VarCorr(model_totalGSL_ind)[2]) # Within-population variance
pst.total.ind <- varpop.obs.ind / (varpop.obs.ind + (2 * varres.obs.ind))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot.ind<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals_ind))
MS.bootpop.ind<-matrix(NA,nrow=nperm,ncol=1)
MS.bootres.ind<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.ind<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs_ind <- do.call(rbind, lapply(split(GSL_totals_ind,GSL_totals_ind$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot_ind <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs_ind)
  
  # Store variances and calculate PST
  variance_components_totalGSLb_ind <- as.data.frame(VarCorr(model_totalGSLboot_ind))
  MS.bootpop.ind[i] <- as.numeric(VarCorr(model_totalGSLboot_ind)[1]) 
  MS.bootres.ind[i] <- as.numeric(VarCorr(model_totalGSLboot_ind)[2])
  tablePSTboot.ind[i] <- MS.bootpop.ind[i] / (MS.bootpop.ind[i] + 2 * MS.bootres.ind[i])
}


# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST_ind <- mean(tablePSTboot.ind, na.rm = TRUE)
se_PST_ind <- sd(tablePSTboot.ind, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.ind)))
pst_ci_ind <- quantile(tablePSTboot.ind, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST_ind, "\n")
cat("Standard Error:", se_PST_ind, "\n")
cat("95% Confidence Interval:", pst_ci_ind, "\n")

# Calculate confidence intervals
quantpop_ind <- quantile(MS.bootpop.ind, c(0.025, 0.975))
MS.bootpop.min.ind <- as.numeric(quantpop_ind[1])
MS.bootpop.max.ind <- as.numeric(quantpop_ind[2])

quantres_ind <- quantile(MS.bootres.ind, c(0.025, 0.975))
MS.bootres.min.ind <- as.numeric(quantres_ind[1])
MS.bootres.max.ind <- as.numeric(quantres_ind[2])

# Calculate min and max Pst
PST.boot.min.ind <- MS.bootpop.min.ind / (MS.bootpop.min.ind + 2 * MS.bootres.min.ind)


# Plot PST as a function of r
PSTobs.ind <- (dist.r * varpop.obs.ind) / ((dist.r * varpop.obs.ind) + (2 * varres.obs.ind))
PSTobs.ord.ind <- sort(PSTobs.ind)
plot(
  dist.r, PSTobs.ord.ind, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)

# PST min and max as functions of r
PSTmin.ind <- (dist.r * min(MS.bootpop.ind)) / ((dist.r * min(MS.bootpop.ind)) + 2 * min(MS.bootres.ind))
PSTmin.ord.ind <- sort(PSTmin.ind)
lines(dist.r, PSTmin.ord.ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.ind <- (dist.r * max(MS.bootpop.ind)) / ((dist.r * max(MS.bootpop.ind)) + 2 * max(MS.bootres.ind))
PSTmax.ord.ind <- sort(PSTmax.ind)
lines(dist.r, PSTmax.ord.ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

# Add horizontal lines
abline(h = 0.03, lty = "dashed", lwd = 2, col = "black")

# Plot histogram of bootstrapped PST distribution
hist(tablePSTboot.ind, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci_ind, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = pst.total.ind, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > without extreme values ----

# log total GSLs 
hist(GSL_totals_ind_ex_rem$logGSL)

model_totalGSL_ex_rem <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals_ind_ex_rem)

# PST for dist.r = 1
varpop.obs.ind.ex.rem <- as.numeric(VarCorr(model_totalGSL_ex_rem)[1]) # Between-population variance
varres.obs.ind.ex.rem <- as.numeric(VarCorr(model_totalGSL_ex_rem)[2]) # Within-population variance
pst_ex_rem <- varpop.obs.ind.ex.rem / (varpop.obs.ind.ex.rem + (2 * varres.obs.ex.rem))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot.ind.ex.rem<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals_ind_ex_rem))
MS.bootpop.ind.ex.rem<-matrix(NA,nrow=nperm,ncol=1)
MS.bootres.ind.ex.rem<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.ind.ex.rem<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs_ind_ex_rem <- do.call(rbind, lapply(split(GSL_totals_ind_ex_rem,GSL_totals_ind_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot_ind_ex_rem <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs_ind_ex_rem)
  
  # Store variances and calculate PST
  variance_components_totalGSLb_ind_ex_rem <- as.data.frame(VarCorr(model_totalGSLboot_ind_ex_rem))
  MS.bootpop.ind.ex.rem[i] <- as.numeric(VarCorr(model_totalGSLboot_ind_ex_rem)[1]) 
  MS.bootres.ind.ex.rem[i] <- as.numeric(VarCorr(model_totalGSLboot_ind_ex_rem)[2])
  tablePSTboot.ind.ex.rem[i] <- MS.bootpop.ind.ex.rem[i] / (MS.bootpop.ind.ex.rem[i] + 2 * MS.bootres.ind.ex.rem[i])
}

# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST_ind_ex_rem <- mean(tablePSTboot.ind.ex.rem, na.rm = TRUE)
se_PST_ind_ex_rem <- sd(tablePSTboot.ind.ex.rem, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.ind.ex.rem)))
pst_ci_ind_ex_rem <- quantile(tablePSTboot.ind.ex.rem, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST_ind_ex_rem, "\n")
cat("Standard Error:", se_PST_ind_ex_rem, "\n")
cat("95% Confidence Interval:", pst_ci_ind_ex_rem, "\n")

#### > Without OH-I3M ----

# log total GSLs 
hist(GSL_totals_ind_rem_ohi3m$logGSL)

model_totalGSL_ind_rem_ohi3m <- lme(logGSL ~ 1, random = ~1|Population, data = GSL_totals_ind_rem_ohi3m)

# PST for dist.r = 1
varpop.obs.ind.rem.ohi3m <- as.numeric(VarCorr(model_totalGSL_ind_rem_ohi3m)[1]) # Between-population variance
varres.obs.ind.rem.ind.ohi3m <- as.numeric(VarCorr(model_totalGSL_ind_rem_ohi3m)[2]) # Within-population variance
pst_ind_rem_ohi3m <- varpop.obs.ind.rem.ohi3m / (varpop.obs.ind.rem.ohi3m + (2 * varres.obs.ind.ex.rem))

# define params set up matrices
nperm<-10000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
GSLboot.ex.ind.rem.ohi3m<- matrix(NA, nrow = nperm, ncol = nrow(GSL_totals_ind_rem_ohi3m))
MS.bootpop.ind.rem.ohi3m <-matrix(NA,nrow=nperm,ncol=1)
MS.bootres.ind.rem.ohi3m <-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.ind.rem.ohi3m<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalGSLs_ind_rem_ohi3m <- do.call(rbind, lapply(split(GSL_totals_ind_rem_ohi3m,GSL_totals_ind_rem_ohi3m$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalGSLboot_ind_rem_ohi3m <- lme(logGSL ~ 1, random = ~1|Population, data = resampled_data_totalGSLs_ind_rem_ohi3m)
  
  # Store variances and calculate PST
  variance_components_totalGSLb_ind_rem_ohi3m <- as.data.frame(VarCorr(model_totalGSLboot_ind_rem_ohi3m))
  MS.bootpop.ind.rem.ohi3m[i] <- as.numeric(VarCorr(model_totalGSLboot_ind_rem_ohi3m)[1]) 
  MS.bootres.ind.rem.ohi3m[i] <- as.numeric(VarCorr(model_totalGSLboot_ind_rem_ohi3m)[2])
  tablePSTboot.ind.rem.ohi3m[i] <- MS.bootpop.ind.rem.ohi3m[i] / (MS.bootpop.ind.rem.ohi3m[i] + 2 * MS.bootres.ind.rem.ohi3m[i])
}

# Calculate mean, standard error, and confidence intervals for bootstrapped PST values
mean_PST_ind_rem_ohi3m <- mean(tablePSTboot.ind.rem.ohi3m, na.rm = TRUE)
se_PST_ind_rem_ohi3m <- sd(tablePSTboot.ind.rem.ohi3m, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.ind.rem.ohi3m)))
pst_ci_ind_rem_ohi3m <- quantile(tablePSTboot.ind.rem.ohi3m, c(0.025, 0.5, 0.975))

cat("Bootstrap Results:\n")
cat("Mean PST:", mean_PST_ind_rem_ohi3m, "\n")
cat("Standard Error:", se_PST_ind_rem_ohi3m, "\n")
cat("95% Confidence Interval:", pst_ind_ci_rem_ohi3m, "\n")


#### total indoles  ----

#### > load data ----

# Read and preprocess data
indole_totals_ind <- read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)

#### > with extreme values ----

# Plot histogram of logindoles
hist(indole_totals_ind$logindoles)

# Boxplot to visualize logindoles by population
ggplot(aes(x = Population, y = logindoles), data = indole_totals_ind) +
  geom_boxplot()

# Fit mixed model for total indoles
model_totalindoles_ind <- lme(logindoles ~ 1, random = ~1 | Population, data = indole_totals_ind)

# Calculate observed variances
varpop.obs.indole.ind <- as.numeric(VarCorr(model_totalindoles_ind)[1])  # Between-population variance
varres.obs.indole.ind <- as.numeric(VarCorr(model_totalindoles_ind)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_indole_ind <- varpop.obs.indole.ind / (varpop.obs.indole.ind + 2 * varres.obs.indole.ind)

# Print PST for verification
print(paste("Observed PST:", pst_observed_indole_ind))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
indole_bootpop_ind <- matrix(NA, nrow = nperm, ncol = 1)
indole_bootres_ind <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_indole_ind <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_indole_ind <- do.call(rbind, lapply(split(indole_totals_ind, indole_totals_ind$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalindoles_boot_ind <- lme(logindoles ~ 1, random = ~1 | Population, data = resampled_data_indole_ind)
  
  # Store variances and calculate PST for the bootstrap sample
  indole_bootpop_ind[i] <- as.numeric(VarCorr(model_totalindoles_boot_ind)[1])
  indole_bootres_ind[i] <- as.numeric(VarCorr(model_totalindoles_boot_ind)[2])
  table_pst_boot_indole_ind[i] <- indole_bootpop_ind[i] / (indole_bootpop_ind[i] + 2 * indole_bootres_ind[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_indole_pst_ind <- mean(table_pst_boot_indole_ind, na.rm = TRUE)
se_indole_pst_ind <- sd(table_pst_boot_indole_ind, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_indole_ind)))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_indole_pst_ind))
print(paste("Standard Error Bootstrapped PST:", se_indole_pst_ind))
pst_ci_indole <- quantile(table_pst_boot_indole_ind, c(0.025, 0.5, 0.975), na.rm = "TRUE")

# Calculate confidence intervals for bootstrapped variances
quantpop_indole_ind <- quantile(indole_bootpop_ind, c(0.025, 0.975))
quantres_indole_ind <- quantile(indole_bootres_ind, c(0.025, 0.975))

indole_bootpop_min_ind <- as.numeric(quantpop_indole_ind[1])
indole_bootpop_max_ind <- as.numeric(quantpop_indole_ind[2])

indole_bootres_min_ind <- as.numeric(quantres_indole_ind[1])
indole_bootres_max_ind <- as.numeric(quantres_indole_ind[2])

# Calculate min and max PST for confidence intervals
pst_boot_min_indole_ind <- indole_bootpop_min_ind / (indole_bootpop_min_ind + 2 * indole_bootres_max_ind)
pst_boot_max_indole_ind <- indole_bootpop_max_ind / (indole_bootpop_max_ind + 2 * indole_bootres_min_ind)

# Plot PST as a function of r
pst_obs_indole_ind <- (dist_r * varpop.obs.indole.ind) / ((dist_r * varpop.obs.indole.ind) + 2 * varres.obs.indole.ind)
pst_obs_ordered_indole_ind <- sort(pst_obs_indole_ind)

plot(
  dist_r, pst_obs_ordered_indole_ind, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 Populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed"), lwd = c(3, 2, 2),
  legend = c("PST", "95% CI (Bootstrap)", "Upper Bound of 95% CI (FST)"),
  col = c("black", "grey", "black")
)

# PST min and max as functions of r
PSTmin_indole_ind <- (dist.r * min(indole_bootpop_ind)) / ((dist.r * min(indole_bootpop_ind)) + 2 * min(indole_bootres_ind))
PSTmin.ord_indole_ind <- sort(PSTmin_indole_ind)
lines(dist.r, PSTmin.ord_indole_ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax_indole_ind <- (dist.r * max(indole_bootpop_ind)) / ((dist.r * max(indole_bootpop_ind)) + 2 * max(indole_bootres_ind))
PSTmax.ord_indole_ind <- sort(PSTmax_indole_ind)
lines(dist.r, PSTmax.ord_indole_ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

# Add horizontal lines
abline(h = 0.03, lty = "dashed", lwd = 2, col = "black")

# Histogram plot of bootstrapped PST values with additional lines and labels
pst_ci_indole_ind <- quantile(table_pst_boot_indole_ind, c(0.025, 0.5, 0.975))
hist(table_pst_boot_indole_ind, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci_indole_ind, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = pst_observed_indole_ind, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### > without extreme values ----

# Plot histogram of logindoles
hist(indole_totals_ind_ex_rem$logindoles)

# Boxplot to visualize logindoles by population
ggplot(aes(x = Population, y = logindoles), data = indole_totals_ex_rem) +
  geom_boxplot()

# Fit mixed model for total indoles
model_totalindoles_ex_rem <- lme(logindoles ~ 1, random = ~1 | Population, data = indole_totals_ex_rem)

# Calculate observed variances
varpop.obs.indole.ex.rem <- as.numeric(VarCorr(model_totalindoles_ex_rem)[1])  # Between-population variance
varres.obs.indole.ex.rem <- as.numeric(VarCorr(model_totalindoles_ex_rem)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_indole_ex_rem <- varpop.obs.indole.ex.rem / (varpop.obs.indole.ex.rem + 2 * varres.obs.indole.ex.rem)

# Print PST for verification
print(paste("Observed PST:", pst_observed_indole_ex_rem))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
indole_bootpop_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
indole_bootres_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_indole_ex_rem <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_indole_ex_rem <- do.call(rbind, lapply(split(indole_totals_ex_rem, indole_totals_ex_rem$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalindoles_boot_ex_rem <- lme(logindoles ~ 1, random = ~1 | Population, data = resampled_data_indole_ex_rem)
  
  # Store variances and calculate PST for the bootstrap sample
  indole_bootpop_ex_rem[i] <- as.numeric(VarCorr(model_totalindoles_boot_ex_rem)[1])
  indole_bootres_ex_rem[i] <- as.numeric(VarCorr(model_totalindoles_boot_ex_rem)[2])
  table_pst_boot_indole_ex_rem[i] <- indole_bootpop_ex_rem[i] / (indole_bootpop_ex_rem[i] + 2 * indole_bootres_ex_rem[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_indole_pst_ex_rem <- mean(table_pst_boot_indole_ex_rem, na.rm = TRUE)
se_indole_pst_ex_rem <- sd(table_pst_boot_indole_ex_rem, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_indole_ex_rem)))
pst_ci_indole_ex_rem <- quantile(table_pst_boot_indole_ex_rem, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_indole_pst))
print(paste("Standard Error Bootstrapped PST:", se_indole_pst))
cat("95% Confidence Interval:", pst_ci, "\n")

#### > without OH-I3M ----

# Plot histogram of logindoles
hist(indole_totals_rem_ohi3m$logindoles)
# compare to regular plot
hist(indole_totals$logindoles)


# Boxplot to visualize logindoles by population
ggplot(aes(x = Population, y = logindoles), data = indole_totals_rem_ohi3m) +
  geom_boxplot()
# compare to regular plot
ggplot(aes(x = Population, y = logindoles), data = indole_totals) +
  geom_boxplot()

# Fit mixed model for total indoles
model_totalindoles_rem_ohi3m <- lme(logindoles ~ 1, random = ~1 | Population, data = indole_totals_rem_ohi3m)

# Calculate observed variances
varpop.obs.indole.rem.ohi3m <- as.numeric(VarCorr(model_totalindoles_rem_ohi3m)[1])  # Between-population variance
varres.obs.indole.rem.ohi3m <- as.numeric(VarCorr(model_totalindoles_rem_ohi3m)[2])  # Within-population variance

# Calculate observed PST for dist.r = 1
pst_observed_indole_rem_ohi3m <- varpop.obs.indole.rem.ohi3m / (varpop.obs.indole.rem.ohi3m + 2 * varres.obs.indole.rem.ohi3m)

# Print PST for verification
print(paste("Observed PST:", pst_observed_indole_rem_ohi3m))

# Set bootstrap parameters
nperm <- 1000
dist_r <- seq(0, 1, by = 0.0001)  # Range for PST calculation

# Initialize matrices for bootstrapping results
indole_bootpop_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)
indole_bootres_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)
table_pst_boot_indole_rem_ohi3m <- matrix(NA, nrow = nperm, ncol = 1)

# Bootstrap loop
for (i in 1:nperm) {
  # Resample data within each population
  resampled_data_indole_rem_ohi3m <- do.call(rbind, lapply(split(indole_totals_rem_ohi3m, indole_totals_rem_ohi3m$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalindoles_boot_rem_ohi3m <- lme(logindoles ~ 1, random = ~1 | Population, data = resampled_data_indole_rem_ohi3m)
  
  # Store variances and calculate PST for the bootstrap sample
  indole_bootpop_rem_ohi3m[i] <- as.numeric(VarCorr(model_totalindoles_boot_rem_ohi3m)[1])
  indole_bootres_rem_ohi3m[i] <- as.numeric(VarCorr(model_totalindoles_boot_rem_ohi3m)[2])
  table_pst_boot_indole_rem_ohi3m[i] <- indole_bootpop_rem_ohi3m[i] / (indole_bootpop_rem_ohi3m[i] + 2 * indole_bootres_rem_ohi3m[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_indole_pst_rem_ohi3m <- mean(table_pst_boot_indole_rem_ohi3m, na.rm = TRUE)
se_indole_pst_rem_ohi3m <- sd(table_pst_boot_indole_rem_ohi3m, na.rm = TRUE) / sqrt(length(na.omit(table_pst_boot_indole_rem_ohi3m)))
pst_ci_indole_rem_ohi3m <- quantile(table_pst_boot_indole_rem_ohi3m, c(0.025, 0.5, 0.975))

# Print bootstrapped PST mean and standard error
print(paste("Mean Bootstrapped PST:", mean_indole_pst))
print(paste("Standard Error Bootstrapped PST:", se_indole_pst))
cat("95% Confidence Interval:", pst_ci, "\n")


#### aliphatics ####

#total aliphatic
aliphatic_totals_ind <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logaliphatics")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# Fit the random effects model

# log total GSLs
hist(aliphatic_totals_ind$logaliphatics)

model_totalaliphatics_ind <- lme(logaliphatics ~ 1, random = ~1|Population, data = aliphatic_totals_ind)

# PST for dist.r = 1
varpop.obs.aliphatic.ind <- as.numeric(VarCorr(model_totalaliphatics_ind)[1]) # Between-population variance
varres.obs.aliphatic.ind <- as.numeric(VarCorr(model_totalaliphatics_ind)[2]) # Within-population variance
aliphatic.pst.ind <- varpop.obs.aliphatic.ind / (varpop.obs.aliphatic.ind + (2 * varres.obs.aliphatic.ind))

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
aliphaticboot.ind <- matrix(NA, nrow = nperm, ncol = nrow(aliphatic_totals_ind))
aliphatic.bootpop.ind <-matrix(NA,nrow=nperm,ncol=1)
aliphatic.bootres.ind <-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.aliphatic.ind <-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalaliphatic_ind <- do.call(rbind, lapply(split(aliphatic_totals_ind,aliphatic_totals_ind$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalaliphatic.boot.ind <- lme(logaliphatics ~ 1, random = ~1|Population, data = resampled_data_totalaliphatic_ind)
  
  # Store variances and calculate PST
  aliphatic.bootpop.ind[i] <- as.numeric(VarCorr(model_totalaliphatic.boot.ind)[1]) 
  aliphatic.bootres.ind[i] <- as.numeric(VarCorr(model_totalaliphatic.boot.ind)[2])
  tablePSTboot.aliphatic.ind[i] <- aliphatic.bootpop.ind[i] / (aliphatic.bootpop.ind[i] + 2 * aliphatic.bootres.ind[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_aliphatic_PST_ind <- mean(tablePSTboot.aliphatic.ind, na.rm = TRUE)
se_aliphatic_PST_ind <- sd(tablePSTboot.aliphatic.ind, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.aliphatic.ind)))

# Calculate confidence intervals
quantpop.aliphatic.ind <- quantile(aliphatic.bootpop.ind, c(0.025, 0.975))
aliphatic.bootpop.min.ind <- as.numeric(quantpop.aliphatic.ind[1])
aliphatic.bootpop.max.ind <- as.numeric(quantpop.aliphatic.ind[2])

quantres.aliphatic.ind <- quantile(aliphatic.bootres.ind, c(0.025, 0.975))
aliphatic.bootres.min.ind <- as.numeric(quantres.aliphatic.ind[1])
aliphatic.bootres.max.ind <- as.numeric(quantres.aliphatic.ind[2])

# Calculate min and max Pst
PST.boot.min.aliphatic.ind <- aliphatic.bootpop.min.ind / (aliphatic.bootpop.min.ind + 2 * aliphatic.bootres.min.ind)
PST.boot.max.aliphatic.ind <- aliphatic.bootpop.max.ind / (aliphatic.bootpop.max.ind + 2 * aliphatic.bootres.max.ind)


# Plot PST as a function of r
PSTobs.aliphatic.ind <- (dist.r * varpop.obs.aliphatic.ind) / ((dist.r * varpop.obs.aliphatic.ind) + (2 * varres.obs.aliphatic.ind))
PSTobs.ord.aliphatic.ind <- sort(PSTobs.aliphatic.ind)
plot(
  dist.r, PSTobs.ord.aliphatic.ind, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.aliphatic.ind <- (dist.r * min(aliphatic.bootpop.ind)) / ((dist.r * min(aliphatic.bootpop.ind)) + 2 * min(aliphatic.bootres.ind))
PSTmin.ord.aliphatic.ind <- sort(PSTmin.aliphatic.ind)
lines(dist.r, PSTmin.ord.aliphatic.ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.aliphatic.ind <- (dist.r * max(aliphatic.bootpop.ind)) / ((dist.r * max(aliphatic.bootpop.ind)) + 2 * max(aliphatic.bootpop.ind))
PSTmax.ord.aliphatic.ind <- sort(PSTmax.aliphatic.ind)
lines(dist.r, PSTmax.ord.aliphatic.ind, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")

# Histogram plot of bootstrapped PST values with additional lines and labels
pst_ci_aliphatic_ind <- quantile(tablePSTboot.aliphatic.ind, c(0.025, 0.5, 0.975))
hist(tablePSTboot.aliphatic.ind, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci_aliphatic_ind, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = aliphatic.pst.ind, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
abline(v = 0.03, col = "black", lwd = 2, lty = "dotted")

#### flavonoids ----

#### > load data ----

#total flavonoid
flavonoid_totals_ind <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logflavonoids")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

flavonoid_totals_ex_rem_ind <- read.csv("./data/mf_means_ex_rem.csv") %>%
  filter(treatment == "CW") %>%  # Filter for control treatment
  select("Population", "mf", "totalflavonoid", "logflavonoids") %>%
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

#### shannon diversity ----

shannon_totals_ind <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "shannon_diversity")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

str(shannon_totals_ind)

#### > with extreme values ----

# Fit the random effects model

model_shannon_ind <- lme(shannon_diversity ~ 1, random = ~1|Population, data = shannon_totals_ind)

# PST for dist.r = 1
varpop.obs.shannon.ind <- as.numeric(VarCorr(model_shannon_ind)[1]) # Between-population variance

varres.obs.shannon.ind <- as.numeric(VarCorr(model_shannon_ind)[2]) # Within-population variance
shannon.pst.ind <- varpop.obs.shannon.ind / (varpop.obs.shannon.ind + (2 * varres.obs.shannon.ind))

print(shannon.pst.ind)

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
shannonboot.ind<- matrix(NA, nrow = nperm, ncol = nrow(shannon_totals))
shannon.bootpop.ind<-matrix(NA,nrow=nperm,ncol=1)
shannon.bootres.ind<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.shannon.ind<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_shannon_ind <- do.call(rbind, lapply(split(shannon_totals_ind, shannon_totals_ind$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_shannon.boot.ind <- lme(shannon_diversity ~ 1, random = ~1|Population, data = resampled_data_shannon_ind)
  
  # Store variances and calculate PST
  shannon.bootpop.ind[i] <- as.numeric(VarCorr(model_shannon.boot.ind)[1]) 
  shannon.bootres.ind[i] <- as.numeric(VarCorr(model_shannon.boot.ind)[2])
  tablePSTboot.shannon.ind[i] <- shannon.bootpop.ind[i] / (shannon.bootpop.ind[i] + 2 * shannon.bootres.ind[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_shannon_PST <- mean(tablePSTboot.shannon, na.rm = TRUE)
se_shannon_PST <- sd(tablePSTboot.shannon, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.flavonoid)))

# Calculate confidence intervals
boot.quant.shannon.ind <- quantile(tablePSTboot.shannon.ind, c(0.025, 0.975))
shannon.bootpop.min.ind <- as.numeric(boot.quant.shannon.ind[1])
shannon.bootpop.max.ind <- as.numeric(boot.quant.shannon.ind[2])

# plot it
varpop.obs.shannon / (varpop.obs.shannon + (2 * varres.obs.shannon))
PSTobs.shannon <- (dist.r * varpop.obs.shannon) / ((dist.r * varpop.obs.shannon) + (2 * varres.obs.shannon))
PSTobs.ord <- sort(PSTobs.shannon)
plot(
  dist.r, PSTobs.shannon, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.shannon <- dist.r*shannon.bootpop.min
PSTmin.ord.shannon <- sort(PSTmin.shannon)
lines(dist.r, PSTmin.ord.shannon, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.shannon <- dist.r * shannon.bootpop.max
PSTmax.ord.shannon <- sort(PSTmax.shannon)
lines(dist.r, PSTmax.ord.shannon, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")


#### > with extreme values ----

# Fit the random effects model

hist(flavonoid_totals_ind$logflavonoids)

model_totalflavonoids_ind <- lme(logflavonoids ~ 1, random = ~1|Population, data = flavonoid_totals_ind)

# PST for dist.r = 1
varpop.obs.flavonoid.ind <- as.numeric(VarCorr(model_totalflavonoids_ind)[1]) # Between-population variance
varres.obs.flavonoid.ind <- as.numeric(VarCorr(model_totalflavonoids_ind)[2]) # Within-population variance
flavonoid.pst.ind <- varpop.obs.flavonoid.ind / (varpop.obs.flavonoid.ind + (2 * varres.obs.flavonoid.ind))

# define params set up matrices
nperm<-1000
dist.r <- seq(0, 1, by = 0.0001)  # see expression for PST
flavonoidboot.ind<- matrix(NA, nrow = nperm, ncol = nrow(flavonoid_totals_ind))
flavonoid.bootpop.ind<-matrix(NA,nrow=nperm,ncol=1)
flavonoid.bootres.ind<-matrix(NA,nrow=nperm,ncol=1)
tablePSTboot.flavonoid.ind<-matrix(NA,nrow=nperm,ncol=1)

# Bootstrap
for (i in 1:nperm) {
  # Resample within each population
  resampled_data_totalflavonoid.ind <- do.call(rbind, lapply(split(flavonoid_totals_ind,flavonoid_totals_ind$Population), function(pop_data) {
    pop_data[sample(nrow(pop_data), nrow(pop_data), replace = TRUE), ]
  }))
  
  # Fit the model to the bootstrap sample
  model_totalflavonoid.boot.ind <- lme(logflavonoids ~ 1, random = ~1|Population, data = resampled_data_totalflavonoid.ind)
  
  # Store variances and calculate PST
  flavonoid.bootpop.ind[i] <- as.numeric(VarCorr(model_totalflavonoid.boot.ind)[1]) 
  flavonoid.bootres.ind[i] <- as.numeric(VarCorr(model_totalflavonoid.boot.ind)[2])
  tablePSTboot.flavonoid.ind[i] <- flavonoid.bootpop.ind[i] / (flavonoid.bootpop.ind[i] + 2 * flavonoid.bootres.ind[i])
}

# Calculate mean and standard error of bootstrapped PST values
mean_flavonoid_PST_ind <- mean(tablePSTboot.flavonoid.ind, na.rm = TRUE)
se_flavonoid_PST_ind <- sd(tablePSTboot.flavonoid.ind, na.rm = TRUE) / sqrt(length(na.omit(tablePSTboot.flavonoid.ind)))

# Calculate confidence intervals
boot.quant.flavonoid.ind <- quantile(tablePSTboot.flavonoid.ind, c(0.025, 0.975))
flavonoid.bootpop.min.ind <- as.numeric(boot.quant.flavonoid.ind[1])
flavonoid.bootpop.max.ind <- as.numeric(boot.quant.flavonoid.ind[2])

# plot it
plot(
  dist.r, PSTobs.ord.flavonoid, type = "l", lty = "solid", lwd = 3, col = "black",
  xlab = "r", ylab = "PST", main = "PST across 13 populations",
  xlim = c(0, 1), ylim = c(0, 1)
)
legend(
  "topleft", lty = c("solid", "dotted", "dashed", "dotted"), lwd = c(3, 2, 2, 2),
  legend = c("PST", "95% CI", "Upper bound of 95% CI of FST"),
  col = c("black", "grey", "black", "black")
)


# PSTmin and PSTmax as functions of r
PSTmin.flavonoid <- dist.r*flavonoid.bootpop.min
PSTmin.ord.flavonoid <- sort(PSTmin.flavonoid)
lines(dist.r, PSTmin.ord.flavonoid, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.flavonoid <- dist.r * flavonoid.bootpop.max
PSTmax.ord.flavonoid <- sort(PSTmax.flavonoid)
lines(dist.r, PSTmax.ord.flavonoid, type = "l", lty = "dotted", lwd = 2, col = "grey")

abline(h = 0.03, lty = "dotted", lwd = 2, col = "black")


### Pairwise Pst ----

#### Total GSLs ----

# Calculate mean and variance within each population
totalGSL_population_stats <- GSL_totals %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(logGSL),
    PopVar = var(logGSL),
    .groups = 'drop'
  )

str(totalGSL_population_stats)

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(totalGSL_population_stats$Population),
                                Pop2 = unique(totalGSL_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

str(population_pairs)

# Initialize a data frame to store pairwise pst results
pairwise_pst_totalGSL <- data.frame(Pop1 = character(), Pop2 = character(), pairwise_pst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate pst
# Loop through each pair to calculate pst

for (i in 1:nrow(population_pairs)) {
  pair_totalGSL <- population_pairs[i, ]
  pop1_stats_totalGSL <- totalGSL_population_stats %>% filter(Population == pair_totalGSL$Pop1)
  pop2_stats_totalGSL <- totalGSL_population_stats %>% filter(Population == pair_totalGSL$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_totalGSL <- mean(c(pop1_stats_totalGSL$PopMean, pop2_stats_totalGSL$PopMean))
  between_pop_var_totalGSL <- sum((c(pop1_stats_totalGSL$PopMean, pop2_stats_totalGSL$PopMean) - mean_of_pop_means_totalGSL)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_totalGSL <- mean(c(pop1_stats_totalGSL$PopVar, pop2_stats_totalGSL$PopVar))
  
  # Calculate pst
  pst_totalGSL <- between_pop_var_totalGSL / (between_pop_var_totalGSL + 2*within_pop_var_totalGSL)
  
  # Store the result
  pairwise_pst_totalGSL <- rbind(pairwise_pst_totalGSL, data.frame(Pop1 = pair_totalGSL$Pop1, Pop2 = pair_totalGSL$Pop2, pairwise_pst = pst_totalGSL))
}

str(pairwise_pst_totalGSL)

# Filter to keep only unique pairs
#unique_pairwise_pst_totalGSL <- pairwise_pst_totalGSL %>%
#  rowwise() %>%
#  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
#  ungroup() %>%
#  distinct(PairID, .keep_all = TRUE) %>%
#  select(-PairID) # Remove the helper column


#### Total indoles ----

#  Calculate mean and variance within each population
indoles_population_stats <- indole_totals %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(logindoles),
    PopVar = var(logindoles),
    .groups = 'drop'
  )

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(indoles_population_stats$Population),
                                Pop2 = unique(indoles_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise pst results
pairwise_pst_indoles <- data.frame(Pop1 = character(), Pop2 = character(), pairwise_pst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate pst
# Loop through each pair to calculate pst
for (i in 1:nrow(population_pairs)) {
  pair_indoles <- population_pairs[i, ]
  pop1_stats_indoles <- indoles_population_stats %>% filter(Population == pair_indoles$Pop1)
  pop2_stats_indoles <- indoles_population_stats %>% filter(Population == pair_indoles$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_indoles <- mean(c(pop1_stats_indoles$PopMean, pop2_stats_indoles$PopMean))
  between_pop_var_indoles <- sum((c(pop1_stats_indoles$PopMean, pop2_stats_indoles$PopMean) - mean_of_pop_means_indoles)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_indoles <- mean(c(pop1_stats_indoles$PopVar, pop2_stats_indoles$PopVar))
  
  # Calculate pst
  pst_indoles <- between_pop_var_indoles / (between_pop_var_indoles + 2*within_pop_var_indoles)
  
  # Store the result
  pairwise_pst_indoles <- rbind(pairwise_pst_indoles, data.frame(Pop1 = pair_indoles$Pop1, Pop2 = pair_indoles$Pop2, pairwise_pst = pst_indoles))
}

str(pairwise_pst_indoles)

#### Total aliphatics ----

#  Calculate mean and variance within each population
aliphatics_population_stats <- aliphatic_totals %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(logaliphatics),
    PopVar = var(logaliphatics),
    .groups = 'drop'
  )

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(aliphatics_population_stats$Population),
                                Pop2 = unique(aliphatics_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise pst results
pairwise_pst_aliphatics <- data.frame(Pop1 = character(), Pop2 = character(), pst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate pst
# Loop through each pair to calculate pst
for (i in 1:nrow(population_pairs)) {
  pair_aliphatics <- population_pairs[i, ]
  pop1_stats_aliphatics <- aliphatics_population_stats %>% filter(Population == pair_aliphatics$Pop1)
  pop2_stats_aliphatics <- aliphatics_population_stats %>% filter(Population == pair_aliphatics$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_aliphatics <- mean(c(pop1_stats_aliphatics$PopMean, pop2_stats_aliphatics$PopMean))
  between_pop_var_aliphatics <- sum((c(pop1_stats_aliphatics$PopMean, pop2_stats_aliphatics$PopMean) - mean_of_pop_means_aliphatics)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_aliphatics <- mean(c(pop1_stats_aliphatics$PopVar, pop2_stats_aliphatics$PopVar))
  
  # Calculate pst
  pst_aliphatics <- between_pop_var_aliphatics / (between_pop_var_aliphatics + 2*within_pop_var_aliphatics)
  
  # Store the result
  pairwise_pst_aliphatics <- rbind(pairwise_pst_aliphatics, data.frame(Pop1 = pair_aliphatics$Pop1, Pop2 = pair_aliphatics$Pop2, pairwise_pst = pst_aliphatics))
}

str(pairwise_pst_aliphatics)

#### Total flavonoids ----

#  Calculate mean and variance within each population
flavonoid_population_stats <- flavonoid_totals %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(logflavonoids),
    PopVar = var(logflavonoids),
    .groups = 'drop'
  )

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(flavonoid_population_stats$Population),
                                Pop2 = unique(flavonoid_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise pst results
pairwise_pst_flavonoids <- data.frame(Pop1 = character(), Pop2 = character(), pairwise_pst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate pst
# Loop through each pair to calculate pst
for (i in 1:nrow(population_pairs)) {
  pair_flavonoids <- population_pairs[i, ]
  pop1_stats_flavonoids <- flavonoid_population_stats %>% filter(Population == pair_flavonoids$Pop1)
  pop2_stats_flavonoids <- flavonoid_population_stats %>% filter(Population == pair_flavonoids$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_flavonoids <- mean(c(pop1_stats_flavonoids$PopMean, pop2_stats_flavonoids$PopMean))
  between_pop_var_flavonoids <- sum((c(pop1_stats_flavonoids$PopMean, pop2_stats_flavonoids$PopMean) - mean_of_pop_means_flavonoids)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_flavonoids <- mean(c(pop1_stats_flavonoids$PopVar, pop2_stats_flavonoids$PopVar))
  
  # Calculate pst
  pst_flavonoids <- between_pop_var_flavonoids / (between_pop_var_flavonoids + 2*within_pop_var_flavonoids)
  
  # Store the result
  pairwise_pst_flavonoids <- rbind(pairwise_pst_flavonoids, data.frame(Pop1 = pair_flavonoids$Pop1, Pop2 = pair_flavonoids$Pop2, pairwise_pst = pst_flavonoids))
}

str(pairwise_pst_flavonoids)

#### Shannon's  diversity ----

#  Calculate mean and variance within each population
shannon_population_stats <- shannon_totals %>%
  group_by(Population) %>%
  summarise(
    PopMean = mean(shannon_diversity),
    PopVar = var(shannon_diversity),
    .groups = 'drop'
  )

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(shannon_population_stats$Population),
                                Pop2 = unique(shannon_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise pst results
pairwise_pst_shannon <- data.frame(Pop1 = character(), Pop2 = character(), pairwise_pst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate pst
# Loop through each pair to calculate pst
for (i in 1:nrow(population_pairs)) {
  pair_shannon <- population_pairs[i, ]
  pop1_stats_shannon <- shannon_population_stats %>% filter(Population == pair_shannon$Pop1)
  pop2_stats_shannon <- shannon_population_stats %>% filter(Population == pair_shannon$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_shannon <- mean(c(pop1_stats_shannon$PopMean, pop2_stats_shannon$PopMean))
  between_pop_var_shannon <- sum((c(pop1_stats_shannon$PopMean, pop2_stats_shannon$PopMean) - mean_of_pop_means_shannon)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_shannon <- mean(c(pop1_stats_shannon$PopVar, pop2_stats_shannon$PopVar))
  
  # Calculate pst
  pst_shannon <- between_pop_var_shannon / (between_pop_var_shannon + 2*within_pop_var_shannon)
  
  # Store the result
  pairwise_pst_shannon <- rbind(pairwise_pst_shannon, data.frame(Pop1 = pair_shannon$Pop1, Pop2 = pair_shannon$Pop2, pairwise_pst = pst_shannon))
}

str(pairwise_pst_shannon)

#### join pairwise pst ----

# total GSLs
combined_data_Pst <- pairwise_pst_totalGSL %>%
  #inner_join(fst, by = c("Pop1", "Pop2")) %>%
  dplyr::select(Pop1, Pop2, total_gsl_pairwise_pst = pairwise_pst) %>%
  left_join(pairwise_pst_indoles, by = c("Pop1", "Pop2"))  %>%
  dplyr::select(Pop1, Pop2, total_gsl_pairwise_pst, indole_pairwise_pst = pairwise_pst) %>%
  left_join(pairwise_pst_aliphatics, by = c("Pop1", "Pop2")) %>%
  dplyr::select(Pop1, Pop2, total_gsl_pairwise_pst, indole_pairwise_pst, aliphatic_pairwise_pst = pairwise_pst) %>%
  left_join(pairwise_pst_flavonoids, by = c("Pop1", "Pop2")) %>%
  dplyr::select(Pop1, Pop2, total_gsl_pairwise_pst, indole_pairwise_pst, aliphatic_pairwise_pst, flavonoid_pairwise_pst = pairwise_pst) %>% 
  left_join(pairwise_pst_shannon, by = c("Pop1", "Pop2")) %>%
  dplyr::select(Pop1, Pop2, total_gsl_pairwise_pst, indole_pairwise_pst, aliphatic_pairwise_pst, flavonoid_pairwise_pst, shannon_pairwise_Pst = pairwise_pst)
  
str(combined_data_Pst)
  
  
#  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-"))


### Join with climate PC data ----

#### load climate data  ----

climate_data <- read.csv(file = "data/pop_means_with_clim.csv") %>%
  select(-X) %>%
  filter(treatment == "C") %>% # filter for the controls
  select(Population, contemporary_PC1, contemporary_PC2, historic_PC1, historic_PC2) %>%
  filter(Population %in% overlap$pop)

climate_data

#### Create all pairwise combinations of populations ----

climate_pop_pairs <- expand.grid(Pop1 = climate_data$Population, Pop2 = climate_data$Population) %>%
  filter(Pop1 != Pop2) 

# Calculate pairwise differences for each climate variable
str(combined_data_Pst)
pairwise_climate_differences <- climate_pop_pairs %>%
  rowwise() %>%
  mutate(
    # contemporary_PC1
    contemporary_PC1_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(contemporary_PC1) -
        climate_data %>% filter(Population == Pop2) %>% pull(contemporary_PC1)
    ),
    # contemporary_PC2
    contemporary_PC2_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(contemporary_PC2) -
        climate_data %>% filter(Population == Pop2) %>% pull(contemporary_PC2)
    ),
    # historic PC1 Diff
    historic_PC1_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(historic_PC1) -
        climate_data %>% filter(Population == Pop2) %>% pull(historic_PC1)
    ),
    # historic PC2 Diff
    historic_PC2_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(historic_PC2) -
        climate_data %>% filter(Population == Pop2) %>% pull(historic_PC2)
    ),
  ) %>%
  ungroup() %>%
  left_join(combined_data_Pst, by = c("Pop1", "Pop2"))

str(pairwise_climate_differences)

### Mantel test ----

# > make mantel dfs ----

# pst

# mantel test 
str(pairwise_climate_differences)

# Get unique population names
populations <- unique(c(pairwise_climate_differences$Pop1, pairwise_climate_differences$Pop2))

# Create an empty matrix with population names as row and column names
gsl_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                      dimnames = list(populations, populations))
aliphatic_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                          dimnames = list(populations, populations))
indole_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                          dimnames = list(populations, populations))
flavonoid_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                          dimnames = list(populations, populations))
shannon_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                          dimnames = list(populations, populations))

# Populate the matrix with pairwise Pst scores

# Defense traits

# > total GSLs
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$total_gsl_pairwise_pst[i]
  
  gsl_dist_matrix[pop1, pop2] <- pst
  gsl_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}
str(gsl_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(gsl_dist_matrix) <- 0

# Print the matrix
print(gsl_dist_matrix)

# > total aliphatics
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$aliphatic_pairwise_pst[i]
  
  aliphatic_dist_matrix[pop1, pop2] <- pst
  aliphatic_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}
str(aliphatic_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(aliphatic_dist_matrix) <- 0

# Print the matrix
print(aliphatic_dist_matrix)

# > total indoles
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$indole_pairwise_pst[i]
  
  indole_dist_matrix[pop1, pop2] <- pst
  indole_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}
str(indole_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(indole_dist_matrix) <- 0

# Print the matrix
print(indole_dist_matrix)

# > total flavonoid
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$flavonoid_pairwise_pst[i]
  
  flavonoid_dist_matrix[pop1, pop2] <- pst
  flavonoid_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}
str(flavonoid_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(flavonoid_dist_matrix) <- 0

# Print the matrix
print(flavonoid_dist_matrix)

# > shannon diversity
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$shannon_pairwise_Pst[i]
  
  shannon_dist_matrix[pop1, pop2] <- pst
  shannon_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}
str(shannon_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(shannon_dist_matrix) <- 0

# Print the matrix
print(aliphatic_dist_matrix)

# Climate

hist_PC1_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                          dimnames = list(populations, populations)) 
hist_PC2_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                               dimnames = list(populations, populations)) 
contemp_PC1_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                               dimnames = list(populations, populations)) 
contemp_PC2_dist_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                               dimnames = list(populations, populations)) 

# historic PC1
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$historic_PC1_Diff[i]
  
  hist_PC1_dist_matrix[pop1, pop2] <- pst
  hist_PC1_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}

str(hist_PC1_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(hist_PC1_dist_matrix) <- 0

# Print the matrix
print(hist_PC1_dist_matrix)

# historic PC2
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$historic_PC2_Diff[i]
  
  hist_PC2_dist_matrix[pop1, pop2] <- pst
  hist_PC2_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}

str(hist_PC2_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(hist_PC2_dist_matrix) <- 0

# Print the matrix
print(hist_PC2_dist_matrix)

# contemporary PC1
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$contemporary_PC1_Diff[i]
  
  contemp_PC1_dist_matrix[pop1, pop2] <- pst
  contemp_PC1_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}

str(contemp_PC1_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(contemp_PC1_dist_matrix) <- 0

# Print the matrix
print(contemp_PC1_dist_matrix)

# contemporary PC2
for (i in 1:nrow(pairwise_climate_differences)) {
  pop1 <- pairwise_climate_differences$Pop1[i]
  pop2 <- pairwise_climate_differences$Pop2[i]
  pst <- pairwise_climate_differences$contemporary_PC2_Diff[i]
  
  contemp_PC2_dist_matrix[pop1, pop2] <- pst
  contemp_PC2_dist_matrix[pop2, pop1] <- pst  # Ensure symmetry
}

str(contemp_PC2_dist_matrix)

# Fill diagonal with zeros (self-comparisons)
diag(contemp_PC2_dist_matrix) <- 0

# Print the matrix
print(contemp_PC2_dist_matrix)

# > run analysis ----

### GSL

#GSL & historic PC1
gsl_mantel_result <- mantel(gsl_dist_matrix, hist_PC1_dist_matrix, method = "pearson", permutations = 999)
print(gsl_mantel_result) # p = 0.37

#GSL & historic PC2
gsl_mantel_result_2 <- mantel(gsl_dist_matrix, hist_PC2_dist_matrix, method = "pearson", permutations = 999)
print(gsl_mantel_result_2) # p = 0.47

#GSL & contemp PC1
gsl_mantel_result_3 <- mantel(gsl_dist_matrix, contemp_PC1_dist_matrix, method = "pearson", permutations = 999)
print(gsl_mantel_result_3) # p = 0.65

#GSL & contemp PC2
gsl_mantel_result_4 <- mantel(gsl_dist_matrix, contemp_PC2_dist_matrix, method = "pearson", permutations = 999)
print(gsl_mantel_result_4) # p = 0.47

### Indoles

#indole & historic PC1
indole_mantel_result_1 <- mantel(indole_dist_matrix, hist_PC1_dist_matrix, method = "pearson", permutations = 999)
print(indole_mantel_result_1)
# p = 0.7

#indole & historic PC2
indole_mantel_result_2 <- mantel(indole_dist_matrix, hist_PC2_dist_matrix, method = "pearson", permutations = 999)
print(indole_mantel_result_2) # p = 0.06

#indole & contemp PC1
indole_mantel_result_3 <- mantel(indole_dist_matrix, contemp_PC1_dist_matrix, method = "pearson", permutations = 999)
print(indole_mantel_result_3) # p = .7

#indole & contemp PC2
indole_mantel_result_4 <- mantel(indole_dist_matrix, contemp_PC2_dist_matrix, method = "pearson", permutations = 999)
print(indole_mantel_result_4) # p = .06                                                                                                                                                                                                                                                                     

### Aliphatics

#aliphatic & historic PC1
aliphatic_mantel_result <- mantel(aliphatic_dist_matrix, hist_PC1_dist_matrix, method = "pearson", permutations = 999)
print(aliphatic_mantel_result) # p  = 0.72

#aliphatic & historic PC2
aliphatic_mantel_result_2 <- mantel(aliphatic_dist_matrix, hist_PC2_dist_matrix, method = "pearson", permutations = 999)
print(aliphatic_mantel_result_2) # p  = 0.39

#aliphatic & contemp PC1
aliphatic_mantel_result_3 <- mantel(aliphatic_dist_matrix, contemp_PC1_dist_matrix, method = "pearson", permutations = 999)
print(aliphatic_mantel_result_3) # p  = 0.65

#aliphatic & contemp PC2
aliphatic_mantel_result_4 <- mantel(aliphatic_dist_matrix, contemp_PC2_dist_matrix, method = "pearson", permutations = 999)
print(aliphatic_mantel_result_4) # p  = 0.38

### Flavonoids

# flavonoid & historic PC1
flavonoid_mantel_result_1 <- mantel(flavonoid_dist_matrix, hist_PC1_dist_matrix, method = "pearson", permutations = 999)
print(flavonoid_mantel_result_1) # p = 0.7

# flavonoid & historic PC2
flavonoid_mantel_result_2 <- mantel(flavonoid_dist_matrix, hist_PC2_dist_matrix, method = "pearson", permutations = 999)
print(flavonoid_mantel_result_2) # p  = 0.36

# flavonoid & contemp PC1
flavonoid_mantel_result_3 <- mantel(flavonoid_dist_matrix, contemp_PC1_dist_matrix, method = "pearson", permutations = 999)
print(flavonoid_mantel_result_3) # p = 0.69

# flavonoid & contemp PC2
flavonoid_mantel_result_4 <- mantel(flavonoid_dist_matrix, contemp_PC2_dist_matrix, method = "pearson", permutations = 999)
print(flavonoid_mantel_result_4) # 0.32

### Shannon 

# shannon & historic PC1
shannon_mantel_result_1 <- mantel(shannon_dist_matrix, hist_PC1_dist_matrix, method = "pearson", permutations = 999)
print(shannon_mantel_result_1) # p = 0.94

# shannon & historic PC2
shannon_mantel_result_2 <- mantel(shannon_dist_matrix, hist_PC2_dist_matrix, method = "pearson", permutations = 999)
print(shannon_mantel_result_2) # p = 0.47

# shannon & contemp PC1
shannon_mantel_result_3 <- mantel(shannon_dist_matrix, contemp_PC1_dist_matrix, method = "pearson", permutations = 999)
print(shannon_mantel_result_3) # p = 0.95

# shannon & contemp PC2
shannon_mantel_result_4 <- mantel(shannon_dist_matrix, contemp_PC2_dist_matrix, method = "pearson", permutations = 999)
print(shannon_mantel_result_4) # p = 0.40

### OLD/ARCHIVED CODE ----

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
rds_data <- rownames_to_column(rds_data, var = "tube_label")

# Merge the data frames based on the "ID" column
merged_data <- merge(rds_data, csv_data, by = "tube_label") %>%
  select(-c("tube_label"))

# filter out pops that are not in my experiment
overlap <- read.csv("data/overlap.csv")
pops_to_keep <- overlap$pop
filtered_merged_data <- merged_data %>%
  filter(pop %in% pops_to_keep)

# move pop and mf to frint of df
filtered_merged_data <- filtered_merged_data %>%
  select(mf, everything()) %>%
  select(pop, everything())

# Calculate Fst

#########################         wc.calc        #########################       
#
# returns the a, b, and c values required to compute Fst according to Weir and Cockerham 1984
#	modified from code previously implemented in hierfstat package 
#	(Jerome Goudet, http://cran.r-project.org/web/packages/hierfstat/index.html)
#
#	ndat		data frame with first column indicating population of origin and following representing loci
#	diploid 	Whether data are diploid

wc.calc <- function(ndat, diploid = TRUE){
  if (!diploid){
    dum <- ndat[, -1]
    nd <- max(dum, na.rm = TRUE)
    modu <- 1000
    if(nd < 10)		modu <- 10
    if(nd < 100)	modu <- 100
    dum <- dum * modu + dum
    ndat <- data.frame(ndat[, 1], dum)
  }
  pop <- ndat[, 1]
  ni <- length(pop)
  dat <- ndat
  loc.names <- names(dat)[-1]
  n <- t(ind.count(dat)) ## NEED in.count to give number of inds genotyped per locus and per pop
  nt <- apply(n, 1, sum, na.rm = TRUE)
  untyped.loc <- which(nt == 0)
  typed.loc <- which(nt != 0)
  if(length(untyped.loc) > 0){
    dat <- dat[, -(untyped.loc + 1)]
    n <- t(ind.count(dat))
    nt <- apply(n, 1, sum, na.rm = TRUE)
  }
  alploc <- nb.alleles(cbind(rep(1, ni), dat[, -1]))
  np <- dim(n)[2]
  npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
  nl <- dim(n)[1]
  p <- pop.freq(dat, diploid)
  pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
  n <- matrix(unlist(n), ncol = np)
  nal <- n[rep(1:nl, alploc), ]
  nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
  ntal <- rep(nt, alploc)
  ncal <- rep(nc, alploc)
  p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
  pb <- matrix(unlist(pb), ncol = 1)
  if(diploid){
    dum <- getal.b(dat[, -1])
    all.loc <- apply(dum, 2, tempfun1 <- function(y) as.numeric(dimnames(table(y))[[1]]))
    hetpl <- apply(dum, 2, fun <- function(z) {
      lapply(as.numeric(dimnames(table(z))[[1]]), who.is.het <- function(y) apply(z == 
                                                                                    y, 1, ind.is.het <- function(x) xor(x[1], x[2])))
    })
    mho <- lapply(hetpl, tempfun2 <- function(x) matrix(unlist(lapply(x, 
                                                                      tempfun3 <- function(y) tapply(y, pop, sum, na.rm = TRUE))), 
                                                        ncol = np))
    mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
    mhom <- (2 * nal * p - mho)/2
  }else{mhom <- nal * p}
  SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
  dum <- nal * (p - 2 * p^2) + mhom
  SSi <- apply(dum, 1, sum, na.rm = TRUE)
  dum1 <- nal * (sweep(p, 1, pb))^2
  SSP <- 2 * apply(dum1, 1, sum, na.rm = TRUE)
  ntalb <- rep(npl, alploc)
  MSG <- SSG/ntal
  MSP <- SSP/(ntalb - 1)
  MSI <- SSi/(ntal - ntalb)
  sigw <- MSG
  sigb <- 0.5 * (MSI - MSG)
  siga <- 1/2/ncal * (MSP - MSI)
  
  abc.mat <- cbind(siga, sigb, sigw)
  # this returns a matrix of a, b, and c values in columns with one allele per row
  return(abc.mat) 
}


#########################       fst.sample      #########################       
#
# returns the Fst value computed on a sample of the alleles in the data
# @param:
#  - obs:  a table containing the components of variance for each locus
#          the table must have one line per allele and at least 3 columns corresponding
#          to the three coefficient a, b, and c as defined in Weir&Cockerham 1984
#
#  - nalleles: the size of the sample (i.e. num of alleles to draw from the table)

fst.sample <- function(obs, nalleles) {
  # choose the alleles to randomly sample:
  allele.smpl <- sample(1: nalleles,size= nalleles,replace=TRUE)
  #select the sampled alleles from the input table:
  dat <- obs[allele.smpl,]
  # Fst = a/(a+b+c); from Weir & Cockerham 1984
  return( sum(dat[,1])/sum(dat[,1]+dat[,2]+dat[,3]) )  
}

#########################       read.fst.input      #########################
#
#function takes from user: filename, number of pops, the number of any 
#	extra columns in front (default is 1), and if there is a header row (default is yes)
#data must be in the form: n columns, columns of all q_hat values by pop and in order, columns 
#	of all q_hat variances by pop in order, any additional columns all saved in a .csv

read.fst.input <- function(q_hat.dat, num.pops, num.extra.columns=1){
  num.pops <- num.pops
  # which columns are the q_hats per population
  q.columns <- num.extra.columns + (1:num.pops)
  # which columns are the corresponding variances in q_hat per population
  last.q.column <- tail(q.columns, n=1)
  first.var.column <- last.q.column+1
  var.columns <- first.var.column:(first.var.column+num.pops-1)	
  
  # extract the data from input file
  dat <- q_hat.dat
  q_hat.matrix <- dat[,q.columns[1]:tail(q.columns, n=1)]
  var.q_hat.matrix <- dat[,var.columns[1]:tail(var.columns, n=1)]
  
  return(list(
    q_hat.matrix, 
    var.q_hat.matrix))
  # this returns a list containing the 2 matrices
  # assign results of the function to an object, then query [1] or [2]
  # 		for the q_hat.matrix and variance.matrix respectively
}



#########################       mean.fst       #########################   
#
# calculate the mean fst from q_hat values for one sample
# this function will then be used in the fst.sample function
# necessary calculations come from Lynch and Milligan 1994, as implemented in AFLP-SURV (Vekemans et al. 2002)
#
#  - dat: a MATRIX of q_hat values per population in columns, per locus in rows (can be produced with AFLP-SURV)
#
#  - vardat: a MATRIX of the variance in q_hat, corresponding to the values in dat (can be produced with AFLP-SURV)
#
#  - num.loci: number of loci individuals are genotyped at
#
#  - num.pops: number of populations sampled

mean.fst <- function(dat, vardat, num.loci, num.pops){
  
  ##### calculate H_j(i)
  H_j_i.matrix <- matrix(NA, nrow=num.loci, ncol=num.pops)
  #	H_j(i) = 2q_j(i)[1-q_j(i)] + 2Var[q_j(i)]
  # j is one population, i is one locus
  for(j in 1:num.pops){ # loop through population columns
    for(i in 1:num.loci){ # loop through loci rows
      q <- dat[i,j]	# q value of one locus at one population
      #x <- xdat[i,j]	# x value of one locus at one population
      #N <- samps[i,j]	# sample size for one locus at one population
      var.q <- vardat[i,j] # variance in q for that locus across all populations
      H_j_i <- 2*q*(1-q) + 2*var.q # measure of gene diversity at one locus
      H_j_i.matrix[i,j] <- H_j_i # put it into a matrix			
    }	
  }
  
  ##### calculate H_j
  # H_j = average of H_j(i) across all loci for that population
  H_j.matrix <- matrix(NA, nrow=1, ncol=num.pops)
  for(j in 1:num.pops){ # get the H_j for each pop
    H_j.temp <- mean(H_j_i.matrix[,j]) # sum over all loci for one pop and divide by number of loci
    # done by taking mean per pop
    H_j.matrix[1,j] <- H_j.temp
  }	
  H_j <- H_j.matrix # this is the mean observed gene diversity in each pop, j
  
  ##### calculate H_jk
  H_jk <- matrix(0, nrow=num.pops, ncol=num.pops) 
  # make a matrix of num.pops x num.pops (j by k)
  for(j in 1:num.pops){
    for(k in 1:num.pops){
      sum.H_jki <- 0
      for(i in 1:num.loci){
        # first need H'_jk(i) where j and k are a pair of populations,  eqn 9a
        Hprime_jki <- dat[i,j] + dat[i,k] - 2*dat[i,j]*dat[i,k]
        # then calculate H_jk(i)  (i.e. not 'prime'),  eqn 10a
        H_jki <- Hprime_jki - (H_j_i.matrix[i,j] + H_j_i.matrix[i,k])/2
        # sum this across all loci in this loop
        sum.H_jki <- sum.H_jki + H_jki
      }
      # then divide by number of loci to get the H_jk for that pair of pops (eqn 12), and store it in a matrix
      H_jk[j,k] <- sum.H_jki/num.loci
    }
  }
  
  ##### calculate H_B, mean between-pop gene diversity
  # H_B = 2/(n(n-1)) * sum(H_jk)  where n is the number of populations  
  first.half.distinct.list <- NULL
  second.half.distinct.list <- NULL
  H_jk.distinct <- NULL
  for(a in 1:(num.pops-1)){
    b <- a+1
    pair <- expand.grid(a, b:num.pops)
    list5 <- pair[,1]
    list6 <- pair[,2]
    first.half.distinct.list <- c(first.half.distinct.list, list5)
    second.half.distinct.list <- c(second.half.distinct.list, list6)
  } # this made 2 lists that when lined up are all possible distinct pairs of populations
  # extract all the H_jk values of distinct population pairs using these lists
  for(z in 1:length(first.half.distinct.list)){
    H_jk.distinct.temp <- H_jk[first.half.distinct.list[z], second.half.distinct.list[z]]
    # print(c(first.half.distinct.list[z], second.half.distinct.list[z])) # check only getting distinct pop pairs, yes.
    H_jk.distinct <- c(H_jk.distinct, H_jk.distinct.temp)
  } # now do eqn 13a
  H_B <- (2 / (num.pops*(num.pops - 1))) * sum(H_jk.distinct)
  
  ##### calculate H_W, mean within-pop diversity
  # H_W = 1/n * sum(H_j)  where n is the number of populations
  H_W <- (1/num.pops)*sum(H_j)
  
  ##### calculate H_T
  # H_T = H_B + H_W
  H_T <- H_B + H_W
  
  ##### calculate Var(H_W)
  var.H_W <- (1 / (num.pops*(num.pops - 1))) * sum( (H_j - H_W)^2 )
  
  ##### calculate V_B
  # modified from Xavier Vekemans C code in AFLP-SURV, found in the Isomain.c file
  # V_B is the variance among diversity measures in non-overlapping pairs of populations - the variance in H_jk
  count <- 0
  sum.H_jk <- 0
  for(k in 1:(num.pops/2)*2){ # a sequence to a decimal only goes to the rounded down decimal, so if odd num.pops still works
    j <- k-1
    temp.H_jk <- H_jk[j,k]
    # sum the non-overlapping pop pairs H_jk's
    sum.H_jk <- sum.H_jk + temp.H_jk
    count <- count+1
    # print(c(j,k))
  }
  mean.H_jk <- sum.H_jk / count # mean non-overlapping pop H_jk
  sum.diffs.squared <- 0
  for(k in 1:(num.pops/2)*2){ # a sequence to a decimal only goes to the rounded down decimal, so if odd num.pops still works
    j <- k-1
    temp.diff.squared <- (H_jk[j,k] - mean.H_jk)^2
    # sum the squared differences to get the variance
    sum.diffs.squared <- sum.diffs.squared + temp.diff.squared
  }
  V_B <- sum.diffs.squared / ((count) * (count-1)) 
  if(num.pops <=3){V_B <- 0} # for 2 and 3 deme cases, V_B and C_B are set to zero 
  
  ##### calculate C_B
  # modified from Xavier Vekemans C code in AFLP-SURV, found in the Isomain.c file
  # C_B is the covariance among diversity measures in overlapping population pairs - the covariance in H_jk
  # covariance is sum over all x minus x_bar times y minus y_bar all divided by n-1
  # overlapping pops are those that share either of the pops in the pair
  count <- 0
  sum.H_jk <- 0
  sum.H_j2k2 <- 0
  for(j in 1:(num.pops-1)){ # loop through pop 1 of X pair
    for(k in (j+1):num.pops){ # loop through pop 2 of X pair
      # get first set of Y pair of populations, those which share pop 1 of X
      j2 <- j # this is first pop of Y pair, shared with first pop of X pair
      k2 <- k+1 # this is second pair of Y pop
      while(k2 <= num.pops){ # loop through first set, second pair of Y pop
        X.H_jk <- H_jk[j, k]
        sum.H_jk <- sum.H_jk + X.H_jk # get the mean of X pairs
        Y.H_jk <- H_jk[j2, k2]
        sum.H_j2k2 <- sum.H_j2k2 + Y.H_jk # get the mean of Y pairs
        # print(c(j,k,"  ",j2,k2)) 
        k2 <- k2+1
        count <- count + 1
      }
      # get second set of Y pair of populations, those which share pop 2 of X
      j2 <- k # this is first pop of Y pair, shared with second pop of X pair
      k2 <- k+1 # this is second pair of Y pop
      while(k2 <= num.pops){ # loop through second set, second pair of Y pop
        X.H_jk <- H_jk[j, k]
        sum.H_jk <- sum.H_jk + X.H_jk
        Y.H_jk <- H_jk[j2, k2]
        sum.H_j2k2 <- sum.H_j2k2 + Y.H_jk
        # print(c(j,k,"  ",j2,k2)) 
        k2 <- k2+1
        count <- count + 1
      }
    }
  }
  mean.H_jk <- sum.H_jk/count
  mean.H_j2k2 <- sum.H_j2k2/count
  cov <- 0
  # do the same as above but now get the covariance, i.e subtracting the mean from each and multiplying together
  for(j in 1:(num.pops-1)){
    for(k in (j+1):num.pops){
      j2 <- j
      k2 <- k+1
      while(k2 <= num.pops){
        cov.temp <- (H_jk[j, k] - mean.H_jk) * (H_jk[j2, k2] - mean.H_j2k2)
        cov <- cov + cov.temp
        k2 <- k2+1
      }
      j2 <- k
      k2 <- k+1
      while(k2 <= num.pops){
        cov.temp <- (H_jk[j, k] - mean.H_jk) * (H_jk[j2, k2] - mean.H_j2k2)
        cov <- cov + cov.temp
        k2 <- k2+1
      }
    }
  }
  C_B <- cov / ((count) * (count-1))
  if(num.pops <=3){C_B <- 0} # for 2 and 3 deme cases, V_B and C_B are set to zero 
  
  ##### calculate Var(H_B)
  var.H_B <- (2 * (V_B + 2*(num.pops-2)*C_B)) / (num.pops * (num.pops - 1))
  
  ##### calculate Cov(H_B, H_W)
  # need the sum of the sums term in there - sum of all H_jk per population j (except where j=k) times H_j
  # multiply the sum of each row of the H_jk matrix by its corresponding H_j from the H_j.matrix excluding pops where j=k in that sum
  sum.H_j.H_jk <- 0
  for(j in 1:num.pops){
    # sum the row of the H_jk matrix, then subtract off the H_jk where j=k
    temp.sum.H_jk <- 0
    for(k in 1:num.pops){
      if(j==k) same.jk <- H_jk[j,k]  #this will be the value to subtract off the sum of H_jk for that row
      temp.sum.H_jk <- temp.sum.H_jk + H_jk[j,k]
      # print(c(j,k))
    }
    sum.H_jk <- temp.sum.H_jk - same.jk
    temp.H_j.H_jk <- H_j[j] * sum.H_jk
    sum.H_j.H_jk <- sum.H_j.H_jk + temp.H_j.H_jk
  }
  cov.H_B.H_W <- ( (sum.H_j.H_jk/(num.pops*(num.pops-1)))-(H_W*H_B) )/num.pops
  
  ##### Fst!
  mean.Fst <- (H_B / H_T) * (1 + (H_B*var.H_W - H_W*var.H_B + (H_B - H_W)*cov.H_B.H_W) / (H_B*H_T^2 ))^(-1)
  
  return(mean.Fst)
}

### Run Fst

fst.dat <- filtered_merged_data

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
  select(pop, everything())

fst.dat[] <- lapply(fst.dat, as.numeric)
fst.dat <- na.omit(fst.dat)

nloci <- ncol(fst.dat)-1
abc.mat <- wc.calc(fst.dat, diploid=TRUE)
nalleles <- nrow(abc.mat)

#the observed Fst:
fst.obs <- sum(abc.mat[,1])/sum(abc.mat[,1] + abc.mat[,2] + abc.mat[,3])

### Pairwise Fst ----

# Create a matrix to store pairwise Fst
num_pops <- length(unique(fst.dat$pop))
pairwise_fst <- matrix(0, nrow = num_pops, ncol = num_pops)

# Get population names
pop_names <- unique(fst.dat$pop)

# Calculate pairwise Fst
for (i in 1:(num_pops - 1)) {
  for (j in (i + 1):num_pops) {
    pop1 <- pop_names[i]
    pop2 <- pop_names[j]
    
    # Subset data for the two populations
    data_subset <- fst.dat[fst.dat$pop %in% c(pop1, pop2), ]
    
    # Apply wc.calc to subset for these populations only
    wc_values <- wc.calc(data_subset, diploid = TRUE)
    
    # Calculate Fst for this pair using the Fst formula
    pairwise_fst[i, j] <- sum(wc_values[, 1]) / sum(wc_values[, 1] + wc_values[, 2] + wc_values[, 3])
  }
}

# Convert matrix to data frame for better readability
pairwise_fst_df <- as.data.frame(pairwise_fst)
rownames(pairwise_fst_df) <- colnames(pairwise_fst_df) <- pop_names

# Print the pairwise Fst matrix
print(pairwise_fst_df)

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
  select(pop, everything())

#### Trait Data (Qst) ----

#load trait data

#total gsl - danielle 11-19-24 test
GSL_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

#total aliphatic
aliphatic_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "totalaliphatic")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

#total indole
indole_totals <-  read.csv("./data/dw.csv") %>%
  filter(treatment == "C") %>% # filter for the controls
  select("Population", "mf", "logtotalGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3", "YO10"))

# replace pops with numbers using key

# fix rownames

# total GSLs
colnames(GSL_totals) = c("pop", "mf", "totalGSL")

# Aliphatics
colnames(aliphatic_totals) = c("pop", "mf", "totalaliphatics")

# totalGSL 
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



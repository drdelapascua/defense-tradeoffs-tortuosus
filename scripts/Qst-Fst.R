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
library(ggdoctheme) # Custom package https://github.com/rishavray/ggdoctheme
library(ggthemes)
library(ggtext)
library(lme4)
library(ggrepel)


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
abline(h = 0.03, lty = "dashed", lwd = 2, col = "black")

# Plot histogram of bootstrapped PST distribution
hist(tablePSTboot, breaks = 30, col = "lightblue", main = "Bootstrap PST Distribution",
     xlab = "PST", border = "black")
abline(v = pst_ci, col = c("red", "green", "red"), lwd = 2, lty = c("dashed", "solid", "dashed"))
abline(v = pst, col = "blue", lwd = 2, lty = "solid")
legend("topright", legend = c("95% CI Bounds", "Observed PST", "Bootstrapped Mean PST", "Observed FST"), col = c("red", "blue", "green", "black"), lty = c("dashed", "solid", "solid", "dotted"))
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
PSTmin_indole <- (dist.r * min(indole_bootpop)) / ((dist.r * min(indole_bootpop)) + 2 * min(indole_bootres))
PSTmin.ord_indole <- sort(PSTmin_indole)
lines(dist.r, PSTmin.ord_indole, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax_indole <- (dist.r * max(indole_bootpop)) / ((dist.r * max(indole_bootpop)) + 2 * max(indole_bootres))
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
PSTmin.aliphatic <- (dist.r * min(aliphatic.bootpop)) / ((dist.r * min(aliphatic.bootpop)) + 2 * min(aliphatic.bootres))
PSTmin.ord.aliphatic <- sort(PSTmin.aliphatic)
lines(dist.r, PSTmin.ord.aliphatic, type = "l", lty = "dotted", lwd = 2, col = "grey")

PSTmax.aliphatic <- (dist.r * max(aliphatic.bootpop)) / ((dist.r * max(aliphatic.bootpop)) + 2 * max(aliphatic.bootpop))
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


# Rinse & repeat for induced defense profile

# Induced defenses ----

#### Total GSLs ----

# load data
GSL_totals_ind <-  read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>% # filter for the controls
  select("Population", "mf", "logGSL")%>% 
  filter(Population %in% c("BH", "IH", "TM2", "CP2", "DPR", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3"))

# Fit the random effects model

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

#### total indoles #### ----

# Read and preprocess data
indole_totals_ind <- read.csv("./data/mf_means.csv") %>%
  filter(treatment == "CW") %>%  # Filter for control treatment
  select("Population", "mf", "totalindole", "logindoles") %>%
  filter(Population %in% c("BH", "IH", "TM2", "KC2", "LV1", "LV2", "SQ1", "WL1", "WL2", "WL3")) %>%
  filter(logindoles != -Inf)

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
lines(dist.r, PSTmin.ord_indole_min, type = "l", lty = "dotted", lwd = 2, col = "grey")

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

# Get pairwise combinations of populations
population_pairs <- expand.grid(Pop1 = unique(totalGSL_population_stats$Population),
                                Pop2 = unique(totalGSL_population_stats$Population)) %>%
  filter(Pop1 != Pop2)

# Initialize a data frame to store pairwise Qst results
pairwise_qst_totalGSL <- data.frame(Pop1 = character(), Pop2 = character(), Qst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate Qst
# Loop through each pair to calculate Qst
for (i in 1:nrow(population_pairs)) {
  pair_totalGSL <- population_pairs[i, ]
  pop1_stats_totalGSL <- totalGSL_population_stats %>% filter(Population == pair_totalGSL$Pop1)
  pop2_stats_totalGSL <- totalGSL_population_stats %>% filter(Population == pair_totalGSL$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_totalGSL <- mean(c(pop1_stats_totalGSL$PopMean, pop2_stats_totalGSL$PopMean))
  between_pop_var_totalGSL <- sum((c(pop1_stats_totalGSL$PopMean, pop2_stats_totalGSL$PopMean) - mean_of_pop_means_totalGSL)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_totalGSL <- mean(c(pop1_stats_totalGSL$PopVar, pop2_stats_totalGSL$PopVar))
  
  # Calculate Qst
  qst_totalGSL <- between_pop_var_totalGSL / (between_pop_var_totalGSL + 2*within_pop_var_totalGSL)
  
  # Store the result
  pairwise_qst_totalGSL <- rbind(pairwise_qst_totalGSL, data.frame(Pop1 = pair_totalGSL$Pop1, Pop2 = pair_totalGSL$Pop2, Qst = qst_totalGSL))
}


# Filter to keep only unique pairs
#unique_pairwise_qst_totalGSL <- pairwise_qst_totalGSL %>%
#  rowwise() %>%
#  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
#  ungroup() %>%
#  distinct(PairID, .keep_all = TRUE) %>%
#  select(-PairID) # Remove the helper column


left_join(indole_pairwise_qst, fst, by = c("Pop1", "Pop2"))


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

# Initialize a data frame to store pairwise Qst results
pairwise_qst_indoles <- data.frame(Pop1 = character(), Pop2 = character(), Qst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate Qst
# Loop through each pair to calculate Qst
for (i in 1:nrow(population_pairs)) {
  pair_indoles <- population_pairs[i, ]
  pop1_stats_indoles <- indoles_population_stats %>% filter(Population == pair_indoles$Pop1)
  pop2_stats_indoles <- indoles_population_stats %>% filter(Population == pair_indoles$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_indoles <- mean(c(pop1_stats_indoles$PopMean, pop2_stats_indoles$PopMean))
  between_pop_var_indoles <- sum((c(pop1_stats_indoles$PopMean, pop2_stats_indoles$PopMean) - mean_of_pop_means_indoles)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_indoles <- mean(c(pop1_stats_indoles$PopVar, pop2_stats_indoles$PopVar))
  
  # Calculate Qst
  qst_indoles <- between_pop_var_indoles / (between_pop_var_indoles + 2*within_pop_var_indoles)
  
  # Store the result
  pairwise_qst_indoles <- rbind(pairwise_qst_indoles, data.frame(Pop1 = pair_indoles$Pop1, Pop2 = pair_indoles$Pop2, Qst = qst_indoles))
}

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

# Initialize a data frame to store pairwise Qst results
pairwise_qst_aliphatics <- data.frame(Pop1 = character(), Pop2 = character(), Qst = numeric(), stringsAsFactors = FALSE)

# Loop through each pair to calculate Qst
# Loop through each pair to calculate Qst
for (i in 1:nrow(population_pairs)) {
  pair_aliphatics <- population_pairs[i, ]
  pop1_stats_aliphatics <- aliphatics_population_stats %>% filter(Population == pair_aliphatics$Pop1)
  pop2_stats_aliphatics <- aliphatics_population_stats %>% filter(Population == pair_aliphatics$Pop2)
  
  # Calculate between-population variance for this pair
  mean_of_pop_means_aliphatics <- mean(c(pop1_stats_aliphatics$PopMean, pop2_stats_aliphatics$PopMean))
  between_pop_var_aliphatics <- sum((c(pop1_stats_aliphatics$PopMean, pop2_stats_aliphatics$PopMean) - mean_of_pop_means_aliphatics)^2) / 1 # 1 degree of freedom for two populations
  
  # Calculate within-population variance
  within_pop_var_aliphatics <- mean(c(pop1_stats_aliphatics$PopVar, pop2_stats_aliphatics$PopVar))
  
  # Calculate Qst
  qst_aliphatics <- between_pop_var_aliphatics / (between_pop_var_aliphatics + 2*within_pop_var_aliphatics)
  
  # Store the result
  pairwise_qst_aliphatics <- rbind(pairwise_qst_aliphatics, data.frame(Pop1 = pair_aliphatics$Pop1, Pop2 = pair_aliphatics$Pop2, Qst = qst_aliphatics))
}

### Genetic data/Fst ----

#### Load Genetic data ----

fst <- read.table("data/fst.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL) %>%
  select(-FST) %>%
  select("Pop1" = row.names, "Pop2" = Pop1, "Fst" = Pop2)

#### Normalize pair order in Qst data  ----

#total GSLs
pairwise_qst_totalGSL <- pairwise_qst_totalGSL %>% 
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

#indoles
pairwise_qst_indoles <- pairwise_qst_indoles %>% 
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

#aliphatics
pairwise_qst_aliphatics <- pairwise_qst_aliphatics %>% 
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

#### Normalize pair order in Fst data ----

fst <- fst %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

#### join pairwise Qst ----

# total GSLs
combined_data_totalGSL <- pairwise_qst_totalGSL %>%
  inner_join(fst, by = "PairID") %>%
  select(-c("Pop1.y", "Pop2.y", "PairID")) %>%
  select("Pop1" = Pop1.x, "Pop2" = Pop2.x, "Qst" = Qst, "Fst" = Fst)

# indoles
combined_data_indoles <- pairwise_qst_indoles %>%
  inner_join(fst, by = "PairID") %>%
  select(-c("Pop1.y", "Pop2.y", "PairID")) %>%
  select("Pop1" = Pop1.x, "Pop2" = Pop2.x, "Qst" = Qst, "Fst" = Fst)

# aliphatics
combined_data_aliphatics <- pairwise_qst_aliphatics %>%
  inner_join(fst, by = "PairID") %>%
  select(-c("Pop1.y", "Pop2.y", "PairID")) %>%
  select("Pop1" = Pop1.x, "Pop2" = Pop2.x, "Qst" = Qst, "Fst" = Fst)

#### Calculate Qst-Fst in all dataframes ----

# total GSLs
combined_data_totalGSL$Qst_Fst <- combined_data_totalGSL$Qst - combined_data_totalGSL$Fst
combined_data_totalGSL$Qst_Fst_abs <- abs(combined_data_totalGSL$Qst - combined_data_totalGSL$Fst)
head(combined_data_totalGSL)
hist(combined_data_totalGSL$Qst_Fst)
hist(combined_data_totalGSL$Qst_Fst_abs)

# indoles
combined_data_indoles$Qst_Fst <- combined_data_indoles$Qst - combined_data_indoles$Fst
combined_data_indoles$Qst_Fst_abs <- abs(combined_data_indoles$Qst - combined_data_indoles$Fst)
head(combined_data_indoles)
hist(combined_data_indoles$Qst_Fst)
hist(combined_data_indoles$Qst_Fst_abs)

# aliphatics
combined_data_aliphatics$Qst_Fst <- combined_data_aliphatics$Qst - combined_data_aliphatics$Fst
combined_data_aliphatics$Qst_Fst_abs <- abs(combined_data_aliphatics$Qst - combined_data_aliphatics$Fst)
head(combined_data_aliphatics)
hist(combined_data_aliphatics$Qst_Fst)
hist(combined_data_aliphatics$Qst_Fst_abs)

### Join with climate PC data ----

#### load climate data  ----
climate_data <- read.csv(file = "data/climate_PC_data.csv") %>%
  select(-X) %>%
  mutate(Population = str_replace(Population, "YOSE10", "YO10")) %>%
  filter(Population %in% overlap$pop)

#### Create all pairwise combinations of populations ----

climate_pop_pairs <- expand.grid(Pop1 = climate_data$Population, Pop2 = climate_data$Population) %>%
  filter(Pop1 != Pop2) 

# Calculate pairwise differences for each climate variable

pairwise_climate_differences <- climate_pop_pairs %>%
  rowwise() %>%
  mutate(
    # climate water deficit
    cwd_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(cwd) -
        climate_data %>% filter(Population == Pop2) %>% pull(cwd)
    ),
    # precipitation
    ppt_mm_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(ppt_mm) -
        climate_data %>% filter(Population == Pop2) %>% pull(ppt_mm)
    ),
    # snowpack
    pck_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(pck) -
        climate_data %>% filter(Population == Pop2) %>% pull(pck)
    ),
    # snowfall
    snw_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(snw) -
        climate_data %>% filter(Population == Pop2) %>% pull(snw)
    ),
    # min annual temp
    tmin_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(tmin) -
        climate_data %>% filter(Population == Pop2) %>% pull(tmin)
    ),
    # max annual temp
    tmax_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(tmax) -
        climate_data %>% filter(Population == Pop2) %>% pull(tmax)
    ),
    # latitude
    latitude_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(latitude) -
        climate_data %>% filter(Population == Pop2) %>% pull(latitude)
    ), 
    # longitude
    longitude_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(longitude) -
        climate_data %>% filter(Population == Pop2) %>% pull(longitude)
    ),
    # elevation
    elevation_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(elevation) -
        climate_data %>% filter(Population == Pop2) %>% pull(elevation)
    ), 
    # PC1
    PC1_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(PC1) -
        climate_data %>% filter(Population == Pop2) %>% pull(PC1)
    ),
    # PC2
    PC2_Diff = abs(
      climate_data %>% filter(Population == Pop1) %>% pull(PC2) -
        climate_data %>% filter(Population == Pop2) %>% pull(PC2)
    )
  ) %>%
  ungroup()

#### Normalize pair order in Fst data ----

pairwise_climate_differences  <- pairwise_climate_differences  %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()


#### Merge with Qst-Fst scripts ----

# total GSLs 

head(combined_data_totalGSL)
combined_data_totalGSL  <- combined_data_totalGSL %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

qst_fst_climate_totalGSL <- combined_data_totalGSL %>%
  inner_join(pairwise_climate_differences, by = "PairID", relationship = "many-to-many") %>%
  select(-c("Pop1.y", "Pop2.y")) %>%
  distinct(PairID, .keep_all = TRUE) %>%
  rename(
    Pop1 = Pop1.x,
    Pop2 = Pop2.x
  )

# indoles

head(combined_data_indoles)
dim(combined_data_indoles)
combined_data_indoles  <- combined_data_indoles %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

qst_fst_climate_indoles <- combined_data_indoles %>%
  inner_join(pairwise_climate_differences, by = "PairID", relationship = "many-to-many") %>%
  select(-c("Pop1.y", "Pop2.y")) %>%
  distinct(PairID, .keep_all = TRUE) %>%
  rename(
    Pop1 = Pop1.x,
    Pop2 = Pop2.x
  )

# aliphatics 

head(combined_data_aliphatics)
dim(combined_data_aliphatics)
combined_data_aliphatics  <- combined_data_aliphatics %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "-")) %>%
  ungroup()

qst_fst_climate_aliphatics <- combined_data_aliphatics %>%
  inner_join(pairwise_climate_differences, by = "PairID", relationship = "many-to-many") %>%
  select(-c("Pop1.y", "Pop2.y")) %>%
  distinct(PairID, .keep_all = TRUE) %>%
  rename(
    Pop1 = Pop1.x,
    Pop2 = Pop2.x
  )

### Models comparing control defense Qst-Fst with PC1 and PC2 ----

#### Total GSLS ----

# PC1
hist(qst_fst_climate_totalGSL$Qst_Fst)
hist(qst_fst_climate_totalGSL$PC1_Diff)
PC1_GSLtotal <- lm(Qst_Fst ~ PC1_Diff, data = qst_fst_climate_totalGSL)

summary(PC1_GSLtotal)  
# slope estimate -0.037, diff than 0
# y-axis is Qst-Fst, smaller/negative values mean Qst is smaller than Fst, or unifying selection (similar genotype campared to phenotype) where larger/positive values are the opposite (genotype is similar but traits have diverged more)
# As PC1 goes down, totalGSL goes up. PC1 associated with prepcip, dryer and hotter values lower wetter cooler values, 
  
# plot the relationship

PC1_GSLtotal_plot <- ggplot(qst_fst_climate_totalGSL, aes(x = PC1_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_label(aes(label = PairID), vjust = -1) + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC1",
    x = "Climate PC1",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()

# PC2
hist(qst_fst_climate_totalGSL$Qst_Fst)
hist(qst_fst_climate_totalGSL$PC2_Diff)
PC2_GSLtotal <- lm(Qst_Fst ~ PC2_Diff, data = qst_fst_climate_totalGSL)

summary(PC2_GSLtotal)  
# similar relationship but not as strong

# Plot the relationship
PC2_GSLtotal_plot <- ggplot(qst_fst_climate_totalGSL, aes(x = PC2_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC2",
    x = "Climate PC2",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()

#### Indoles ----

# PC1
hist(qst_fst_climate_indoles$Qst_Fst)
hist(qst_fst_climate_indoles$PC1_Diff)
PC1_indoles <- lm(Qst_Fst ~ PC1_Diff, data = qst_fst_climate_indoles)

summary(PC1_indoles)  
# slope estimate -0.036, sig diff from 0, higher R2 valu

# plot the relationship

PC1_indole_plot <- ggplot(qst_fst_climate_indoles, aes(x = PC1_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_label(aes(label = PairID), vjust = -1) + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC1",
    x = "Climate PC1 Difference",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()

# PC2
hist(qst_fst_climate_indoles$Qst_Fst)
hist(qst_fst_climate_indoles$PC2_Diff)
PC2_indoles <- lm(Qst_Fst ~ PC2_Diff, data = qst_fst_climate_indoles)

summary(PC2_indoles)  
# similar relationship but not as strong, sign diff tho!

# Plot the relationship
PC2_indoles_plot <- ggplot(qst_fst_climate_indoles, aes(x = PC2_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC2",
    x = "Climate PC2 Difference",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()


#### Aliphatics ----

# PC1
hist(qst_fst_climate_aliphatics$Qst_Fst)
hist(qst_fst_climate_aliphatics$PC1_Diff)
PC1_aliphatics <- lm(Qst_Fst ~ PC1_Diff, data = qst_fst_climate_aliphatics)

summary(PC1_aliphatics)  
# slope estimate -0.036, sig diff from 0, higher R2 valu

# plot the relationship

PC1_aliphatics_plot <- ggplot(qst_fst_climate_aliphatics, aes(x = PC1_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_label(aes(label = PairID), vjust = -1) + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC1",
    x = "Climate PC1 Difference",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()

PC1_aliphatics_plot

# PC2
hist(qst_fst_climate_aliphatics$Qst_Fst)
hist(qst_fst_climate_aliphatics$PC2_Diff)

PC2_aliphatics <- lm(Qst_Fst ~ PC2_Diff, data = qst_fst_climate_aliphatics)

summary(PC2_indoles)  
# similar relationship but not as strong

# Plot the relationship
PC2_aliphatics_plot <- ggplot(qst_fst_climate_aliphatics, aes(x = PC2_Diff, y = Qst_Fst)) +
  geom_point(color = "blue") + 
  geom_smooth(method = "lm", color = "red") + 
  labs(
    title = "Relationship Between Pairwise Qst-Fst and Climate PC2",
    x = "Climate PC2 Difference",
    y = "Pairwise Qst-Fst"
  ) +
  theme_minimal()

### Create correlational matrix ----

#### total GSLs ----

subset_totalGSL <- qst_fst_climate_totalGSL[, c("Qst", "Fst", "Qst_Fst", "Qst_Fst_abs", "cwd_Diff", "ppt_mm_Diff", "pck_Diff", "snw_Diff", "tmin_Diff", "tmax_Diff", "latitude_Diff", "longitude_Diff", "PC1_Diff", "PC2_Diff")]

# Initialize empty matrices to store results
n_vars_totalGSL <- ncol(subset_totalGSL)
cor_matrix_totalGSL <- matrix(NA, nrow = n_vars_totalGSL, ncol = n_vars_totalGSL, dimnames = list(names(subset_totalGSL), names(subset_totalGSL)))
pval_matrix_totalGSL <- matrix(NA, nrow = n_vars_totalGSL, ncol = n_vars_totalGSL, dimnames = list(names(subset_totalGSL), names(subset_totalGSL)))
r2_matrix_totalGSL <- matrix(NA, nrow = n_vars_totalGSL, ncol = n_vars_totalGSL, dimnames = list(names(subset_totalGSL), names(subset_totalGSL)))

# Calculate correlations, p-values, and R2 values
for (i in 1:n_vars_totalGSL) {
  for (j in i:n_vars_totalGSL) {
    if (i != j) {
      test <- cor.test(subset_totalGSL[[i]], subset_totalGSL[[j]])  # Perform correlation test
      cor_matrix_totalGSL[i, j] <- test$estimate  # Correlation coefficient
      cor_matrix_totalGSL[j, i] <- test$estimate  # Symmetric entry
      pval_matrix_totalGSL[i, j] <- test$p.value  # P-value
      pval_matrix_totalGSL[j, i] <- test$p.value  # Symmetric entry
      r2_matrix_totalGSL[i, j] <- test$estimate^2  # R-squared
      r2_matrix_totalGSL[j, i] <- test$estimate^2  # Symmetric entry
    } else {
      cor_matrix_totalGSL[i, j] <- 1  # Diagonal entries
      r2_matrix_totalGSL[i, j] <- 1
      pval_matrix_totalGSL[i, j] <- NA  # No p-value for diagonal
    }
  }
}

results_totalGSL_mats <- list(
  Correlation = cor_matrix_totalGSL,
  PValue = pval_matrix_totalGSL,
  R2 = r2_matrix_totalGSL
)

results_totalGSL_mats$Correlation

# save tables

write.csv(cor_matrix_totalGSL, "data/corr_mat_totalGSL_r.csv")
write.csv(r2_matrix_totalGSL, "data/corr_mat_totalGSL_r2.csv")
write.csv(pval_matrix_totalGSL, "data/corr_mat_totalGSL_p.csv")

#### Indoles ----

subset_indoles <- qst_fst_climate_indoles[, c("Qst", "Fst", "Qst_Fst", "Qst_Fst_abs", "cwd_Diff", "ppt_mm_Diff", "pck_Diff", "snw_Diff", "tmin_Diff", "tmax_Diff", "latitude_Diff", "longitude_Diff", "PC1_Diff", "PC2_Diff")]

# Initialize empty matrices to store results
n_vars_indoles <- ncol(subset_indoles)
cor_matrix_indoles <- matrix(NA, nrow = n_vars_indoles, ncol = n_vars_indoles, dimnames = list(names(subset_indoles), names(subset_indoles)))
pval_matrix_indoles <- matrix(NA, nrow = n_vars_indoles, ncol = n_vars_indoles, dimnames = list(names(subset_indoles), names(subset_indoles)))
r2_matrix_indoles <- matrix(NA, nrow = n_vars_indoles, ncol = n_vars_indoles, dimnames = list(names(subset_indoles), names(subset_indoles)))

# Calculate correlations, p-values, and R2 values
for (i in 1:n_vars_indoles) {
  for (j in i:n_vars_indoles) {
    if (i != j) {
      test <- cor.test(subset_indoles[[i]], subset_indoles[[j]])  # Perform correlation test
      cor_matrix_indoles[i, j] <- test$estimate  # Correlation coefficient
      cor_matrix_indoles[j, i] <- test$estimate  # Symmetric entry
      pval_matrix_indoles[i, j] <- test$p.value  # P-value
      pval_matrix_indoles[j, i] <- test$p.value  # Symmetric entry
      r2_matrix_indoles[i, j] <- test$estimate^2  # R-squared
      r2_matrix_indoles[j, i] <- test$estimate^2  # Symmetric entry
    } else {
      cor_matrix_indoles[i, j] <- 1  # Diagonal entries
      r2_matrix_indoles[i, j] <- 1
      pval_matrix_indoles[i, j] <- NA  # No p-value for diagonal
    }
  }
}

results_indoles_mats <- list(
  Correlation = cor_matrix_indoles,
  PValue = pval_matrix_indoles,
  R2 = r2_matrix_indoles
)

results_indoles_mats$Correlation

# save tables

write.csv(cor_matrix_indoles, "data/corr_mat_indoles_r.csv")
write.csv(r2_matrix_indoles, "data/corr_mat_indoles_r2.csv")
write.csv(pval_matrix_indoles, "data/corr_mat_indoles_p.csv")

#### Aliphatics ----

subset_aliphatics <- qst_fst_climate_aliphatics[, c("Qst", "Fst", "Qst_Fst", "Qst_Fst_abs", "cwd_Diff", "ppt_mm_Diff", "pck_Diff", "snw_Diff", "tmin_Diff", "tmax_Diff", "latitude_Diff", "longitude_Diff", "PC1_Diff", "PC2_Diff")]

# Initialize empty matrices to store results
n_vars_aliphatics <- ncol(subset_aliphatics)
cor_matrix_aliphatics <- matrix(NA, nrow = n_vars_aliphatics, ncol = n_vars_aliphatics, dimnames = list(names(subset_aliphatics), names(subset_aliphatics)))
pval_matrix_aliphatics <- matrix(NA, nrow = n_vars_aliphatics, ncol = n_vars_aliphatics, dimnames = list(names(subset_aliphatics), names(subset_aliphatics)))
r2_matrix_aliphatics <- matrix(NA, nrow = n_vars_aliphatics, ncol = n_vars_aliphatics, dimnames = list(names(subset_aliphatics), names(subset_aliphatics)))

# Calculate correlations, p-values, and R2 values
for (i in 1:n_vars_aliphatics) {
  for (j in i:n_vars_aliphatics) {
    if (i != j) {
      test <- cor.test(subset_aliphatics[[i]], subset_aliphatics[[j]])  # Perform correlation test
      cor_matrix_aliphatics[i, j] <- test$estimate  # Correlation coefficient
      cor_matrix_aliphatics[j, i] <- test$estimate  # Symmetric entry
      pval_matrix_aliphatics[i, j] <- test$p.value  # P-value
      pval_matrix_aliphatics[j, i] <- test$p.value  # Symmetric entry
      r2_matrix_aliphatics[i, j] <- test$estimate^2  # R-squared
      r2_matrix_aliphatics[j, i] <- test$estimate^2  # Symmetric entry
    } else {
      cor_matrix_aliphatics[i, j] <- 1  # Diagonal entries
      r2_matrix_aliphatics[i, j] <- 1
      pval_matrix_aliphatics[i, j] <- NA  # No p-value for diagonal
    }
  }
}

results_aliphatics_mats <- list(
  Correlation = cor_matrix_aliphatics,
  PValue = pval_matrix_aliphatics,
  R2 = r2_matrix_aliphatics
)

results_aliphatics_mats$Correlation

# save tables

write.csv(cor_matrix_aliphatics, "data/corr_mat_aliphatics_r.csv")
write.csv(r2_matrix_aliphatics, "data/corr_mat_aliphatics_r2.csv")
write.csv(pval_matrix_aliphatics, "data/corr_mat_aliphatics_p.csv")

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
  qst <- between_pop_var / (between_pop_var + 2*within_pop_var)
  
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

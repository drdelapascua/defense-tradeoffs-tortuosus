### NMDS Script ----
### Danielle De La Pascua
### 7-17-24

### libraries
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ade4)
library(funspace)
library(installr)
library(tibble)
library(mclust)
library(ggrepel)


### load data ----
dw <-  read.csv("./data/dw.csv") %>%
  select("Population", "mf", "rep", "Elevation", "treatment", "X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")
  #whatever order you put elect in is the order of the df, pop mf trt first

mf_means_nmds <-  read.csv("./data/mf_means.csv") %>%
  select("Population", "mf", "treatment", "X3MSO", "OHAlkenyl", "X4MSO", "Allyl", "X5MSO", "Butenyl", "X3MT", "MSOO", "OHI3M", "X4MT", "Flavonol16", "I3M", "Flavonol17", "Flavonol18", "Indole")
#whatever order you put elect in is the order of the df, pop mf trt first


### center and standardize data ----
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0
head(mf_means_nmds)
mf_means_nmds_scaling <- as.data.frame(scale(mf_means_nmds[4:18]))
mf_means_nmds_scaled <- mf_means_nmds[1:3]
mf_means_nmds_scaled[4:18] <- mf_means_nmds_scaling # merge two dataframes

# create data frame that is just columns with scaled compounds and labels as colnames
data_scaled2 = mf_means_nmds_scaled %>%
  # make and ID variable
  mutate(ID = paste(Population, treatment, mf, sep = "_")) %>%
  # make rownames combo of pop, tmt, mf, and rep
  column_to_rownames(var = "ID")  %>%
  select(-Population, -treatment, -mf)


# Perform NMDS (all data) ----

# bray - needs to be positive, could use gower
set.seed(123)
nmds_result2 <- metaMDS(data_scaled2, k = 2, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3 <- metaMDS(data_scaled2, k = 3, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4 <- metaMDS(data_scaled2, k = 4, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5 <- metaMDS(data_scaled2, k = 5, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result6 <- metaMDS(data_scaled2, k = 6, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result7 <- metaMDS(data_scaled2, k = 7, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)

#look at results, stress level should be under 0.1
nmds_result2 # 0.19
nmds_result3 # 0.14
nmds_result4 # 0.10 
nmds_result5 # 0.09
nmds_result6 # 0.10
nmds_result7 # 0.09 - stress dips below 0.1

nmds_result
stressplot(nmds_result2)
stressplot(nmds_result3)
stressplot(nmds_result4) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5)
stressplot(nmds_result7)

# Plot NMDS results
plot(nmds_result5, type = "points", display = "sites") # 'sites' for samples/observations


# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores <- scores(nmds_result5, display = "sites")  # Extract NMDS scores
scores # gives by each individual ID (rownames)

plot(nmds_result5)
envfit

# Calculate correlations between original variables and NMDS scores
correlations1 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 1]))
correlations2 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 2]))
correlations3 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 3]))
correlations4 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 4]))
correlations5 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 5]))
#correlations_6 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 6]))
#correlations_7 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 7]))

#when error - ask what is it doing/trying to do?

correlations1 # strongest is negative butenyl, -.8. Allyl is .5, 4MSO -.6, 3MSO is .65
correlations_2 # pos associated with Flavonol-18 & 3MT, neg associated w/ 5MSO
correlations_3 # pos associated with OH-Alkenyl, neg associated with Allyl, 5MSO & Flavonol 18
correlations_4 # pos associated with i3M and 4MSO, neg associated with indole
correlations_5 # pos associated with indole, neg associated with flavonol-16 and OH-Alkenyl

# Combine all correlations into a data frame
correlation_table <- data.frame(
  Variable = colnames(data_scaled2),
  NMDS1 = correlations1,
  NMDS2 = correlations2,
  NMDS3 = correlations3,
  NMDS4 = correlations4,
  NMDS5 = correlations5
)

# save dataframe
write.csv(correlation_table, file = "output/NMDS_axes_correlation_table.csv")

#plot correlations
barplot(correlations_1, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")


barplot(correlations_2, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_3, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_4, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_5, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_6, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_7, names.arg = colnames(data_scaled2),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

####visualize NMDS ----

elevation <- read.csv(file = "data/elevation.csv") %>%
  mutate(Population = ifelse(Population == "YOSE10", "YO10", Population))

mf_means_nmds_scaled_withel <- mf_means_nmds_scaled %>%
  left_join(elevation, by = "Population") %>%
  select(-c("Lat", "Long", "Seed.year"))

#create grouping variable
all_groups = mf_means_nmds_scaled_withel %>%
  select(Population, treatment, mf, Elevation)

nmds_dist_scores = as.data.frame(scores(nmds_result5, "sites")) %>%
  mutate(population = all_groups$Population, treatment = all_groups$treatment, mf = all_groups$mf, elevation = all_groups$Elevation)

# plot it 
nmds_dist_plot_trt = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt

nmds_dist_plot_el = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = elevation, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  facet_wrap(~treatment) +
  labs(color = "Elevation", shape= "treatment") 

nmds_dist_plot_el # need to add el back in 12-12-24

nmds_dist_plot_pop = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = elevation, shape = treatment), size = 5, alpha = 0.75) + 
  geom_text_repel(aes(label = all_groups$Population), size = 3, max.overlaps = 10) +  # Add repelled labels
  theme_bw() +
  facet_wrap(~treatment) +
  labs(color = "Elevation", shape = "Treatment") +
  theme(legend.position = "right")  +
  scale_color_gradient(low="orange", high="blue")


nmds_dist_plot_pop


#zeros are creating the straight lines - maybe 0-0 or something between groups 

# axes 2 & 3
nmds_dist_plot_pop_2_3 = ggplot(data = nmds_dist_scores, aes(x = NMDS2, y= NMDS3)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 
nmds_dist_plot_pop_2_3

# plot only trt
nmds_dist_plot_trt23 = ggplot(data = nmds_dist_scores, aes(x = NMDS2, y= NMDS3)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt23


#axes 3 & 4

#population & treatment
nmds_dist_plot_pop_3_4 = ggplot(data = nmds_dist_scores, aes(x = NMDS3, y= NMDS4)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

nmds_dist_plot_pop_3_4

#just treatment
nmds_dist_plot_trt34 = ggplot(data = nmds_dist_scores, aes(x = NMDS3, y= NMDS4)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt34

#axes 4 & 5

#population & treatment
nmds_dist_plot_pop_4_5 = ggplot(data = nmds_dist_scores, aes(x = NMDS4, y= NMDS5)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

nmds_dist_plot_pop_4_5

#just treatment
nmds_dist_plot_trt45 = ggplot(data = nmds_dist_scores, aes(x = NMDS4, y= NMDS5)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt45


#sleuthing out zeroes 

# Count number of values greater than zero 
count_5MSO_zero <- sum(data_scaled2$X3MSO_5.2 > 0) # 105

count_OH_Alkenyl_zero <- sum(data_scaled2$OH.Alkenyl_6 > 0) #30 - potentially exclude, likely introducing lots of zeroes?

count_4MSO_zero <- sum(data_scaled2$X4MSO_7.1 > 0) #100

count_Allyl_zero <- sum(data_scaled2$Allyl_7.4 > 0) #85

count_5MSO_zero <- sum(data_scaled2$X5MSO_10.2 > 0) # 90

count_Butenyl_zero <- sum(data_scaled2$Butenyl_12.1 > 0) #104

count_3MT_zero <- sum(data_scaled2$X3MT_13.6 > 0) # 98

count_MSOO_zero <- sum(data_scaled2$MSOO_13.8 > 0) #73

count_OH_I3M_zero <- sum(data_scaled2$OH.I3M_15.1 > 0) # 43 - potentially exclude, likely introducing lots of zeroes?

count_4MT_zero <- sum(data_scaled2$X4MT._15.5 > 0) # 107

count_Flavonol16_zero <- sum(data_scaled2$Flavonol_16.1 > 0) # 120

count_I3M_zero <- sum(data_scaled2$I3M_16.7 > 0) # 83

count_Flavonol17_zero <- sum(data_scaled2$Flavonol_17.5 > 0) # 81

count_Flavonol18_zero <- sum(data_scaled2$Flavonol_18.5 > 0) # 109

count_Indole_zero <- sum(data_scaled2$Indole_18.8 > 0) # 93

#### Permanova ----

summary(data_scaled2)

all_gower_rel = vegdist(data_scaled2, method = 'gower')
permanova_gower = adonis2(all_bray_rel ~ Population*treatment, perm= 999, data = all_groups)
summary(permanova_gower)
permanova_gower

# pops are different, treatments are different, but the way theyre responding across pops is the same

# cluster analysis - heirarchical clustering
nmds_scores <- nmds_result7$points
hc <- hclust(dist(nmds_scores))  # Compute distances between NMDS scores
plot(hc)  # Plot dendrogram
clusters <- cutree(hc, k = 4)  # Cut dendrogram into 3 clusters (adjust k as needed)
plot(clusters)
clusters

# plotting clusters on NMDS plot
plot(nmds_result7, display = "sites", type = "n")  # Plot NMDS configuration
points(nmds_result7, display = "sites", col = clusters)  # Add points colored by cluster membership

#using paper methods
gmm_model <- Mclust(nmds_scores)

# Extract cluster assignments
clusters <- gmm_model$classification
clusters
# Example scatter plot of NMDS scores with clusters colored
plot(nmds_scores, col = clusters, pch = 16, main = "Cluster Analysis on NMDS Scores")
legend("topright", legend = unique(clusters), col = unique(clusters), pch = 16, title = "Cluster")

# Example scatter plot of NMDS scores with GMM clusters colored - 9 clusters
plot(nmds_scores, col = clusters, pch = 16, main = "Gaussian Mixture Model Clustering on NMDS Scores")
legend("topright", legend = unique(clusters), col = unique(clusters), pch = 16, title = "Cluster")

# Perform NMDS (remove extreme values) ----

# load data 
mf_means_nmds_ex_rem <-  read.csv("./data/mf_means_ex_rem.csv") %>%
  select("Population", "mf", "treatment", "X3MSO", "OHAlkenyl", "X4MSO", "Allyl", "X5MSO", "Butenyl", "X3MT", "MSOO", "OHI3M", "X4MT", "Flavonol16", "I3M", "Flavonol17", "Flavonol18", "Indole")

#### center and standardize data ----
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0
head(mf_means_nmds_ex_rem)
mf_means_nmds_scaling_ex_rem <- as.data.frame(scale(mf_means_nmds_ex_rem[4:18]))
mf_means_nmds_scaled_ex_rem <- mf_means_nmds_ex_rem[1:3]
mf_means_nmds_scaled_ex_rem[4:18] <- mf_means_nmds_scaling_ex_rem # merge two dataframes

# create data frame that is just columns with scaled compounds and labels as colnames
data_scaled2_ex_rem = mf_means_nmds_scaled_ex_rem %>%
  # make and ID variable
  mutate(ID = paste(Population, treatment, mf, sep = "_")) %>%
  # make rownames combo of pop, tmt, mf, and rep
  column_to_rownames(var = "ID")  %>%
  select(-Population, -treatment, -mf)



# bray - needs to be positive, could use gower
set.seed(123)
nmds_result2_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 2, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 3, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 4, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 5, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result6_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 6, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result7_ex_rem <- metaMDS(data_scaled2_ex_rem, k = 7, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)

#look at results, stress level should be under 0.1
nmds_result2_ex_rem # 0.22
nmds_result3_ex_rem # 0.17
nmds_result4_ex_rem # 0.14 
nmds_result5_ex_rem # best sol
nmds_result6_ex_rem # 0.10
nmds_result7_ex_rem # 0.09 - stress dips below 0.1

nmds_result
stressplot(nmds_result2_ex_rem)
stressplot(nmds_result3_ex_rem)
stressplot(nmds_result4_ex_rem) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5_ex_rem)
stressplot(nmds_result7_ex_rem)

# Plot NMDS results
plot(nmds_result5_ex_rem, type = "points", display = "sites")  # 'sites' for samples/observations


# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores_ex_rem <- scores(nmds_result5_ex_rem, display = "sites")  # Extract NMDS scores
scores_ex_rem # gives by each individual ID (rownames)

# Calculate correlations between original variables and NMDS scores
correlations1_ex_rem <- apply(data_scaled2_ex_rem, 2, function(x) cor(x, scores_ex_rem[, 1]))
correlations2_ex_rem <- apply(data_scaled2_ex_rem, 2, function(x) cor(x, scores_ex_rem[, 2]))
correlations3_ex_rem <- apply(data_scaled2_ex_rem, 2, function(x) cor(x, scores_ex_rem[, 3]))
correlations4_ex_rem <- apply(data_scaled2_ex_rem, 2, function(x) cor(x, scores_ex_rem[, 4]))
correlations5_ex_rem <- apply(data_scaled2_ex_rem, 2, function(x) cor(x, scores_ex_rem[, 5]))
#correlations_6 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 6]))
#correlations_7 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 7]))

# Combine all correlations into a data frame
correlation_table_ex_rem <- data.frame(
  Variable = colnames(data_scaled2_ex_rem),
  NMDS1 = correlations1_ex_rem,
  NMDS2 = correlations2_ex_rem,
  NMDS3 = correlations3_ex_rem,
  NMDS4 = correlations4_ex_rem,
  NMDS5 = correlations5_ex_rem
)

# save dataframe
write.csv(correlation_table_ex_rem, file = "output/NMDS_axes_correlation_table_ex_rem.csv")

### Perform NMDS (without OH-I3M)

# Perform NMDS (remove extreme values) ----

# load data 
mf_means_nmds_rem_ohi3m <-  read.csv("./data/mf_means_rem_ohi3m.csv") %>%
  select("Population", "mf", "treatment", "X3MSO", "OHAlkenyl", "X4MSO", "Allyl", "X5MSO", "Butenyl", "X3MT", "MSOO", "X4MT", "Flavonol16", "I3M", "Flavonol17", "Flavonol18", "Indole")

#### center and standardize data ----
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0
head(mf_means_nmds_rem_ohi3m)
dim(mf_means_nmds_rem_ohi3m)
mf_means_nmds_scaling_rem_ohi3m <- as.data.frame(scale(mf_means_nmds_rem_ohi3m[4:17]))
mf_means_nmds_scaled_rem_ohi3m <- mf_means_nmds_rem_ohi3m[1:3]
mf_means_nmds_scaled_rem_ohi3m[4:17] <- mf_means_nmds_scaling_rem_ohi3m # merge two dataframes

# create data frame that is just columns with scaled compounds and labels as colnames
data_scaled2_rem_ohi3m = mf_means_nmds_scaled_rem_ohi3m %>%
  # make and ID variable
  mutate(ID = paste(Population, treatment, mf, sep = "_")) %>%
  # make rownames combo of pop, tmt, mf, and rep
  column_to_rownames(var = "ID")  %>%
  select(-Population, -treatment, -mf)


# bray - needs to be positive, could use gower
set.seed(123)
nmds_result2_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 2, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 3, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 4, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 5, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result6_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 6, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result7_rem_ohi3m <- metaMDS(data_scaled2_rem_ohi3m, k = 7, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)

#look at results, stress level should be under 0.1
nmds_result2_rem_ohi3m # 0.22
nmds_result3_rem_ohi3m # 0.17
nmds_result4_rem_ohi3m # 0.14 
nmds_result5_rem_ohi3m # best sol
nmds_result6_rem_ohi3m # 0.10
nmds_result7_rem_ohi3m # 0.09 - stress dips below 0.1

nmds_result
stressplot(nmds_result2_rem_ohi3m)
stressplot(nmds_result3_rem_ohi3m)
stressplot(nmds_result4_rem_ohi3m) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5_rem_ohi3m)
stressplot(nmds_result7_rem_ohi3m)

# Plot NMDS results
plot(nmds_result5_rem_ohi3m, type = "points", display = "sites")  # 'sites' for samples/observations


# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores_rem_ohi3m <- scores(nmds_result5_rem_ohi3m, display = "sites")  # Extract NMDS scores
scores_rem_ohi3m # gives by each individual ID (rownames)

# Calculate correlations between original variables and NMDS scores
correlations1_rem_ohi3m <- apply(data_scaled2_rem_ohi3m, 2, function(x) cor(x, scores_rem_ohi3m[, 1]))
correlations2_rem_ohi3m <- apply(data_scaled2_rem_ohi3m, 2, function(x) cor(x, scores_rem_ohi3m[, 2]))
correlations3_rem_ohi3m <- apply(data_scaled2_rem_ohi3m, 2, function(x) cor(x, scores_rem_ohi3m[, 3]))
correlations4_rem_ohi3m <- apply(data_scaled2_rem_ohi3m, 2, function(x) cor(x, scores_rem_ohi3m[, 4]))
correlations5_rem_ohi3m <- apply(data_scaled2_rem_ohi3m, 2, function(x) cor(x, scores_rem_ohi3m[, 5]))
#correlations_6 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 6]))
#correlations_7 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 7]))

# Combine all correlations into a data frame
correlation_table_rem_ohi3m <- data.frame(
  Variable = colnames(data_scaled2_rem_ohi3m),
  NMDS1 = correlations1_rem_ohi3m,
  NMDS2 = correlations2_rem_ohi3m,
  NMDS3 = correlations3_rem_ohi3m,
  NMDS4 = correlations4_rem_ohi3m,
  NMDS5 = correlations5_rem_ohi3m
)

# save dataframe
write.csv(correlation_table_rem_ohi3m, file = "output/NMDS_axes_correlation_table_rem_ohi3m.csv")


### visualize NMDS ----

elevation <- read.csv(file = "data/elevation.csv") %>%
  mutate(Population = ifelse(Population == "YOSE10", "YO10", Population))

mf_means_nmds_scaled_withel <- mf_means_nmds_scaled %>%
  left_join(elevation, by = "Population") %>%
  select(-c("Lat", "Long", "Seed.year"))

#create grouping variable
all_groups = mf_means_nmds_scaled_withel %>%
  select(Population, treatment, mf, Elevation)

nmds_dist_scores = as.data.frame(scores(nmds_result5, "sites")) %>%
  mutate(population = all_groups$Population, treatment = all_groups$treatment, mf = all_groups$mf, elevation = all_groups$Elevation)

# plot it 
nmds_dist_plot_trt = ggplot(data = nmds_dist_scores_, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt

nmds_dist_plot_el = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = elevation, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  facet_wrap(~treatment) +
  labs(color = "Elevation", shape= "treatment") 

nmds_dist_plot_el # need to add el back in 12-12-24

nmds_dist_plot_pop = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = elevation, shape = treatment), size = 5, alpha = 0.75) + 
  geom_text_repel(aes(label = all_groups$Population), size = 3, max.overlaps = 10) +  # Add repelled labels
  theme_bw() +
  facet_wrap(~treatment) +
  labs(color = "Elevation", shape = "Treatment") +
  theme(legend.position = "right")  +
  scale_color_gradient(low="orange", high="blue")


nmds_dist_plot_pop


#zeros are creating the straight lines - maybe 0-0 or something between groups 

# axes 2 & 3
nmds_dist_plot_pop_2_3 = ggplot(data = nmds_dist_scores, aes(x = NMDS2, y= NMDS3)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 
nmds_dist_plot_pop_2_3

# plot only trt
nmds_dist_plot_trt23 = ggplot(data = nmds_dist_scores, aes(x = NMDS2, y= NMDS3)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt23


#axes 3 & 4

#population & treatment
nmds_dist_plot_pop_3_4 = ggplot(data = nmds_dist_scores, aes(x = NMDS3, y= NMDS4)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

nmds_dist_plot_pop_3_4

#just treatment
nmds_dist_plot_trt34 = ggplot(data = nmds_dist_scores, aes(x = NMDS3, y= NMDS4)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt34

#axes 4 & 5

#population & treatment
nmds_dist_plot_pop_4_5 = ggplot(data = nmds_dist_scores, aes(x = NMDS4, y= NMDS5)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

nmds_dist_plot_pop_4_5

#just treatment
nmds_dist_plot_trt45 = ggplot(data = nmds_dist_scores, aes(x = NMDS4, y= NMDS5)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt45


#sleuthing out zeroes 

# Count number of values greater than zero 
count_5MSO_zero <- sum(data_scaled2$X3MSO_5.2 > 0) # 105

count_OH_Alkenyl_zero <- sum(data_scaled2$OH.Alkenyl_6 > 0) #30 - potentially exclude, likely introducing lots of zeroes?

count_4MSO_zero <- sum(data_scaled2$X4MSO_7.1 > 0) #100

count_Allyl_zero <- sum(data_scaled2$Allyl_7.4 > 0) #85

count_5MSO_zero <- sum(data_scaled2$X5MSO_10.2 > 0) # 90

count_Butenyl_zero <- sum(data_scaled2$Butenyl_12.1 > 0) #104

count_3MT_zero <- sum(data_scaled2$X3MT_13.6 > 0) # 98

count_MSOO_zero <- sum(data_scaled2$MSOO_13.8 > 0) #73

count_OH_I3M_zero <- sum(data_scaled2$OH.I3M_15.1 > 0) # 43 - potentially exclude, likely introducing lots of zeroes?

count_4MT_zero <- sum(data_scaled2$X4MT._15.5 > 0) # 107

count_Flavonol16_zero <- sum(data_scaled2$Flavonol_16.1 > 0) # 120

count_I3M_zero <- sum(data_scaled2$I3M_16.7 > 0) # 83

count_Flavonol17_zero <- sum(data_scaled2$Flavonol_17.5 > 0) # 81

count_Flavonol18_zero <- sum(data_scaled2$Flavonol_18.5 > 0) # 109

count_Indole_zero <- sum(data_scaled2$Indole_18.8 > 0) # 93

### Permanova ----

summary(data_scaled2)

all_gower_rel = vegdist(data_scaled2, method = 'gower')
permanova_gower = adonis2(all_gower_rel ~ Population + treatment + Population:treatment, perm= 999, data = all_groups, by = "terms")
summary(permanova_gower)
permanova_gower

### Old code ----

# filter for only controls
#control_pops <- filter(pop_means, treatment == "C")
#induced_pops <- filter(pop_means, treatment == "CW")

### center and standardize data
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0

# pop means not c
data_scaled <- as.data.frame(scale(control_pops[4:18]))
data_scaled_labels <- control_pops[2] # vector with id info for above
rownames(data_scaled) = data_scaled_labels

data_scaled = control_pops %>%
  # make rownames combo of pop and treatment
  
  # Perform NMDS
  set.seed(123)
nmds_result <- metaMDS(data_scaled, k = 2)  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result2 <- metaMDS(data_scaled, k = 2, weakties = FALSE, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3 <- metaMDS(data_scaled, k = 3, weakties = FALSE, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4 <- metaMDS(data_scaled, k = 4, weakties = FALSE, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5 <- metaMDS(data_scaled, k = 5, weakties = FALSE, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
?metaMDS
nmds_result
nmds_result2
nmds_result3
nmds_result4

nmds_result
stressplot(nmds_result2)
stressplot(nmds_result3)
stressplot(nmds_result4) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5)

# Plot NMDS results
plot(nmds_result4, type = "points", display = "sites")  # 'sites' for samples/observations

# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores <- scores(nmds_result4, display = "sites")  # Extract NMDS scores
scores
nmds_result4
# Calculate correlations between original variables and NMDS scores
correlations_1 <- apply(data_scaled, 2, function(x) cor(x, scores[, 1]))
correlations_2 <- apply(data_scaled, 2, function(x) cor(x, scores[, 2]))
correlations_3 <- apply(data_scaled, 2, function(x) cor(x, scores[, 3]))
correlations_4 <- apply(data_scaled, 2, function(x) cor(x, scores[, 4]))
#when error - ask what is it doing/trying to do?

nmds_result

correlations_1 # Flavonol 16, Allyl, 5MSO are loading negatively, Flavonol 18 and indole are loading positively
correlations_2 # Butenyl positive corr, Flavonol 16 slightly negatively corr, Allyl, Flavonol 18 & Indole neg corr
correlations_3 # Butenyl positive corr, Flavonol 16 slightly negatively corr, Allyl, Flavonol 18 & Indole neg corr
correlations_4 # Butenyl positive corr, Flavonol 16 slightly negatively corr, Allyl, Flavonol 18 & Indole neg corr

# for both loadings, only Allyl, Butenyl, Flavonol 16, Flavonol 18, and Indole are showing up
# One observation, these are among the most abundant - maybe the others dont have enough data?

#plot correlations
barplot(correlations_1, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_2, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

# plot with labels
data <- data_scaled
data[16] <- data_scaled_labels
points <- nmds_result$points
points <- as.data.frame(points) 
points[3] <- data_scaled_labels

ggplot(points, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point() +  # Plot points
  geom_text(nudge_x = 0.05) +  # Add text labels with a slight offset
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", title = "NMDS Plot with Labels")
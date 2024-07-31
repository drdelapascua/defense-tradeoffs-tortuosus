### NMDS Script
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

### load data
dl <-  read.csv("./data/dw.csv") %>%
  select("Population", "mf", "rep", "treatment", "X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4", "X5MSO_10.2", "Butenyl_12.1", "X3MT_13.6", "MSOO_13.8", "OH.I3M_15.1", "X4MT._15.5", "Flavonol_16.1", "I3M_16.7", "Flavonol_17.5", "Flavonol_18.5", "Indole_18.8")
  #whatever order you put elect in is the order of the df, pop mf trt first


### center and standardize data
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0
dl_scaling <- as.data.frame(scale(dl[5:19]))
dl_scaled <- dl[1:4]
dl_scaled[5:19] <- dl_scaling # merge two dataframes

# create data frame that is just columns with scaled compounds and labels as colnames
data_scaled2 = dl_scaled %>%
  # make and ID variable
  mutate(ID = paste(Population, treatment, mf, rep, sep = "_")) %>%
  # make rownames combo of pop, tmt, mf, and rep
  column_to_rownames(var = "ID")  %>%
  select(-Population, -treatment, -mf, -rep)
  
# Perform NMDS

# bray - needs to be positive, could use gower
set.seed(123)
nmds_result2 <- metaMDS(data_scaled2, k = 2, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3 <- metaMDS(data_scaled2, k = 3, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4 <- metaMDS(data_scaled2, k = 4, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5 <- metaMDS(data_scaled2, k = 5, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result6 <- metaMDS(data_scaled2, k = 6, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result7 <- metaMDS(data_scaled2, k = 7, distance = "gower")  # Adjust 'k' as needed (usually 2 or 3 dimensions)

#look at results, stress level should be under 0.1
nmds_result2 # 0.22
nmds_result3 # 0.17
nmds_result4 # 0.14 
nmds_result5 # best sol
nmds_result6 # 0.10
nmds_result7 # 0.09 - stress dips below 0.1

nmds_result
stressplot(nmds_result2)
stressplot(nmds_result3)
stressplot(nmds_result4) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5)
stressplot(nmds_result7)

# Plot NMDS results
plot(nmds_result5, type = "points", display = "sites")  # 'sites' for samples/observations


# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores <- scores(nmds_result5, display = "sites")  # Extract NMDS scores
scores # gives by each individual ID (rownames)

# Calculate correlations between original variables and NMDS scores
correlations1 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 1]))
correlations2 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 2]))
correlations3 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 3]))
correlations4 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 4]))
correlations5 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 5]))
#correlations_6 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 6]))
#correlations_7 <- apply(data_scaled2, 2, function(x) cor(x, scores[, 7]))

#when error - ask what is it doing/trying to do?

correlations_1 # strongest is negative butenyl, -.8. Allyl is .5, 4MSO -.6, 3MSO is .65
correlations_2 # 
correlations_3
correlations_4 
correlations_5 
correlations_6
correlations_7

df <- data.frame(
  Axis_1 = correlations_1,
  Axis_2 = correlations_2,
  Axis_3 = correlations_3,
  Axis_4 = correlations_4,
  Axis_5 = correlations_5
)

df

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

#visualize NMDS

#create grouping variable
all_groups = dl_scaled %>%
  select(Population, treatment, mf, rep)

nmds_dist_scores = as.data.frame(scores(nmds_result5, "sites")) %>%
  mutate(population = all_groups$Population, treatment = all_groups$treatment, mf = all_groups$mf, rep = all_groups$rep)

# plot it 
nmds_dist_plot_trt = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  scale_color_manual(values = c("firebrick4" , "darkolivegreen4", "dodgerblue4")) +
  scale_shape_manual(values = c(1,18))+
  labs(color = "treatment", shape= "population") 

nmds_dist_plot_trt

nmds_dist_plot_pop = ggplot(data = nmds_dist_scores, aes(x = NMDS1, y= NMDS2)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

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

# Permanova

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

############################ OLD CODE ##########################################


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
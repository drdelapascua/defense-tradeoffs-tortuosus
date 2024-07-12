### NMDS Script
### Danielle De La Pascua
### 6-30-24

### libraries
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ade4)
library(funspace)

### load data
dl <-  read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv") %>%
  select(-starts_with("Junk")) %>% #remove any columns that start with Junk 

dl2 <- dl %>%
  select(c("X3MSO_5.2", "OH.Alkenyl_6", "X4MSO_7.1", "Allyl_7.4")) #whatever order you put elect in is the order of the df, pop mf trt first

dl2 = dl[,c(11:28)]

#make NAs = 0 (NA means compound wasn't present in the sample)
dl3 = dl2 %>%
  mutate(across(c(1:15), ~replace_na(.x, 0))) 
  
dim(dl2)
names(dl)

# Danielle clean up above and below :)

#dl <- dl[-17] #getting rid of junk
#dl <- dl[-25] #getting rid of junk

# filter for only controls
#control_pops <- filter(pop_means, treatment == "C")
#induced_pops <- filter(pop_means, treatment == "CW")

### center and standardize data
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0
dl_scaled <- as.data.frame(scale(dl3[1:15]))
dl_labels <- dl3[16:18]
dl_scaled[16:18] <- dl_labels # merge two dataframes

### exported df above and added reps for below - add to code 

dl_scaled_2 <-  read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/manual_reps.csv")


dl_scaled_2 <- dl_scaled_2 %>%
  mutate(ID = paste(Population, treatment, mf, rep, sep = "_"))

row.names(dl_scaled) = dl_scaled_2$ID

# create data frame that is just columns with scaled compounds and labels as colnames
data_scaled = dl_scaled %>%
  # make rownames combo of pop, tmt, mf, and rep
  column_to_rownames(var = "ID")  %>%
  select(-Population, -treatment, -mf, -rep) %>%
  
# Perform NMDS
set.seed(123)
nmds_result <- metaMDS(data_scaled, k = 2)  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result2 <- metaMDS(data_scaled, k = 2, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result3 <- metaMDS(data_scaled, k = 3, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result4 <- metaMDS(data_scaled, k = 4, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result5 <- metaMDS(data_scaled, k = 5, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result6 <- metaMDS(data_scaled, k = 6, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)
nmds_result7 <- metaMDS(data_scaled, k = 7, distance = "bray")  # Adjust 'k' as needed (usually 2 or 3 dimensions)

#look at results, stress level should be under 0.1
nmds_result
nmds_result2
nmds_result3
nmds_result4
nmds_result5
nmds_result6
nmds_result7 # stress dips below 0.1

nmds_result
stressplot(nmds_result2)
stressplot(nmds_result3)
stressplot(nmds_result4) # this seems the best, going up to k = 5 not change much
stressplot(nmds_result5)
stressplot(nmds_result7)

# Plot NMDS results
plot(nmds_result7, type = "points", display = "sites")  # 'sites' for samples/observations

# Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores <- scores(nmds_result7, display = "sites")  # Extract NMDS scores
scores
nmds_result4
# Calculate correlations between original variables and NMDS scores
correlations_1 <- apply(data_scaled, 2, function(x) cor(x, scores[, 1]))
correlations_2 <- apply(data_scaled, 2, function(x) cor(x, scores[, 2]))
correlations_3 <- apply(data_scaled, 2, function(x) cor(x, scores[, 3]))
correlations_4 <- apply(data_scaled, 2, function(x) cor(x, scores[, 4]))
correlations_5 <- apply(data_scaled, 2, function(x) cor(x, scores[, 5]))
correlations_6 <- apply(data_scaled, 2, function(x) cor(x, scores[, 6]))
correlations_7 <- apply(data_scaled, 2, function(x) cor(x, scores[, 7]))

#when error - ask what is it doing/trying to do?

correlations_1 # potentially describing a t-o between indoles and alkenyl compounds?
correlations_2 # 
correlations_3
correlations_4 
correlations_5 
correlations_6
correlations_7

# for both loadings, only Allyl, Butenyl, Flavonol 16, Flavonol 18, and Indole are showing up
# One observation, these are among the most abundant - maybe the others dont have enough data?

#plot correlations
barplot(correlations_1, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")


barplot(correlations_2, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_3, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_4, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_5, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_6, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

barplot(correlations_7, names.arg = colnames(data_scaled),
        main = "Correlations with NMDS Axis 1",
        xlab = "Variables", ylab = "Correlation")

#visualize NMDS

#create grouping variable
all_groups = dl_scaled %>%
  select(Population, treatment, mf, rep)

nmds_dist_scores = as.data.frame(scores(nmds_result7, "sites")) %>%
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

#zeros are creating the straight lines - maybe 0-0 or something between groups 

#sleuthing out zeroes 

# Count number of values greater than zero 
count_5MSO_zero <- sum(data_scaled$X3MSO_5.2 > 0)
count_5MSO_zero

count_OH_Alkenyl_zero <- sum(data_scaled$OH.Alkenyl_6 > 0)
count_OH_Alkenyl_zero

count_4MSO_zero <- sum(data_scaled$X4MSO_7.1 > 0)
count_4MSO_zero

count_Allyl_zero <- sum(data_scaled$Allyl_7.4 > 0)
count_Allyl_zero


count_5MSO_zero <- sum(data_scaled$X5MSO_10.2 > 0)
count_5MSO_zero

count_Butenyl_zero <- sum(data_scaled$Butenyl_12.1 > 0)
count_Butenyl_zero

count_3MT_zero <- sum(data_scaled$X3MT_13.6 > 0)
count_3MT_zero

count_MSOO_zero <- sum(data_scaled$MSOO_13.8 > 0)
count_MSOO_zero

count_OH_I3M_zero <- sum(data_scaled$OH.I3M_15.1 > 0)
count_OH_I3M_zero

count_4MT_zero <- sum(data_scaled$X4MT._15.5 > 0)
count_4MT_zero

count_Flavonol16_zero <- sum(data_scaled$Flavonol_16.1 > 0)
count_Flavonol16_zero 

count_I3M_zero <- sum(data_scaled$I3M_16.7 > 0)
count_I3M_zero

count_Flavonol17_zero <- sum(data_scaled$Flavonol_17.5 > 0)
count_Flavonol17_zero


count_Flavonol18_zero <- sum(data_scaled$Flavonol_18.5 > 0)
count_Flavonol18_zero

nmds_dist_plot_pop

correlations_1 
correlations_2

# axes 2 & 3
nmds_dist_plot_pop_2_3 = ggplot(data = nmds_dist_scores, aes(x = NMDS2, y= NMDS3)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 
nmds_dist_plot_pop_2_3
#so strange!!!


#axes 3 & 4
nmds_dist_plot_pop_3_4 = ggplot(data = nmds_dist_scores, aes(x = NMDS3, y= NMDS4)) +
  geom_point(aes(colour = population, shape = treatment), size= 5, alpha = 0.75)   + theme_bw() +
  labs(color = "population", shape= "treatment") 

nmds_dist_plot_pop_3_4

# fun space
nmds_coords <- nmds_result7$points
nmds_data <- data.frame(X = nmds_coords[,1], Y = nmds_coords[,2], Group = all_groups$treatment)
# Create NMDS plot with groups
funplot(nmds_data[,1:2], group = nmds_data$Group, cex = 1.5) # getting an error here - why?

# Permanova

all_bray_rel = vegdist(data_scaled, method = 'bray')
permanova_bray = adonis2(all_bray_rel ~ Population*treatment, perm= 999, data = all_groups)
summary(permanova_bray)
permanova_bray

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
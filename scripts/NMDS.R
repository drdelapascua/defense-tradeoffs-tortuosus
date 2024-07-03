### NMDS Script
### Danielle De La Pascua
### 6-30-24

### libraries
library(vegan)
library(dplyr)
library(ggplot2)

### load data
pop_means <-  read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/pop_means.csv")

# filter for only controls
control_pops <- filter(pop_means, treatment == "C")
induced_pops <- filter(pop_means, treatment == "CW")

### center and standardize data
# standardizing the data (all will have variance = 1)
# centering gives all a mean of 0

data_scaled <- as.data.frame(scale(control_pops[4:18]))
data_scaled_labels <- control_pops[2] # vector with id info for above

# Calculate dissimilarities (for example, using Euclidean distance)
distances <- dist(data_scaled, method = "euclidean")

# Perform NMDS
set.seed(123)
nmds_result <- metaMDS(distances, k = 2)  # Adjust 'k' as needed (usually 2 or 3 dimensions)

# Plot NMDS results
plot(nmds_result, type = "points", display = "sites")  # 'sites' for samples/observations

# # Assuming 'nmds_result' is your NMDS result object from 'metaMDS'
scores <- scores(nmds_result, display = "sites")  # Extract NMDS scores

# Calculate correlations between original variables and NMDS scores
correlations_1 <- apply(data_scaled, 2, function(x) cor(x, scores[, 1]))
correlations_2 <- apply(data_scaled, 2, function(x) cor(x, scores[, 2]))

nmds_result

correlations_1 # Flavonol 16, Allyl, 5MSO are loading negatively, Flavonol 18 and indole are loading positively
correlations_2 # Butenyl positive corr, Flavonol 16 slightly negatively corr, Allyl, Flavonol 18 & Indole neg corr
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

# Load necessary libraries
library(dplyr)
library(tidyverse)

# populations need to be numbers not names/code? not sure

# some pops only in elenas study, some only in mine, we just want overlap

# make genetic dummy data

# Define row names and column names
#row_names <- c("ST13L16", "ST1L25", "ST5L2D5", "ST1L43", "ST4L24", "ST15L33", "ST2L28")
#col_names <- c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4")
# Generate random values
#set.seed(123)  # Set seed for reproducibility (optional)
#values <- sample(c("11", "12", "21", "22"), length(row_names) * length(col_names), replace = TRUE)
# Create the data frame
#dummy_df <- matrix(values, nrow = length(row_names), ncol = length(col_names),
#                   dimnames = list(row_names, col_names))
# Convert matrix to data frame
#dummy_df <- as.data.frame(dummy_df)

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

# Save the merged data as a new CSV file
write.csv(filtered_merged_data, "data/merged_data.csv", row.names = FALSE)


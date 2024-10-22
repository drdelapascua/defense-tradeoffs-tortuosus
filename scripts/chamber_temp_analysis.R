# temp 

# load packages
library(dplyr)
library(tidyverse)

# load data 
gc422 <- read.csv("./data/GC_422-Rformatted.csv") %>%
  filter(Date != c("9/19/24", "9/20/24", "9/21/24", "9/22/24", "9/23/24"))

sliced_gc422 <- gc422 %>% 
  slice(-1:-93)

round1 <- read.csv("./data/2021_round3_Rformatted.csv")



colnames(sliced_gc422) <- c("Date", "Time", "temp_2024")
colnames(round1) <- c("Date", "Time", "Unit", "temp_2021")

# merge data 
data <- full_join(x = sliced_gc422, y = round1)
dim(gc422)
dim(round1)
dim(data)
head(data)

# pivot longer
dl <- data %>%
  pivot_longer(cols = c(temp_2021, temp_2024), names_to = "temp_year", values_to = "temp")

dl$Time <- hms(dl$Time)
summary(dl)

ggplot(data = dl, aes(x = Time, y = temp, group = temp_year, color = temp_year)) + 
  geom_line() +
  geom_point() +  # Optional: add points for better visibility
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Temperature Comparison: 2021 vs 2024",
       x = "Time",
       y = "Temperature (°C)",
       color = "Year") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(maps)
library(rworldmap)
install.packages("dplyr")
install.packages("")
setwd("~/GitHub/defense-tradeoffs-tortuosus")
updateR()

#### Map of populations for defense trade-off project

#load points
pops <- read.csv("defense_tradeoff_tortuosus_populations.csv", header = T) 
plot(x =pops$Long , y =pops$Lat )

#preliminary data plotted
prelim_plot <- ggplot(pops, aes(x = Long, y = Lat, 
    colour = Visited)) +
    geom_point()
prelim_plot

#get world map data
world <- getMap(resolution = "low")
world

# Make a vector of country names
usa <- c("United States of America")

# Call the vector in `borders()`
world_usa <- world[world@data$ADMIN %in% usa, ]

with_world_usa <- ggplot() +
  geom_polygon(data = world_usa, 
               aes(x = long, y = lat, group = group),
               fill = NA, colour = "black") + 
  geom_point(data = pops,  # Add and plot species data
             aes(x = Long, y = Lat, 
                 colour = Visited)) +
  coord_quickmap() +  # Prevents stretching when resizing
  theme_classic() +  # Remove ugly grey background
  xlab("Longitude") +
  ylab("Latitude") + 
  guides(colour=guide_legend(title="Visited?"))
with_world_usa

#too big, need california and a few other states

MainStates <- map_data("state")

#plot all states with ggplot2, using black borders and light blue fill
ggplot() + 
  geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="black", fill="lightblue" ) +
  geom_point(data = pops,  # Add and plot species data
           aes(x = Long, y = Lat, 
               colour = Visited)) +
  coord_quickmap() + # Prevents stretching when resizing
  theme_classic() +  # Remove ugly grey background
  xlab("Longitude") +
  ylab("Latitude") + 
  guides(colour=guide_legend(title="Visited?")) + 
  coord_cartesian(xlim=c(-130, -115), ylim = c(34, 43))
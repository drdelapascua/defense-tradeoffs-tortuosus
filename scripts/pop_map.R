library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(maps)
library(rworldmap)

setwd("~/GitHub/defense-tradeoffs-tortuosus")

#### Map of populations for defense trade-off project

#load points
pops <- read.csv("defense_tradeoff_tortuosus_populations.csv", header = T) 

#points for herb drought
pops <- read.csv("herb_drought_populations.csv")

plot(x =pops$long , y =pops$lat )

# need california and a few other states

MainStates <- map_data("state")

#plot all states with ggplot2, using black borders and light blue fill

map <- ggplot() + 
  geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="black", fill="lightblue" ) +
  geom_point(data = pops,  # Add and plot species data
             aes(x = long, y = lat,
               colour = elevation), size = 3) +
  coord_quickmap() + # Prevents stretching when resizing
  theme_classic() +  # Remove ugly grey background
  xlab("Longitude") +
  ylab("Latitude") + 
  coord_cartesian(xlim=c(-125, -117.5), ylim = c(35.5, 43)) +
  scale_color_gradient(low="orange", high="blue")


map

map + 
  geom_label_repel(aes(label = pops$site_abrev, x=pops$long, y=pops$lat),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50')

head(pops)

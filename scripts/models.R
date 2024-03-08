### Danielle De La Pascua
### Ch2 Testing Intraspecifc Defense Trade-offs
### 1-28-24

### Linear Mixed Models

#libraries
library(lme4)
library(nlme)
library(dplyr)
library(tidyr)

#pull data
data <- read.csv("~/GitHub/defense-tradeoffs-tortuosus/data/dl-induced.csv")


#Question: do growth and defense traits trade off in S. tortuosus?
q1 <- lmer(biomass ~ X3MSO_5.2 + OH.Alkenyl_6 + X4MSO_7.1 + X5MSO_10.2 + MSOO_13.8 + OH.I3M_15.1 + X4MT._15.5 + Flavonol_16.1 + I3M_16.7 + Flavonol_17.5 + Flavonol_18.5 + Indole_18.8 + (1 | Elevation),data = data,  na.action = na.exclude)

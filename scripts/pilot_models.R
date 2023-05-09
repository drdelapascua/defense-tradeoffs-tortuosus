#Danielle De La Pascua
#Pilot tissue analysis
#5-3-23

#libraries
library(ggplot2)
library(tidyr)
library(dplyr)

#data
dat <- read.csv("data/pilot_GSL_data_clean.csv")

#data visualization


#linear models

# 3MSO
m1 <- lm(X3MSO ~ trt + Population + sample_time, data = dat)
summary(m1)

# 4MSO
lm(X4MSO ~ trt + Population + sample_time, data = dat )

#Allyl
m2 <- lm(Allyl ~ trt + Population + sample_time, data = dat)
summary(m2)

#but-3-enyl
lm(but.3.enyl ~ trt + Population + sample_time, data = dat)


#X3C Modified
lm(X3C.modified ~ trt + Population + sample_time, data = dat)

#X8MSOO
lm(X8MSOO ~ trt + Population + sample_time, data = dat)

#indol
lm(indol ~ trt + Population + sample_time, data = dat)

#X4C modified
lm(X4C.modified ~ trt + Population + sample_time, data = dat)

#flavanol sulfate
lm(flavonol.sulfate.or.f.gsl ~ trt + Population + sample_time, data = dat)

#indol 1
lm(indol.1 ~ trt + Population + sample_time, data = dat)

#flavonol
lm(flavonol ~ trt + Population + sample_time, data = dat)

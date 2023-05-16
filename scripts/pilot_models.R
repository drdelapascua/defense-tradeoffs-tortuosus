#Danielle De La Pascua
#Pilot tissue analysis
#5-3-23

#libraries ----
library(ggplot2)
library(tidyr)
library(dplyr)

#data ----
dat <- read.csv("data/pilot_GSL_data_clean.csv")

# run sample time as a factor
dat$sample_time <- as.factor(dat$sample_time)

#data visualization ----

#3MSO
plot.3MSO <- ggplot(data = dat, aes(x = sample_time, y = X3MSO)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.3MSO

#4MSO
plot.4MSO <- ggplot(data = dat, aes(x = sample_time, y = X4MSO)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.4MSO

#Allyl
plot.allyl <- ggplot(data = dat, aes(x = sample_time, y = Allyl)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.allyl

#but-3-enyl 
plot.but.3.enyl <- ggplot(data = dat, aes(x = sample_time, y = but.3.enyl)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.but.3.enyl

#3c modified
plot.3C.modified <- ggplot(data = dat, aes(x = sample_time, y = X3C.modified)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.3C.modified

#X8MSOO
plot.8MSOO <- ggplot(data = dat, aes(x = sample_time, y = X8MSOO)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.8MSOO

#indol
plot.indol <- ggplot(data = dat, aes(x = sample_time, y = indol)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.indol

#4C modified
plot.4C.modified <- ggplot(data = dat, aes(x = sample_time, y = X4C.modified)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.4C.modified

#flavonol.sulfate.or.f.gsl
plot.f.gsl <- ggplot(data = dat, aes(x = sample_time, y = flavonol.sulfate.or.f.gsl)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.f.gsl

#indol.1
plot.indol.1 <- ggplot(data = dat, aes(x = sample_time, y = indol.1)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.indol.1

#flavonol
plot.flavonol <- ggplot(data = dat, aes(x = sample_time, y = flavonol)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.flavonol

#Gen-cov aka pop differences (egentic covariance calculated by allyl/sum of allyl and butenyl)
plot.gencov <- ggplot(data = dat, aes(x = sample_time, y = Gen_cov)) +
  geom_boxplot(aes(fill = trt)) +
  facet_wrap(~Population)
plot.gencov
#Ben Hur has more allyls, table moutian has more butanyl

#linear models ----

# > simple linear models ----

# 3MSO
m1 <- lm(X3MSO ~ trt + Population + sample_time, data = dat)
summary(m1)


# 4MSO
m2 <- lm(X4MSO ~ trt + Population + sample_time, data = dat )
summary(m2)

#Allyl
m3 <- lm(Allyl ~ trt + Population + sample_time, data = dat)
summary(m3)

#but-3-enyl
m4 <- lm(but.3.enyl ~ trt + Population + sample_time, data = dat)
summary(m4)

# 3C Modified
m5 <- lm(X3C.modified ~ trt + Population + sample_time, data = dat)
summary(m5)

# 8MSOO
m6 <- lm(X8MSOO ~ trt + Population + sample_time, data = dat)
summary(m6)

#indol
m7 <- lm(indol ~ trt + Population + sample_time, data = dat)
summary(m7)

# 4C modified
m8 <- lm(X4C.modified ~ trt + Population + sample_time, data = dat)
summary(m8)

#flavanol sulfate
m9 <- lm(flavonol.sulfate.or.f.gsl ~ trt + Population + sample_time, data = dat)
summary(m9)

#indol 1
m10 <- lm(indol.1 ~ trt + Population + sample_time, data = dat)
summary(m10)

#flavonol
m11 <- lm(flavonol ~ trt + Population + sample_time, data = dat)
summary(m11)


# > only table mountain ----
TM2dat <- read.csv("data/pilot_GSL_data_TM2.csv")

# 3MSO
TM2m1 <- lm(X3MSO ~ trt + sample_time + trt*sample_time , data = TM2dat)
summary(TM2m1)

# 4MSO
TM2m2 <- lm(X4MSO ~ trt + sample_time + trt*sample_time, data = TM2dat )
summary(TM2m2)

#Allyl
TM2m3 <- lm(Allyl ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m3)

#but-3-enyl
TM2m4 <- lm(but.3.enyl ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m4)

# 3C Modified
TM2m5 <- lm(X3C.modified ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m5)

# 8MSOO
TM2m6 <- lm(X8MSOO ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m6)

#indol
TM2m7 <- lm(indol ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m7)

# 4C modified
TM2m8 <- lm(X4C.modified ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m8) 

#flavanol sulfate
TM2m9 <- lm(flavonol.sulfate.or.f.gsl ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m9)

#indol 1
TM2m10 <- lm(indol.1 ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m10)

#flavonol
TM2m11 <- lm(flavonol ~ trt + sample_time + trt*sample_time, data = TM2dat)
summary(TM2m11)

# > only Ben Hur ----

BHdat <- read.csv("data/pilot_GSL_data_BH.csv")

# 3MSO
BHm1 <- lm(X3MSO ~ trt + sample_time + trt*sample_time , data = BHdat)
summary(BHm1)

# 4MSO
BHm2 <- lm(X4MSO ~ trt + sample_time + trt*sample_time, data = BHdat )
summary(BHm2)

#Allyl
BHm3 <- lm(Allyl ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm3)

#but-3-enyl
BHm4 <- lm(but.3.enyl ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm4)

# 3C Modified
BHm5 <- lm(X3C.modified ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm5)

# 8MSOO
BHm6 <- lm(X8MSOO ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm6)

#indol
BHm7 <- lm(indol ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm7)

# 4C modified
BHm8 <- lm(X4C.modified ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm8) 

#flavanol sulfate
BHm9 <- lm(flavonol.sulfate.or.f.gsl ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm9)

#indol 1
BHm10 <- lm(indol.1 ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm10)

#flavonol
BHm11 <- lm(flavonol ~ trt + sample_time + trt*sample_time, data = BHdat)
summary(BHm11)

# > big model with everything ---- 

# 3MSO
bigm1 <- lm(X3MSO ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm1)

#4MSO
bigm2 <- lm(X4MSO ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm2)

#allyl
bigm3 <- lm(Allyl ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm3)

#but-3-enyl
bigm4 <- lm(but.3.enyl ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm4)

# 3C modified
bigm5 <- lm(X3C.modified ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm5)

# 8MSOO
bigm6 <- lm(X8MSOO ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm6)

# indol
bigm7 <- lm(indol ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm7)

# 4C modified
bigm8 <- lm(X4C.modified ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm8)

# favonol sulfates
bigm9 <- lm(flavonol.sulfate.or.f.gsl ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm9)

# indol 1
bigm10 <- lm(indol.1 ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm10)

# flavonol
bigm11 <- lm(flavonol ~ trt + Population + sample_time + trt*Population + trt*sample_time + sample_time*Population + trt*Population*sample_time, data = dat)
summary(bigm11)

# posthoc tests ----

# use EM means package to get the estimates at the different - estimated marginal means - these are estimates given other factors given the model. gives estimates marginal means, confidence limits etc. 
# emmeans also can do posthoc comparissons to look at estimate of differences. know which are highest, which are different (post hoc!)

#emmeans(bigm1, pairwise ~)

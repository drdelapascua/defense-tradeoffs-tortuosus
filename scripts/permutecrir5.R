# Randomization test for non-random association between measures of constitutive and induced resistance

# Uses two measures of association:
# 1) correlation between induced resistance (damaged minus control family means) and constitutive 
#    resistance (control family mean)
# 2) slope of a regression of damaged vs. control family means

# See Morris, W.F., et al. Oikos (2006) 112: 102-110.

###########################  DATA FILE REQUIREMENTS  ###################################
# 														                                                         #
# The data file must be a csv (comma-separated values) file with 4 columns,		         #
# all of which must contain numeric data:								                               #
# 1) family, genotype, or species index, going from 1 to the number of groups;	       #
#		numbers should not be skipped								                                       #
# 2) treatment: 1 = control 2 = damaged treatment						                           #
# 3) replicate plant number										                                         #
# 4) measure of resistance (assumed to be a non-negative number)				               #
#														                                                           #	
# Column headings in datafile must be "fam", "trt", "rep", "res" (without the quotes)  #
#														                                                           #
########################################################################################

# libararies
library(tidyverse)

# prep script for above format

#total gsl
GSL_totals <-  read.csv("./data/mf_means.csv") %>%
  select(fam = "Population", trt = "treatment",  rep = "mf", res = "totalGSL") %>%
  # filter out pop with no replication
  filter(!fam %in% c("MtSH", "YO10")) %>% 
  mutate(fam = as.numeric(as.integer(as.factor(fam)))) %>%
  # make C = 1 , CW = 2 
  mutate(trt = case_when(
    trt == "C" ~ 1,
    trt == "CW" ~ 2,
  )) %>%
  # make each unique mf its own number chronologically
  group_by(fam) %>%
  mutate(rep_numeric = as.integer(factor(rep))) %>%
  ungroup() %>%
  mutate(rep = as.numeric(rep_numeric)) %>%
  select(-rep_numeric)

str(GSL_totals)
  
#total aliphatic
aliphatic_totals <-  read.csv("./data/mf_means.csv") %>%
  select(fam = "Population", trt = "treatment",  rep = "mf", res = "totalaliphatic") %>%
  # filter out pop with no replication
  filter(!fam %in% c("MtSH", "YO10")) %>% 
  mutate(fam = as.numeric(as.integer(as.factor(fam)))) %>%
  # make C = 1 , CW = 2 
  mutate(trt = case_when(
    trt == "C" ~ 1,
    trt == "CW" ~ 2,
  )) %>%
  # make each unique mf its own number chronologically
  group_by(fam) %>%
  mutate(rep_numeric = as.integer(factor(rep))) %>%
  ungroup() %>%
  mutate(rep = as.numeric(rep_numeric)) %>%
  select(-rep_numeric)

str(aliphatic_totals)

#total indole
indole_totals <-  read.csv("./data/mf_means.csv") %>%
  select(fam = "Population", trt = "treatment",  rep = "mf", res = "totalindole") %>%
  # filter out pop with no replication
  filter(!fam %in% c("MtSH", "YO10")) %>% 
  mutate(fam = as.numeric(as.integer(as.factor(fam)))) %>%
  # make C = 1 , CW = 2 
  mutate(trt = case_when(
    trt == "C" ~ 1,
    trt == "CW" ~ 2,
  )) %>%
  # make each unique mf its own number chronologically
  group_by(fam) %>%
  mutate(rep_numeric = as.integer(factor(rep))) %>%
  ungroup() %>%
  mutate(rep = as.numeric(rep_numeric)) %>%
  select(-rep_numeric)

str(indole_totals)

# This script has a toggle (see "gonormal" below) for two options in cases
# in which negative resistances are likely to arise in the randomization
# procedure.

# Written by W.F. Morris
# Last revised: 6 Dec 2009

#rm(list=ls(all=TRUE)) # clear memory
#graphics.off() # close open figures

# *************************** USER DEFINED PARAMETERS *******************************

# number of bootstrap replicates to perform
maxreps=10000  

# number of families, genotypes, or species in the data set
numfams= 19     

# Chose one of the following options for what to do if, for a family, Imean, the randomly 
# generated mean of induced resistance, is negative and less in absolute value than Cmin, 
# the lowest randomly generated constitutive resistance (which would frequently lead to 
# negative resistances in the randomization):
#   gonormal = 1 -> generate normally distributed induced resistances but
#       discard values less than -Cmin
#   gonormal = 0 -> try generating new family mean constitutive and induced
#       resistances until Imean > -Cmin; this will be slower.
gonormal=1 


# ***********************************************************************************

# set some parameters
q=numfams-1
tcrit=qt(0.975,numfams)           # critical value of t statistic with numfams degrees of freedom
chi2critlo=qchisq(.975,numfams-1) # critical lower value of chi-square statistic with numfams-1 degrees of freedom     
chi2critup=qchisq(.025,numfams-1) # critical upper value of chi-square statistic with numfams-1 degrees of freedom   

data = GSL_totals
attach(data)  # column headings in datafile must be fam, trt, rep, res

Nc=Nd=Cest=Dest=Cc=Cd=matrix(0,numfams,1) # making matrices for each variable (see below for what they are)

# Calculate population-level means and coeffs of variation of each treatment
for (f in 1:numfams) {
  # loop through each population,   
    I=which(fam==f & trt==1) 	# find control plants for fam f
    Nc[f]=length(I)           # Nc = number of reps in each population in the control trt
    x=res[I]
    Cest[f]=mean(x)           # Cest = population mean resistances in control plants
    Cc[f]=sd(x)/Cest[f]      	# Cc =  population coefficient of variation in control plants 

    I=which(fam==f & trt==2) 	# repeat for damage trt.
    Nd[f]=length(I)           # Nd = number of reps in each population in the damage trt
    x=res[I]
    Dest[f]=mean(x)           # Dest = population mean resistance in damaged plants
    Cd[f]=sd(x)/Dest[f]       # Cd = population coefficient of variation in damaged plants
    
}

# Calculate other parameters
Iest=Dest-Cest      # difference measure of induced resistance
CVs=(Nc*Cc+Nd*Cd)/(Nc+Nd)     # weighted mean CV across trts., by pop
CV=mean(CVs)                  # mean weighted CV across pops 
sigma2=log(1+CV^2)                  
sigma=sigma2^.5               # SD of exponent in lognormal resistance 

# Plot the data to make sure it looks correct

split.screen( matrix(c(.25,.75,.5,1,.25,.75,0,.5),2,4,byrow=TRUE) )

screen(1)
xx=seq(0,2,.01)
plot( Cest,Iest,type="p",
	xlab="Constitutive Resistance",
	ylab="Induced Resistance",
	xlim=c(min(Cest)-.01,max(Cest)+.01),
	ylim=c(min(Iest)-.01,max(Iest)+.01)
)
lines(c(min(Cest)-.01,max(Cest)+.01),c(0,0),lty="dashed")

screen(2)
xx=seq(min(min(Dest),min(Cest)),max(max(Dest),max(Cest)),.1)
plot(Cest,Dest,type="p",
	xlab="Control Resistance",
	ylab="Damaged Resistance",
	xlim=c(min(Cest)-.01,max(Cest)+.01),
	ylim=c(min(Dest)-.01,max(Dest)+.01)
)
lines(xx,xx)

# compute confidence interval for species mean and among-pop variance in constitutive resistance
Cmeanhat=mean(Cest)   
Cvarhat=var(Cest)
SE.Cmean=(Cvarhat/numfams)^.5
CI.Cmean=c(Cmeanhat-tcrit*SE.Cmean, Cmeanhat+tcrit*SE.Cmean)
CI.Cvar=c(Cvarhat*q/chi2critlo, Cvarhat*q/chi2critup)

# compute confidence interval for species mean and among-pop variance in induced resistance
Imeanhat=mean(Iest) # uses Iest, which is the difference from Damaged - Control         
Ivarhat=var(Iest)
SE.Imean=(Ivarhat/numfams)^.5
CI.Imean=c(Imeanhat-tcrit*SE.Imean, Imeanhat+tcrit*SE.Imean)
CI.Ivar=c(Ivarhat*q/chi2critlo, Ivarhat*q/chi2critup)

# observed tradeoff measures
corrobs=cor(Iest,Cest)          # correlation measure of trade-off
bobs=cov(Dest,Cest)/var(Cest)   # slope measure of trade-off

print("Observed correlation between induced (damage minus control)and constitutive (control) resistance",quote=FALSE)
corrobs
print("Observed slope of regression of damaged vs. control fam means",quote=FALSE)
bobs

correst=best=matrix(0,1,maxreps)

# generate pseudodata with no correlation between C and I  using Monte Carlo sims

rep=0

while(rep<maxreps) {
	
	done=0
    
    # generate grand pop mean and var for C and I using their confidence intervals 
    
    Cmean=Inf  
    while( Cmean<CI.Cmean[1]||Cmean>CI.Cmean[2] ) {
        Cmean=Cmeanhat+SE.Cmean*rnorm(1)
    }
	
	Cvar=Inf  
    while( Cvar<CI.Cvar[1]||Cvar>CI.Cvar[2] ) {
     	Cvar=Cvarhat*rchisq(1,df=q)/q
    }
    
    Imean=Inf  
    while( Imean<CI.Imean[1]||Imean>CI.Imean[2] ) {
     	Imean=Imeanhat+SE.Imean*rnorm(1)
    }

 	Ivar=Inf  
 	while( Ivar<CI.Ivar[1]||Ivar>CI.Ivar[2] ) {
     		Ivar=Ivarhat*rchisq(1,df=q)/q
    }        
 	Isd=Ivar^.5
           
 	# generate (lognormal) resistance measures for each pop in control and damage treatments
 	sig2=log(1+Cvar/Cmean^2)
 	C=exp( rnorm(numfams,log(Cmean)-0.5*sig2,sig2^.5) ) # generate lognormal C's                 
 	Cmin=min(C)
 	if(Imean>-Cmin) {  # abs(Imean)>Cmin, so mean of the following lognormal is positive
		done=1
      	rep=rep+1
      	sig2=log(1+Ivar/(Imean+Cmin)^2)
      	I=exp( rnorm(numfams,log(Imean+Cmin)-0.5*sig2,sig2^.5) ) - Cmin  # generate I's as shifted lognormal, with Imin=-Cmin  
 	} else {
		if(gonormal) {  # Imean<=-Cmin, so mean of preceding lognormal would be neg; instead draw normal I's and discard if < -Cmin
        		done=1
        		rep=rep+1
        		nf=0
        		while(nf<numfams) {
            		Itemp=rnorm(1,Imean,Isd)
            		if(Itemp>-Cmin) {
                			nf=nf+1
                			I[nf]=Itemp
            		}
			}
		}
 	}
 	# NOTE: The preceding "else" statement will cause the mean and
 	# variance of I to differ from Imean and Ivar
        
    if(done) {

		D=C+I
		Cest=Dest=matrix(0,numfams,1)
		for(f in 1:numfams) { # generate lognormal resistance measures for each plant, and take means 
            	Cest[f]=mean( exp( rnorm( Nc[f],log(C[f])-.5*sigma2,sigma ) ) ) # use of sigma2 and sigma assumes fixed CV across families and trts.
            	Dest[f]=mean( exp( rnorm( Nd[f],log(D[f])-.5*sigma2,sigma ) ) ) 
        	}
        	Iest=Dest-Cest
        	correst[rep]=cor(Iest,Cest)
		best[rep]=cov(Dest,Cest)/var(Cest) # slope of regression of Dest on Cest

	} 
} # while(rep<maxreps)

corrlow=quantile(correst,.05)
slopelow=quantile(best,.05)

print("Lower 5th percentile of the bootstrap distribution of the correlation coefficient",quote=FALSE)
corrlow
print("Lower 5th percentile of the bootstrap distribution of the regression slope",quote=FALSE)
slopelow

# for plotting, set lower limit at largest value divisible by .2 below estimates
ccmin=round(min(correst),1)
if(ccmin>min(correst)) ccmin=ccmin-.1
if( !( ccmin %% .2 < 1e-10) ) ccmin=ccmin-.1 
# for plotting, set upper limit at smallest value divisible by .2 above estimates
ccmax=round(max(correst),1)
if(ccmax<max(correst)) ccmax=ccmax+.1
if( !(ccmax %% .2 < 1e-10) ) ccmax=ccmax+.1 
dcc=(ccmax-ccmin)/50
edgecc=seq(ccmin,ccmax,dcc)
cctks=seq(ccmin,ccmax,.2)

yy=hist(correst,breaks=edgecc,plot=FALSE)
yc=yy$counts
lowcc=min(which(yc>0))
highcc=max(which(yc>0))
yccmax=1.1*max(yc)
Fcc=cumsum(yc)/maxreps
cci=min(which(Fcc>=0.05))    # index of critical correlation coeff.
cco=min(which(edgecc>=as.numeric(corrobs)))   # index of observed correlation coeff.
pccobs=Fcc[cco]
print("Approximate probability of observed or smaller correlation",quote=FALSE)
pccobs

bmin=floor(min(best))
bmax=ceiling(max(best))
db=(bmax-bmin)/50
edgeb=seq(bmin,bmax,db)
btks=seq(bmin,bmax,1)

yy=hist(best,breaks=edgeb,plot=FALSE)
lowb=min(which(yy$counts>0))
highb=max(which(yy$counts>0))
ybmax=1.1*max(yy$counts)
Fb=cumsum(yy$counts)/maxreps
bi=min(which(Fb>=0.05))    # index of critical regression slope
bo=min(which(edgeb>=as.numeric(bobs)))   # index of observed regression slope
pbobs=Fb[bo]
print("Approximate probability of observed or smaller regression slope",quote=FALSE)
pbobs


# Plot the bootstrap distribution of the correlation coefficient;
# plot the lower 5th percentile as a vertical dashed line
# plot the observed value as a vertical solid line

windows()

hist(correst,breaks=edgecc,ylab="Number of replicates",xlab="Correlation coefficient",
	main="Solid:Observed; Dashed:CumProb=0.05",
	xlim=c(ccmin,ccmax),ylim=c(0,yccmax),lwd=1.5,ps=16,cex.main=1,xaxt = "n")
axis(1,cctks)
lines(edgecc[cci]*c(1,1),c(0,yccmax),lwd=3,lty="dashed")
lines(edgecc[cco]*c(1,1),c(0,yccmax),lwd=3)


# Plot the bootstrap distributions of the correlation coefficient and regression slope,
# as well as the cumulative distributions of the two measures;
# plot the lower 5th percentile as a vertical dashed line;
# plot the observed value as a vertical solid line

windows()

split.screen(c(2,2))

screen(1)
hist(correst,breaks=edgecc,ylab="Number of replicates",xlab="Correlation coefficient",
	main="Solid: Observed; Dashed: CumProb=0.05",
	xlim=c(ccmin,ccmax),ylim=c(0,yccmax),lwd=1.5,ps=16,cex.main=.9,xaxt = "n")
axis(1,cctks)
lines(edgecc[cci]*c(1,1),c(0,yccmax),lwd=3,lty="dashed")
lines(edgecc[cco]*c(1,1),c(0,yccmax),lwd=3)
   
screen(2)
hist(best,breaks=edgeb,ylab="Number of replicates",xlab="Regression slope",
	main="Solid: Observed; Dashed: CumProb=0.05",
	xlim=c(bmin,bmax),ylim=c(0,ybmax),lwd=1.5,ps=16,cex.main=.9,
	xaxp=c(bmin,bmax,bmax-bmin)	)
lines(edgeb[bi]*c(1,1),c(0,ybmax),lwd=3,lty="dashed")
lines(edgeb[bo]*c(1,1),c(0,ybmax),lwd=3)

screen(3)
plot(edgecc[-length(edgecc)],Fcc,type="l",ylab="Cumulative probability",xlab="Correlation coefficient",
	main="Horiz. line shows Pr(observed correlation)",
	xlim=c(ccmin,ccmax),ylim=c(0,1),lwd=1.5,ps=16,cex.main=.9,xaxt = "n")
axis(1,cctks)
lines(edgecc[cco]*c(1,1),c(0,Fcc[cco]),lwd=2)
lines(c(ccmin,edgecc[cco]),Fcc[cco]*c(1,1),lwd=2)

screen(4)
plot(edgeb[-length(edgeb)],Fb,type="l",ylab="Cumulative probability",xlab="Regression slope",
	main="Horiz. line shows Pr(observed slope)",
	xlim=c(bmin,bmax),ylim=c(0,1),lwd=1.5,ps=16,cex.main=.9,xaxt = "n")
axis(1,btks)
lines(edgeb[bo]*c(1,1),c(0,Fb[bo]),lwd=2)
lines(c(bmin,edgeb[bo]),Fb[bo]*c(1,1),lwd=2)
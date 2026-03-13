#*******************************************************************************************
#**                         II. HYPOTHESIS TESTING in R                                   **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#*******************************************************************************************


setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data")
dir()

library(lattice)
library(RColorBrewer)

# Loading data (HUD census tract data on vacancies, Q3 2010)
     
census <-read.csv("HUDUser/DA2_Census_TractsLarge.csv", h=T, na.strings="#N/A")
summary(census)

length(census$res_vac)
length(census$res_vac[census$metro>0])

hist(census$metro)
table(census$metro)

densityplot(census$res_vac[census$metro==1])


# Using factors to visualise different subsamples
fac         <- factor(census$metro, labels=c("Nation","DeWaLi", "Ann Arbor MSA", "Flint MSA", "Monroe MSA"))
summary(fac)

fac.int     <- as.integer(fac)
summary(fac.int)

nfac <- max(fac.int)

# ColorBrewer Palettes (http://colorbrewer2.org/)
display.brewer.all(n=11, exact.n=FALSE)

plotclr     <- brewer.pal(nfac,"Spectral")

boxplot(census$res_vac ~ fac, horizontal = T, col=plotclr, ylab="Cluster", xlab="Number of vacant residential properties", cex.lab=1 ,cex.axis=0.8)

sub.census  <- subset(census, census$metro>0) 
fac.1       <- factor(sub.census$metro, labels=c("DeWaLi", "Ann Arbor MSA", "Flint MSA", "Monroe MSA"))
nfac.1      <- max(as.integer(fac.1))
plotclr.1   <- brewer.pal(nfac.1,"Spectral")

boxplot(sub.census$res_vac ~ fac.1, horizontal = T, col=plotclr.1, ylab="Cluster", xlab="Number of vacant residential properties", cex.lab=1 ,cex.axis=0.8)

# Hypothesis Tesing
mu.all      <- mean(census$res_vac)
mu.da2      <- mean(census$res_vac[census$metro>0])
mu.dwl      <- mean(census$res_vac[census$metro==1])
mu.a2       <- mean(census$res_vac[census$metro==2])
mu.flint    <- mean(census$res_vac[census$metro==3])

sd.all      <- sd(census$res_vac)
sd.da2      <- sd(census$res_vac[census$metro>0])
sd.dwl      <- sd(census$res_vac[census$metro==1])
sd.a2       <- sd(census$res_vac[census$metro==2])
sd.flint    <- sd(census$res_vac[census$metro==3])

mu.all
mu.da2
mu.dwl
mu.a2
mu.flint

t.test(census$res_vac[census$metro==1], mu=mu.all, alternative=c("two.sided"))

t.test(census$res_vac[census$metro==1], mu=mu.dwl, alternative=c("two.sided"))
t.test(census$res_vac[census$metro==1], mu=mu.flint, alternative=c("two.sided"))
t.test(census$res_vac[census$metro==1], mu=mu.da2, alternative=c("two.sided"))

t.result <- t.test(census$res_vac[census$metro==2], mu=mu.dwl, alternative=c("two.sided"))
summary(t.result)

t.test(census$res_vac[census$metro==2], mu=mu.dwl, alternative=c("two.sided"))

1 - pnorm(1.96,mean=0,sd=1)

2*(1 - pnorm(1.96,0,1))
1 - pf(3.6,4,43)
1 - pt(0.8328,1287)

d.all       <- density(sub.census$res_vac)
d.dwl       <- density(sub.census$res_vac[sub.census$metro==1])
d.a2        <- density(sub.census$res_vac[sub.census$metro==2])

plot(d.all)
lines(d.dwl, col="red")
lines(d.a2, col="blue")

hist(sub.census$bus_vac[sub.census$bus_vac<200], breaks=55)
hist(sub.census$res_vac[sub.census$res_vac<200], breaks=55, add=T, col="gray")

# Replacing values of variables
summary(sub.census$structure_age)

is.na(sub.census$structure_age) <- sub.census$structure_age < 0
is.na(sub.census$hh_median)     <- sub.census$hh_median < 0
is.na(sub.census$structure_rooms) <- sub.census$structure_rooms < 0
is.na(sub.census$structure_value) <- sub.census$structure_value < 0

summary(sub.census$structure_age)
# do NOT use census$structure_age <- replace(census$structure_age,census$structure_age<0,0) [why?]

boxplot(sub.census$structure_age ~ fac.1, horizontal = T, col=plotclr.1, ylab="Cluster", xlab="Median house age", cex.lab=1 ,cex.axis=0.8)

sub.census$val <- sub.census$structure_value / 1000
boxplot(sub.census$val ~ fac.1, horizontal = T, col=plotclr.1, ylab="Cluster", xlab="House value ($000s)", cex.lab=1 ,cex.axis=0.8)

boxplot(sub.census$hh_median ~ fac.1, horizontal = T, col=plotclr.1, ylab="Cluster", xlab="Median annual household income", cex.lab=1 ,cex.axis=0.8)

library(corrgram)
corrgram(census[,2:4], order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt,
  main="Census tract-level vacancies")


# The importance of data visualisation (Testing distributional properties of observational data)

hist(census$res_vac[census$metro==1], freq=F, breaks=25)

nobs             <- length(census$res_vac[census$metro==1])

sample.norm      <- rnorm(nobs, mean=mu.dwl, sd=sd.dwl)
sample.weibull   <- rweibull(nobs, shape=0.9, scale=mu.dwl)  # Weibull shape parameter indicates rate of change of failure rate. k=1 constant failure rate.

lines(density(sample.norm), col="red", lwd=2)
lines(density(sample.weibull), col="blue", lwd=2)


# Standardisation induces normality even in finite samples (Levy-Lindberg Central Limit Theorem)
z.res_vac        <- (census$res_vac[census$metro==1]-mean(census$res_vac[census$metro==1],na.rm=T))/sd(census$res_vac[census$metro==1],na.rm=T)
hist(z.res_vac, freq=F)
curve(dnorm, add=T, col="red", lwd="2")

#Linear hypothesis testing (using  , see also http://www.youtube.com/watch?v=pwZFqdTPcB4)
library(car)

linearHypothesis(lm.object, restriction)

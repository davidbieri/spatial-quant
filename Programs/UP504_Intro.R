#*******************************************************************************************
#**                         I. INTRO TO R: BASIC CONCEPTS                                 **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#*******************************************************************************************


setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data")

dir()

## Clear the workspace

rm(list=ls())

## Chapter 1 ISwR -- Basics

library(ISwR)

height    <- c(1.75, 1.80, 1.65, 1.90, 1.74, 1.94, 1.81, 1.53, 1.65, 1.67)
weight    <- c(60, 72, 57, 90, 95, 83, 80, 49, 66, 62)
sex       <- c("Female", "Female", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Female")

data.entry(height)
fix(height)

demo      <- cbind(height, weight, sex)
demo.df   <- data.frame(height, weight, sex)

summary(demo)
summary(demo.df)

names(demo)
names(demo.df)

head(demo)
head(demo.df)

attributes(demo)
attributes(demo.df)

class(demo)
class(demo.df)


plot(weight)
plot(weight,height)

plot(demo)
plot(demo.df)

plot(demo.df$weight[demo.df$sex=="Male"], demo.df$height[demo.df$sex=="Male"], pch="M", col="red")
points(demo.df$weight[demo.df$sex=="Female"], demo.df$height[demo.df$sex=="Female"], pch="F", col="blue")


X         <- rnorm(30, 75, 15.44)           # Demonstrate Law of Large Numbers by changing N=3,000; N=300,000)
plot(X)
hist(X)
summary(X)

x         <- as.data.frame(X)
fix(X)

# Random number generation and visualisation
x         <- rnorm(10000)
mean(x)
plot(x)

y         <- dnorm(x)
plot(x,y)

Fn        <- ecdf(x)
plot(knots(Fn),Fn(knots(Fn)))


# Loading data (HUD census tract data on vacancies, Q3 2010)
library(lattice)
library(RColorBrewer)

census          <- read.csv("HUDUser/DA2_Census_TractsLarge.csv", h=T, na.strings="#N/A")

summary(census)
head(census)
names(census)

length(census$res_vac)
length(census$res_vac[census$metro>0])

densityplot(census$res_vac[census$metro==1])

# Factoring the "metro" indicator variable to use as levels a in boxplot
fac             <- factor(census$metro, labels=c("Nation","DeLiWa MSA", "Ann Arbor MSA", "Flint MSA", "Monroe MSA"))
summary(fac)

fac.int         <- as.integer(fac)
summary(fac.int)
nfac            <- max(as.integer(fac))
plotclr         <- brewer.pal(nfac,"Spectral")  #Selecting the appropriate plotting colours from www.ColorBrewer2.org

boxplot(census$res_vac ~ fac, horizontal = T, col=plotclr, ylab="Cluster", xlab="Number of vacant properties", cex.lab=1 ,cex.axis=0.8)


# Producing the same boxplot as above for the Detroit region only
sub.census      <- subset(census, census$metro>0) 
fac.1           <- factor(sub.census$metro, labels=c("DeLiWa MSA", "Ann Arbor MSA", "Flint MSA", "Monroe MSA"))
nfac.1          <- max(as.integer(fac.1))
plotcl.1        <- brewer.pal(nfac.1,"Spectral")

boxplot(sub.census$res_vac ~ fac.1, horizontal = T, col=plotclr, ylab="Cluster", xlab="Number of vacant properties", cex.lab=1 ,cex.axis=0.8)


# Some simple density plots
d.all <- density(census$res_vac)
d.de <- density(census$res_vac[vac$metro==1])
d.a2 <- density(census$res_vac[vac$metro==2])
d.fl <- density(census$res_vac[vac$metro==3])

plot(d.a2, lwd=2, main="Density plots")
lines(d.de, col="red", lwd=2)
lines(d.fl, col="blue", lwd=2)

hist(census$bus_vac[census$bus_vac<200], breaks=55)
hist(census$res_vac[census$res_vac<200], breaks=55, add=T, col="gray")

library(corrgram)
corrgram(census[,2:4], order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie, text.panel=panel.txt,
         main="Census tract-level vacancies")


#Hypothesis Testing
mu.all <- mean(census$res_vac)
mu.da2 <- mean(census$res_vac[census$metro>0])
mu.dlw <- mean(census$res_vac[census$metro==1])
mu.a2  <- mean(census$res_vac[census$metro==2])

mu.all
mu.da2
mu.dlw
mu.a2

t.test(census$res_vac[census$metro==1], mu=mu.all)

1 - pnorm(1.96,mean=0,sd=1)

2*(1 - pnorm(1.96,0,1))
1 - pf(3.6,4,43)


# Accessing data in other common formats (.csv, .txt. .xls, .dta)
library(XLConnect)
library(lattice)


# Data visualisation
library(foreign)
library(googleVis)          ## interface between R and the Google Visualisation API
library(rgl)                ## 3-D Visualisation

mi              <- read.dta("Earnings_MSA_MI.dta")
names(mi)
MI              <- gvisMotionChart(mi, idvar="area", timevar="time")
plot(MI)
print(MI, file="MI_GoogleVisChart.html")



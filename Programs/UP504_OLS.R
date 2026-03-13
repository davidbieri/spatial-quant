#*******************************************************************************************
#**                   III. INTRODUCTION TO LINEAR REGRESSION in R                         **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#*******************************************************************************************

setwd("C:/Documents and Settings/[your directory]/Data")
dir()

library(lattice)
library(RColorBrewer)


# Loading data (HUD census tract data on vacancies, Q3 2010)
     
census      <- read.csv("HUDUser/DA2_Census_Tracts.csv", h=T, na.strings="#N/A")
fac         <- factor(census$metro, labels=c("DeWaLi", "Ann Arbor MSA", "Flint MSA", "Monroe MSA"))

nfac        <- max(as.integer(fac))
plotclr     <- brewer.pal(nfac,"Spectral")

# ColorBrewer Palettes (http://colorbrewer2.org/)
display.brewer.all(n=11, exact.n=FALSE)


regvars     <- c("rooms","val","age","li_units","qct")
regdata     <- census[regvars]
length(regdata[,1])

# Standard OLS calculation using canned lm function
reg.0       <- lm(census$val ~ census$li_units, data=census)
reg.1       <- lm(census$val ~ census$li_units + census$rooms + census$age, data=census)
reg.2       <- lm(census$val ~ census$li_units + census$rooms + census$age + census$vacrate + census$hh_median + census$edu + census$black + census$ooh + census$assist, data=census)

output      <- summary(reg.0)
output
names(output)
summary(output)

output      <- summary(reg.2)
output

coef(reg.2)
confint(reg.2, level = 0.95)

rss         <- deviance(reg.2)
r2          <- output$r.squared
vb          <- vcov(reg.2)

reg.3       <- lm(census$val ~ census$rooms, data=census)
plot(census$rooms, census$val)
abline(reg.3, col="dark green", lwd=2)


# Quantile regressions
library(quantreg)
library(latticeExtra)

qreg.0       <- rq(census$val ~ census$li_units, tau= 1:49/50, data=census)
qreg.1       <- rq(census$val ~ census$li_units + census$rooms + census$age, tau= 1:49/50, data=census)
qreg.2       <- rq(census$val ~ census$li_units + census$rooms + census$age + census$vacrate + census$hh_median + census$edu + census$black + census$ooh + census$assist, tau= 1:49/50, data=census)

#Plot parameter estimates across quantiles
plot(summary(qreg.1))

# Visualise outliers: http://quickfacts.census.gov/qfd/states/26/26161.html
plot(census$rooms, census$val, pch="")
text(census$rooms, census$val, labels=census$geoid, cex=.6)

mu.qct     <- mean(census$val[census$qct==1], na.rm=T)   
mu.other   <- mean(census$val[census$qct==0], na.rm=T)   

t.result <- t.test(census$val[census$qct==0], y=census$val[census$qct==1], alternative=c("two.sided"))
t.result

fac.1       <- factor(census$qct, labels=c("Not qualified", "LIHTC qualified census tract"))
nfac        <- max(as.integer(fac.1))
boxplot(census$val ~ fac.1, horizontal = T, ylab="Policy treatment", xlab="Median house value ($000s)", cex.lab=1 ,cex.axis=0.9, col=c("green","red"))

reg.resid   <- resid(reg)
reg.stdres  <- rstandard(reg)


#Estimation with FGLS
sig2        <- t(e.hat)%*%e.hat / n
Omega.H0    <- diag(sig2[1,1] * rep(1,n))

M           <- cbind(M, e.hat)
colnames(M)[6] <- 'e.hat'
e.dda       <- M[M[,"dda"]==1,"e.hat"] 
e.norm      <- M[M[,"dda"]==0,"e.hat"]

sig2.dda    <- (t(e.dda)%*%e.dda) / length(e.dda)
sig2.norm   <- (t(e.norm)%*%e.norm) / length(e.norm)

Omega.H1    <- diag(ifelse(M[,"dda"]==1,sig2.dda,sig2.norm))

b.gls1      <- solve(t(X)%*%solve(Omega.H0)%*%X)%*%t(X)%*%solve(Omega.H0)%*%y
b.gls2      <- solve(t(X)%*%solve(Omega.H1)%*%X)%*%t(X)%*%solve(Omega.H1)%*%y

Vb1         <- solve(t(X)%*%solve(Omega.H0)%*%X)
Vb2         <- solve(t(X)%*%solve(Omega.H1)%*%X)

se1         <- sqrt(diag(Vb1))
se2         <- sqrt(diag(Vb2))

t.val1      <- b.gls1/se1
t.val2      <- b.gls2/se2

reg         <- lm(census$val ~ census$structure_rooms + census$age, data=census)
output      <- summary(reg)
names(output)
summary(output)

coef(reg)

rss         <- deviance(reg)
r2          <- output$r.squared
vb          <- vcov(reg)

Omega.H0    <- model.he$residuals^2                        
Omega.H1    <- model.he$residuals^2*(n/(n-k))


qqnorm(reg.stdres, ylab="Standardized Residuals", xlab="Normal Scores", main="")
qqline(reg.stdres)

par(mfrow = c(3, 2))
plot(reg, which=1:6)
par(mfrow = c(1, 1))

# Basic matrix algebra review

a       <- matrix(seq(1,6), ncol=1)
a1      <- matrix(seq(1,6), ncol=2)

i       <- matrix(rep(1,6), ncol=1)

b       <- matrix(seq(1,6), nrow=1)
b1      <- matrix(seq(1,6), nrow=2)

A       <- a%*%t(a)
B       <- b%*%t(b)

sum(a)
t(i)%*%a            # convenient summation of vector elements
t(i)%*%i


a%*%b
b%*%a

SS.a    <- t(a)%*%a

diag(A)
sum(diag(A))

# Centering  of a variable

n           <- length(y)
i           <- matrix(rep(1,n), ncol=1)

Pi          <- i%*%solve(t(i)%*%i)%*%t(i)    # or write directly as (i%*%t(i))/n to save resources
I           <- diag(n)
Mi          <- I - Pi

z           <- scale(y, center = T, scale = F)
summary(z)

z.alt       <- Mi%*%y


Z           <- cbind(z, z.alt)
colnames(Z) <- c("z","z.alt")

dotchart(A, gcol=c("red", "blue"))

summary(z.alt)          


# Plotting coefficients
# http://www.r-statistics.com/2010/07/visualization-of-regression-coefficients-in-r/
## 
## set up the work environment
##
## packages
require(RevoScaleR)
require(coefplot)
## directories
sampleDataDir   <- rxGetOption("sampleDataDir")
workDir         <- getwd()
## get loan default dataset (from RevoScaleR)
mortgagesXdf    <- file.path(workDir, "Mortgages.xdf")
xdfAppend <- "none"
xdfOverwrite    <- TRUE
for (myYear in 2000:2009) {
csvFileName <- file.path(sampleDataDir, paste("mortDefaultSmall", myYear, ".csv", sep=""))
rxTextToXdf(inFile=csvFileName, outFile=mortgagesXdf, append = xdfAppend, overwrite = xdfOverwrite)
xdfAppend <- "rows"
xdfOverwrite <- FALSE
}

rxGetInfoXdf(file=mortgagesXdf, getVarInfo=TRUE)
model1          <- rxLinMod(default ~F(year) + yearsEmploy + ccDebt + creditScore, 
cube=TRUE, data=mortgagesXdf)
summary(model1)
coefplot(model1)
coefplot(model1, factors=c("F_year"))
coefplot(model1, factors=c("F_year"), numeric=TRUE)



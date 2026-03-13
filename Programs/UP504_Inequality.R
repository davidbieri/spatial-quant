#*******************************************************************************************
#**                               VI. ECONOMIC ANALYSIS                                   **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#*******************************************************************************************

setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data")
census      <- read.csv("HUDUser/DA2_Census_Tracts.csv", h=T, na.strings="#N/A")


# The "ineq" Package
# http://cran.r-project.org/web/packages/ineq/ineq.pdf
library(ineq)

# Programming the Gini coefficient as a function (not needed for HW6)
gini            <- function(x, unbiased = TRUE, na.rm = FALSE){
    if (!is.numeric(x)){
        warning("'x' is not numeric; returning NA")
        return(NA)
    }
    if (!na.rm && any(na.ind <- is.na(x)))
        stop("'x' contain NAs")
    if (na.rm)
        x <- x[!na.ind]
    n <- length(x)
    mu <- mean(x)
    N <- if (unbiased) n * (n - 1) else n * n
    ox <- x[order(x)]
    dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
    dsum / (mu * N)
}

gini(c(100,0,0,0))

# Lorenz curves based on census tract median household incomes
inc             <- census[, c("hh_median", "metro", "res_vac", "res_tot", "bus_vac", "bus_tot", "ooh", "black")]

inc$res         <- inc$res_vac/inc$res_tot
inc$bus         <- inc$bus_vac/inc$bus_tot

inc             <- na.omit(inc)

# Preparing for Lorenz Curves (income; .1= residential vacancy rate; .2= business vacancy rate; .3= owner-occupied housing; .4= African-American households)
Lc.d            <- Lc(inc$hh_median[inc$metro==1])
Lc.d1           <- Lc(inc$res[inc$metro==1])
Lc.d2           <- Lc(inc$bus[inc$metro==1])
Lc.d3           <- Lc(inc$ooh[inc$metro==1])
Lc.d4           <- Lc(inc$black[inc$metro==1])

Lc.a            <- Lc(inc$hh_median[inc$metro==2])
Lc.a1           <- Lc(inc$res[inc$metro==2])
Lc.a2           <- Lc(inc$bus[inc$metro==2])
Lc.a3           <- Lc(inc$ooh[inc$metro==2])
Lc.a4           <- Lc(inc$black[inc$metro==2])
 
Lc.f            <- Lc(inc$hh_median[inc$metro==3])
Lc.f1           <- Lc(inc$res[inc$metro==3])
Lc.f2           <- Lc(inc$bus[inc$metro==3])
Lc.f3           <- Lc(inc$ooh[inc$metro==3])
Lc.f4           <- Lc(inc$black[inc$metro==3])

Lc.m            <- Lc(inc$hh_median[inc$metro==4])
Lc.m1           <- Lc(inc$res[inc$metro==4])
Lc.m2           <- Lc(inc$bus[inc$metro==4])
Lc.m3           <- Lc(inc$ooh[inc$metro==4])
Lc.m4           <- Lc(inc$black[inc$metro==4])

#Plot Lorenz curves (graphs LC only for one variable)
plot(Lc.d, main="MSA-level household income inequality")
lines(Lc.a, col=2)
lines(Lc.f, col=3)
lines(Lc.m, col=4)
legend(x=c("topleft"),c("Detroit","Ann Arbor", "Flint", "Monroe"), col=c("black","blue", "green", "red"), 
       lwd=c(2,2,2,2), box.lty=0, cex=0.85, ncol=3, horiz=F)


plot(Lc.d1, main="MSA-level inequality in residential vacancy rates")
lines(Lc.a1, col=2)
lines(Lc.f1, col=3)
lines(Lc.m1, col=4)
legend(x=c("topleft"),c("Detroit","Ann Arbor", "Flint", "Monroe"), col=c("black","blue", "green", "red"), 
        lwd=c(2,2,2,2), box.lty=0, cex=0.85, ncol=3, horiz=F)

# Theil indicies 
# Two equivalent ways of applying a function to a data frame grouped by a condition
by(inc$hh_median, inc$metro, mean)
by(inc$hh_median, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))
tapply(inc$hh_median, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))


tapply(inc$res_vac, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))
tapply(inc$bus_vac, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))
tapply(inc$ooh, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))
tapply(inc$black, inc$metro, function(x) ineq(x, type=c("Theil"), parameter=0))

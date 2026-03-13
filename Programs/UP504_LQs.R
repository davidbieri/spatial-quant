#*******************************************************************************************
#**                               VI. ECONOMIC ANALYSIS                                   **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#*******************************************************************************************

library(foreign)
library(googleVis)          ## interface between R and the Google Visualisation API
library(rgl)                ## 3-D Visualisation
library(Rcmdr)              ## 3-D Visualisation (opens up addtional GUI "R Commander", do not close)
library(RColorBrewer)       ## ColorBrewer implementation for R

setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data")


# Data from County Business Patterns (http://www.census.gov/econ/cbp/)
national        <- read.csv("LQs_CBP_National.csv", h=T, na.strings="#N/A") # Reads in national CBP data for 2-digit NAICS industries  
msa             <- read.csv("LQs_CBP_AnnArbor.csv", h=T, na.strings="#N/A")

summary(national)

# Show total employment for nation and Ann Arbor; NAICS code 0 indicates summary accross all 2-digit categories
national$Empl08[national$NAICS==0]
msa$Empl08[msa$NAICS==0]


#*******************************************************************************************
#**    1) Measures of economic activity and spatial concentration                         **
#*******************************************************************************************

#*******************************************************************************************
#**    1.1) Location quotients                                                            **
#*******************************************************************************************

msa$lq98        <- (msa$Empl98/msa$Empl98[msa$NAICS==0])/(national$Empl98/national$Empl98[national$NAICS==0])
msa$lq08        <- (msa$Empl08/msa$Empl08[msa$NAICS==0])/(national$Empl08/national$Empl08[national$NAICS==0])


#*******************************************************************************************
#**    1.2) Shift-share analysis                                                          **
#*******************************************************************************************
#  Defining the different growth rates; national growth rate (gn), national industry differential (gni),
#  local industry differential (gl)

msa$gn          <- (national$Empl08[msa$NAICS==0]/national$Empl98[national$NAICS==0]) - 1
msa$gni         <- (national$Empl08/national$Empl98) - 1
msa$gl          <- (msa$Empl08/msa$Empl98) - 1

# Shift share components; national share (NS), industry mix (IM), local factors (LF)
msa$ns          <- msa$Empl98*(1 + msa$gn)
msa$im          <- msa$Empl98*(msa$gni - msa$gn)
msa$lf          <- msa$Empl98*(msa$gl - msa$gni)

msa$ns + msa$im + msa$lf
(msa$ns + msa$im + msa$lf)[1]   # Shift-share total for Ann Arbor in 2008 (should be equal to msa$Empl08[msa$NAICS==0])

msa$lf[msa$NAICS==54]   # Sci-Tech
msa$lf[msa$NAICS==62]   # Health Care


#*******************************************************************************************
#**    1.3) Growth rate analysis                                                          **
#*******************************************************************************************

## Wage and employment growth analysis 
msa$dlq         <- msa$lq08 - msa$lq98
msa$wage08      <- msa$AnnPay08 / msa$Empl08
msa$wage98      <- msa$AnnPay98 / msa$Empl98
msa$size        <- msa$Empl08 / msa$Est08

msa$w_growth    <- (((msa$wage08/msa$wage98)^(1/11)) - 1)*100
msa$e_growth    <- (((msa$Empl08/msa$Empl98)^(1/11)) - 1)*100
 

#*******************************************************************************************
#**    2) Visualising economic data                                                       **
#*******************************************************************************************

plot.msa        <- subset(msa,(msa$NAICS>0 & !msa$NAICS==95))   # Drop NAICS category 95 (Auxiliary)
plotclr         <- brewer.pal(11,"Spectral")   # setting the colour scheme from ColorBrewer (http://colorbrewer2.org/)


#*******************************************************************************************
#**    2.1) Bubble charts                                                                 **
#*******************************************************************************************

symbols(plot.msa$e_growth, plot.msa$w_growth, circles=plot.msa$wage08, xlab="Annual employment growth (%, 1998-2008)", 
    ylab="Annual nominal wage growth (%, 1998-2008)", fg="white", bg=plotclr, inches=0.75)
text(plot.msa$e_growth, plot.msa$w_growth, plot.msa$Des_short, cex=0.8)

symbols(plot.msa$wage08, plot.msa$w_growth, circles=plot.msa$size, xlab="Average annual wage ($000s, 2008)", 
    ylab="Annual nominal wage growth (%, 1998-2008)", fg="white", bg=plotclr, inches=0.75)
text(plot.msa$wage08, plot.msa$w_growth, plot.msa$Des_short, cex=0.8)

symbols(plot.msa$size, plot.msa$w_growth, circles=plot.msa$wage08, xlab="Average firm size (2008)", 
    ylab="Annual nominal wage growth (%, 1998-2008)", fg="white", bg=plotclr, inches=0.75)
text(plot.msa$size, plot.msa$w_growth, plot.msa$Des_short, cex=0.8)

## Bubble chart for HW5
[your code]

abline(v=1, col="gray", lty=2, lwd=2) # adds vertical line at x=1
abline(h=0, col="gray", lty=2, lwd=2) # adds horizontal line at y=0

[more code]

## Alternative bubble chart for HW5 using browser web interface
lq.web          <- gvisBubbleChart(data=plot.msa, idvar="Des_short", xvar="lq08", yvar="dlq", sizevar="Est08", options=list(height=800, width=1200))
plot(lq.web)


#*******************************************************************************************
#**    2.2) 3D plots of the shift-share analysis                                          **
#*******************************************************************************************

attach(msa)

# 3D scatter plot with labels
plot3d(ns[!NAICS==0],im[!NAICS==0],lf[!NAICS==0], col="red", xlab="National share", ylab="Industry mix", zlab="Local factors",
       size=3, type="s", radius=plot.msa$Empl08/10, top=T)

texts3d(ns[!NAICS==0], im[!NAICS==0], lf[!NAICS==0], Des, adj = 0.9, cex=0.9, col="darkred")

# Alternative version of 3D scatter plot with trend surface
scatter3d(ns[!NAICS==0],im[!NAICS==0],lf[!NAICS==0], xlab="National share", ylab="Industry mix", zlab="Local factors",
fogtype="exp2", grid="True", model.summary=F)

detach(msa)


#*******************************************************************************************
#**    2.3) Dynamic visualisation via Google API                                          **
#*******************************************************************************************
#  Frist data set reads in GDP data for all MSAs

gdp             <- read.csv("GDP_Metro_ALL.csv", h=T, na.strings="#N/A")

# Define the GDP deflator as a measure of inflation (Nominal GDP / Real GDP)
gdp$gdp_defl    <- (gdp$gdp_n/gdp$gdp_r)*100

# Define the LQ for manufacturing, FIPS code "998" indicates national data
gdp$lq_manuf    <- ((gdp$gdp_manuf/gdp$gdp_n))/(gdp$gdp_manuf[gdp$fips==998]/gdp$gdp_n[gdp$fips==998])

names(gdp)

# Plot empirical distribution of LQs for manufacturing
hist(gdp$lq_manuf, main="Distribution of LQs", sub="(All US cities at 2-digit NAICS)", cex.sub=0.8, breaks=25, 
     xlab="Relative specialisation of manufacturing (LQs)", col="lightgrey")


# Plot dynamic chart with changes in manufacturing specialisation
gdp.mi          <- gdp[gdp$state=="MI",]
MI              <- gvisMotionChart(gdp.mi, idvar="metropolitanarea", timevar="time", options=list(height=800, width=1200))
plot(MI)

# Plot spatial distribution of LQs for manufacturing across MSAs
MSA             <- gvisGeoChart(gdp, locationvar="metro_code", colorvar="lq_manuf", options=list(region="US", resolution="metros", height=800, width=1200))
plot(MSA)


# Second data set reads long panel data of MSA-level personal incomes, earnings by industry and industry-specific earning shares
ind.shares      <- read.csv("Earnings_MSA_MI.csv", header=T, na.strings="#N/A")
names(ind.shares)

Shares          <- gvisMotionChart(ind.shares, idvar="area", timevar="time", options=list(height=800, width=1200))
plot(Shares)



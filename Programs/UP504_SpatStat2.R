#*******************************************************************************************
#**                         VIII. APPLIED SPATIAL METHODS (Part 2)                        **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#**  Overview of spatial packages in R: http://cran.r-project.org/web/views/Spatial.html  **
#**                                                                                       **
#**                                                                                       **
#*******************************************************************************************



## Load spatial packages
library(maps)         ## Projections and basic maps
library(maptools)     ## Map manipulation
library(sp)           ## Spatial data manipulation
library(spdep)        ## Spatial statistics
library(gstat)        ## Geostatistics
library(splancs)      ## Kernel density modelling
library(spatstat)     ## Geostatistics
library(pgirmess)     ## Spatial autocorrelation
library(RColorBrewer) ## ColorBrewer implementation for R
library(classInt)     ## Class intervals
library(spgwr)        ## Geographically-weighted regression
library(RgoogleMaps)  ## GoogleMaps API interface
library(rgdal)        ## Bindings for the Geospatial Data Abstraction Library (GDAL)
library(fields)       ## Distance measurements

setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data/")



#*******************************************************************************************
#**   5)  Spatial flow data analysis                                                      **
#*******************************************************************************************

# http://flowingdata.com/2011/05/11/how-to-map-connections-with-great-circles/
# See also 

# Rae, A. (2011) "Flow-data Analysis with Geographical Information Systems: A Visual Approach"
#         Environment and Planning B, 2011, 38, 776-794



#*******************************************************************************************
#**   6)  Spatial DGPs and exploratory spatial data analysis                              **
#*******************************************************************************************

## From Contiguity neighbours (C) to weights matrices (W)
counties.MI         <- readShapePoly("GIS/county_MI.shp")
centroids.MI        <- coordinates(counties.MI)


tracts.DA2          <- readShapePoly("GIS/DA2F_Tracts.shp")
centroids.DA2       <- coordinates(tracts.DA2)
#water               <- tracts.DA2[(tracts.DA2$AWATER00/(tracts.DA2$ALAND00+tracts.DA2$AWATER00))>0.60,]
#tracts.DA2          <- tracts.DA2[!(tracts.DA2$AWATER00/(tracts.DA2$ALAND00+tracts.DA2$AWATER00))>0.60,]

C.DA2               <- poly2nb(tracts.DA2, queen=T)
C.DA2.rook          <- poly2nb(tracts.DA2, queen=F)

W.DA2               <- nb2listw(C.DA2, style="W", zero.policy=TRUE)
W.DA2.rook          <- nb2listw(C.DA2.rook, style="W", zero.policy=TRUE)
summary(W.DA2)

C.MI                <- poly2nb(counties.MI, queen=T)
C.MI.rook           <- poly2nb(counties.MI, queen=F)

W.MI                <- nb2listw(C.MI, style="W", zero.policy=TRUE)
W.MI.rook           <- nb2listw(C.MI.rook, style="W", zero.policy=TRUE)
summary(W.MI)

## Defining the Michigan State GeoRef projection (http://www.michigan.gov/documents/DNR_Map_Proj_and_MI_Georef_Info_20889_7.pdf)
proj.MI             <- CRS("+proj=omerc +lat_0=45.309166667 +lonc=-86.0 +alpha=337.255555556 
                       +k=0.9996 +x_0=499839.8337 +y_0=528600.2398 +ellps=GRS80 +datum=NAD83 +units=m")

proj.longlat        <- CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")   # used by USGS and NOOA


coords.MI           <- SpatialPoints(coordinates(counties.MI), proj4string=proj.MI)
coords.MI.new       <- spTransform(coords.MI, CRS("+proj=longlat"))

proj4string(counties.MI) <- proj.MI
counties.MI.new     <- spTransform(counties.MI, CRS("+proj=longlat"))

## Adding population weighted centroids
c.pop               <- cbind(counties.MI$POPLON, counties.MI$POPLAT)
colnames(c.pop)     <- c("LON", "LAT")
centroids.pop       <- SpatialPoints(coords=c.pop, proj4string=proj.longlat)


plot(counties.MI.new, border="lightslategrey", axes=T)
plot(coords.MI.new, add=T, pch=16, cex=0.5)
plot(centroids.pop, add=T, col="red", pch=16, cex=0.5)

## Plot connections for queen vs. rook contiguities
par(mar=rep(0,4))
plot(W.DA2,coords=coordinates(tracts.DA2), pch=21, cex=0.5, col="lightslategray")
plot(W.DA2.rook,coords=coordinates(tracts.DA2),pch=21, cex=0.5, col="red", add=T)

par(mar=rep(0,4))
plot(W.MI,coords=centroids.MI, pch=21, cex=0.75, col="lightslategray")
plot(W.MI.rook,coords=centroids.MI, pch=21, cex=0.75, col="red", add=T)

## Visualising the W matrix
z                   <- t(listw2mat(W.MI))
image(1:length(counties.MI), 1:length(counties.MI), z[,ncol(z):1], col=c("black", "orange", "orange", "orange", "yellow", "yellow", "yellow", "yellow", "yellow", "white", "white"), 
      main="B style", axes=F)

z                   <- t(listw2mat(W.DA2))
image(1:length(tracts.DA2), 1:length(tracts.DA2), z[,ncol(z):1], col=brewer.pal(9,"RdYlBu"), main="B style", axes=F)

image(1:length(tracts.DA2), 1:length(tracts.DA2), z[,ncol(z):1], col=c("black", "orange", "orange", "orange", "yellow", "yellow", "yellow", "yellow", "yellow", "white", "white"), 
      main="Spatial weights matrix", xlab="Observations (i)", ylab="Neighbours (j)")

dev.off()

## Testing for Spatial Dependence: Global and Local Autocorrelation (Moran's I)
lag.W.DA2           <- lag.listw(W.DA2,tracts.DA2$VAL, zero.policy=T)
cor.test(tracts.DA2$VAL, lag.W.DA2)

moran.plot(tracts.DA2$VAL, W.DA2, zero.policy=T)
obs                 <- moran.plot(tracts.DA2$VAL, W.DA2, zero.policy=T)
summary(obs)

LM                  <- localmoran(tracts.DA2$VAL, listw=W.DA2, zero.policy=T)
tracts.DA2$LM       <- abs(LM[,4])                 ## Extract z-scores for local Moran's I

LM.palette          <- colorRampPalette(c("white", "orange", "red", "red2", "red3"), space = "rgb")
spplot(tracts.DA2, zcol="LM", col.regions=LM.palette(40), main="Local Moran's I (|z| scores)", pretty=T)


#*******************************************************************************************
#**   7)  Spatial regression models                                                       **
#*******************************************************************************************

## OLS baseline model
mod.ols             <- lm(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2)
res.ols             <- mod.ols$residuals
summary(mod.ols)

classes_fx          <- classIntervals(res.ols, n=7, style="fixed", fixedBreaks=c(-400,-150,-30,10,75,250,750,1000), rtimes = 1)
pal                 <- brewer.pal(7,"RdYlGn")
cols                <- findColours(classes_fx, pal)

plot(tracts.DA2,col=cols, border="white", axes=T)
legend(x=-83.00, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from OLS Model", ncol=2)


## Residual Autocorrelation
moran.test(res.ols, listw=W.DA2, zero.policy=T)


## SAR Model (Careful: Model is computationally intensive)
mod.sar             <- lagsarlm(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2, 
                                listw=W.DA2, zero.policy=T, tol.solve=1e-18)
res.sar             <- mod.sar$residuals
summary(mod.sar)

classes_fx          <- classIntervals(res.sar, n=7, style="fixed", fixedBreaks=c(-400,-150,-30,10,75,250,750,1000), rtimes = 1)
pal                 <- brewer.pal(7,"RdYlGn")
cols                <- findColours(classes_fx, pal)

plot(tracts.DA2, col=cols, border="white", axes=T)
legend(x=-83.00, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SAR Model", ncol=2)

moran.test(res.sar, listw=W.DA2, zero.policy=T)

## SEM Model (Careful: Model is computationally intensive)
mod.sem             <- errorsarlm(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2, 
                                  listw=W.DA2, zero.policy=T, tol.solve=1e-18)
res.sem             <- mod.sem$residuals
summary(mod.sem)

classes_fx          <- classIntervals(res.sem, n=7, style="fixed", fixedBreaks=c(-400,-150,-30,10,75,250,750,1000), rtimes = 1)
pal                 <- brewer.pal(7,"RdYlGn")
cols                <- findColours(classes_fx, pal)

plot(tracts.DA2, col=cols, border="white", axes=T)
legend(x=-83.00, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SEM Model", ncol=2)

moran.test(res.sem, listw=W.DA2, zero.policy=T)


## SDM Model (Careful: Model is computationally intensive)
mod.sdm             <- lagsarlm(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2, type="mixed",
                                listw=W.DA2, zero.policy=T, tol.solve=1e-18)
res.sdm             <- mod.sdm$residuals
summary(mod.sdm)

classes_fx          <- classIntervals(res.sdm, n=7, style="fixed", fixedBreaks=c(-400,-150,-30,10,75,250,750,1000), rtimes = 1)
pal                 <- brewer.pal(7,"RdYlGn")
cols                <- findColours(classes_fx, pal)

plot(tracts.DA2, col=cols, border="white", axes=T)
legend(x=-83.00, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from SDM Model", ncol=2)

moran.test(res.sdm, listw=W.DA2, zero.policy=T)

## Comparing the residuals from all four models
RES                 <- cbind(res.ols, res.sar, res.sem, res.sdm)

summary(RES)
RES.trunc           <- RES[(res.sdm>-1.877e+01 & res.sdm<1.583e+01),]
boxplot(RES.trunc, col=c("grey","orange","green", "blue"))
abline(h=0, col="red", lwd=2)


## Geographically weighted regression, using baseline specification of HW2(c)
bwG                 <- gwr.sel(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2, gweight=gwr.Gauss, verbose=F)
mod.gwr             <- gwr(VAL ~ LI_UNITS + AGE + ROOMS + VACRATE + OOH + HH_MEDIAN + EDU + BLACK, data = tracts.DA2, bandwidth=bwG, gweight=gwr.Gauss)


res.gwr             <- mod.gwr$SDF$gwr.e
plot(res.gwr)

classes_fx          <- classIntervals(res.gwr, n=7, style="fixed", fixedBreaks=c(-300,-200,-100,0,100,250,500,1000), rtimes = 1)
#res.palette        <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
#pal                <- res.palette(7)
pal                 <- brewer.pal(7,"RdYlGn")
cols                <- findColours(classes_fx,pal)

plot(tracts.DA2,col=cols, border="white", axes=T)
legend(x=-83.00, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from GWR Model", ncol=2)

moran.test(res.gwr, listw=W.DA2, zero.policy=T)

coef.li_units       <- mod.gwr$SDF$LI_UNITS
coef.age            <- mod.gwr$SDF$AGE
coef.edu            <- mod.gwr$SDF$EDU
coef.black          <- mod.gwr$SDF$BLACK

plot(coef.li_units)
abline(h=0, col="red", lwd=2)

classes_qt          <- classIntervals(round(coef.li_units, digits=2), n=9, style="quantile", rtimes = 1)
cols                <- findColours(classes_qt, pal)

plot(tracts.DA2,col=cols, border="white", axes=T)
legend(x=-83.15, y=42.125, cex=0.75,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Local Coefficient Estimates (LIHCT Units)",ncol=2)

 
spplot(mod.gwr$SDF, "AGE", col.regions=pal, cuts=9, main="Local Coefficient Estimates (Age of House)")





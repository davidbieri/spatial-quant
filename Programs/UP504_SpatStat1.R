#*******************************************************************************************
#**                         VIII. APPLIED SPATIAL METHODS (Part 1)                        **
#**                                                                                       **
#**                     UP 504 -- Quantitative Methods W 2013                             **
#**                                  David Bieri                                          **
#**                                                                                       **
#**  Overview of spatial packages in R: http://cran.r-project.org/web/views/Spatial.html  **
#**                                                                                       **
#**                                                                                       **
#*******************************************************************************************



## Load spatial packages
library(foreign)
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

setwd("C:/Documents and Settings/cbbadbi/My Documents/Teaching/UP_504/Data")


#*******************************************************************************************
#**    1) Loading and Manipulating Spatial Data (using Points and Polygons)               **
#*******************************************************************************************

## SPATIAL POINTS
## Load LIHTC project locations (original data from http://lihtc.huduser.org/)
lihtc.data          <- read.dta("HUDUser/LIHTC_Projects_MI.dta")
summary(lihtc.data)

sp.vars             <- c("n_units", "li_units", "yr_pis", "yr_alloc", "type", "bond", "latitude", "longitud")
lihtc               <- lihtc.data[sp.vars]
lihtc               <- na.omit(lihtc)
length(lihtc[,1])

## Display LIHTC project locations
sp.points           <- cbind(lihtc$longitud, lihtc$latitude)
colnames(sp.points) <- c("LON", "LAT")
proj.utm            <- CRS("+proj=utm +zone=17 +datum=WGS84")                                   # Universal Transverse Mercator [60 zones, 6' band lat], World Geodesic Datum, used by DoD
proj.longlat        <- CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")   # used by USGS and NOOA
data.sp             <- SpatialPointsDataFrame(coords=sp.points, data=lihtc, proj4string=proj.utm)
data.sp2            <- SpatialPointsDataFrame(coords=sp.points, data=lihtc, proj4string=proj.longlat)

## Bounding box of data points
bbox(data.sp)
plot(data.sp, pch=1, axes=T, col="red")
plot(data.sp2, pch=1, axes=T, col="blue")


## SPATIAL POLYGONS
## Load spatial polygons with attributes as SpatialPolygonsDataframe
tracts.DA2          <- readShapePoly("GIS/DA2F_Tracts.shp")
centroids.DA2       <- coordinates(tracts.DA2)
summary(tracts.DA2)

plot(tracts.DA2, axes=T)

## Separate shapefiles for census tracts where  >60% surface area is water
#water               <- tracts.DA2[tracts.DA2$ALAND00==0,]
#tracts.DA2          <- tracts.DA2[!tracts.DA2$ALAND00==0,]

water               <- tracts.DA2[(tracts.DA2$AWATER00/(tracts.DA2$ALAND00+tracts.DA2$AWATER00))>0.60,]
tracts.DA2          <- tracts.DA2[!(tracts.DA2$AWATER00/(tracts.DA2$ALAND00+tracts.DA2$AWATER00))>0.60,]

plot(water, axes=T, xlim=bbox(tracts.DA2)[1,], ylim=bbox(tracts.DA2)[2,], col="lightblue")
plot(tracts.DA2, add=T, col="tan3")

## UTM projected data vs. long-lat plot 
par(mar=c(2,2,2,1)+0.1, mfrow=c(1,2))
plot(data.sp, pch=16, cex=.75, axes=T, col="blue")
plot(tracts.DA2, add=T)

plot(data.sp2, pch=16, cex=.75, axes=T, col="red")
plot(tracts.DA2, add=T)

par(mar=c(2,2,2,1)+0.1, mfrow=c(1,2))
plot(tracts.DA2,xlim=c(-83.5,-82.5),ylim=c(42,42.5),col="tan3", axes=T)    #City of Detroit
plot(data.sp2,pch=16, cex=.75,add=T, col="blue")

plot(tracts.DA2,xlim=c(-83.8,-83.4),ylim=c(42.1,42.4),col="tan3", axes=T)  #City of Ann Arbor and Ypsilanti
plot(data.sp2,pch=16, cex=.95,add=T, col="darkgreen")
dev.off()


#*******************************************************************************************
#**   1b)  Exporting maps to Google Earth                                                 **
#*******************************************************************************************

# Code modified from http://spatial-analyst.net/wiki/index.php?title=Export_maps_to_GE#Polygon_maps

proj4string(tracts.DA2) <- proj.longlat
tracts.GE   <- spTransform(tracts.DA2, CRS("+proj=longlat +datum=WGS84"))
writeOGR(tracts.GE["VACRATE"], "vacrate.kml", "VACRATE", "KML") 

varname     <- "VACRATE"  # variable name
maxvar      <- max(tracts.GE[varname]@data)  # maximum value
filename    <- file(paste(varname, "_bubble.kml", sep=""), "w",  blocking=FALSE)
write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", filename)
write("<kml xmlns=\"http://earth.google.com/kml/2.2\">", filename, append = TRUE)
write("<Document>", filename, append = TRUE)
write(paste("<name>", varname, "</name>", sep=" "), filename, append = TRUE)
write("<open>1</open>", filename, append = TRUE)
for (i in 1:length(tracts.GE@data[[1]])) {
   write(paste('  <Style id="','pnt', i,'">',sep=""), filename, append = TRUE)
   write("    <LabelStyle>", filename, append = TRUE)      
   write("     <scale>0.7</scale>", filename, append = TRUE)
   write("    </LabelStyle>", filename, append = TRUE)
   write("      <IconStyle>", filename, append = TRUE)
   write("  <color>ff0000ff</color>", filename, append = TRUE)
#  write("  <colorMode>random</colorMode>", filename, append = TRUE)
   write(paste("			<scale>", tracts.GE[i,varname]@data[[1]]/maxvar*2+0.3, "</scale>", sep=""), filename, append = TRUE)
   write("			<Icon>", filename, append = TRUE)
   write("				<href>http://maps.google.com/mapfiles/kml/shapes/realestate.png</href>", filename, append = TRUE)
   write("			</Icon>", filename, append = TRUE)
   write("		</IconStyle>", filename, append = TRUE)
   write("	</Style>", filename, append = TRUE)
 }
write("<Folder>", filename, append = TRUE)
write(paste("<name>Donut icon for", varname,"</name>"), filename, append = TRUE)
for (i in 1:length(tracts.GE@data[[1]])) {
   write("  <Placemark>", filename, append = TRUE)
   write(paste("  <name>", tracts.GE[i,varname]@data[[1]],"</name>", sep=""), filename, append = TRUE)
   write(paste("  <styleUrl>#pnt",i,"</styleUrl>", sep=""), filename, append = TRUE)
   write("    <Point>", filename, append = TRUE)
   write(paste("      <coordinates>",coordinates(tracts.GE)[[i,1]],",", coordinates(tracts.GE)[[i,2]],",10</coordinates>", sep=""), filename, append = TRUE)
   write("    </Point>", filename, append = TRUE)
   write("    </Placemark>", filename, append = TRUE)
 }
write("</Folder>", filename, append = TRUE)
write("</Document>", filename, append = TRUE)
write("</kml>", filename, append = TRUE)
close(filename)



#*******************************************************************************************
#**   2)  Mapping and geocoding spatial data (using GoogleMaps and OpenStreetMap)         **
#*******************************************************************************************

## Overlay maps from OSM (www.openstreetmap.org) and GoogleMaps

bb                  <- qbbox(bbox(sp.points)[1,],bbox(sp.points)[2,])

map.DA2.g           <- GetMap.bbox(bb$latR, bb$lonR, zoom=9, destfile = "DetroitGoogle.png", maptype="hybrid")
#map.DA2.g           <- GetMap.bbox(bb$latR, bb$lonR, zoom=10, destfile = "DetroitGoogle.png", maptype="hybrid")
PlotOnStaticMap(map.DA2.g, lon=data.sp$longitud, lat=data.sp$latitude, pch=21, col="red", bg="white", axes=T)

map.DA2.osm         <- GetMap.OSM(lonR=c(-83.947, -82.783), latR=c(42.084, 42.789), scale=200000, destfile = "DetroitOSM.png")
PlotOnStaticMap(map.DA2.osm, lon=data.sp$longitud, lat=data.sp$latitude, NEWMAP=T)

## More on mapping using new OSM library in R
#  http://www.r-bloggers.com/plot-maps-like-a-boss/?utm_source=feedburner&utm_medium=email&utm_campaign=Feed%3A+RBloggers+%28R+bloggers%29

library(plotGoogleMaps)

data(meuse)
coordinates(meuse)<-~x+y # convert to SPDF
proj4string(meuse) <- CRS('+init=epsg:28992')

m <-segmentGoogleMaps(meuse, zcol=c('zinc','dist.m'), mapTypeId='ROADMAP', filename='myMap4.htm', colPalette=c('#E41A1C','#377EB8'), strokeColor='black')


#*******************************************************************************************
#**   2a)  Analysing the spatial distribution of bus ridership in AA                      **
#*******************************************************************************************

## AATA Rideship Data (AATA Line 4, Washtenaw Ave, Sept.- Dec. 2011)
aata                <- read.csv("AATA/Ridership_F11/fall2011_aataridership.csv", header=T, na.strings="#N/A")
summary(aata)
aata                <- aata[ which(aata$Longitude < 0), ]

aata.A2             <- aata[aata$Direction=="Inbound", ]
aata.Yp             <- aata[aata$Direction=="Outbound", ]

## Display ridership locations
sp.points           <- as.data.frame(cbind(aata$Longitud, aata$Latitude, aata$Direction, aata$Time))
colnames(sp.points) <- c("Longitude", "Latitude", "Direction", "Time")

sp.points.A2        <- sp.points[sp.points$Direction==1, ]      # 1=Inbound, i.e. to Ann Arbor   
sp.points.Yp        <- sp.points[sp.points$Direction==2, ]

proj.utm            <- CRS("+proj=utm +zone=17 +datum=WGS84")                                   # World Geodesic Datum, used by DoD
proj.longlat        <- CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0")   # used by USGS and NOOA
data.sp             <- SpatialPointsDataFrame(coords=sp.points, data=aata, proj4string=proj.utm)
data.sp2            <- SpatialPointsDataFrame(coords=sp.points, data=aata, proj4string=proj.longlat)


## Bounding box of data points
bbox(data.sp)
plot(data.sp, pch=1, axes=T, col="red")


## Overlay maps from OSM (www.openstreetmap.org) and GoogleMaps

bb                  <- qbbox(bbox(data.sp)[1,],bbox(data.sp)[2,])

map.AATA.g          <- GetMap.bbox(bb$latR, bb$lonR, zoom=12, destfile = "A2Google.png", maptype="hybrid")
PlotOnStaticMap(map.AATA.g, lon=aata$Longitude, lat=aata$Latitude, col="red")


# Reading in bus stop data
AATALine4.A2        <- readShapePoints("AATA/stop_to_A2.shp")
stops.A2            <- as.data.frame(coordinates(AATALine4.A2))
colnames(stops.A2)  <- c("Longitude", "Latitude")
stops.A2$StopID     <- AATALine4.A2$StopID
stops.A2$StopName   <- AATALine4.A2$StopName

AATALine4.Yp        <- readShapePoints("AATA/stop_to_Ypsi.shp")
stops.Yp            <- as.data.frame(coordinates(AATALine4.Yp))
colnames(stops.Yp)  <- c("Longitude", "Latitude")
stops.Yp$StopID     <- AATALine4.Yp$StopID
stops.Yp$StopName   <- AATALine4.Yp$StopName


# Distance between bus stops and passenger entries
dist.A2          <- as.data.frame(rdist.earth(sp.points.A2[, c("Longitude", "Latitude")], stops.A2[, c("Longitude", "Latitude")], miles=T))
dist.Yp          <- as.data.frame(rdist.earth(sp.points.Yp[, c("Longitude", "Latitude")], stops.Yp[, c("Longitude", "Latitude")], miles=T))


# Assigning passenger records to nearest bus stop
d.A2             <- dist.A2
d.Yp             <- dist.Yp

cutoff           <- 0.025       # Distance cut-off per bus stop in miles, 50 yards ~ 0.025mi

d.A2[(d.A2 > cutoff)] = 0 
d.A2[(d.A2 <= cutoff) & (d.A2>0)] = 1
d.Yp[(d.Yp > cutoff)] = 0 
d.Yp[(d.Yp <= cutoff) & (d.Yp>0)] = 1

stops.A2$Passengers      <- colSums(d.A2)
stops.Yp$Passengers      <- colSums(d.Yp)

radius.A2        <- sqrt(stops.A2$Passengers/pi)
radius.Yp        <- sqrt(stops.Yp$Passengers/pi)

PlotOnStaticMap(map.AATA.g, lon=stops.A2$Longitude, lat=stops.A2$Latitude, col="green", lwd=3, pch=21, cex=radius.A2/15, add=T)
PlotOnStaticMap(map.AATA.g, lon=stops.Yp$Longitude, lat=stops.Yp$Latitude, col="yellow", lwd=3, pch=21, cex=radius.Yp/15, add=T)


#*******************************************************************************************
#**   2b)  Spatial point analysis with geocoding                                          **
#*******************************************************************************************

library("RJSONIO") # Serialize R objects to JSON, JavaScript Object Notation

getGeoCode <- function(gcStr)  {
  gcStr <- gsub(' ','%20',gcStr) #Encode URL Parameters
  #Open Connection
  connectStr <- paste('http://maps.google.com/maps/api/geocode/json?sensor=false&address=',gcStr, sep="") 
  con <- url(connectStr)
  data.json <- fromJSON(paste(readLines(con), collapse=""))
  close(con)
  #Flatten the received JSON
  data.json <- unlist(data.json)
  if(data.json["status"]=="OK")   {
    lat <- data.json["results.geometry.location.lat"]
    lng <- data.json["results.geometry.location.lng"]
    gcodes <- c(lng, lat)
    names(gcodes) <- c("Lon", "Lat")
    return (gcodes)
  }
}

reverseGeoCode <- function(latlng) {
  latlngStr <-  gsub(' ','%20', paste(latlng, collapse=","))#Collapse and Encode URL Parameters
  library("RJSONIO") #Load Library
  #Open Connection
  connectStr <- paste('http://maps.google.com/maps/api/geocode/json?sensor=false&latlng=',latlngStr, sep="")
  con <- url(connectStr)
  data.json <- fromJSON(paste(readLines(con), collapse=""))
  close(con)
  #Flatten the received JSON
  data.json <- unlist(data.json)
  if(data.json["status"]=="OK")
    address <- data.json["results.formatted_address"]
  return (address)
}


##     Geocoding of TCAUP staff/faculty addresses using GoogleMaps 

addresses           <- c("515 W. Madison Street, Ann Arbor MI", "2224 Applewood Court, Ann Arbor MI", "505 N. Division Street, Ann Arbor MI", 
                        "1010 Fairmount Drive, Ann Arbor MI", "2200 Fuller Court, Ann Arbor MI", "1689 Butterweed Court, Ann Arbor MI", 
                        "808 Duncan Street, Ann Arbor MI", "805 W. First Street, Ann Arbor MI", "223 East Ann Street, Ann Arbor MI", 
                        "1415 E. Stadium Boulevard, Ann Arbor MI", "456 Hilldale Drive, Ann Arbor MI", "1525 Westfield Avenue, Ann Arbor MI", 
                        "2118 Harding Avenue, Ypsilanti MI", "715 Spring Valley Road, Ann Arbor MI", "2939 Aurora Street, Ann Arbor MI", 
                        "6 Westbury Court, Ann Arbor MI", "632 W. Davis Avenue, Ann Arbor MI", "655 N. 5th Avenue, Ann Arbor MI")
                        

coords              <- lapply(addresses, function(val){getGeoCode(val)})
coords              <- matrix(unlist(coords), ncol=2, byrow=TRUE)
class(coords)       <- "numeric"
colnames(coords)    <- c("long","lat")

## Coordinates of TCAUP at 2000 Bonisteel Boulevard, A2
TCAUP               <- matrix(nrow=1, ncol=2)
colnames(TCAUP)     <- c("long", "lat")
TCAUP[1,1]          <- -83.71743
TCAUP[1,2]          <- 42.28898
tcaup               <- SpatialPoints(coords=coordinates(TCAUP))

## Calculate distances from TCAUP
dist.TCAUP          <- rdist.earth(TCAUP, coordinates(coords), miles=T)
attributes          <- as.data.frame(t(dist.TCAUP))
colnames(attributes)<- c("distance")
rownames(attributes)<- c("Bieri", "Allen", "Arquero de Alcaron", "Campbell", "Charles", "Deng", "Grengs", "Junghans",
                        "Kelbaugh", "Larsen", "Levine", "Norton", "Tietjen", "Ponce de Leon", "Shatkin", "Thomas", "Hoey", "Fishman")

data.urp            <- SpatialPointsDataFrame(coords, attributes, proj4string=proj.longlat)

bb                  <- qbbox(bbox(data.urp)[1,],bbox(data.urp)[2,])
map.A2.g            <- GetMap.bbox(bb$latR, bb$lonR, zoom=12, destfile = "A2Google.png", maptype="hybrid")

PlotOnStaticMap(map.A2.g, lon=coordinates(data.urp)[,1], lat=coordinates(data.urp)[,2], pch=21, col="blue", bg="gold", 
                cex=1.2, axes=T)
PlotOnStaticMap(map.A2.g, lon=coordinates(data.urp)[,1], lat=coordinates(data.urp)[,2], labels=rownames(data.urp@data), pos=3,
                font=2, col="white", cex=0.9, FUN=text, add=T)
PlotOnStaticMap(map.A2.g, lon=coordinates(tcaup)[,1], lat=coordinates(tcaup)[,2], col="gold", pch="M", font=2, cex=1.5, axes=T, add=T)


## Scaling names by inverse distances to TCAUP
PlotOnStaticMap(map.A2.g, lon=coordinates(tcaup)[,1], lat=coordinates(tcaup)[,2], col="gold", pch="M", font=2, cex=1.5, axes=T)
PlotOnStaticMap(map.A2.g, lon=coordinates(data.urp)[,1], lat=coordinates(data.urp)[,2], labels=rownames(data.urp@data), font=2, col="white", 
                cex=(1/(data.urp$distance)^(1/2)), FUN=text, add=T)


#*******************************************************************************************
#**   3)  Spatial point patterns: Identifying spatial neighbours                          **
#*******************************************************************************************

## k-nearest neighbours
k                   <- 4
nn                  <- knearneigh(coordinates(data.urp), k, longlat=T)
nn.urp              <- knn2nb(nn)

plot(data.urp, axes=T, pch=21, col="blue", bg="gold")
text(coordinates(data.urp), rownames(data.urp@data), font=2, col="black", pos=3, cex=0.9)
plot(tcaup, add=T, col="gold", pch="M", cex=1.5)
text(coordinates(tcaup), label="TCAUP", font=2, col="gold", pos=3, cex=0.9)

plot(nn.urp, coordinates(data.urp), col="green", lwd=2, add=T)


## Distance-based neighbors
d                   <- 2.5
dn.urp              <- dnearneigh(coordinates(data.urp), 0, d, longlat=TRUE)
plot(dn.urp, coordinates(data.urp),add=T, lwd=2, col="red")


## Spatial lags
lag.urp             <- nblag(nn.urp, 6)
plot(lag.urp[[2]], coordinates(data.urp), lty=2, add=T, lwd=2, col="orange")
plot(lag.urp[[6]], coordinates(data.urp), lty=2, add=T, lwd=2, col="blue")


#*******************************************************************************************
#**   4)  Visualising spatial patterns: Chloropeth maps and data-driven colour schemes    **
#*******************************************************************************************

#library(rgdal)      Alternative method for loading shape files using OGR
#tracts.DA2.alt      <- readOGR(dsn="[... path]/Teaching/UP_504/Data/GIS", layer="DA2F_Tracts")

plot(tracts.DA2, border="gray")
plot(tracts.DA2, col=tracts.DA2$COUNTYFP00, border="white")
plot(tracts.DA2, col=tracts.DA2$METRO, border="white")
plot(water, col="lightblue", add=T)
points(centroids.DA2, col="orange", pch=20)

#Defining colour scheme on the basis of county names
tracts.DA2$fac      <- factor(tracts.DA2@data$COUNTYFP00,label=c("Genesee", "Lapeer", "Livingstone", "Macomb", "Monroe", "Oakland", "St. Clair", "Washtenaw", "Wayne"))
tracts.DA2$fac.int  <- as.integer(tracts.DA2$fac)
nfac                <- max(as.integer(tracts.DA2$fac))
plotclr             <- brewer.pal(nfac,"Spectral")


## Chloropeth graphs using spplot
spplot(tracts.DA2, zcol="fac", col.regions=plotclr, scales=list(draw = TRUE), colorkey=TRUE)

trellis.par.set(sp.theme())
spplot(tracts.DA2, zcol="VAL", scales=list(draw = TRUE), colorkey=TRUE)
spplot(tracts.DA2, zcol="EDU", scales=list(draw = TRUE), colorkey=TRUE)
spplot(tracts.DA2, zcol="HH_MEDIAN", scales=list(draw = TRUE), colorkey=TRUE)

## Create different colour palettes
col.palette         <- colorRampPalette(c("red","green"), space = "rgb")
col.palette(n=5)

pal                 <- col.palette(n=5)
var                 <- tracts.DA2$VAL

classes_fx <- classIntervals(var, n=5, style="fixed", fixedBreaks=c(0, 50, 100, 250, 500, 1000), rtimes = 1)
classes_sd <- classIntervals(var, n=5, style = "sd", rtimes = 1)
classes_fi <- classIntervals(var, n=5, style = "fisher", rtimes = 3)
classes_eq <- classIntervals(var, n=5, style = "equal", rtimes = 1)
classes_km <- classIntervals(var, n=5, style = "kmeans", rtimes = 1)
classes_qt <- classIntervals(var, n=5, style = "quantile", rtimes = 1)

## Compare classes
par(mar=c(2,2,2,1)+0.1, mfrow=c(2,3))
plot(classes_fx, pal=pal, main="Fixed Intervals", xlab="", ylab="")
plot(classes_sd, pal=pal, main="Standard Deviation", xlab="", ylab="")
plot(classes_fi, pal=pal, main="Fisher-Jenks", xlab="", ylab="")
plot(classes_km, pal=pal, main="K Means", xlab="", ylab="")
plot(classes_eq, pal=pal, main="Equal Interval", xlab="", ylab="")
plot(classes_qt, pal=pal, main="Quantile", xlab="", ylab="")

dev.off()

cols_fx <- findColours(classes_fx, pal)
cols_sd <- findColours(classes_sd, pal)
cols_fi <- findColours(classes_fi, pal)
cols_km <- findColours(classes_km, pal)
cols_eq <- findColours(classes_eq, pal)
cols_qt <- findColours(classes_qt, pal)


par(mar=c(2,2,2,1)+0.1, mfrow=c(2,2))
plot(tracts.DA2, col=cols_fx, axes=F)
legend(x=-83.10, y=42.0,cex=.7,fill=attr(cols_fx,"palette"),bty="n",legend=names(attr(cols_fx, "table")),title="Median house value (2008, $000s)",ncol=2)

plot(tracts.DA2, col=cols_fi, axes=F)
legend(x=-83.10, y=42.0,cex=.7,fill=attr(cols_fi,"palette"),bty="n",legend=names(attr(cols_fi, "table")),title="Median house value (2008, $000s)",ncol=2)

plot(tracts.DA2, col=cols_eq, axes=F)
legend(x=-83.10, y=42.0,cex=.7,fill=attr(cols_eq,"palette"),bty="n",legend=names(attr(cols_eq, "table")),title="Median house value (2008, $000s)",ncol=2)

plot(tracts.DA2, col=cols_qt, axes=F)
legend(x=-83.10, y=42.0,cex=.7,fill=attr(cols_qt,"palette"),bty="n",legend=names(attr(cols_qt, "table")),title="Median house value (2008, $000s)",ncol=2)
dev.off()

## Useful code for finding colours that contain "grey" in name: colours()[grep("grey",colours())]



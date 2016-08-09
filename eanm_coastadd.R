rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("~/dropbox/testing todd/")

source("matching/eanm_functions.R")

library("sp")
library("maptools")
library("rgeos")
library("foreign")
library("stringi")
library("geosphere")
library("parallel")
data(wrld_simpl)

mdn = maptools::readShapeSpatial("matching/dat/mdn_withtrielev", proj4string=wgs)

coast110 = maptools::readShapeSpatial("~/downloads/data/ne_110m_coastline/ne_110m_coastline.shp", proj4string=wgs)
coast50 = maptools::readShapeSpatial("~/downloads/data/50m_physical/ne_50m_coastline.shp", proj4string=wgs)
coast10 = maptools::readShapeSpatial("~/downloads/data/ne_10m_coastline/ne_10m_coastline.shp", proj4string=wgs)

ethn_centres = rgeos::gCentroid(mdn, byid=TRUE)
plot(ethn_centres)

system.time(geosphere::dist2Line(ethn_centres[1], coast110))
system.time(geosphere::dist2Line(ethn_centres[1], coast50))
system.time(geosphere::dist2Line(ethn_centres[1], coast10))
# 10m: 15 sec -> 50h
# 50m: 4 sec -> 13h
# 110m: .4 sec -> ~ 1h

# parrallize
plot(c(8, 50, 100), c(34.938, 225.852, 425.139), type='b')
lines(c(8, 50, 100), c(18.620, 108.425, 200.247), type='b')

ethnlist = list()
for (i in 1:nrow(mdn)){
    ethnlist[[i]] = ethn_centres[i, ] 
}

ncores = parallel::detectCores() - 4
cl = parallel::makeCluster(ncores)
parallel::clusterExport(cl=cl, "coast50")
distances_50m <- parallel::parLapply(cl=cl, X=ethnlist, fun=function(x)geosphere::dist2Line(x, coast50))
stopCluster(cl)
write.csv(do.call(rbind, distances_50m), 'results_50m.csv')

cl = parallel::makeCluster(ncores)
parallel::clusterExport(cl=cl, "coast10")
system.time(distances_10m <- parallel::parLapply(cl=cl, X=ethnlist, fun=function(x)geosphere::dist2Line(x, coast10)))
stopCluster(cl)
write.csv(do.call(rbind, distances_10m), 'results_10m.csv')


curacao = lapply(ethnlist[1], geosphere::dist2Line, coast50)
plot(mdn[1, ])
plot(coast110, add=T, col='blue') # not plotted
plot(coast50, add=T, col='red')
plot(coast10, add=T, col='green')
points(ethn_centres[1, ])
points(curacao[[1]][,'lon'], curacao[[1]][,'lat'])

tahiti = lapply(ethnlist[c(2330, 2331, 2332)], geosphere::dist2Line, coast50)
plot(mdn[2330:2332, ])
plot(coast110, add=T, col='blue')
plot(coast50, add=T, col='red')
plot(coast10, add=T, col='green')
points(ethn_centres[2330:2332, ])
points(do.call(rbind, tahiti)[, c('lon', 'lat')])
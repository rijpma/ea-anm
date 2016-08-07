rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("~/dropbox/testing todd/")

source("matching/eanm_functions.R")

library("sp")
library("rgeos")
library("raster")
library("foreign")
library("stringi")
library("maptools")
data(wrld_simpl)

world_elev = raster::raster("~/downloads/gtopo_all.tif")

mdn = maptools::readShapeSpatial("matching/dat/murocknarodov.shp", proj4str=wgs)
names(mdn)

# ruggedness
tri = raster::terrain(world_elev, opt="TRI")
area = raster::area(world_elev)
areaXtri = area * tri
area = area * !is.na(tri) # na in tri -> area = 0
# patience!
wx = raster::extract(areaXtri, mdn, sum, na.rm=T)
w = raster::extract(area, mdn, sum)

mdn@data$tri = c(wx / w)

mdn@data = factor2char(mdn@data)
maptools::writeSpatialShape(mdn, 'matching/dat/mdn_withtri')
# mdn = maptools::readShapeSpatial("matching/dat/mdn_withtri.shp")

area = raster::area(world_elev)
areaXelev = area * world_elev
area = area * !is.na(world_elev)
wx = raster::extract(areaXelev, mdn, sum, na.rm=TRUE)
w = raster::extract(area, mdn, sum)
mdn@data$elev = c(wx / w)
maptools::writeSpatialShape(mdn, "matching/dat/mdn_withtrielev")

# quick inspection
pdf('matching/elev_by_eth.pdf', width=10, height=7)
plot(mdn, col=to_col(mdn$elev))
to_col(mdn$elev))
dev.off()

pdf('matching/figs/rgd_by_eth.pdf', width=10, height=7)
plot(mdn, col=to_col(mdn$tri), lwd=0.2)
add_legend(mdn$tri)
dev.off()

list.files(dirname(rasterTmpFile()))
removeTmpFiles(h=0)

aggregate(tri ~ v8, data=mdn, mean)
summary(lm(tri ~ v8, data=mdn))
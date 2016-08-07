rm(list=ls())

options(stringsAsFactors=FALSE)
setwd("~/dropbox/testing todd/")

source("matching/eanm_functions.R")

# gdistance w/o igraph dependencies, works fine
# install.packages("raster")
# install.packages("igraph")
# install.packages("NMF", dependencies=FALSE)
# remove.packages(c('dichromat', 'munsell', 'labeling', 'gtable', 'scales', 'ggplot2'))

library("sp")
library("rgeos")
library("raster")
library("gdistance")
library("maptools")
data(wrld_simpl)

raster::showTmpFiles()
raster::removeTmpFiles(h=0)

mdn = maptools::readShapeSpatial("matching/dat/mdn_withtrielev.shp")

# From Price and Bar-Yosef
middle_east = rect2pol(33.75, 32.6208701832, 44.296875, 38.376115424)
india = rect2pol(68.1113739014, 19.3111433551, 86.30859375, 31.6533813997)
china1 = rect2pol(106.0400390625, 34.5970415161, 111.8408203125, 39.1982053489)
china2 = rect2pol(110.0830078125, 29.4969875965, 123.3764648438, 33.02708758)
new_guinea = rect2pol(141.1083984375, -7.373362481, 145.0634765625, -4.2916356326)
africa = rect2pol(-17.05078125, 12.2970682929, 36.2109375, 20.6739052647)
usa = rect2pol(-88.0883789063, 35.9424357526, -82.529296875, 40.2292181887)
mexico = rect2pol(-101.2730164375, 16.5393294943, -97.2719054375, 19.4893735171)
guianas = rect2pol(-81.474609375, -5.3972734077, -55.1953125, 7.0354756524)
altiplano=rect2pol(-71.19140625, -16.3412256192, -65.56640625, -10.9627642564)
brazil_bolivia = Polygon(rbind(c(-65.6103515625, -9.167178733),
                c(-59.58984375, -18.958246486),
                c(-55.986328125, -15.982453523),
                c(-62.490234375, -7.4278365287)))
neo = pol2spdf(pol=list(middle_east, india, china1, china2, new_guinea, africa, usa,   mexico, guianas, altiplano, brazil_bolivia), 
    dat=data.frame(year=c(10e3,     4.5e3, 8e3,    8e3,    7e3,        3e3,    4.5e3, 7.3e3,  5.9e3,   6e3,       7e3),
                    region=c("middle east", "india", "inland china", "yangtze china", "new guinea", "africa", "usa", "mexico", 'guianas/amazon', 'peru/bolivia', 'brazil/bolivia'),
                    crop=c('wheat/barley', 'millet/pulses', 'millet', 'rice/foxnut', 'yam/banana/taro', 'rice/millet/sorghum', 'squash/sunflower', 'squash/maize/bean', 'yam/cotton/potatoe', 'potato/quinoa', 'peanut/manioc/chile')))
proj4string(neo) = wgs

world = raster::raster("~/downloads/gtopo_all.tif")
world_poly = rgeos::gUnaryUnion(wrld_simpl)
world = aggregate(world, 50)
world_nosea = raster::mask(world, world_poly)
world_nosea_nb = raster::focal(world_nosea, matrix(1, nrow=3, ncol=3), function(x) ifelse(is.na(x[-5]), x[5], x[-5]))
plot(world, xlim=c(-5, 5), ylim=c(50, 60))
plot(world_nosea, xlim=c(-5, 5), ylim=c(50, 60))
plot(world_nosea_nb, xlim=c(-5, 5), ylim=c(50, 60))

pc = coordinates(neo)
geo_dist = raster::pointDistance(pc, coordinates(mdn), longlat=TRUE)

world_nosea_nb = reclassify(world_nosea_nb, cbind(NA, 0))
world_nosea_nb[world_nosea_nb != 0] = 1

# still avoiding himalayas
tr = gdistance::transition(world_nosea_nb, mean, directions=8)
afghan_routes = gdistance::shortestPath(tr, coordinates(mdn[5, ]), pc[c(1:4, 6), ], output="SpatialLines")

png('matching/fig/distance_by_land_exmpl.png', width=720)
plot(world_nosea_nb, xlim=c(0, 130), ylim=c(0, 70), col=c(0, 'lightgray'))
points(pc)
points(coordinates(mdn[5, ]), col=2)
plot(afghan_routes, add=T, col=2)
dev.off()

tr_c_noscale = gdistance::geoCorrection(tr, "c") # cost
plot(raster(tr_c_noscale), xlim=c(-5, 5), ylim=c(50, 60))
least_cost_distance = gdistance::costDistance(tr_c_noscale, coordinates(mdn), pc)

dist2nearestneo = apply(least_cost_distance, 1, min)
neo_id = apply(least_cost_distance, 1, which.min)
neo_id[apply(least_cost_distance, 1, function(x) all(is.infinite(x)))] = NA

mdn@data$dist2neo_crow = apply(geo_dist, 2, min)
mdn@data[c('neo_crow_year_bp', 'neo_crow_region', 'neo_crow_crops')] = neo@data[apply(geo_dist, 2, which.min), ]
mdn@data$dist2neo_land = dist2nearestneo
mdn@data[c('neo_by_land_year_bp', 'neo_by_land_region', 'neo_by_land_crops')] = neo@data[neo_id, ]

maptools::writeSpatialShape(mdn, "matching/dat/mdn_withtrielevneo")

pdf('matching/neoregion.pdf', width=10)
plot(mdn, col=as.factor(mdn$neo_crow_region), lwd=0.2)
dev.off()

pdf('matching/neoregion_byland.pdf', width=10)
plot(mdn, col=as.factor(mdn$neo_by_land_region), lwd=0.2)
dev.off()

pdf('matching/dist2neo_land.pdf', width=10)
plot(mdn[!is.infinite(mdn$dist2neo_land), ], col=to_col(mdn$dist2neo_land[!is.infinite(mdn$dist2neo_land)]), lwd=0.2)
dev.off()

pdf('matching/dist2neo_crow.pdf', width=10)
plot(mdn, col=to_col(mdn$dist2neo_crow), lwd=0.2)
dev.off()
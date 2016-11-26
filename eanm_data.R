rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("~/dropbox/testing todd/")

source("matching/eanm_functions.R")

library("maptools")
library("foreign")
library("stringi")

# -------
# data in
# -------
greg2md = read.csv("matching/matchdat_mnl.csv")
names(greg2md)[5] = "greg_g1slname"

head(greg2md[, 1:10])

table(greg2md$greg_g1slname=='')

greg = maptools::readShapeSpatial('GREG/greg.shp')
greg@data = as.data.frame(sapply(greg@data, iconv, 'latin1', 'utf8'))


# fix encodings and transcription errors
# original encodings
greg@data[grep('Guar', greg@data$G1SHORTNAM), ]
greg@data[c(1, 150, 260, 262, 462), "G1SHORTNAM"]

enc_wrong = which(!stringi::stri_enc_isutf8(greg2md$greg_g1slname))

greg2md[enc_wrong, "greg_g1slname"]
greg2md[enc_wrong[1], "greg_g1slname"] = "Guarañoco [Guarañoco]"
greg2md[enc_wrong[2], "greg_g1slname"] = "Tehuelche and Óna [Tehuelche and Óna]"
greg2md["greg_g1sname"] = gsub(' \\[.*', '', greg2md$greg_g1slname)
greg2md$murdock_name[grep('NAIL', greg2md$murdock_name)] = "NAIL. . ."

ea = foreign::read.spss('old files/EthnographicAtlasWCRevisedByWorldCultures.sav')
ea = as.data.frame(ea)

# jutta's people list
asia <- read.csv('murdock_raw/asia-processed_if.csv')
africa <- read.csv('murdock_raw/africa-processed_if.csv')
europe <- read.csv("murdock_raw/Europe-processed_IF.csv")
latam <- read.csv("murdock_raw/LatinAmericaCarribean-processed_IF.csv")
northam <- read.csv("murdock_raw/North America-processed_IF.csv")
oceanea <- read.csv("murdock_raw/Oceanea-processed_IF.csv")
# socialism <- read.csv("murdock_raw/socialism-processed_IF.csv")
westasia <- read.csv("murdock_raw/WesternAsia-processed_IF.csv")
# nms <- c("ANM.name", "Society.Name")
murdock <- rbind(asia, africa, europe, latam, northam, oceanea, westasia)

# check greg and md names
greg2md$greg_g1sname[!greg2md$greg_g1sname %in% greg@data$G1SHORTNAM]
# only the missings!
greg2md$murdock_name[!tolower(gsub('[^A-z]', '', greg2md$murdock_name)) %in% 
                     tolower(gsub('[^A-z]', '', ea$v107))]
# only the missings and palestinians which aren't in the EA

dim(unique(greg2md[ , c('greg_g1sname', 'murdock_name')]))
dim(unique(greg2md[ , c('greg_g1sname', 'anm_name')]))
dim(unique(greg2md[ , c('greg_g1sname', 'murdock_name', 'anm_name')]))

greg2md_unq = unique(greg2md[ , c('greg_g1sname', 'murdock_name')])
ea$murdock_name = tolower(gsub('[^A-z]', '', ea$v107))
greg2md_unq$murdock_name = tolower(gsub('[^A-z]', '', greg2md_unq$murdock_name))

greg2md_unq[duplicated(greg2md_unq$murdock_name), ]
greg2md_unq[duplicated(greg2md_unq$greg_g1sname), ]

ea2map = merge(ea, greg2md_unq, by='murdock_name', all=TRUE)
ea2map[grep("NEWENGLAN", ea2map$v107), c("v107", "greg_g1sname")]
sum(duplicated(ea2map$greg_g1sname[ea2map$greg_g1sname!='' & !is.na(ea2map$greg_g1sname)]))
# watch these duplicates when working with the shapefile

# make suffienct number of extra polys for each
# give them new poly IDs
# and give both a numbered g1name

ea2map = ea2map[!is.na(ea2map$greg_g1sname) & ea2map$greg_g1sname!='', ]
# ea2map$greg_g1sname_unq[duplicated(ea2map$greg_g1sname)] = paste0(ea2map$greg_g1sname[duplicated(ea2map$greg_g1sname)], 1:sum(duplicated(ea2map$greg_g1sname)))
# ea2map$greg_g1sname_unq[!duplicated(ea2map$greg_g1sname)] = ea2map$greg_g1sname[!duplicated(ea2map$greg_g1sname)]
ethnicgrps = ea2map$greg_g1sname

head(sapply(slot(greg, "polygons"), function(x) slot(x, "ID")))
head(rownames(as(greg, "data.frame"))) # this one begins counting at 1
                                       # might need fixing   
greg_exp = greg
rownames(greg_exp@data) = as.numeric(rownames(greg_exp@data)) - 1
# sapply(slot(greg_exp, "polygons"), function(x) slot(x, "ID")) ==
# rownames(as(greg_exp, "data.frame")) # this one begins counting at 1
all.equal(sapply(slot(greg_exp, "polygons"), function(x) slot(x, "ID")), 
    rownames(as(greg_exp, "data.frame")))

grps_dpl = ethnicgrps[duplicated(ethnicgrps)]
for (i in 1:length(grps_dpl)){
    cat(i, '\n')
    group = grps_dpl[i]
    cat(group, '\n', ethnicgrps[duplicated(ethnicgrps)][i])

    tobind = greg[greg$G1SHORTNAM==group, ]
    tobind = spChFIDs(tobind, paste0(sapply(slot(tobind, "polygons"), function(x) slot(x, "ID")), group, i))
    rownames(tobind@data) = sapply(slot(tobind, "polygons"), function(x) slot(x, "ID"))

    grps_dpl[i] = paste0(group, i)
    tobind@data$G1SHORTNAM = paste0(group, i)
    greg_exp = maptools::spRbind(greg_exp, tobind)
}
head(rownames(greg_exp@data))
tail(rownames(greg_exp@data))

# place modified group names to match greg's
ea2map$greg_g1sname[duplicated(ea2map$greg_g1sname)] = grps_dpl

x = data.frame(greg_exp@data, ea2map[match(greg_exp$G1SHORTNAM, ea2map$greg_g1sname), ])
greg_exp@data = x
unique(greg_exp@data[grep("NEWENGLAN", greg_exp$v107), c("v107", "greg_g1sname")])

# what happens to the duplicated polygons in making neighbours
x = greg_exp[grep("Akan", greg_exp$greg_g1sname), ]
nb_exp = spdep::poly2nb(x)
nb_orig = spdep::poly2nb(x[grep("Akan$", x$greg_g1sname), ])
plot(nb_orig, jitter(coordinates(x[grep("Akan$", x$greg_g1sname), ]), factor=100))
plot(nb_exp, jitter(coordinates(x), factor=100))
# each is neighboured with itself
# which is ok
# but each self is neighboured with each self in each neighour
# have to average each poly maybe?

x = unique(cbind(coordinates(greg_exp), greg_exp$v104, greg_exp$v106))

# these multiple values are probably ok because  
# some people (eskimo, han chinese) live in many polys over wide area
plot(x[, 1], x[,2 ])
plot(x[, 4], x[,3 ])
plot(x[, 1], x[, 4])
plot(x[, 2], x[, 3])
max(as.numeric(x[, 3]), na.rm=T)
plot(greg_exp[greg_exp$v104=="  71" & !is.na(greg_exp$v104), ])
plot(greg_exp[greg_exp$v104=="  24" & !is.na(greg_exp$v104), ])
axis(1)
axis(2)


maptools::writeSpatialShape(greg_exp, 'matching/dat/murocknarodov')

table(!is.na(greg_exp$greg_g1sname))
sort(table(greg_exp$G1SHORTNAM))

png('matching/fig/maptest.png', width=480*3, height=480*2)
plot(greg_exp, col=greg_exp$v8, lwd=0.1)
legend('bottomleft', fill=unique(greg_exp$v8), legend=unique(greg_exp$v8))
dev.off()

pdf('matching/maptest.pdf', width=10, height=7)
plot(greg_exp, col=greg_exp$v8, lwd=0.1)
legend('bottomleft', fill=unique(greg_exp$v8), legend=unique(greg_exp$v8))
dev.off()
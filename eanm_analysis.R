rm(list=ls())
options(stringsAsFactors=FALSE, digits=3)
setwd("~/dropbox/testing todd")

source('matching/eanm_functions.r')

library("sp")
library("spdep")
library("maptools")
library("psych")
library("lmtest")
library("texreg")
library("countrycode")

oecdregions = read.csv("/Users/auke/Dropbox/cliodata/oecdregions.csv")

mdn = maptools::readShapeSpatial("matching/dat/mdn_withtrielevneo.shp")
dist2coast = read.csv('matching/dat/results_50m.csv')
mdn@data[c('dist2coast', 'lon_coast', 'lat_coast', 'coast_id')] = dist2coast[, -1]

names(mdn)[names(mdn)=='v102'] = 'year'

head(mdn@data$G1SHORTNAM)
length(unique(mdn@data$greg_g1sna)) / length(unique(mdn@data$G1SHORTNAM))

table(mdn$v74)
mdn$patrinherit = mdn$v74 %in% c("other patrilineal heirs", "patrilineal (sons)")
table(mdn$v7)
mdn$tofambride = mdn$v7 %in% c("bride price or wealth, to bride's family", "bride Service, to bride's family")
table(mdn$v12)
mdn$marresman = mdn$v12 %in% c("Patrilocal", "Virilocal")
table(mdn$v43)
mdn$patridescent = mdn$v43 %in% "patrilineal"
table(mdn$v15)
mdn$clans = mdn$v15 %in% c("Clan communities, or clan barrios", "Segmented communities")
table(mdn$v25)
mdn$cousinmar = mdn$v25 %in% c("Duolateral, matrilateral preference", 
    "Duolateral, patrilineal preference", 
    "Duolateral, symmetrical preference", 
    "Duolateral, with maternal cousins only, MoBrDa", 
    "Matrilineal cross-cousin, MoBrDa only", 
    "Patrilineal cross-cousin, FaSiDa only", 
    "Quadrilateral, FaSiDa preferred", 
    "Quadrilateral, matrilineal preference", 
    "Quadrilateral, symmetrical preference", 
    "Trilateral with bilateral preference")
table(mdn$v8)
mdn$polygamous = mdn$v8 %in% c("Independent polyandrous families",
    "Polygynous: unusual co-wives",
    "Polygynous: usual co-wives")
table(mdn$v8)
mdn$extended = mdn$v8 %in% c("Large extended families", "Small extended families")

mdn$patrinherit[mdn$v74=="missing data" | is.na(mdn$v74)] = NA
mdn$tofambride[mdn$v7=="Missing data" | is.na(mdn$v7)] = NA
mdn$marresman[mdn$v12=="missing data" | is.na(mdn$v12)] = NA
mdn$patridescent[mdn$v43=="missing data" | is.na(mdn$v43)] = NA
mdn$clans[mdn$v15=="Missing data" | is.na(mdn$v15)] = NA
mdn$cousinmar[mdn$v25=="missing data" | is.na(mdn$v25)] = NA
mdn$polygamous[mdn$v8=="Missing data" | is.na(mdn$v8)] = NA
mdn$extended[mdn$v8=="Missing data" | is.na(mdn$v8)] = NA
vrbs <- c('patrinherit',
          'tofambride',
          # 'd_bothfamilies',  # no burden, no trade
          "marresman",
          'patridescent',
          'clans',
          # 'd_patrilineages',   # 'd_matrilineages', # 'd_bilatlineages',  # no hierarchy
          # 'd_patrikinterms', # no hierarchy
          'cousinmar',
          'polygamous',
          'extended'
          )
fullvrbs <- c(patrinherit='Patr. inheritance',
              tofambride='Marr. payment: to family bride',
              marresman='Residence near family groom',
              patridescent='Patr. descent',
              clans='Clans present',
              # d_patrilineages='Patri. lineages',
              cousinmar='Consang. marriage',
              polygamous='Polygamous',
              extended='Extended household'
          # 'd_bothfamilies',
          # 'd_patrikinterms',
          )

mat = mdn@data[, vrbs]
mat[, ] = lapply(mat, as.numeric)
mat = mat[complete.cases(mat), ]
psych::alpha(mat)

fan = factanal(mat, factors=1, scores="regression")
fap = psych::fa(mat, poly=TRUE, cor="poly")
prc = prcomp(mat)


mdn@data$fa = fan$scores[match(rownames(mdn@data), rownames(fan$scores))]
mdn@data$fapoly = fap$scores[match(rownames(mdn@data), rownames(fan$scores))]
mdn@data$pc1 = prc$x[match(rownames(mdn@data), rownames(fan$scores)), "PC1"]
mdn@data$pc1 = max(mdn@data$pc1, na.rm=T) - mdn@data$pc1 / max(mdn@data$pc1, na.rm=T)

mdn@data = factor2char(mdn@data)

mdn_a = aggregate(mdn[, !sapply(mdn@data, is.character)], by=list(mdn$FeatureID), FUN=function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm=T)))

# no work atm
# twowaycol = RColorBrewer::brewer.pal(3, 'Set3')[1:2]
# for (vrb in vrbs){
#     png(paste0('matching/fig/fvrbs_by_eth', vrb, '.png'), 
#         width=1080, height=600)
#     plot(mdn_a, col=to_col(mdn_a@data[, vrb]), lwd=0.2)
#     add_legend(mdn_a@data[, vrb])
#     title(main=fullvrbs[vrb])
#     dev.off()
# }

# x = unique(data.frame(mdn_a@data[, vrb], to_col(mdn_a@data[, vrb])))
# x = x[order(x[, 1]), ]
# plot(x[, 1], 1:13, col=x[,2])

# commented out for speed
# png('matching/fig/fcowa_fa.png', width=1080, height=600)
# plot(mdn_a, col=to_col(mdn_a@data$fa), lwd=0.2, main='Fam. constraints on women’s agency: FA')
# add_legend(x=mdn_a@data$fa)
# dev.off()

# png('matching/fig/fcowa_fapolych.png', width=1080, height=600)
# plot(mdn_a, col=to_col(mdn_a@data$fapoly), lwd=0.2, main='Fam. constraints on women’s agency: polych. FA')
# add_legend(mdn_a@data$fapoly)
# dev.off()

# png('matching/fig/fcowa_prcomp.png', width=1080, height=600)
# plot(mdn_a, col=to_col(mdn_a@data$pc1), lwd=0.2, main='Fam. constraints on women’s agency: 1st principal component')
# add_legend(mdn_a@data$pc1)
# dev.off()

# png('matching/fig/dist2neo_crow.png', width=1080, height=600)
# plot(mdn, col=to_col(mdn@data$dist2neo_c), lwd=0.2, main='Distance to nearest neolitihic revolution site, “as crow flies”')
# add_legend(mdn@data$dist2neo_c)
# dev.off()

# png('matching/fig/dist2neo_land.png', width=1080, height=600)
# plot(mdn[!is.infinite(mdn@data$dist2neo_l), ], col=to_col(mdn@data$dist2neo_l[!is.infinite(mdn@data$dist2neo_l)]), lwd=0.2, main='Distance to nearest neolitihic revolution site, “as wolf runs”')
# add_legend(mdn@data$dist2neo_l[!is.infinite(mdn@data$dist2neo_l)])
# dev.off()

# png('matching/fig/neoregion_land.png', width=1080, height=600)
# plot(mdn, col=as.factor(mdn@data$neo_by_lan.1), lwd=0.2, main='Neolothic region (by land)')
# dev.off()

# png('matching/fig/dist2coast.png', width=1080, height=600)
# plot(mdn, col=to_col(mdn@data$dist2coast), lwd=0.2, main='Distance to coast, “as crow flies” from poly. centre')
# add_legend(mdn@data$dist2coast)
# dev.off()

# png('matching/fig/elevation.png', width=1080, height=600)
# plot(mdn, col=to_col(mdn@data$elev), lwd=0.2, main='Mean elevation (m)')
# add_legend(mdn@data$elev)
# dev.off()

# png('matching/fig/ruggedness.png', width=1080, height=600)
# plot(mdn, col=to_col(mdn@data$elev), lwd=0.2, main='Mean ruggedness (mean of the absolute differences between cell its 8 neighbours')
# add_legend(mdn@data$elev)
# dev.off()

mdn@data$dist2coast_trunc = ifelse(mdn@data$dist2coast < 0, 0, mdn@data$dist2coast)
mdn@data$iso3c = countrycode::countrycode(mdn@data$COW, 'cown', 'iso3c')
mdn@data$iso3n = countrycode::countrycode(mdn@data$COW, 'cown', 'iso3n')
mdn@data$region = countrycode::countrycode(mdn@data$COW, 'cown', 'region')
rownames(oecdregions) = oecdregions$ccode
mdn@data$oecdregion = oecdregions["region"][as.character(mdn@data$iso3n), 'region']
# or v91?

m_neo_l = lm(fapoly ~ log1p(dist2neo_l) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_neo_lt = lm(fapoly ~ log1p(dist2neo_l / neo_by_lan) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_neo_c = lm(fapoly ~ log1p(dist2neo_c) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_ele = lm(fapoly ~ log1p(elev) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_cst = lm(fapoly ~ log1p(dist2coast_trunc) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_rug = lm(fapoly ~ log1p(tri) + year - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])

range(exp(m_ele$model[, 2]))
range(exp(m_cst$model[, 2]))

m_neo_l_wrg = lm(fapoly ~ log1p(dist2neo_l) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_neo_lt_wrg = lm(fapoly ~ log1p(dist2neo_l / neo_by_lan) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_neo_c_wrg = lm(fapoly ~ log1p(dist2neo_c) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_ele_wrg = lm(fapoly ~ log1p(elev) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_cst_wrg = lm(fapoly ~ log1p(dist2coast_trunc) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])
m_rug_wrg = lm(fapoly ~ log1p(tri) + year + oecdregion - G1SHORTNAM - FeatureID, data=mdn@data[!is.infinite(mdn@data$dist2neo_l), ])

reglist = list(m_neo_l, m_neo_lt, m_neo_c, m_ele, m_cst, m_rug)
reglist_wrg = list(m_neo_l_wrg, m_neo_lt_wrg, m_neo_c_wrg, m_ele_wrg, m_cst_wrg, m_rug_wrg)

cfs = lapply(reglist, function(m) lmtest::coeftest(m, vcov=vcovCL(m, 'G1SHORTNAM')))
selist = lapply(cfs, function(x) x[, 2])
pvlist = lapply(cfs, function(x) x[, 4])
htmlreg(reglist, override.se=selist, override.pval=pvlist,
    digits=3, stars=c(0.1, 0.05, 0.01),
    file='matching/tab/regs.html')

cfs = lapply(reglist_wrg, function(m) lmtest::coeftest(m, vcov=vcovCL(m, 'G1SHORTNAM')))
selist = lapply(cfs, function(x) x[, 2])
pvlist = lapply(cfs, function(x) x[, 4])
htmlreg(reglist_wrg, override.se=selist, override.pval=pvlist,
    digits=3, stars=c(0.1, 0.05, 0.01),
    file='matching/tab/regs_wrgn.html')
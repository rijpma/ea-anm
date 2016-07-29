rect2pol = function(xleft, ybottom, xright, ytop){
    out = Polygon(matrix(c(xleft, xleft, xright, xright, xleft,
               ybottom, ytop, ytop, ybottom, ybottom), ncol=2))
    return(out)
}
pol2spdf = function(pol, dat){
    pols = list()
    for (i in 1:length(pol)){
        pols[[i]] = Polygons(list(pol[[i]]), i)
    }
    spols = SpatialPolygons(pols, )
    spdf = SpatialPolygonsDataFrame(spols, data=dat, match.ID=F)
    return(spdf)
}
to_col = function(x, cuts=9, pal='RdPu'){
    plt <- RColorBrewer::brewer.pal(cuts, pal)
    col <- plt[cut(x, 9)]
    return(col)
}
factor2char <- function(dat){
    factors <- sapply(dat, class) == 'factor'
    dat[factors] <- sapply(dat[factors], as.character)
    return(dat)
}
lat2km <- function(lat){
    110.574 / lat
}
lon2km <- function(lon, lat){
    (111.320 * cos(lat * (pi/180))) / lon
}
vcovCL <- function(fit, cluster1name, cluster2name=NULL){
    # there is also mfx:::clusterVCV but it requires a dataset + formula 
    # rather than a fitted model
    library(sandwich)
    cluster1 <- fit$model[, cluster1name]
    cluster2 <- fit$model[, cluster2name]
    cluster12 <- paste0(fit$model[, cluster1name], fit$model[, cluster2name])

    N <- length(cluster1)
    K <- fit$rank
    
    M1 <- length(unique(cluster1))
    dfc1 <- (M1 / (M1 - 1)) * ((N - 1) / (N - K))
    u1j <- apply(estfun(fit), 2, function(x) tapply(x, cluster1, sum))
    vc1 <- dfc1 * sandwich(fit, meat=crossprod(u1j) / N)
    out <- vc1
    if (!is.null(cluster2name)){
        M2 <- length(unique(cluster2))
        M12 <- length(unique(cluster12))
        dfc2 <- (M2 / (M2 - 1)) * ((N - 1) / (N - K))
        dfc12 <- (M12 / (M12 - 1)) * ((N - 1) / (N - K))
        u2j <- apply(estfun(fit), 2, function(x) tapply(x, cluster2, sum))
        u12j <- apply(estfun(fit), 2, function(x) tapply(x, cluster12, sum))
        vc2 <- dfc2 * sandwich(fit, meat=crossprod(u2j) / N)
        vc12 <- dfc12 * sandwich(fit, meat=crossprod(u12j) / N)
        out <- vc1 + vc2 - vc12
    }
    return(out)
}
factor2char <- function(dat){
    factors <- sapply(dat, class) == 'factor'
    dat[factors] <- sapply(dat[factors], as.character)
    return(dat)
}
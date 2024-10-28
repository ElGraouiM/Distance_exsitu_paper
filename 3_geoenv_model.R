#### 3 - fitting models to compute geo-environmental distance
## packages used: terra

wd <- "XS_vigna"

setwd(wd)
# get species data
sv <- readRDS("vigna_occurrences.rds")
sv <- terra::vect(sv, c("lon", "lat"), crs="+proj=longlat")

# get environmental data
env <- terra::rast("data/intermediate/wc.tif")[[c("bio_1", "bio_12")]]
names(env) <- c("tmp", "prc")

##### create envdist function 
cells <- terra::cellFromXY(terra::rast(), terra::crds(sv))
smp <- sv[!duplicated(cells), ]
gd <- round(terra::distance(smp, unit="km"))
x <- terra::extract(env, smp, ID=F)
ed1 <- dist(x[,1])
ed2 <- dist(x[,2])
g <- data.frame(g=gd, tmp=round(ed1, 1), prc=round(ed2))

a1 <- aggregate(g~tmp, data=g, median)
a2 <- aggregate(g~prc, data=g, median)

write.csv(a1, "data/intermediate/geo_temp.csv", row.names=FALSE)
write.csv(a2, "data/intermediate/geo_prec.csv", row.names=FALSE)

b1 <- aggregate(tmp~g, data=g, median)
b2 <- aggregate(prc~g, data=g, median)

# fit regression model
mtmp <- loess(g~tmp, data=b1)
p <- predict(mtmp)
e <- b1$g - p
itmp <- e < 250
mtmp <- loess(g~tmp, data=b1[itmp,], span=0.25)

mprc <-  loess(g~prc, data=b2)
p <- predict(mprc)
e <-  b1$g - p
iprc <-  e < 250
mprc <-  loess(g~prc, data=b2[iprc,], span=.25)

# export models
saveRDS(mprc, "data/intermediate/m_prc.rds")
saveRDS(mtmp, "data/intermediate/m_tmp.rds")

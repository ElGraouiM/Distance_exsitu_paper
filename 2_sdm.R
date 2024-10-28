## 2- species distribution modeling and range definition
## packages used: terra, predicts, XSitu

wd <- "XC_vigna"

setwd(wd)
wc <- terra::rast("data/intermediate/wc.tif")
afr <- terra::readRDS("data/intermediate/africa.rds")
spp <- XSdata::vigna_occ()


options(java.parameters = "-Xmx8g")

#run model
set.seed(123)
bg <- predicts::backgroundSample(wc, 10000)

usp <- sort(unique(spp$species))
eva <- list()

for (sp in usp) {
  cat(sp, "\n"); flush.console()
  obs <- unique(spp[spp$species==sp, c("lon", "lat")])
  obs <- terra::vect(obs, crs="+proj=longlat")
  
  nfld <- min(5, nrow(obs))
  k <- predicts::folds(obs, nfld)
  keva <- list()
  for (i in 1:nfld) {
    if (nfld == 1) {
      m <- predicts::MaxEnt(wc, obs)
      p <- obs		
    } else {
      m <- predicts::MaxEnt(wc, obs[k!=i,])
      p <- obs[k==i,]
    }
    nbg <- max(100, min(nrow(bg), nrow(p) * 5))
    b <- bg[sample(nbg), ]
    pa <- suppressWarnings(predicts::pa_evaluate(p, b, m, wc))
    keva[[i]] <- cbind(pa@stats, pa@thresholds)
  }
  eva[[sp]] <- do.call(rbind, keva) |> colMeans()
  
  f1 <- paste0("data/intermediate/sdm/raw/raw_", sp, ".tif")
  #	if (file.exists(f1)) next
  
  m <- predicts::MaxEnt(wc, obs)
  p <- terra::predict(wc, m, filename=f1, overwrite=TRUE)
  
  f2 <- paste0("data/intermediate/sdm/range/range_", sp, ".tif")
  th <- eva[[sp]]["equal_sens_spec"]
  p_range <- p >= th
  p_range <- terra::subst(p_range, FALSE, NA, filename=f2, overwrite=TRUE)
  
  f3 <- paste0("data/intermediate/sdm/adj/adj_", sp, ".tif")
  
  p_adj <- XSitu::adjust_range(p_range, obs, afr, CAmin=100000, CAmax=250000)
  
  p_adj <- terra::writeRaster(p_adj, f3, overwrite=TRUE)
  
}

#write evaluation csv
eva <- do.call(rbind, eva)
write.csv(eva, "data/intermediate/sdm/evaluation.csv")


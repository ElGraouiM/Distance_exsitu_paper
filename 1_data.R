## 1- data preparation
## packages used: terra, geodata

wd <- "XC_vigna"

dir.create(wd, FALSE, TRUE)
setwd(wd)

#prepare directories
dir.create("data/intermediate/sdm/raw", FALSE, TRUE)
dir.create("data/intermediate/sdm/range", FALSE)
dir.create("data/intermediate/sdm/adj", FALSE)

#prepare data and export to directories
##download vigna_occurrences.rds and put it in your "vigna" folder

afext <- terra::ext(-20, 54, -36, 40)

w <- geodata::world(path = "data/raw")
cc <- geodata::country_codes()
cc <- cc[cc$continent == "Africa", ]
w <- w[w$GID_0 %in% cc$ISO3]
afext <- c(-20, 54, -36, 40)
w <- terra::crop(w, afext)
terra::saveRDS(w, "data/intermediate/africa.rds")

wc <- geodata::worldclim_global("bio", res=10, path="data/raw", version="2.1")
wc <- wc[[c(1,8,9,12,15,16,17)]]
names(wc) <- gsub("wc2.1_10m_", "", names(wc))
wc <- terra::crop(wc, afext)
wc <- terra::mask(wc, w, filename="data/intermediate/wc.tif", overwrite=TRUE)


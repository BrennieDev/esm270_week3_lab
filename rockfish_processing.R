library(here)
library(tidyverse)
library(raster)

# Set our probability threshold for a species to be considered "present"
presence_threshold = 0.75

# Get list of all the sebastes files
files <- list.files(path=here("sebastes"), pattern="*.tif", full.names=TRUE, recursive=FALSE)

# Read in cumulative threats raster
threats_all <- raster(here("full_modelnv.tif"))

#raster_stack <- stack(files)
#presence_stack <- raster_stack %>% 
# reclassify(rcl = c(-Inf, 0.75, 0, 0.75,1,1))
#sebastes_richness <- raster::calc(presence_stack, fun = sum)

# Read in first sebastes raster
base_raster <- raster(files[1]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))

# Read in each of the other sebastes rasters
for (i in c(2:length(files))) {
  # Reclassify as present (1) or not present (0) based on 0.75 probability threshold
  temp <- raster(files[i]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  
  # Add this raster to the cumulative sum of other sebastes presence rasters
  base_raster <- sum(base_raster, temp)
}

# base_raster now contains a raster with the number of rockfish species present in each cell

#### Function to output a binary raster based on a user-given quantile (default is top 20%) ###
reclassify_topx <- function(rast,quant=0.8) {
  topx <- quantile(rast,quant) #find the 80% quantile of the raster values
  maxVal <- cellStats(rast,max) #find the maximum
  rcl <- c(-Inf,topx,0,
           topx,maxVal,1) # reclassify matrix (see help file for ?reclassify)
  out <- reclassify(rast,rcl=rcl)
  return(out) # returns the new binary raster
}

# Resample sebastes richness to higher resolution
base_raster <- resample(base_raster,threats_all,method='ngb',progress='text')

# Get top 20% of rockfish richness values
rockfish_top_20 <- reclassify_topx(base_raster, quant=0.8)

# Get top 20% of impact values
impact_top_20 <- reclassify_topx(threats_all, quant=0.8)

# Overlay rasters to get regions that are in both top 20% of impact and top 20% of rockfish diversity
hotspots <- overlay(rockfish_top_20,impact_top_20,fun=function(x,y){x*y})

# Set any values of 0 to NA
hotspots <- reclassify(hotspots, c(-Inf, 0, NA), right=TRUE)

# write files to view in gis
writeRaster(rockfish_top_20, "sebastes_top20.tif", overwrite = TRUE)
writeRaster(impact_top_20, "threats_top20.tif", overwrite = TRUE)
writeRaster(hotspots, "hotspots_cumulative.tif", overwrite = TRUE)


# Perform analysis for dem destr fishing and recreational fishing threats
threats_ddf <- raster(here("impact_dem_d.tif"))
impact_ddf_top_20 <- reclassify_topx(threats_ddf, quant=0.8)
hotspots_ddf <- overlay(rockfish_top_20,impact_ddf_top_20,fun=function(x,y){x*y})
hotspots_ddf <- reclassify(hotspots_ddf, c(-Inf, 0, NA), right=TRUE)
writeRaster(hotspots_ddf, "hotspots_dem_des_fish.tif", overwrite = TRUE)

threats_rf <- raster(here("impact_rec_fish.tif"))
impact_rf_top_20 <- reclassify_topx(threats_rf, quant=0.8)
hotspots_rf <- overlay(rockfish_top_20,impact_rf_top_20,fun=function(x,y){x*y})
hotspots_rf <- reclassify(hotspots_rf, c(-Inf, 0, NA), right=TRUE)
writeRaster(hotspots_rf, "hotspots_rec_fish.tif", overwrite = TRUE)






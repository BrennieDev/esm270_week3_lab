library(here)
library(tidyverse)
library(raster)


## Function to output a binary raster based on a user-given quantile (default is top 20%) ###
reclassify_topx <- function(rast,quant=0.8) {
  topx <- quantile(rast,quant) #find the 80% quantile of the raster values
  maxVal <- cellStats(rast,max) #find the maximum
  rcl <- c(-Inf,topx,0,
           topx,maxVal,1) # reclassify matrix (see help file for ?reclassify)
  out <- reclassify(rast,rcl=rcl)
  return(out) # returns the new binary raster
}

## Function to output the intersected hotspots raster
calculate_hotspots <- function(richness_raster, impact_raster, save_with = NULL) {
  
  # Get top 20% of richness values
  richness_raster[richness_raster <= 0] <- NA
  richness_top_20 <- reclassify_topx(richness_raster, quant=0.8)
  
  # Get top 20% of impact values
  impact_raster[impact_raster <= 0] <- NA
  impact_top_20 <- reclassify_topx(impact_raster, quant=0.8)
  
  # Resample richness to higher resolution
  richness_top_20 <- resample(richness_top_20,impact_raster,method='ngb',progress='text')
  
  # Overlay rasters to get regions that are in both top 20% of impact and top 20% richness
  hotspots <- overlay(richness_top_20,impact_top_20,fun=function(x,y){x*y})
  
  # Set any values of 0 to NA
  hotspots <- reclassify(hotspots, c(-Inf, 0, NA), right=TRUE)
  
  if (!is.null(save_with)) {
    writeRaster(richness_top_20, paste(save_with,"_richness_20.tif"), overwrite = TRUE)
    writeRaster(impact_top_20, paste(save_with,"_impact_20.tif"), overwrite = TRUE)
  }
  return(hotspots)
}


### Cumulative impact + total richness Analysis ###

# Read in cumulative threats raster
threats_all <- raster(here("full_modelnv.tif"))

# Read in total species richness data
total_richness <- raster(here("ca_curr_sp_rich.tif"))

# Find hotspots
cumulative_hotspots <- calculate_hotspots(total_richness, threats_all,"total")
writeRaster(cumulative_hotspots, "hotspots_cumulative.tif", overwrite = TRUE)


### Sebastes Analysis ###

# Set our probability threshold for a species to be considered "present"
presence_threshold = 0.75

# Get list of all the sebastes files
files <- list.files(path=here("sebastes"), pattern="*.tif", full.names=TRUE, recursive=FALSE)

# Read in first sebastes raster
rockfish_raster <- raster(files[1]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))

# Read in each of the other sebastes rasters
for (i in c(2:length(files))) {
  # Reclassify as present (1) or not present (0) based on 0.75 probability threshold
  temp <- raster(files[i]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  
  # Add this raster to the cumulative sum of other sebastes presence rasters
  rockfish_raster <- sum(rockfish_raster, temp)
}

#writeRaster(rockfish_raster, "rockfish_richness.tif", overwrite = TRUE)

# Calculate rockfish+commercial fishing hotspots
files <- list.files(path=here("fishing"), pattern="*.tif", full.names=TRUE, recursive=FALSE)

# Read in first commercial fishing raster
fishing_raster <- raster(files[1])
fishing_raster[is.na(fishing_raster[])] <- 0
fishing_raster <- reclassify(fishing_raster, c(-Inf, 0, 0), right=TRUE)

# Read in each of the other commercial fishing rasters
for (i in c(2:length(files))) {
  temp <- raster(files[i])
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  temp <- reclassify(temp, c(-Inf, 0, 0), right=TRUE)
  
  
  # Add this raster to the cumulative sum of other commercial fishing threat rasters
  fishing_raster <- sum(fishing_raster, temp)
}

#writeRaster(fishing_raster, "commercial_fishing_threat.tif", overwrite = TRUE)

rockfish_commercial_hotspots <- calculate_hotspots(rockfish_raster, fishing_raster)
writeRaster(rockfish_commercial_hotspots, "hotspots_commercial_fish.tif", overwrite = TRUE)

# Calculate rockfish+dem destr fishing hotspots
threats_rf <- raster(here("impact_rec_fish.tif"))
rockfish_rf_hotspots <- calculate_hotspots(rockfish_raster, threats_rf)
writeRaster(rockfish_rf_hotspots, "hotspots_rec_fish.tif", overwrite = TRUE)

# Calculate rockfish+dem destr fishing hotspots
rockfish_cumulative_hotspots <- calculate_hotspots(rockfish_raster, threats_all, "cumlt")
writeRaster(rockfish_cumulative_hotspots, "hotspots_rockfish_cumulative.tif", overwrite = TRUE)



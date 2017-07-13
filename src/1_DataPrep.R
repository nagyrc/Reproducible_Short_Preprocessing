# This script is the first step in the WUI project. 
# Here we import, project, intersect, organize data layers
# Key layers are the Short ignitions, Radeloff WUI product, MTBS data

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(scales)
library(gridExtra)
library(grid)
library(raster)
library(sf)
library(lubridate)
library(ncdf4)
library(doParallel)
library(foreach) 

# Import Shapefiles from Geodatabase -------------------------------------------------------

dir <- "/Users/NateM/Dropbox/Professional/RScripts/"

# An Albers equal area projection for the lower 48
proj <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 "

# Read the feature class
shrt_fire = st_read(dsn = paste0(dir, "Short_Update/data/Short_9215/Data/FPA_FOD_20170508.gdb"),
                    layer = "Fires", quiet= TRUE) %>%
  subset(., STATE != "AK" & STATE != "PR" & STATE != "HI" & FIRE_SIZE >= 0.01) %>%
  select(FPA_ID, ICS_209_INCIDENT_NUMBER, ICS_209_NAME, MTBS_ID, MTBS_FIRE_NAME,
         FIRE_YEAR, DISCOVERY_DOY, STAT_CAUSE_DESCR, FIRE_SIZE, STATE) %>%
  mutate(IGNITION = ifelse(STAT_CAUSE_DESCR == "Lightning", "Lightning", "Human"))

#Import the USA States layer
usa_shp <- st_read(dsn = paste0(dir,  "Short_Update/data/cb_2016_us_state_20m"),
                   layer = "cb_2016_us_state_20m", quiet= TRUE) %>%
  st_transform(., proj) %>%
  subset(., NAME != "Alaska" &
           NAME != "Hawaii" &
           NAME != "Puerto Rico") %>%
  mutate(area_m2 = as.numeric(st_area(geometry)),
         StArea_km2 = area_m2/1000000,
         group = 1) %>%
  st_simplify(., preserveTopology = TRUE)
plot(usa_shp[5])

# Dissolve to the USA Boundary
conus <- usa_shp %>%
  group_by(group) %>%
  st_union()
plot(conus)

# Import the Level 3 Ecoregions
eco = paste0(dir, "Short_Update/data/us_eco_l3")
ecoreg <- st_read(dsn = eco, layer = "us_eco_l3", quiet= TRUE) %>%
  st_transform(., proj) %>%
  st_simplify(., preserveTopology = TRUE, dTolerance = 1000) %>%
  mutate(area_m2 = as.numeric(st_area(geometry)),
         EcoArea_km2 = area_m2/1000000)
plot(ecoreg[2])

# Intersects the region 
state_eco <- st_intersection(usa_shp, ecoreg) %>%
  dplyr::select(STUSPS, NAME, StArea_km2, US_L3CODE, US_L3NAME, NA_L2NAME, NA_L1CODE, NA_L1NAME, geometry)

ign_state_eco <- st_intersection(shrt_fire, state_eco)


# Import Wind -------------------------------------------------------------

wind_dl <- list.files("/Users/NateM/Dropbox/Professional/RScripts/Mietkiewicz_etal_HumanIgnProb/data/raw/historical_gridmet/", pattern = "nc", full.names = TRUE)

netcdf_to_tif_numabove95 <- function(file, mask, outfolder) {
  file_split <- wind_dl %>%
    basename %>%
    strsplit(split = "_") %>%
    unlist
  var <- file_split[1]
  year <- substr(file_split[2], start = 1, stop = 4)
  
  start_date <- as.Date(paste(year, "01", "01", sep = "-"))
  end_date <- as.Date(paste(year, "12", "31", sep = "-"))
  date_seq <- seq(start_date, end_date, by = "1 day")
  monthly_seq <- seq(start_date, end_date, by = "1 month")
  month_seq <- month(date_seq)
  
  nc1 <- nc_open(wind_dl[1])
  nc_att <- attributes(nc1$var)$names
  ncvar <- ncvar_get(nc1, nc_att)
  tvar <- aperm(ncvar, c(1,2,3))
  proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  rbrck1 <- brick(tvar, crs= proj)
  extent(rbrck1) <- c(-124.772163391113, -67.06383005778, 25.0626894632975, 49.3960227966309) #from metadata in Arc or IDRISI, when you open .nc file
  
  nc1 <- nc_open(wind_dl[2])
  nc_att <- attributes(nc1$var)$names
  ncvar <- ncvar_get(nc1, nc_att)
  tvar <- aperm(ncvar, c(1,2,3))
  proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  rbrck2 <- brick(tvar, crs= proj)
  extent(rbrck2) <- c(-124.772163391113, -67.06383005778, 25.0626894632975, 49.3960227966309) #from metadata in Arc or IDRISI, when you open .nc file
  
  mst_stack <- stack(rbrck1, rbrck2)
  
  res <- stackApply(mst_stack, indices = month_seq, fun = mean)
  
  names(res) <- paste(var, unique(month(monthly_seq)),
                      sep = "_")
  masked_res <- mask(res, mask)
  dir.create(paste0(dirc, dir_proc, var), showWarnings = FALSE)
  dir.create(paste0(dirc, dir_proc, var, "/", outfolder), showWarnings = FALSE)
  out <- paste0(dirc, dir_proc, var,  "/", outfolder, "/")
  writeRaster(masked_res, filename = paste0(out,names(masked_res)),
              format = "GTiff", bylayer=TRUE, overwrite = TRUE)
  return(paste("File", names(masked_res), "written"))
}

UseCores <- detectCores() -1
cl       <- makeCluster(UseCores)

ffwi_95 <- foreach(i = 1:length(ffwi_dl)) %dopar% {
  netcdf_to_tif_numabove95(ffwi_dl[i], usa_shp, "NumDaysAbove95th")}

stopCluster(cl)

# Create Fishnets ---------------------------------------------------------
# make an object the size and shape of the output you want
globe_bb <- matrix(c(-2236857, -1691466, 2127027,  1328720), nrow=2) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(., crs = 4326)

# generate grid of 20 x 10 tiles (200 in total, each 18 x 18 degrees)
globe_grid_18x18 <- st_make_grid(globe_bb, n = c(20, 10), 
                                 crs = 4326, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))



r10k <- raster(extent(matrix( c(-2236857, -1691466, 2127027,  1328720), nrow=2)), 
               res = 10000, 
               crs = proj)            
r10k[] <- 1:ncell(r10k)

fish10k <- rasterToPolygons(r10k, fun=NULL, dissolve=TRUE)
fish10k <- SpatialPolygonsDataFrame(fish10k, fish10k@data)
fish10k$fishid10k <- fish10k$layer
fish10k$Area_SqKm <- area(fish10k) / 1000000
fish10k_conus <- intersect(fish10k, CONUS)


r25k <- raster(extent(matrix( c(-2236857, -1691466, 2127027,  1328720), nrow=2)), 
               res = 25000, 
               crs = proj)            
r25k[] <- 1:ncell(r25k)

fish25k <- rasterToPolygons(r25k, fun=NULL, dissolve=TRUE)
fish25k <- SpatialPolygonsDataFrame(fish25k, fish25k@data)
fish25k$fishid25k <- fish25k$layer
fish25k$Area_SqKm <- area(fish25k) / 1000000
fish25k_conus <- intersect(fish25k, CONUS)


r50k <- raster(extent(matrix( c(-2236857, -1691466, 2127027,  1328720), nrow=2)), 
               res = 50000, 
               crs = proj)            
r50k[] <- 1:ncell(r50k)

fish50k <- rasterToPolygons(r50k, fun=NULL, dissolve=TRUE)
fish50k <- SpatialPolygonsDataFrame(fish50k, fish50k@data)
fish50k$fishid50k <- fish50k$layer
fish50k$Area_SqKm <- area(fish50k) / 1000000
fish50k_conus <- intersect(fish50k, CONUS)

StEcoReg_Fish10 <- intersect(StEcoReg, fish10k)
StEcoReg_Fish10_25 <- intersect(StEcoReg_Fish10, fish25k)
StEcoReg_Fish10_25_50 <- intersect(StEcoReg_Fish10_25, fish50k)

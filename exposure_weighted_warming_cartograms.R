library(ncdf4); library(raster); library(cartogram); library(rgdal); library(rgeos)
# weighted median function from spatstat package
source(paste0(wdir, "fun_WeightedMedian.R"))

# INPUT #######################################################################
wdir <- "" #enter working directory

# temperature 
period <- "preind" # OR "18501900"
median.temp <- raster(paste0(wdir, "5yravg.2014-2018.wrt_", period, "_median.asc"))
# total population count
popc.nc <- raster::stack(paste0(wdir, "popc_05x05.nc"))
popc.nc <- stackApply(popc.nc[[165:169]], 1, median, na.rm = TRUE) # take the median population of the last 5 years
# rural population
popc_r <- raster::stack(paste0(wdir, "rurc_05x05.nc"))
popc_r <- stackApply(popc_r[[165:169]], 1, median, na.rm = TRUE)
popc_r[is.na(popc_r)] <- 0
# urban population
popc_u <- raster::stack(paste0(wdir, "urbc_05x05.nc"))
popc_u <- stackApply(popc_u[[165:169]], 1, median, na.rm = TRUE)
popc_u[is.na(popc_u)] <- 0
# urban heat island 
uhi.r <- raster(paste0(wdir, "uhi.asc"))
uhi.r[is.na(uhi.r)] <- 0
# country shapefile
cntr <- shapefile(paste0(wdir, "NE_WORLD_POP2016_Alaska.shp"))


## Create temperature - population table #######################################
median.temp <- mask(median.temp, popc.nc)

popc.v <- rasterToPolygons(popc.nc)
names(popc.v) <- "pop"
median.temp.v <- rasterToPolygons(median.temp)
median.temp.v$popc_median <- popc.v$pop # add population to temperature spatial df

df.temp <- data.frame(temp=median.temp.v[paste0("X5yravg.2014.2018.wrt_", 
                                                period, "_median")], 
                      pop=median.temp.v$popc_median)
df.temp <- df.temp[df.temp$pop!=0,]

write.csv(df.temp, paste0(wdir, "temp_pop_", period, ".csv"), row.names = FALSE)

## create temperature means for each country weighted by population and area ####
# UHI weighted 5-year median temperature
#  T_(to colour country by) = sum[ median_2016*rural_popn + (median_2016+UHI)*urban_popn ] / sum [rural_popn+urban_popn]
temp.uhi <- (median.temp * popc_r + (median.temp+uhi.r) * popc_u) / (popc_r + popc_u)
temp.uhi[is.na(temp.uhi)] <- median.temp[is.na(temp.uhi)] # fill gaps in UHI with non-UHI temperature

# load country shp
cntr <- cntr[!is.na(cntr@data$POP2016),]
cntr <- cntr[cntr@data$POP2016 > 0, ]
cntr$FID <- 1:length(cntr)

# temperature weighted by population for each country country
for (i in 1:length(cntr)){
  mask.poly <- SpatialPolygons(list(cntr@polygons[[i]]))
  p <- extract(x = popc_r + popc_u, mask.poly)
  m <- extract(x = temp.uhi, mask.poly)
  m <- m[[1]]
  p <- p[[1]]
  if (!is.null(p)){
    # calculate weighted median and add to shapefile
    cntr$wmedian_t[cntr$FID==i] <- WeightedMedian(m, p)
  } else {
    m[is.null(m)] <- NA
    # not covered by population - use the value of the grid cell
    cntr$wmedian_t[cntr$FID==i] <- m
  }
}

# area weighted temperature by country
for (i in 1:length(cntr)){
  mask.poly <- SpatialPolygons(list(cntr@polygons[[i]]))
  p <- extract(x = area(median.temp), mask.poly)
  m <- extract(x = median.temp, mask.poly)
  m <- m[[1]]
  p <- p[[1]]
  if (!is.null(p)){
    # calculate weighted mean and add to shapefile
    cntr$warea_t[cntr$FID==i] <- weighted.mean(m, p)
  } else {
    m[is.null(m)] <- NA
    # not covered by population - use the value of the grid cell
    cntr$warea_t[cntr$FID==i] <- m
  }
}

write.csv(cntr@data[,c(9, 97:99)], paste0(wdir, "CNTR_TempArea_", period, ".csv"))

## make cartogram based on population ##########################################
worldR <- readOGR(wdir, layer = "NE_WORLD_POP2016_Alaska")
worldR <- worldR[!is.na(worldR@data$POP2016),]
worldR <- worldR[worldR@data$POP2016 > 0, ]
worldR$FID <- 1:length(worldR)

# update Alaska and remainder US population
worldR@data$POP2016[worldR@data$ADM0_A3=="ALA"] <- 736.341 # in 1k, http://worldpopulationreview.com/states/alaska-population/
worldR@data$POP2016[worldR@data$ADM0_A3=="USA"] <- worldR@data$POP2016[worldR@data$ADM0_A3=="USA"] - 736.341

# project
world.carto <- spTransform(worldR, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "))

# create cartogram
world.carto <- cartogram_cont(x = world.carto, weight = "POP2016", itermax = 3) 

writeOGR(world.carto, dsn = wdir, layer = paste0("cartogram_", period, "_uhi"), 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)


# make cartogram based on country size #########################################
worldR <- readOGR(wdir, layer = "NE_WORLD_POP2016_Alaska")
worldR <- worldR[!is.na(worldR@data$POP2016),]
worldR <- worldR[worldR@data$POP2016 > 0, ]
worldR$FID <- 1:length(worldR)
# project 
worldR <- spTransform(worldR, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "))

# calculate areas
worldR@data$AREA_sqkm <- NA
for(i in 1:length(worldR@data[,1])){
  worldR@data$AREA_sqkm[i] <- gArea(worldR[i,]) / 1000
}

# create cartogram
worldR <- cartogram_cont(x = worldR, weight = "AREA_sqkm", itermax = 3)

writeOGR(worldR, dsn = wdir, layer = "area_cartogram", driver = "ESRI Shapefile", 
         overwrite_layer = TRUE)

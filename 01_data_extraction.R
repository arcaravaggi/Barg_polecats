setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data prep for polecat distribution modelling ####

# Library
library(sp)
library(elevatr)
library(raster)
library(spatstat)
library(rgeos)
library(rgdal)
library(dplyr)
library(tidyr)
library(lme4)
library(reshape2)
library(ggplot2)
library(lubridate)
library(maptools)
library(mapdata)
library(ggmap)
library(maps)

setwd(".")

proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                +towgs84=0,0,0")

# All polecat mortality data are available via Local Environmental Records Centres (LERCs). 
# These centres often track requests and the destination of data, hence we have chosen not
# to provide this data here as information on direct requests to LERCs would benefit the
# biodiverity records community.

roadkill <- readOGR(dsn = ".", layer = "mortalities")
roadkill_buffer <- readOGR(dsn = ".", layer = "mortalitybuffer")
str(roadkill) #structure
head(roadkill)

#---------------------------------Habitat---------------------------------------
# Load Corine raster file
corine <- raster("../Shapefiles/LCM/lcm15_wales_WGS84_geo.tif")
plot(corine)
plot(roadkill, add = T) 

# Combine similar habitat classes 
corine[corine %in% c(6:7)] <- 5 # semi-nat grasslands = neutral, calcareous, acid grasslands and pastures
corine[corine %in% 21] <- 20 # urban cover = urban & suburban
corine[corine %in% 11] <- 8 # wetlands = fen & bog
corine[corine %in% 10] <-9 # heath - heather & heather grassland
corine[corine %in% c(16:19)] <-15 # coastal

# Extract habitat types
hab_pt <- raster::extract(corine, roadkill)

#add hab_pt column to roadkill
roadkill$habitat_pt <- hab_pt
head(roadkill)

# Extract landcover for buffered areas
hab_buffer <- raster::extract(corine, roadkill_buffer, df = TRUE)

# Calculate proportional habitat within buffers
prop.lc <- hab_buffer %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  complete(ID, nesting(lc_type), 
           fill = list(pland = 0)) %>%             # fill in implicit landcover 0s
  spread(lc_type, pland)                  # convert to long format

# Add percent cover columns for the specific habitat types I'm interested in
roadkill$pc_broad <- as.numeric(prop.lc$'1')*100  
roadkill$pc_conif <- as.numeric(prop.lc$'2')*100
roadkill$pc_ag <- as.numeric(prop.lc$'3')*100
roadkill$pc_igrass <- as.numeric(prop.lc$'4') * 100
roadkill$pc_sgrass <- as.numeric(prop.lc$'5') * 100
roadkill$pc_wet <- as.numeric(prop.lc$'8')*100
roadkill$pc_heath <- as.numeric(prop.lc$'9')*100
roadkill$pc_urban <- as.numeric(prop.lc$'20')*100

# Time to figure out number of habitat patches
prop.lc2 <- prop.lc # Make a copy cuz I'm sure I'll screw this up

# Change all values greater than 0 to a 1, excluding the ID column
prop.lc2[ ,2:13][prop.lc2[ ,2:13] > "0.000000000"] <- 1
head(prop.lc2)

# Gonna go ahead and add the row totals as a new column in prop.lc dataframe and save that as a 
# csv so I have it forever because it took me a long time to make and I'm scared to screw it up again.
prop.lc$patches <- rowSums(prop.lc2[ ,2:13])
head(prop.lc)

# Add number of habitat patches to roadkill dataframe
roadkill$patches <- rowSums(prop.lc2[ ,2:13])
head(roadkill)

#---------------------------------Elevation-------------------------------------
# Calculate elevation with Elevatr
roadkill.e <- get_elev_point(roadkill, src = "aws")
roadkill$elevation <- roadkill.e$elevation
head(roadkill)

#---------------------------------Roads-----------------------------------------

# Read OSM road shapefile 
road <- readOGR(dsn = ".", layer = "roadlayer")
plot(road)

# Convert to spatstat object (psp)
road_psp <- as.psp.SpatialLines(road)

# Convert roadkill to spatstat spatial points object (ppp)
roadkill$X <- as.numeric(roadkill$X)
roadkill$Y <- as.numeric(roadkill$Y)
roadkill_ppp <- with(roadkill, ppp(roadkill$X, roadkill$Y, 
                                 owin(range(roadkill$X), range(roadkill$Y))))

# Calculate minimum distance between individual points and nearest road line
road.d <- project2segment(roadkill_ppp, road_psp)
road.distance <- road.d$d

# Add distance to roads to the roadkill dataframe
roadkill$d_road <- road.distance
head(roadkill) 
head(road)

#---------------------------------Road density----------------------------------


# Calculate total line length and density for a given set of polygons
#
# Intersects SpatialLine and SpatialPolygon objects, calculates individual line length,
# appends to dataframe extracted from intersect object and summarises by polygon ID.
# If density is not required, do not specify area parameter `a`.
#
# The function assumes a spatial projection system where distance is given in metres.
#
# E.g.
# sp1 = object of class SpatialLine/SpatialLineDataFrame
# sp2 = object of class SpatialPolygon/SpatialPolygonDataFrame
# id = name of  polygon column for summary. Currently requires a-priori insertion into the function.
# a = area of polygon (default = 0)
#
# E.g.
# df <- len.dens(lines, polygon, a = pi*4.58^2)
spLine.ld  <- function(sp1, sp2, a = 0){
  r <- raster::intersect(sp1, sp2)
  y <- gLength(r, byid = TRUE)
  t <- r@data
  t["r_length"] <- y/1000
  if(a == 0){
    dat <- ddply(t,.(ORIG_FID),summarize,r_length=sum(r_length))
  } else{
    dat <- ddply(t,.(ORIG_FID),summarize,r_length=sum(r_length)) 
    dat["r_dens"] <- dat$r_length / a
  }
  dat
}

crs(roadkill_buffer)
roadkill_buffer <- spTransform(roadkill_buffer, proj.crs)
road <- spTransform(road, proj.crs)

t <- spLine.ld(road, roadkill_buffer, a = pi*1500^2)

# Add density to roadkill
roadkill$road_dens <- t

#---------------------------------Rivers----------------------------------------

# Read OSM waterways shapefile and reproject
water <- readOGR(dsn = ".", layer = "waterlayer")

# Convert roads to spatstat object (PSP = spatstat versionb of SpatialLines)
water_psp <- as.psp.SpatialLines(water)

# Calculate minimum distance between individual points and nearest waterway
water.d <- project2segment(roadkill_ppp, water_psp)

# Add distances to roadkill
roadkill$d_water <- water.d$d
head(roadkill)

#---------------------------------Buildings-------------------------------------

# Read OSM buildings shapefile and reproject
buildings <- readOGR(dsn = ".", layer = "buildingslayer")
head(buildings)
unique(buildings$type)
crs(buildings)

# Calculate distance from buildings
buildings.d <- data.frame((gDistance(buildings, roadkill, byid=TRUE)))

# Calculate MINIMUM distance from buildings and add to alldat
buildings.min.d <- do.call(pmin, buildings.d)
roadkill$d_building <- buildings.min.d
head(roadkill)

# Add an ID column just to help keep myself straight
roadkill$ID <- c(1:171)

# Create 1500m buffer as object psuedobuffer

#---------------------------------Pseudoabsences--------------------------------

# Load data
roadkill <- readOGR(dsn = ".", layer = "mortalitybuffer")
plot(roadkill)
roadkill$ID <- "pseudo"
y <- gUnaryUnion(roadkill, id = roadkill@data$ID)
plot(y)

boundaries <- readOGR(dsn = ".", layer = "boundaries")
boundaries$ct <- "wales"
wales <- gUnaryUnion(boundaries, id = boundaries@data$ct)
plot(wales)
plot(y, add = T)

wales_hole <- gDifference(wales, y)
plot(wales_hole)

samp1 <- spsample(wales_hole, 1200, type = 'random', iter = 25) 
plot(samp1)
plot(y, add = T)

# Create new dataframe for coordinates and merge to create SpatialPointsDataFrame
df <- data.frame(lat = coordinates(samp1)[,1], lon =  coordinates(samp1)[,2])
pseudo <- SpatialPointsDataFrame(samp1, df)


#---------------------------------Habitat&Elevation-----------------------------
# habitat point data
ps_hab_pt <- raster::extract(corine, pseudo)
pseudo$habitat_pt <- ps_hab_pt
head(pseudo)

# Extract landcover for buffered areas
ps_hab_buffer <- raster::extract(corine, pseudo_buffer, df = TRUE) 
ps_hab_buffer <- as.data.frame(ps_hab_buffer) 

# Calculate proportional habitat within buffers
ps.prop.lc = ps_hab_buffer %>%
  setNames(c("ID", "lc_type")) %>%        # rename for ease
  group_by(ID, lc_type) %>%               # group by point (ID) and lc class 
  summarise(n = n()) %>%                  # count the number of occurences of each class
  mutate(pland = n / sum(n)) %>%          # calculate percentage
  ungroup() %>%                           # convert back to original form
  dplyr::select(ID, lc_type, pland) %>%   # keep only these vars
  complete(ID, nesting(lc_type), 
           fill = list(pland = 0)) %>%             # fill in implicit landcover 0s
  spread(lc_type, pland)                  # convert to long format

# Add percent cover columns for the specific habitat types I'm interested in
pseudo$pc_broad <- as.numeric(ps.prop.lc$'1')*100  
pseudo$pc_conif <- as.numeric(ps.prop.lc$'2')*100
pseudo$pc_ag <- as.numeric(ps.prop.lc$'3')*100
pseudo$pc_igrass <- as.numeric(ps.prop.lc$'4') * 100
pseudo$pc_sgrass <- as.numeric(ps.prop.lc$'5') * 100
pseudo$pc_wet <- as.numeric(ps.prop.lc$'8')*100
pseudo$pc_heath <- as.numeric(ps.prop.lc$'9')*100
pseudo$pc_urban <- as.numeric(ps.prop.lc$'20')*100

# Time to figure out number of habitat patches
ps.prop.lc2 <- ps.prop.lc # still gonna make a copy even though it might be unnecessary at this point

# Change all values greater than 0 to a 1, excluding the ID column
ps.prop.lc2[ ,2:13][ps.prop.lc2[ ,2:13] > "0.000000000"] <- 1
head(ps.prop.lc2)

# Add number of habitat patches to roadkill dataframe
pseudo$patches <- rowSums(ps.prop.lc2[ ,2:13])
head(pseudo)

# Calculate elevation with Elevatr
pseudo.e <- get_elev_point(pseudo, src = "aws")
pseudo$elevation <- pseudo.e$elevation
head(pseudo)

#---------------------------------Roads&Rivers$Buildings------------------------
# Convert pseudo to spatstat spatial points object (ppp)
pseudo_ppp <- with(pseudo, ppp(pseudo$lat, pseudo$lon, 
                                   owin(range(pseudo$lat), range(pseudo$lon))))

# Calculate minimum distance between individual points and nearest road line
ps.road.d <- project2segment(pseudo_ppp, road_psp)
ps.road.distance. <- ps.road.d$d

# Add distance to roads to the alldat dataframe
pseudo$d_road <- ps.road.distance

# Road density
pseudo_buffer <- spTransform(pseudo_buffer, proj.crs)
t2 <- spLine.ld(road, pseudo_buffer, a = pi*1500^2)

# Add density to pseudo
pseudo$road_dens <- t2

# Calculate minimum distance between individual points and nearest waterway
ps.water.d <- project2segment(pseudo_ppp, water_psp)

# Add distances to pseudo
pseudo$d_water <- ps.water.d$d
head(pseudo)

# Calculate distance from buildings
ps.buildings.d <- data.frame((gDistance(buildings, pseudo, byid=TRUE)))

# Calculate MINIMUM distance from buildings and add to alldat
ps.buildings.min.d <- do.call(pmin, ps.buildings.d)

head(pseudo)

#--------------------------------Rebuild for modelling--------------------------
# Reorganise pseudoabsences 
pseudo <- pseudo[, c("X", "Y", "habitat_pt", "pc_grass", "pc_wet", "pc_ag", "pc_urban", "pc_conif",
                     "pc_broad", "pc_mixed", "pc_humans", "patches", "elevation", "d_road", "d_water", "d_building", 
                     "fclass", "all_mt_", "presence", "moran")]

# Do a quick reorganise on my roadkill data 
roadkill$X <- roadkill$lat
roadkill$Y <- roadkill$long

roadkill <- roadkill[, c("X", "Y", "habitat_pt", "pc_grass", "pc_wet", "pc_ag", "pc_urban", "pc_conif",
                         "pc_broad", "pc_mixed", "pc_humans", "patches", "elevation", "d_road", "d_water", "d_building", 
                         "fclass", "all_mt_", "presence", "moran")]

all.data <- rbind(roadkill, pseudo)
summary(all.data)

write.csv(all.data, "polecat_model_data.csv")

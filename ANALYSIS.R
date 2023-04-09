################################################################################
################# FULL SCRIPT FOR FIELD ET AL. 2021 ############################
################################################################################

################################################################################
####################### SETTING UP FOR ANLAYSIS ################################
################################################################################
#### IMPORT NECESSARY LIBRARIES 
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(gdistance)
library(leastcostpath)
library(devtools)
library(sf)
library(spatialEco)
library(magrittr)

# set working directory and set master projection if re-projecting
#setwd()

master.projection <- sp::CRS('+proj=utm +datum=NAD83 +zone=13')


################################################################################
######################## ANALYSIS FOR CORRIDOR #1 ##############################
################################################################################

###########################################################
################## STEP #1: IMPORT DATA ###################
###########################################################
# Import data and 
# aggregate DEM to reduce computational costs 
origins <- readOGR(dsn="./DATA/origin locations", layer="origin locations")
# NOTE: Insert your own destinations and waypoints for each corridor
sites <- readOGR(dsn="./DATA/Sites",layer="USER_SPECIFIC_SITE_DATA")
destination1 <- subset(sites,Name == "Corridor 2 destination")
destination2 <- subset(sites,Name == "Corridor 3 destination")

dem <- raster('./DATA/DEM/DEM_corridor1.tif',na.rm=T)

# data downloaded through repository has already been aggregated via following steps 
#dem.aggregate <- 5
#dem <- aggregate(dem, fact=dem.aggregate, fun=mean)

###########################################################
####### STEP #2.1: CREATE CONDUCTANCE SURFACES  & #########
########### CALCULATE TERRAIN RUGGEDNESS ##################
###########################################################
# Calculate terrain ruggedness index
# this calculates mean difference between center cell and near-by cells
tri_dem <- spatialEco::tri(dem)

# List of time functions available in "leastcostpath" R-package (Lewis 2021)
cfs <- c("tobler", "tobler offpath", "irmischer-clarke male", 
         "irmischer-clarke offpath male", "irmischer-clarke female", 
         "irmischer-clarke offpath female","modified tobler",
          "wheeled transport", "herzog", "llobera-sluckin", "campbell 2019")

# Set the number of neighbors considered in creating transition matrix
# The greater the number, the more representative of human movement, but 
# requires more computational power
neigh <- 16

#### Create time-based transition layers 
# NOTE: Transition layers can be transformed to a raster layer with "raster" function 
# and exported with "writeRaster" 
tobler_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "tobler", neighbours = neigh)
ic_off_m_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "irmischer-clarke offpath male", neighbours = neigh)
campbell_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "campbell 2019", neighbours = neigh)

#### Create energy-based transition layers 
herzog_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "herzog", neighbours = neigh)
ls_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "llobera-sluckin", neighbours = neigh)

###########################################################
########### STEP #2.2: CREATE PANDOLF SURFACE #############
###########################################################
# creating a transition layer with Pandolf et al. (1977) cost-function
# requires computation beyond Lewis (2021) script 

# Define inputs for Pandolf et al. (1977) function
W <-63.5 # body mass
L <-20 # load mass (kg)
V <- .35 # velocity (m/s)
N <- 1.2 # terrain

# Create conductance surface using Van Etten (2017) methodology
# First, create function to calculate the attitudinal difference between adjacent 
# cells and apply it to the elevation data. 
altDiff <- function(x){x[2] - x[1]}
hd <- gdistance::transition(dem, altDiff, directions=neigh, symm=FALSE)

# Use the geoCorrection function to divide the attitudinal difference by the 
# distance between cells (i.e., calculating slope as rise over run)
slope <- gdistance::geoCorrection(hd)

# Values between non-adjacent cells is 0, but the slope between these cells is not 0! 
# so therefore, we need to restrict the calculation to adjacent cells. 
# This is done by creating an index for adjacent cells (`adj`) with the 
# function `adjacent` (Van Etten 2017)
adj <- raster::adjacent(dem, cells=1:ncell(dem), pairs=TRUE, directions=neigh)

# Replicate slope for duplicated use in following steps
cost <- slope

# Create Pandolf et al. (1977) function and apply this function to the entire 
# surface  to create cost function 
cost_function_pan <- function(x){ 1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100))) }
cost[adj] <- cost_function_pan(slope)
Conductance_1 <- gdistance::geoCorrection(cost)

# NOTE: Importantly, Pandolf et al. (1977) function creates negative values 
# on negative slopes, so the function needs to be applied differently 
# depending on the slope (Herzog 2010, 2014; White 2012; White and Barber 2012).
# Therefore two conductance surfaces are going to be created, delimited depending 
# on the slope they represent, and then combined 

# Create a second surface using Santee et al. (2001) function, which fixes issue
# of using Pandolf et al. (1977) over negative slopes 
cost_function_santee <- function(x){ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))- N * ((V * (abs(x[adj]) * 100)* (W + L) / 3.5) - ((W + L) * ((abs(x[adj])*100) + 6)^2 / W) + (25-(V^2)))) }
cost[adj] <- cost_function_santee(slope)
Conductance_2 <- gdistance::geoCorrection(cost)

# Determine which portions of DEM are positive and negative slope by classifying 
# rasters. This will produce two rasters where cells equal 1 (positive slopes) and
# 0 (negative slopes), OR vice versa. 
# NOTE: ensure "0" of classified rasters equal N/A
slope_r <- raster(slope)
slope_r[is.na(slope_r)] <- maxValue(slope_r)
positiveslope <- slope_r >=0
negativeslope <- slope_r < 0
negativeslope[negativeslope == 0] <- NA
positiveslope[positiveslope==0] <- NA

# Create barrier out of positive and negative raster using function from Lewis (2021)
positiveslopesonly <- leastcostpath::create_barrier_cs(raster = dem, barrier = negativeslope, neighbours = neigh, field = 0, background = 1)
negativeslopeonly <- leastcostpath::create_barrier_cs(raster = dem, barrier = positiveslope, neighbours = neigh, field = 0, background = 1)

# Multiply Pandolf et al. (1977) derived conductance surface by barrier raster 
# and multiply Santee et al. (2001) derived conductance surface by other barrier 
# raster. 
# NOTE: Make sure positive and negative slopes cells are represented by 1 
# to ensure you are note distorting the conductance values.
pandolfpositive_conductance <- Conductance_1*positiveslopesonly
pandolfnegative_conductance <- Conductance_2*negativeslopeonly

# Merge the two delimited conductance surfaces 
pandolf_cs <- pandolfpositive_conductance+pandolfnegative_conductance

#### remove unnecessary variables 
rm(hd)
rm(slope)
rm(adj)
rm(cost)
rm(pandolfnegative_conductance)
rm(pandolfpositive_conductance)
rm(positiveslope)
rm(negativeslope)
rm(positiveslopesonly)
rm(negativeslopeonly)
rm(Conductance_1)
rm(Conductance_2)
rm(L)
rm(N)
rm(V)
rm(W)



###########################################################
########### STEP #3: CALCULATE LEAST COST PATHS ###########
########### FROM SERIES OF ORIGIN POINTS TO  ##############
########### SEPARATE DESTINATIONS AND COMPUTE #############
########### DETERMINING VARIABLES (TERRAIN, DISTANCE) #####
###########################################################
# SEE EXAMPLE: "Corridor #1" case study Field et al. (2021)

# Define length of loop by number of unique origin locations 
n <- length(origins)
# Build table to store results
result_destination1 <-as.data.frame(matrix(NA,nrow=n,ncol=9))
colnames(result_destination1)<- c("time (sq.m)","energy(sq.m)","overlap(sq.m)","mean time and energy (sq.m)", "overlap as percent of average",
                                  "mean terrain ruggedness within 1 km of time lcps","median terrain ruggedness within 1 km of time lcps",
                                  "mean terrain ruggedness within 1 km of energy lcps","median terrain ruggedness within 1 km of energy lcps" )
# Execute loop that 
# 1: Calculates LCPs from origin to destination
# 2: Calculates distance between all paths and nearby terrain for time- and 
# energy=based LCP corridors
# 3: Calculates summary statistics for each iteration

for(i in 1:n){
  origin <- origins[i,]
  proj4string(origin) <- sp::CRS("+proj=utm +datum=NAD83 +zone=13")
  # Calculate least cost paths 
  
  t_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)
  ic_lcp <- leastcostpath::create_lcp(cost_surface = ic_off_m_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)
  c_lcp <- leastcostpath::create_lcp(cost_surface = campbell_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)
  h_lcp <- leastcostpath::create_lcp(cost_surface = herzog_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)
  ls_lcp <- leastcostpath::create_lcp(cost_surface = ls_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)
  p_lcp <- leastcostpath::create_lcp(cost_surface = pandolf_cs, origin = origin, destination = destination1, directional = T,cost_distance=T)

  # Create points on time-based LCPs at 100 m intervals
  points <- gLength(t_lcp)/100
  t_p<-spsample(t_lcp, n = points, type = "regular")
  points <- gLength(ic_lcp)/100
  ic_p<-spsample(ic_lcp, n = points, type = "regular")
  points <- gLength(c_lcp)/100
  c_p<-spsample(c_lcp, n = points, type = "regular")
  
  # Find common minimum number of points created in previous step and transform 
  # 100 m points in df that is limited to that minimum number (this will remove 
  # some data but will make graphing  and second loop possible)
  length <- data.frame(length(c_p),length(ic_p),length(t_p))
  l <- min(length)
  t_p <- t_p[1:l,]
  ic_p <- ic_p[1:l,]
  c_p <- c_p[1:l,]
  
  # Create table to store values of 
  # 1: distance between each time-based LCP to each energy-based LCP at every 
  # 100 m of travel along time-based LCP
  # 2: distance between LCPs at every 100 m point of travel and distance from 
  # every 100 m point of travel and the final destination
  result_difference <-as.data.frame(matrix(NA,nrow=l,ncol=18))
  colnames(result_difference)<- c("tobler_pandolf","tobler_lloberasluckin","tobler_herzog",
                                  "irmischer_pandolf","irmischer_herzog","irmischer_lloberasluckin",
                                  "campbell_pandolf","campbell_lloberasluckin","campbell_herzog",
                                  "tp_pctdistancetodestination","tls_pctdistancetodestination","th_pctdistancetodestination",
                                  "ip_pctdistancetodestination","ils_pctdistancetodestination","ih_pctdistancetodestination",
                                  "cp_pctdistancetodestination","cls_pctdistancetodestination","ch_pctdistancetodestination")
  
  # Run embedded loop to populate dataframe based on minimum number
  # of 100 m intervals  
   for(j in 1:l){
      result_difference[j,1]<-gDistance(t_p[j,], p_lcp)
      result_difference[j,2]<-gDistance(t_p[j,],ls_lcp)
      result_difference[j,3]<-gDistance(t_p[j,],h_lcp)
      result_difference[j,4]<-gDistance(ic_p[j,], p_lcp)
      result_difference[j,5]<-gDistance(ic_p[j,],ls_lcp)
      result_difference[j,6]<-gDistance(ic_p[j,],h_lcp)
      result_difference[j,7]<-gDistance(c_p[j,], p_lcp)
      result_difference[j,8]<-gDistance(c_p[j,],ls_lcp)
      result_difference[j,9]<-gDistance(c_p[j,],h_lcp)
      result_difference[j,10:12]<-(abs(gDistance(t_p[j,], destination1)-gLength(t_lcp))/gLength(t_lcp))
      result_difference[j,13:15]<-(abs(gDistance(ic_p[j,], destination1)-gLength(ic_lcp))/gLength(ic_lcp))
      result_difference[j,16:18]<-(abs(gDistance(c_p[j,], destination1)-gLength(c_lcp))/gLength(c_lcp))
      }

  # Export table 
  names <- paste0("result_difference_1_", i)
  write.csv(result_difference,paste0("./OUTPUT",names,".csv"), row.names = FALSE)
  
  # Bind all time-based and energy-based LCPs together 
  t_all <- rbind(t_lcp,ic_lcp,c_lcp)
  e_all <- rbind(h_lcp,ls_lcp,p_lcp)
  
  # Create a corridor of time-efficient and energy-efficient travel by applying
  # a 125 m buffer to the binded paths (equal to a 500 m area around path).
  # Export the paths as a spatialpolygondataframe
  t_all_buff <- gBuffer(t_all,width=125)
  e_all_buff <- gBuffer(e_all,width=125)
  time_id <- sapply(slot(t_all_buff, "polygons"), function(x) slot(x, "ID"))
  t_all_df <- data.frame( ID=1:length(t_all_buff), row.names = time_id)
  t_buffer_spdf <- SpatialPolygonsDataFrame(t_all_buff,t_all_df)
  names<- paste0("all_time_buff_origin",i)
  writeOGR(t_buffer_spdf,"./OUTPUT",layer=names,driver="ESRI Shapefile")
  energy_id <- sapply(slot(e_all_buff, "polygons"), function(x) slot(x, "ID"))
  e_all_df <- data.frame( ID=1:length(e_all_buff), row.names = energy_id)
  e_buffer_spdf <- SpatialPolygonsDataFrame(e_all_buff,e_all_df)
  names<- paste0("all_energy_buff_origin",i)
  writeOGR(e_buffer_spdf,"./OUTPUT",layer=names,driver="ESRI Shapefile")
  
  # Calculate the overlap between the time- and energy-efficient travel corridors
  # adn export as a spatial polygondataframe
  t_e_overlap <- gIntersection(t_all_buff, e_all_buff)
  overlap_id <- sapply(slot(t_e_overlap, "polygons"), function(x) slot(x, "ID"))
  overlap_df <- data.frame( ID=1:length(t_e_overlap), row.names = overlap_id)
  overlap_spdf <- SpatialPolygonsDataFrame(t_e_overlap,overlap_df)
  names<- paste0("overlap_buff_origin",i)
  writeOGR(overlap_spdf,"./OUTPUT",layer=names,driver="ESRI Shapefile")
  
  # Extract all terrain ruggedness index values within the LCP corridors
  t_buff_tri_extract <- extract(tri_dem,t_buff_tri)
  e_buff_tri_extract <- extract(tri_dem,e_buff_tri)
  
  # Bind all values into a table and export table 
  e_tri <- as.data.frame(e_buff_tri_extract)%>%
    dplyr::mutate(cost="energy")%>%
    magrittr::set_names(c("value", "cost"))
  t_tri <- as.data.frame(t_buff_tri_extract)%>%
    dplyr::mutate(cost="time")%>%
    magrittr::set_names(c("value", "cost"))
  all_tri <- rbind(t_tri,e_tri)
  names <- paste0("result_tri_1_", i)
  write.csv(all_tri,paste0("./OUTPUT",names,".csv"), row.names = FALSE)
  
  # Calculate general statistics for each loop
  result_destination1[i,1]<-area(t_all_buff)
  result_destination1[i,2]<-area(e_all_buff)
  result_destination1[i,3]<-area(t_e_overlap)
  result_destination1[i,4]<-(area(t_all_buff)+area(e_all_buff))/2
  result_destination1[i,5]<-area(t_e_overlap)/((area(t_all_buff)+area(e_all_buff))/2)
  result_destination1[i,6]<- mean(t_buff_tri_extract[[1]])
  result_destination1[i,7]<- median(t_buff_tri_extract[[1]])
  result_destination1[i,8]<- mean(e_buff_tri_extract[[1]])
  result_destination1[i,9]<- median(e_buff_tri_extract[[1]])
}

# Export the final table after loop is complete
write.csv(result_destination1,"./OUTPUT/results_destination_1.csv", row.names = FALSE)

# Same loop can be ran for multiple destinations
# NOTE: Be sure to recall the variable "destination" depending on 
# your new destination point, and change location or name of outputs


# The same set of analyses can be used to determine LCPs and compare differences 
# between them in contexts where you are modelling movement from one site to the 
# next in a stepwise manner until reaching a final destination. 
# SEE EXAMPLE: "Corridor #2" case study in Field et al. (2021)

# For example, you might model a time- based LCP from an "origin" to an "ultimate destination" 
origin <- subset(sites,Name == "Corridor 2 origin")
destination <- subset(sites,Name == "Corridor 2 destination")
waypoint1 <- subset(sites,Name == "Corridor 2 waypoint 1")
waypoint2 <- subset(sites,Name == "Corridor 2 waypoint 2")
waypoint3 <- subset(sites,Name == "Corridor 2 waypoint 3")


# EXAMPLE BELOW:
t_o_d_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = origin, destination = destination, directional = T,cost_distance=T)

# And compare it to other time-based models where movement is determined from  
# the origin to another site, then from that site to a third site, and so on,  
# until reaching the destination. 
# EXAMPLE BELOW: 
t_o_wp1_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = origin, destination = waypoint1, directional = T,cost_distance=T)
t_wp1_wp2_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = waypoint1, destination = waypoint2, directional = T,cost_distance=T)
t_wp2_wp3_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = waypoint2, destination = waypoint3, directional = T,cost_distance=T)
t_wp3_d_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = waypoint3, destination = destination, directional = T,cost_distance=T)

# This process can be repeated for an energy-based LCP and then the two can be 
# compared using the same processes [lines 265-315] applied in the looped portion of "CORRIDOR #1". 


################################################################################
################# STEP 4: INTEGRATING ROADS AS CONDUITS ########################
################################################################################
# Parts of the built environment, such as roads, can also be integrated into the 
# conductance surfaces to simulate increased ease of travel over or across those 
# features. 
# SEE EXAMPLE: "Corridor #3" case study in Field et al. (2021)

origin <- subset(sites,Name == "Corridor 3 origin")
destination <- subset(sites,Name == "Corridor 3 destination")
waypoint1 <- subset(sites,Name == "Corridor 3 waypoint 1")
waypoint2 <- subset(sites,Name == "Corridor 3 waypoint 2")


# Import the raster surface that represents conductive features 
conduit <- raster('./DATA/DEM/DEM_corridor3_conduit.tif',na.rm=T)
dem <- raster('./DATA/DEM/DEM_corridor3.tif',na.rm=T)

# Make sure the raster representing the conduit is the same resolution as the DEM 
# you are using to create the conductance surface 
conduit <- resample(conduit,dem,method='bilinear')

# Ensure raster values represent conduit surface (i.e. value of 1 for the 
# location of the feature and a value of 0 for all places where the feature is not 
# located).  Then, subset the conduit raster so only the conductive surface is 
# represented with numerical values. 
# Note: Change min and max depending on the raster values in the surface 
# labeled "conduit"
min <- 1.5
max <- 1.9999
conduit <- conduit > min & conduit < max
conduit[conduit == 0] <- NA

# Create a conductance surface
tobler_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "tobler", neighbours = neigh)

# Convert surface to raster and calculate range of values in conductance surface.
# Divide range by the amount you wish to increase conductance by. In this 
# instance, conductance will be increased by ten percent. To increase conductance 
# by 50 %, divide by two, etc.  
consurf<-raster(tobler_cs)
max <- maxValue(consurf)
min <- minValue(consurf)
con <- (max - min)/10

# Create conductance surface where ten percent increase in conductance range 
# is applied to the road surface, and nothing to all other locations. 
conduit_cs <- leastcostpath::create_barrier_cs(raster = conduit, barrier = conduit, neighbours = neigh, field = con, background = 0)

# Add the conduit surface to the base conductance surface to create combined  
# surface that incorporates cost across all of the surface, and a reduced cost 
#across the road surface. 
t_cs_withconduit <- tobler_cs+conduit_cs

# Use this surface to model your LCPs.
# EXAMPLE BELOW: 
t_lcp_ex <- leastcostpath::create_lcp(cost_surface = tobler_cs_withconduit, origin = origin, destination = waypoint1, directional = T,cost_distance=T)


################################################################################
################# STEP 5: MODELING LCPs WITH ALTERNATING TIME OR ###############
################# ENERGY CONSIDERATIONS ########################################
################################################################################
# Import data and create conductance surfaces
# Note: The example below demonstrates the method with multiple time- and 
# energy- based LCPs.
destination <- subset(sites,Name == "Corridor 2 destination")
origin <- origins[1,]
dem <- raster('./DATA/DEM/DEM_corridor1.tif',na.rm=T)

# Create an accumulated cost surfaces for each time-based conductance surface 
# that you have created. An accumulated cost surface represents the cumulative 
# cost of travel away from the origin. Even though this method uses conductance, 
# values representing time and energy cost are still preserved in the matrices, 
# so the "accCOst" function is tabulating the cost of travel and not the conductance.

# Note: All cost functions used her calculate time in seconds, however not all cost 
# functions do, so ensure the tabulated costs are constant between surfaces.
acc_tobler_cs<- accCost(tobler_cs,origin)
acc_ic_cs<- accCost(ic_off_m_cs,origin)
acc_campbell_cs<- accCost(campbell_cs,origin)

# Artificial thresholds are used here to determine when movers transition from 
# time to energy-based decision making when moving. See Field et al. (2021) for 
# discussion on these thresholds.

# Limit each time-based surface to the time threshold (30 min/1800 sec) when fatigue
# first sets in (i.e., first fatigue threshold). 
tobler_limit<- acc_tobler_cs<1800
ic_limit<- acc_ic_cs<1800
campbell_limit<- acc_campbell_cs<1800

# Create a contour around the limited surfaces to define the perimeter of travel
# within the first fatigue threshold. 
toblerlimit_l<-rasterToContour(tobler_limit)
toblerlimit_l <- toblerlimit_l[1,]
iclimit_l<-rasterToContour(ic_limit)
iclimit_l <- iclimit_l[1,]
campbelllimit_l<-rasterToContour(campbell_limit)
campbelllimit_l <- campbelllimit_l[1,]

# Calculate LCP from origin to final destination
t_lcp <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = origin, destination = destination, directional = T,cost_distance=T)
ic_lcp <- leastcostpath::create_lcp(cost_surface = ic_off_m_cs, origin = origin, destination = destination, directional = T,cost_distance=T)
c_lcp <- leastcostpath::create_lcp(cost_surface = campbell_cs, origin = origin, destination = destination, directional = T,cost_distance=T)

# Find out when the full LCP intersects with the first fatigue threshold
tobler_end1 <-rgeos::gIntersection(toblerlimit_l, t_lcp)
ic_end1 <- rgeos::gIntersection(iclimit_l, ic_lcp)
c_end1 <- rgeos::gIntersection(campbelllimit_l, c_lcp)

# If creating multiple paths from different conductance surfaces, the first fatigue 
# threshold and the point of intersection will likely be different for each.
# Therefore to limit the computational costs in the following steps, the average 
# point of intersection was identified and used.
coordinates1 <- coordinates(tobler_end1)
coordinates2 <- coordinates(ic_end1)
coordinates3 <- coordinates(c_end1)
end1 <- c(coordinates1[1,1],coordinates2[1,1],coordinates3[1,1])
end2 <- c(coordinates1[1,2],coordinates2[1,2],coordinates3[1,2])
end1<-mean(end1)
end2<-mean(end2)
end_first <- merge(end1,end2)
end_first <- SpatialPointsDataFrame(coords = end_first, data = end_first)
proj4string(end_first) <- CRS("+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs")


# Re-calculate time-based LCPs from origin to average point of intersection
t_lcp_first <- leastcostpath::create_lcp(cost_surface = tobler_cs, origin = origin, destination = end_first, directional = T,cost_distance=T)
ic_lcp_first <- leastcostpath::create_lcp(cost_surface = ic_off_m_cs, origin = origin, destination = end_first, directional = T,cost_distance=T)
c_lcp_first <- leastcostpath::create_lcp(cost_surface = campbell_cs, origin = origin, destination = end_first, directional = T,cost_distance=T)


# Calculate first set of energy paths from average point of intersection to 
# final destination 
h_lcp <- leastcostpath::create_lcp(cost_surface = herzog_cs, origin = end_first, destination = destination, directional = T,cost_distance=T)
ls_lcp <- leastcostpath::create_lcp(cost_surface = ls_cs, origin = end_first, destination = destination, directional = T,cost_distance=T)
p_lcp <- leastcostpath::create_lcp(cost_surface = pandolf_cs, origin = end_first, destination = destination, directional = T,cost_distance=T)


# Limit accumulated cost surfaces (calculated from the origin) 
# to the second time threshold (90 min/5400 sec) of fatigue
tobler_limit<- acc_tobler_cs<5400
ic_limit<- acc_ic_cs<5400
campbell_limit<- acc_campbell_cs<5400

# Create a contour around the limited surfaces to define the perimeter of travel
# within the second fatigue threshold.
toblerlimit_l<-rasterToContour(tobler_limit)
toblerlimit_l <- toblerlimit_l[1,]
iclimit_l<-rasterToContour(ic_limit)
iclimit_l <- iclimit_l[1,]
campbelllimit_l<-rasterToContour(campbell_limit)
campbelllimit_l <- campbelllimit_l[1,]

# Find out when the energy-based full LCP intersects with the second fatigue 
# threshold.
t_end2_1 <-rgeos::gIntersection(toblerlimit_l, h_lcp)
t_end2_2 <-rgeos::gIntersection(toblerlimit_l, ls_lcp)
t_end2_3 <-rgeos::gIntersection(toblerlimit_l, p_lcp)
ic_end2_1 <- rgeos::gIntersection(iclimit_l, h_lcp)
ic_end2_2 <- rgeos::gIntersection(iclimit_l, ls_lcp)
ic_end2_3 <- rgeos::gIntersection(iclimit_l, p_lcp)
c_end2_1 <- rgeos::gIntersection(campbelllimit_l, h_lcp)
c_end2_2 <- rgeos::gIntersection(campbelllimit_l, ls_lcp)
c_end2_3 <- rgeos::gIntersection(campbelllimit_l, p_lcp)


# Calculate average point of intersection 
coordinates1 <- coordinates(t_end2_1)
coordinates2 <- coordinates(t_end2_2)
coordinates3 <- coordinates(t_end2_3)
coordinates4 <- coordinates(ic_end2_1)
coordinates5 <- coordinates(ic_end2_2)
coordinates6 <- coordinates(ic_end2_3)
coordinates7 <- coordinates(c_end2_1)
coordinates8 <- coordinates(c_end2_2)
coordinates9 <- coordinates(c_end2_3)
end1 <- c(coordinates1[1,1],coordinates2[1,1],coordinates3[1,1],
          coordinates4[1,1],coordinates5[1,1],coordinates6[1,1],
          coordinates7[1,1],coordinates8[1,1],coordinates9[1,1])
end2 <- c(coordinates1[1,2],coordinates2[1,2],coordinates3[1,2],
          coordinates4[1,2],coordinates5[1,2],coordinates6[1,2],
          coordinates7[1,2],coordinates8[1,2],coordinates9[1,2])
end1<-mean(end1)
end2<-mean(end2)
end_second <- merge(end1,end2)
end_second <- SpatialPointsDataFrame(coords = end_second, data = end_second)
proj4string(end_second) <- CRS("+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs")


# Re-calculate energy-based LCPs from first the point representing the average first
# fatigue threshold to point representing the average second fatigue threshold
h_lcp_second <- leastcostpath::create_lcp(cost_surface = herzog_cs, origin = end_first, destination = end_second, directional = T,cost_distance=T)
ls_lcp_second <- leastcostpath::create_lcp(cost_surface = ls_cs, origin = end_first, destination = end_second, directional = T,cost_distance=T)
p_lcp_second <- leastcostpath::create_lcp(cost_surface = pandolf_cs, origin = end_first, destination = end_second, directional = T,cost_distance=T)

# This process can be reiterated for as many series as necessary to get to the 
# final destination. 
# NOTE: After completing the first two legs of travel new accumulated cost surfaces 
# must be constructed, so that the next two fatigue thresholds are representative 
# of cost away from the point of travel and not from the original point of origin.

# EXAMPLE BELOW: For calculating the next phases of travel, calculate a new accumulated 
# surface using the same conductance surface created at the beginning of the process 
# and the point representing the average second fatigue threshold. 
acc_tobler_cs<- accCost(tobler_cs,end_second)



################################################################################
############################## CITATIONS #######################################
################################################################################
# Field, Sean, Donna M. Glowacki, and Lee Gettler 2021. Energetics in Least Cost
# Analysis. Submitted for review in the Journal of Archaeological Method and Theory

# Herzog, Irmela. 2010. Theory and Practice of Cost Functions. Paper presented at the 38th Annual
# Conference on Computer Applications and Quantitative Methods in Archaeology. Granada,
# Spain.

# Herzog, Irmela. 2014. A Review of Case Studies in Archaeological Least-Cost Analysis.
# Archaeologia e Calcolatori, 25: 223-239.

# Lewis, J. 2021a. leastcostpath: Modelling Pathways and Movement Potential Within a Landscape (version 1.8.2). 
# Available at: https://cran.r-project.org/web/packages/leastcostpath/index.html

# Pandolf, Kent B., B. Givoni and R.F. Goldman. 1977. Predicting energy expenditure with loads
# while standing or walking very slowly. Journal of Applied Physiology, 43: 577-581.

# Santee, W.R., W.F. Allison, L.A. Blanchard, and M.G. Small. 2001. A Proposed Model for Load 
# Carriage on Sloped Terrain. Aviation, Space and Environmental Medicine 72(6):562-566.

# Van Etten, Jacob. 2017. R package gdistance: distances and routes on geographical grids. 
# Journal of Statistical Software 76(1):1-21.

# White, Devin A. 2012. Prehistoric Trail Networks of the Western Papagueria: A Multifaceted
# Least Cost Graph Theory Analysis. In: Least Cost Analysis of Social Landscapes: Archaeological
# Case Studies, edited by Devin A. White and Sarah Surface-Evans, 188-208. Salt Lake City, UT:
# University of Utah Press.

# White, Devin A. and Sarah B. Barber. 2012. Geospatial modeling of pedestrian transportation
# networks: a case study from pre-Columbian Oaxaca, Mexico. Journal of Archaeological Science
# 39(8):2684-2696.


library(lidR)
library(sf)
library(terra)
library(viridis)
library(dplyr)
 
###########################
###  First Steps ###########
###########################
 
#Import the city districts 
boundary <- sf::st_read("C:/Users/LaraO/EAGLE_Master/2_Semester/LiDAR-Abgabe/AOIs/Boundary_Lindleinsmuehle.gpkg", quiet = TRUE)
 
#Import the .laz files for the city districts using the LAScatalog processing engine 
  #Since the city districts are consisting of several LAS files the LAS catalog is used to extract the region of interest without 
  #invalidity between the tiles (r-lidar-github)
 
ctg <- readLAScatalog("C:/Users/LaraO/EAGLE_Master/2_Semester/LiDAR-Abgabe/LiDAR_Daten/lindleinsmuehle") 
 
#Set the coordinate Reference System 
projection (ctg) <- 25832
print(ctg) #checking if it worked
 
#Clip the catalog to the boundaries of the city districts 
clipped <-clip_roi(ctg, boundary)
plot(clipped)
print(clipped)
 
 
#Tabulate Classification 
table(clipped$Classification)
plot(clipped, color="Classification")
 
 
###########################################
#### Height Normalization ###########
###########################################
 
# In the first step of the analysis the elevations of the point clouds have to be transformed from absolute 
# elevations above sea level to heights relative to the ground. The Height normalization is archieved with the
#  Triangular irregular network (tin()) algorithm. 
 
nlas <- normalize_height(clipped, tin())
 
# Check if the point cloud normalization worked out (all ground points should be exactly 0)
hist(filter_ground(nlas)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
          #They are! 
 
#################################################################
##### Creation of a Canopy Height Model without buildings ####
#################################################################
 
# Since the ROIs are located in urban settings the point clouds have to be cleared from 
# interfering factors as f.ex. buildings to create a canopy height model just consisting of vegetation
# for the following tree detection 

# Filter point clouds for vegetation and buildings
 
veg      <- filter_poi(nlas, Classification == 20 & Z >= 3 & Z <= 35)
building <- filter_poi(nlas, Classification == 6)
 
# Creating a CHM of the vegetation with the p2r() algorithm 
chm_raw <- rasterize_canopy(
  veg, 
  res = 0.5, 
  algorithm = p2r(subcircle = 0.15), 
  pkg = "terra"
)
 
# # Create a building mask with a buffer of 50cm
building_mask <- rasterize_canopy(
  building, 
  res = chm_raw, 
  algorithm = p2r(), 
  pkg = "terra"
)
 
building_buffered <- focal(
  building_mask, 
  w = matrix(1, 3, 3), 
  fun = max, 
  na.rm = TRUE
)
 
building_buffered <- crop(
  extend(building_buffered, chm_raw), 
  chm_raw
)
 
# Smoothing of the CHM 
 
gf <- focalMat(chm_raw, 0.5, "Gauss")
chm_smooth <- focal(
  chm_raw, 
  w = gf, 
  fun = sum, 
  na.rm = TRUE
)
 
# Masking and Clipping
 
# Mask out the buildings 
chm_final <- mask(
  chm_smooth, 
  building_buffered, 
  inverse = TRUE
)
# Crop to the city borders
chm_final <- crop(chm_final, boundary)
chm_final <- mask(chm_final, boundary)
 
plot(
  chm_final, 
  col = viridis(100), 
  main = "CHM without buildings"
)
 
 
 
#####################################################
########## Individual Tree Detection ###############
#####################################################
 
f <- function(x) {x * 0.2 + 5}
heights <- seq(0,30,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

#Locate tree tops (filter for trees taller than 3m)

ttops <- locate_trees(chm_final, lmf(f))
 ttops <- ttops[ttops$Z >3,] 
 
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
 
########################################################
######### Segmentation of the point cloud ##############
########################################################
 
algo <- dalponte2016(chm_final, ttops, th_tree=3) # treshold of 3m to be detected 
veg_las_segmented <- segment_trees(veg, algo)
 
plot(veg_las_segmented, color ="treeID")
 
 
########################################################
############ Crown Area ################################
########################################################

# Crown Metrics and Filtering 

# Define a custom function for maximum crown height 
f_metrics <- function(z) {
  return(list(max_z = max(z)))
}
 
# Calculate concave crown polygons with metrics
tree_metrics <- crown_metrics(
  las  = veg_las_segmented,
  func = ~f_metrics(Z),
  geom = "concave"
)
 
# Convert to sf object and calculate crown area and radius
crowns_sf <- st_as_sf(tree_metrics)
crowns_sf$crown_area <- as.numeric(st_area(crowns_sf))
crowns_sf$crown_radius <- sqrt(crowns_sf$crown_area / pi)

summary(crowns_sf$crown_area)
summary(crowns_sf$crown_radius)
 
# Filter: Every tree with a crown radius below 1m are being removed 
crowns_filtered <- crowns_sf[crowns_sf$crown_radius >= 1.0, ]


######################################################
############# Tree Count & Density Calculation #######
######################################################

n_trees <- nrow(crowns_filtered)
print(paste("Amount of trees:", n_trees))
 
# Compute and map the density of trees with a 10m resolution
ttops_filtered <- ttops[ttops$treeID %in% crowns_filtered$treeID, ]

r <- terra::rast(x = ttops_filtered)
terra::res(r) <- 10
r <- terra::rasterize(x = ttops_filtered, y = r, "treeID", fun = 'count')
plot(r, col = viridis(20), main = "Tree Density(10m)")
 
# Calculate the tree density per ha for the AOI 
scene_area_m2 <- sum(sf::st_area(boundary))
scene_area_ha <- as.numeric(scene_area_m2) / 10000
 
tree_density <- n_trees / scene_area_ha
print(paste("Tree density:", round(tree_density, 1), "Trees / ha"))


#####################################
#### Height data 
##################################

tt_height <- ttops_filtered$Z

height_stats <- data.frame(
  Total_Trees   = length(tt_height),
  Min_Height    = min(tt_height, na.rm = TRUE),
  Max_Height    = max(tt_height, na.rm = TRUE),
  Mean_Height   = mean(tt_height, na.rm = TRUE),
  Median_Height = median(tt_height, na.rm = TRUE),
  Std_Dev       = sd(tt_height, na.rm = TRUE),
  Q25           = quantile(tt_height, 0.25, na.rm = TRUE),
  Q75           = quantile(tt_height, 0.75, na.rm = TRUE)
)

print(height_stats)

###################################################
########### Crown Metrics and Explorative Diagrams
# ##############################
# 
 
# Calculation of crown area (m²) und and crown radius (m) 
crowns_sf$crown_area <- as.numeric(st_area(crowns_sf))
crowns_sf$crown_radius <- sqrt(crowns_sf$crown_area / pi)
 
 
# A) Histogram of crown radius
hist(crowns_filtered$crown_radius, 
     main = "Histogram of crown radius", 
     xlab = "Crown radius (m)", 
     col = "forestgreen", 
     border = "white")
 
# B) Histogram of the tree heights 
hist(ttops_filtered$Z, 
     main = "Histogram of tree height", 
     xlab = "Tree height (m)", 
     col = "brown", 
     border = "white")
 
# C) Scatter plot: Height vs. crown radius
plot(crowns_filtered$crown_radius ~ crowns_filtered$max_z, 
    main = "Crown Radius vs. Tree Height",
     xlab = "Tree height (m)", 
     ylab = "Crown radius (m)", 
     pch = 16, 
     col = "#00000033")
 




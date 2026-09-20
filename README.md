# Urban Green Canopy Structure Across the City Districts of Würzburg
Final Assignement for the course Novel Image Analysis Methods/ LIDAR Remote sensing of the M.Sc. Applied Earth Observation and Geoanalysis(Matr.3248301)

## Introduction 


## Data 
The project used airbone LiDAR data from the Bayerische Vermessungsverwaltung. The data is freely available and can be downloaded at the [Bayerische Vermessungsverwaltung](https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=laserdaten). The data is is downloaded in tiles with a resolution of 1x1km. The tiles covering the areas studies were downloaded accordingly and had previously been cropped in QGIS using the Point Cloud data management tool: Clip. 

## Workflow 

### Getting the project started
#### Needed Packages
```R
library(lidR)
library(terra)  
library(sf)    
library(dplyr)
library(ggplot2)
library(viridis)

```
#### Adjust the Region of Interest (ROI)

To start with the analyis the .laz files have to be clipped accordingly to the city districts being analyzed. Herefore the boundaries of the city districs as well as the .laz files are needed. The boundaries have been previously extracted using an OSM query in QGIS. The boundaries are available as a Geopackage. Since the city districts are consisting of several LAS files the LAS catalog is used to extract the region of interest without invalidity or edge artifacts between the tiles. 

```R
#Import the city districts 
boundary <- sf::st_read("C:/Users/LaraO/EAGLE_Master/2_Semester/LiDAR-Abgabe/AOIs/Boundary_Lindleinsmuehle.gpkg", quiet = TRUE)

#Import the .laz files for the city districts using the LAScatalog processing engine 
ctg <- readLAScatalog("C:/Users/LaraO/EAGLE_Master/2_Semester/LiDAR-Abgabe/LiDAR_Daten/lindleinsmuehle")

#Set the coordinate Reference System 
projection (ctg) <- 25832
print(ctg) #checking if it worked

#Clip the catalog to the boundaries of the city districts 
clipped <-clip_roi(ctg, boundary)
plot(clipped)
print(clipped)

```
#### Adjust the Region of Interest (ROI)


```R
#Tabulate Classification 
table(clipped$Classification)
plot(clipped, color="Classification")
```
### Height Normalization 

In the first step of the analysis the elevations of the point clouds have to be transformed from absolute elevations above sea level to heights relative to the ground. The Height normalization is archieved with the Triangular irregular network (tin()) algorithm.

```R
nlas <- normalize_height(clipped, tin())

# Check if the point cloud normalization worked out (all ground points should be exactly 0)
hist(filter_ground(nlas)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")
          #They are! 
```
### Creation of a Canopy Height Model without buildings
Since the ROIs are located in urban settings the point clouds have to be cleared from interfering factors as f.ex. buildings to create a canopy height model just consisting of vegetation for the following tree detection 

```R
# Filter Points 
veg      <- filter_poi(nlas, Classification == 20 & Z >= 2 & Z <= 35) 
building <- filter_poi(nlas, Classification == 6)

#Create a CHM of the vegetation with the p2r() algorithm
chm_raw <- rasterize_canopy(
  veg, 
  res = 0.5, 
  algorithm = p2r(subcircle = 0.15), 
  pkg = "terra"
)

# Create a building mask with a buffer of 50cm
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

#Masking and Clipping
chm_final <- mask(
  chm_smooth, 
  building_buffered, 
  inverse = TRUE
)

chm_final <- crop(chm_final, boundary)
chm_final <- mask(chm_final, boundary)

plot(
  chm_final, 
  col = viridis(100), 
  main = "CHM without buildings (NA)"
)

```

### Individual Tree Detection 

```R
f <- function(x) {x * 0.1 + 3}
heights <- seq(0,30,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))
 
ttops <- locate_trees(chm_final, lmf(f))
 
plot(chm_final, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

```

### Segmentation of the point cloud 

```R
algo <- dalponte2016(chm_final, ttops)
veg_las_segmented <- segment_trees(veg, algo)
 
plot(veg_las_segmented, color ="treeID")
 ```
### Crown Area 

Computing the hulls of each tree 

```R
crowns <- crown_metrics(veg_las_segmented, func = .stdtreemetrics, geom = "convex")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
```
### Getting Individual tree metrics 

```R
f <- function(z) {
  return(list(max_z = max(z)))
}
 
tree_metrics <- crown_metrics(
  las  = veg_las_segmented,
  func = ~f(Z),
  geom = "concave"
)
 
tree_metrics
 
plot(tree_metrics["max_z"])


# Count the amount of trees
n_trees <- nrow(ttops)
print(n_trees)
```
### Compute and map the density of the trees 
```R
# Calcuate the tree density per ha for the AOI

scene_area_m2 <- prod(res(chm_final)[1:2]) * prod(dim(chm_final)[1:2])
scene_area_ha <- scene_area_m2 / 10000
 
n_trees <- nrow(ttops)
 
tree_density <- n_trees / scene_area_ha
print(paste("Tree density:", round(tree_density, 1), "Trees / ha"))

# Calculate and map the density of trees with a 10m resolution
r <- terra::rast(x=ttops)
terra::res(r) <-10
r <-terra::rasterize(x=ttops, y=r, "treeID", fun ='count')
plot (r, col =viridis(20))

```

### Crown metrics and explorative diagrams 
```R
# 1. Calculation of crown area (m²) und and crown radius (m) 
crowns_sf <- st_as_sf(crowns)
crowns_sf$crown_area <- as.numeric(st_area(crowns_sf))
crowns_sf$crown_radius <- sqrt(crowns_sf$crown_area / pi)
 
# 2. Make a 1x3 panel plot
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1)) # 1 Zeile, 3 Spalten
 
# A) Histogram of crown radius
hist(crowns_sf$crown_radius, 
     main = "Histogram of crown radius", 
     xlab = "Crown radius (m)", 
     col = "darkseagreen3", 
     border = "white")
 
# B) Histogramm: Treeheights
hist(ttops$Z, 
     main = "Histogram of tree height", 
     xlab = "Tree height (m)", 
     col = "skyblue3", 
     border = "white")
 
# C) Streudiagramm: Höhe vs. Kronenradius
plot(crowns_sf$crown_radius ~ ttops$Z[1:nrow(crowns_sf)], 
     xlab = "Tree height (m)", 
     ylab = "Crown radius (m)", 
     pch = 16, 
     col = "#00000033")
 


```

## Results

### Tutorials used
- Cannon (n.d.): LiDAR Tutorials: Getting Started with LiDAR Series (https://lab.jonesctr.org/cannon/resources/)
-  Campbell (2024): Working with Lidar Data in R (https://rstudio-pubs-static.s3.amazonaws.com/1218724_8a65c98174fd49ed9289474c4b23fb86.html)
- Goodbody/Roussel (2023): lidR: (A workshop for) Airborne LiDAR Data Manipulation and Visualization for Forestry Application (https://tgoodbody.github.io/lidRtutorial/)
- Roussel et al. (2026): The lidR package: A guide to the lidR package (https://r-lidar.github.io/lidRbook/) 

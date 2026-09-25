# Urban Green Canopy Structure Across the City Districts of Würzburg
Final Assignement for the course Novel Image Analysis Methods/ LIDAR Remote sensing of the M.Sc. Applied Earth Observation and Geoanalysis (Matr.3248301) by Dr. Julia Rieder.

## Introduction 
Urban Green Spaces are key elements of sustainable urban development. Urban green spaces offer benefits in various areas, e.g. ecological benefits in the form of cooling and air purification, health-related benefits such as improved mental health, and economic benefits, for example through higher property values; at the same time, they often serve as meeting places (Zhang/Quian 2024). Urban Green Spaces are also especially important in terms of adapting to the effects of climate change because of the potential of urban heat island mitigation and the enhancement of the resilience of cities to extreme weather events. 

However the effectiveness of urban green depends on the structure of the urban canopy structure. Tree registries such as the [Baumkataster Würzburg](https://opendata.wuerzburg.de/explore/assets/baumkataster_stadt_wuerzburg/) provide an indication of the amount and distribution of trees within the city. Nevertheless, they have a major methodological weakness because only trees on public land are included. Trees located on private land, such as in gardens or allotments, are not included, which results in a significant discrepancy from the actual number of trees. 

In the following study the Urban Green Canopy Structure across the three city districts: Altstadt, Sanderau and Lindleinsmühle of Würzburg are being analyzed. 

The Sanderau (1st), the Altstadt (2nd) and Lindleinsmühle (3rd) are the most populous districts of the city of Würzburg. However, the districts differ primarily in terms of their size and structure.
The Old Town is the oldest district of Würzburg, but it was largely destroyed during the Second World War and quickly rebuilt. The Old Town is characterised primarily by commercial activity, as it is home to the city’s main shopping streets. However, the Old Town area also includes the River Main, the Ringpark and the Hofgarten, as well as the Mainviertel on the left bank of the Main, which also includes the fortress. 
The Sanderau covers an area of ca. 1.62 km². As the most populous district, it is characterised primarily by its role as a residential area. In the northern part, the town is characterised primarily by the ‘Gründerzeit’ architecture of the 19th century, whilst in the south, post-war buildings dominate.
The Lindleinsmühle district is the smallest district in Würzburg (ca. 0.94 km²) and was developed in the 1960s as a result of rapid population growth and to compensate for the destroyed city centre. Consequently, the district is primarily characterised by blocks of flats and high-rise buildings, although areas with detached and terraced houses have also been developed (Würzburg Wiki). 

Given the varying architectural styles across the districts, a difference in urban green space is to be expected. To which extent there is a difference is being analyzed in this study using airborne LiDAR data.


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

The code was executed equally for all districts (The following code block shows an example for Lindleinsmühle).

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

<img width="3884" height="1676" alt="Clipped to city districts" src="https://github.com/user-attachments/assets/74de9d09-11dd-471e-8061-586fdec23528" />


#### Look at the classification 
The LiDAR data are preclassified. The description can be found at the [Landesamt für Digitalisierung, Breitband und Vermessung](https://www.ldbv.bayern.de/mam/ldbv/dateien/laserdaten_punktklassenbeschreibung.pdf). Relevant for this project are the classes 6 -Building Point and 20 - Object Point (e.g. vegetation).

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
# Filter the points the points for vegetatin and buildings
veg      <- filter_poi(nlas, Classification == 20 & Z >= 3 & Z <= 35) 
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
<img width="3878" height="1557" alt="Canopy Height Model" src="https://github.com/user-attachments/assets/f8af1451-81ad-4c95-b977-e1cdf9761749" />



### Individual Tree Detection (ITD) 

The ITD is done using a Local Maximum Filter with variable windows size to adapt to the variant sizes of the trees.
```R
f <- function(x) {x * 0.2 + 5}
heights <- seq(0,30,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

#Locate tree tops (filter for trees taller than 3m)
ttops <- locate_trees(chm_final, lmf(f))
 ttops <- ttops[ttops$Z >3,] #just count trees higher than 3m 
 
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
```

### Segmentation of the point cloud 
Performing a tree segmentation based on the CHM using the Dalponte2016 Algorithm

```R
algo <- dalponte2016(chm_final, ttops)
veg_las_segmented <- segment_trees(veg, algo)
 
plot(veg_las_segmented, color ="treeID")
 ```
### Crown Area 

The crown metrics such as crown area and radius are calculated from concave tree crown polygons 

```R
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
 
# Filter: Every tree with a crown radius below 1m are being removed to eliminate artifacts that got falsely  segmented

crowns_filtered <- crowns_sf[crowns_sf$crown_radius >= 1.0, ]

```


### Tree Count & Density Calculation 


```R
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

```
<img width="1851" height="832" alt="image" src="https://github.com/user-attachments/assets/ae4a2f4f-530f-4c4e-a2cd-62577a936e6e" />


### Height data 
To make precise statements about the vertical structure of the tree canopy, a statistical summary of the (filtered) tree heights was calculated. 


```R
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
```

### Creation of a Explorative Diagrams

For a better understanding of the structural characteristics the urban trees, explorative diagrams were generated including a histogram of the crown radius, the tree heights and a scatterplot of the relationship between tree heights and crown radius was calculated. 

```R
Calculation of crown area (m²) und and crown radius (m) 
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

```


## Results

<img width="1920" height="1080" alt="1" src="https://github.com/user-attachments/assets/bd6cb860-e63a-4ee9-91a8-5504df08c627" />


<img width="1858" height="465" alt="image" src="https://github.com/user-attachments/assets/43a8a8f5-8b54-42c0-989d-5afd0b67d0d4" />





### Tutorials used
- Cannon (n.d.): LiDAR Tutorials: Getting Started with LiDAR Series (https://lab.jonesctr.org/cannon/resources/)
-  Campbell (2024): Working with Lidar Data in R (https://rstudio-pubs-static.s3.amazonaws.com/1218724_8a65c98174fd49ed9289474c4b23fb86.html)
- Goodbody/Roussel (2023): lidR: (A workshop for) Airborne LiDAR Data Manipulation and Visualization for Forestry Application (https://tgoodbody.github.io/lidRtutorial/)
- Roussel et al. (2026): The lidR package: A guide to the lidR package (https://r-lidar.github.io/lidRbook/)

### Literature 
- Würzburgwiki: Verwaltungsgliederunh der Stadt Würzburg (https://wuerzburgwiki.de/wiki/Verwaltungsgliederung_der_Stadt_W%C3%BCrzburg)
-  Zhang/Qian (2024): A comprehensive review of the environmental benefits of urban green spaces (https://doi.org/10.1016/j.envres.2024.118837 )

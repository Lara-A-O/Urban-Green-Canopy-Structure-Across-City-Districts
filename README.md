# Urban Green Canopy Structure Across City Districts
Final Assignement for the course Novel Image Analysis Methods/ LIDAR Remote sensing of the M.Sc. Applied Earth Observation and Geoanalysis(Matr.3248301)

## Introduction 


## Data 
The project used airbone LiDAR data from the Bayerische Vermessungsverwaltung. The data is freely available and can be downloaded at the [Bayerische Vermessungsverwaltung](https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=laserdaten). The data is is downloaded in tiles with a resolution of 1x1km. The tiles covering the areas studies were downloaded accordingly and had previously been cropped in QGIS using the Point Cloud data management tool: Clip. 

## Workflow 

### Getting the project started


```R
library(lidR)
library(terra)  
library(sf)    
library(dplyr)
library(ggplot2)
library(viridis)

```







## Results

### Tutorials used
- Cannon (n.d.): LiDAR Tutorials: Getting Started with LiDAR Series (https://lab.jonesctr.org/cannon/resources/)
-  Campbell (2024): Working with Lidar Data in R (https://rstudio-pubs-static.s3.amazonaws.com/1218724_8a65c98174fd49ed9289474c4b23fb86.html)
- Goodbody/Roussel (2023): lidR: (A workshop for) Airborne LiDAR Data Manipulation and Visualization for Forestry Application (https://tgoodbody.github.io/lidRtutorial/)
- Roussel et al. (2026): The lidR package: A guide to the lidR package (https://r-lidar.github.io/lidRbook/) 

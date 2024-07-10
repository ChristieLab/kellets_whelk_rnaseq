## script written by andy lee on 02/15/2023; contact at andymuanlee@gmail.com
## modified code from: https://blog.devgenius.io/making-publication-quality-inset-maps-in-r-using-ggplot2-85359c492ac8 
## Rversion:

library(sf)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial") # scale bar and north arrow 
library("cowplot") # to use inset maps
library("ggrepel")
library("ggforce")


manual_col <- c("#000004FF", "#D1426FFF","#FEB77EFF")
manual_col4 <- c("#000004FF", "#972C80FF", "#D1426FFF", "#FEB77EFF")

world <- ne_countries(scale = "large", returnclass = "sf")

inset_us <- ggplot(data = world) +
  geom_sf() + 
  coord_sf(xlim=c(-130, -70), ylim=c(53,23), expand=FALSE) + 
  geom_rect(aes(xmin = -123, xmax = -116, ymin = 31, ymax = 37), color = "black", fill = NA) + 
  labs(x = NULL, y = NULL) +
  theme_test() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "lightblue")) 

# Manually input our collection sites 
coordinates <- data.frame(rbind(c("Monterey","-121.8970","36.61817", "Expanded"), c("Naples", "-119.9523", "34.42197", "Historical"), c("Point Loma", "-117.2615", "32.66533", "Historical")))
# c("Diablo Canyon", "-120.8775", "35.22445","Expanded")
colnames(coordinates) <- c("location", "lon", "lat", "range")
coordinates$lon <- as.numeric(coordinates$lon)
coordinates$lat <- as.numeric(coordinates$lat)

coordinates

# create the study map 
# study <- ggplot(data = world) +
#   geom_sf() + 
#   coord_sf(xlim=c(-123.75, -116), ylim=c(37,31), expand=FALSE) + 
#   annotation_scale(location = "br") +
#   annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering) +
#   geom_label_repel(data= coordinates, aes(lon,lat, label=location), nudge_x = -.5, nudge_y =-.5) +
#   geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill=manual_col, shape=23, size=3.5) +
#   theme_test(base_size = 12) + 
#   theme(panel.background = element_rect(fill = "lightblue")) + 
#   xlab("Longitude") + 
#   ylab("Latitude")  + 
#   geom_segment(aes(x = -123.75, y = 34.4486, xend = -120.4716, yend = 34.4486), color = "black", linetype = "dashed") 


# maps joining 
tiff("map.tiff", width=5, height=5, units = "in", res=300)
ggdraw() +
  draw_plot(study) + 
  draw_plot(inset_us, x  = .703, y = .724, width = 0.28, height = 0.25) 
dev.off()

### Create Map with 4 populations
manual_col <- c("#000004FF", "#972C80FF", "#D1426FFF", "#FEB77EFF")

coordinates <- data.frame(rbind(c("Monterey","-121.8970","36.61817", "Expanded"), c("Naples", "-119.9523", "34.42197", "Historical"), c("Point Loma", "-117.2615", "32.66533", "Historical"), c("Diablo Canyon", "-120.8775", "35.22445","Expanded")))
colnames(coordinates) <- c("location", "lon", "lat", "range")
coordinates$lon <- as.numeric(coordinates$lon)
coordinates$lat <- as.numeric(coordinates$lat)
coordinates$loc_id <- c("MON","NAP", "POL", "DIC")

coordinates <- coordinates[order(coordinates$lat, decreasing = TRUE),]
# create the study map 
study_map <- ggplot(data = world) +
  geom_sf() + 
  coord_sf(xlim=c(-123.75, -116), ylim=c(37,31), expand=FALSE) +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering) +
  geom_label_repel(data= coordinates, aes(lon,lat, label=location), nudge_x = -.5, nudge_y =-.5)  +
  geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill=manual_col4, shape=23, size=5) +
  theme_test(base_size = 12) + 
  theme(panel.background = element_rect(fill = "lightblue")) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  geom_segment(aes(x = -123.75, y = 34.4486, xend = -120.4716, yend = 34.4486), color = "black", linetype = "dashed")

# map with inset
tiff("map.tiff", width=5, height=5, units = "in", res=300)
ggdraw() +
  draw_plot(study_map) + 
  draw_plot(inset_us, x = .703, y = .724, width = 0.28, height = 0.25) 
dev.off()


study_map# geom_segment(aes(x = -123.75, y = 34.4486, xend = -120.4716, yend = 34.4486), color = "black", linetype = "dashed") 


##############################################################################################################################################################################  

## Using the ggOceanMaps package to get ocean depth data 

# load packages and set WD
# devtools::install_github("MikkoVihtakari/ggOceanMapsData")
library(ggOceanMapsData)
library(ggOceanMaps)
library(ggrepel)

coordinates <- data.frame(rbind(c("Monterey","-121.8970","36.61817", "Expanded"), c("Naples", "-119.9523", "34.42197", "Historical"), c("Point Loma", "-117.2615", "32.66533", "Historical")))
# c("Diablo Canyon", "-120.8775", "35.22445","Expanded")
colnames(coordinates) <- c("location", "lon", "lat", "range")

coordinates$lon <- as.numeric(coordinates$lon)
coordinates$lat <- as.numeric(coordinates$lat)
# add bathymetry? , bathymetry = TRUE
basemap(limits=c(-123.75, -116, 37, 31))

#   base <- basemap("ArcticStereographic", bathymetry = TRUE, limits=c(-123.75, -116, 37, 31), rotate = TRUE)
#   
#   base + 
#     geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill="orange", shape=23, size=3.5) +
#     labs( x = " Latitude", y = "Longitude") + 
#     annotation_scale(location = "bl") + 
#     annotation_north_arrow(location = "tr", which_north = "true") 

manual_col <- c("#FEB77EFF", "#972C80FF", "#000004FF")

basemap(limits=c(-123.75, -116, 37, 31), grid.col = NA) + 
  geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill=manual_col, shape=23, size=3.5) +
  labs( x = " Latitude", y = "Longitude") + 
  annotation_scale(location = "br") + 
  annotation_north_arrow(location = "bl", which_north = "true") + 
  theme(panel.background = element_rect(fill = "lightblue"), panel.ontop = FALSE) +
  geom_label(data= coordinates, aes(lon,lat, label=location), nudge_x = -1)


# Full Range Map 
full <- basemap("ArcticStereographic", limits=c(-123.75,-111,37,26), bathymetry=TRUE, rotate = TRUE)  
full + 
  geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill="orange", shape=23, size=3.5) + 
  labs( x = " Latitude", y = "Longitude") + 
  # annotation_scale(location = "br") + 
  annotation_north_arrow(location = "tr", which_north = "true") 
  
full_nobath <-  basemap("ArcticStereographic", limits=c(-123.75,-111,37,26), rotate = TRUE)  
full_nobath + 
  geom_spatial_point(data = coordinates, aes(x = lon, y = lat), color = "black", fill="orange", shape=23, size=3.5) + 
  labs( x = " Latitude", y = "Longitude") + 
  # annotation_scale(location = "br") + 
  annotation_north_arrow(location = "tr", which_north = "true") 



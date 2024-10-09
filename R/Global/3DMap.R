
install.packages("RcppArmadillo", type = "binary")
install.packages("rayimage")
install.packages("rayvertex")
install.packages("rayrender", type = "binary")
install.packages("rayshader", type = "binary")
install.packages("rgl", type = "binary")
install.packages("ggplot2")
install.packages("sf")
install.packages("maps")

# Load libraries
library(ggplot2)
library(sf)
library(rayshader)
library(maps)
library(dplyr)
library(tidyverse) 

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Extract sample site locations -------------------------------------------
# Average tPCB per sample site
tpcb.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = wdc, mean)

# Define the latitude and longitude boundaries for the Great Lakes region
lat_min <- 41.5
lat_max <- 47.5
long_min <- -88
long_max <- -78

# Filter the tpcb.ave data for the Great Lakes region
great_lakes_data <- tpcb.ave %>%
  filter(Latitude >= lat_min, Latitude <= lat_max, 
         Longitude >= long_min, Longitude <= long_max)

# Check for duplicates and summarize
duplicates <- great_lakes_data %>%
  group_by(Latitude, Longitude) %>%
  summarise(tPCB = mean(tPCB, na.rm = TRUE), .groups = 'drop')

# Create a grid for rendering
tPCB_matrix <- duplicates %>%
  select(Longitude, Latitude, tPCB) %>%
  pivot_wider(names_from = Longitude, values_from = tPCB) %>%
  as.data.frame()

# Print the filtered data dimensions to ensure it has been processed correctly
print(dim(tPCB_matrix))

# Check for missing values in the matrix
missing_values <- any(is.na(tPCB_matrix))
if (missing_values) {
  warning("There are missing values in the matrix, please check the data.")
}

# Create the matrix from the filtered data, ensuring it is properly filled
reshaped_tPCB_matrix <- as.matrix(tPCB_matrix[,-1])  # Remove Latitude for the matrix creation

# Now ensure that reshaped_tPCB_matrix has the correct dimensions
print(dim(reshaped_tPCB_matrix))

# Create the 3D plot for the Great Lakes region
reshaped_tPCB_matrix %>%
  sphere_shade(texture = "imhof1") %>%
  plot_3d(reshaped_tPCB_matrix, zscale = 0.1, fov = 0, theta = 45, phi = 30, 
          windowsize = c(800, 800), zoom = 0.75)

# Render depth
render_depth()


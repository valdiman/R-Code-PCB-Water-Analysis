# Install packages
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
{
  library(ggplot2)
  library(rayshader)
  library(viridis)
  library(dplyr)
  library(maps)
  library(rgl)
}

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Extract sample site locations -------------------------------------------
# Average tPCB per sample site
tpcb.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = wdc, mean)

# Define the latitude and longitude boundaries for the entire USA
lat_min <- 24.396308
lat_max <- 49.384358
long_min <- -125.0
long_max <- -66.93457

# Filter the tpcb.ave data for the entire USA
usa_data <- tpcb.ave %>%
  filter(Latitude >= lat_min, Latitude <= lat_max, 
         Longitude >= long_min, Longitude <= long_max) %>%
  mutate(tPCB_log10 = log10(tPCB))  # Apply log10 transformation to tPCB

# Get map data for the USA and states
usa_map <- map_data("usa")
states_map <- map_data("state")  # Get state boundaries
lakes_map <- map_data("lakes")

gg_map_usa <- ggplot() +
  # Fill USA land with a light gray color
  geom_polygon(data = usa_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = NA) +  
  # Add state lines using the state boundaries
  geom_path(data = states_map, aes(x = long, y = lat, group = group), 
            color = "black", size = 0.1) +  # Draw state lines
  # Fill lakes with a light blue color
  geom_polygon(data = lakes_map, aes(x = long, y = lat, group = group), 
               fill = "lightblue", color = NA, alpha = 0.5) +  # Set lakes to be a bit transparent
  # Plot log10-transformed data points, with transparency added using alpha
  geom_point(data = usa_data, aes(x = Longitude, y = Latitude, color = tPCB_log10), 
             size = 2, alpha = 0.7, show.legend = FALSE) +  # Remove legend
  scale_color_viridis_c(option = "D", name = "log10(tPCB)", direction = -1) +  # Reverse color scale
  coord_map(xlim = c(long_min, long_max), ylim = c(lat_min, lat_max)) +  # Limit plot to specific latitude/longitude bounds
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +  # Minimal theme with clean background
  theme(panel.background = element_rect(fill = "white"),  # Set the panel background to white
        plot.background = element_rect(fill = "white"))  

# Plot the ggplot object for the USA
print(gg_map_usa)

# Convert the ggplot into a 3D plot using rayshader without the legend
plot_gg(gg_map_usa, width = 5, height = 4, scale = 300, multicore = TRUE, windowsize = c(1000, 800))

# Adjust the camera settings
render_camera(fov = 70, zoom = 0.5, theta = 130, phi = 35)

# Render the snapshot without the legend
Sys.sleep(0.2)
render_snapshot(clear = TRUE)  # No keep_user_par argument

# Save the 3D plot as an image file
rgl::snapshot3d("Output/Maps/Global/USA_tPCB_3D_map.png")


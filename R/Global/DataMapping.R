## Water PCB concentrations mapping.

# Install packages
install.packages("ggplot2")
install.packages("devtools")
install.packages("dplyr")
install.packages("stringr")
install.packages("maps")
install.packages("mapdata")
install.packages("ggmap")
install.packages("usethis")
install.packages("GISTools")
install.packages("rgeos")
install.packages("ggsn")
install.packages("ggrepel")
install.packages("ggpp")
install.packages("scales")
install.packages("viridis")
install.packages("tidyr")

# Load libraries
{
  library(dplyr)
  library(usethis)
  library(devtools)
  library(ggplot2)
  library(ggmap) # function map_data
  library(maps)
  library(leaflet)
  library(rgeos)
  library(ggsn)
  library(ggrepel)
  library(reshape2)
  library(ggpmisc)
  library(scales) # add commas in legend in maps
  library(cowplot)
  library(viridis) # customize color legend
  library(tidyr)
}

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Extract sample site locations -------------------------------------------
# Average tPCB per sample site
tpcb.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = wdc, mean)

# Find number of samples per StateSampled
wdc.2 <- wdc %>%
  group_by(StateSampled) %>%
  summarise(n = n())

# Convert wdc.2 to long format
wdc.3 <- gather(wdc.2, key = "Variable", value = "Value", -StateSampled)

# Find number of samples per LocationName
wdc.4 <- wdc %>%
  group_by(LocationName) %>%
  summarise(n = n())

# Convert wdc.4 to long format
wdc.5 <- gather(wdc.4, key = "Variable", value = "Value", -LocationName)

# USA/State maps -------------------------------------------------------------
us <- map_data("usa")
states <- map_data("state")
counties <- map_data("county")

# (1) Map of US with locations
maploc <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "transparent") +
  coord_fixed(1.3, xlim = c(min(us$long) - 10, max(us$long) + 2)) +
  theme_void() +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = tpcb.ave, aes(x = Longitude, y = Latitude),
             color = "black",
             size = 2.2, shape = 20) +
  geom_text(data = wdc.3, aes(x = -66.9, y = 50.2 - seq_along(StateSampled),
                              label = paste(StateSampled, Value), hjust = 0,
                              vjust = 1), size = 3) +
  geom_text(data = wdc.5, aes(x = -138, y = 50.2 - seq_along(LocationName),
                              label = paste(LocationName, Value), hjust = 0,
                              vjust = 1), size = 3) +
  geom_text(aes(x = -69.8, y = 50.5, label = "States/Samples", hjust = 0,
                vjust = 1), size = 3, fontface = "bold") +
  geom_text(aes(x = -138, y = 50.5, label = "Location/Samples", hjust = 0,
                vjust = 1), size = 3, fontface = "bold")

print(maploc)

# Save map in folder
ggsave("Output/Maps/Global/maploc.pdf", plot = maploc,
       width = 14, height = 6, dpi = 300)

# (2) Map + tPCB
maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = tpcb.ave, aes(x = Longitude, y = Latitude,
                                  fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(1, 3500000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(maptPCB)  # Print the plot

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCB.pdf", plot = maptPCB,
       width = 14, height = 4)

# Individual PCB Maps -----------------------------------------------------
# PCB11
pcb11 <- wdc[wdc$PCB11 != 0, ]

# Average PCB11 per sample site
pcb.11.ave <- aggregate(PCB11 ~ SiteID + Latitude + Longitude,
                      data = pcb11, mean)

# Plot
mapPCB11 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.11.ave, aes(x = Longitude, y = Latitude,
                                  fill = PCB11), alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCB 11")),
    limits = c(0.1, 2000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(1.1, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )
        
print(mapPCB11)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB11.pdf", plot = mapPCB11,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB 20+21+28+31+33+50+53
pcb20 <- wdc[wdc$PCB20.21.28.31.33.50.53 != 0, ]

# Average PCB11 per sample site
pcb.20.ave <- aggregate(PCB20.21.28.31.33.50.53 ~ SiteID + Latitude + Longitude,
                      data = pcb20, mean)

# Plot
mapPCB20 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.20.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB20.21.28.31.33.50.53),
             alpha = 1, color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(name = expression(bold(atop("PCBs 20+21+28",
                                                   paste("+31+33+50+53")))),
    limits = c(1, 1000000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(1.15, 0.5),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB20)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB20.pdf", plot = mapPCB20,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB44+47+65
pcb44 <- wdc[wdc$PCB44.47.65 != 0, ]

# Average PCB11 per sample site
pcb.44.ave <- aggregate(PCB44.47.65 ~ SiteID + Latitude + Longitude,
                      data = pcb44, mean)

# Plot
mapPCB44 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.44.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB44.47.65), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCBs 44+47+65")),
    limits = c(0.5, 200000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(1.12, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB44)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB44.pdf", plot = mapPCB44,
       width = 14, height = 4)

# Filter out rows with 0 values for PCB67
pcb67 <- wdc[wdc$PCB67 != 0, ]

# Average PCB67 per sample site
pcb.67.ave <- aggregate(PCB67 ~ SiteID + Latitude + Longitude,
                      data = pcb67, mean)

# Plot
mapPCB67 <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = pcb.67.ave, aes(x = Longitude, y = Latitude,
                                    fill = PCB67), alpha = 1,
             color  = "black",
             shape = 21, size = 2, stroke = 0.75) +
  scale_fill_viridis_c(
    name = expression(bold("PCB 67")),
    limits = c(1, 500),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(1.05, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18)  # Adjust the size of legend title
  )

print(mapPCB67)  # Print the plot

# Save map in folder
ggsave("Output/Maps/Global/mapPCB67.pdf", plot = mapPCB67,
       width = 14, height = 4)

# Specific locations ------------------------------------------------------
# Version 1
# Midwest
mid.west <- subset(wdc, LocationName %in% c("DEQ (Michigan)", "Fox River", "Kalamazoo River",
                                            "Lake Michigan Mass Balance",
                                            "Indiana Harbor and Ship Canal",
                                            "USGS (Wisconsin)"))

mid.west.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = mid.west, mean)

maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = mid.west.ave, aes(x = Longitude, y = Latitude,
                                      fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 3, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(20, 730000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3, xlim = c(-89.1, -82.4), ylim = c(41.4, 47.2)) +  # Adjust these values accordingly
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position.inside = c(0.5, 0.55),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                              margin = margin(b = 10))  # Adjust title position
  ) +
  ggtitle("Lake Michigan and nearby locations") +
  geom_text(aes(x = -87, y = 41.3, label = "IHSC"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -84.9, y = 42.6, label = "KR"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -88.7, y = 44.7, label = "FR"), size = 6,
            color = "black", fontface = "bold")

print(maptPCB)

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBMidWest.pdf", plot = maptPCB,
       width = 14, height = 4)

# North east coast
nec <- subset(wdc, LocationName %in% c("Housatonic River", "Hudson River",
                                       "Passaic River", "Chesapeake Bay",
                                       "New Bedford Harbor", "Newtown Creek",
                                       "Richardson Hill Road Landfill",
                                       "Anacostia River", "Susquehanna River"))

nec.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                     data = nec, mean)

maptPCB <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_polygon(data = counties, aes(x = long, y = lat, group = group),
               color = "grey", fill = NA) +  # County boundaries
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "black") +
  geom_point(data = nec.ave, aes(x = Longitude, y = Latitude,
                                 fill = tPCB), alpha = 1, color  = "black",
             shape = 21, size = 3, stroke = 0.75) +
  scale_fill_viridis_c(
    name = element_blank(),
    limits = c(10, 3500000),
    trans = "log10",
    labels = scales::comma,
    begin = 1,  # Adjust the starting color (lower value)
    end = 0.001     # Adjust the ending color (higher value)
  ) +
  coord_fixed(1.3, xlim = c(-78, -70.25), ylim = c(37.5, 43.2)) +  # Adjust these values accordingly
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(3, "lines"),
    legend.position = c(1.3, 0.6),  # Adjust the legend position
    legend.text = element_text(size = 18),  # Adjust the size of legend labels
    legend.title = element_text(size = 18),
    legend.background = element_rect(color = "white"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                              margin = margin(b = 10))  # Adjust title position
  ) +
  ggtitle("Northeast Coast Locations") +
  geom_text(aes(x = -70.9, y = 42, label = "NBH"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -72.5, y = 42.6, label = "HoR"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -74.4, y = 43.2, label = "HuR"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -75.8, y = 42.6, label = "RHRL"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -74.5, y = 41.2, label = "PR"), size = 6,
            color = "black", fontface = "bold") +
  geom_text(aes(x = -77.4, y = 40.3, label = "AR"), size = 6,
            color = "black", fontface = "bold") +
  geom_segment(aes(x = -77.4, y = 40.1, xend = -76.9, yend = 39.1),
               arrow = arrow(length = unit(0.3, "cm")), color = "black",
               linewidth = 1.5) +
  geom_text(aes(x = -73.2, y = 40.2, label = "NC"), size = 6,
            color = "black", fontface = "bold") +
  geom_segment(aes(x = -73.2, y = 40.4, xend = -73.7, yend = 40.7),
               arrow = arrow(length = unit(0.3, "cm")), color = "black",
               linewidth = 1.5) +
  geom_text(aes(x = -76.2, y = 40.4, label = "SR"), size = 6,
            color = "black", fontface = "bold") +
  geom_segment(aes(x = -76.2, y = 40.2, xend = -76.2, yend = 39.8),
               arrow = arrow(length = unit(0.3, "cm")), color = "black",
               linewidth = 1.5) +
  geom_text(aes(x = -74.2, y = 38.7, label = "CB"), size = 6,
            color = "black", fontface = "bold") +
  geom_segment(aes(x = -74.7, y = 38.7, xend = -76.4, yend = 38.7),
               arrow = arrow(length = unit(0.3, "cm")), color = "black",
               linewidth = 1.5)

print(maptPCB)

# Save the plot as PDF
ggsave("Output/Maps/Global/maptPCBNorthEastCost.pdf", plot = maptPCB,
       width = 14, height = 4)

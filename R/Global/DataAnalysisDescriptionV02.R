## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of the sum of PCB.

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("zoo")
install.packages("reshape")
install.packages("sf")
install.packages("sfheaders")
install.packages("viridis")
install.packages("lubridate")
install.packages("cowplot")

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(zoo) # yields seasons
  library(reshape)
  library(sf)
  library(sfheaders) # Create file to be used in Google Earth
  library(viridis) # colors for plots
  library(lubridate)
  library(cowplot) # to create figure with many plots
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/Water_Congener_Aroclor.csv")

# General information -----------------------------------------------------
# Number of locations and number of site per location
location_count <- wdc %>%
  group_by(Location..Location.Name.) %>%
  summarise(count = n())

print(location_count)
# Median amount of samples per location
median(location_count$count)

# Number of locations per states
state_count <- wdc %>%
  group_by(Location..State.Sampled.) %>%
  summarise(count = n())

print(state_count)

# Number of sites and number of replicates per site
site_count <- wdc %>%
  group_by(Sample.ID) %>%
  summarise(count = n())

print(site_count)

# Media of number of site available
median(site_count$count)

# Number of replicates per sites from the same day
site_repli_count <- wdc %>%
  group_by(Sample.ID) %>%
  summarise(count = n())

print(site_repli_count)

# Media of number of replicates per sites from the same day
median(site_repli_count$count)

# Find the Sample.ID with the highest count
max_count_sample <- site_repli_count %>%
  filter(count == max(count)) %>%
  pull(Sample.ID)

filtered_wdc <- wdc %>%
  filter(Sample.ID == max_count_sample)

# Extract Site Name
extracted_string <- filtered_wdc %>%
  select(Site..Site.Name.) %>%
  pull()

# Display Site Name
print(extracted_string[1])

# Aroclor Congener Summary ------------------------------------------------
# Create a new data frame with 'Year' and 'AroclorCongener' columns
wdc.aroclorcongener <- wdc %>%
  mutate(Year = as.integer(format(as.Date(Date.Time..Sample.Date., format = "%Y-%m-%d"), "%Y"))) %>%
  filter(Comment..Aroclor.Congener. %in% c("Congener", "Aroclor")) %>%
  select(Year, Comment..Aroclor.Congener.)

# Calculate counts for 'Congener' and 'Aroclor' for each year
yearly_counts <- wdc.aroclorcongener %>%
  group_by(Year, Comment..Aroclor.Congener.) %>%
  summarize(Count = n())

# Calculate percentages for 'Congener' and 'Aroclor' for each year
percentages <- yearly_counts %>%
  group_by(Year) %>%
  mutate(Percent = Count / sum(Count))

# Precompute the percentages for 'Congener' and 'Aroclor' for each year
precomputed_percentages <- percentages %>%
  group_by(Year, Comment..Aroclor.Congener.) %>%
  summarize(Percent = sum(Percent))

# Calculate the number of samples per year for 'Congener' and 'Aroclor'
sample_counts <- yearly_counts %>%
  group_by(Year) %>%
  summarize(TotalSamples = sum(Count))

# Convert Year to character
precomputed_percentages$Year <- as.character(precomputed_percentages$Year)

# Create a stacked bar plot with percentages as percentages (0 to 100)
plot.aroclor.congener <- ggplot(precomputed_percentages,
                                aes(x = factor(Year,
                                               levels = unique(precomputed_percentages$Year)))) +
  geom_bar(aes(y = Percent * 100, fill = AroclorCongener), stat = "identity") +
  geom_text(data = sample_counts, aes(label = TotalSamples, y = 100), size = 2.5,
            fontface = "bold") +
  theme_classic() +
  labs(title = "Percentage of Congener and Aroclor Over the Years",
       x = "Year",
       y = "Percentage") +
  scale_fill_manual(values = c("Congener" = "grey50", "Aroclor" = "grey90")) +
  scale_y_continuous(limits = c(0, 100)) +  # Set the y-axis limits
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11)) +
  coord_cartesian(ylim = c(0, 100))

# See the plot
print(plot.aroclor.congener)

# Save plot in folder
ggsave("Output/Plots/Global/AroclorCongener.png", plot = plot.aroclor.congener,
       width = 10, height = 3, dpi = 300)

# Aroclor summary ---------------------------------------------------------
# Number of samples analyzed using Aroclor method
count_Aroclor <- sum(wdc$Comment..Aroclor.Congener. == "Aroclor")
total_samples <- length(wdc[,1])
percent_aroclor <- count_Aroclor/total_samples*100

# Calculate sample % for each Aroclor mixtures
aroclors <- c('PCB.water..pg.l...A1016.', 'PCB.water..pg.l...A1221.',
              'PCB.water..pg.l...A1232.', 'PCB.water..pg.l...A1242.',
              'PCB.water..pg.l...A1248.', 'PCB.water..pg.l...A1254.',
              'PCB.water..pg.l...A1260.')

# Calculate the number of non-NA values (frequency of numbers) in each Aroclor
frequency_aroclors <- lapply(wdc[aroclors], function(column) {
  length(na.omit(column))
})

# Percentage of Aroclor mixtures in relation to all Aroclors
for (i in seq_along(aroclors)) {
  column_name <- aroclors[i]
  print(paste(column_name))
  print(frequency_aroclors[[i]]/count_Aroclor*100)
}

# Congener frequency ------------------------------------------------------
{
  # Remove metadata
  wdc.nmd <- subset(wdc, select = -c(Reference:Comment..Aroclor.Congener.)) #nmd = no meta data
  # Remove Aroclor data
  wdc.nmd <- subset(wdc.nmd, select = -c(PCB.water..pg.l...A1016.:PCB.water..pg.l...A1260.))
  # (2) Only consider congener data
  wdc.cong <- subset(wdc, Comment..Aroclor.Congener. == "Congener")
  # Remove metadata
  wdc.cong.1 <- subset(wdc.cong, select = -c(Reference:Comment..Aroclor.Congener.))
  # Remove Aroclor data and tPCB
  wdc.cong.1 <- subset(wdc.cong.1, select = -c(PCB.water..pg.l...A1016.:PCB.water..pg.l...tPCB.))
}

# Create a frequency detection plot
{
  wdc.cong.freq <- colSums(! is.na(wdc.cong.1) & (wdc.cong.1 !=0))/nrow(wdc.cong.1)
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  colnames(wdc.cong.freq) <- c("PCB.frequency")
  congener <- row.names(wdc.cong.freq)
  wdc.cong.freq <- cbind(congener, wdc.cong.freq$PCB.frequency)
  colnames(wdc.cong.freq) <- c("congener", "PCB.frequency")
  wdc.cong.freq <- data.frame(wdc.cong.freq)
  wdc.cong.freq$congener <- as.character(wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('PCB\\.water\\.\\.pg\\.l\\.\\.\\.', '', wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('(PCB[0-9]+)\\.([0-9]+)', '\\1+\\2', wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('\\.$', '', wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('\\+$', '', wdc.cong.freq$congener)
  wdc.cong.freq$congener <- gsub('\\.(?=[0-9])', '+', wdc.cong.freq$congener, perl = TRUE)
  wdc.cong.freq$congener <- gsub('\\.$', '', wdc.cong.freq$congener)
  # Convert PCB.frequency to numeric
  wdc.cong.freq$PCB.frequency <- as.numeric(as.character(wdc.cong.freq$PCB.frequency))
  # Change the order of the factor levels
  wdc.cong.freq$congener <- factor(wdc.cong.freq$congener,
                                   levels = rev(wdc.cong.freq$congener))
}

# Summary statistic of frequency of detection
summary(wdc.cong.freq$PCB.frequency)

# Frequency detection plot
plot.cong.freq <- ggplot(wdc.cong.freq, aes(x = 100*PCB.frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "white", color = "black") +
  geom_vline(xintercept = 100*mean(wdc.cong.freq$PCB.frequency),
             color = "red") +
  ylab("") +
  theme_bw() +
  xlim(c(0,100)) +
  theme(aspect.ratio = 20/5) +
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.text.y = element_text(face = "bold", size = 7))

# See plot
print(plot.cong.freq)

# Save map in folder
ggsave("Output/Plots/Global/FreqPCB.png", plot = plot.cong.freq,
       width = 5, height = 10, dpi = 300)

# Total PCB description ---------------------------------------------------
# Data preparation
{
  # Change date format
  SampleDate <- as.Date(wdc$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb <- cbind(factor(wdc$SiteID), SampleDate,
                wdc$Latitude, wdc$Longitude, wdc$tPCB,
                data.frame(time.day), season.s)
  # Add column names
  colnames(tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                      "tPCB", "time", "season")
}

# Get coordinates per site to plot in Google Earth
location <- tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                      data = location, mean)
# Create an sf data frame
sf_location <- st_as_sf(location, coords = c("Longitude", "Latitude"))
# Set the CRS to WGS 84 (EPSG:4326)
sf_location <- st_set_crs(sf_location, 4326)
# Define the full file path for the KML file
kmlFilePath <- "Output/Data/Global/PCBSampleLocations.kml"
# Write the KML file to the specified directory
st_write(sf_location, kmlFilePath, driver = "kml", append = FALSE)

# Summary statistic of total PCB (congeners + Aroclor) in pg/L not including 0s
summary(wdc$tPCB)

# Highest sample location and time
max_sample <- wdc[which.max(wdc$tPCB), ]
max_sample <- max_sample[c("LocationName", "SampleDate", "SiteName")]
print(max_sample)

# Individual congeners description
summary(wdc.cong.1, na.rm = T, zero = T)
# Get the max value for each congener
cong.max <-as.numeric(sub('.*:', '',
                          summary(wdc.cong.1, na.rm = T,
                                  zero = T)[6,]))
# Add congener
cong.max <- cbind(congener, data.frame(cong.max))

# Obtain the median for each individual congener
cong.median <- as.numeric(sub('.*:',
                              '', summary(wdc.cong.1, na.rm = T,
                                          zero = T)[3,]))
# Add congener
cong.median <- cbind(congener, data.frame(cong.median))
# Min
print(min(cong.median$cong.median))
#Max
print(max(cong.median$cong.median))
# Mean
print(mean(cong.median$cong.median))

# Global plots ------------------------------------------------------------
# (1) Histogram
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))

## (2) Total PCBs in 1 box plot
## include 64 and 640 pg/L from EPA
## Select data
selected_data_1 <- wdc[, c("SampleDate", "tPCB")]
# Convert SampleDate to year
selected_data_1$Year <- lubridate::year(as.Date(selected_data_1$SampleDate,
                                                format = "%Y-%m-%d"))

# Group years into intervals of 7 years
selected_data_1$Year_Group <- cut(selected_data_1$Year,
                                  breaks = seq(1979, 2021, by = 7),
                                  include.lowest = TRUE)

# Create ggplot2 plot with viridis color palette
tpcb_points_ggplot_1 <- ggplot(selected_data_1, aes(x = "", y = tPCB,
                                                    color = factor(Year_Group),
                                                    fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3),
             size = 1.5, alpha = 0.6) +
  geom_boxplot(lwd = 0.8, width = 0.7, outlier.shape = NA, alpha = 0.7) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.5, 10^7.5)) +
  xlab(expression(bold("(n = 5132)"))) +
  labs(y = expression(bold("Water Concentration " *Sigma*"PCB (pg/L)")),
       shape = "PCB Method") +
  theme_bw() +
  theme(aspect.ratio = 5/1,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11,
                                    angle = 60, hjust = 0, vjust = 0.3),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1979-1985", "1986-1992", "1993-1999",
                                 "2000-2006", "2007-2013", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1979-1985", "1986-1992", "1993-1999",
                                "2000-2006", "2007-2013", "2014-2020")) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 640, color = "#9999CC",
             linewidth = 0.8) +
  geom_hline(yintercept = 64, color = "#CC6666",
             linewidth = 0.8)

tpcb_points_ggplot_1  # Display the plot

# Save plot in folder
ggsave("Output/Plots/Global/tPCBBoxPlotV02.jpg", plot = tpcb_points_ggplot_1,
       width = 10, height = 4.38, dpi = 500)

# Calculate % samples above both EPA thresholds
EPA640 <- sum(wdc$tPCB > 640)/nrow(wdc)*100
EPA64 <- sum(wdc$tPCB > 64)/nrow(wdc)*100

# (2.1) Total PCBs for selected locations
# Selecte data
selected_data_2 <- wdc[, c("LocationName", "SampleDate", "tPCB")]

sites_to_include <- c("Anacostia River", "Housatonic River",
                      "New Bedford Harbor", "Passaic River",
                      "Lake Michigan Mass Balance", "Hudson River",
                      "Kalamazoo River", "Fox River", "Portland Harbor",
                      "Lake Michigan Mass Balance", "Spokane River",
                      "Chesapeake Bay", "Bannister Fed Complex", "DEQ (Michigan)")

# Filter the data to include only the specified sites
filtered_data_2 <- selected_data_2 %>%
  filter(LocationName %in% sites_to_include)

# Convert SampleDate to year
filtered_data_2$Year <- lubridate::year(as.Date(filtered_data_2$SampleDate,
                                              format = "%Y-%m-%d"))

# Group years into intervals of 7 years
filtered_data_2$Year_Group <- cut(filtered_data_2$Year,
                                breaks = seq(1979, 2021, by = 7),
                                include.lowest = TRUE)

# Calculate the number of samples in each LocationName group
sample_counts <- table(filtered_data_2$LocationName)

# Create ggplot2 plot with viridis color palette
tpcb_points_ggplot_2 <- ggplot(filtered_data_2, aes(x = LocationName, y = tPCB,
                                                color = factor(Year_Group),
                                                fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.5, 10^7.5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration " *Sigma*"PCB (pg/L)")),
       shape = "PCB Method") +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1979-1985", "1986-1992", "1993-1999",
                                 "2000-2006", "2007-2013", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1979-1985", "1986-1992", "1993-1999",
                                "2000-2006", "2007-2013", "2014-2020")) +
  annotation_logticks(sides = "l")

# Add the number of samples to the x-axis labels
tpcb_points_ggplot_2 <- tpcb_points_ggplot_2 + 
  scale_x_discrete(labels = function(x) paste(x, "(n=", sample_counts[x], ")"))

tpcb_points_ggplot_2  # Display the plot

# Save plot in folder
ggsave("Output/Plots/Global/tPCBSiteV03.jpg", plot = tpcb_points_ggplot_2,
       width = 10, height = 6, dpi = 500)

# (3) Box plot for individual PCBs
PCBi_boxplot <- ggplot(stack(wdc.cong.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, shape = 21, outlier.fill = "white",
               fill = "white", outlier.shape = 21) +
  scale_x_discrete(labels = wdc.cong.freq$congener) +
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# See plot
print(PCBi_boxplot)

# Save map in folder
ggsave("Output/Plots/Global/PCBiBoxPlot.png", plot = PCBi_boxplot,
       width = 10, height = 5, dpi = 300)

# Individual PCBs per locations -------------------------------------------
# Filter out rows with NA and 0 values in the 'PCBi' column
# (5.1) PCB5+8 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb5 <- wdc %>%
  filter(!is.na(PCB5.8), !(PCB5.8 == 0)) %>%
  select(LocationName, SampleDate, PCB5.8) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb5_points_ggplot <- ggplot(filtered_data_pcb5, aes(x = LocationName, y = PCB5.8,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^7)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 5+8 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007", "2007-2014", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007", "2007-2014", "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb5_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB5Site.jpg", plot = pcb5_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.2) PCB11 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb11 <- wdc %>%
  filter(!is.na(PCB11), !(PCB11 == 0)) %>%
  select(LocationName, SampleDate, PCB11) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb11_points_ggplot <- ggplot(filtered_data_pcb11, aes(x = LocationName, y = PCB11,
                                                    color = factor(Year_Group),
                                                    fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^4)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCB 11 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("2000-2007", "2007-2014", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("2000-2007", "2007-2014", "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb11_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB11Site.jpg", plot = pcb11_points_ggplot,
       width = 6, height = 10, dpi = 500)

# (5.3) PCB20.21.28.31.33.50.53 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb20 <- wdc %>%
  filter(!is.na(PCB20.21.28.31.33.50.53),
         !(PCB20.21.28.31.33.50.53 == 0)) %>%
  select(LocationName, SampleDate, PCB20.21.28.31.33.50.53) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb20_points_ggplot <- ggplot(filtered_data_pcb20, aes(x = LocationName,
                                                       y = PCB20.21.28.31.33.50.53,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.1, 10^6.2)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 20+21+28+31+33+50+53 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007","2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007","2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb20_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB20Site.jpg", plot = pcb20_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.4) PCB35
# Filter the data to remove NAs and 0s
filtered_data_pcb35 <- wdc %>%
  filter(!is.na(PCB35), !(PCB35 == 0)) %>%
  select(LocationName, SampleDate, PCB35) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7), include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb35_points_ggplot <- ggplot(filtered_data_pcb35, aes(x = LocationName, y = PCB35,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^3.5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCB 35 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("2000-2007", "2007-2014", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("2000-2007", "2007-2014", "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb35_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB35Site.jpg", plot = pcb35_points_ggplot,
       width = 6, height = 10, dpi = 500)

# (5.5) PCB43+49+52+69+73
# Filter the data to remove NAs and 0s
filtered_data_pcb43 <- wdc %>%
  filter(!is.na(PCB43.49.52.69.73), !(PCB43.49.52.69.73 == 0)) %>%
  select(LocationName, SampleDate, PCB43.49.52.69.73) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb43_points_ggplot <- ggplot(filtered_data_pcb43, aes(x = LocationName, y = PCB43.49.52.69.73,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^6.5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 43+49+52+69+73 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007", "2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007", "2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb43_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB43Site.jpg", plot = pcb43_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.6) PCBB44+47+65 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb44 <- wdc %>%
  filter(!is.na(PCB44.47.65), !(PCB44.47.65 == 0)) %>%
  select(LocationName, SampleDate, PCB44.47.65) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb44_points_ggplot <- ggplot(filtered_data_pcb44, aes(x = LocationName,
                                                       y = PCB44.47.65,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.1, 10^6.2)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 44+47+65 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007","2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007","2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb44_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB44Site.jpg", plot = pcb44_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.7) PCBB45.51 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb45 <- wdc %>%
  filter(!is.na(PCB45.51), !(PCB45.51 == 0)) %>%
  select(LocationName, SampleDate, PCB45.51) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb45_points_ggplot <- ggplot(filtered_data_pcb45, aes(x = LocationName,
                                                       y = PCB45.51,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^4.5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 45+51 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007","2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007","2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb45_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB45Site.jpg", plot = pcb45_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.8) PCBB67 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb67 <- wdc %>%
  filter(!is.na(PCB67), !(PCB67 == 0)) %>%
  select(LocationName, SampleDate, PCB67) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb67_points_ggplot <- ggplot(filtered_data_pcb67, aes(x = LocationName,
                                                       y = PCB67,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^4)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCB 67 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("2000-2007","2007-2014", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("2000-2007","2007-2014", "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb67_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB67Site.jpg", plot = pcb67_points_ggplot,
       width = 6, height = 10, dpi = 500)

# (5.9) PCBB68 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb68 <- wdc %>%
  filter(!is.na(PCB68), !(PCB68 == 0)) %>%
  select(LocationName, SampleDate, PCB68) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb68_points_ggplot <- ggplot(filtered_data_pcb68, aes(x = LocationName,
                                                       y = PCB68,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^3.5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCB 68 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("2000-2007","2007-2014", "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("2000-2007","2007-2014", "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb68_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB68Site.jpg", plot = pcb68_points_ggplot,
       width = 6, height = 10, dpi = 500)

# (5.10) PCB 77+85+110+111+115+116+117 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb77 <- wdc %>%
  filter(!is.na(PCB77.85.110.111.115.116.117),
         !(PCB77.85.110.111.115.116.117 == 0)) %>%
  select(LocationName, SampleDate, PCB77.85.110.111.115.116.117) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb77_points_ggplot <- ggplot(filtered_data_pcb77, aes(x = LocationName,
                                                       y = PCB77.85.110.111.115.116.117,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^6)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCBs 77+85+110+111+115+116+117 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007","2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007","2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb77_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB77Site.jpg", plot = pcb77_points_ggplot,
       width = 8, height = 10, dpi = 500)

# (5.11) PCB209 plot
# Filter the data to remove NAs and 0s
filtered_data_pcb209 <- wdc %>%
  filter(!is.na(PCB209), !(PCB209 == 0)) %>%
  select(LocationName, SampleDate, PCB209) %>%
  mutate(Year = year(as.Date(SampleDate, format = "%Y-%m-%d")),
         Year_Group = cut(Year, breaks = seq(1979, 2021, by = 7),
                          include.lowest = TRUE))

# Create ggplot2 plot with viridis color palette
pcb209_points_ggplot <- ggplot(filtered_data_pcb209, aes(x = LocationName, y = PCB209,
                                                       color = factor(Year_Group),
                                                       fill = factor(Year_Group))) +
  geom_point(shape = 21, color = "black", position = position_jitter(0.3), size = 1.5) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.01, 10^5)) +
  labs(x = NULL,
       y = expression(bold("Water Concentration PCB 209 (pg/L)"))) +
  theme_bw() +
  theme(aspect.ratio = 10/15,
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.text.x = element_text(face = "bold", size = 11,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8),
        axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  scale_color_viridis(discrete = TRUE, name = "Year Group",
                      labels = c("1993-2000", "2000-2007", "2007-2014",
                                 "2014-2020")) +
  scale_fill_viridis(discrete = TRUE, name = "Year Group",
                     labels = c("1993-2000", "2000-2007", "2007-2014",
                                "2014-2020")) +
  annotation_logticks(sides = "l")

# Display the plot
pcb209_points_ggplot

# Save plot in folder
ggsave("Output/Plots/Global/PCB209Site.jpg", plot = pcb209_points_ggplot,
       width = 6, height = 10, dpi = 500)

# Time trend plots --------------------------------------------------------
# Create a data frame with tPCB
tpcb_2 <- data.frame(
  location = wdc$LocationName,
  date = SampleDate,
  tpcb = wdc$tPCB
)

plot.time.tPCB <- ggplot(tpcb_2, aes(y = tpcb,
                                   x = format(date,'%Y'),
                                   color = factor(location),
                                   fill = factor(location))) +
  geom_point(shape = 21, color = "black", size = 1.5) +
  xlab("") +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB (pg/L)"))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 10/20) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.height = unit(0.1, "cm"), # Adjust the height of the legend key
        legend.spacing.y = unit(0.1, "cm")) + 
  annotation_logticks(sides = "l") +
  scale_fill_viridis(discrete = TRUE, name = "Location") +
  guides(fill = guide_legend(ncol = 1))

# See plot
plot.time.tPCB

# Save plot in folder
ggsave("Output/Plots/Global/tPCBTime.png", plot = plot.time.tPCB,
       width = 10, height = 8, dpi = 300)

# Individual PCB trend plots
# Select the PCB column you want plot
i <- "PCB5.8"
i <- "PCB11"
i <- "PCB20.21.28.31.33.50.53"
i <- "PCB35"
i <- "PCB43.49.52.69.73"
i <- "PCB44.47.65"
i <- "PCB45.51"
i <- "PCB67"
i <- "PCB68"
i <- "PCB77.85.110.111.115.116.117"

# Create a data frame with the selected PCB column
pcbi <- data.frame(
  location = wdc$LocationName,
  date = SampleDate,
  pcb_value = wdc[[i]]  # Renamed the column to pcb_value
)

# Remove rows with NA values or values equal to 0 in the selected PCB column
pcbi <- pcbi[complete.cases(pcbi$pcb_value) & pcbi$pcb_value != 0, ]

# Plotting
plot.time.pcbi <- ggplot(pcbi, aes(y = pcb_value,
                                   x = format(date,'%Y'),
                                   color = factor(location),
                                   fill = factor(location))) +
  geom_point(shape = 21, color = "black", size = 1.5) +
  xlab("") +
  ylab(bquote(bold(.(i) ~ "(pg/L)"))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 10/10) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.height = unit(0.1, "cm"), # Adjust the height of the legend key
        legend.spacing.y = unit(0.1, "cm")) + 
  annotation_logticks(sides = "l") +
  scale_fill_viridis(discrete = TRUE, name = "Location") +
  guides(fill = guide_legend(ncol = 1)) 

# See plot
plot.time.pcbi

# Change name of plot.time.pcbi to the name of the congener in "i"
plot.time.pcb77 <- plot.time.pcbi

# plot all together
plot_combined <- plot_grid(plot.time.pcb5.8, plot.time.pcb11, plot.time.pcb20,
                           plot.time.pcb35, plot.time.pcb43, plot.time.pcb44,
                           plot.time.pcb45, plot.time.pcb67, plot.time.pcb68,
                           plot.time.pcb77,
                           ncol = 2, align = 'hv')

plot_combined

# Plot Aroclor congeners
plot_combined_1 <- plot_grid(plot.time.pcb5.8, plot.time.pcb20,
                             plot.time.pcb45, plot.time.pcb77,
                             ncol = 2, align = 'hv')

# See plot
plot_combined_1

# Save plot in folder
ggsave("Output/Plots/Global/aroclorTime.png", plot = plot_combined_1,
       width = 12, height = 7, dpi = 300)

# Plot non-Aroclor congeners
plot_combined_2 <- plot_grid(plot.time.pcb11, plot.time.pcb35,
                             plot.time.pcb67, plot.time.pcb68,
                             ncol = 2, align = 'hv')

# See plot
plot_combined_2

# Save plot in folder
ggsave("Output/Plots/Global/nonaroclorTime.png", plot = plot_combined_2,
       width = 12, height = 7, dpi = 300)


## Water PCB concentrations data analysis per site
## Anacostia River
## data only tPCB
## Random forest model

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("zoo")
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("tidyr")
install.packages('patchwork')
install.packages("scales")
install.packages("sf")
install.packages("units")
install.packages("sfheaders")
install.packages('ranger')
install.packages('caret')

# Load libraries
{
  library(ggplot2)
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(zoo) # yields seasons
  library(dataRetrieval) # read data from USGS
  library(reshape)
  library(tidyr) # function gather
  library(patchwork) # combine plots
  library(sf) # Create file to be used in Google Earth
  library(units)
  library(ranger) # Random Forest functions
  library(caret) # For cross-validation
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select anrsatonic River data ---------------------------------------------------
anr <- wdc[str_detect(wdc$LocationName, 'Anacostia River'),]

# Located northern location & calculate distance to other locations -------
{
  # Identify the northern sample based on the maximum Longitude
  index_of_northern_sample <- which.max(anr$Longitude)
  # Extract coordinates for the northern sample
  northern_sample_latitude <- anr$Latitude[index_of_northern_sample]
  northern_sample_longitude <- anr$Longitude[index_of_northern_sample]
  # Define source coordinates for the northern sample
  northern_source <- c(Latitude = northern_sample_latitude,
                       Longitude = northern_sample_longitude)
  # Create an sf point for the northern source
  northern_source_sf <- st_sfc(st_point(c(northern_source["Longitude"],
                                          northern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(northern_source_sf) <- 4326
  # Transform northern_source_sf to UTM Zone 18N (EPSG:32618)
  northern_source_sf_utm <- st_transform(northern_source_sf, 32618)
  # Convert the data frame to an sf object for the northern sample (anr)
  sf_anr <- st_as_sf(anr, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_anr <- st_set_crs(sf_anr, 4326)
  # Transform to UTM Zone 18N
  sf_anr_utm <- st_transform(sf_anr, 32618)
  # Calculate distances in meters from each location to northern source
  distances_meters_anr <- st_distance(sf_anr_utm, northern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_anr <- units::set_units(distances_meters_anr, "km")
  # Extract numeric values and assign to the DistanceToNorthernSource column
  anr$DistanceToNorthernLocation <- as.numeric(distances_km_anr[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  anr$SampleDate <- as.Date(anr$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(anr$SampleDate),
                                  min(as.Date(anr$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(anr$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  anr.tpcb <- cbind(factor(anr$SiteID), anr$SampleDate, as.matrix(anr$tPCB),
                    data.frame(time.day), season.s, anr$DistanceToNorthernLocation)
  # Add column names
  colnames(anr.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToNorthernLocation")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # https://maps.waterdata.usgs.gov/mapper/index.html
  # Include flow data from USGS station Anacostia River
  siteanrN1 <- "01649500" # flow @ NORTHEAST BRANCH ANACOSTIA RIVER AT RIVERDALE, MD
  siteanrN2 <- "01651000" # flow @ NORTHWEST BR ANACOSTIA RIVER NR HYATTSVILLE, MD
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteanrN1, paramflow,
                     min(anr.tpcb$date), max(anr.tpcb$date))
  flow.2 <- readNWISdv(siteanrN2, paramflow,
                       min(anr.tpcb$date), max(anr.tpcb$date))
  temp.1 <- readNWISdv(siteanrN1, paramtemp,
                     min(anr.tpcb$date), max(anr.tpcb$date))
  # Add USGS data to anr.tpcb.2, matching dates, conversion to m3/s
  anr.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(anr.tpcb$date, flow.1$Date)]
  anr.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(anr.tpcb$date, flow.2$Date)]
  anr.tpcb$temp.1 <- 273.15 + temp.1$X_00010_00003[match(anr.tpcb$date,
                                                         temp.1$Date)]
}

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
anr.tpcb <- select(anr.tpcb, -c(date, flow.1))

# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(anr.tpcb), 0.8 * nrow(anr.tpcb))
train_data <- anr.tpcb[train_indices, ]
test_data <- anr.tpcb[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(5, 10, 20)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
# Use flow.2
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.2 + temp.1 +
    DistanceToNorthernLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.2 + temp.1 +
    DistanceToNorthernLocation,
  data = train_data,
  num.trees = 1000, # need to manualy modify this parameter
  mtry = best_mtry,
  importance = 'permutation',
  seed = 123
)

# Get predictions on the test data
predictions <- predict(final_rf_model, data = test_data)$predictions

# Evaluate model performance
mse <- mean((predictions - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)

# Pearson correlation
correlation <- cor(log10(test_data$tPCB), predictions)

# Calculate Factor2
compare_df <- data.frame(observed = test_data$tPCB,
                         predicted = 10^predictions)
compare_df$factor2 <- compare_df$observed / compare_df$predicted
factor2_percentage <- nrow(compare_df[compare_df$factor2 > 0.5 & compare_df$factor2 < 2, ]) / nrow(compare_df) * 100

# Extract parameters from the random forest model
rf_parameters <- data.frame(
  NumTrees = final_rf_model$num.trees,
  Mtry = final_rf_model$mtry,
  SplitRule = final_rf_model$splitrule,
  MinNodeSize = final_rf_model$min.node.size
)

# Combine everything into a single data frame
results_rf_tPCB <- data.frame(
  Observation = test_data$tPCB,
  Predicted = 10^predictions,
  RMSE = rmse,
  Pearson = correlation,
  Factor2_Percentage = factor2_percentage,
  rf_parameters
)

results_rf_tPCB <- cbind(location = "Anacostia River", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRF)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/AnacostiaRiver/AnacostiaRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)


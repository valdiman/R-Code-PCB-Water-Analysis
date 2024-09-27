## Water PCB concentrations data analysis per site
## Housatonic River
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @ https://semspub.epa.gov/work/01/574882.pdf
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
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source1 <- c(Latitude = 42.456479, Longitude = -73.217587) # GE Pittsfield
  # Create an sf point for the source
  source1_sf <- st_sfc(st_point(c(source1["Longitude"], source1["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source1_sf) <- 4326
  # Transform source1_sf to UTM Zone 18N (EPSG:32618)
  source1_sf_utm <- st_transform(source1_sf, 32618)
  # Convert the data frame to an sf object
  sf_hou <- st_as_sf(hou, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_hou <- st_set_crs(sf_hou, 4326)
  # Transform to UTM Zone 18N
  sf_hou_utm <- st_transform(sf_hou, 32618)
  # Calculate distances in meters from each location to the source
  distances_meters <- st_distance(sf_hou_utm, source1_sf_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hou$DistanceToSource <- as.numeric(distances_km[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  hou$SampleDate <- as.Date(hou$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(hou$SampleDate),
                                  min(as.Date(hou$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$SampleDate, as.matrix(hou$tPCB),
                    data.frame(time.day), season.s, hou$DistanceToSource)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToSource")
}

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
# Use flow.1 (remove flow.2)
hou.tpcb <- select(hou.tpcb, -c(date, flow.2))

# Train-Test Split
set.seed(123)

train_indices <- sample(1:nrow(hou.tpcb), 0.8 * nrow(hou.tpcb))
train_data <- hou.tpcb[train_indices, ]
test_data <- hou.tpcb[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(5, 10, 20)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.1 + DistanceToSource,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.1 + DistanceToSource,
  data = train_data,
  num.trees = 1000, # need to manualy modify this parameter
  mtry = best_mtry,
  importance = 'permutation',
  seed = 123
)

# Get predictions on the test data
predictions <- predict(final_rf_model, data = test_data)$predictions

# Calculate RMSE
mse <- mean((predictions - log10(test_data$tPCB))^2)
rmse <- sqrt(mse)

# Calculate Pearson correlation
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

results_rf_tPCB <- cbind(location = "Housatonic River", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot using ggplot2
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(10^3, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^3, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/HousatonicRiver/HousatonicRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

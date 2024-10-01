## Water PCB concentrations data analysis per site
## Richardson Hill Road Landfill
## Arolclor method
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
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select Richardson Hill Road Landfill data ---------------------------------------------------
rhr <- wdc[str_detect(wdc$LocationName, 'Richardson Hill Road Landfill'),]

# Located southern location & calculate distance to other locations -------
{
  # Identify the northern sample based on the maximum Longitude
  index_of_northern_sample <- which.min(rhr$Longitude)
  # Extract coordinates for the northern sample
  northern_sample_latitude <- rhr$Latitude[index_of_northern_sample]
  northern_sample_longitude <- rhr$Longitude[index_of_northern_sample]
  # Define source coordinates for the northern sample
  northern_source <- c(Latitude = northern_sample_latitude,
                       Longitude = northern_sample_longitude)
  # Create an sf point for the northern source
  northern_source_sf <- st_sfc(st_point(c(northern_source["Longitude"],
                                          northern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(northern_source_sf) <- 4326
  # Transform northern_source_sf to UTM Zone 18N (EPSG:32610)
  northern_source_sf_utm <- st_transform(northern_source_sf, 32618)
  # Convert the data frame to an sf object for the northern sample (rhr)
  sf_rhr <- st_as_sf(rhr, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_rhr <- st_set_crs(sf_rhr, 4326)
  # Transform to UTM Zone 18N
  sf_rhr_utm <- st_transform(sf_rhr, 32618)
  # Calculate distances in meters from each location to northern source
  distances_meters_rhr <- st_distance(sf_rhr_utm, northern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_rhr <- units::set_units(distances_meters_rhr, "km")
  # Extract numeric values and assign to the DistanceToSouthernSource column
  rhr$DistanceToSouthernLocation <- as.numeric(distances_km_rhr[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  rhr$SampleDate <- as.Date(rhr$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(rhr$SampleDate),
                                  min(as.Date(rhr$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(rhr$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  rhr.tpcb <- cbind(factor(rhr$SiteID), as.matrix(rhr$tPCB), data.frame(time.day),
                    season.s, rhr$DistanceToSouthernLocation)
  # Add column names
  colnames(rhr.tpcb) <- c("SiteID", "tPCB", "time", "season",
                          "DistanceToSouthernLocation")
}

# Random Forest Model tPCB ------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(rhr.tpcb), 0.8 * nrow(rhr.tpcb))
train_data <- rhr.tpcb[train_indices, ]
test_data <- rhr.tpcb[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(2, 3, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + DistanceToSouthernLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + DistanceToSouthernLocation,
  data = train_data,
  num.trees = 100, # need to manually modify this parameter
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

results_rf_tPCB <- cbind(location = "Richardson Hill Road Landfill", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/Richardson/RichardsonRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(100, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(100, 10^7),
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
ggsave("Output/Plots/Sites/ObsPred/Richardson/RichardsonRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)


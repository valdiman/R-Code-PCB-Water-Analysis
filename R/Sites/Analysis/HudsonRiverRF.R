## Water PCB concentrations data analysis per site
## Hudson River

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
install.packages('viridis')

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
  library(viridis) # For color blind
}

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select Hudson River data ---------------------------------------------------
hud <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Add distance to source --------------------------------------------------
{
  # Define source coordinates
  source1 <- c(Latitude = 43.295369, Longitude = -73.590631)  # GE Hudson Falls Plant
  source2 <- c(Latitude = 43.28639, Longitude = -73.588380)  # GE Fort Edward Plant
  # Create an sf point for the source
  source1_sf <- st_sfc(st_point(c(source1["Longitude"], source1["Latitude"])))
  source2_sf <- st_sfc(st_point(c(source2["Longitude"], source1["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(source1_sf) <- 4326
  st_crs(source2_sf) <- 4326
  # Transform source1_sf to UTM Zone 18N (EPSG:32618)
  source1_sf_utm <- st_transform(source1_sf, 32618)
  # Transform source2_sf to UTM Zone 18N (EPSG:32618)
  source2_sf_utm <- st_transform(source2_sf, 32618)
  # Convert the data frame to an sf object
  sf_hud <- st_as_sf(hud, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_hud <- st_set_crs(sf_hud, 4326)
  # Transform to UTM Zone 18N
  sf_hud_utm <- st_transform(sf_hud, 32618)
  # Calculate distances in meters from each location to source1
  distances_meters1 <- st_distance(sf_hud_utm, source1_sf_utm)
  # Convert distances to kilometers
  distances_km1 <- units::set_units(distances_meters1, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hud$DistanceToSource1 <- as.numeric(distances_km1[, 1])
  # Calculate distances in meters from each location to source2
  distances_meters2 <- st_distance(sf_hud_utm, source2_sf_utm)
  # Convert distances to kilometers
  distances_km2 <- units::set_units(distances_meters2, "km")
  # Extract numeric values and assign to the DistanceToSource column
  hud$DistanceToSource2 <- as.numeric(distances_km2[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  hud$SampleDate <- as.Date(hud$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(hud$SampleDate),
                                  min(as.Date(hud$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hud.tpcb <- cbind(factor(hud$SiteID), hud$SampleDate, as.matrix(hud$tPCB),
                    data.frame(time.day), season.s, hud$DistanceToSource1,
                    hud$DistanceToSource2)
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceSource1", "DistanceSource2")
}

# Remove site -------------------------------------------------------------
## Remove site Bakers Falls. Upstream source
## North Bakers Falls = WCPCB-HUD006 and
## South Bakers Falls = WCPCB-HUD006.
hud.tpcb.1 <- subset(hud.tpcb, SiteID != c("WCPCB-HUD006"))
hud.tpcb.1 <- subset(hud.tpcb.1, SiteID != c("WCPCB-HUD010"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Hudson River
{
  sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY No temp!
  sitehudN2 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY, no temp!
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN4 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitehudN1, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.3 <- readNWISdv(sitehudN3, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.4 <- readNWISdv(sitehudN4, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  
  # Add USGS data to hud.tpcb.2, matching dates
  hud.tpcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.1$Date)]
  hud.tpcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.2$Date)]
  hud.tpcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.3$Date)]
  hud.tpcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.4$Date)]
  hud.tpcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.tpcb.1$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  hud.tpcb.2 <- na.omit(hud.tpcb.1)
}

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
# Use flow.2 and DistanceSource1
hud.tpcb.2 <- select(hud.tpcb.2, -c(date, flow.1, flow.3, flow.4,
                                    DistanceSource2))

# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(hud.tpcb.2), 0.8 * nrow(hud.tpcb.2))
train_data <- hud.tpcb.2[train_indices, ]
test_data <- hud.tpcb.2[-train_indices, ]

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
  log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceSource1,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceSource1,
  data = train_data,
  num.trees = 5000, # need to manualy modify this parameter
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

results_rf_tPCB <- cbind(location = "Hudson River", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(10^2, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2, 10^6),
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
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -----------------------------------
{
  hud.pcb <- subset(hud, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  hud.pcb <- subset(hud.pcb, select = -c(A1016:DistanceToSource2))
  # Log10 individual PCBs 
  hud.pcb <- log10(hud.pcb)
  # Replace -inf to NA
  hud.pcb <- do.call(data.frame,
                     lapply(hud.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  hud.pcb.1 <- hud.pcb[,
                       -which(colSums(is.na(hud.pcb))/nrow(hud.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(hud$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to source
  DistanceToSource1 <- hud$DistanceToSource1
  # Add date and time to hud.pcb.1
  hud.pcb.1 <- cbind(hud.pcb.1, as.factor(hud$SiteID), SampleDate,
                     data.frame(time.day), season.s, DistanceToSource1)
  # Remove site Bakers Falls. Upstream source
  # North Bakers Falls = WCPCB-HUD006 and
  # South Bakers Falls = WCPCB-HUD010.
  # Remove rows with SiteID equal to "WCPCB-HUD006" or "WCPCB-HUD010"
  hud.pcb.1 <- hud.pcb.1[!(hud.pcb.1$`as.factor(hud$SiteID)` %in% c("WCPCB-HUD006",
                                                              "WCPCB-HUD010")), ]
  # Include flow data from USGS station Hudson River
  # sitehudN3 for flow and sitehudN5 for water temperature
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Retrieve USGS data
  # Flow (ft3/s)
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  
  # Add USGS data to hud.tpcb.1 matching dates
  hud.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                      flow.2$Date)]
  hud.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with flow.2 = NA
  hud.pcb.2 <- hud.pcb.1[!is.na(hud.pcb.1$flow.2), ]
  # Remove metadata not use in the random forest
  hud.pcb.2 <- hud.pcb.2[, !(names(hud.pcb.2) %in% c("SampleDate"))]
}

# Function to perform grid search for ranger model hyperparameters
perform_grid_search <- function(congener_name, seed = 123) {
  # Extract data for the specified PCB congener
  congener_data <- hud.pcb.2[[congener_name]]
  
  # Get indices of non-missing values in the specified PCB congener
  non_missing_indices <- which(!is.na(congener_data))
  
  # Extract covariate columns (excluding columns starting with "PCB")
  covariate_columns <- names(hud.pcb.2)[!grepl("^PCB", names(hud.pcb.2))]
  
  # Select rows with non-missing values in the specified PCB congener and extract covariates
  congener_covariate_data <- hud.pcb.2[non_missing_indices, covariate_columns]
  
  # Include the specified PCB congener column in the covariate data
  congener_covariate_data[[congener_name]] <- congener_data[non_missing_indices]
  
  # Count the number of non-missing observations in the specified PCB congener data
  num_obs_congener <- nrow(congener_covariate_data)
  
  # Calculate the number of observations for training and testing data
  num_train <- round(0.8 * num_obs_congener)
  num_test <- num_obs_congener - num_train
  
  # Sample indices for training and testing data
  set.seed(seed)
  train_indices <- sample(1:num_obs_congener, num_train)
  test_indices <- setdiff(1:num_obs_congener, train_indices)
  
  # Extract training and testing data
  train_data <- congener_covariate_data[train_indices, ]
  test_data <- congener_covariate_data[test_indices, ]
  
  # Define hyperparameter grid
  hyperparameters <- expand.grid(
    num.trees = c(1000, 4000, 5000),  # Example values, adjust as needed
    mtry = c(4, 6, min(6, ncol(train_data))),  # Ensure mtry is less than or equal to the number of variables
    min.node.size = c(3, 4, 5)      # Example values, adjust as needed
  )
  
  # Initialize variables to store best hyperparameters and performance metrics
  best_hyperparameters <- NULL
  best_performance <- Inf
  
  # Iterate over each combination of hyperparameters
  for (i in 1:nrow(hyperparameters)) {
    # Train the ranger model with current hyperparameters
    ranger_model <- ranger(
      dependent.variable.name = congener_name,
      data = train_data,
      num.trees = hyperparameters$num.trees[i],
      mtry = hyperparameters$mtry[i],
      min.node.size = hyperparameters$min.node.size[i],
      seed = seed
    )
    
    # Predict on the test set
    predictions <- predict(ranger_model, data = test_data)$predictions
    
    # Calculate RMSE
    mse <- mean((predictions - test_data[[congener_name]])^2, na.rm = TRUE)
    
    # Update best hyperparameters and performance metrics if current model performs better
    if (mse < best_performance) {
      best_performance <- mse
      best_hyperparameters <- hyperparameters[i, ]
    }
  }
  
  # Predict on the test set using the best hyperparameters
  best_ranger_model <- ranger(
    dependent.variable.name = congener_name,
    data = train_data,
    num.trees = best_hyperparameters$num.trees,
    mtry = best_hyperparameters$mtry,
    min.node.size = best_hyperparameters$min.node.size,
    seed = seed
  )
  
  best_predictions <- predict(best_ranger_model, data = test_data)$predictions
  
  # Create dataframe for predictions and actual values
  results_rf_PCBi <- data.frame(
    Location = rep("Hudson River", length(test_indices)),
    Congener = rep(congener_name, length(test_indices)),
    Test_Data = test_data[[congener_name]],
    Predicted_Data = best_predictions,
    RMSE = rep(sqrt(best_performance), length(test_indices)),
    Correlation = rep(cor(test_data[[congener_name]], best_predictions), length(test_indices)),
    NumTrees = rep(best_hyperparameters$num.trees, length(test_indices)),
    Mtry = rep(best_hyperparameters$mtry, length(test_indices)),
    SplitRule = rep("variance", length(test_indices)), # Example value, adjust as needed
    MinNodeSize = rep(best_hyperparameters$min.node.size, length(test_indices)),
    Factor2 = numeric(length(test_indices)), # Placeholder, to be calculated later
    stringsAsFactors = FALSE
  )
  
  # Calculate Factor2 for each Congener in results_rf_PCBi
  results_rf_PCBi$Factor2 <- with(results_rf_PCBi, {
    factor2 <- 10^(Test_Data) / 10^(Predicted_Data)
    as.numeric(factor2 > 0.5 & factor2 < 2)
  })
  
  # Calculate percentage of Factor2 for each Congener
  results_rf_PCBi$Factor2 <- aggregate(Factor2 ~ Congener,
                                       data = results_rf_PCBi, FUN = function(x) {
                                         sum(x == 1) / length(x) * 100
                                       })$Factor2
  
  # Return the dataframe
  return(results_rf_PCBi)
}

# Initialize an empty dataframe to store the results
results_rf_PCBi <- data.frame(
  Location = character(),
  Congener = character(),
  Test_Data = numeric(),
  Predicted_Data = numeric(),
  RMSE = numeric(),
  Correlation = numeric(),
  NumTrees = numeric(),
  Mtry = numeric(),
  SplitRule = character(),
  MinNodeSize = numeric(),
  Factor2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each PCB congener
for (congener_name in colnames(hud.pcb.2)) {
  if (grepl("^PCB", congener_name)) {
    # Perform grid search and get predictions
    result <- perform_grid_search(congener_name, seed = 123)
    
    # Append the results to the results_rf_PCBi dataframe
    results_rf_PCBi <- rbind(results_rf_PCBi, result)
  }
}

# Export results
write.csv(results_rf_PCBi,
          file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(results_rf_PCBi, aes(x = 10^(Test_Data), y = 10^(Predicted_Data))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/HudsonRiver/HudsonRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


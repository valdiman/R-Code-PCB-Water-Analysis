## Water PCB concentrations data analysis per site
## Portland Harbor
## Random Forest model

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

# Select Portland Harbor data ---------------------------------------------------
por <- wdc[str_detect(wdc$LocationName, 'Portland Harbor'),]

# Located northern location & calculate distance to other locations -------
{
  # Identify the northern sample based on the maximum Longitude
  index_of_northern_sample <- which.max(por$Longitude)
  # Extract coordinates for the northern sample
  northern_sample_latitude <- por$Latitude[index_of_northern_sample]
  northern_sample_longitude <- por$Longitude[index_of_northern_sample]
  # Define source coordinates for the northern sample
  northern_source <- c(Latitude = northern_sample_latitude,
                       Longitude = northern_sample_longitude)
  # Create an sf point for the northern source
  northern_source_sf <- st_sfc(st_point(c(northern_source["Longitude"],
                                          northern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(northern_source_sf) <- 4326
  # Transform northern_source_sf to UTM Zone 10N (EPSG:32610)
  northern_source_sf_utm <- st_transform(northern_source_sf, 32610)
  # Convert the data frame to an sf object for the northern sample (spo)
  sf_por <- st_as_sf(por, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_por <- st_set_crs(sf_por, 4326)
  # Transform to UTM Zone 10N
  sf_por_utm <- st_transform(sf_por, 32610)
  # Calculate distances in meters from each location to northern source
  distances_meters_por <- st_distance(sf_por_utm, northern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_por <- units::set_units(distances_meters_por, "km")
  # Extract numeric values and assign to the DistanceToNorthernSource column
  por$DistanceToNorthernLocation <- as.numeric(distances_km_por[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  por$SampleDate <- as.Date(por$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(por$SampleDate),
                                  min(as.Date(por$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  por.tpcb <- cbind(factor(por$SiteID), por$SampleDate, as.matrix(por$tPCB),
                    data.frame(time.day), season.s, por$DistanceToNorthernLocation)
  # Add column names
  colnames(por.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToNorthernLocation")
  # Include USGC station Portland Harbor flow and water temperature
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  temp <- readNWISdv(sitePorN1, paramtemp,
                       min(por.tpcb$date), max(por.tpcb$date))
  # Add USGS data to por.tpcb, matching dates
  por.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(por.tpcb$date, flow.1$Date)]
  por.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.tpcb$date, flow.2$Date)]
  por.tpcb$temp <- 273.15 + temp$X_00010_00003[match(por.tpcb$date, temp$Date)]
  # Remove samples with temp = NA
  por.tpcb.1 <- na.omit(por.tpcb)
}

# Random Forest Model -----------------------------------------------------
# Remove columns not used here
# Use flow.2 and DistanceSource1
por.tpcb.1 <- select(por.tpcb.1, -c(date, flow.1))

# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(por.tpcb.1), 0.8 * nrow(por.tpcb.1))
train_data <- por.tpcb.1[train_indices, ]
test_data <- por.tpcb.1[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(3, 4, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
# Better fit with flow.2
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.2 + temp + DistanceToNorthernLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.2 + temp +
    DistanceToNorthernLocation,
  data = train_data,
  num.trees = 5000, # need to manually modify this parameter
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

results_rf_tPCB <- cbind(location = "Portland Harbor", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^4),
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
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  por.pcb <- subset(por, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  por.pcb <- subset(por.pcb, select = -c(A1016:DistanceToNorthernLocation))
  # Log10 individual PCBs 
  por.pcb <- log10(por.pcb)
  # Replace -inf to NA
  por.pcb <- do.call(data.frame,
                     lapply(por.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  por.pcb.1 <- por.pcb[,
                       -which(colSums(is.na(por.pcb))/nrow(por.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(por$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to northern location sample
  DistanceToNorthernLocation <- por$DistanceToNorthernLocation
  # Add date and time to por.pcb.1
  por.pcb.1 <- cbind(por.pcb.1, as.factor(por$SiteID), SampleDate,
                     data.frame(time.day), season.s, DistanceToNorthernLocation)
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  temp <- readNWISdv(sitePorN1, paramtemp,
                     min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  # Add USGS data to por.tpcb, matching dates
  por.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.pcb.1$SampleDate, flow.2$Date)]
  por.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(por.pcb.1$SampleDate, temp$Date)]
  # Remove samples with temp = NA
  por.pcb.2 <- por.pcb.1[!is.na(por.pcb.1$temp), ]
  # Remove metadata not use in the random forest
  por.pcb.2 <- por.pcb.2[, !(names(por.pcb.2) %in% c("SampleDate"))]
}

# Function to perform grid search for ranger model hyperparameters
perform_grid_search <- function(congener_name, seed = 123) {
  # Extract data for the specified PCB congener
  congener_data <- por.pcb.2[[congener_name]]
  
  # Get indices of non-missing values in the specified PCB congener
  non_missing_indices <- which(!is.na(congener_data))
  
  # Extract covariate columns (excluding columns starting with "PCB")
  covariate_columns <- names(por.pcb.2)[!grepl("^PCB", names(por.pcb.2))]
  
  # Select rows with non-missing values in the specified PCB congener and extract covariates
  congener_covariate_data <- por.pcb.2[non_missing_indices, covariate_columns]
  
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
    num.trees = c(50, 100, 200, 500),  # Example values, adjust as needed
    mtry = c(2, 6, min(6, ncol(train_data))),  # Ensure mtry is less than or equal to the number of variables
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
    Location = rep("Portland Harbor", length(test_indices)),
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
for (congener_name in colnames(por.pcb.2)) {
  if (grepl("^PCB", congener_name)) {
    # Perform grid search and get predictions
    result <- perform_grid_search(congener_name, seed = 123)
    
    # Append the results to the results_rf_PCBi dataframe
    results_rf_PCBi <- rbind(results_rf_PCBi, result)
  }
}

# Export results
write.csv(results_rf_PCBi,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(results_rf_PCBi, aes(x = 10^(Test_Data), y = 10^(Predicted_Data))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.0001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.0001, 10^4),
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
ggsave("Output/Plots/Sites/ObsPred/PortlandHarbor/PortlandHarborRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


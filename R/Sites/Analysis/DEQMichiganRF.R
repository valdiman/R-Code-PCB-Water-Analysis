## Water PCB concentrations data analysis per site
## Data from DEQ dmihigan
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
# Data in pg/L
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Select DEQ Michigan data ---------------------------------------------------
# issue with ( ).
dmi <- wdc[str_detect(wdc$LocationName, fixed("DEQ (Michigan)")), ]

# Calculate central location. ---------------------------------------------
{
  # Calculate the mean latitude and longitude
  center_lat <- mean(dmi$Latitude)
  center_lon <- mean(dmi$Longitude)
  # Create a data frame for the center
  center_df <- data.frame(
    SiteID = "Center",
    Latitude = center_lat,
    Longitude = center_lon
  )
  # Convert the center data frame to an sf object
  sf_center <- st_as_sf(center_df, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326) for the center
  sf_center <- st_set_crs(sf_center, 4326)
  # Transform to UTM Zone 16 for the center
  sf_center_utm <- st_transform(sf_center, 32616)
  # Convert the data frame to an sf object
  sf_dmi <- st_as_sf(dmi, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_dmi <- st_set_crs(sf_dmi, 4326)
  # Transform to UTM Zone 16
  sf_dmi <- st_transform(sf_dmi, 32616)
  # Calculate distances in meters from each location to the center
  distances_meters <- st_distance(sf_dmi, sf_center_utm)
  # Convert distances to kilometers
  distances_km <- units::set_units(distances_meters, "km")
  # Extract numeric values and assign to the DistanceToCentroid column
  dmi$DistanceToCentroid <- as.numeric(distances_km[, 1])
}

# tPCB data preparation ---------------------------------------------------
{
  # Change date format
  dmi$SampleDate <- as.Date(dmi$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(dmi$SampleDate),
                                  min(as.Date(dmi$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(dmi$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  dmi.tpcb <- cbind(factor(dmi$SiteID), as.matrix(dmi$tPCB),
                    data.frame(time.day), season.s, dmi$DistanceToCentroid)
  
  # Add column names
  colnames(dmi.tpcb) <- c("SiteID", "tPCB", "time", "season",
                          "DistanceToCentroid")
}

# Random Forest Model tPCB ------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Train-test split
train_indices <- sample(1:nrow(dmi.tpcb), 0.8 * nrow(dmi.tpcb))
train_data <- dmi.tpcb[train_indices, ]
test_data <- dmi.tpcb[-train_indices, ]

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
  log10(tPCB) ~ time + SiteID + season + DistanceToCentroid,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + DistanceToCentroid,
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

results_rf_tPCB <- cbind(location = "DEQ MI", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/DEQMichigan/DEQMIRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(1, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^6),
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
ggsave("Output/Plots/Sites/ObsPred/DEQMichigan/DEQMIRFObsPredtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
{
  dmi.pcb <- subset(dmi, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  dmi.pcb <- subset(dmi.pcb, select = -c(A1016:DistanceToCentroid))
  # Log10 individual PCBs 
  dmi.pcb <- log10(dmi.pcb)
  # Replace -inf to NA
  dmi.pcb <- do.call(data.frame,
                     lapply(dmi.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  dmi.pcb.1 <- dmi.pcb[,
                       -which(colSums(is.na(dmi.pcb))/nrow(dmi.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(dmi$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Add distance to the centroid
  centroid <- dmi$DistanceToCentroid
  # Include season
  yq.s <- as.yearqtr(as.yearmon(dmi$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to dmi.pcb.1
  dmi.pcb.1 <- cbind(dmi.pcb.1, as.factor(dmi$SiteID), data.frame(time.day),
                     season.s, centroid)
}

# Function to perform grid search for ranger model hyperparameters
perform_grid_search <- function(congener_name, seed = 123) {
  # Extract data for the specified PCB congener
  congener_data <- dmi.pcb.1[[congener_name]]
  
  # Get indices of non-missing values in the specified PCB congener
  non_missing_indices <- which(!is.na(congener_data))
  
  # Extract covariate columns (excluding columns starting with "PCB")
  covariate_columns <- names(dmi.pcb.1)[!grepl("^PCB", names(dmi.pcb.1))]
  
  # Select rows with non-missing values in the specified PCB congener and extract covariates
  congener_covariate_data <- dmi.pcb.1[non_missing_indices, covariate_columns]
  
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
    num.trees = c(500, 700, 1000),  # Example values, adjust as needed
    mtry = c(2, 3, min(3, ncol(train_data))),  # Ensure mtry is less than or equal to the number of variables
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
    Location = rep("DEQ MI", length(test_indices)),
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
for (congener_name in colnames(dmi.pcb.1)) {
  if (grepl("^PCB", congener_name)) {
    # Perform grid search and get predictions
    result <- perform_grid_search(congener_name, seed = 123)
    
    # Append the results to the results_rf_PCBi dataframe
    results_rf_PCBi <- rbind(results_rf_PCBi, result)
  }
}

# Export results
write.csv(results_rf_PCBi,
          file = "Output/Data/Sites/csv/DEQMichigan/DEQMIRFPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(results_rf_PCBi, aes(x = 10^(Test_Data), y = 10^(Predicted_Data))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.001, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^6),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) + # 1:1 line
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) +   # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/DEQMichigan/DEQMIRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


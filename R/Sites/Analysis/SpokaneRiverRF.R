## Water PCB concentrations data analysis per site
## Spokane River
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
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select Spokane River data ---------------------------------------------------
spo <- wdc[str_detect(wdc$LocationName, 'Spokane River'),]

# Located eastern location & calculate distance to other locations ---------
{
  # Identify the eastern sample based on the maximum Latitude
  index_of_eastern_sample <- which.max(spo$Latitude)
  # Extract coordinates for the eastern sample
  eastern_sample_latitude <- spo$Latitude[index_of_eastern_sample]
  eastern_sample_longitude <- spo$Longitude[index_of_eastern_sample]
  # Define source coordinates for the eastern sample
  eastern_source <- c(Latitude = eastern_sample_latitude,
                      Longitude = eastern_sample_longitude)
  # Create an sf point for the eastern source
  eastern_source_sf <- st_sfc(st_point(c(eastern_source["Longitude"],
                                         eastern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(eastern_source_sf) <- 4326
  # Transform eastern_source_sf to UTM Zone 10N (EPSG:32610)
  eastern_source_sf_utm <- st_transform(eastern_source_sf, 32610)
  # Convert the data frame to an sf object for the eastern sample (spo)
  sf_spo <- st_as_sf(spo, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_spo <- st_set_crs(sf_spo, 4326)
  # Transform to UTM Zone 10N
  sf_spo_utm <- st_transform(sf_spo, 32610)
  # Calculate distances in meters from each location to eastern source
  distances_meters_spo <- st_distance(sf_spo_utm, eastern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_spo <- units::set_units(distances_meters_spo, "km")
  # Extract numeric values and assign to the DistanceToEasternSource column
  spo$DistanceToEasternLocation <- as.numeric(distances_km_spo[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  spo$SampleDate <- as.Date(spo$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(spo$SampleDate) - min(as.Date(spo$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  spo.tpcb <- cbind(factor(spo$SiteID), spo$SampleDate,
                    as.matrix(spo$tPCB), data.frame(time.day), season.s,
                    spo$DistanceToEasternLocation)
  # Add column names
  colnames(spo.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToEasternLocation")
  # Include flow data from USGS station Spokane River
  siteSpoN1 <- "12417650" # SPOKANE RIVER BLW BLACKWELL NR COEUR D ALENE ID
  siteSpoN2 <- "12419000" # Spokane River near Post Falls, ID
  siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
  siteSpoN4 <- "12424000" # Hangman Creek at Spokane, WA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  #paramtemp <- "00010" # water temperature, C No data
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteSpoN1, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.2 <- readNWISdv(siteSpoN2, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.3 <- readNWISdv(siteSpoN3, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  flow.4 <- readNWISdv(siteSpoN4, paramflow,
                       min(spo.tpcb$date), max(spo.tpcb$date))
  
  # Add USGS data to spo.tpcb, matching dates
  spo.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(spo.tpcb$date,
                                                     flow.1$Date)]
  spo.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(spo.tpcb$date,
                                                     flow.2$Date)]
  spo.tpcb$flow.3 <- 0.03*flow.3$X_00060_00003[match(spo.tpcb$date,
                                                     flow.3$Date)]
  spo.tpcb$flow.4 <- 0.03*flow.4$X_00060_00003[match(spo.tpcb$date,
                                                     flow.4$Date)]
}

# Remove site -------------------------------------------------------------
## Sample sites not located at the Spokane River
{
  spo.tpcb.1 <- subset(spo.tpcb, SiteID != c("WCPCB-SPR002")) # City of Spokane WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR005")) # Regional WRF
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR006")) # Inland Empire paper
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR008")) # Kaiser Aluminum
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR010")) # Liberty Lake sewer
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR011")) # Post Falls WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR013")) # Coeur d'Alene WWTP
  spo.tpcb.1 <- subset(spo.tpcb.1, SiteID != c("WCPCB-SPR015")) # Hangman Creek
}

# Random Forest Model tPCB  -----------------------------------------------
# Remove columns not used here
# Using flow.3
spo.tpcb.1 <- select(spo.tpcb.1, -c(date, flow.1, flow.2, flow.4))

# Set seed for reproducibility
set.seed(123)

# Train-Test Split
train_indices <- sample(1:nrow(spo.tpcb.1), 0.8 * nrow(spo.tpcb.1))
train_data <- spo.tpcb.1[train_indices, ]
test_data <- spo.tpcb.1[-train_indices, ]

# Define hyperparameter grid
param_grid <- expand.grid(
  mtry = seq(1, ncol(train_data) - 1),  # Adjust mtry values based on your data
  splitrule = c("gini", "extratrees"),
  min.node.size = c(3, 4, 5)
)

# Prepare training control
ctrl <- trainControl(method = "cv", number = 5, search = "grid")

# Perform grid search with cross-validation using ranger
rf_model <- train(
  log10(tPCB) ~ time + SiteID + season + flow.3 + DistanceToEasternLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow.3 +
    DistanceToEasternLocation,
  data = train_data,
  num.trees = 1000, # need to manually modify this parameter
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

results_rf_tPCB <- cbind(location = "Spokane River", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
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
ggsave("Output/Plots/Sites/ObsPred/SpokaneRiver/SpokaneRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)

# Random Forest Model individual PCBs -------------------------------------
{
  # Remove metadata
  spo.pcb <- subset(spo, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  spo.pcb <- subset(spo.pcb, select = -c(A1016:DistanceToEasternLocation))
  # Log10 individual PCBs 
  spo.pcb <- log10(spo.pcb)
  # Replace -inf to NA
  spo.pcb <- do.call(data.frame,
                     lapply(spo.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  spo.pcb.1 <- spo.pcb[, colSums(is.na(spo.pcb))/nrow(spo.pcb) <= 0.7]
  # Change date format
  SampleDate <- as.Date(spo$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(spo$SampleDate) - min(as.Date(spo$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to eastern location
  eastern <- spo$DistanceToEasternLocation
  # Add date and time to spo.pcb.1
  spo.pcb.1 <- cbind(spo.pcb.1, as.factor(spo$SiteID), SampleDate,
                     data.frame(time.day), season.s, eastern)
  # Include flow data from USGS station Spokane River
  siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.3 <- readNWISdv(siteSpoN3, paramflow,
                       min(spo.pcb.1$SampleDate), max(spo.pcb.1$SampleDate))
  # Add USGS data to spo.tpcb, matching dates
  spo.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(spo.pcb.1$SampleDate,
                                                     flow.3$Date)]
  # Sample sites not located at the Spokane River
  spo.pcb.2 <- subset(spo.pcb.1, !(as.factor(spo$SiteID) %in% c("WCPCB-SPR002",
                                                                "WCPCB-SPR005",
                                                                "WCPCB-SPR006",
                                                                "WCPCB-SPR008",
                                                                "WCPCB-SPR010",
                                                                "WCPCB-SPR011",
                                                                "WCPCB-SPR013",
                                                                "WCPCB-SPR015")))
  # Remove metadata not use in the random forest
  spo.pcb.2 <- spo.pcb.2[, !(names(spo.pcb.2) %in% c("SampleDate"))]
}

# Function to perform grid search for ranger model hyperparameters
perform_grid_search <- function(congener_name, seed = 123) {
  # Extract data for the specified PCB congener
  congener_data <- spo.pcb.2[[congener_name]]
  
  # Get indices of non-missing values in the specified PCB congener
  non_missing_indices <- which(!is.na(congener_data))
  
  # Extract covariate columns (excluding columns starting with "PCB")
  covariate_columns <- names(spo.pcb.2)[!grepl("^PCB", names(spo.pcb.2))]
  
  # Select rows with non-missing values in the specified PCB congener and extract covariates
  congener_covariate_data <- spo.pcb.2[non_missing_indices, covariate_columns]
  
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
    num.trees = c(200, 300, 500),  # Example values, adjust as needed
    mtry = c(2, 5, min(5, ncol(train_data))),  # Ensure mtry is less than or equal to the number of variables
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
    Location = rep("Spokane River", length(test_indices)),
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
for (congener_name in colnames(spo.pcb.2)) {
  if (grepl("^PCB", congener_name)) {
    # Perform grid search and get predictions
    result <- perform_grid_search(congener_name, seed = 123)
    
    # Append the results to the results_rf_PCBi dataframe
    results_rf_PCBi <- rbind(results_rf_PCBi, result)
  }
}

# Export results
write.csv(results_rf_PCBi,
          file = "Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(results_rf_PCBi, aes(x = 10^(Test_Data), y = 10^(Predicted_Data))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.01, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.01, 10^3),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted concentration PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Print the plot
print(plotRFPCBi)

# Save plot in folder
ggsave("Output/Plots/Sites/ObsPred/SpokaneRiver/SpokaneRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


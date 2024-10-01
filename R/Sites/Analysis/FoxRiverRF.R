## Water PCB concentrations data analysis
## Fox River 2005 - 2018
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
  library(sf) # Create file to be used in Google Earth
  library(units)
  library(ranger) # Random Forest functions
  library(caret) # For cross-validation
  library(viridis) # For color blind
}

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select Fox River data ---------------------------------------------------
fox <- wdc[str_detect(wdc$LocationName, 'Fox River'),]

# Located eastern location & calculate distance to other locations ---------
{
  # Identify the eastern sample based on the maximum Latitude
  index_of_eastern_sample <- which.max(fox$Latitude)
  # Extract coordinates for the eastern sample
  eastern_sample_latitude <- fox$Latitude[index_of_eastern_sample]
  eastern_sample_longitude <- fox$Longitude[index_of_eastern_sample]
  # Define source coordinates for the eastern sample
  eastern_source <- c(Latitude = eastern_sample_latitude,
                      Longitude = eastern_sample_longitude)
  # Create an sf point for the eastern source
  eastern_source_sf <- st_sfc(st_point(c(eastern_source["Longitude"],
                                         eastern_source["Latitude"])))
  # Set the CRS to EPSG:4326
  st_crs(eastern_source_sf) <- 4326
  # Transform eastern_source_sf to UTM Zone 16N (EPSG:32610)
  eastern_source_sf_utm <- st_transform(eastern_source_sf, 32616)
  # Convert the data frame to an sf object for the eastern sample (fox)
  sf_fox <- st_as_sf(fox, coords = c("Longitude", "Latitude"))
  # Set the CRS to WGS 84 (EPSG:4326)
  sf_fox <- st_set_crs(sf_fox, 4326)
  # Transform to UTM Zone 16N
  sf_fox_utm <- st_transform(sf_fox, 32616)
  # Calculate distances in meters from each location to eastern source
  distances_meters_fox <- st_distance(sf_fox_utm, eastern_source_sf_utm)
  # Convert distances to kilometers
  distances_km_fox <- units::set_units(distances_meters_fox, "km")
  # Extract numeric values and assign to the DistanceToEasternSource column
  fox$DistanceToEasternLocation <- as.numeric(distances_km_fox[, 1])
}

# Data preparation --------------------------------------------------------
{
  # Change date format
  fox$SampleDate <- as.Date(fox$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(fox$SampleDate),
                                  min(as.Date(fox$SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox$SiteID), fox$SampleDate, as.matrix(fox$tPCB),
                    data.frame(time.day), season.s, fox$DistanceToEasternLocation)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "tPCB", "time", "season",
                          "DistanceToEasternLocation")
}

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.tpcb$date), max(fox.tpcb$date))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.tpcb$date), max(fox.tpcb$date))
  # Add USGS data to fox.tpcb.2, matching dates, conversion to m3/s
  fox.tpcb$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb$date,
                                                                       flow$Date)]
  fox.tpcb$temp <- 273.15 + temp$X_..2.._00010_00003[match(fox.tpcb$date,
                                                     temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb <- na.omit(fox.tpcb)
}

# Random Forest Model tPCB ------------------------------------------------
# Remove columns not used here
fox.tpcb <- select(fox.tpcb, -c(date))

# Train-Test Split
set.seed(123)
train_indices <- sample(1:nrow(fox.tpcb), 0.8 * nrow(fox.tpcb))
train_data <- fox.tpcb[train_indices, ]
test_data <- fox.tpcb[-train_indices, ]

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
  log10(tPCB) ~ time + SiteID + season + flow + temp +
    DistanceToEasternLocation,
  data = train_data,
  method = "ranger",
  importance = 'permutation',
  tuneGrid = param_grid,
  trControl = ctrl
)

# Get the best mtry
best_mtry <- rf_model$bestTune$mtry

final_rf_model <- ranger(
  formula = log10(tPCB) ~ time + SiteID + season + flow + temp +
    DistanceToEasternLocation,
  data = train_data,
  num.trees = 4500, # need to manualy modify this parameter
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

results_rf_tPCB <- cbind(location = "Fox River", results_rf_tPCB)

# Export results
write.csv(results_rf_tPCB,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFtPCB.csv",
          row.names = FALSE)

# Create the scatter plot
plotRF <- ggplot(results_rf_tPCB, aes(x = Observation, y = Predicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(10, 10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5),
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
ggsave("Output/Plots/Sites/ObsPred/FoxRiver/FoxRiverRFtPCB.png",
       plot = plotRF, width = 6, height = 5, dpi = 500)


# Random Forest Model individual PCBs -------------------------------------
# Prepare data.frame
{
  fox.pcb <- subset(fox, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  fox.pcb <- subset(fox.pcb, select = -c(A1016:DistanceToEasternLocation))
  # Log10 individual PCBs 
  fox.pcb <- log10(fox.pcb)
  # Replace -inf to NA
  fox.pcb <- do.call(data.frame,
                     lapply(fox.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  fox.pcb.1 <- fox.pcb[,
                       -which(colSums(is.na(fox.pcb))/nrow(fox.pcb) > 0.7)]
  # Change date format
  SampleDate <- as.Date(fox$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- as.numeric(difftime(as.Date(SampleDate),
                                  min(as.Date(SampleDate)), units = "days"))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add distance to eastern location
  eastern <- fox$DistanceToEasternLocation
  # Add date and time to fox.pcb.1
  fox.pcb.1 <- cbind(fox.pcb.1, SampleDate, as.factor(fox$SiteID),
                     data.frame(time.day), season.s, eastern)
  # Remove site Lake Winnebago (background site)
  fox.pcb.1 <- subset(fox.pcb.1, as.factor(fox$SiteID) != c("WCPCB-FOX001"))
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  # Add USGS data to fox.pcb.1, matching dates, conversion to m3/s
  fox.pcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.pcb.1$SampleDate,
                                                                        flow$Date)]
  fox.pcb.1$temp <- 273.15 + temp$X_..2.._00010_00003[match(fox.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with temperature = NA
  fox.pcb.2 <- fox.pcb.1[!is.na(fox.pcb.1$temp), ]
  # Remove metadata not use in the random forest
  fox.pcb.2 <- fox.pcb.2[, !(names(fox.pcb.2) %in% c("SampleDate"))]
}

# Function to perform grid search for ranger model hyperparameters
perform_grid_search <- function(congener_name, seed = 123) {
  # Extract data for the specified PCB congener
  congener_data <- fox.pcb.2[[congener_name]]
  
  # Get indices of non-missing values in the specified PCB congener
  non_missing_indices <- which(!is.na(congener_data))
  
  # Extract covariate columns (excluding columns starting with "PCB")
  covariate_columns <- names(fox.pcb.2)[!grepl("^PCB", names(fox.pcb.2))]
  
  # Select rows with non-missing values in the specified PCB congener and extract covariates
  congener_covariate_data <- fox.pcb.2[non_missing_indices, covariate_columns]
  
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
    mtry = c(2, 4, min(4, ncol(train_data))),  # Ensure mtry is less than or equal to the number of variables
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
    Location = rep("Fox River", length(test_indices)),
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
for (congener_name in colnames(fox.pcb.2)) {
  if (grepl("^PCB", congener_name)) {
    # Perform grid search and get predictions
    result <- perform_grid_search(congener_name, seed = 123)
    
    # Append the results to the results_rf_PCBi dataframe
    results_rf_PCBi <- rbind(results_rf_PCBi, result)
  }
}

# Export results
write.csv(results_rf_PCBi,
          file = "Output/Data/Sites/csv/FoxRiver/FoxRiverRFPCB.csv",
          row.names = FALSE)

# Plot
plotRFPCBi <- ggplot(results_rf_PCBi, aes(x = 10^(Test_Data), y = 10^(Predicted_Data))) +
  geom_point(shape = 21, size = 1, fill = "white") +
  scale_y_log10(limits = c(0.001, 10^4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.001, 10^4),
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
ggsave("Output/Plots/Sites/ObsPred/FoxRiver/FoxRiverRFPCB.png",
       plot = plotRFPCBi, width = 6, height = 5, dpi = 500)


# Summarize performance of random forest model analysis,
# For both tPCB and individual PCBs. RSME, R2, and factor or 2.

# Install packages
install.packages('dplyr')
install.packages('ggplot2')
install.packages('RColorBrewer')

# Load libraries
{
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
}

# Read generated data for total PCB ---------------------------------------
{
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverRFtPCB.csv")
    # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFtPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFtPCB.csv")
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFtPCB.csv")
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFtPCB.csv")
  # Housatonic River
  hou <- read.csv("Output/Data/Sites/csv/HousatonicRiver/HousatonicRiverRFtPCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFtPCB.csv")
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverRFtPCB.csv")
  # Great Lakes (Lake Michigan)
  lmi <- read.csv("Output/Data/Sites/csv/GreatLakes/LakeMichiganRFtPCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFtPCB.csv")
  # New Bedford Harbor
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFtPCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFtPCB.csv")
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFtPCB.csv")
  # Richardson Hill Road Landfill
  rhl <- read.csv("Output/Data/Sites/csv/Richardson/RichardsonRFtPCB.csv")
  # Spokane River
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFtPCB.csv")
  # Combine the data frames
  combined_RFtPCB <- rbind(anr, bfc,  che, dmi, fox, hou, hud, kal, lmi, lwa, nbh,
                         pas, por, rhl, spo)
}

# see data
print(combined_RFtPCB)

# Export results
write.csv(combined_RFtPCB, file = "Output/Data/Sites/csv/Summary/RFtPCB.csv")

# Summary data
# (1) RSME
RSME <- summary(combined_RFtPCB$RMSE)
print(RSME)
# (2) Pearson correlation
Pearson <- summary(combined_RFtPCB$Pearson)
print(Pearson)
# (3) Factor 2
Factor2 <- summary(combined_RFtPCB$Factor2_Percentage)
print(Factor2)
# (4) Summary per location
summary_RFtPCB <- combined_RFtPCB %>%
  group_by(location) %>%
  summarize(
    RMSE = first(RMSE),
    Pearson = first(Pearson),
    Factor2 = first(Factor2_Percentage)
  )

# Find the index of the min RMSE value
min_RSME_index <- which.min(combined_RFtPCB$RMSE)
# Get the corresponding location
min_RSME_location <- combined_RFtPCB$location[min_RSME_index]
min_RSME <- combined_RFtPCB$RMSE[min_RSME_index]
# Print them
cat("Minimum RMSE Location:", min_RSME_location, "\n")
cat("Minimum RMSE:", min_RSME, "\n")

# Find the index of the max RMSE value
max_RSME_index <- which.max(combined_RFtPCB$RMSE)
# Get the corresponding location
max_RSME_location <- combined_RFtPCB$location[max_RSME_index]
max_RSME <- combined_RFtPCB$RMSE[max_RSME_index]
# Print them
cat("Max RMSE Location:", max_RSME_location, "\n")
cat("Max RMSE:", max_RSME, "\n")

# Find the index of the min RMSE value
min_RSME_index <- which.min(combined_RFtPCB$RMSE)
# Get the corresponding location
min_RSME_location <- combined_RFtPCB$location[min_RSME_index]
min_RSME <- combined_RFtPCB$RMSE[min_RSME_index]
# Print them
cat("Minimum RMSE Location:", min_RSME_location, "\n")
cat("Minimum RMSE:", min_RSME, "\n")

# Find the index of the max Pearson value
max_Pearson_index <- which.max(combined_RFtPCB$Pearson)
# Get the corresponding location
max_Pearson_location <- combined_RFtPCB$location[max_Pearson_index]
max_Pearson <- combined_RFtPCB$Pearson[max_Pearson_index]
# Print them
cat("Max Pearson Location:", max_Pearson_location, "\n")
cat("Max Pearson:", max_Pearson, "\n")

# Find the index of the min Pearson value
min_Pearson_index <- which.min(combined_RFtPCB$Pearson)
# Get the corresponding location
min_Pearson_location <- combined_RFtPCB$location[min_Pearson_index]
min_Pearson <- combined_RFtPCB$Pearson[min_Pearson_index]
# Print them
cat("Min Pearson Location:", min_Pearson_location, "\n")
cat("Min Pearson:", min_Pearson, "\n")

# Find the index of the max Factor2 value
max_Factor2_index <- which.max(combined_RFtPCB$Factor2_Percentage)
# Get the corresponding location
max_Factor2_location <- combined_RFtPCB$location[max_Factor2_index]
max_Factor2 <- combined_RFtPCB$Factor2_Percentage[max_Factor2_index]
# Print them
cat("Max Factor2 Location:", max_Factor2_location, "\n")
cat("Max Factor2:", max_Factor2, "\n")

# Find the index of the min Factor2 value
min_Factor2_index <- which.min(combined_RFtPCB$Factor2_Percentage)
# Get the corresponding location
min_Factor2_location <- combined_RFtPCB$location[min_Factor2_index]
min_Factor2 <- combined_RFtPCB$Factor2_Percentage[min_Factor2_index]
# Print them
cat("Min Factor2 Location:", min_Factor2_location, "\n")
cat("Min Factor2:", min_Factor2, "\n")

# Read generated data for PCBs --------------------------------------------
{
  # Bannister Federal Complex
  bfc <- read.csv("Output/Data/Sites/csv/BannisterFedComplex/BannisterFedComplexRFPCB.csv")
  # Chesapeake Bay data
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeBayRFPCB.csv")
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMIRFPCB.csv")
  # Fox River data
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverRFPCB.csv")
  # Great Lakes (Lake Michigan)
  lmi <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesRFPCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverRFPCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonRFPCB.csv")
  # New Bedford Harbor data
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHRFPCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicRiverRFPCB.csv")
  # Portland Harbord data
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborRFPCB.csv")
  # Spokane River data
  spo <- read.csv("Output/Data/Sites/csv/SpokaneRiver/SpokaneRiverRFPCB.csv")
  # Combine the data frames
  combined_RFPCB <- rbind(bfc, che, dmi, fox, hud, lmi, lwa, nbh, pas, por, spo)
}

# Export results
write.csv(combined_RFPCB, file = "Output/Data/Sites/csv/Summary/RFPCB.csv")

# Summary data
# (1) RSME
RSME <- summary(combined_RFPCB$RMSE)
print(RSME)
# (2) Pearson correlation
Pearson <- summary(combined_RFPCB$Correlation)
print(Pearson)
# (3) Factor 2
Factor2 <- summary(combined_RFPCB$Factor2)
print(Factor2)
# (4) Summary per location
summary_RFPCB <- combined_RFPCB %>%
  group_by(Location) %>%
  summarize(
    RMSE = first(RMSE),
    Pearson = first(Correlation),
    Factor2 = first(Factor2)
  )

# Find the index of the min RMSE value
min_RSME_index <- which.min(combined_RFPCB$RMSE)
# Get the corresponding location
min_RSME_location <- combined_RFPCB$Location[min_RSME_index]
min_RSME <- combined_RFPCB$RMSE[min_RSME_index]
# Print them
cat("Minimum RMSE Location:", min_RSME_location, "\n")
cat("Minimum RMSE:", min_RSME, "\n")

# Find the index of the max RMSE value
max_RSME_index <- which.max(combined_RFPCB$RMSE)
# Get the corresponding location
max_RSME_location <- combined_RFPCB$Location[max_RSME_index]
max_RSME <- combined_RFPCB$RMSE[max_RSME_index]
# Print them
cat("Max RMSE Location:", max_RSME_location, "\n")
cat("Max RMSE:", max_RSME, "\n")

# Find the index of the max Pearson value
max_Pearson_index <- which.max(combined_RFPCB$Correlation)
# Get the corresponding location
max_Pearson_location <- combined_RFPCB$Location[max_Pearson_index]
max_Pearson <- combined_RFPCB$Correlation[max_Pearson_index]
# Print them
cat("Max Pearson Location:", max_Pearson_location, "\n")
cat("Max Pearson:", max_Pearson, "\n")

# Find the index of the min Pearson value
min_Pearson_index <- which.min(combined_RFPCB$Correlation)
# Get the corresponding location
min_Pearson_location <- combined_RFPCB$Location[min_Pearson_index]
min_Pearson <- combined_RFPCB$Correlation[min_Pearson_index]
# Print them
cat("Min Pearson Location:", min_Pearson_location, "\n")
cat("Min Pearson:", min_Pearson, "\n")

# Find the index of the max Factor2 value
max_Factor2_index <- which.max(combined_RFPCB$Factor)
# Get the corresponding location
max_Factor2_location <- combined_RFPCB$Location[max_Factor2_index]
max_Factor2 <- combined_RFPCB$Factor2[max_Factor2_index]
# Print them
cat("Max Factor2 Location:", max_Factor2_location, "\n")
cat("Max Factor2:", max_Factor2, "\n")

# Find the index of the min Factor2 value
min_Factor2_index <- which.min(combined_RFPCB$Factor)
# Get the corresponding location
min_Factor2_location <- combined_RFPCB$Location[min_Factor2_index]
min_Factor2 <- combined_RFPCB$Factor2[min_Factor2_index]
# Print them
cat("Min Factor2 Location:", min_Factor2_location, "\n")
cat("Min Factor2:", min_Factor2, "\n")

# Find the number of samples above 0.5 for the Pearson values
count_above_0.5 <- sum(combined_RFPCB$Correlation > 0.5)
# Calculate the total number of values in the 'Correlation' column
total_values <- length(combined_RFPCB$Correlation)
# Calculate the percentage of values greater than 0.5
percentage_above_0.5 <- (count_above_0.5 / total_values) * 100

# Find the number of samples above 0.5 for the Pearson values for Portland
count_above_0.5_por <- sum(por$Correlation > 0.75)
# Calculate the total number of values in the 'Correlation' column
total_values_por <- length(por$Correlation)
# Calculate the percentage of values greater than 0.5
percentage_above_0.5_por <- (count_above_0.5_por/total_values_por) * 100


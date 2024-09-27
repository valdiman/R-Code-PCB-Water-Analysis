# Code to combine all Tributary River and Open Lake data from LMMB

# Install packages
install.packages("dplyr")

# Load libraries
{
  library(dplyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  lmmb.1 <- read.csv("Data/LMMB/Tributaries/Tributaries.csv")
  lmmb.2 <- read.csv("Data/LMMB/OpenLake/OpenLake.csv")
}

# Combine all the dataset in one.
merged_lmmb <- rbind(lmmb.1, lmmb.2)

# Delete the first column
merged_lmmb <- merged_lmmb[, -1]

# Create a data frame with unique Latitude and Longitude combinations
unique_combinations <- merged_lmmb %>%
  distinct(Latitude, Longitude)

# Function to generate SiteID based on LATITUDE and LONGITUDE
generate_SiteID <- function(lat, lon) {
  index <- which(unique_combinations$Latitude == lat & unique_combinations$Longitude == lon)
  SiteID <- paste("WCPCB-LMM", sprintf("%03d", index), sep = "")
  return(SiteID)
}

# Add a SiteID column based on LATITUDE and LONGITUDE
merged_lmmb <- merged_lmmb %>%
  mutate(SiteID = mapply(generate_SiteID, Latitude, Longitude))

# Create a new column "SampleDate_yyyymmdd" with the date in yyyymmdd format
merged_lmmb$SampleDate_yyyymmdd <- format(as.Date(merged_lmmb$SampleDate,
                                                 format = "%m/%d/%y"), format = "%Y%m%d")

# Create SampleID by concatenating SiteID and SampleDate_yyyymmdd
merged_lmmb$SampleID <- paste(merged_lmmb$SiteID,
                             merged_lmmb$SampleDate_yyyymmdd, sep = "-")

# Create a new column "SampleCount" to count samples with the same SampleID
merged_lmmb <- merged_lmmb %>%
  group_by(SampleID) %>%
  mutate(SampleCount = row_number())

# Add a dot and SampleCount to SampleID if SampleCount is greater than 1
for (i in unique(merged_lmmb$SampleID)) {
  subset_df <- merged_lmmb[merged_lmmb$SampleID == i, ]
  num_samples <- nrow(subset_df)
  if (num_samples > 1) {
    for (j in 1:num_samples) {
      subset_df$SampleID[j] <- paste(subset_df$SampleID[j], ".", j, sep = "")
    }
    merged_lmmb[merged_lmmb$SampleID == i, ] <- subset_df
  } else {
    merged_lmmb[merged_lmmb$SampleID == i, "SampleID"] <- paste(i, ".1", sep = "")
  }
}

# Remove the SampleCount and SampleDate_yyyymmdd columns if no longer needed
merged_lmmb$SampleCount <- NULL
merged_lmmb$SampleDate_yyyymmdd <- NULL

# Remove samples with tPCB = 0
merged_lmmb <- merged_lmmb[!merged_lmmb$tPCB == 0, ]

# Export results
write.csv(merged_lmmb, file = "Data/LMMB/lmmb.csv")


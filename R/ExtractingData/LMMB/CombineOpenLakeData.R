# Code to combine all opebutary River data from LMMB

# Install packages
install.packages("dplyr")

# Load libraries
{
  library(dplyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  ope.1 <- read.csv("Data/LMMB/OpenLake/1994SOL.csv")
  ope.2 <- read.csv("Data/LMMB/OpenLake/1994SuOL.csv")
  ope.3 <- read.csv("Data/LMMB/OpenLake/1995SOL.csv")
  ope.4 <- read.csv("Data/LMMB/OpenLake/1995SuOL.csv")
}

# Combine all the dataset in one.
merged_ope <- rbind(ope.1, ope.2, ope.3, ope.4)

# Delete the first column
merged_ope <- merged_ope[, -1]

# Transform SampleDate into date format
merged_ope$SampleDate <- as.Date(merged_ope$SampleDate, format = "%Y/%m/%d")

# Convert it back to character with the desired format
merged_ope$SampleDate <- format(merged_ope$SampleDate, format = "%m/%d/%y")

# Rename SAMPLE_ID to SiteName and change the name of the SAMPLE_ID
names(merged_ope)[names(merged_ope) == "SAMPLE_ID"] <- "SiteName"
merged_ope$SiteName <- sub("^[^-]*-", "", merged_ope$SiteName)

# Names and values for the new columns
new_col_names <- c("SampleID", "EPARegion", "StateSampled", "LocationName",
                   "SiteID")
new_col_values <- c("SampleIDValue", "R5", NA, "Lake Michigan Mass Balance",
                    "SiteIDValue")

# Add new columns at the beginning (from column 1)
merged_ope <- cbind(setNames(data.frame(matrix(NA, nrow = nrow(merged_ope),
                                               ncol = length(new_col_names))),
                             new_col_names), merged_ope)

# Fill the new columns with values
for (i in 1:length(new_col_names)) {
  col_name <- new_col_names[i]
  col_value <- new_col_values[i]
  merged_ope[, col_name] <- col_value
}

# Organize column order
desired_column_order <- c(
  "SampleID",
  "EPARegion",
  "StateSampled",
  "LocationName",
  "SiteName",
  "SiteID",
  "SampleDate",
  "Latitude",
  "Longitude",
  "Units"
)

merged_ope <- merged_ope %>%
  select(desired_column_order, everything())

# Fill StatedSampled column
state_mapping <- list(
  IL = c(1, "1-7m", "1-8m", 5, "5-17m", "5-25m", "5-3m", "5-8m", "MB9", "MB9-20m", "MB9-7m"),
  WI = c("17-10m", "17-50m", "17-74m", "180-13m", "180-34m", "180-50m", "180-58m", "180-5m",
         "180-7m", "24-25m", "240-30m", "240-37m", "240-7m", "280-10m", "280-42m", "280-42m-285L", "280-42m-380L",
         "280-54m", "280-5m", "280-62m", "31-13m", "31-6m", "40M-10m", "40M-110m", "40M-113m", "40M-116m",
         "40M-116m", "40M-15m", "40M-25m", "40M-5m", "40M-82m", "45-14m", "45-15m", "45-5m", "47M", "47M-110m",
         "47M-110m", "47M-112m", "47M-119m", "47M-30m", "47M-33m", "47M-5m",  "47M-8m", "47M-96m", "GB100M",
         "GB100M-10m", "GB100M-34m", "GB100M-44m", "GB100M-45m", "GB100M-50m", "GB100M-8m", "GB17-11m",
         "GB17-11m-190L-Full", "GB17-11m-190L-Half", "GB17-11m-285L-Half", "GB17-17m", "GB17-7m",
         "GM100M-10m", "GM100M-34m", "GM100M-44m", "GM100M-45m", "GM100M-50m", "GM100M-8m", "MB21-24m",
         "MB21-5m", "MB25-15m", "MB25-22m", "MB25-5.5m", "MB25-6m", "MB38-10m", "MB38-11.5m", "MB38-17m",
         "MB38-5.5m")
)

# Create a function to map SiteName to StateSampled
map_site_to_state <- function(site_name) {
  for (state in c("IN", "IL", "WI")) {  # Change the order as needed
    if (site_name %in% state_mapping[[state]]) {
      return(state)
    }
  }
  return("MI")  # Default to "MI" if no match is found
}

# Apply the mapping function to create the StateSampled column
merged_ope$StateSampled <- sapply(merged_ope$SiteName, map_site_to_state)

# Define the values for the new columns
PhaseMeasuredValue <- "SurfaceWater"
EPAMethodValue <- "M1668"
AroclorCongenerValue <- "Congener"

# Insert the new columns at position 11
merged_ope <- merged_ope %>%
  mutate(PhaseMeasured = PhaseMeasuredValue, EPAMethod = EPAMethodValue,
         AroclorCongener = AroclorCongenerValue) %>%
  select(1:10, PhaseMeasured, EPAMethod, AroclorCongener, everything())

# Define the column names
column_names <- c("A1016", "A1221", "A1232", "A1242", "A1248", "A1254", "A1260")

# Initialize these columns with "NA"
merged_ope[, column_names] <- NA

# Insert the columns at position 118
merged_ope <- merged_ope %>%
  select(1:117, all_of(column_names), everything())

# Export results
write.csv(merged_ope, file = "Data/LMMB/OpenLake/OpenLake.csv")

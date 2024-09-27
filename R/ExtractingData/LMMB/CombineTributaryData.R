# Code to combine all Tributary River data from LMMB

# Install packages
install.packages("dplyr")

# Load libraries
{
  library(dplyr)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
{
  tri.1 <- read.csv("Data/LMMB/Tributaries/1994T.csv")
  tri.2 <- read.csv("Data/LMMB/Tributaries/1995T.csv")
}

# Combine all the dataset in one.
merged_tri <- rbind(tri.1, tri.2)

# Delete the first column
merged_tri <- merged_tri[, -1]

# Transform SampleDate into date format
merged_tri$SampleDate <- as.Date(merged_tri$SampleDate, format = "%Y/%m/%d")

# Convert it back to character with the desired format
merged_tri$SampleDate <- format(merged_tri$SampleDate, format = "%m/%d/%y")

# Names and values for the new columns
new_col_names <- c("SampleID", "EPARegion", "StateSampled", "LocationName",
                   "SiteID")
new_col_values <- c("SampleIDValue", "R5", NA, "Lake Michigan Mass Balance",
                    "SiteIDValue")

# Add new columns at the beginning (from column 1)
merged_tri <- cbind(setNames(data.frame(matrix(NA, nrow = nrow(merged_tri),
                                               ncol = length(new_col_names))),
                             new_col_names), merged_tri)

# Fill the new columns with values
for (i in 1:length(new_col_names)) {
  col_name <- new_col_names[i]
  col_value <- new_col_values[i]
  merged_tri[, col_name] <- col_value
}

# Create SiteName column
# Define the code_to_name mapping
code_to_name <- c(
  "TFOXRB" = "Tributary Fox River", #WI
  "TGRANH" = "Tributary Grand River", #MI
  "TIHCAE" = "Tributary Indiana Harbor Canal", #IN
  "TKALAG" = "Tributary Kalamazoo River", #MI
  "TMANIS" = "Tributary Manistee River", #MI
  "TMANIK" = "Tributary Manistique River", #MI
  "TMENOA" = "Tributary Menominee River", #WI
  "TMILWD" = "Tributary Milwaukee River", #WI
  "TMUSKI" = "Tributary Muskegon River", #MI
  "TPEREJ" = "Tributary Pere Marquette River", #MI
  "TSHEBC" = "Tributary Sheboygan River", #WI
  "TSTJOF" = "Tributary St. Joseph River" #MI
)

# Extract the relevant part of SAMPLE_ID and create the SiteName columns
merged_tri$SiteName <- sapply(merged_tri$SAMPLE_ID, function(sample_id) {
  match_key <- gsub("[0-9]", "", sample_id)
  return(code_to_name[match_key])
})

# Define the code_to_name mapping (site names to state abbreviations)
code_to_state <- c(
  "Tributary Fox River" = "WI",
  "Tributary Grand River" = "MI",
  "Tributary Indiana Harbor Canal" = "IN",
  "Tributary Kalamazoo River" = "MI",
  "Tributary Manistee River" = "MI",
  "Tributary Manistique River" = "MI",
  "Tributary Menominee River" = "WI",
  "Tributary Milwaukee River" = "WI",
  "Tributary Muskegon River" = "MI",
  "Tributary Pere Marquette River" = "MI",
  "Tributary Sheboygan River" = "WI",
  "Tributary St. Joseph River" = "MI"
)

# Create the StateSampled column based on the SiteName column
merged_tri$StateSampled <- sapply(merged_tri$SiteName, function(site_name) {
  return(code_to_state[site_name])
})

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

merged_tri <- merged_tri %>%
  select(desired_column_order, everything())

# Remove SAMPLE_ID column
merged_tri$SAMPLE_ID <- NULL

# Define the values for the new columns
PhaseMeasuredValue <- "SurfaceWater"
EPAMethodValue <- "M1668"
AroclorCongenerValue <- "Congener"

# Insert the new columns at position 11
merged_tri <- merged_tri %>%
  mutate(PhaseMeasured = PhaseMeasuredValue, EPAMethod = EPAMethodValue,
         AroclorCongener = AroclorCongenerValue) %>%
  select(1:10, PhaseMeasured, EPAMethod, AroclorCongener, everything())

# Define the column names
column_names <- c("A1016", "A1221", "A1232", "A1242", "A1248", "A1254", "A1260")

# Initialize these columns with "NA"
merged_tri[, column_names] <- NA

# Insert the columns at position 118
merged_tri <- merged_tri %>%
  select(1:117, all_of(column_names), everything())

# Export results
write.csv(merged_tri, file = "Data/LMMB/Tributaries/Tributaries.csv")



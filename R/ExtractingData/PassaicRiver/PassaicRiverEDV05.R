
# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
}

# Read the CSV file into a data frame
PS_data <- read_csv("Data/PassaicRiver/2012 CPG CWCM Sampling - Round 5.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
PS_data <- PS_data %>% arrange(SAMPLE_NAME, SAMPLE_DATE)

# Filter and keep only the columns you need
PS_data <- PS_data %>%
  select(SAMPLE_NAME, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE, ANALYTIC_METHOD,
         Y_COORD, X_COORD, REPORT_RESULT_UNIT, CHEMICAL_NAME, RESULT_NUMERIC)

# Filter rows based on the ANALYTIC_METHOD condition
PS_data <- PS_data %>%
  filter(ANALYTIC_METHOD == "E1668A")

# Modify the CHEMICAL_NAME column to extract "PCB X" where X is the number
PS_data <- PS_data %>%
  mutate(CHEMICAL_NAME = sub(".*\\(PCB (\\d+)\\).*", "PCB \\1", CHEMICAL_NAME))

# Filter rows where CHEMICAL_NAME contains "PCB" and not "Total"
PS_data <- PS_data %>%
  filter(grepl("PCB", CHEMICAL_NAME,
               fixed = TRUE) & !grepl("Total", CHEMICAL_NAME, fixed = TRUE))

# Create a new data frame with transposed values
transposed_data <- PS_data %>%
  pivot_wider(
    id_cols = c(SAMPLE_NAME, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE,
                ANALYTIC_METHOD, Y_COORD, X_COORD, REPORT_RESULT_UNIT),
    names_from = CHEMICAL_NAME,
    values_from = RESULT_NUMERIC
  )

# Read the JSON file with new congener list from code NewPCBList.R
pcb_groups <- read_json("Data/pcb_groups.json")

# Create an empty data frame to store the grouped data
grouped_data <- PS_data %>%
  distinct(SAMPLE_NAME, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE,
           ANALYTIC_METHOD, Y_COORD, X_COORD, REPORT_RESULT_UNIT)

# Remove spaces and special characters from column names in transposed_data
colnames(transposed_data) <- gsub("[^[:alnum:]]", "", colnames(transposed_data))

# Iterate through the pcb_groups list and sum columns
for (group_name in names(pcb_groups)) {
  group_columns <- pcb_groups[[group_name]]
  group_columns <- gsub("[^[:alnum:]]", "", group_columns)  # Clean group column names
  grouped_data[group_name] <- rowSums(transposed_data[group_columns], na.rm = TRUE)
}

# Remove SAMPLE_TYPE_CODE and MATRIX_CODE columns from grouped_data
grouped_data <- grouped_data %>%
  select(-SAMPLE_TYPE_CODE, -MATRIX_CODE, -ANALYTIC_METHOD)

# Multiply PCB columns by 1000
pcb_columns <- names(grouped_data)[grep("^PCB\\d+", names(grouped_data))]
grouped_data[, pcb_columns] <- grouped_data[, pcb_columns] * 1000

# Update REPORT_RESULT_UNIT to "pg/l"
grouped_data$REPORT_RESULT_UNIT <- "pg/l"

# Create a new column named "tPCB" that sums columns 6 to 109
grouped_data <- grouped_data %>%
  mutate(tPCB = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Replace 0s with NA for columns 6 to 109
grouped_data <- grouped_data %>%
  mutate_at(vars(starts_with("PCB")), ~ ifelse(. == 0, NA, .))

# Change SAMPLE_NAME
grouped_data <- grouped_data %>%
  mutate(SAMPLE_NAME = substr(SAMPLE_NAME, 1, nchar(SAMPLE_NAME) - 3))

# Remove SAMPLE_NAME with no name (NA)
grouped_data <- grouped_data %>%
  filter(!is.na(SAMPLE_NAME))

# Change the name of columns to be consistent
colnames(grouped_data)[colnames(grouped_data) == "SAMPLE_DATE"] <- "SampleDate"
colnames(grouped_data)[colnames(grouped_data) == "Y_COORD"] <- "Latitude"
colnames(grouped_data)[colnames(grouped_data) == "X_COORD"] <- "Longitude"
colnames(grouped_data)[colnames(grouped_data) == "REPORT_RESULT_UNIT"] <- "Units"

# Convert SampleDate to the desired format "m/d/yy"
grouped_data$SampleDate <- as.Date(grouped_data$SampleDate, format = "%m/%d/%y")

# Convert it back to character with the desired format
grouped_data$SampleDate <- format(grouped_data$SampleDate, format = "%m/%d/%y")

# Export results
write.csv(grouped_data, file = "Data/PassaicRiver/pass05.csv")


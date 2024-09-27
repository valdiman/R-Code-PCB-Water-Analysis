
# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("jsonlite")
install.packages("stringr")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(stringr)
}

# Read the CSV file into a data frame
PS_data <- read_csv("Data/PassaicRiver/2007-08 USEPA-MPI.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
PS_data <- PS_data %>% arrange(sys_sample_code, sample_date)

# Filter and keep only the columns you need
PS_data <- PS_data %>%
  select(sys_sample_code, sample_date, analytic_method, latitude,
         longitude, chemical_name2, result_unit, result_numeric)

# Filter rows based on the ANALYTIC_METHOD condition & just water samples (i.e., pg/l)
# and initial analysis
PS_data <- PS_data %>%
  filter(analytic_method == "AXYS MLA013COEXTRACT", result_unit == "pg/l")

# Modify the CHEMICAL_NAME column to extract "PCB X" where X is the number
PS_data <- PS_data %>%
  mutate(chemical_name2 = sub(".*\\(PCB (\\d+)\\).*", "PCB \\1", chemical_name2))

# Filter rows where CHEMICAL_NAME contains "PCB" and not "Total"
# and "Polychlorinated Biphenyl (PCB)"
PS_data <- PS_data %>%
  filter(grepl("PCB", chemical_name2, fixed = TRUE) &
           !grepl("Total", chemical_name2, fixed = TRUE) &
           !grepl("Polychlorinated Biphenyl", chemical_name2))

# Create a new data frame with transposed values
transposed_data <- PS_data %>%
  pivot_wider(
    id_cols = c(sys_sample_code, sample_date, analytic_method, latitude,
                longitude),
    names_from = chemical_name2,
    values_from = result_numeric,
    values_fn = list  # Treat duplicate values as lists
  )

# Read the JSON file with new congener list from code NewPCBList.R
pcb_groups <- read_json("Data/pcb_groups.json")

# Create an empty data frame to store the grouped data
grouped_data <- PS_data %>%
  distinct(sys_sample_code, sample_date, analytic_method, latitude,
           longitude, result_unit)

# Remove spaces and special characters from column names in transposed_data
colnames(transposed_data) <- gsub("[^[:alnum:]]", "", colnames(transposed_data))

# Function to safely coerce to numeric, replacing non-numeric values with NAs
safe_as_numeric <- function(x) {
  as.numeric(as.character(x))
}

# Iterate through the pcb_groups list and sum columns
for (group_name in names(pcb_groups)) {
  group_columns <- pcb_groups[[group_name]]
  group_columns <- gsub("[^[:alnum:]]", "", group_columns)  # Clean group column names
  
  # Coerce group_columns to numeric
  transposed_data[group_columns] <- lapply(transposed_data[group_columns], safe_as_numeric)
  
  # Calculate row sums, replacing NA with 0
  grouped_data[group_name] <- rowSums(transposed_data[group_columns], na.rm = TRUE)
  
  # Clean up: replace 0 with NA
  grouped_data[group_name][grouped_data[group_name] == 0] <- NA
}

# Create a new column named "tPCB" that sums columns 6 to 109
grouped_data <- grouped_data %>%
  mutate(tPCB = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Replace 0s with NA for columns 6 to 109
grouped_data <- grouped_data %>%
  mutate_at(vars(starts_with("PCB")), ~ ifelse(. == 0, NA, .))

# Concatenate new_column and last_three_digits
grouped_data <- grouped_data %>%
  mutate(new_column = gsub("^(.*[a-zA-Z]).*$", "\\1", sys_sample_code),
         last_three_digits = str_extract(sys_sample_code, "\\d{3}$"),
         concatenated_column = paste(new_column, last_three_digits, sep = "-"))

# Move concatenated_column to sys_sample_code
# and eliminate new_column and last_three_digits
# Concatenate new_column and last_three_digits
grouped_data <- grouped_data %>%
  mutate(sys_sample_code = paste(new_column, last_three_digits, sep = "-"))

# Remove new_column, last_three_digits, and concatenated_column
grouped_data <- grouped_data %>%
  select(-new_column, -last_three_digits, -concatenated_column, -analytic_method)

# Remove SAMPLE_NAME with no name (NA)
grouped_data <- grouped_data %>%
  filter(!is.na(sys_sample_code))

# Change the name of columns to be consistent
colnames(grouped_data)[colnames(grouped_data) == "latitude"] <- "Latitude"
colnames(grouped_data)[colnames(grouped_data) == "longitude"] <- "Longitude"
colnames(grouped_data)[colnames(grouped_data) == "sample_date"] <- "SampleDate"
colnames(grouped_data)[colnames(grouped_data) == "sys_sample_code"] <- "SAMPLE_NAME"
colnames(grouped_data)[colnames(grouped_data) == "result_unit"] <- "Units"

# Convert SampleDate to the desired format "m/d/yy"
grouped_data$SampleDate <- as.Date(grouped_data$SampleDate, format = "%m/%d/%y")

# Convert it back to character with the desired format
grouped_data$SampleDate <- format(grouped_data$SampleDate, format = "%m/%d/%y")

# Export results
write.csv(grouped_data, file = "Data/PassaicRiver/pass10.csv")


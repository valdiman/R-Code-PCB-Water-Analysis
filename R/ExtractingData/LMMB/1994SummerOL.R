# https://www.epa.gov/greatlakes/lake-michigan-mass-balance-results-and-publications
# Data: https://cdxapps.epa.gov/cdx-glenda/action/querytool/querySystem
# Welcome to the Great Lakes Environmental Database System (GLENDA)
# GLENDA Query System
# Selected Query: lake Michigan Mass Balance Results
# Project Code: BALP: LAKE PCB
# Year: 1994
# Season: Summer
# Medium: surface water

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("jsonlite")
install.packages("lubridate")
install.packages("strings")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(lubridate)
  library(stringr)
}

# Read the CSV file into a data frame
LMMB_data <- read_csv("Data/LMMB/OpenLake/1994SummerOL.csv")

# Function to clean and convert VALUE_ columns to numeric
clean_and_convert_value <- function(x) {
  # Remove asterisks and any non-numeric characters except '.'
  cleaned_value <- gsub("[^0-9.*-]", "", x)
  # Replace asterisk with an empty string
  cleaned_value <- gsub("\\*", "", cleaned_value)
  # Handle empty cells by converting them to NA
  cleaned_value[cleaned_value == ""] <- NA
  # Convert to numeric
  as.numeric(cleaned_value)
}

# Identify the VALUE_ columns
value_columns <- grep("^VALUE_", names(LMMB_data), value = TRUE)

# Clean and convert each VALUE_ column
LMMB_data[value_columns] <- lapply(LMMB_data[value_columns],
                                   clean_and_convert_value)

# Identify the range of columns for chemical data
start_column <- which(names(LMMB_data) == "ANL_CODE_1")
pivot_columns <- names(LMMB_data)[start_column:length(LMMB_data)]

# Pivot the data while keeping metadata columns and maintaining the "ANAL_CODE" column
LMMB_data_long <- LMMB_data %>%
  pivot_longer(
    cols = all_of(pivot_columns),
    names_to = c(".value", "Chemical"),
    names_pattern = "^(ANL_CODE|ANALYTE|VALUE|UNITS|FRACTION|RESULT_REMARK)_([0-9]+)$"
  ) %>%
  filter(!is.na(VALUE))

# Select samples
LMMB_data_long <- LMMB_data_long %>%
  filter(FRACTION == "Filtrate", UNITS == "ng/l",
         QC_TYPE == "routine field sample",
         SAMPLE_TYPE == "Individual",
         grepl("^PCB", ANL_CODE))

# Remove the metadata column not necessaries
LMMB_data_long <- LMMB_data_long %>%
  select(-Row, -PROJECT, -PROJ_CODE, -YEAR, -MONTH, -SEASON,
         -CRUISE_ID, -VISIT_ID, -STATION_ID, -TIME_ZONE,
         -STN_DEPTH_M, -SAMPLE_DEPTH_M, -SAMPLE_TYPE, -QC_TYPE,
         -Chemical, -FRACTION, -RESULT_REMARK)

# Remove PCB congeners with 2 or more measurements (both measurements are the same)
LMMB_data_long <- LMMB_data_long %>%
  group_by(LATITUDE, LONGITUDE, SAMPLING_DATE, MEDIUM, SAMPLE_ID, ANL_CODE, ANALYTE) %>%
  mutate(VALUE = ifelse(row_number() == 1, VALUE, NA)) %>%
  ungroup() %>%
  na.omit()

# Change date format
LMMB_data_long$SAMPLING_DATE <- sub(" \\d+:\\d+", "",
                                    LMMB_data_long$SAMPLING_DATE)

# Function to change PCB names
transform_ANL_CODE <- function(anl_code) {
  # Replace "PCB-" with "PCB"
  anl_code <- gsub("^PCB-", "PCB", anl_code)
  
  # Replace "PCB" with "PCB"
  anl_code <- gsub("^PCB", "PCB", anl_code)
  
  # Extract numbers using regular expressions
  numbers <- str_extract_all(anl_code, "\\d+")
  
  if (length(numbers) > 0) {
    # Sort numbers in ascending order
    sorted_numbers <- sort(as.integer(unlist(numbers)))
    
    # Replace dots with plus signs and combine sorted numbers
    transformed_code <- paste("PCB", gsub("\\.", "+", paste(sorted_numbers, collapse = ".")), sep = "")
  } else {
    transformed_code <- anl_code
  }
  
  return(transformed_code)
}

# Apply function to change the PCB names
LMMB_data_long$Transformed_ANL_CODE <- sapply(LMMB_data_long$ANL_CODE,
                                              transform_ANL_CODE)

# Reshape the dataset
transposed_data <- LMMB_data_long %>%
  pivot_wider(
    id_cols = c(LATITUDE, LONGITUDE, SAMPLING_DATE, SAMPLE_ID, UNITS),
    names_from = Transformed_ANL_CODE,
    values_from = VALUE
  )

# Multiply PCB columns by 1000
pcb_columns <- names(transposed_data)[grep("^PCB\\d+", names(transposed_data))]
transposed_data[, pcb_columns] <- transposed_data[, pcb_columns] * 1000

# Multiply all values in the "PCB" column by 1000
transposed_data$PCB <- transposed_data$PCB * 1000

# Update RESULT_UNIT to "pg/l"
transposed_data$UNITS <- "pg/l"



write.csv(transposed_data, file = "Data/LMMB/OpenLake/1994SuOL.csv")

# Fix and create new PCB congener list
# List of column renaming rules
column_renaming_rules <- list(
  "PCB18" = "PCB18+30",
  "PCB44" = "PCB44+47+65",
  "PCB47+48" = "PCB48+59+62+75",
  "PCB91" = "PCB88+91",
  "PCB101" = "PCB90+101+113",
  "PCB118" = "PCB106+118",
  "PCB105+132+153" = "PCB132+153+161+168", # Issue with the original PCB coelution
  "PCB134" = "PCB134+143"
)

# Apply column renaming rules
for (old_name in names(column_renaming_rules)) {
  new_name <- column_renaming_rules[[old_name]]
  colnames(transposed_data)[colnames(transposed_data) == old_name] <- new_name
}

# Change NA to 0s
transposed_data[is.na(transposed_data)] <- 0

# Function to combine columns and remove original columns
combine_and_remove <- function(data, new_column, columns_to_combine) {
  data[[new_column]] <- rowSums(data[columns_to_combine], na.rm = TRUE)
  data <- data[, !colnames(data) %in% columns_to_combine]
  return(data)
}

# Specify the combinations and call the function for each
combinations <- list(
  list("PCB12+13", c("PCB12", "PCB13")),
  list("PCB16+32", c("PCB16", "PCB32")),
  list("PCB20+21+28+31+33+50+53", c("PCB21", "PCB28+31", "PCB33", "PCB53")),
  list("PCB24+27", c("PCB24", "PCB27")),
  list("PCB26+29", c("PCB26", "PCB29")),
  list("PCB40+41+64+71+72", c("PCB40", "PCB41+71", "PCB64")),
  list("PCB43+49+52+69+73", c("PCB43", "PCB49", "PCB52")),
  list("PCB45+51", c("PCB45", "PCB51")),
  list("PCB61+66+70+74+76+93+95+98+100+102", c("PCB66", "PCB70+76", "PCB74", "PCB95", "PCB100")),
  list("PCB77+85+110+111+115+116+117", c("PCB77+110", "PCB85")),
  list("PCB81+86+87+97+107+108+109+112+119+124+125", c("PCB81", "PCB87", "PCB97", "PCB107")),
  list("PCB82+135+144+151+154", c("PCB82", "PCB135+144", "PCB151")),
  list("PCB83+99", c("PCB83", "PCB99")),
  list("PCB114+122+131+133+142+146+165", c("PCB114+131", "PCB146")),
  list("PCB123+139+140+147+149", c("PCB123+149", "PCB124+147")), # issue with the original coelution.
  list("PCB128+162+166+167", c("PCB128", "PCB167")),
  list("PCB129+137+138+158+160+163+164+176+178", c("PCB129", "PCB137+176", "PCB138+163", "PCB158", "PCB178")),
  list("PCB156+157+172+197+200", c("PCB156", "PCB157+200", "PCB172+197")),
  list("PCB171+173+202", c("PCB171+202", "PCB173")),
  list("PCB180+193", c("PCB180", "PCB193")),
  list("PCB183+185", c("PCB183", "PCB185")),
  list("PCB194+205", c("PCB194", "PCB205")),
  list("PCB196+203", c("PCB196", "PCB203")),
  list("PCB198+199+201", c("PCB198", "PCB199", "PCB201"))
)

# Apply the function for each combination
for (combo in combinations) {
  transposed_data <- combine_and_remove(transposed_data, combo[[1]], combo[[2]])
}

# Especial column cases
{
  # Create new columns newPCB15 and newPCB17
  transposed_data$PCB15 <- transposed_data$'PCB15+17' * 0.5
  transposed_data$PCB17 <- transposed_data$'PCB15+17' * 0.5
  
  # Remove the original PCB15+17 and newPCB17 columns
  transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("PCB15+17")]
}

# Define the desired order of columns
desired_column_order <- c(
  "LATITUDE", "LONGITUDE", "SAMPLING_DATE", "SAMPLE_ID", "UNITS", "PCB", "PCB1", "PCB2", "PCB3",
  "PCB4+10", "PCB5+8", "PCB6", "PCB7+9", "PCB11", "PCB12+13",
  "PCB14", "PCB15", "PCB16+32", "PCB17", "PCB18+30", "PCB19", "PCB20+21+28+31+33+50+53",
  "PCB22", "PCB23", "PCB24+27", "PCB25", "PCB26+29", "PCB34", "PCB35", "PCB36", "PCB37+42",
  "PCB38", "PCB39", "PCB40+41+64+71+72", "PCB43+49+52+69+73", "PCB44+47+65", "PCB45+51", "PCB46",
  "PCB48+59+62+75", "PCB54", "PCB55", "PCB56+60", "PCB57", "PCB58", "PCB61+66+70+74+76+93+95+98+100+102",
  "PCB63", "PCB67", "PCB68", "PCB77+85+110+111+115+116+117", "PCB78", "PCB79", "PCB80",
  "PCB81+86+87+97+107+108+109+112+119+124+125", "PCB82+135+144+151+154", "PCB83+99", "PCB84+92",
  "PCB88+91", "PCB89", "PCB90+101+113", "PCB94", "PCB96", "PCB103", "PCB104", "PCB105",
  "PCB106+118", "PCB114+122+131+133+142+146+165", "PCB120", "PCB121", "PCB123+139+140+147+149",
  "PCB126", "PCB127", "PCB128+162+166+167", "PCB129+137+138+158+160+163+164+176+178",
  "PCB130", "PCB132+153+161+168", "PCB134+143", "PCB136", "PCB141", "PCB145", "PCB148", "PCB150",
  "PCB152", "PCB155", "PCB156+157+172+197+200", "PCB159", "PCB169", "PCB170+190", "PCB171+173+202",
  "PCB174", "PCB175", "PCB177", "PCB179", "PCB180+193", "PCB181", "PCB182+187", "PCB183+185",
  "PCB184", "PCB186", "PCB188", "PCB189", "PCB191", "PCB192", "PCB194+205", "PCB195+208",
  "PCB196+203", "PCB198+199+201", "PCB204", "PCB206", "PCB207", "PCB209"
)

# Check and add missing columns
missing_columns <- setdiff(desired_column_order, colnames(transposed_data))
transposed_data[missing_columns] <- NA

# Reorder the columns based on the desired order
transposed_data <- transposed_data[desired_column_order]

# Change NA to 0s
transposed_data[is.na(transposed_data)] <- 0

# Rename the "PCB" column to "tPCB"
colnames(transposed_data)[colnames(transposed_data) == "PCB"] <- "tPCB"

# Reorder the columns so that "tPCB" is at the end
transposed_data <- transposed_data %>%
  select(-tPCB, everything())

# Create a new column named "tPCB.2" that sums columns 5 to 109 to check original tPCB value
transposed_data <- transposed_data %>%
  mutate(tPCB.2 = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Remove tPCB.2
transposed_data <- transposed_data[, !colnames(transposed_data) %in% c("tPCB.2")]

# Change the name of columns to be consistent
colnames(transposed_data)[colnames(transposed_data) == "SAMPLING_DATE"] <- "SampleDate"
colnames(transposed_data)[colnames(transposed_data) == "LATITUDE"] <- "Latitude"
colnames(transposed_data)[colnames(transposed_data) == "LONGITUDE"] <- "Longitude"
colnames(transposed_data)[colnames(transposed_data) == "UNITS"] <- "Units"

# Export results
write.csv(transposed_data, file = "Data/LMMB/OpenLake/1994SuOL.csv")


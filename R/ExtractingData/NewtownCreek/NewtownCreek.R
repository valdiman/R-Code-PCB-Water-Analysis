# Data downloaded from: https://newtowncreek.com/resources/ri-report/#elf_l1_Lw

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")
install.packages("stringr")
install.packages("sp")
install.packages("sf")

# Load libraries
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(sp)
  library(sf)
}

# Read the CSV file into a data frame
NC_data <- read_csv("Data/NewtownCreek/Bi-B10-8_Analytical-Results.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
NC_data <- NC_data %>% arrange(sys_sample_code, sample_date)

# Filter and keep only the columns you need
NC_data <- NC_data %>%
  select(sys_sample_code, subfacility_code, sample_date, matrix_code, fraction,
         analytic_method, y_coord, x_coord, target_unit, chemical_name,
         result_value, dilution_factor)

# Filter rows based on the analytic_method and fraction condition
NC_data <- NC_data %>%
  filter(analytic_method == "E1668A",
         fraction == "D")

# Modify the CHEMICAL_NAME column to extract "PCB X" where X is the number
NC_data <- NC_data %>%
  mutate(chemical_name = sub(".*\\(PCB (\\d-)\\).*", "PCB \\1", chemical_name))

transform_chemical_name <- function(chem_code) {
  # Replace "PCB-" with "PCB"
  chem_code <- gsub("PCB-", "PCB", chem_code)
  
  # Extract numbers using regular expressions
  numbers <- str_extract_all(chem_code, "\\d+")
  
  if (length(numbers) > 0) {
    # Sort numbers in ascending order
    sorted_numbers <- sort(as.integer(unlist(numbers)))
    
    # Combine sorted numbers with "+"
    transformed_code <- paste("PCB", paste(sorted_numbers, collapse = "+"), sep = "")
  } else {
    transformed_code <- chem_code
  }
  
  return(transformed_code)
}

# Apply function to change the PCB names
NC_data$Transformed_chemical_name <- sapply(NC_data$chemical_name,
                                              transform_chemical_name)

# Create a new data frame with transposed values
transposed_data <- NC_data %>%
  pivot_wider(
    id_cols = c(sys_sample_code, subfacility_code, sample_date, matrix_code,
                fraction, analytic_method, y_coord, x_coord, target_unit),
    names_from = Transformed_chemical_name,
    values_from = result_value
  )

# Multiply PCB columns by 1000
pcb_columns <- names(transposed_data)[grep("^PCB\\d+", names(transposed_data))]
transposed_data[, pcb_columns] <- transposed_data[, pcb_columns] * 1000

# Update target_unit to "pg/l"
transposed_data$target_unit <- "pg/l"

# Change NA to 0s
transposed_data[is.na(transposed_data)] <- 0

# Function to combine PCB columns and remove original columns
combine_and_remove <- function(data, new_column, columns_to_combine) {
  data[[new_column]] <- rowSums(data[columns_to_combine], na.rm = TRUE)
  data <- data[, !colnames(data) %in% columns_to_combine]
  return(data)
}

# Specify the combinations and call the function for each
combinations <- list(
  list("PCB4+10", c("PCB4", "PCB10")),
  list("PCB5+8", c("PCB5", "PCB8")),
  list("PCB7+9", c("PCB7", "PCB9")),
  list("PCB16+32", c("PCB16", "PCB32")),
  list("PCB20+21+28+31+33+50+53", c("PCB20+28", "PCB21+33", "PCB31", "PCB50+53")),
  list("PCB24+27", c("PCB24", "PCB27")),
  list("PCB37+42", c("PCB37", "PCB42")),
  list("PCB40+41+64+71+72", c("PCB40+71", "PCB41", "PCB64", "PCB72")),
  list("PCB43+49+52+69+73", c("PCB43", "PCB49+69", "PCB52", "PCB73")),
  list("PCB45+51", c("PCB45", "PCB51")),
  list("PCB48+59+62+75", c("PCB48", "PCB59+62+75")),
  list("PCB56+60", c("PCB56", "PCB60")),
  list("PCB61+66+70+74+76+93+95+98+100+102", c("PCB61+70+74+76", "PCB66",
                                               "PCB93+100", "PCB95", "PCB98",
                                               "PCB102")),
  list("PCB77+85+110+111+115+116+117", c("PCB77", "PCB85+116", "PCB110", "PCB111",
                                         "PCB115", "PCB117")),
  list("PCB81+86+87+97+107+108+109+112+119+124+125", c("PCB81", "PCB86+87+97+108+119+125",
                                                       "PCB107+124", "PCB109",
                                                       "PCB112")),
  list("PCB82+135+144+151+154", c("PCB82", "PCB135+151", "PCB144", "PCB154")),
  list("PCB83+99", c("PCB83", "PCB99")),
  list("PCB84+92", c("PCB84", "PCB92")),
  list("PCB88+91", c("PCB88", "PCB91")),
  list("PCB106+118", c("PCB106", "PCB118")),
  list("PCB114+122+131+133+142+146+165", c("PCB114", "PCB122", "PCB131", "PCB133",
                                           "PCB142", "PCB146", "PCB165")),
  list("PCB123+139+140+147+149", c("PCB123", "PCB139+140", "PCB147+149")),
  list("PCB128+162+166+167", c("PCB128+166", "PCB162", "PCB167")),
  list("PCB129+137+138+158+160+163+164+176+178", c("PCB129+138+163", "PCB137",
                                                   "PCB158", "PCB160", "PCB164",
                                                   "PCB176", "PCB178")),
  list("PCB132+153+161+168", c("PCB132", "PCB153+168", "PCB161")),
  list("PCB134+143", c("PCB134", "PCB143")),
  list("PCB156+157+172+197+200", c("PCB156+157", "PCB172", "PCB197", "PCB200")),
  list("PCB170+190", c("PCB170", "PCB190")),
  list("PCB171+173+202", c("PCB171+173", "PCB202")),
  list("PCB182+187", c("PCB182", "PCB187")),
  list("PCB183+185", c("PCB183", "PCB185")),
  list("PCB194+205", c("PCB194", "PCB205")),
  list("PCB195+208", c("PCB195", "PCB208")),
  list("PCB196+203", c("PCB196", "PCB203")),
  list("PCB198+199+201", c("PCB198+199", "PCB201"))
)

# Apply the function for each combination
for (combo in combinations) {
  transposed_data <- combine_and_remove(transposed_data, combo[[1]], combo[[2]])
}

# Define the desired order of columns
desired_column_order <- c(
  "sys_sample_code", "subfacility_code", "sample_date", "y_coord", "x_coord", "target_unit", "analytic_method",
  "PCB1", "PCB2", "PCB3", "PCB4+10", "PCB5+8", "PCB6", "PCB7+9", "PCB11", "PCB12+13",
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

# Reorder the columns based on the desired order
transposed_data <- transposed_data[desired_column_order]

# Change the name of columns to be consistent
{
  colnames(transposed_data)[colnames(transposed_data) == "sys_sample_code"] <- "SampleID"
  colnames(transposed_data)[colnames(transposed_data) == "subfacility_code"] <- "SiteName"
  colnames(transposed_data)[colnames(transposed_data) == "sample_date"] <- "SampleDate"
  colnames(transposed_data)[colnames(transposed_data) == "y_coord"] <- "Latitude"
  colnames(transposed_data)[colnames(transposed_data) == "x_coord"] <- "Longitude"
  colnames(transposed_data)[colnames(transposed_data) == "target_unit"] <- "Units"
  colnames(transposed_data)[colnames(transposed_data) == "analytic_method"] <- "EPAMethod"
  }

# Remove "-NYC" from SiteID
transposed_data$SiteName <- gsub("-NYC", "", transposed_data$SiteName)

# Add a digit to the SiteName
transposed_data <- transposed_data %>%
  group_by(SiteName) %>%
  mutate(
    SiteName = if(n_distinct(Latitude) == 1) paste0(SiteName, "_1") else paste0(SiteName, "_", row_number())
  ) %>%
  ungroup()

# Create a SiteID column with three digits at the end
transposed_data <- transposed_data %>%
  mutate(SiteID = sprintf("WCPCB-NTC%03d", as.integer(factor(SiteName))))

transposed_data <- transposed_data %>%
  relocate(SiteID, .before = SampleDate)

# Format SampleID
# Extract last 8 characters from SampleID
transposed_data$SampleID <- str_sub(transposed_data$SampleID, start = -8)

# Combine SiteID and SampleID
transposed_data$SampleID <- paste(transposed_data$SiteID,
                                  transposed_data$SampleID, sep = "-")

transposed_data <- transposed_data %>%
  group_by(SampleID) %>%
  mutate(SampleID = ifelse(n() > 1, paste0(SampleID, ".", row_number()),
                           paste0(SampleID, ".1"))) %>%
  ungroup()

# Add new columns with metadata
# Names and values for the new columns
new_col_names <- c("EPARegion", "StateSampled", "LocationName")
new_col_values <- c("R2", "NY", "Newtown Creek")

# Create empty columns for the new columns
for (col_name in new_col_names) {
  transposed_data[[col_name]] <- NA_character_
}

# Insert new columns between SampleID and SiteName
transposed_data <- transposed_data %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character type
  mutate(across(all_of(new_col_names), ~NA_character_)) %>%  # Create empty columns
  relocate(SampleID, .after = "SampleID") %>%  # Move SampleID column to the beginning
  relocate(all_of(new_col_names), .after = "SampleID")  # Move new columns after SampleID

# Assign values to the new columns
for (i in seq_along(new_col_names)) {
  transposed_data[[new_col_names[i]]] <- new_col_values[i]
}

# Name and value for the new column
new_col_name <- c("PhaseMeasured")
new_col_value <- c("SurfaceWater")

# Create empty columns for the new columns
for (col_name in new_col_name) {
  transposed_data[[col_name]] <- NA_character_
}

# Insert new columns in desired positions
transposed_data <- transposed_data %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character type
  mutate(across(all_of(new_col_name), ~NA_character_)) %>%  # Create empty columns
  relocate(SampleID, .after = "SampleID") %>%  # Move SampleID column to the beginning
  relocate(all_of(new_col_name), .before = "EPAMethod")  # Move new columns before EPAMethod

# Assign values to the new columns
transposed_data$PhaseMeasured <- new_col_value[1]

# Name and value for the new column
new_col_name <- c("AroclorCongener")
new_col_value <- c("Congener")

# Create empty columns for the new columns
for (col_name in new_col_name) {
  transposed_data[[col_name]] <- NA_character_
}

# Insert new columns in desired positions
transposed_data <- transposed_data %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character type
  mutate(across(all_of(new_col_name), ~NA_character_)) %>%  # Create empty columns
  relocate(SampleID, .after = "SampleID") %>%  # Move SampleID column to the beginning
  relocate(all_of(new_col_name), .before = "PCB1")  # Move new columns before EPAMethod

transposed_data$AroclorCongener <- new_col_value[1]

# Name and value for the new column
transposed_data[, c("A1016", "A1221", "A1232", "A1242", "A1248",
                    "A1254", "A1260")] <- NA

# Transfor PCBs into numerical values
pcbs <- grep("^PCB", names(transposed_data), value = TRUE)
transposed_data[pcbs] <- lapply(transposed_data[pcbs], as.numeric)

# Create a new column named "tPCB" that sums columns 14 to 117
transposed_data <- transposed_data %>%
  mutate(tPCB = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Update the string in the "EPAMethod" column
transposed_data <- transposed_data %>%
  mutate(EPAMethod = ifelse(EPAMethod == "E1668A", "M1668", EPAMethod))

# Change coordinate system
# Modify 'Latitude' and 'Longitude' columns to be numeric
transposed_data$Latitude <- as.numeric(transposed_data$Latitude)
transposed_data$Longitude <- as.numeric(transposed_data$Longitude)

# Convert the data frame to an sf object
# First, create a simple feature object from the data frame by specifying the coordinates
transposed_sf <- st_as_sf(transposed_data, coords = c("Longitude", "Latitude"), crs = 2263) # Assuming your data is initially in WGS84

# Transform the coordinates to the target CRS
target_crs <- 4326 # EPSG code for WGS84
transformed_sf <- st_transform(transposed_sf, crs = target_crs)

# Transform coordinates back into a regular data frame
transposed_data_transformed <- as.data.frame(transformed_sf)

# Extract them into separate columns:
transposed_data$Longitude <- st_coordinates(transformed_sf)[, 1]
transposed_data$Latitude <- st_coordinates(transformed_sf)[, 2]

# Export results
write.csv(transposed_data, file = "Data/NewtownCreek/ntcV01.csv")


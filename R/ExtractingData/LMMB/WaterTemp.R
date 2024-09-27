# Code to read and change format of water temperatures for Lake Michigan 93-95
# https://www.ndbc.noaa.gov/station_history.php?station=45002

# Read the specific columns (YY, MM, DD, WTMP)
# Temperature in C
{
  WT1993 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1993.txt',
                       header = TRUE,
                       colClasses = c("numeric", "numeric",
                                      "numeric",
                                      "numeric"))[, c("YY", "MM", "DD", "WTMP")]
  
  WT1994 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1994.txt',
                       header = TRUE,
                       colClasses = c("numeric", "numeric",
                                      "numeric",
                                      "numeric"))[, c("YY", "MM", "DD", "WTMP")]
  
  WT1995 <- read.table('Data/LMMB/OpenLake/WaterTemp/45002h1995.txt',
                       header = TRUE,
                       colClasses = c("numeric", "numeric",
                                      "numeric",
                                      "numeric"))[, c("YY", "MM", "DD", "WTMP")]
}

# Combine all the datasets
LakeMichiganWT <- rbind(WT1993, WT1994, WT1995)

# Daily average, not including 999s
LakeMichiganWT <- aggregate(WTMP ~ YY + MM + DD, data = LakeMichiganWT,
                            mean, subset = WTMP != 999)

# Transforming temperature from C to K
LakeMichiganWT$WTMP_K <- LakeMichiganWT$WTMP + 273.15

# Create a new column with date
LakeMichiganWT$Date <- as.Date(paste(LakeMichiganWT$MM, LakeMichiganWT$DD,
                                     LakeMichiganWT$YY, sep = "/"),
                               format = "%m/%d/%y")

# Export results
write.csv(LakeMichiganWT,
          file = "Output/Data/Sites/csv/GreatLakes/WaterTemp/LakeMichiganWT.csv")


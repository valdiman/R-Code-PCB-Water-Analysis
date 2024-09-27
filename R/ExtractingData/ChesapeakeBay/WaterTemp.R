# Code to read and change format of water temperatures for Chesapeake Bay
# https://www.ndbc.noaa.gov/station_history.php?station=tplm2
# Read the specific columns (YY, MM, DD, WTMP)
# Temperature in C
{
  WT2001 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2001.txt',
                              header = TRUE), select = c("YYYY", "MM", "DD",
                                                         "WTMP"))
  WT2003 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2003.txt',
                              header = TRUE), select = c("YYYY", "MM", "DD",
                                                         "WTMP"))
  WT2006 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2006.txt',
                              header = TRUE), select = c("YYYY", "MM", "DD",
                                                         "WTMP"))
  WT2008 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2008.txt',
                              header = TRUE), select = c("X2008", "X01", "X01.1",
                                                         "X5.7"))
  # Rename the columns
  colnames(WT2008) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2010 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2010.txt',
                              header = TRUE), select = c("X2010", "X01", "X01.1",
                                                         "X3.6"))
  # Rename the columns
  colnames(WT2010) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2011 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2011.txt',
                              header = TRUE), select = c("X2011", "X01", "X01.1",
                                                         "X2.4"))
  # Rename the columns
  colnames(WT2011) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2012 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2012.txt',
                              header = TRUE), select = c("X2012", "X01", "X01.1",
                                                         "X7.2"))
  # Rename the columns
  colnames(WT2012) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2013 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2013.txt',
                              header = TRUE), select = c("X2013", "X01", "X01.1",
                                                         "X5.6"))
  # Rename the columns
  colnames(WT2013) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2014 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2014.txt',
                              header = TRUE), select = c("X2014", "X01", "X01.1",
                                                         "X4.4"))
  # Rename the columns
  colnames(WT2014) <- c("YYYY", "MM", "DD", "WTMP")
  
  WT2015 <- subset(read.table('Data/ChesapeakeBay/WaterTemp/tplm2h2015.txt',
                              header = TRUE), select = c("X2015", "X01", "X01.1",
                                                         "X5.2"))
  # Rename the columns
  colnames(WT2015) <- c("YYYY", "MM", "DD", "WTMP")
}

# Combine all the datasets
ChesapeakeBayWT <- rbind(WT2001, WT2003, WT2006, WT2008, WT2010, WT2011, WT2012, WT2013,
                         WT2014, WT2015)

# Daily average, not including 999s
ChesapeakeBayWT <- aggregate(WTMP ~ YYYY + MM + DD, data = ChesapeakeBayWT, mean,
                      subset = WTMP != 999)

ChesapeakeBayWT$Date <- as.Date(paste(ChesapeakeBayWT$MM, ChesapeakeBayWT$DD,
                                      ChesapeakeBayWT$YYYY, sep = "/"),
                                format = "%m/%d/%Y")

# Transforming temperature from C to K
ChesapeakeBayWT$WTMP_K <- ChesapeakeBayWT$WTMP + 273.15

# Export results
write.csv(ChesapeakeBayWT,
          file = "Output/Data/Sites/csv/ChesapeakeBay/WaterTemp/ChesapeakeBayWT.csv")



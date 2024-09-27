## Water PCB concentrations data analysis
## Housatonic River
## Only Linear Mixed-Effects Model (lme)
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @https://semspub.epa.gov/work/01/574882.pdf

# Install packages
install.packages("tidyverse")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")
install.packages("zoo")
install.packages("dataRetrieval")
install.packages("reshape")
install.packages("tidyr")
install.packages("scales")

# Load libraries
{
  library(scales) # function trans_breaks
  library(stringr) # str_detect
  library(robustbase) # function colMedians
  library(dplyr) # performs %>%
  library(tibble) # adds a column
  library(lme4) # performs lme
  library(MuMIn) # gets Rs from lme
  library(lmerTest) # gets the p-value from lme
  library(zoo) # yields seasons
  library(dataRetrieval) # read data from USGS
  library(reshape)
  library(tidyr) # function gather
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Select Housatonic River data ---------------------------------------------------
hou <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  hou$SampleDate <- as.Date(hou$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hou$SampleDate) - min(as.Date(hou$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou$SiteID), hou$SampleDate, as.matrix(hou$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- hou.tpcb$tPCB
time <- hou.tpcb$time
flow <- hou.tpcb$flow.1
site <- hou.tpcb$SiteID
season <- hou.tpcb$season

lme.hou.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb)

# Shapiro test
shapiro.test(resid(lme.hou.tpcb)) # p-value < 0.05

# Lme does not provide a good model for the data.

# Perform Linear Mixed-Effects Model (lme) w/ quadratic flow.
lme.hou.tpcb.2 <- lmer(log10(tpcb) ~ 1 + time + poly(flow, 2) + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb.2)

# Shapiro test
shapiro.test(resid(lme.hou.tpcb.2)) # p-value < 0.05


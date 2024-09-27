## Water PCB concentrations data analysis
## Hudson River
## Only Linear Mixed-Effects Model (lme)

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

# Select Hudson River data ---------------------------------------------------
hud <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Data preparation --------------------------------------------------------
{
  # Change date format
  hud$SampleDate <- as.Date(hud$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hud$SampleDate) - min(as.Date(hud$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hud.tpcb <- cbind(factor(hud$SiteID), hud$SampleDate,
                    hud$Latitude, hud$Longitude, as.matrix(hud$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "season")
}

# Remove site -------------------------------------------------------------
## Remove site Bakers Falls. Upstream source
## North Bakers Falls = WCPCB-HUD006 and
## South Bakers Falls = WCPCB-HUD006.
hud.tpcb.1 <- subset(hud.tpcb, SiteID != c("WCPCB-HUD006"))
hud.tpcb.1 <- subset(hud.tpcb.1, SiteID != c("WCPCB-HUD010"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Hudson River
{
  sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY No temp!
  sitehudN2 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY, no temp!
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN4 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitehudN1, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.3 <- readNWISdv(sitehudN3, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  flow.4 <- readNWISdv(sitehudN4, paramflow,
                       min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.tpcb.1$date), max(hud.tpcb.1$date))
  
  # Add USGS data to hud.tpcb.2, matching dates
  hud.tpcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.1$Date)]
  hud.tpcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.2$Date)]
  hud.tpcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.3$Date)]
  hud.tpcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.tpcb.1$date,
                                                       flow.4$Date)]
  hud.tpcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.tpcb.1$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  hud.tpcb.2 <- na.omit(hud.tpcb.1)
}

# LME Model tPCB --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- hud.tpcb.2$tPCB
time <- hud.tpcb.2$time
site <- hud.tpcb.2$SiteID
season <- hud.tpcb.2$season
flow <- hud.tpcb.2$flow.2 # flow.2
wtemp <- hud.tpcb.2$temp
# tPCB vs. time + season + flow + temp + site
lme.hud.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + wtemp + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"),
                     na.action = na.exclude)

# See results
summary(lme.hud.tpcb)

# Normality test
shapiro.test(resid(lme.hud.tpcb))$p.value # It doesn't work.

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  hud.pcb <- subset(hud, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  hud.pcb <- subset(hud.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  hud.pcb <- log10(hud.pcb)
  # Replace -inf to NA
  hud.pcb <- do.call(data.frame,
                     lapply(hud.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  hud.pcb.1 <- hud.pcb[,
                       -which(colSums(is.na(hud.pcb))/nrow(hud.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(hud$SiteID)
  # Change date format
  SampleDate <- as.Date(hud$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                             labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to hud.pcb.1
  hud.pcb.1 <- cbind(hud.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove site Bakers Falls. Upstream source
  # North Bakers Falls = WCPCB-HUD006 and
  # South Bakers Falls = WCPCB-HUD006.
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD006"))
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD010"))
  # Include flow data from USGS station Hudson River
  sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY No temp!
  sitehudN2 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY, no temp!
  sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
  sitehudN4 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY, no temp!
  sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Retrieve USGS data
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitehudN1, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  flow.2 <- readNWISdv(sitehudN2, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  flow.3 <- readNWISdv(sitehudN3, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  flow.4 <- readNWISdv(sitehudN4, paramflow,
                       min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  # Water temperature in Celsius
  temp <- readNWISdv(sitehudN5, paramtemp,
                     min(hud.pcb.1$SampleDate), max(hud.pcb.1$SampleDate))
  
  # Add USGS data to hud.tpcb.1 matching dates
  hud.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.1$Date)]
  hud.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.2$Date)]
  hud.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.3$Date)]
  hud.pcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.4$Date)]
  hud.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(hud.pcb.1$SampleDate,
                                                       temp$Date)]
  # Remove samples with flow.3 = NA
  hud.pcb.2 <- hud.pcb.1[!is.na(hud.pcb.1$flow.3), ]
  # Remove metadata for both hud.pcb.1 and hud.pcb.2
  hud.pcb.3 <- subset(hud.pcb.1, select = -c(SiteID:temp))
  # Flow.3 = NA was removed here. When using Flow.3, use only hud.pcb.4!
  hud.pcb.4 <- subset(hud.pcb.2, select = -c(SiteID:temp))
}

# Get covariates
time <- hud.pcb.2$time
flow <- hud.pcb.2$flow.3 # For flow.3, use hud.pcb.4 in lme (see line 235), not here
wtemp <- hud.pcb.2$temp
season <- hud.pcb.2$season
site <- hud.pcb.2$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(hud.pcb.4[1,]), ncol = 25)

# Perform LME
for (i in 1:length(hud.pcb.4[1,])) {
  fit <- lmer(hud.pcb.4[,i] ~ 1 + time + flow + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
  lme.pcb[i,1] <- fixef(fit)[1] # intercept
  lme.pcb[i,2] <- summary(fit)$coef[1,"Std. Error"] # intercept error
  lme.pcb[i,3] <- summary(fit)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.pcb[i,4] <- fixef(fit)[2] # time
  lme.pcb[i,5] <- summary(fit)$coef[2,"Std. Error"] # time error
  lme.pcb[i,6] <- summary(fit)$coef[2,"Pr(>|t|)"] # time p-value
  lme.pcb[i,7] <- fixef(fit)[3] # flow
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # flow error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.pcb[i,10] <- fixef(fit)[4] # temperature
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # temperature error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.pcb[i,13] <- fixef(fit)[5] # season 2
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 2 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,16] <- fixef(fit)[6] # season 3
  lme.pcb[i,17] <- summary(fit)$coef[6,"Std. Error"] # season 3 error
  lme.pcb[i,18] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,19] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,20] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,21] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,22] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,23] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,24] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- hud.pcb.4[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 25] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(hud.pcb.4[,1]),
                      ncol = length(hud.pcb.4[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(hud.pcb.4[1,]))

for (i in 1:length(hud.pcb.4[1,])) {
  fit <- lmer(hud.pcb.4[,i] ~ 1 + time + flow + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(hud.pcb.4[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(hud.pcb.4)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "temperature",
                       "temperature.error", "temperature.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Hudson River", nrow(lme.pcb)), lme.pcb)

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]
# Select only congeners with significant time coefficients
lme.pcb.t <- lme.pcb[lme.pcb$time.pv < 0.05, ]
# Select relevant columns
lme.pcb.t <- lme.pcb.t[, c("LocationName", "Congeners", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.pcb.t, file = "Output/Data/Sites/csv/HudsonRiver/HudsonRiverLmePCB.csv",
          row.names = FALSE)


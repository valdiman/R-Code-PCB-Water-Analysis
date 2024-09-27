## Water PCB concentrations data analysis per site
## Passaic River
## Lme using a quadratic function for the flow

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

# Select Passaic River data ---------------------------------------------------
pas <- wdc[str_detect(wdc$LocationName, 'Passaic River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  pas$SampleDate <- as.Date(pas$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(pas$SampleDate) - min(as.Date(pas$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  pas.tpcb <- cbind(factor(pas$SiteID), pas$SampleDate, as.matrix(pas$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(pas.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN2 <- "01379500" # No temp
  sitepasN3 <- "01389005" # No flow
  sitepasN4 <- "01389010" # No temp
  sitepasN5 <- "01389500" # No temp
  sitepasN6 <- "01389890" # No temp

  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.2 <- readNWISdv(sitepasN2, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.3 <- readNWISdv(sitepasN4, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.4 <- readNWISdv(sitepasN5, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  flow.5 <- readNWISdv(sitepasN6, paramflow,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.tpcb$date), max(pas.tpcb$date))
  
  # Add USGS data to pass.tpcb.2, matching dates, conversion to m3/s
  pas.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.1$Date)]
  pas.tpcb$flow.2 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.2$Date)]
  pas.tpcb$flow.3 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.3$Date)]
  pas.tpcb$flow.4 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.4$Date)]
  pas.tpcb$flow.5 <- 0.03*flow.1$X_00060_00003[match(pas.tpcb$date,
                                                      flow.5$Date)]
  pas.tpcb$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.tpcb$date,
                                                       temp$Date)]
}

# Remove site -------------------------------------------------------------
# Remove site located in the ocean. Possible typo in original coordinates.
pas.tpcb.1 <- subset(pas.tpcb, SiteID != c("WCPCB-PAS022"))

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- pas.tpcb$tPCB
time <- pas.tpcb$time
site <- pas.tpcb$SiteID
season <- pas.tpcb$season
flow <- pas.tpcb$flow.5
tem <- pas.tpcb$temp
# tPCB vs. time + season + flow + temp + site
lme.pas.tpcb <- lmer(log10(tpcb) ~ 1 + time + poly(flow, 2) + tem + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lme.pas.tpcb)

# Look at residuals
{
  res.pas.tpcb <- resid(lme.pas.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/PassaicRiverQ-QtPCB.pdf")
  qqnorm(res.pas.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.pas.tpcb)
  dev.off()
}
# Shapiro test
shapiro.test(resid(lme.pas.tpcb)) # Lme doesn't work.

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove site located in the ocean. Possible typo in original coordinates.
  pas.pcb <- subset(pas, SiteID != c("WCPCB-PAS022"))
  # Remove metadata
  pas.pcb.1 <- subset(pas.pcb, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  pas.pcb.1 <- subset(pas.pcb.1, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  pas.pcb.1 <- log10(pas.pcb.1)
  # Replace -inf to NA
  pas.pcb.1 <- do.call(data.frame,
                     lapply(pas.pcb.1,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  pas.pcb.2 <- pas.pcb.1[,
                       -which(colSums(is.na(pas.pcb.1))/nrow(pas.pcb.1) > 0.7)]
  # Add site ID
  SiteID <- factor(pas.pcb$SiteID)
  # Change date format
  SampleDate <- as.Date(pas.pcb$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(pas.pcb$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to pass.pcb
  pas.pcb.2 <- cbind(pas.pcb.2, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Include flow data from USGS station Passaic River
  sitepasN1 <- "01381900" # No temp
  sitepasN2 <- "01379500" # No temp
  sitepasN3 <- "01389005" # No flow
  sitepasN4 <- "01389010" # No temp
  sitepasN5 <- "01389500" # No temp
  sitepasN6 <- "01389890" # No temp
  
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(sitepasN1, paramflow,
                       min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  flow.2 <- readNWISdv(sitepasN2, paramflow,
                       min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  flow.3 <- readNWISdv(sitepasN4, paramflow,
                       min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  flow.4 <- readNWISdv(sitepasN5, paramflow,
                       min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  flow.5 <- readNWISdv(sitepasN6, paramflow,
                       min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  temp <- readNWISdv(sitepasN3, paramtemp,
                     min(pas.pcb.2$SampleDate), max(pas.pcb.2$SampleDate))
  
  # Add USGS data to pass.tpcb.1, matching dates, conversion to m3/s
  pas.pcb.2$flow.1 <- 0.03*flow.1$X_00060_00003[match(pas.pcb.2$SampleDate,
                                                      flow.1$Date)]
  pas.pcb.2$flow.2 <- 0.03*flow.2$X_00060_00003[match(pas.pcb.2$SampleDate,
                                                      flow.2$Date)]
  pas.pcb.2$flow.3 <- 0.03*flow.3$X_00060_00003[match(pas.pcb.2$SampleDate,
                                                      flow.3$Date)]
  pas.pcb.2$flow.4 <- 0.03*flow.4$X_00060_00003[match(pas.pcb.2$SampleDate,
                                                      flow.4$Date)]
  pas.pcb.2$flow.5 <- 0.03*flow.5$X_00060_00003[match(pas.pcb.2$SampleDate,
                                                      flow.5$Date)]
  pas.pcb.2$temp <- 273.15 + temp$X_.from.middle.intake_00010_00003[match(pas.pcb.2$SampleDate,
                                                                          temp$Date)]
  # Remove metadata
  pas.pcb.3 <- subset(pas.pcb.2, select = -c(SiteID:temp))
  # Remove samples with flow.3 = NA
  pas.pcb.4 <- pas.pcb.2[!is.na(pas.pcb.2$flow.3), ]
  # Remove metadata
  pas.pcb.5 <- subset(pas.pcb.4, select = -c(SiteID:temp))
}

# Get covariates
time <- pas.pcb.4$time
season <- pas.pcb.4$season
site <- pas.pcb.4$SiteID
flow <- pas.pcb.4$flow.3 # flow.3 yiels 6 PCB congeners
tem <- pas.pcb.4$temp

# Create matrix to store results
lme.pcb <- matrix(nrow = length(pas.pcb.5[1,]), ncol = 28)

# Perform LME
# Need to remove season 3, not enough data.
for (i in 1:length(pas.pcb.5[1,])) {
  fit <- lmer(pas.pcb.5[,i] ~ 1 + time + poly(flow, 2) + tem + season + (1|site),
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
  lme.pcb[i,10] <- fixef(fit)[4] # flow quad
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # flow error quad
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # flow p-value quad
  lme.pcb[i,13] <- fixef(fit)[4] # temp
  lme.pcb[i,14] <- summary(fit)$coef[4,"Std. Error"] # temp error
  lme.pcb[i,15] <- summary(fit)$coef[4,"Pr(>|t|)"] # temp p-value
  lme.pcb[i,16] <- fixef(fit)[4] # # season 1
  lme.pcb[i,17] <- summary(fit)$coef[5,"Std. Error"] # season 1 error
  lme.pcb[i,18] <- summary(fit)$coef[5,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,19] <- fixef(fit)[6] # season 2
  lme.pcb[i,20] <- summary(fit)$coef[6,"Std. Error"] # season 2 error
  lme.pcb[i,21] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,22] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,23] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,24] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,25] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,26] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,27] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- pas.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 28] <- sqrt(mean(residuals[non_na_indices]^2))
  
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(pas.pcb.5[,1]),
                      ncol = length(pas.pcb.5[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(pas.pcb.5[1,]))

for (i in 1:length(pas.pcb.5[1,])) {
  fit <- lmer(pas.pcb.5[,i] ~ 1 + time + poly(flow, 2) + tem + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(pas.pcb.5[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(pas.pcb.5)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "flow2", "flow2.error",
                       "flow2.pv", "temp", "temp.error", "temp.pv", "season1",
                       "season1.error", "season1.pv", "season2",
                       "season2.error", "season2, pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Passaic River Q", nrow(lme.pcb)), lme.pcb)

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
write.csv(lme.pcb.t, file = "Output/Data/Sites/csv/PassaicRiver/Quadratic/PassaicLmePCB.csv",
          row.names = FALSE)

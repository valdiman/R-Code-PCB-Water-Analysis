## Water PCB concentrations data analysis per site
## Portland Harbor

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

# Select Portland Harbor data ---------------------------------------------------
por <- wdc[str_detect(wdc$LocationName, 'Portland Harbor'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  por$SampleDate <- as.Date(por$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(por$SampleDate) - min(as.Date(por$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  por.tpcb <- cbind(factor(por$SiteID), por$SampleDate, as.matrix(por$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(por.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Include USGS flow data --------------------------------------------------
{
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR No!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  temp.1 <- readNWISdv(sitePorN1, paramtemp,
                       min(por.tpcb$date), max(por.tpcb$date))
  # Add USGS data to por.tpcb, matching dates
  por.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(por.tpcb$date, flow.1$Date)]
  por.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.tpcb$date, flow.2$Date)]
  por.tpcb$temp.1 <- 273.15 + temp.1$X_00010_00003[match(por.tpcb$date, temp.1$Date)]
  # Remove samples with temp = NA
  por.tpcb.2 <- na.omit(por.tpcb)
}

# LME Model tPCB --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- por.tpcb.2$tPCB
time <- por.tpcb.2$time
site <- por.tpcb.2$SiteID
season <- por.tpcb.2$season
flow <- por.tpcb.2$flow.1 # use 1
tem <- por.tpcb.2$temp
# tPCB vs. time + season + flow + temp + site
lme.por.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + tem + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.por.tpcb)

# Look at residuals
{
  res.por.tpcb <- resid(lme.por.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QPortlandHarbortPCB.pdf")
  qqnorm(res.por.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.por.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmePortlandHarbortPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(10^(predict(lme.por.tpcb)), resid(lme.por.tpcb),
       points(10^(predict(lme.por.tpcb)), resid(lme.por.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 2000, 500), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Shapiro test
shapiro.test(resid(lme.por.tpcb)) # p-value = 0.9216

# Create matrix to store results from lme analysis
{
  lme.tpcb <- matrix(nrow = 1, ncol = 28)
  lme.tpcb[1] <- fixef(lme.por.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.por.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.por.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.por.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.por.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.por.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.por.tpcb)[3] # flow
  lme.tpcb[8] <- summary(lme.por.tpcb)$coef[3,"Std. Error"] # flow error
  lme.tpcb[9] <- summary(lme.por.tpcb)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb[10] <- fixef(lme.por.tpcb)[4] # temperature
  lme.tpcb[11] <- summary(lme.por.tpcb)$coef[4,"Std. Error"] # temperature error
  lme.tpcb[12] <- summary(lme.por.tpcb)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.tpcb[13] <- fixef(lme.por.tpcb)[5] # season 1
  lme.tpcb[14] <- summary(lme.por.tpcb)$coef[5,"Std. Error"] # season 1 error
  lme.tpcb[15] <- summary(lme.por.tpcb)$coef[5,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[16] <- fixef(lme.por.tpcb)[6] # season 2
  lme.tpcb[17] <- summary(lme.por.tpcb)$coef[6,"Std. Error"] # season 2 error
  lme.tpcb[18] <- summary(lme.por.tpcb)$coef[6,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[19] <- fixef(lme.por.tpcb)[7] # season 3
  lme.tpcb[20] <- summary(lme.por.tpcb)$coef[7,"Std. Error"] # season 3 error
  lme.tpcb[21] <- summary(lme.por.tpcb)$coef[7,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[22] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[23] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[24] <- as.data.frame(VarCorr(lme.por.tpcb))[1,'sdcor']
  lme.tpcb[25] <- as.data.frame(r.squaredGLMM(lme.por.tpcb))[1, 'R2m']
  lme.tpcb[26] <- as.data.frame(r.squaredGLMM(lme.por.tpcb))[1, 'R2c']
  lme.tpcb[27] <- shapiro.test(resid(lme.por.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.por.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[28] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.por.tpcb <- as.data.frame(fitted(lme.por.tpcb))
# Add column name
colnames(fit.lme.values.por.tpcb) <- c("predicted")
# Add predicted values to data.frame
por.tpcb.2$predicted <- 10^(fit.lme.values.por.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
por.tpcb.2$factor2 <- por.tpcb.2$tPCB/por.tpcb.2$predicted
factor2.tpcb <- nrow(por.tpcb.2[por.tpcb.2$factor2 > 0.5 & por.tpcb.2$factor2 < 2,
                                ])/length(por.tpcb.2[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "flow", "flow.error", "flow.pv", "temperature",
                        "temperature.error", "temperature.pv", "season1",
                        "season1.error", "season1.pv", "season2",
                        "season2.error", "season2.pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                        "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Portland Harbor",
                                     nrow(lme.tpcb)), lme.tpcb)
# No significance on time coefficient.

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  por.pcb <- subset(por, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  por.pcb <- subset(por.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  por.pcb <- log10(por.pcb)
  # Replace -inf to NA
  por.pcb <- do.call(data.frame,
                     lapply(por.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  por.pcb.1 <- por.pcb[,
                       -which(colSums(is.na(por.pcb))/nrow(por.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(por$SiteID)
  # Change date format
  SampleDate <- as.Date(por$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to por.pcb.1
  por.pcb.1 <- cbind(por.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR No!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  temp <- readNWISdv(sitePorN1, paramtemp,
                       min(por.pcb.1$SampleDate), max(por.pcb.1$SampleDate))
  # Add USGS data to por.pcb.1, matching dates
  por.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(por.pcb.1$SampleDate,
                                                     flow.1$Date)]
  por.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.pcb.1$SampleDate,
                                                     flow.2$Date)]
  por.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(por.pcb.1$SampleDate,
                                                         temp$Date)]
  # Remove samples with temp = NA
  por.pcb.2 <- por.pcb.1[!is.na(por.pcb.1$temp), ]
  # Remove metadata
  por.pcb.3 <- subset(por.pcb.2, select = -c(SiteID:temp))
}

# Get covariates
time <- por.pcb.2$time
flow <- por.pcb.2$flow.1
temp <- por.pcb.2$temp
season <- por.pcb.2$season
site <- por.pcb.2$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(por.pcb.3[1,]), ncol = 25)

# Perform LME for each column
for (i in 1:length(por.pcb.3[1,])) {
  fit <- lmer(por.pcb.3[,i] ~ 1 + time + flow + temp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
    lme.pcb[i, 1] <- fixef(fit)[1] # intercept
    lme.pcb[i, 2] <- summary(fit)$coef[1, "Std. Error"] # intercept error
    lme.pcb[i, 3] <- summary(fit)$coef[1, "Pr(>|t|)"] # intercept p-value
    lme.pcb[i, 4] <- fixef(fit)[2] # time
    lme.pcb[i, 5] <- summary(fit)$coef[2, "Std. Error"] # time error
    lme.pcb[i, 6] <- summary(fit)$coef[2, "Pr(>|t|)"] # time p-value
    lme.pcb[i, 7] <- fixef(fit)[3] # flow
    lme.pcb[i, 8] <- summary(fit)$coef[3, "Std. Error"] # flow error
    lme.pcb[i, 9] <- summary(fit)$coef[3, "Pr(>|t|)"] # flow p-value
    lme.pcb[i, 10] <- fixef(fit)[4] # temperature
    lme.pcb[i, 11] <- summary(fit)$coef[4, "Std. Error"] # temperature error
    lme.pcb[i, 12] <- summary(fit)$coef[4, "Pr(>|t|)"] # temperature p-value
    lme.pcb[i, 13] <- fixef(fit)[5] # season 2
    lme.pcb[i, 14] <- summary(fit)$coef[5, "Std. Error"] # season 2 error
    lme.pcb[i, 15] <- summary(fit)$coef[5, "Pr(>|t|)"] # season 2 p-value
    lme.pcb[i, 16] <- fixef(fit)[6] # season 3
    lme.pcb[i, 17] <- summary(fit)$coef[6, "Std. Error"] # season 3 error
    lme.pcb[i, 18] <- summary(fit)$coef[6, "Pr(>|t|)"] # season 3 p-value
    lme.pcb[i, 19] <- -log(2) / lme.pcb[i, 4] / 365 # t0.5
    lme.pcb[i, 20] <- abs(-log(2) / lme.pcb[i, 4] / 365) * lme.pcb[i, 5] / abs(lme.pcb[i, 4]) # t0.5 error
    lme.pcb[i, 21] <- as.data.frame(VarCorr(fit))[1, 'sdcor']
    lme.pcb[i, 22] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
    lme.pcb[i, 23] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
    lme.pcb[i, 24] <- shapiro.test(resid(fit))$p.value
    # Calculate RMSE
    # Predictions
    predictions <- predict(fit)
    # Calculate residuals and RMSE
    residuals <- por.pcb.3[, i] - predictions
    non_na_indices <- !is.na(residuals)
    lme.pcb[i, 25] <- sqrt(mean(residuals[non_na_indices]^2))
    
  }

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(por.pcb.3[,1]),
                      ncol = length(por.pcb.3[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(por.pcb.3[1,]))

for (i in 1:length(por.pcb.3[1,])) {
  fit <- lmer(por.pcb.3[,i] ~ 1 + time + flow + temp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(por.pcb.3[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(por.pcb.3)
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
lme.pcb <- cbind(LocationName = rep("Portland Harbor", nrow(lme.pcb)), lme.pcb)

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
write.csv(lme.pcb.t,
          file = "Output/Data/Sites/csv/PortlandHarbor/PortlandHarborLmePCB.csv",
          row.names = FALSE)


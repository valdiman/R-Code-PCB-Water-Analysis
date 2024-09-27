## Water PCB concentrations data analysis
## Fox River 2005 - 2018
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
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Select Fox River data ---------------------------------------------------
fox <- wdc[str_detect(wdc$LocationName, 'Fox River'),]
# Lake Winnebago is a background site.
# Data preparation --------------------------------------------------------
{
  # Change date format
  fox$SampleDate <- as.Date(fox$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(fox$SampleDate) - min(as.Date(fox$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox$SiteID), fox$SampleDate, as.matrix(fox$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb.1 <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

# Include USGS flow and temperature data --------------------------------------------------
{
  # https://maps.waterdata.usgs.gov/mapper/index.html
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.tpcb.1$date), max(fox.tpcb.1$date))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.tpcb.1$date), max(fox.tpcb.1$date))
  # Add USGS data to fox.tpcb.2, matching dates, conversion to m3/s
  fox.tpcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb.1$date,
                                                                         flow$Date)]
  fox.tpcb.1$temp <- 273.15 + temp$X_00010_00003[match(fox.tpcb.1$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb.1 <- na.omit(fox.tpcb.1)
}

# LME Model tPCB ----------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- fox.tpcb.1$tPCB
time <- fox.tpcb.1$time
site <- fox.tpcb.1$SiteID
season <- fox.tpcb.1$season
flow <- fox.tpcb.1$flow
wtemp <- fox.tpcb.1$temp
# tPCB vs. time + season + flow + temp + site
lme.fox.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + wtemp + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.fox.tpcb)

# Look at residuals
{
  res.fox.tpcb <- resid(lme.fox.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QFoxRivertPCB.pdf")
  qqnorm(res.fox.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.fox.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeFoxRivertPCB.png", width = 800,
      height = 600)
  # Create plot
  plot(10^(predict(lme.fox.tpcb)), resid(lme.fox.tpcb),
       points(10^(predict(lme.fox.tpcb)), resid(lme.fox.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(0, 5000, 1000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Shapiro-Wilk normality  test
shapiro.test(resid(lme.fox.tpcb)) # p-value = 0.04757

# Create matrix to store results from lme analysis
{
  lme.tpcb <- matrix(nrow = 1, ncol = 25)
  lme.tpcb[1] <- fixef(lme.fox.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.fox.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.fox.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.fox.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.fox.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.fox.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.fox.tpcb)[3] # flow
  lme.tpcb[8] <- summary(lme.fox.tpcb)$coef[3,"Std. Error"] # flow error
  lme.tpcb[9] <- summary(lme.fox.tpcb)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb[10] <- fixef(lme.fox.tpcb)[4] # temperature
  lme.tpcb[11] <- summary(lme.fox.tpcb)$coef[4,"Std. Error"] # temperature error
  lme.tpcb[12] <- summary(lme.fox.tpcb)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.tpcb[13] <- fixef(lme.fox.tpcb)[5] # season 2
  lme.tpcb[14] <- summary(lme.fox.tpcb)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb[15] <- summary(lme.fox.tpcb)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[16] <- fixef(lme.fox.tpcb)[6] # season 3
  lme.tpcb[17] <- summary(lme.fox.tpcb)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb[18] <- summary(lme.fox.tpcb)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[19] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[20] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[21] <- as.data.frame(VarCorr(lme.fox.tpcb))[1,'sdcor']
  lme.tpcb[22] <- as.data.frame(r.squaredGLMM(lme.fox.tpcb))[1, 'R2m']
  lme.tpcb[23] <- as.data.frame(r.squaredGLMM(lme.fox.tpcb))[1, 'R2c']
  lme.tpcb[24] <- shapiro.test(resid(lme.fox.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.fox.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[25] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.fox.tpcb <- as.data.frame(fitted(lme.fox.tpcb))
# Add column name
colnames(fit.lme.values.fox.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.tpcb.1$predicted <- 10^(fit.lme.values.fox.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
fox.tpcb.1$factor2 <- fox.tpcb.1$tPCB/fox.tpcb.1$predicted
factor2.tpcb <- nrow(fox.tpcb.1[fox.tpcb.1$factor2 > 0.5 & fox.tpcb.1$factor2 < 2,
                                ])/length(fox.tpcb.1[,1])*100

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
                       "temperature.error", "temperature.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Fox River",
                                     nrow(lme.tpcb)), lme.tpcb)
# Select relevant columns
lme.tpcb.t <- lme.tpcb[, c("LocationName", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.tpcb.t, file = "Output/Data/Sites/csv/FoxRiver/FoxRiverLmetPCB.csv",
          row.names = FALSE)

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  fox.pcb <- subset(fox, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  fox.pcb <- subset(fox.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  fox.pcb <- log10(fox.pcb)
  # Replace -inf to NA
  fox.pcb <- do.call(data.frame,
                     lapply(fox.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  fox.pcb.1 <- fox.pcb[,
                       -which(colSums(is.na(fox.pcb))/nrow(fox.pcb) > 0.7)]
  
  # Add site ID
  SiteID <- factor(fox$SiteID)
  # Change date format
  SampleDate <- as.Date(fox$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  fox.pcb.1 <- cbind(fox.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove site Lake Winnebago (background site)
  fox.pcb.1 <- subset(fox.pcb.1, SiteID != c("WCPCB-FOX001"))
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  # Add USGS data to fox.pcb.1, matching dates, conversion to m3/s
  fox.pcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.pcb.1$SampleDate,
                                                                        flow$Date)]
  fox.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(fox.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with temperature = NA
  fox.pcb.2 <- fox.pcb.1[!is.na(fox.pcb.1$temp), ]
  # Remove metadata
  fox.pcb.3 <- subset(fox.pcb.2, select = -c(SiteID:temp))
}

# Get covariates
time <- fox.pcb.2$time
flow <- fox.pcb.2$flow
wtemp <- fox.pcb.2$temp
season <- fox.pcb.2$season
site <- fox.pcb.2$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(fox.pcb.3[1,]), ncol = 25)

# Perform LME
for (i in 1:length(fox.pcb.3[1,])) {
    fit <- lmer(fox.pcb.3[,i] ~ 1 + time + flow + wtemp + season + (1|site),
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
    residuals <- fox.pcb.3[, i] - predictions
    non_na_indices <- !is.na(residuals)
    lme.pcb[i, 25] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(fox.pcb.3[,1]),
                      ncol = length(fox.pcb.3[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(fox.pcb.3[1,]))

for (i in 1:length(fox.pcb.3[1,])) {
  fit <- lmer(fox.pcb.3[,i] ~ 1 + time + flow + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(fox.pcb.3[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(fox.pcb.3)
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
lme.pcb <- cbind(LocationName = rep("Fox River", nrow(lme.pcb)), lme.pcb)

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
write.csv(lme.pcb.t, file = "Output/Data/Sites/csv/FoxRiver/FoxRiverLmePCB.csv",
          row.names = FALSE)


## Water PCB concentrations data analysis per site
## Chesapeake Bay & Delaware Canal

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

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  che$SampleDate <- as.Date(che$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(che$SampleDate) - min(as.Date(che$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  che.tpcb <- cbind(factor(che$SiteID), che$SampleDate, as.matrix(che$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Add water temperature data ----------------------------------------------
# See code: R/ExtractingData/ChesapeakeBay/WaterTemp.R
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/WaterTemp/ChesapeakeBayWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  che.tpcb$temp <- wtp$WTMP_K[match(che.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  che.tpcb <- na.omit(che.tpcb)
}

# LME Model tPCB ----------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- che.tpcb$tPCB
time <- che.tpcb$time
wtemp <- che.tpcb$temp
site <- che.tpcb$SiteID
season <- che.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.che.tpcb <- lmer(log10(tpcb) ~ 1 + time + wtemp + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.che.tpcb)

# Look at residuals
{
  res.che.tpcb <- resid(lme.che.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QChesapeakeBaytPCB.pdf")
  qqnorm(res.che.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.che.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeChesapeakeBaytPCB.png", width = 800,
      height = 600)
  # Create your plot
  plot(10^(predict(lme.che.tpcb)), resid(lme.che.tpcb),
       points(10^(predict(lme.che.tpcb)), resid(lme.che.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = c(-2, 2), col = "grey")
  abline(v = seq(0, 30000, 5000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Shapiro test
shapiro.test(resid(lme.che.tpcb)) # p-value = 0.0731

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 25)
  lme.tpcb[1] <- fixef(lme.che.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.che.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.che.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.che.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.che.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.che.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.che.tpcb)[3] # water temperature
  lme.tpcb[8] <- summary(lme.che.tpcb)$coef[3,"Std. Error"] # water temperature error
  lme.tpcb[9] <- summary(lme.che.tpcb)$coef[3,"Pr(>|t|)"] # water temperature p-value
  lme.tpcb[10] <- fixef(lme.che.tpcb)[4] # season 1
  lme.tpcb[11] <- summary(lme.che.tpcb)$coef[4,"Std. Error"] # season 1 error
  lme.tpcb[12] <- summary(lme.che.tpcb)$coef[4,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[13] <- fixef(lme.che.tpcb)[5] # season 2
  lme.tpcb[14] <- summary(lme.che.tpcb)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb[15] <- summary(lme.che.tpcb)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[16] <- fixef(lme.che.tpcb)[6] # season 3
  lme.tpcb[17] <- summary(lme.che.tpcb)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb[18] <- summary(lme.che.tpcb)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[19] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[20] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[21] <- as.data.frame(VarCorr(lme.che.tpcb))[1,'sdcor']
  lme.tpcb[22] <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2m']
  lme.tpcb[23] <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2c']
  lme.tpcb[24] <- shapiro.test(resid(lme.che.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.che.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[25] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.che.tpcb <- as.data.frame(fitted(lme.che.tpcb))
# Add column name
colnames(fit.lme.values.che.tpcb) <- c("predicted")
# Add predicted values to data.frame
che.tpcb$predicted <- 10^(fit.lme.values.che.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
che.tpcb$factor2 <- che.tpcb$tPCB/che.tpcb$predicted
factor2.tpcb <- nrow(che.tpcb[che.tpcb$factor2 > 0.5 & che.tpcb$factor2 < 2,
                              ])/length(che.tpcb[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "temperature", "temperature.error", "temperature.pv",
                        "season1", "season1.error", "season1.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                        "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Chesapeake Bay",
                                     nrow(lme.tpcb)), lme.tpcb)
# p-value of time coefficient = 0.0533

# Select relevant columns
lme.tpcb.t <- lme.tpcb[, c("LocationName", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.tpcb.t, file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmetPCB.csv",
          row.names = FALSE)

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  che.pcb <- subset(che, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  che.pcb <- subset(che.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  che.pcb <- log10(che.pcb)
  # Replace -inf to NA
  che.pcb <- do.call(data.frame,
                     lapply(che.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  che.pcb.1 <- che.pcb[,
                       -which(colSums(is.na(che.pcb))/nrow(che.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(che$SiteID)
  # Change date format
  SampleDate <- as.Date(che$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  che.pcb.1 <- cbind(che.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Add water temperature to grl.pcb.1
  che.pcb.1$wtemp <- wtp$WTMP_K[match(che.pcb.1$SampleDate, wtp$Date)]
  # Remove metadata
  che.pcb.2 <- subset(che.pcb.1, select = -c(SiteID:wtemp))
}

# Get covariates
time <- che.pcb.1$time
wtemp <- che.pcb.1$wtemp
season <- che.pcb.1$season
site <- che.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(che.pcb.2[1,]), ncol = 25)

# Perform LME
for (i in 1:length(che.pcb.2[1,])) {
  fit <- lmer(che.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # water temperature
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # water temperature error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # water temperature p-value
  lme.pcb[i,10] <- fixef(fit)[4] # # season 1
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 1 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # # season 1 p-value
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
  residuals <- che.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 25] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(che.pcb.2[,1]),
                      ncol = length(che.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(che.pcb.2[1,]))

for (i in 1:length(che.pcb.2[1,])) {
  fit <- lmer(che.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(che.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(che.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "temperature", "temperature. error", "temperature.pv",
                       "season1", "season1.error", "season1.pv", "season2",
                       "season2.error", "season3, pv", "season3",
                       "season3.error", "season3, pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Chesapeake Bay", nrow(lme.pcb)), lme.pcb)

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
          file = "Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmePCB.csv",
          row.names = FALSE)

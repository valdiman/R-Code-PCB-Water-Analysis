## Water PCB concentrations data analysis per site
## Data from Lake Washington

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

# Select Lake Washington data ---------------------------------------------------
lwa <- wdc[str_detect(wdc$LocationName, 'Lake Washington'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  lwa$SampleDate <- as.Date(lwa$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(lwa$SampleDate) - min(as.Date(lwa$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(lwa$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  lwa.tpcb <- cbind(factor(lwa$SiteID), lwa$SampleDate, as.matrix(lwa$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(lwa.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# LME Model tPCB --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- lwa.tpcb$tPCB
time <- lwa.tpcb$time
site <- lwa.tpcb$SiteID
season <- lwa.tpcb$season
# tPCB vs. time + season + site
lme.lwa.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.lwa.tpcb)

# Look at residuals
{
  res.lwa.tpcb <- resid(lme.lwa.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QLakeWashingtontPCB.pdf")
  qqnorm(res.lwa.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.lwa.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeLakewashingtontPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(10^(predict(lme.lwa.tpcb)), resid(lme.lwa.tpcb),
       points(10^(predict(lme.lwa.tpcb)), resid(lme.lwa.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 10^6, 10^4), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Shapiro test
shapiro.test(resid(lme.lwa.tpcb)) # p-value = 0.2424

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 22)
  lme.tpcb[1] <- fixef(lme.lwa.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.lwa.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.lwa.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.lwa.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.lwa.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.lwa.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.lwa.tpcb)[3] # season 1
  lme.tpcb[8] <- summary(lme.lwa.tpcb)$coef[3,"Std. Error"] # season 1 error
  lme.tpcb[9] <- summary(lme.lwa.tpcb)$coef[3,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[10] <- fixef(lme.lwa.tpcb)[4] # season 2
  lme.tpcb[11] <- summary(lme.lwa.tpcb)$coef[4,"Std. Error"] # season 2 error
  lme.tpcb[12] <- summary(lme.lwa.tpcb)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[13] <- fixef(lme.lwa.tpcb)[5] # season 3
  lme.tpcb[14] <- summary(lme.lwa.tpcb)$coef[5,"Std. Error"] # season 3 error
  lme.tpcb[15] <- summary(lme.lwa.tpcb)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[16] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[17] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[18] <- as.data.frame(VarCorr(lme.lwa.tpcb))[1,'sdcor']
  lme.tpcb[19] <- as.data.frame(r.squaredGLMM(lme.lwa.tpcb))[1, 'R2m']
  lme.tpcb[20] <- as.data.frame(r.squaredGLMM(lme.lwa.tpcb))[1, 'R2c']
  lme.tpcb[21] <- shapiro.test(resid(lme.lwa.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.lwa.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.lwa.tpcb <- as.data.frame(fitted(lme.lwa.tpcb))
# Add column name
colnames(fit.lme.values.lwa.tpcb) <- c("predicted")
# Add predicted values to data.frame
lwa.tpcb$predicted <- 10^(fit.lme.values.lwa.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
lwa.tpcb$factor2 <- lwa.tpcb$tPCB/lwa.tpcb$predicted
factor2.tpcb <- nrow(lwa.tpcb[lwa.tpcb$factor2 > 0.5 & lwa.tpcb$factor2 < 2,
                              ])/length(lwa.tpcb[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season1", "season1.error", "season1.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                        "RMSE", "Factor2")


# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Lake Washington",
                                     nrow(lme.tpcb)), lme.tpcb)
# No significance on time coefficient

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  lwa.pcb <- subset(lwa, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  lwa.pcb <- subset(lwa.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  lwa.pcb <- log10(lwa.pcb)
  # Replace -inf to NA
  lwa.pcb <- do.call(data.frame,
                     lapply(lwa.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  lwa.pcb.1 <- lwa.pcb[,
                       -which(colSums(is.na(lwa.pcb))/nrow(lwa.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(lwa$SiteID)
  # Change date format
  SampleDate <- as.Date(lwa$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(lwa$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  lwa.pcb.1 <- cbind(lwa.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove metadata
  lwa.pcb.2 <- subset(lwa.pcb.1, select = -c(SiteID:season.s))
}

# Get covariates
time <- lwa.pcb.1$time
season <- lwa.pcb.1$season
site <- lwa.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(lwa.pcb.2[1,]), ncol = 22)

# Perform LME
for (i in 1:length(lwa.pcb.2[1,])) {
  fit <- lmer(lwa.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # season 1
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 1 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,10] <- fixef(fit)[4] # season 2
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 2 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,13] <- fixef(fit)[5] # season 3
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 3 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,16] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,17] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,18] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,19] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,20] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,21] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- lwa.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(lwa.pcb.2[,1]),
                      ncol = length(lwa.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(lwa.pcb.2[1,]))

for (i in 1:length(lwa.pcb.2[1,])) {
  fit <- lmer(lwa.pcb.2[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(lwa.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(lwa.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season1", "season1.error", "season1.pv",
                       "season2", "season2.error", "season2.pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Lake Washington", nrow(lme.pcb)), lme.pcb)

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
          file = "Output/Data/Sites/csv/LakeWashington/LakeWashingtonLmePCB.csv",
          row.names = FALSE)


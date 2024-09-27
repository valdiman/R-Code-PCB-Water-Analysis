## Water PCB concentrations data analysis per site
## Blue River or Bannister Fed Complex

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
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
install.packages('patchwork')
install.packages("scales")

# Load libraries
{
  library(ggplot2)
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
  library(patchwork) # combine plots
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Select Bannister Federal Complex data ---------------------------------------------------
bfc <- wdc[str_detect(wdc$LocationName, 'Bannister Fed Complex'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  bfc$SampleDate <- as.Date(bfc$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(bfc$SampleDate) - min(as.Date(bfc$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(bfc$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  bfc.tpcb <- cbind(factor(bfc$SiteID), bfc$SampleDate, as.matrix(bfc$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(bfc.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# LME Model tPCB --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- bfc.tpcb$tPCB
time <- bfc.tpcb$time
site <- bfc.tpcb$SiteID
season <- bfc.tpcb$season
lme.bfc.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.bfc.tpcb)

# Look at residuals
{
  res.bfc.tpcb <- resid(lme.bfc.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QBannisterFedComplextPCB.pdf")
  qqnorm(res.bfc.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.bfc.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeBannisterFedComplextPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(10^(predict(lme.bfc.tpcb)), resid(lme.bfc.tpcb),
       points(10^(predict(lme.bfc.tpcb)), resid(lme.bfc.tpcb), pch = 16, 
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
shapiro.test(resid(lme.bfc.tpcb)) # p-value = 0.4844

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 19)
  lme.tpcb[1] <- fixef(lme.bfc.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.bfc.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.bfc.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.bfc.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.bfc.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.bfc.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.bfc.tpcb)[3] # season 2
  lme.tpcb[8] <- summary(lme.bfc.tpcb)$coef[3,"Std. Error"] # season 2 error
  lme.tpcb[9] <- summary(lme.bfc.tpcb)$coef[3,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[10] <- fixef(lme.bfc.tpcb)[4] # season 3
  lme.tpcb[11] <- summary(lme.bfc.tpcb)$coef[4,"Std. Error"] # season 3 error
  lme.tpcb[12] <- summary(lme.bfc.tpcb)$coef[4,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[13] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[14] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[15] <- as.data.frame(VarCorr(lme.bfc.tpcb))[1,'sdcor']
  lme.tpcb[16] <- as.data.frame(r.squaredGLMM(lme.bfc.tpcb))[1, 'R2m']
  lme.tpcb[17] <- as.data.frame(r.squaredGLMM(lme.bfc.tpcb))[1, 'R2c']
  lme.tpcb[18] <- shapiro.test(resid(lme.bfc.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.bfc.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[19] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.bfc.tpcb <- as.data.frame(fitted(lme.bfc.tpcb))
# Add column name
colnames(fit.lme.values.bfc.tpcb) <- c("predicted")
# Add predicted values to data.frame
bfc.tpcb$predicted <- 10^(fit.lme.values.bfc.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
bfc.tpcb$factor2 <- bfc.tpcb$tPCB/bfc.tpcb$predicted
factor2.tpcb <- nrow(bfc.tpcb[bfc.tpcb$factor2 > 0.5 & bfc.tpcb$factor2 < 2,
                              ])/length(bfc.tpcb[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "season2", "season2.error", "season2.pv",
                        "season3", "season3.error", "season3.pv", "t05",
                        "t05.error", "RandonEffectSiteStdDev", "R2nR", "R2R",
                        "Normality", "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Bannister Fed Complex",
                                     nrow(lme.tpcb)), lme.tpcb)

# Time coefficients not significant.

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  bfc.pcb <- subset(bfc, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  bfc.pcb <- subset(bfc.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  bfc.pcb <- log10(bfc.pcb)
  # Replace -inf to NA
  bfc.pcb <- do.call(data.frame,
                     lapply(bfc.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  bfc.pcb.1 <- bfc.pcb[,
                       -which(colSums(is.na(bfc.pcb))/nrow(bfc.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(bfc$SiteID)
  # Change date format
  SampleDate <- as.Date(bfc$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(bfc$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to bfc.pcb.1
  bfc.pcb.1 <- cbind(bfc.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove metadata
  bfc.pcb.2 <- subset(bfc.pcb.1, select = -c(SiteID:season.s))
}

# Get covariates
time <- bfc.pcb.1$time
season <- bfc.pcb.1$season
site <- bfc.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(bfc.pcb.2[1,]), ncol = 19)

# Perform LME
for (i in 1:length(bfc.pcb.2[1,])) {
  fit <- lmer(bfc.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # # season 2
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 2 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 2 p-value
  lme.pcb[i,10] <- fixef(fit)[4] # season 3
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 3 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,13] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,14] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,15] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,16] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,17] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,18] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- bfc.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 19] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(bfc.pcb.2[,1]),
                      ncol = length(bfc.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(bfc.pcb.2[1,]))

for (i in 1:length(bfc.pcb.2[1,])) {
  fit <- lmer(bfc.pcb.2[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(bfc.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(bfc.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season2", "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Bannister Fed Complex",
                                    nrow(lme.tpcb)), lme.pcb)

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.05, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]
# Select only congeners with significant time coefficients
lme.pcb.t <- lme.pcb[lme.pcb$time.pv < 0.05, ]
# No congeners!

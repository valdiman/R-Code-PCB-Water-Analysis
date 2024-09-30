## Water PCB concentrations data analysis
## Lake Michigan Mass Balance & Great Lakes
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
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select LMMB and Great Lakes data ---------------------------------------------------
grl <- wdc[str_detect(wdc$LocationName, 'Lake Michigan Mass Balance|Great Lakes'), ]

# (1) Just get lake data, remove data from tributaries
grl <- grl[!grepl("^Tributary", grl$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Change date format
  grl$SampleDate <- as.Date(grl$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(grl$SampleDate) - min(as.Date(grl$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  grl.tpcb <- cbind(factor(grl$SiteID), grl$SampleDate, as.matrix(grl$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(grl.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Add water temperature data ----------------------------------------------
# See code: R/ExtractingData/LMMB/WaterTemp.R
{
  # Read water temperature
  wtp <- read.csv("Output/Data/Sites/csv/GreatLakes/WaterTemp/LakeMichiganWT.csv")
  # Convert date columns to Date format
  wtp$Date <- as.Date(wtp$Date)
  # Add water temperature to grl.tpcb
  grl.tpcb$temp <- wtp$WTMP_K[match(grl.tpcb$date, wtp$Date)]
  # Remove samples with temp = NA
  grl.tpcb <- na.omit(grl.tpcb)
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- grl.tpcb$tPCB
time <- grl.tpcb$time
wtemp <- grl.tpcb$temp
site <- grl.tpcb$SiteID
season <- grl.tpcb$season
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + wtemp + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)

# Shapiro test
shapiro.test(resid(lme.grl.tpcb))  # Lme doesn't work, p-value < 0.05

# (2) Samples from Lake Michigan
grl.tpcb.1 <- subset(grl.tpcb, grepl("LMM", SiteID))

# Using grl.tpcb.1
# Get variables
tpcb <- grl.tpcb.1$tPCB
time <- grl.tpcb.1$time
wtemp <- grl.tpcb.1$temp
site <- grl.tpcb.1$SiteID
season <- grl.tpcb.1$season
lme.grl.tpcb <- lmer(log10(tpcb) ~ 1 + time + wtemp + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.grl.tpcb)

# Shapiro test
shapiro.test(resid(lme.grl.tpcb))  # Lme doesn't work, p-value < 0.05

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Select only Lake Michigan data
  grl.pcb.0 <- subset(grl, grepl("LMM", SiteID))
  # Remove metadata
  grl.pcb <- subset(grl.pcb.0, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  grl.pcb <- subset(grl.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  grl.pcb <- log10(grl.pcb)
  # Replace -inf to NA
  grl.pcb <- do.call(data.frame,
                     lapply(grl.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  grl.pcb.1 <- grl.pcb[,
                       -which(colSums(is.na(grl.pcb))/nrow(grl.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(grl.pcb.0$SiteID)
  # Change date format
  SampleDate <- as.Date(grl.pcb.0$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(grl.pcb.0$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to grl.pcb.1
  grl.pcb.1 <- cbind(grl.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Add water temperature to grl.tpcb
  grl.pcb.1$temp <- wtp$WTMP_K[match(grl.pcb.1$SampleDate, wtp$Date)]
  # Remove samples with temp = NA
  grl.pcb.1 <- grl.pcb.1[!is.na(grl.pcb.1$temp), ]
  # Remove metadata
  grl.pcb.2 <- subset(grl.pcb.1, select = -c(SiteID:temp))
}

# Get covariates
time <- grl.pcb.1$time
wtemp <- grl.pcb.1$temp
season <- grl.pcb.1$season
site <- grl.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(grl.pcb.2[1,]), ncol = 22)

# Perform LME
for (i in 1:length(grl.pcb.2[1,])) {
  fit <- lmer(grl.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[2] # water temperature
  lme.pcb[i,8] <- summary(fit)$coef[2,"Std. Error"] # water temperature error
  lme.pcb[i,9] <- summary(fit)$coef[2,"Pr(>|t|)"] # water temperature p-value
  lme.pcb[i,10] <- fixef(fit)[3] # # season 2
  lme.pcb[i,11] <- summary(fit)$coef[3,"Std. Error"] # season 2 error
  lme.pcb[i,12] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 2 p-value
  lme.pcb[i,13] <- fixef(fit)[4] # season 3
  lme.pcb[i,14] <- summary(fit)$coef[4,"Std. Error"] # season 3 error
  lme.pcb[i,15] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 3 p-value
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
  residuals <- grl.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(grl.pcb.2[,1]),
                      ncol = length(grl.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(grl.pcb.2[1,]))

for (i in 1:length(grl.pcb.2[1,])) {
  fit <- lmer(grl.pcb.2[,i] ~ 1 + time + wtemp + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(grl.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(grl.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "temperature", "temperature.error", "temperature.pv",
                       "season2", "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("Lake Michigan", nrow(lme.pcb)), lme.pcb)

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
write.csv(lme.pcb.t, file = "Output/Data/Sites/csv/GreatLakes/GreatLakesLmePCB.csv",
          row.names = FALSE)


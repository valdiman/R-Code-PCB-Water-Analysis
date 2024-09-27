## Water PCB concentrations data analysis per site
## Tributaries from Lake Michigan Mass Balance

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

# Select tributary data from LMMB  data ---------------------------------------------------
# Just tributaries data
trb <- wdc[grepl("^Tributary", wdc$SiteName), ]

# Data preparation --------------------------------------------------------
{
  # Change date format
  trb$SampleDate <- as.Date(trb$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(trb$SampleDate) - min(as.Date(trb$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(trb$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  trb.tpcb <- cbind(factor(trb$SiteID), as.matrix(trb$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(trb.tpcb) <- c("SiteID", "tPCB", "time", "season")
}

# LME Model tPCB ----------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- trb.tpcb$tPCB
time <- trb.tpcb$time
site <- trb.tpcb$SiteID
season <- trb.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.trb.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.trb.tpcb)

# Look at residuals
{
  res.trb.tpcb <- resid(lme.trb.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QTributariestPCB.pdf")
  qqnorm(res.trb.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.trb.tpcb)
  dev.off()
}

# Shapiro test
shapiro.test(resid(lme.trb.tpcb)) # p-value < 0.05

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  trb.pcb <- subset(trb, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  trb.pcb <- subset(trb.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  trb.pcb <- log10(trb.pcb)
  # Replace -inf to NA
  trb.pcb <- do.call(data.frame,
                     lapply(trb.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  trb.pcb.1 <- trb.pcb[,
                       -which(colSums(is.na(trb.pcb))/nrow(trb.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(trb$SiteID)
  # Change date format
  SampleDate <- as.Date(trb$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(trb$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  trb.pcb.1 <- cbind(trb.pcb.1, SiteID, SampleDate, data.frame(time.day),
                     season.s)
  # Remove metadata
  trb.pcb.2 <- subset(trb.pcb.1, select = -c(SiteID:season.s))
}

# Get covariates
time <- trb.pcb.1$time
season <- trb.pcb.1$season
site <- trb.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(trb.pcb.2[1,]), ncol = 22)

# Perform LME
for (i in 1:length(trb.pcb.2[1,])) {
  fit <- lmer(trb.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # # season 1
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
  residuals <- trb.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(trb.pcb.2[,1]),
                      ncol = length(trb.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(trb.pcb.2[1,]))

for (i in 1:length(trb.pcb.2[1,])) {
  fit <- lmer(trb.pcb.2[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(trb.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(trb.pcb.2)
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
lme.pcb <- cbind(LocationName = rep("Lake Michigan tributaries",
                                    nrow(lme.pcb)), lme.pcb)

# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb$Normality <- as.numeric(lme.pcb$Normality)
# Get the congeners that are not showing normality
lme.pcb.out <- lme.pcb[lme.pcb$Normality < 0.045, ]
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.045, ]
# Select only congeners with significant time coefficients
lme.pcb.t <- lme.pcb[lme.pcb$time.pv < 0.05, ]
# Select relevant columns
lme.pcb.t <- lme.pcb.t[, c("LocationName", "Congeners", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.pcb.t,
          file = "Output/Data/Sites/csv/GreatLakes/Tributaries/TributariesLmePCB.csv",
          row.names = FALSE)


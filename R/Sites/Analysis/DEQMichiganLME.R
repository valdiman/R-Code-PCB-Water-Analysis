## LME analysis
## Data from DEQ dmihigan

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
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select DEQ Michigan data ---------------------------------------------------
dmi <- wdc[str_detect(wdc$LocationName, fixed("DEQ (Michigan)")),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  dmi$SampleDate <- as.Date(dmi$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(dmi$SampleDate) - min(as.Date(dmi$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(dmi$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  dmi.tpcb <- cbind(factor(dmi$SiteID), as.matrix(dmi$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(dmi.tpcb) <- c("SiteID", "tPCB", "time", "season")
}

# LME Model tPCB ----------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- dmi.tpcb$tPCB
time <- dmi.tpcb$time
site <- dmi.tpcb$SiteID
season <- dmi.tpcb$season
# tPCB vs. time + season + site
lme.dmi.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.dmi.tpcb)

# Shapiro test
shapiro.test(resid(lme.dmi.tpcb))
# lme doesn't work here. p-value = 5e-13

# LME for individual PCBs -------------------------------------------------
# Prepare data.frame
{
  # Remove metadata
  dmi.pcb <- subset(dmi, select = -c(Source:AroclorCongener))
  # Remove Aroclor data
  dmi.pcb <- subset(dmi.pcb, select = -c(A1016:tPCB))
  # Log10 individual PCBs 
  dmi.pcb <- log10(dmi.pcb)
  # Replace -inf to NA
  dmi.pcb <- do.call(data.frame,
                     lapply(dmi.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  dmi.pcb.1 <- dmi.pcb[,
                       -which(colSums(is.na(dmi.pcb))/nrow(dmi.pcb) > 0.7)]
  # Add site ID
  SiteID <- factor(dmi$SiteID)
  # Change date format
  SampleDate <- as.Date(dmi$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Change name time.day to time
  colnames(time.day) <- "time"
  # Include season
  yq.s <- as.yearqtr(as.yearmon(dmi$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Add date and time to fox.pcb.1
  dmi.pcb.1 <- cbind(dmi.pcb.1, SiteID, data.frame(time.day), season.s)
  # Remove metadata
  dmi.pcb.2 <- subset(dmi.pcb.1, select = -c(SiteID:season.s))
}

# Get covariates
time <- dmi.pcb.1$time
season <- dmi.pcb.1$season
site <- dmi.pcb.1$SiteID

# Create matrix to store results
lme.pcb <- matrix(nrow = length(dmi.pcb.2[1,]), ncol = 19)

# Perform LME
for (i in 1:length(dmi.pcb.2[1,])) {
  fit <- lmer(dmi.pcb.2[, i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i, 7] <- fixef(fit)[3] # season 2
  lme.pcb[i, 8] <- summary(fit)$coef[3, "Std. Error"] # season 2 error
  lme.pcb[i, 9] <- summary(fit)$coef[3, "Pr(>|t|)"] # # season 2 p-value
  lme.pcb[i, 10] <- fixef(fit)[4] # season 3
  lme.pcb[i, 11] <- summary(fit)$coef[4, "Std. Error"] # season 3 error
  lme.pcb[i, 12] <- summary(fit)$coef[4, "Pr(>|t|)"] # season 3 p-value
  lme.pcb[i, 13] <- -log(2)/lme.pcb[i, 4]/365 # t0.5
  lme.pcb[i, 14] <- abs(-log(2)/lme.pcb[i, 4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i, 15] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i, 16] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i, 17] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i, 18] <- shapiro.test(resid(fit))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(fit)
  # Calculate residuals and RMSE
  residuals <- dmi.pcb.2[, i] - predictions
  non_na_indices <- !is.na(residuals)
  lme.pcb[i, 19] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Transform result to data.frame so factor 2 can be included
lme.pcb <- as.data.frame(lme.pcb)

# Add factor of 2
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(dmi.pcb.2[,1]),
                      ncol = length(dmi.pcb.2[1,]))

# Create a vector to store factor 2 for each congener
factor2_vector <- numeric(length = length(dmi.pcb.2[1,]))

for (i in 1:length(dmi.pcb.2[1,])) {
  fit <- lmer(dmi.pcb.2[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  
  lme.fit.pcb[,i] <- fitted(fit)
  
  # Calculate factor2 for each congener
  factor2 <- 10^(lme.fit.pcb[, i])/10^(dmi.pcb.2[, i])
  factor2_vector[i] <- sum(factor2 > 0.5 & factor2 < 2,
                           na.rm = TRUE) / (sum(!is.na(factor2))) * 100
}

# Add factor 2 to lme.pcb data.frame
lme.pcb$factor2 <- factor2_vector

# Change number format of factor 2 to 3 significant figures
lme.pcb$factor2 <- formatC(signif(lme.pcb$factor2, digits = 3))

# Add congener names
congeners <- colnames(dmi.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))

# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season2", "season2.error", "season2.pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                       "RMSE", "Factor2")

# Add Location Name
lme.pcb <- cbind(LocationName = rep("DEQ MI", nrow(lme.pcb)), lme.pcb)

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
write.csv(lme.pcb.t, file = "Output/Data/Sites/csv/DEQMichigan/DEQMILmePCB.csv",
          row.names = FALSE)

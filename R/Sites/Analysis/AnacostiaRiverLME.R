## LME analysis per site
## Anacostia River
## data only tPCB

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
}

# Read data ---------------------------------------------------------------
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("Data/USAWaterPCB.csv")

# Select Anacostia River data ---------------------------------------------------
anr <- wdc[str_detect(wdc$LocationName, 'Anacostia River'),]

# Data preparation --------------------------------------------------------
{
  # Change date format
  anr$SampleDate <- as.Date(anr$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(anr$SampleDate) - min(as.Date(anr$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(anr$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  anr.tpcb <- cbind(factor(anr$SiteID), anr$SampleDate, as.matrix(anr$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(anr.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Include USGS flow and temperature data --------------------------------------------------
{
  # https://maps.waterdata.usgs.gov/mapper/index.html
  # Include flow data from USGS station Anacostia River
  siteanrN1 <- "01649500" # flow @ NORTHEAST BRANCH ANACOSTIA RIVER AT RIVERDALE, MD
  siteanrN2 <- "01651000" # flow @ NORTHWEST BR ANACOSTIA RIVER NR HYATTSVILLE, MD
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteanrN1, paramflow,
                     min(anr.tpcb$date), max(anr.tpcb$date))
  flow.2 <- readNWISdv(siteanrN2, paramflow,
                       min(anr.tpcb$date), max(anr.tpcb$date))
  temp.1 <- readNWISdv(siteanrN1, paramtemp,
                     min(anr.tpcb$date), max(anr.tpcb$date))
  # Dont use temp.2
  # temp.2 <- readNWISdv(siteanrN2, paramtemp,
  #                   min(anr.tpcb$date), max(anr.tpcb$date))
  # Add USGS data to anr.tpcb.2, matching dates, conversion to m3/s
  anr.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(anr.tpcb$date, flow.1$Date)]
  anr.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(anr.tpcb$date, flow.2$Date)]
  anr.tpcb$temp.1 <- 273.15 + temp.1$X_00010_00003[match(anr.tpcb$date,
                                                         temp.1$Date)]
  #anr.tpcb$temp.2 <- 273.15 + temp.2$X_Center.of.flow_00010_00003[match(anr.tpcb$date,
  #                                                                    temp.2$Date)]
}

# LME Model tPCB ----------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- anr.tpcb$tPCB
time <- anr.tpcb$time
flow.1 <- anr.tpcb$flow.1 # use this one
# flow.2 <- anr.tpcb$flow.2
wtemp <- anr.tpcb$temp.1
site <- anr.tpcb$SiteID
season <- anr.tpcb$season

# tPCB vs. time + flow + season + site
lme.anr.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow.1 + wtemp + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.anr.tpcb)

# Look at residuals
{
  res.anr.tpcb <- resid(lme.anr.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QAnacostiaRivertPCB.pdf")
  qqnorm(res.anr.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.anr.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeAnacostiaRivertPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(10^(predict(lme.anr.tpcb)), resid(lme.anr.tpcb),
       points(10^(predict(lme.anr.tpcb)), resid(lme.anr.tpcb), pch = 16, 
              col = "white"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 10000, 2000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Shapiro test
shapiro.test(resid(lme.anr.tpcb)) # p-value = 0.1738

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 28)
  lme.tpcb[1] <- fixef(lme.anr.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.anr.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.anr.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.anr.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.anr.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.anr.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.anr.tpcb)[3] # flow
  lme.tpcb[8] <- summary(lme.anr.tpcb)$coef[3,"Std. Error"] # flow error
  lme.tpcb[9] <- summary(lme.anr.tpcb)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb[10] <- fixef(lme.anr.tpcb)[4] # temperature
  lme.tpcb[11] <- summary(lme.anr.tpcb)$coef[4,"Std. Error"] # temperature error
  lme.tpcb[12] <- summary(lme.anr.tpcb)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.tpcb[13] <- fixef(lme.anr.tpcb)[5] # season 1
  lme.tpcb[14] <- summary(lme.anr.tpcb)$coef[5,"Std. Error"] # season 1 error
  lme.tpcb[15] <- summary(lme.anr.tpcb)$coef[5,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb[16] <- fixef(lme.anr.tpcb)[6] # season 2
  lme.tpcb[17] <- summary(lme.anr.tpcb)$coef[6,"Std. Error"] # season 2 error
  lme.tpcb[18] <- summary(lme.anr.tpcb)$coef[6,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[19] <- fixef(lme.anr.tpcb)[7] # season 3
  lme.tpcb[20] <- summary(lme.anr.tpcb)$coef[7,"Std. Error"] # season 3 error
  lme.tpcb[21] <- summary(lme.anr.tpcb)$coef[7,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[22] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[23] <- (-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[24] <- as.data.frame(VarCorr(lme.anr.tpcb))[1,'sdcor']
  lme.tpcb[25] <- as.data.frame(r.squaredGLMM(lme.anr.tpcb))[1, 'R2m']
  lme.tpcb[26] <- as.data.frame(r.squaredGLMM(lme.anr.tpcb))[1, 'R2c']
  lme.tpcb[27] <- shapiro.test(resid(lme.anr.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.anr.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[28] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.anr.tpcb <- as.data.frame(fitted(lme.anr.tpcb))
# Add column name
colnames(fit.lme.values.anr.tpcb) <- c("predicted")
# Add predicted values to data.frame
anr.tpcb$predicted <- 10^(fit.lme.values.anr.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
anr.tpcb$factor2 <- anr.tpcb$tPCB/anr.tpcb$predicted
factor2.tpcb <- nrow(anr.tpcb[anr.tpcb$factor2 > 0.5 & anr.tpcb$factor2 < 2,
                              ])/length(anr.tpcb[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "flow", "flow.error", "flow.pv",
                        "temp", "temp.error", "temp.pv",
                        "season1", "season1.error", "season1.pv",
                        "season2", "season2.error", "season2.pv",
                        "season3", "season3.error", "season3.pv", "t05",
                        "t05.error", "RandonEffectSiteStdDev", "R2nR", "R2R",
                        "Normality", "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Anacostia River",
                                     nrow(lme.tpcb)), lme.tpcb)
# Select relevant columns
lme.tpcb.t <- lme.tpcb[, c("LocationName", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.tpcb.t,
          file = "Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverLmetPCB.csv",
          row.names = FALSE)

# Modeling plot
# Plot prediction vs. observations, 1:1 line
p <- ggplot(anr.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "white") +
  scale_y_log10(limits = c(1, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# See plot
print(p)

# Save plot
ggsave("Output/Plots/Sites/ObsPred/AnacostiaRiver/AnacostiaRiverLmeObsPredtPCB.png",
       plot = p, width = 8, height = 8, dpi = 500)


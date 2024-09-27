## Water PCB concentrations data analysis per site
## Kalamazoo River
## Source: Allied Paper, Inc. Recycling of carbonless copy paper
## Aroclors 1242, 1254 and 1260, no congener analysis

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

# Select Kalamazoo data ---------------------------------------------------
kal <- wdc[str_detect(wdc$LocationName, 'Kalamazoo River'),]
# Superfund site from Morrow Dam (Kalamazoo River) to Lake Michigan
# and 30 miles of Portage Creek (south), Cork St and Portage Creek Cork St sites
# Dredging occurred at Plainwell Dam site.

# Data preparation --------------------------------------------------------
{
  # Change date format
  kal$SampleDate <- as.Date(kal$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(kal$SampleDate) - min(as.Date(kal$SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(kal$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  kal.tpcb <- cbind(factor(kal$SiteID), kal$SampleDate, as.matrix(kal$tPCB),
                    data.frame(time.day), season.s)
  # Add column names
  colnames(kal.tpcb) <- c("SiteID", "date", "tPCB", "time", "season")
}

# Remove site -------------------------------------------------------------
# Remove site PlainwellDam. Dredging = WCPCB-KAL023
kal.tpcb.1 <- subset(kal.tpcb, SiteID != c("WCPCB-KAL023"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Kalamazoo River, no temperature available
{
  siteKalN1 <- "04108660" # KALAMAZOO RIVER AT NEW RICHMOND, MI
  siteKalN2 <- "04106000" # KALAMAZOO RIVER AT COMSTOCK, MI
  # Codes to retrieve data
  paramflow <- "00060" # discharge
  # Flow (ft3/s)
  flow.1 <- readNWISdv(siteKalN1, paramflow,
                       min(kal.tpcb.1$date), max(kal.tpcb.1$date))
  flow.2 <- readNWISdv(siteKalN2, paramflow,
                       min(kal.tpcb.1$date), max(kal.tpcb.1$date))
  kal.tpcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(kal.tpcb.1$date,
                                                       flow.1$Date)]
  kal.tpcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(kal.tpcb.1$date,
                                                       flow.2$Date)]
  # Create flow, flow.3
  kal.tpcb.1$flow.3 <- ifelse(is.na(kal.tpcb.1$flow.1) == TRUE,
                              kal.tpcb.1$flow.2/0.46, kal.tpcb.1$flow.1)
  # Remove samples with flow.1 = NA
  kal.tpcb.2 <- na.omit(kal.tpcb.1)
}

# LME Model tPCB -------------------------------------------------------------
# Perform Linear Mixed-Effects Model (LMEM)
# Use kal.tpcb.2
tpcb <- kal.tpcb.2$tPCB
time <- kal.tpcb.2$time
site <- kal.tpcb.2$SiteID
season <- kal.tpcb.2$season
flow <- kal.tpcb.2$flow.1

lme.kal.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.kal.tpcb)

# Look at residuals
{
  res.kal.tpcb <- resid(lme.kal.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/Q-QKalamazooRivertPCB.pdf")
  qqnorm(res.kal.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.kal.tpcb)
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Sites/Residual/res_plotlmeKalamazooRivertPCB.png",
      width = 800, height = 600)
  # Create your plot
  plot(10^(predict(lme.kal.tpcb)), resid(lme.kal.tpcb),
       points(10^(predict(lme.kal.tpcb)), resid(lme.kal.tpcb), pch = 16, 
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

# Shapiro-Wilk normatily test
shapiro.test(resid(lme.kal.tpcb)) # p-value = 0.1079

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 22)
  lme.tpcb[1] <- fixef(lme.kal.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.kal.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.kal.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.kal.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.kal.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.kal.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.kal.tpcb)[3] # flow
  lme.tpcb[8] <- summary(lme.kal.tpcb)$coef[3,"Std. Error"] # flow error
  lme.tpcb[9] <- summary(lme.kal.tpcb)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb[10] <- fixef(lme.kal.tpcb)[4] # season 2
  lme.tpcb[11] <- summary(lme.kal.tpcb)$coef[4,"Std. Error"] # season 2 error
  lme.tpcb[12] <- summary(lme.kal.tpcb)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[13] <- fixef(lme.kal.tpcb)[5] # season 3
  lme.tpcb[14] <- summary(lme.kal.tpcb)$coef[5,"Std. Error"] # season 3 error
  lme.tpcb[15] <- summary(lme.kal.tpcb)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[16] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[17] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[18] <- as.data.frame(VarCorr(lme.kal.tpcb))[1,'sdcor']
  lme.tpcb[19] <- as.data.frame(r.squaredGLMM(lme.kal.tpcb))[1, 'R2m']
  lme.tpcb[20] <- as.data.frame(r.squaredGLMM(lme.kal.tpcb))[1, 'R2c']
  lme.tpcb[21] <- shapiro.test(resid(lme.kal.tpcb))$p.value
  # Calculate RMSE
  # Predictions
  predictions <- predict(lme.kal.tpcb)
  # Calculate residuals and RMSE
  residuals <- log10(tpcb) - predictions
  non_na_indices <- !is.na(residuals)
  lme.tpcb[22] <- sqrt(mean(residuals[non_na_indices]^2))
}

# Obtain observations and predictions
# Get predicted values tpcb
fit.lme.values.kal.tpcb <- as.data.frame(fitted(lme.kal.tpcb))
# Add column name
colnames(fit.lme.values.kal.tpcb) <- c("predicted")
# Add predicted values to data.frame
kal.tpcb.2$predicted <- 10^(fit.lme.values.kal.tpcb$predicted)
# Estimate a factor of 2 between observations and predictions
kal.tpcb.2$factor2 <- kal.tpcb.2$tPCB/kal.tpcb.2$predicted
factor2.tpcb <- nrow(kal.tpcb.2[kal.tpcb.2$factor2 > 0.5 & kal.tpcb.2$factor2 < 2,
                                ])/length(kal.tpcb.2[,1])*100

# Transform lme.tpcb to data.frame so factor 2 can be included
lme.tpcb <- as.data.frame(lme.tpcb)

# Add factor 2 to lme.pcb data.frame
lme.tpcb$factor2 <- factor2.tpcb

# Change number format of factor 2 to 3 significant figures
lme.tpcb$factor2 <- formatC(signif(lme.tpcb$factor2, digits = 3))

# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "flow", "flow.error", "flow.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality",
                        "RMSE", "Factor2")

# Add Location Name
lme.tpcb <- cbind(LocationName = rep("Kalamazoo River",
                                     nrow(lme.tpcb)), lme.tpcb)

# Select relevant columns
lme.tpcb.t <- lme.tpcb[, c("LocationName", "t05", "t05.error",
                           "R2R", "RMSE", "Factor2")]

# Export results
write.csv(lme.tpcb.t,
          file = "Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverLmetPCB.csv",
          row.names = FALSE)


## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Using log10 of the sum of PCB.

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
install.packages("reshape")
install.packages("sf")

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
  library(reshape)
  library(sf)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataPangaea20240606.csv")

# Total Concentration Analysis --------------------------------------------
# Data preparation
{
  # Change date format
  SampleDate <- as.Date(wdc$SampleDate, format = "%Y-%m-%d")
  # Calculate sampling time
  time.day <- data.frame(as.Date(SampleDate) - min(as.Date(SampleDate)))
  # Include season
  yq.s <- as.yearqtr(as.yearmon(wdc$SampleDate, "%Y-%m-%d") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  tpcb <- cbind(factor(wdc$SiteID), SampleDate,
                wdc$Latitude, wdc$Longitude, wdc$tPCB,
                data.frame(time.day), season.s)
  # Add column names
  colnames(tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                      "tPCB", "time", "season")
}

# Total PCB Regressions ---------------------------------------------------
# Get variables
tPCB <- tpcb$tPCB
time <- tpcb$time
site <- tpcb$SiteID
season <- tpcb$season

# Perform Linear Mixed-Effects Model (lme)
lme.tpcb <- lmer(log10(tPCB) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.tpcb)

#Create a Q-Q plot and save it.
{
  # Create a new PNG graphics device
  png("Output/Plots/Global/Q-QlmetPCB.png", width = 600, height = 500)
  res <- resid(lme.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res)
  # Close the PNG device
  dev.off()
}

# Plot residuals vs. predictions
{
  # Open a PNG graphics device
  png("Output/Plots/Global/res_plotlmetPCB.png", width = 600, height = 500)
  # Create your plot
  plot(10^(predict(lme.tpcb)), resid(lme.tpcb),
       points(10^(predict(lme.tpcb)), resid(lme.tpcb), pch = 16, col = "white"),
       ylim = c(-4, 4),
       xlim = c(1, 10^6.1),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  # Add lines to the plot
  abline(0, 0)
  abline(h = seq(-4, 4, 1), col = "grey")
  abline(v = seq(0, 1200000, 200000), col = "grey")
  # Close the PNG graphics device
  dev.off()
}

# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.tpcb))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.tpcb <- as.data.frame(fitted(lme.tpcb))
# Add column name
colnames(fit.values.tpcb) <- c("lme.predicted")
# Add predicted values to data.frame
tpcb$lmepredicted <- 10^(fit.values.tpcb$lme.predicted)

# Plot prediction vs. observations, 1:1 line
tPCBObsPred <- ggplot(tpcb, aes(x = tPCB, y = lmepredicted)) +
  geom_point(shape = 21, size = 2, fill = "white") +
  scale_y_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue",
              linewidth = 0.7) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue",
              linewidth = 0.7) + # 2:1 line (factor of 2)
  theme_bw() +
  annotation_logticks(sides = "bl")

# See plot
print(tPCBObsPred)

# Save plot in folder
ggsave("Output/Plots/Global/tPCBObsPred.png", plot = tPCBObsPred,
       width = 6, height = 5, dpi = 500)

# Individual PCB Regressions ----------------------------------------------
# PCB5.8
pcb5.8 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB5.8,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb5.8) <- c("SiteID", "date", "PCB5.8", "time", "season")
# Remove 0s and NA values
pcb5.8 <- pcb5.8[complete.cases(pcb5.8$PCB5.8) & pcb5.8$PCB5.8 != 0, ]

# Get variables
PCBi <- pcb5.8$PCB5.8
time <- pcb5.8$time
site <- pcb5.8$SiteID
season <- pcb5.8$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb5.8 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb5.8)

# Shapiro test
shapiro.test(resid(lme.pcb5.8)) # p-value <<< 0.5

# PCB11
pcb11 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB11,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb11) <- c("SiteID", "date", "PCB11", "time", "season")
# Remove 0s and NA values
pcb11 <- pcb11[complete.cases(pcb11$PCB11) & pcb11$PCB11 != 0, ]

# Get variables
PCBi <- pcb11$PCB11
time <- pcb11$time
site <- pcb11$SiteID
season <- pcb11$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb11 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb11)

# Shapiro test
shapiro.test(resid(lme.pcb11)) # p-value <<< 0.5

# PCB18.30
pcb18.30 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB18.30,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb18.30) <- c("SiteID", "date", "PCB18.30", "time", "season")
# Remove 0s and NA values
pcb18.30 <- pcb18.30[complete.cases(pcb18.30$PCB18.30) & pcb18.30$PCB18.30 != 0, ]

# Get variables
PCBi <- pcb18.30$PCB18.30
time <- pcb18.30$time
site <- pcb18.30$SiteID
season <- pcb18.30$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb18.30 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb18.30)

# Shapiro test
shapiro.test(resid(lme.pcb18.30)) # p-value <<< 0.5

# PCB20.21.28.31.33.50.53
pcb20 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB20.21.28.31.33.50.53,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb20) <- c("SiteID", "date", "PCB20", "time", "season")
# Remove 0s and NA values
pcb20 <- pcb20[complete.cases(pcb20$PCB20) & pcb20$PCB20 != 0, ]

# Get variables
PCBi <- pcb20$PCB20
time <- pcb20$time
site <- pcb20$SiteID
season <- pcb20$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb20 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb20)

# Shapiro test
shapiro.test(resid(lme.pcb20)) # p-value <<< 0.5

# PCB44+47+65
pcb44 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB44.47.65,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb44) <- c("SiteID", "date", "PCB44", "time", "season")
# Remove 0s and NA values
pcb44 <- pcb44[complete.cases(pcb44$PCB44) & pcb44$PCB44 != 0, ]

# Get variables
PCBi <- pcb44$PCB44
time <- pcb44$time
site <- pcb44$SiteID
season <- pcb44$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb44 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb44)

# Shapiro test
shapiro.test(resid(lme.pcb44)) # p-value <<< 0.5

# PCB 67
pcb67 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB67,
               data.frame(time.day), season.s)
# Add column names
colnames(pcb67) <- c("SiteID", "date", "PCB67", "time", "season")
# Remove 0s and NA values
pcb67 <- pcb67[complete.cases(pcb67$PCB67) & pcb67$PCB67 != 0, ]

# Get variables
PCBi <- pcb67$PCB67
time <- pcb67$time
site <- pcb67$SiteID
season <- pcb67$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb67 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb67)

# Shapiro test
shapiro.test(resid(lme.pcb67)) # p-value <<< 0.5

# PCB 106+118
pcb106.118 <- cbind(factor(wdc$SiteID), SampleDate, wdc$PCB106.118,
                    data.frame(time.day), season.s)
# Add column names
colnames(pcb106.118) <- c("SiteID", "date", "PCB106.118", "time", "season")
# Remove 0s and NA values
pcb106.118 <- pcb106.118[complete.cases(pcb106.118$PCB106.118) & pcb106.118$PCB106.118 != 0, ]

# Get variables
PCBi <- pcb106.118$PCB106.118
time <- pcb106.118$time
site <- pcb106.118$SiteID
season <- pcb106.118$season

# Perform Linear Mixed-Effects Model (lme)
lme.pcb106.118 <- lmer(log10(PCBi) ~ 1 + time + season + (1|site),
                       REML = FALSE,
                       control = lmerControl(check.nobs.vs.nlev = "ignore",
                                             check.nobs.vs.rankZ = "ignore",
                                             check.nobs.vs.nRE = "ignore"))

# See results
summary(lme.pcb106.118)

# Shapiro test
shapiro.test(resid(lme.pcb106.118)) # p-value <<< 0.5

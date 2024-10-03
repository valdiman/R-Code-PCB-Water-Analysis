
# Install packages
install.packages("ggplot2")
install.packages("scales")
install.packages("viridis")
install.packages("dplyr")

# Load libraries
{
  library(ggplot2)
  library(scales)
  library(viridis)
  library(dplyr)
}

# Read data RFtPCB.csv ----------------------------------------------------
rf_tPCB <- read.csv("Output/Data/Sites/csv/Summary/RFtPCB.csv")

# Format data for plotting
# Round values
rf_tPCB$Pearson <- round(rf_tPCB$Pearson, 2)
rf_tPCB$RMSE <- round(rf_tPCB$RMSE, 2)
rf_tPCB$Factor2_Percentage <- round(rf_tPCB$Factor2_Percentage, 0)
# Calculate correlation and # of samples for each location
correlation_by_location <- aggregate(rf_tPCB$Pearson, by = list(rf_tPCB$location),
                                     FUN = function(x) round(mean(x), 2))
# Rename columns for clarity
colnames(correlation_by_location) <- c("location", "corr_RMSE")
# Merge with original data
rf_tPCB <- merge(rf_tPCB, correlation_by_location, by = "location", all.x = TRUE)
# Create site_number column with correlation information
rf_tPCB$corr_RMSE <- paste0(rf_tPCB$location, " (", rf_tPCB$corr_RMSE, ")")
# Remove closing parenthesis from correlation2
rf_tPCB$corr_RMSE <- gsub("\\)$", "", rf_tPCB$corr_RMSE)
# Combine correlation2_no_parenthesis and RMSE into a single column with the desired format
rf_tPCB$corr_RMSE <- paste0("(", rf_tPCB$corr_RMSE, ", ",
                                        rf_tPCB$RMSE, ")")
# Remove the initial parenthesis from correlation2_combined
rf_tPCB$corr_RMSE <- gsub("^\\(", "", rf_tPCB$corr_RMSE)
# Convert correlation2 to a factor
rf_tPCB$corr_RMSE <- factor(rf_tPCB$corr_RMSE)
# Reverse the levels of the factor
rf_tPCB$corr_RMSE <- factor(rf_tPCB$corr_RMSE,
                                        levels = rev(levels(rf_tPCB$corr_RMSE)))

# Need to change Richardson Hill Road Landfill to Richardson Landfill to fit in the plot
rf_tPCB$location <- gsub("Richardson Hill Road Landfill",
                         "Richardson Landfill", rf_tPCB$location)

# Create the new column
rf_tPCB <- rf_tPCB %>%
  group_by(location) %>%
  mutate(summary = paste0(location, " (",
                          round(mean(RMSE), 2), ", ",
                          round(mean(Pearson), 2), ", ",
                          unique(Factor2_Percentage), ")"))

CombinePredObstPCBV2 <- ggplot(rf_tPCB, aes(x = Observation, y = Predicted, fill = summary)) +
  geom_point(shape = 21, size = 4.5, alpha = 0.5, color = "black") +  # Fill aesthetic already specified in aes()
  scale_fill_viridis_d(name = NULL, option = "plasma") + # Proper discrete color scale for 'corr_RMSE'
  xlab(expression(bold("Observed " *Sigma* "PCB (pg/L)"))) +
  ylab(expression(bold("Predicted " *Sigma* "PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) + # 1:1 line
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line
  theme_bw() +
  theme(aspect.ratio = 1,
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    axis.title = element_text(size = 24),
    legend.position = c(0.29, 0.785),
    legend.text = element_text(size = 19),
    legend.title = element_blank(),  # Remove legend title
    legend.background = element_rect(fill = "transparent")) +
  scale_y_log10(limits = c(1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(1, 10^8),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb")

CombinePredObstPCBV2

# Save plot
ggsave("Output/Plots/Sites/ObsPred/Summary/CombineRFObsPredtPCB.png",
       plot = CombinePredObstPCBV2, width = 10, height = 10, dpi = 1300)

# Read data RFPCB.csv ----------------------------------------------------
rf_PCB <- read.csv("Output/Data/Sites/csv/Summary/RFPCB.csv")

# Calculate count of individual congeners per location
congener_counts <- rf_PCB %>%
  group_by(Location) %>%
  summarize(congener_count = n_distinct(Congener))

# Merge congener counts with main dataset
rf_PCB <- left_join(rf_PCB, congener_counts, by = "Location")

# Create legend text including congener count
rf_PCB$Legend_Text <- paste(rf_PCB$Location,
                            " (n=", rf_PCB$congener_count, ")", sep = "")

# Plot
CombinePredObsPCBV1 <- ggplot(rf_PCB, aes(x = 10^(Test_Data), y = 10^(Predicted_Data), fill = Legend_Text)) +
  geom_point(shape = 21, size = 3, alpha = 0.5, color = "black") +  # Fill aesthetic already specified in aes()
  scale_fill_viridis_d(name = NULL, option = "plasma") + # Proper discrete color scale for 'Location'
  xlab(expression(bold("Observed PCBi (pg/L)"))) +
  ylab(expression(bold("Predicted PCBi (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "black", linewidth = 0.7) + # 1:1 line
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.7) + # 1:2 line
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.7) + # 2:1 line
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    axis.title = element_text(size = 24),
    legend.position = c(0.24, 0.83),
    legend.text = element_text(size = 19),
    legend.title = element_blank(),  # Remove legend title
    legend.key.height = unit(0.5, "lines"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.background = element_rect(fill = "transparent")) +
  scale_y_log10(limits = c(10^-4, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^-4, 10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "lb")

CombinePredObsPCBV1

# Save plot
ggsave("Output/Plots/Sites/ObsPred/Summary/CombineRFObsPredPCB.png",
       plot = CombinePredObsPCBV1, width = 10, height = 10, dpi = 1300)

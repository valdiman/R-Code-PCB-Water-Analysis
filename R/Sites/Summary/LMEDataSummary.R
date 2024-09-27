# Results from the LME model, including t05, its error,
# RSME, R2 and factor of 2.
# Data only normality (p-value) > 0.05,
# and time coefficient is significant (p-value < 0.05)

# Read generated data for total PCB ---------------------------------------
{
  # Anacostia River
  anr <- read.csv("Output/Data/Sites/csv/AnacostiaRiver/AnacostiaRiverLmetPCB.csv")
  # Chesapeake Bay
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmetPCB.csv")
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverLmetPCB.csv")
  # Kalamazoo River
  kal <- read.csv("Output/Data/Sites/csv/KalamazooRiver/KalamazooRiverLmetPCB.csv")
  # Kalamazoo River (quadratic flow function)
  kal.2 <- read.csv("Output/Data/Sites/csv/KalamazooRiver/Quadratic/KalamazooRiverLmetPCB.csv")
  # New Bedford Harbor
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHLmetPCB.csv")
  # Combine the data frames
  combined_tPCB <- rbind(anr, che,  fox, kal, kal.2, nbh)
  colnames(combined_tPCB) <- c("LocationName", "t05", "t05.error", "R2",
                               "RMSE", "Factor2")
}

print(combined_tPCB)

# Export results
write.csv(combined_tPCB, file = "Output/Data/Sites/csv/Summary/AllLmetPCB.csv",
          row.names = FALSE)

# Read generated data for PCB ---------------------------------------------
{
  # DEQ MI
  dmi <- read.csv("Output/Data/Sites/csv/DEQMichigan/DEQMILmePCB.csv")
  # Chesapeake Bay
  che <- read.csv("Output/Data/Sites/csv/ChesapeakeBay/ChesapeakeLmePCB.csv")
  # Fox River
  fox <- read.csv("Output/Data/Sites/csv/FoxRiver/FoxRiverLmePCB.csv")
  # Fox River (Quadratic flow function)
  fox.2 <- read.csv("Output/Data/Sites/csv/FoxRiver/Quadratic/FoxRiverLmeQuadPCB.csv")
  # Lake Michigan
  grl <- read.csv("Output/Data/Sites/csv/GreatLakes/GreatLakesLmePCB.csv")
  # Hudson River
  hud <- read.csv("Output/Data/Sites/csv/HudsonRiver/HudsonRiverLmePCB.csv")
  # Lake Washington
  lwa <- read.csv("Output/Data/Sites/csv/LakeWashington/LakeWashingtonLmePCB.csv")
  # New Bedford Harbor
  nbh <- read.csv("Output/Data/Sites/csv/NewBedfordHarbor/NBHLmePCB.csv")
  # Passaic River
  pas <- read.csv("Output/Data/Sites/csv/PassaicRiver/PassaicLmePCB.csv")
  # Passaic River (Quadratic flow function)
  pas.2 <- read.csv("Output/Data/Sites/csv/PassaicRiver/Quadratic/PassaicLmePCB.csv")
  # Portland Harbor
  por <- read.csv("Output/Data/Sites/csv/PortlandHarbor/PortlandHarborLmePCB.csv")
  # Portland Harbor (Quadratic flow function)
  por.2 <- read.csv("Output/Data/Sites/csv/PortlandHarbor/Quadratic/PortlandHarborLmeQuadPCB.csv")
  # Spokane River (Quadratic flow function)
  spo.2 <- read.csv("Output/Data/Sites/csv/SpokaneRiver/Quadratic/SpokaneRiverLmeQuadPCB.csv")
  # Lake Michigan tributaries
  trb <- read.csv("Output/Data/Sites/csv/GreatLakes/Tributaries/TributariesLmePCB.csv")
  # Combine the data frames
  combined_PCB <- rbind(dmi, che, fox, fox.2, grl, hud, lwa, nbh, pas, pas.2,
                     por, por.2, spo.2, trb)
  colnames(combined_PCB) <- c("LocationName", "Congeners", "t05", "t05.error", "R2",
                               "RMSE", "Factor2")
}

print(combined_PCB)

# Export results
write.csv(combined_PCB, file = "Output/Data/Sites/csv/Summary/AllLmePCB.csv",
          row.names = FALSE)



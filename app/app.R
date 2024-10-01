# Load the data (adjust the path as necessary)
setwd("app")
# Data (pg/L) downloaded from Pangaea using code: R/Pangaea/PangaeaDownloadDataset.R
wdc <- read.csv("../Data/USAWaterPCB.csv")

# Source the new server and UI files
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)

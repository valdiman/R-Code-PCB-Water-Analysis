# Load the data (adjust the path as necessary)
setwd("app")
wdc <- read.csv("WaterDataCongenerAroclor09072023.csv")

# Source the new server and UI files
source("uiV02.R")
source("serverV02.R")

# Run the application
shinyApp(ui = ui, server = server)

# Install required packages if they are not already installed
packages <- c("shiny", "leaflet", "ggplot2", "data.table", "rsconnect")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# Load necessary libraries
library(shiny)
library(leaflet)
library(ggplot2)
library(data.table)
library(rsconnect)

# Need to indicate where wdc is
setwd("app")

# Set account information for ShinyApps.io
rsconnect::setAccountInfo(name='valdiman', token='0F7C710D1AF5503E411EF82C60E7F695',
                          secret='eecFcwuY6FM97Qvet85zDZJtjF9pth/JsTk6OhxQ')

# Deploy the app to ShinyApps.io
rsconnect::deployApp()

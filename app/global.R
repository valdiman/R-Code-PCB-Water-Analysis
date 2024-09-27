# global.R
library(shiny)
library(leaflet)
library(ggplot2) # Assuming you're using ggplot2 in your server logic

# Load data
wdc <- read.csv("WaterDataCongenerAroclor09072023.csv")

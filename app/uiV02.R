# uiV02.R
ui <- fluidPage(
  titlePanel("PCB Water Concentration Data Visualization"),
  selectInput("location_select", "Select Location:", 
              choices = c("All", unique(wdc$LocationName))),  # Include "All" option
  leafletOutput("map"),
  splitLayout(
    tableOutput("data"),
    plotOutput("plot", height = "300px", width = "700px"),
    cellWidths = c("35%", "65%")
  ),
  verbatimTextOutput("plot_text")
)

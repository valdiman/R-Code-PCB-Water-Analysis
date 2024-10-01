library(shiny)
library(leaflet)
library(ggplot2)
library(data.table)

# Load the dataset from the RDS file
#wdc <- readRDS("WaterDataCongenerAroclor09072023.rds")

# Define the Shiny server
# Define the server
server <- function(input, output, session) {
  # Convert numeric columns
  wdc$Latitude <- as.numeric(wdc$Latitude)
  wdc$Longitude <- as.numeric(wdc$Longitude)
  wdc$tPCB <- as.numeric(wdc$tPCB)
  
  # Filter the data based on the selected location
  filtered_data <- reactive({
    if (input$location_select == "All") {
      wdc
    } else {
      subset(wdc, LocationName == input$location_select)
    }
  })
  
  # Render the map
  output$map <- renderLeaflet({
    leaflet(filtered_data()) %>%
      addTiles() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude
      )
  })
  
  # Render the table
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      filtered_data_subset <- subset(wdc, SiteID == siteid)[, c("SiteID", "SampleDate", "tPCB")]
      
      # Sort the data by date
      filtered_data_subset <- filtered_data_subset[order(as.Date(filtered_data_subset$SampleDate)), ]
      
      # Get the number of samples
      num_samples <- nrow(filtered_data_subset)
      
      # Rename the SiteID column
      colnames(filtered_data_subset)[1] <- paste("SiteID (n =", num_samples, ")")
      
      # Ensure that the renaming is reflected correctly in the output
      return(filtered_data_subset)
    } else {
      return(data.frame())
    }
  })
  
  # Render the plot
  output$plot <- renderPlot({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      site_data <- subset(filtered_data(), SiteID == siteid)
      
      # Remove NAs
      site_data <- site_data[!is.na(site_data$SampleDate) & !is.na(site_data$tPCB), ]  # Filter out NAs
      
      if (nrow(site_data) == 0) {
        return(NULL)
      }
      
      site_data$week <- cut(as.Date(site_data$SampleDate), breaks = "week")
      data_agg <- aggregate(tPCB ~ week, data = site_data, mean, na.rm = TRUE)
      
      num_values <- nrow(data_agg)
      width <- ifelse(num_values > 5, 0.8, 0.2 + (num_values * 0.1))
      
      p <- ggplot(data_agg, aes(x = week, y = tPCB)) +
        geom_col(fill = "steelblue", width = width) +
        labs(x = NULL, y = paste("\u03A3", "PCB (pg/L)", sep = "")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 8)
        )
      
      max_value <- max(data_agg$tPCB, na.rm = TRUE)
      if (max_value >= 50000) {
        p <- p + scale_y_log10()
      }
      
      return(p)
    } else {
      return(NULL)
    }
  })
  
  output$plot_text <- renderPrint({
    if (!is.null(input$map_marker_click)) {
      cat("Plots are showing the mean data per week.\n")
      cat("If the maximum tPCB is too large (>50,000 pg/L), the y-axis changes to log10 scale.\n")
      cat("Source:")
    }
  })
  
  observe({
    leafletProxy("map", data = filtered_data()) %>%
      clearMarkers() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        layerId = ~SiteID,
        popup = ~paste(
          "Location Name: ", LocationName, "<br>",
          "SiteID: ", SiteID, "<br>",
          "Latitude: ", Latitude, "<br>",
          "Longitude: ", Longitude, "<br>",
          "Number of Samples: ", as.character(table(filtered_data()$SiteID)[as.character(SiteID)])
        )
      )
  })
}

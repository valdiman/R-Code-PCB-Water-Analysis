# Define the Shiny server
server <- function(input, output, session) {
  
  # Data.frame
  tPCB <- data.frame(
    SiteID = with(wdc, SiteID),
    LocationName = with(wdc, LocationName),
    date = with(wdc, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(wdc, as.numeric(Latitude)),
    Longitude = with(wdc, as.numeric(Longitude)),
    tPCB = with(wdc, as.numeric(tPCB))
  )
  
  # Filter the data based on the selected location
  filtered_data <- reactive({
    if (input$location_select == "All") {
      tPCB  # Return all data
    } else {
      subset(tPCB, LocationName == input$location_select)  # Filter by selected location
    }
  })
  
  # Render the map
  output$map <- renderLeaflet({
    leaflet(filtered_data()) %>%
      addTiles() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        label = ~as.character(SiteID)
      )
  })
  
  # Render the table
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      filtered_data <- subset(tPCB, SiteID == siteid)[, c("SiteID", "date", "tPCB")]
      filtered_data$date <- format(as.Date(filtered_data$date, format = "%m/%d/%y"), "%m-%d-%Y")
      colnames(filtered_data)[3] <- paste("\u03A3", "PCB ", "(pg/L)", sep = "")
      
      # Sort the data by date
      filtered_data <- filtered_data[order(as.Date(filtered_data$date, format = "%m-%d-%Y")), ]
      
      # Get the number of samples
      num_samples <- nrow(filtered_data)
      
      colnames(filtered_data)[1] <- paste("SiteID (n =", num_samples, ")")
      
      return(filtered_data)
    } else {
      return(data.frame())
    }
  })
  
  # Render the plot
  output$plot <- renderPlot({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      site_data <- subset(filtered_data(), SiteID == siteid)
      
      if (nrow(site_data) == 0) {
        # No data available for the selected SiteID
        return(NULL)
      }
      
      # Aggregate data by week and calculate the average of PCB values
      site_data$week <- cut(site_data$date, breaks = "week")
      data_agg <- aggregate(tPCB ~ week, data = site_data, mean)
      
      num_values <- nrow(data_agg)
      width <- ifelse(num_values > 5, 0.8, 0.2 + (num_values * 0.1))
      
      p <- ggplot(data_agg, aes(x = week, y = tPCB)) +
        geom_col(fill = "steelblue", width = width) +
        labs(x = NULL, y = paste("\u03A3", "PCB (pg/L)", sep = "")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12),  # Increase the font size of y-axis numbers
          axis.title.y = element_text(size = 14),  # Increase the font size of y-axis title
          axis.text = element_text(size = 8)
        )
      
      # Check if y-values are too large for linear scale
      max_value <- max(data_agg$tPCB)
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

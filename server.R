library(shiny)
library(dplyr)
library(ggplot2)
library(cgdsr)

source("utility_functions.R")

ggplot2::theme_set(theme_classic() +
    theme(axis.line.x = element_blank()) + 
    theme(axis.line.y = element_blank()))

colmutcat <- c("(germline)" = "black", "mutated" = "#1070b8")
alphamutcat <- c("(germline)" = 0.5, "mutated" = 1)
shapemutcat <- c("(germline)" = 1, "mutated" = 16)

function(input, output) {
  
  retrieved_tcga_data <- reactive({
    input$retrieve_button
    ids <- split_query_str(isolate(input$query_str))
    retrieve_tcga_data(ids)
  })
  # retrieved_tcga_data <- eventReactive(input$retrieve_button, {
  #   ids <- split_query_str(isolate(input$query_str))
  #   retrieve_tcga_data(ids)
  # })
  
  output$retrieved_genes <- renderText({
    paste(
      "Data retrieved for genes:", 
      paste(retrieved_tcga_data()$ids, collapse = ", "))
  })
  
  output$var_y_ui = renderUI({
    ids <- retrieved_tcga_data()$ids
    selectInput("var_y", "Gene on vertical axis", 
      choices = ids, selected = ids[1])
  })
  
  output$var_x_ui = renderUI({
    ids <- retrieved_tcga_data()$ids
    selectInput("var_x", "Gene on horizontal axes", 
      choices = ids, selected = ids[min(2, length(ids))])
  })
  
  assembled_graphics_data <- reactive({
    graphics_data <- retrieved_tcga_data()$data %>%
      mutate_(
        x_mut = paste0(input$var_x, "_mutations"), 
        x_gistic = paste0(input$var_x, "_gistic"), 
        x_rna = paste0(input$var_x, "_rna"), 
        y = paste0(input$var_y, "_rna")) %>%
      mutate(
        x_mutcat = 
          factor(x_mut == "(germline)",
            levels = c(TRUE, FALSE),
            labels = c("(germline)", "mutated"))) %>%
      select(x_mut, x_mutcat, x_gistic, x_rna, y)
    graphics_data
  })
  
  output$tab1 <- renderTable({
    tab1 <- assembled_graphics_data() %>%
      filter(!is.na(x_mut) & !is.na(y)) %>%
      select(x_mut) %>%
      table() %>%
      as.data.frame.table()
    names(tab1) <- c(paste0(input$var_x, ", AA change(s)"), "n")
    tab1
  })
  
  output$plot1 <- renderPlot({
    if (input$show_mut) {
      gg1 <- assembled_graphics_data() %>%
        filter(!is.na(x_mut) & !is.na(y)) %>%
        ggplot(aes(x = x_mut, y = y))
    } else {
      gg1 <- assembled_graphics_data() %>%
        filter(!is.na(x_mut) & !is.na(y)) %>%
        ggplot(aes(x = x_mutcat, y = y))
    }
    if (input$mark_mut) {
      gg1 <- gg1 +
        geom_point(aes(col = x_mutcat, alpha = x_mutcat, shape = x_mutcat), 
          position = position_jitter(h = 0,  w = 0.1)) + 
        geom_boxplot(col = "darkred", varwidth = TRUE,
          fill = "transparent", outlier.colour = "transparent") +
        scale_colour_manual(values = colmutcat, na.value = "black", guide = FALSE) + 
        scale_alpha_manual(values = alphamutcat, na.value = 1, guide = FALSE) + 
        scale_shape_manual(values = shapemutcat, na.value = 4, guide = FALSE) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(
          x = paste0(input$var_x, ", somatic point mutations"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"))
    } else {
      gg1 <- gg1 +
        geom_point(shape = 1, alpha = 0.5, 
          position = position_jitter(h = 0,  w = 0.1)) + 
        geom_boxplot(col = "darkred", varwidth = TRUE,
          fill = "transparent", outlier.colour = "transparent") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(
          x = paste0(input$var_x, ", somatic point mutations"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"))
    }
    plot(gg1)
  })
  
  output$plot2 <- renderPlot({
    gg2 <- assembled_graphics_data() %>%
      filter(!is.na(x_gistic) & !is.na(y)) %>%
      ggplot(aes(x = x_gistic, y = y)) 
    if (input$mark_mut) {
      gg2 <- gg2 +
        geom_point(aes(col = x_mutcat, alpha = x_mutcat, shape = x_mutcat), 
          position = position_jitter(h = 0,  w = 0.1)) + 
        geom_boxplot(col = "darkred", varwidth = TRUE,
          fill = "transparent", outlier.colour = "transparent") +
        scale_colour_manual(values = colmutcat, na.value = "black") + 
        scale_alpha_manual(values = alphamutcat, na.value = 1) + 
        scale_shape_manual(values = shapemutcat, na.value = 4) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(
          x = paste0(input$var_x, ", putative CNA (GISTIC)"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"),
          col = input$var_x, alpha = input$var_x, shape = input$var_x)
    } else {
      gg2 <- gg2 +
        geom_point(shape = 1, alpha = 0.5, 
          position = position_jitter(h = 0,  w = 0.1)) + 
        geom_boxplot(col = "darkred", varwidth = TRUE,
          fill = "transparent", outlier.colour = "transparent") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        labs(
          x = paste0(input$var_x, ", putative CNA (GISTIC)"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"))
    }
    plot(gg2)
  })
  
  output$plot3 <- renderPlot({
    gg3 <- assembled_graphics_data() %>%
      filter(!is.na(x_rna) & !is.na(y)) %>%
      ggplot(aes(x = x_rna, y = y)) 
    if (input$mark_mut) {
      gg3 <- gg3 +
        geom_point(aes(col = x_mutcat, alpha = x_mutcat, shape = x_mutcat)) + 
        scale_colour_manual(values = colmutcat, na.value = "black") + 
        scale_alpha_manual(values = alphamutcat, na.value = 1) + 
        scale_shape_manual(values = shapemutcat, na.value = 4) + 
        labs(
          x = paste0(input$var_x, ", mRNA expression (log2 RNA-seq)"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"), 
          col = input$var_x, alpha = input$var_x, shape = input$var_x)
    } else {
      gg3 <- gg3 +
        geom_point(shape = 1, alpha = 0.5) + 
        labs(
          x = paste0(input$var_x, ", mRNA expression (log2 RNA-seq)"), 
          y = paste0(input$var_y, ", mRNA expression (log2 RNA-seq)"))
    }
    if (input$smooth_method3 != "(none)")
      gg3 <- gg3 + geom_smooth(col = "darkred", method = input$smooth_method3)
    plot(gg3)
  })
  
  output$text3 <- renderText({ 
    r <- cor(
      assembled_graphics_data()$x_rna, 
      assembled_graphics_data()$y, 
      use = "complete.obs", method = "spearman")
    paste("Spearman's rank correlation coefficient:", format(r, digits = 2))
  })
  
}
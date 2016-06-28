# eQTL Browser
# Server.R
#
# User interface to browse eQTL data.
# Nick Burns
# June 2016

library(shiny)
library(RMySQL)
library(glida)
library(ggplot2)
library(data.table)
source("eqtl_library.R")


# iris: for test purposes only
data(iris)

shinyServer(function(input, output) {
    
    observeEvent(input$btn_query, {
        
        lcl_gene <- input$txt_query
        data <- browse(lcl_gene, input$radio_query)
        
        output$plt_panel_main <- renderPlot({
            ggplot(data, aes(x = build_37_pos, y = -log10(pvalue+1e-20))) +
                geom_point(colour = "dodgerblue", alpha = 0.8) +
                theme_minimal()
        })
        
        output$plt_panel_thumb <- renderPlot({
                ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) + geom_point()
        })
        
        output$tbl_eqtls <- renderDataTable({data})
    })
})
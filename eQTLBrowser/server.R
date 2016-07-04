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
            ggplot(dummy_data(lcl_gene, data_type = "qtls"), aes(x = build_37_pos, y = -log10(pvalue+1e-20))) +
                geom_point(colour = "dodgerblue", alpha = 0.8) +
                theme_minimal()
        })
        
        output$plt_panel_thumb <- renderPlot({
                
            ggplot(dummy_data(lcl_gene), aes(x = genotype, y = expression, group = genotype)) +
                geom_boxplot(colour = "dodgerblue", fill = "dodgerblue", alpha = 0.6) +
                geom_jitter(colour = "darkgrey", alpha = 0.2) +
                theme_minimal()
        })
        
        output$plt_panel_bottom <- renderPlot({
            
            ggplot(browse_expression(lcl_gene), aes(x = SMTSD, y = rpkm, group = SMTSD)) +
                geom_boxplot(aes(colour = smts, fill = smts), alpha = 0.5) +
                theme_minimal() +
                xlab("") +
                theme(axis.text.x = element_text(angle = 60, vjust = 0.7))
            
        })
        output$tbl_eqtls <- renderDataTable({data})
    })
})
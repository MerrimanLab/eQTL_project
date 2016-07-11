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


shinyServer(function(input, output) {
    
    dimTissue <- lookup_tissues()
    
    display_results <- function (lcl_gene, data_) {
        output$plt_panel_main <- renderPlot({
            display_eqtls(data_)
        })
        
        output$plt_panel_thumb <- renderPlot({
            
            ggplot(dummy_data(lcl_gene), aes(x = genotype, y = expression, group = genotype)) +
                geom_boxplot(colour = "darkgrey", fill = "darkgrey", alpha = 0.1) +
                geom_jitter(colour = "darkgrey", alpha = 0.1) +
                theme_minimal()
            
            
        })
        
        output$plt_panel_bottom <- renderPlot({
            
            display_expression(lcl_gene)
            
        })
        #output$tbl_eqtls <- renderDataTable({data_[order(pvalue, decreasing = FALSE)]})
        output$tbl_eqtls <- renderDataTable({
            snp_info <- if (!is.null(input$plot_click)) {
                nearPoints(data_, input$plot_click)
            } else NULL
            
            if (! is.null(snp_info)) {
                snp_info <- all_snp_info(snp_info)    
            }
            
            return (snp_info)
        })
    }
    
    observeEvent(input$btn_browse, {
        
        lcl_gene <- input$txt_gene_query
        data_ <- browse_qtls(lcl_gene, dimTissue)
        
        display_results(lcl_gene, data_)
        
    })
    
    observeEvent(input$btn_snp, {
        
        lcl_snp <- input$txt_snp_query
        snp_coords <- glida::queryUCSC(glida::updatePositions(lcl_snp))
        genes <- browse_snps(snp_coords)
        
        gene_list <- genes$gene_symbol
        names(gene_list) <- genes$gene_symbol
        
        output$ui_message <- renderUI({
            p(sprintf("The following genes display eQTLs for positions nearby %s:", lcl_snp))
        })
        output$ui_gene_list <- renderUI({
                radioButtons("rad_genes", label = "", choices = gene_list)
        })
    })
    
    observeEvent(input$rad_genes, {
        lcl_gene <- input$rad_genes
        data_ <- browse_qtls(lcl_gene, dimTissue)
        
        display_results(lcl_gene, data_)
    })
})
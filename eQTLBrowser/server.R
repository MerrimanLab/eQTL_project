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
library(gridExtra)
source("eqtl_library.R")

# for testing only:
options(shiny.maxRequestSize=100*1024^2) 

shinyServer(function(input, output) {
    
    dimTissue <- lookup_tissues()
    
    
    #### QTL-relevant code ####
    display_qtl_results <- function (lcl_gene, data_) {
        
        output$plt_panel_main <- renderPlot({
            display_eqtls(data_)
        })
        
        output$plt_panel_bottom <- renderPlot({
            
            display_expression(lcl_gene)
            
        })
        output$tbl_eqtls <- renderDataTable({data_[order(pvalue, decreasing = FALSE)]})
        
    }
    
    observeEvent(input$btn_browse, {
        
        lcl_gene <- input$txt_gene_qtl
        data_ <- browse_qtls(lcl_gene, dimTissue)

        display_qtl_results(lcl_gene, data_) 
        
    })
    
    observeEvent(input$btn_snp, {
        
        lcl_snp <- input$txt_snp_query
        snp_coords <- glida::queryUCSC(glida::updatePositions(lcl_snp))
        genes <- browse_by_snp(snp_coords)
        
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
        
        display_qtl_results(lcl_gene, data_)
    })
    
    #### GWAS : QTL functions ####
    observeEvent(input$btn_gwas, {
        
        file_ <- input$gwas_file
        
        chr_ <- as.integer(input$gwas_chr)
        start_ <- as.integer(input$gwas_start)
        end_ <- as.integer(input$gwas_end)
        
        gwas_data_ <- extract_gwas(file_$datapath, chr_, start_, end_)
        
        long_range_qtls <- qtl_network(chr_, start_, end_)
        
        output$plt_gwas <- renderPlot({
            display_gwas(gwas_data_) + xlim(start_ / 1000000, end_ / 1000000)
        })
        output$plt_qtl_network <- renderPlot({
            viz <- display_eqtls(long_range_qtls[-log10(pvalue) > input$sld_qtl_threshold], 
                          show_genes = FALSE, show_tissues = FALSE, alpha_pvalues = FALSE, show_title = FALSE)
            
            display_qtl_network(viz, long_range_qtls[-log10(pvalue) > input$sld_qtl_threshold], show_endpoint = FALSE) +
                xlim(start_ / 1000000, end_ / 1000000)
        })

    })
    
    # select GWAS peak and re-display the QTL network
    observeEvent(input$gwas_peak, {
        
        start_ <- as.integer(input$gwas_start)
        end_ <- as.integer(input$gwas_end)
        chr_ <- as.integer(input$gwas_chr)
        long_range_qtls <- qtl_network(chr_, start_, end_)
        
        viz <- display_eqtls(long_range_qtls[((-log10(pvalue) > input$sld_qtl_threshold) &
                                                  (build_37_pos > (input$gwas_peak$xmin * 1000000)) &
                                                  (build_37_pos < (input$gwas_peak$xmax * 1000000)))], 
                             show_genes = FALSE, show_tissues = FALSE, alpha_pvalues = FALSE, show_title = FALSE)
        
        output$plt_qtl_network <- renderPlot({
            display_qtl_network(viz, long_range_qtls[((-log10(pvalue) > input$sld_qtl_threshold) &
                                                          (build_37_pos > (input$gwas_peak$xmin * 1000000)) &
                                                          (build_37_pos < (input$gwas_peak$xmax * 1000000)))],
                                show_endpoint = FALSE) + xlim(start_ / 1000000, end_ / 1000000)
        })
    })
})
# eQTL Browser
# UI.R
#
# User interface to browse eQTL data.
# Nick Burns
# June 2016

library(shiny)

shinyUI(fluidPage(
    
    theme = "interface_styles.css",
    
    headerPanel(""),
    sidebarPanel(
        h2("eQTL Browser", class = "heading"),
        br(),
        hr(),
        
        p("Gene search", class = "boldtext"),
        p("Insert a gene name below to browse eQTLs and gene expression data.", class = "standardtext"),
        textInput("txt_gene_query", label = "", placeholder = "example: ABCG2"),
        br(),
        br(),
        actionButton("btn_browse", label = "Search", class = "button"),
        br(),
        hr(),
        
        conditionalPanel(
            condition = "input.conditionedPanels == 1",
            p("SNP search", class = "boldtext"),
            p("Insert a SNP name to identify the genes for which there are relevant eQTL results."),
            textInput("txt_snp_query", label = "", placeholder = "example: rs9930506"),
            br(),
            actionButton("btn_snp", label = "Search", class = "button"),
            br(),
            br(),
            uiOutput("ui_message"),
            uiOutput("ui_gene_list")
        ),
        conditionalPanel(
            condition = "input.conditionedPanels == 2",
            p("GWAS Summary Dataset", class = "boldtext")
        )
        
    ),
    mainPanel(
        tabsetPanel(
            
            # eQTL and gene expression panel
            tabPanel(
                
                h4("eQTL and Gene Expression"),
                br(),
                br(),
                tags$div(
                    tags$div(
                        plotOutput("plt_panel_main", click = "plot_click"),
                        class = "plot_main"
                    ),
                    tags$div(
                        uiOutput("sld_pval")
                    ),
                    tags$div(
                        plotOutput("plt_panel_bottom"),
                        class = "plot_bottom"
                    ), class = "plot_panel"
                ),
                br(),
                dataTableOutput("tbl_eqtls"),
                
                value = 1    # this is for the conditional sidebar panels
            ),
            
            # GWAS : eQTL panel
            tabPanel(
                h4("eQTL and GWAS"),
                
                value = 2
            ),
            
            # download datasets
            tabPanel(
                h4("Download Datasets")
            ),
            
            id = "conditionedPanels"    # this is the main id for the conditional sidebar panels
        )
    )
))
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
        
        conditionalPanel(
            condition = "input.conditionedPanels == 1",
            
            p("Gene search", class = "boldtext"),
            p("Insert a gene name below to browse eQTLs and gene expression data.", class = "standardtext"),
            textInput("txt_gene_qtl", label = "", placeholder = "example: ABCG2"),
            br(),
            br(),
            actionButton("btn_browse", label = "Search", class = "button"),
            br(),
            hr(),
            
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
            p("GWAS Summary Dataset", class = "boldtext"),
            br(),
            p("GWAS summary datasets can be uploaded below. Each file should contain the folowing
              columns: (CHR, POS, P).", class = "standardtext"),
            br(),
            fileInput("gwas_file", p("GWAS File: ", class = "boldtext")),
            hr(),
            p("QTL Threshold", class = "boldtext"),
            p("Adjust the -log10( pvalue ) threshold below. QTLs below this threshold will not be displayed.
              Note: setting this to zero will hide the QTL network.", class = "standardtext"),
            sliderInput("sld_qtl_threshold", label = "", min = 0, max = 10, value = 4),
            br(),
            hr(),
            p("GWAS / eQTL Search: ", class = "boldtext"),
            br(),
            p("Chromosome: ", class = "standardtext"),
            textInput("gwas_chr", label = "", placeholder = "example: 12"),
            br(),
            p("Region start: ", class = "standardtext"),
            textInput("gwas_start", label = "", placeholder = "example: 1000000"),
            br(),
            p("Region end: ", class = "standardtext"),
            textInput("gwas_end", label = "", placeholder = "example: 2000000"),
            br(),
            actionButton("btn_gwas", label = "Search", class = "button")
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
                br(),
                tags$div(
                    plotOutput("plt_gwas"),# brush = "gwas_peak"),
                    class = "plot_main"
                ),
                tags$div(
                    plotOutput("plt_qtl_network"),
                    class = "plot_bottom"
                ),
                
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
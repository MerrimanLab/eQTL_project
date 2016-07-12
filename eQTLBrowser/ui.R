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
        hr(),
        p("Some instructions will go here", class = "standardtext"),
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
    mainPanel(
        # TO DO: add hover functionality to this
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
        dataTableOutput("tbl_eqtls")
    )
))
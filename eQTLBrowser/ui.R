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
        
        p("Search eQTLs for:", class = "boldtext"),
        textInput("txt_query", label = "", placeholder = "example: ABCG2"),
        radioButtons("btn_query", label = "", choices = list(snp = "snp", gene = "gene")),
        br(),
        
        actionButton("btn_query", label = "GO", class = "button")
    ),
    mainPanel()
))
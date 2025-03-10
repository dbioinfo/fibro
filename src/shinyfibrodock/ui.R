library(tidyverse)
library(shiny)
library(shinythemes)
library(DT)
library(shinyjs)
library(plotly)

#setwd("Postdoc/Hindle/dylan/fibro/")

# Define UI
ui <- fluidPage(navbarPage(title = "Multispecies Expression Analysis Dashboard",
  # Apply the flatly theme
  theme = shinytheme("flatly"),
  
  # Create tabsetPanel with 4 tabs
  tabsetPanel(
    # ---- Summary Tab ----
    tabPanel("Summary", 
             sidebarLayout(
               sidebarPanel(
                 # Button to load data
                 actionButton("loadData", "Load Data")
               ),
               mainPanel(
                 h3("Differential Expression Summary"),
                 DTOutput("summaryTable")
               )
             )
    ),
    
    # ---- Expression Plotter Tab ----
    tabPanel("Expression Plotter",
             sidebarLayout(
               sidebarPanel(
                 # Selectize input for species
                 uiOutput("speciesSelectExp"),
                 
                 # Selectize input for gene
                 uiOutput("geneSelect"),
                 
                 # Radio buttons for experimental groups
                 radioButtons("iexpGroup", "Experimental Group:",
                              choices = c("Glucose"="glucose", "Hypoxia"="hypoxia", "Temperature"="temperature")),
                 
                 # X axis groupings 
                 radioButtons("iexpX", "Choose X Axis",
                              choices=c("Individual"='individual', "Treatment"='treatment'), 
                              selected = 'individual'),
                 
                 # Refresh button
                 actionButton("refreshExpressionPlot", "Refresh Plot")
               ),
               mainPanel(
                 h3("Gene Expression Plot"),
                 plotlyOutput("expressionPlot", height = "600px")
               )
             )
    ),
    
    # ---- DE Genes Tab ----
    tabPanel("DE Genes",
             sidebarLayout(
               sidebarPanel(
                
                #select species 
                uiOutput("speciesSelectDE"),
                
                 
                #select one of the three experimental groups
                radioButtons('ideGroup', "Choose Experimental Group",
                             choices= c("Glucose"="glucose", "Hypoxia"="hypoxia","Temperature"="temperature")),
                 
                #conditional panel that chooses contrast based on previous value
                conditionalPanel(
                  condition = "input.ideGroup == 'glucose'",
                  radioButtons("glucoseContrast", 
                               "Select Contrast", 
                               choices = c("2.5mM vs 8mM" = "2.5v8", 
                                           "30mM vs 8mM" = "30v8", 
                                           "2.5mM vs 30mM" = "2.5v30"))
                ),
                conditionalPanel(
                  condition = "input.ideGroup == 'hypoxia'",
                  radioButtons("hypoxiaContrast", 
                               "Select Contrast", 
                               choices = c("6H vs 8mM" = "6v8", 
                                           "24H vs 8mM" = "24v8", 
                                           "6H vs 24H" = "6v24"))
                ),
                conditionalPanel(
                  condition = "input.ideGroup == 'temperature'",
                  radioButtons("tempContrast", 
                               "Select Contrast", 
                               choices = c("32°C vs 8mM" = "32v8", 
                                           "41°C vs 8mM" = "41v8", 
                                           "32°C vs 41°C" = "32v41"))
                ),
                 
                
                #show guides?
                radioButtons("deGuides","Show cutoffs?",
                             choices=c("Yes"=T, "No"=F),
                             selected=T),
                
                #color?
                radioButtons("deColor", "Choose Color Mapping",
                             choices = c("Expression Quartile"="expr", "Significance"='sig'),
                             selected='sig'),
                
                # species for table of DE genes
                 
                 uiOutput("speciesSelectDETable"),
                  
                # Refresh button
                actionButton("refreshDEPlot", "Refresh Plot")

               ),
               mainPanel(
                 h3("DE Visualization"),
                 plotlyOutput("dePlot", height = "600px"),
                 h4("DE Table"),
                 DTOutput("deGenesTable")
               )
             )
    ),
    
    # ---- KEGG Pathway Tab ----
    tabPanel("KEGG Pathway",
             sidebarLayout(
               sidebarPanel(
                 p("KEGG pathway analysis options will be added in future updates.")
               ),
               mainPanel(
                 h3("KEGG Pathway Analysis"),
                 p("This section displays enriched KEGG pathways based on differentially expressed genes."),
                 DTOutput("keggTable")
               )
             )
    )
  )
)
)



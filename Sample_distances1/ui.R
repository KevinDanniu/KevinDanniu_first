####library R packages####
library(shiny)
library(plotly)
library(shinythemes)
library(shinycssloaders)

####File's size####
options(shiny.maxRequestSize=90*1024^2)

####UI####
#----UI----#
shinyUI(fluidPage(
  #use a beautiful shiny theme
  theme = shinytheme("united"),
  
  #Application title
  titlePanel("(DeSeq2)Heatmap Analysis"),
  
  #Sidebar Layout
  sidebarLayout(
    #Sidebar
    sidebarPanel(
      width = 2,
      #Upload file
      HTML('<h3 style="margin-top: 0px;">Data upload</h3>'),
      'Warnings:',br(),    
      '1) The number of groups in the sample was less than 11.',br(),
      '2) The attribution "treat" in gse must exist. ("treat" means group.)',
      fileInput("upload", NULL, multiple = FALSE, accept = c(".Rds",".RDS")),
      
      #More information
      HTML('<h3 style="margin-top:0 px;">More information</h3>'),
      "There are currently",verbatimTextOutput("count"),
      "session(s) connected to this app.",
      hr(),
      
      #Preprocess Select
      HTML('<h3 style="margin-top: 0px;">Preprocess data</h3>'),
      selectInput("preprocess",
                  label = "Preprocess data",
                  choices = c("log2(exp + 1)", "VST", "rlog"),
                  selected = "log2(exp + 1)"),
      
      #Heatmap Select
      HTML('<h3 style="margin-top: 0px;">Distance measure</h3>'),
      selectInput("distancemeasure",
                  label = "Distance measure",
                  choices = c("Euclidean", "Poisson"),
                  selected = "Euclidean"),
      
      #Heatmap Size Select
      HTML('<h3 style="margin-top: 0px;">Heatmap size</h3>'),
      "Default: initial Width is 800 px; Height is 600 px.",
      textInput("heatmapwidth",
                label = "Heatmap width (px)",
                value = 800),
      textInput("heatmapheight",
                label = "Heatmap height (px)",
                value = 600),
      
      #Dimensionality-reduction Select
      HTML('<h3 style="margin-top: 0px;">Dimensionality-reduction measure</h3>'),
      "Warings: PCA based on normalized data; GLM-PCA based on count data;",
      "MDS based on normalized data.",
      selectInput("drmeasure",
                  label = "Dimensionality-reduction measure",
                  choices = c("PCA", "GLM-PCA", "MDS"),
                  selected = "PCA"),
      
      #Dimensionality-reduction Plot Size Select
      HTML('<h3 style="margin-top: 0px;">Dimensionality-reduction plot size</h3>'),
      "Default: initial Width is 600 px; Height is 400 px.",
      textInput("drwidth",
                label = "PCA plot width (px)",
                value = 600),
      textInput("drheight",
                label = "PCA plot height (px)",
                value = 400),
      radioButtons('format', 'Plot output format',
                   choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'),
                   inline = T)
    ),
    
    #Main Panel
    mainPanel(
      #separate panel
      tabsetPanel(
        id = "tabs",
        tabPanel("sample Table", h4("Experiment table of Samples"),
                 dataTableOutput("sampletable")),
        tabPanel("Normalized data Table", h4("Table of Normalized Data"),
                 dataTableOutput("Normalizetable"),
                 absolutePanel(
                   top = 900,
                   downloadButton('dataDown',label="Download Table of Normalized Data")
                 )),
        tabPanel("Heatmap of Samples", h4("Similarity analysis for Samples"),
                 # dataTableOutput("Distancetable"),
                 plotOutput("heatmap"),
                 absolutePanel(
                   top = 900,
                   downloadButton('heatmapDown',label="Download Plot")
                 )),
        tabPanel("PCA plot of Samples", h4("Similarity analysis for Samples"),
                 plotOutput("drplot"),
                 absolutePanel(
                   top = 900,
                   downloadButton('PCADown',label="Download Plot")
                 ))
      )
    )
  )
))
library(rsconnect)
library(DT)
library(shiny)
library(shinyWidgets)
library(matrixStats)
library(plotly)
library(openxlsx)
library(rlist)
library(shinyEffects)
library(shinycssloaders)
library(scales)
library(tidyr)
library(plyr)
library(EnvStats)
library(stringr)
library(sortable)
library(shinybusy)
library(reshape2)

params <- list(allsc_genes        = "Files/sc_P_berghei_averaged.csv",
               allk_genes         = "Files/bulk_H_sapiens_averaged.csv",
               dotplot_genes      = "Files/sc_P_berghei_dotplot.csv",
               UMAP               = "Files/sc_P_berghei_UMAP.csv",
               Counts             = "Files/sc_P_berghei_counts.csv"
)

source("helper_module.R")

sc_genes <- read.csv2(params$allsc_genes, stringsAsFactors = FALSE)

human_genes <- read.csv2(params$allk_genes, stringsAsFactors = FALSE)

sc_dot_plot <- read.csv2(params$dotplot_genes, stringsAsFactors = FALSE)

UMAP_sc <- read.csv2(params$UMAP, stringsAsFactors = FALSE)

Sc_counts <- as.matrix(read.csv2(params$Counts, sep = ";", stringsAsFactors = FALSE))

################################################################################
##                                   UI                                       ##
################################################################################

ui <- navbarPage(
  span(HTML(paste0(tags$sup(tags$i("X"), style = "font-size : 18px; color:#383D3B ; font-weight:bold"),"SPOT")), style = " color: #b22222; font-size : 32px; font-weight: bold"),
  id = "navbar", selected = "value", 
  tags$style(HTML("  *{font-family: Verdana;}
                      .navbar { background-color: #f1f1f1}
                      .navbar-default .navbar-nav > li > a {color:#c00000;}
                      .navbar-default .navbar-nav > .active > a,
                      .navbar-default .navbar-nav > .active > a:focus,
                      .navbar-default .navbar-nav > .active > a:hover {color: white;background-color: #c00000;}
                      .navbar-default .navbar-nav > li > a:hover {color: white;background-color:#c00000;text-decoration:underline;text-color:white}
                      .button_DEA .bttn-danger{background-color: #c00000; border-color: black; margin-bottom: 10px }
                      .button_dd .bttn-danger{color: white; border-color: grey; background-color: #c00000; margin-bottom: 10px}
                      .text_about {color:black; font-size: 14px;text-align: justify; margin-top: 10px; margin-bottom: 20px}
                      .text_cite {color:black; font-size: 14px; font-weight: bold ;text-align: justify}
                      ")),
  tabPanel(span("SPOT expression profiles", style = "font-size : 20px; font-weight: bold "), value = "value",
           
           fluidPage(
             setShadow(id = "Component2a"),
             setShadow(id = "Component2b"),
             
             conditionalPanel(
               condition = "input.dataSwitch2 == 'Plasmodium (single cell)'",
               withSpinner(plotlyOutput("Component2a"), type = "2", color.background = "white", color = "#c00000"),
             ),
             
             conditionalPanel(
               condition = "input.dataSwitch2 == 'Human organs (bulk)'",
               withSpinner(plotlyOutput("Component2b"), type = "2", color.background = "white", color = "#c00000"),
               
             ),
             chooseSliderSkin("Flat"),
             
             setSliderColor(rep("#c00000", 100), c(1:100)),
             
             hr(),
             
             fluidRow(
               
               column(2, 
                      
                      radioGroupButtons(
                        inputId = "Algos",
                        label = "Choose algorithm",
                        choices = c("SPOT", "Correlation"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      conditionalPanel(
                        condition = "input.dataSwitch2 == 'Human organs (bulk)' ",
                        tags$h5("Subset samples", style = "font-weight: bold "),
                        #br(),
                        div(class = "button_dd",
                            dropdown(
                              
                              bucket_list(
                                header = "Choose developmental timepoints per Drag n Drop.\n wpc -> week post conception",
                                add_rank_list(input_id = "bucket_in",
                                              text = "Drag from here",
                                              labels = c("4wpc", "10wpc", "20wpc", "infant", "toddler", "school", 
                                                         "teenager", "youngAdult", "youngMidAge", "olderMidAge", "senior")
                                              ),
                                add_rank_list( input_id = "bucket_out",
                                               text = "to here",
                                               labels = c("newborn")
                                               ),
                                orientation = "horizontal"
                              ),
                              
                              style = "bordered", icon = icon("bars"), label = "Subset",
                              status = "danger", width = "1000px", 
                              height = "200px",
                              animate = animateOptions(
                                enter = animations$fading_entrances$fadeInLeftBig,
                                exit = animations$fading_exits$fadeOutRightBig
                              )
                            )),
                      ),
                      
                      pickerInput(
                        inputId = "dataSwitch2",
                        label = "Select Dataset", 
                        choices = c("Plasmodium (single cell)", "Human organs (bulk)")
                      ),
                    
               ),
               
               conditionalPanel(
                 condition = "input.dataSwitch2 == 'Plasmodium (single cell)' ",
                 slidersUI("1", colnames(sc_genes[4:13])),
               ),
               
               conditionalPanel(
                 condition = "input.dataSwitch2 == 'Human organs (bulk)'",
                 uiOutput("sliders_K"),
               ),
               
               column(2,
                      
                      conditionalPanel(
                        condition = "input.dataSwitch2 == 'Plasmodium (single cell)' ",
                        selectInput("spec_gene2", "Pick genes to visualize:", choices = c("a", "b"), NULL, multiple = TRUE), #sc_genes$PB_ID
                      ),
                      
                      radioGroupButtons(
                        inputId = "Radio2",
                        label = "Alter Visualization",
                        choices = c("Table", "Bar chart", "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      downloadButton('downloadData1', 'Download table as .xlsx')
                      
               )))),
  
  tabPanel(span("Differential Expression Analysis", style = "font-size : 20px; font-weight: bold"),
           fluidPage(
             
             plotlyOutput("Component3"),
             
             hr(),
             
             fluidRow(
               
               column(2, 
                      
                      radioGroupButtons(
                        inputId = "algorithm",
                        label = "Choose algorithm",
                        choices = c("Wilcox", "MAST", "DESeq2"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
                      downloadButton('downloadData2', 'Download table as .xlsx')
                      
               ),
               slidersUI("5", colnames(sc_genes[4:13])),
               column(2,
                      
                      div(class = "button_DEA",
                          actionBttn("action_DEA", label = "Start DEA!", icon = NULL, style = "material-flat", color = "danger",
                                     size = "lg", block = FALSE, no_outline = TRUE
                          )),
                      
                      radioGroupButtons(
                        inputId = "Radio3",
                        label = "Alter Visualization",
                        choices = c("Table", 
                                    "Bar chart", "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      
               ),
             ))),
  
  tabPanel(span("Compare profiles", style = " font-size : 20px; font-weight: bold "), 
           
           fluidPage(
             setShadow(id = "Component1"),
             setShadow("tab_gen"),
             withSpinner(plotlyOutput("Component1"), type = "2", color.background = "white", color = "#c00000"),
             fluidRow(
               column(12,
                      DTOutput("tab_gen"))
             )
           )
  ),
  
  tabPanel(span("Upload your own data", style = "font-size : 20px; font-weight: bold "),
           fluidPage( 
             
             withSpinner(plotlyOutput("Component4"), type = 2, color = "#c00000", color.background = "white"),
             
             hr(),
             
             fluidRow(
               column(2,
                      
                      fileInput('file1', 'Choose file to upload',
                                accept = c('text/csv','text/comma-separated-values',
                                           'text/tab-separated-values', 'text/plain',
                                           '.csv','.tsv', 'text/xlsx','.xlsx'),
                                ),
                      ),
               uiOutput("reactive_sliders"),
               column(2,
                      radioGroupButtons(
                        inputId = "Radio5",
                        label = "Alter Visualization",
                        choices = c("Table", 
                                    "Bar chart",
                                    "Dot plot"),
                        individual = TRUE,
                        checkIcon = list(
                          yes = tags$i(class = "fa fa-circle", 
                                       style = "color: #c00000"),
                          no = tags$i(class = "fa fa-circle-o", 
                                      style = "color: #c00000"))
                      ),
                      downloadButton('downloadData3', 'Download table as .xlsx')
               )
             ))),
  
  tabPanel(span("About", style = "font-size : 20px; font-weight: bold "),
           setShadow(id = "UMAP"),
           fluidRow(
             column(4, offset = 1,
                    div(class = "text_about",
                        "The increasing number of single cell and bulk RNAseq data sets describing 
                        complex gene expression profiles in different organisms, organs or cell types
                        calls for an intuitive tool allowing rapid comparative analysis. Here we present 
                        Swift Profiling Of Transcriptomics (SPOT) as a web tool that allows not only 
                        differential expression analysis but also fast selection of genes fitting transcription
                        profiles of interest. Based on a heuristic approach the spot algorithm ranks the genes 
                        according to their proximity to the user-defined gene expression profile of interest. 
                        The best hits are visualized as a table, bar chart or dot plot and can be exported as 
                        an Excel file. While the tool is generally applicable, we tested it on RNAseq data from
                        malaria parasites that undergo multiple stage transformations during their complex life 
                        cycle as well as on data from multiple human organs during development. SPOT should enable
                        non-bioinformaticians to easily analyse their own and any available dataset. "
                    ),
                    div(class = "text_cite",
                    "Elias Farr, Julia M Sattler, Friedrich Frischknecht, SPOT: a web-tool enabling Swift Profiling Of Transcriptomics"
                    )
             ),
             column(5, offset = 1,
                    plotlyOutput("UMAP")
                    )
             )
           )
  )
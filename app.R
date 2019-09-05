#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)

# load data
dataTab <- readRDS("data/smallHumanCellTypesDf.rds")
speciesDataTab <- readRDS("data/speciesMeltomeDf.rds")

# Define UI for application that draws a histogram
ui <- navbarPage("Meltome Atlas",
                 tabPanel("Home",
                 
                          titlePanel("Welcome to Meltome Atlas!"),
                          
                          sidebarLayout(
                            
                            sidebarPanel(
                              "This R-Shiny App allows interactive exploration of the dataset published by Jarzab et al., (2019). 
                              
                              Have fun!"
                              ),
                            
                            mainPanel(
                              "The different tabs lead you to exploring the either the comparison of thermal proteome profiles across various species or across different human cell lines and cell types."
                            )
                          )
                          ),
                 tabPanel("Cross-species Meltome",
                          titlePanel("Plot protein of interest in different species"),
                          
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Select species and a protein"),
                              
                              selectInput(inputId = "selectSpecies1", label = NULL, 
                                          choices = c("Escherichia coli lysate", 
                                                      "Escherichia coli cells", 
                                                      "Saccharomyces cerevisiae lysate",
                                                      "Arabidopsis thaliana seedling lysate",
                                                      "Drosophila melanogaster SII lysate",
                                                      "Caenorhabditis elegans lysate",
                                                      "Danio rerio Zenodo lysate",
                                                      "Homo sapiens A549 lysate",
                                                      "Homo sapiens Jurkat lysate",
                                                      "Mus musculus BMDC lysate"), 
                                          selected = "Drosophila melanogaster SII lysate"),
                              
                              selectInput("selectSpecies2", NULL, 
                                          choices = c("Escherichia coli lysate", 
                                                      "Escherichia coli cells", 
                                                      "Saccharomyces cerevisiae lysate",
                                                      "Arabidopsis thaliana seedling lysate",
                                                      "Drosophila melanogaster SII lysate",
                                                      "Caenorhabditis elegans lysate",
                                                      "Danio rerio Zenodo lysate",
                                                      "Homo sapiens A549 lysate",
                                                      "Homo sapiens Jurkat lysate",
                                                      "Mus musculus BMDC lysate",
                                                      ""), 
                                          selected = "Mus musculus BMDC lysate"),
                              
                              selectInput("selectSpecies3", NULL, 
                                          choices = c("Escherichia coli lysate", 
                                                      "Escherichia coli cells", 
                                                      "Saccharomyces cerevisiae lysate",
                                                      "Arabidopsis thaliana seedling lysate",
                                                      "Drosophila melanogaster SII lysate",
                                                      "Caenorhabditis elegans lysate",
                                                      "Danio rerio Zenodo lysate",
                                                      "Homo sapiens A549 lysate",
                                                      "Homo sapiens Jurkat lysate",
                                                      "Mus musculus BMDC lysate"), 
                                          selected = "Homo sapiens Jurkat lysate"),
                              
                              textInput("gene_name1", 
                                        label = NULL,
                                        value = "MTOR"),
                                        # placeholder = 'Gene name (e.g. MTOR)'),
                              
                              textInput("gene_name2", 
                                        label = NULL,
                                        placeholder = 'Gene name (e.g. CDK1)'),
                              
                              textInput("gene_name3", 
                                        label = NULL,
                                        placeholder = 'Gene name (e.g. CCNB2)')
                              
                            ),
                            
                            mainPanel(
                              plotOutput("species_melt_curve_plot1")
                            )
                          )),
                 tabPanel("Human cell line Meltome",
                          titlePanel("Plot protein of interest in different human cell lines"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Select cell lines and a protein"),
      
      selectInput(inputId = "selectCell1", NULL, 
                  choices = c("K562", "Jurkat", "HepG2", 
                              "HaCaT", "U937", "K562", 
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"), 
                  selected = "Jurkat"),
      
      selectInput("selectCell2", NULL, 
                  choices = c("K562", "Jurkat", "HepG2", 
                              "HaCaT", "U937", "K562", 
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"), 
                  selected = "U937"),
      
      selectInput("selectCell3", NULL, 
                  choices = c("K562", "Jurkat", "HepG2", 
                              "HaCaT", "U937", "K562", 
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"), 
                  selected = "pTcells"),
      
      textInput("protein", 
                label = NULL,
                value = "JAK1")
                # placeholder = 'Gene name (e.g. ABL1)')
      
    ),
    
    mainPanel(
      plotOutput("melt_curve_plot1")
    )
  ),
  
  titlePanel("Plot different proteins of interest in a desired cell line"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Select proteins and a cell line"),
      
      selectInput("selectCell", "Select a cell line", 
                  choices = c("K562", "Jurkat", "HepG2", 
                              "HaCaT", "U937", "K562", 
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"), 
                  selected = "Jurkat"),
      
      
      textInput("protein1", 
                label = NULL,
                #placeholder = 'Protein (e.g. ABL1)',
                value = "ABL1"),
      
      textInput("protein2", 
                label = NULL,
                value = "BCR"),
                # placeholder = 'Gene name (e.g. ABL1)'),
      
      textInput("protein3", 
                label = NULL,
                placeholder = 'Gene name (e.g. JAK1)')
    ),
    
    mainPanel(
      plotOutput("melt_curve_plot2")
    )
  )
), 
# tabPanel("Human body fluid Meltome"),
tabPanel("Download",
         sidebarLayout(
           
           # Sidebar panel for inputs ----
           sidebarPanel(
             
             # Input: Choose dataset ----
             selectInput("dataset", "Choose a dataset:",
                         choices = c("cross-species", "human")),
             
             # Button
             downloadButton("downloadData", "Download")
             
           ),
           
           # Main panel for displaying outputs ----
           mainPanel(
             
             # tableOutput("table")
             
           ))
)
)

# Server logic ----
server <- function(input, output) {
  
  output$melt_curve_plot1 <- renderPlot({
    
    ggplot(filter(dataTab, cellLine %in% c(input$selectCell1, input$selectCell2, 
                                           input$selectCell3), 
                  Protein_ID == input$protein), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = cellLine), size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "fraction non-denatured") +
      stat_summary(aes(color = cellLine), geom="line", fun.y="median") +
      theme_bw() + 
      theme(text = element_text(size = 20))
    
  })
  
  output$melt_curve_plot2 <- renderPlot({
    
    ggplot(filter(dataTab, cellLine == input$selectCell, 
                  Protein_ID %in% c(input$protein1, input$protein2,
                                    input$protein3)), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = Protein_ID), size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "fraction non-denatured") +
      stat_summary(aes(color = Protein_ID), geom="line", fun.y="median") +
      theme_bw() + 
      theme(text = element_text(size = 20))
    
  })
  
  output$species_melt_curve_plot1 <- renderPlot({
    
    ggplot(filter(speciesDataTab, run_name %in% c(input$selectSpecies1, input$selectSpecies2, 
                                           input$selectSpecies3), 
                  protein %in% c(input$gene_name1, input$gene_name2,
                                 input$gene_name3)), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = run_name, shape = protein), 
                 size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "fraction non-denatured") +
      stat_summary(aes(color = run_name, linetype = protein), geom="line", fun.y="median") +
      scale_color_discrete("species") +
      theme_bw() + 
      theme(text = element_text(size = 20))
    
  })
  
  # Reactive value for selected dataset ----
  datasetInput <- reactive({
    switch(input$dataset,
           "cross-species" = speciesDataTab,
           "human" = dataTab)
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
}

# Run app ----
shinyApp(ui, server)

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
library(ggplot2)

# load data
humanDataTab <- reactive(readRDS("/Users/kurzawa/repos/meltomeAtlasApp/data/humanCellsUpdated.rds"))
speciesDataTab <- reactive(readRDS("/Users/kurzawa/repos/meltomeAtlasApp/data/fullSpeciesDat.rds"))
humanCellsTmSpread <- reactive(readRDS("/Users/kurzawa/repos/meltomeAtlasApp/data/humanCellsTmSpread.rds"))
all_human_gene_names <- readRDS("/Users/kurzawa/repos/meltomeAtlasApp/data/uniqueHumanGeneNames.rds")
#all_species_gene_names <- readRDS("/home/nils/projects/meltomeatlasapp/data/uniqueSpeciesGeneNames.rds")
  
# Define UI for application that draws a histogram
ui <- navbarPage("Meltome Atlas",
                 tabPanel("Home",
                 
                          titlePanel("Welcome to Meltome Atlas!"),
                          
                          sidebarLayout(
                            
                            sidebarPanel(
                              tagList(
                                fluidRow("This R-Shiny App allows interactive exploration of the dataset published by", a("Jarzab et al., (2020), Nat. Methods.", href="https://www.nature.com/articles/s41592-020-0801-4")),
                                tags$br(),
                                fluidRow("The different tabs lead you to exploring either the comparison of thermal proteome profiles across various species or across different human cell lines and cell types."),
                                tags$br(),
                                fluidRow("Have fun!"))
                              ),
                            
                            mainPanel(
                              img(
                                src = "meltome_atlas_fig1a.png",
                                align = "center",
                                width = "743px", height = "500px"
                              )
                            )
                          )
                          ),
                 tabPanel("Cross-species Meltome",
                          titlePanel("Plot proteins of interest in different species"),
                          
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
                                                      "Mus musculus liver lysate",
                                                      "Geobacillus stearothermophilus NCA26 lysate",
                                                      "Thermus thermophilus HB27 lysate",
                                                      "Thermus thermophilus HB27 cells", 
                                                      "Picrophilus torridus DSM9790 lysate",
                                                      "Homo sapiens Jurkat NP40 lysate",
                                                      "Homo sapiens K562 PBS lysate",
                                                      "Homo sapiens Jurkat cells",
                                                      "Homo sapiens K562 cells", 
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
                                                      "Mus musculus liver lysate",
                                                      "Geobacillus stearothermophilus NCA26 lysate",
                                                      "Thermus thermophilus HB27 lysate",
                                                      "Thermus thermophilus HB27 cells", 
                                                      "Picrophilus torridus DSM9790 lysate",
                                                      "Homo sapiens Jurkat NP40 lysate",
                                                      "Homo sapiens K562 PBS lysate",
                                                      "Homo sapiens Jurkat cells",
                                                      "Homo sapiens K562 cells", 
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
                                                      "Mus musculus liver lysate",
                                                      "Geobacillus stearothermophilus NCA26 lysate",
                                                      "Thermus thermophilus HB27 lysate",
                                                      "Thermus thermophilus HB27 cells", 
                                                      "Picrophilus torridus DSM9790 lysate",
                                                      "Homo sapiens Jurkat NP40 lysate",
                                                      "Homo sapiens K562 PBS lysate",
                                                      "Homo sapiens Jurkat cells",
                                                      "Homo sapiens K562 cells", 
                                                      "Mus musculus BMDC lysate"), 
                                          selected = "Homo sapiens Jurkat cells"),
                              
                              
                              textInput("gene_name1",
                                        label = NULL,
                                        value = "MTOR"),
                                        # placeholder = 'Gene symbol (e.g. Mtor)'),
                              # selectizeInput(
                              #   "gene_name1", label = NULL, choices = all_species_gene_names,
                              #   options = list(create = TRUE),
                              #   selected = "Mtor"
                              # ),
                              textInput("gene_name2",
                                        label = NULL,
                                        placeholder = 'Gene symbol (e.g. CDK1)'),
                              
                              textInput("gene_name3",
                                        label = NULL,
                                        placeholder = 'Gene symbol (e.g. CCNB2)')
                              
                              # selectizeInput(
                              #   "gene_name2", label = NULL, choices = all_species_gene_names,
                              #   options = list(create = TRUE),
                              #   selected = ""
                              # ),
                              # 
                              # selectizeInput(
                              #   "gene_name3", label = NULL, choices = all_species_gene_names,
                              #   options = list(create = TRUE),
                              #   selected = ""
                              # )
                              
                            ),
                            
                            mainPanel(
                              plotOutput("species_melt_curve_plot1"),
                              fluidRow(
                                column(12, 
                                       tableOutput('table')
                                )
                            )
                          ))),
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
      
      
      selectizeInput(
        "protein", label = NULL, choices = all_human_gene_names,
        options = list(create = TRUE),
        selected = "JAK1"
      )
    #   textInput("protein",
    #             label = NULL,
    #             value = "JAK1")
    #             # placeholder = 'Gene symbol (e.g. ABL1)')
    ),
    
    mainPanel(
      plotOutput("melt_curve_plot1"),
      fluidRow(
        column(12, 
               tableOutput('table2')
        )
      )
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
      
      
      selectizeInput(
        "protein1", label = NULL, choices = all_human_gene_names,
        options = list(create = TRUE),
        selected = "AAAS"
      ),
      # textInput("protein1",
      #           label = NULL,
      #           #placeholder = 'Protein (e.g. ABL1)',
      #           value = "ABL1"),
      
      # textInput("protein2",
      #           label = NULL,
      #           value = "BCR"),
      #           # placeholder = 'Gene symbol (e.g. ABL1)'),
      selectizeInput(
        "protein2", label = NULL, choices = all_human_gene_names,
        options = list(create = TRUE),
        selected = "ABL1"
      ),
      
      # textInput("protein3",
      #           label = NULL,
      #           placeholder = 'Gene symbol (e.g. JAK1)')
      
      selectizeInput(
        "protein3", label = NULL, choices = all_human_gene_names,
        options = list(create = TRUE),
        selected = "BCR"
      )
      
    ),
    
    mainPanel(
      plotOutput("melt_curve_plot2"),
      fluidRow(
        column(12, 
               tableOutput('table3')
        )
      )
    )
  ),
  titlePanel("Scatterplot of protein melting points in desired cell lines"),

  sidebarLayout(
    sidebarPanel(
      helpText("Select cell lines"),

      selectInput("selectCellTm1", "Select a cell line",
                  choices = c("K562", "Jurkat", "HepG2",
                              "HaCaT", "U937", "K562",
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"),
                  selected = "Jurkat"),


      selectInput("selectCellTm2", "Select a cell line",
                  choices = c("K562", "Jurkat", "HepG2",
                              "HaCaT", "U937", "K562",
                              "HAOEC", "pTcells", "HEK293T",
                              "HL60", "colon_cancer_spheroids"),
                  selected = "HepG2")

    ),

    mainPanel(
      plotOutput("tm_scatter")
        )
      )

), 
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
             
           ))
         
    )

)

# Server logic ----
server <- function(input, output) {
  
  output$melt_curve_plot1 <- renderPlot({
    
    ggplot(filter(humanDataTab(), cell_line_or_type %in% c(input$selectCell1, input$selectCell2, 
                                           input$selectCell3), 
                  toupper(gene_name) == toupper(input$protein)), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = cell_line_or_type), size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "Fraction non-denatured") +
      stat_summary(aes(color = cell_line_or_type), geom="line", fun.y="median") +
      scale_color_discrete("Cell line/type") +
      theme_bw() + 
      theme(text = element_text(size = 20))
    
  })
  
  output$table2 <- renderTable(filter(humanDataTab(), cell_line_or_type %in% c(input$selectCell, input$selectCell2, 
                                                                                   input$selectCell3), 
                                          toupper(gene_name) == toupper(input$protein)) %>% 
                                     group_by(cell_line_or_type) %>% 
                                     summarize(Tm = round(median(meltPoint, na.rm = TRUE), 2)) %>% 
                                     dplyr::select_("Cell line/type" = "cell_line_or_type", "Tm"))
  
  output$melt_curve_plot2 <- renderPlot({
    
    ggplot(filter(humanDataTab(), cell_line_or_type == input$selectCell, 
                  toupper(gene_name) %in% toupper(c(input$protein1, input$protein2,
                                    input$protein3))), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = gene_name), size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "Fraction non-denatured") +
      stat_summary(aes(color = gene_name), geom="line", fun.y="median") +
      scale_color_discrete("Gene symbol") +
      theme_bw() + 
      theme(text = element_text(size = 20))
    
  })
  
  output$table3 <- renderTable(filter(humanDataTab(), cell_line_or_type == input$selectCell,
                                      toupper(gene_name) %in% toupper(c(input$protein1, input$protein2,
                                                                        input$protein3))) %>%
                                     group_by(gene_name) %>%
                                     summarize(Tm = round(median(meltPoint, na.rm = TRUE), 2)) %>%
                                     dplyr::select_("Gene symbol" = "gene_name", "Tm"))
  
  output$tm_scatter <- renderPlot({
    
    ggplot(humanCellsTmSpread(), aes_string(input$selectCellTm1, input$selectCellTm2)) + 
      geom_point() +
      coord_fixed() +
      labs(x = bquote(.(input$selectCellTm1) ~ Tm ~ degree*C),
           y = bquote(.(input$selectCellTm2) ~ Tm ~ degree*C)) +
      theme_bw()  + 
      theme(text = element_text(size = 20))
    
  })
  
  output$species_melt_curve_plot1 <- renderPlot({
    
    ggplot(filter(speciesDataTab(), run_name %in% c(input$selectSpecies1, input$selectSpecies2, 
                                           input$selectSpecies3), 
                  toupper(gene_name) %in% toupper(c(input$gene_name1, input$gene_name2,
                                 input$gene_name3))), 
           aes(temperature, as.numeric(fold_change))) +
      geom_point(aes(color = run_name, shape = gene_name), 
                 size = 3, alpha = 0.75) +
      labs(x = expression('Temperature ('*~degree*C*')'), y = "Fraction non-denatured") +
      stat_summary(aes(color = run_name, linetype = gene_name), geom="line", fun.y="median") +
      scale_color_discrete("Species") +
      scale_linetype_discrete("Gene symbol") +
      scale_shape_discrete("Gene symbol") +
      theme_bw() + 
      theme(text = element_text(size = 20),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vertical")
    
  })
  
  output$table <- renderTable(filter(speciesDataTab(), run_name %in% c(input$selectSpecies1, input$selectSpecies2, 
                                                                       input$selectSpecies3), 
                                     toupper(gene_name) %in% toupper(c(input$gene_name1, input$gene_name2,
                                                                       input$gene_name3))) %>% 
                                 group_by(run_name, gene_name) %>% 
                                 summarize(Tm = round(median(meltPoint, na.rm = TRUE), 2)) %>% 
                                 dplyr::select_("Species" = "run_name", "Gene symbol" = "gene_name", "Tm"))
  
  
  # Reactive value for selected dataset ----
  datasetInput <- reactive({
    switch(input$dataset,
           "cross-species" = speciesDataTab(),
           "human" = humanDataTab())
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

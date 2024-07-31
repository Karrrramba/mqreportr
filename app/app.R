#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(bslib)
library(shiny)
library(sortable)
library(shinyFiles)





# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = bs_theme(version = 5, 
                     bootswatch = "flatly", 
                     primary = "#96ded1",
                     base_font = font_google("Fira Sans")),
    navbarPage(
        "MaxQuant ReportR",
        tabPanel("Parameters",
                 sidebarPanel(
                     shinyDirButton('dir', 'Select MaxQuant txt output folder', 'Please select a folder', FALSE)
                 ),
                 mainPanel(
                     fluidRow(
                         column(
                             tags$b("Input selection"),
                             width = 12,
                             bucket_list(
                                 header = "Select tables for analysis",
                                 group_name = "bucket_list_group",
                                 orientation = "horizontal",
                                 add_rank_list(
                                     text = "Drag from here",
                                     labels = list(
                                         "one",
                                         "two",
                                         "three",
                                         htmltools::tags$div(
                                             htmltools::em("Complex"), " html tag without a name"
                                         ),
                                         "five" = htmltools::tags$div(
                                             htmltools::em("Complex"), " html tag with name: 'five'"
                                         )
                                     ),
                                     input_id = "rank_list_1"
                                 ),
                                 add_rank_list(
                                     text = "to here",
                                     labels = NULL,
                                     input_id = "rank_list_2"
                                 )
                             )
                         )
                     ),
                     fluidRow(
                         column(
                             width = 12,
                             tags$b("Result"),
                             column(
                                 width = 12,
                                 
                                 tags$p("input$rank_list_1"),
                                 verbatimTextOutput("results_1"),
                                 
                                 tags$p("input$rank_list_2"),
                                 verbatimTextOutput("results_2"),
                                 
                                 tags$p("input$bucket_list_group"),
                                 verbatimTextOutput("results_3")
                             )
                         )
                     )
                 )
        ),
        tabPanel("Summary stats"),
        tabPanel("Differential expression analysis")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output){
    volumes = getVolumes()() # this makes the directory at the base of your computer.
    observe({
        shinyDirChoose(input, 'dir', roots = volumes, filetypes = c('txt'))
        print(input$folder)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



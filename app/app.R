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
library(shinyFiles)

# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = bs_theme(version = 4, base_font = font_google("Fira Sans")),
    shinyDirButton('folder', 'Select MaxQuant txt output folder', 'Please select a folder', FALSE)
)

# Define server logic required to draw a histogram
server <- function(input, output){
    volumes = getVolumes()() # this makes the directory at the base of your computer.
    observe({
        shinyDirChoose(input, 'folder', roots = volumes, filetypes = c('', 'txt'))
        print(input$folder)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



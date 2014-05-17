#Run this example using runApp("../examples/egshiny/")

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Interactive Unimodal Demo"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("pi0",
                  "pi0:",
                  min = 0,
                  max = 1,
                  value = 0),
      
      helpText("The histogram shows a distribution of z scores, some null (blue) and some non-null (cyan).", 
               "You control what proportion are colored as being from the null by changing pi0.",
    "The remainder are implicitly colored as being non-null.",
    "As you increase pi0, the distribution of the non-null z becomes less unimodal.",
"ash tries to find the largest pi0 that leaves the cyan distribution unimodal.",
"qvalue tries to find the largest pi0 that leaves the cyan distribution non-negative.")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
library(shiny)

#plot a histogram of z scores, highlighting the alternative distribution
#of z scores that is implied by pi0
nullalthist_z = function(z,pi0,nullcol="blue",altcol="cyan",...){
  h=hist(z, freq=FALSE,col=nullcol,nclass=40,ylim=c(-0.5,0.5),...)
  nb = length(h$breaks)
  nullexp = pi0 * (pnorm(h$breaks)[-1] - pnorm(h$breaks[-nb]))/(h$breaks[-1]-h$breaks[-nb])
  h$density = h$density - nullexp
  plot(h,add=TRUE,col=altcol,freq=FALSE)
}
set.seed(100)
z = c(rnorm(5000,0,2),rnorm(5000,0,1))

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$distPlot <- renderPlot({
    # draw the histogram with the specified number of bins
    nullalthist_z(z, input$pi0)
  })
})

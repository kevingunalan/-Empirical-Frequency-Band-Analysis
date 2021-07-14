library(shiny)
library(momentchi2)
library(fields)
library(viridis)
library(signal)
library(fossil)


ui <- fluidPage(
  
  sliderInput(inputId = "Time", label = "Choose total length of time series", 
              value = 50000, min = 1000, max = 100000),
  
  fileInput("Source", label = "Upload the source file", accept = ".R"),
  
  sliderInput(inputId = "Num", label = "Choose number of observations per approximately stationary block", 
              value = 500, min = 100, max = 1000),
  
  sliderInput(inputId = "Tapers", label = "Choose number of tapers to use in multitaper spectral estimator", 
              value = 15, min = 0, max = 50),
  
  radioButtons(inputId = "Signi", label = "Choose significance level to use for testing partition points using FRESH statistic", 
              c("0.1" = 0.1, "0.5" = 0.5, "0.01" = 0.01)),
  
  radioButtons(inputId = "TF", label = "Should the variance of each stationary block be set to one across all blocks?", 
               c("True" = TRUE, "False" = FALSE)),
  plotOutput("Image_Plot")
)

server <- function(input,output) {
  output$Image_Plot <- renderPlot({
    
    source(input$Source)
    X <- eba.simdata(T=input$Time);
    ebaout.wn <- eba.search(X=X$wn,N=input$Num,K=input$Tapers,std=input$TF,alpha=input$Signi);
    ebaout.bL <- eba.search(X=X$bL,N=input$Num,K=input$Tapers,std=input$TF,alpha=input$Signi);
    ebaout.bS <- eba.search(X=X$bS,N=input$Num,K=input$Tapers,std=input$TF,alpha=input$Signi);
    
    image.plot(x=ebaout.wn$mtspec$t,y=ebaout.wn$mtspec$f,z=t(ebaout.wn$mtspec$mtspec), 
               axes = TRUE, col = inferno(256),zlim=c(0,10), 
               xlab='Time',ylab='Hz',xaxs="i");
  })
}

shinyApp(ui = ui, server = server)
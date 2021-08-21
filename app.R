library(shiny)
library(momentchi2)
library(fields)
library(viridis)
library(signal)
library(fossil)

source("EBA_functions.R")

ui <- fluidPage(
  tags$h1("Empirical Frequency Band Analysis"),
  tags$a(tags$strong(tags$em("Source: Empirical Frequency Band Analysis of Nonstationary Time Series")), 
         href = "https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1671199"),
  tags$p(tags$em("Research reported in this publication was supported by the National Institute of General Medical Sciences 
         of the National Institutes of Health under Award Number R01GM140476. The content is solely the responsibility
         of the authors and does not necessarily represent the official views of the National Institutes of Health")),
  tags$p(tags$strong("Authors: Kevin Gunalan Antony Michael Raj, Tejaswini Prakash")), (tags$em(tags$u("Under the guidance of Professor Scott Bruce"))),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      selectInput("Simsetting", "Simulation Setting",
                  c("White Noise" = "W",
                    "Linear" = "L",
                    "Sinusoidal" = "S")),
      sliderInput(inputId = "Time", label = "Choose total length of time series", 
                  value = 50000, min = 1000, max = 50000, step = 1000),
      
      sliderInput(inputId = "Num", label = "Choose number of observations per approximately stationary block", 
                  value = 500, min = 100, max = 1000, step = 100),
      
      sliderInput(inputId = "Tapers", label = "Choose number of tapers to use in multitaper spectral estimator", 
                  value = 15, min = 0, max = 50, step = 5),
      
      radioButtons(inputId = "Signi", label = "Choose significance level", 
                   c("0.05" = 0.05, "0.01" = 0.01), selected = "0.05"),
      
      radioButtons(inputId = "TF", label = "Standardize", 
                   c("True" = TRUE, "False" = FALSE), selected = FALSE),
      actionButton("go", label = "Run")
    ),
    mainPanel(
      plotOutput("Image_Plot")
    )
  ))

server <- function(input,output) {
  plot.list <- eventReactive(input$go, ignoreNULL = FALSE, {
    X = eba.simdata(T=input$Time)
    
  
  if (input$Simsetting == "W"){
    ebaout.wn <- eba.search(X=X$wn,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    plot.x = ebaout.wn$mtspec$t
    plot.y = ebaout.wn$mtspec$f
    plot.z = t(ebaout.wn$mtspec$mtspec)
    plot.main = "Multitaper Spectrogram for White Noise Setting"
    plot.h = ebaout.wn$part.final[c(-1,-length(ebaout.wn$part.final))]
    plot.data = X$wn  
    
    
    
  } else if (input$Simsetting == "L") {
    ebaout.bL <- eba.search(X=X$bL,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    plot.x = ebaout.bL$mtspec$t
    plot.y = ebaout.bL$mtspec$f
    plot.z = t(ebaout.bL$mtspec$mtspec)
    plot.main = "Multitaper Spectrogram for Linear Setting"
    plot.h = ebaout.bL$part.final[c(-1,-length(ebaout.bL$part.final))]
    plot.data = X$bL
    
    
  } else if (input$Simsetting == "S") {
    ebaout.bS <- eba.search(X=X$bS,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    plot.x = ebaout.bS$mtspec$t
    plot.y = ebaout.bS$mtspec$f
    plot.z = t(ebaout.bS$mtspec$mtspec)
    plot.main = "Multitaper Spectrogram for Sinusoidal Setting"
    plot.h = ebaout.bS$part.final[c(-1,-length(ebaout.bS$part.final))]
    plot.data = X$bS
    
  }
    
  list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
       plot.main = plot.main, plot.h = plot.h, plot.data = plot.data)
    
  });
  

  output$Image_Plot <- renderPlot({
    par(mfrow=c(2,1))
    plot(x = seq(0,1,length.out = length(plot.list()[[6]])), y = plot.list()[[6]], type = "l", xlab = "Time", ylab = "")
    image.plot(x=plot.list()[[1]], y=plot.list()[[2]], z=plot.list()[[3]], 
               axes = TRUE, col = inferno(256), 
               xlab='Time',ylab='Hz',xaxs="i", main = plot.list()[[4]]); 
    abline(h=plot.list()[[5]], col = "green");
    
    
        
  })
}

shinyApp(ui = ui, server = server)
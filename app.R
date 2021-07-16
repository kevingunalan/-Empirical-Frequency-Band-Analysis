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
  tags$p(tags$em("Research reported in thos publication was supported by the National Insitute Of General Medical Sciences 
         of the National Institutes of Health under Award Number R01GM140476. The content is solely the responsbility
         of the authors and does not necessarily represent the official views of the National Institutes of Health")),
  tags$p(tags$strong("Authors: Tejaswini Prakash, Kevin Gunalan Antony Michael Raj")), (tags$em(tags$u("Under the guidance of Professor Scott Bruce"))),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "Time", label = "Choose total length of time series", 
              value = 50000, min = 1000, max = 100000),
  
      sliderInput(inputId = "Num", label = "Choose number of observations per approximately stationary block", 
              value = 500, min = 100, max = 1000),
  
      sliderInput(inputId = "Tapers", label = "Choose number of tapers to use in multitaper spectral estimator", 
              value = 15, min = 0, max = 50),
  
      radioButtons(inputId = "Signi", label = "Choose significance level to use for testing partition points using FRESH statistic", 
              c("0.5" = 0.5, "0.1" = 0.1, "0.05" = 0.05, "0.01" = 0.01), selected = "0.05"),
  
      radioButtons(inputId = "TF", label = "Should the variance of each stationary block be set to one across all blocks?", 
             c("True" = TRUE, "False" = FALSE), selected = FALSE)
    ),
    mainPanel(
    plotOutput("Image_Plot")
  )
))

server <- function(input,output) {
  output$Image_Plot <- renderPlot({
  
    X <- eba.simdata(T=input$Time);
    
    ebaout.wn <- eba.search(X=X$wn,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    ebaout.bL <- eba.search(X=X$bL,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    ebaout.bS <- eba.search(X=X$bS,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
    
    par(mfrow=c(3,1),oma=c(0,0,2,0),mar=c(4,4,1.5,2)+0.1);
    
    image.plot(x=ebaout.wn$mtspec$t,y=ebaout.wn$mtspec$f,z=t(ebaout.wn$mtspec$mtspec), 
                   axes = TRUE, col = inferno(256),zlim=c(0,10), 
                   xlab='Time',ylab='Hz',xaxs="i"); 
    abline(h=ebaout.wn$part.final[c(-1,-length(ebaout.wn$part.final))],col='green');
    
    image.plot(x=ebaout.bL$mtspec$t,y=ebaout.bL$mtspec$f,z=t(ebaout.bL$mtspec$mtspec), 
               axes = TRUE, col = inferno(256), zlim=c(0,10), 
               xlab='Time',ylab='Hz',xaxs="i");
    abline(h=ebaout.bL$part.final[c(-1,-length(ebaout.bL$part.final))],col='green');
    
    image.plot(x=ebaout.bS$mtspec$t,y=ebaout.bS$mtspec$f,z=t(ebaout.bS$mtspec$mtspec), 
               axes = TRUE, col = inferno(256), zlim=c(0,10), 
               xlab='Time',ylab='Hz',xaxs="i");
    abline(h=ebaout.bS$part.final[c(-1,-length(ebaout.bS$part.final))],col='green');
    mtext("Multitaper Spectrogram for White Noise(top) Linear(middle) and Sinusoidal(bottom)", outer = TRUE, cex = 1);
  
    
    })
  
}

shinyApp(ui = ui, server = server)
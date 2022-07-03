library(shiny)
library(momentchi2)
library(fields)
library(viridis)
library(signal)
library(fossil)

wd = "C:/Users/Administrator/Desktop/Summer Project/"
setwd(wd)

source("EBA_functions.R")
df <- read.csv("xtsample.csv")
df <- df[[1]]

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
      conditionalPanel(condition = "input.tabselected==1",
      
      selectInput("Simsetting", "Simulation Setting",
                  c("White Noise" = "W",
                    "Linear" = "L",
                    "Sinusoidal" = "S")),
      sliderInput(inputId = "Time", label = "Choose total length of time series", 
                  value = 50000, min = 1000, max = 50000, step = 1000),
      
      sliderInput(inputId = "Num", label = "Choose number of observations per approximately stationary block", 
                  value = 500, min = 50, max = 1000, step = 100),
      
      sliderInput(inputId = "Tapers", label = "Choose number of tapers to use in multitaper spectral estimator", 
                  value = 15, min = 0, max = 50, step = 5),
      
      radioButtons(inputId = "Signi", label = "Choose significance level", 
                   c("0.05" = 0.05, "0.01" = 0.01), selected = "0.05"),
      
      radioButtons(inputId = "TF", label = "Standardize", 
                   c("True" = TRUE, "False" = FALSE), selected = FALSE),
      actionButton("go", label = "Run")),
      
      conditionalPanel(condition = "input.tabselected==2",
      fileInput("file_csv", "Choose CSV File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      checkboxInput("header", "Header", TRUE),
      sliderInput(inputId = "Num2", label = "Choose number of observations per approximately stationary block", 
                  value = max(floor(length(df)/6),50), min = 50, max = length(df), step = 1),
      
      sliderInput(inputId = "Tapers2", label = "Choose number of tapers to use in multitaper spectral estimator", 
                  value = 15, min = 0, max = floor(2*length(df)*.15-1), step = 1),
      
      radioButtons(inputId = "Signi2", label = "Choose significance level", 
                   c("0.05" = 0.05, "0.01" = 0.01), selected = "0.05"),
      
      radioButtons(inputId = "TF2", label = "Standardize", 
                   c("True" = TRUE, "False" = FALSE), selected = FALSE),
      actionButton("go2", label = "Run"))
    ),
    mainPanel(
      tabsetPanel(type = "tabs", id = 'tabselected', selected = 1,
                  tabPanel("Simulation Setting", value = 1),
                  tabPanel("File Upload", value = 2)),
      conditionalPanel(condition = "input.tabselected==1",
      plotOutput("Image_Plot"),
      radioButtons(inputId = "downloadType", label = "Select download type", choices = list("png","pdf")),
      downloadButton('downloadData','Download the plot')
      ),
      conditionalPanel(condition = "input.tabselected==2",
                       plotOutput("Image_Plot2"),
                       ))
  ))

server <- function(input,output) {
  plot.list2 <- eventReactive(input$go2, ignoreNULL = FALSE, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    dataf <- read.csv(file$datapath, header = input$header)
    dataf <- dataf[[1]]
    ebaoutfu <- eba.search(X=dataf,N= as.numeric(input$Num2),K=as.numeric(input$Tapers2),std=input$TF2,alpha=as.numeric(input$Signi2))
    plot.x = ebaoutfu$mtspec$t
    plot.y = ebaoutfu$mtspec$f
    plot.z = t(ebaoutfu$mtspec$mtspec)
    plot.main = "Multitaper Spectrogram"
    plot.h = ebaoutfu$part.final[c(-1,-length(ebaoutfu$part.final))]
    plot.data = dataf
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.h = plot.h, plot.data = plot.data)})
    
    
    output$Image_Plot2 <- renderPlot({
      par(mfrow=c(2,1))
      plot.ts(plot.list2()[[6]])
      image.plot(x=plot.list2()[[1]], y=plot.list2()[[2]], z=plot.list2()[[3]], 
                 axes = TRUE, col = inferno(256), 
                 xlab='Time',ylab='Hz',xaxs="i", main = plot.list2()[[4]]); 
      abline(h=plot.list2()[[5]], col = "green");
      })
    
    
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
    
    output$downloadData <- downloadHandler(
      filename = function(){
        paste("White_noise",input$downloadType,sep = ".") 
      },
      content = function(file){
        if(input$downloadType == "png")
          png(file)
        else
          pdf(file)
        par(mfrow=c(2,1))
        plot(x = seq(0,1,length.out = length(plot.list()[[6]])), y = plot.list()[[6]], type = "l", xlab = "Time", ylab = "")
        image.plot(x=plot.list()[[1]], y=plot.list()[[2]], z=plot.list()[[3]], 
                   axes = TRUE, col = inferno(256), 
                   xlab='Time',ylab='Hz',xaxs="i", main = plot.list()[[4]]); 
        abline(h=plot.list()[[5]], col = "green");
        
        dev.off()
      }
    )
    
    
        
  })
}

shinyApp(ui = ui, server = server)
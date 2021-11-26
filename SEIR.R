##Shiny Dashboard for SEIR Models in app.R
#Name: SIR and SEIR Models for Infectious Disease Epidemiology

#load the required libraries
library(shiny)
library(shinydashboard)
library(shinyBS)
library(magrittr)
library(tidyverse)
library(rmarkdown)
library(plotly)
library(DT)
library(shinyWidgets)
library(shinycssloaders)
library(deSolve)
library(fontawesome)


ui <- dashboardPage(
  dashboardHeader(title="Simulating Diseases Dynamics"),
  
  dashboardSidebar(
    hr(),
    sidebarMenu(id = "menu",
                menuItem("Simulate Disease Dynamics", tabName = "exp-model", icon = icon("chart-line"), selected = TRUE),
                menuItem("About", tabName = "readme", icon = icon("id-card")),
                menuItem("Code",  icon = icon("code"))
    ),
    
   conditionalPanel(condition = 'input.menu == "exp-model"',
                     actionButton("go", "Run Model Simulation"), 
                     splitLayout(cellWidths = c("50%", "50%"),
                     numericInput(inputId = "beta", label = "Infection Rate", value = 0.8, min = 0, max = 50),
                     numericInput(inputId = "gamma", label = "Recovery parameter", value = 0.3, min = 0)),
                    
                     splitLayout(cellWidths = c("50%", "50%"),
                     numericInput(inputId = "sigma", label = "Rate from exposure to infection", value = 1/3, min = 0),
                     numericInput(inputId = "S", label = "Susceptible population (S)", value = 1-1e-6, min = 0)),
                    
                     splitLayout(cellWidths = c("50%", "50%"),
                     numericInput(inputId = "E", label = "Exposed population (E)", value = 0, min = 0),
                     numericInput(inputId ="I", label = "Number infected (I)", value = 1e-5, min = 0)),
                     
                     sliderInput("R", 
                                 "Reproduction number:",
                                 value = 0,
                                 min = 0, 
                                 max = 50,
                                 step = 0.1),
                     numericInput(inputId = "tfinal", label = "Final simulation time (in days):", value = 100, min = 0, max = 730)),
            hr(),
            helpText("Developed by ", 
                              a("Amos Okutse", href = "http://brown.edu"), 
                              "Powered by",
                              a("shinyapps", href = "http://shinyapps.com"),
                              style = "padding-left:1em; padding-right:1em;position:absolute;")
    ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "exp-model",
              fluidPage(
                
                fluidRow(
                  tabBox(
                    width = 12,
                    height = "800px",
                    title = "Simulated Infectious Disease Model Dynamics", 
                    side = "right",
                    tabPanel(title = "Trajectories",
                             plotlyOutput("plot", width = "95%", height = "650px") %>% 
                               withSpinner()
                    )
                    
                  )),
                
                fluidRow(
                  tabBox(
                    width = 12,
                    title = "Simulated Data Table", 
                    side = "right",
                    #tabPanel(title = "Summary",
                    #tableOutput("model_sum_tab") %>% 
                    #withSpinner()),
                    tabPanel(title = "Trajectories",
                             DT::dataTableOutput("simulation_results") %>% 
                               withSpinner()),
                    tabPanel(title = "Download Data", 
                            downloadLink("downloadData", "Download"))
                  )
                )
              )
            ) #,
      #add the code for the readme file here, and update the About menu tab to have details and pictures of the model employed in the simulation
      
      #tabItem(tabName = "readme",
      #        includeMarkdown("README.md")
      #),
      #tabItem(tabName = "ui",
              #box( width = NULL, status = "primary", solidHeader = TRUE, title = "UI",
                   #downloadButton('downloadData2', 'Download'),
                   #br(),br(),
                   #pre(includeText("ui.R"))
                  
    ) 
    
  ),
  
  skin="purple"
)




server <- function(input, output) { 
  
  #create the SEIR function to be used in solving ODEs
  seir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- - beta*S*I
      dE <- beta*S*I - sigma*E
      dI <- sigma*E - gamma*I
      dR <- gamma*I
      return(list(c(dS, dE, dI, dR)))
    })
  }
  
  observeEvent(input$go, {
   
  output$plot <- renderPlotly({
    #use user inputs as parameters for the simulation model
    parameters <- c(beta=input$beta, gamma=input$gamma, sigma=input$sigma)
    
    #initial values for S, E, I, R
    y<-c(S=input$S, E=input$E, I=input$I, R=input$R)
    
    #create a vector of the simulation times
    times=seq(0,input$tfinal,by=1)
    
    #solve the system of differentials using the deSolve package
    results <- ode(y=y, times=times, func=seir, parms=parameters)
    
    #create a data frame of the outputs from the model
    df<-as.data.frame(results)  
    
        plot_ly(data.frame(df), x = ~df[,1], y = ~df[,2], name = 'S', type = 'scatter', mode = 'lines', color=I("blue")) %>%
        add_trace(y = ~df[,3], name = 'E', mode = 'lines', color=I("orange")) %>%
        add_trace(y = ~df[,4], name = 'I', mode = 'lines', color=I("red")) %>%
        add_trace(y = ~df[,5], name = 'R', mode = 'lines', color=I("purple")) %>%
        layout(xaxis = list(title = "Simulation Period in Days"), yaxis=list(title = "Total Susceptible Population"))
    }) 
  })
  
  #Simulated Data Table
  observeEvent(input$go, {
    output$simulation_results<-DT::renderDataTable({
      parameters <- c(beta=input$beta, gamma=input$gamma, sigma=input$sigma)
      y<-c(S=input$S, E=input$E, I=input$I, R=input$R)
      times=seq(0,input$tfinal,by=1)
      results <- DT::datatable(ode(y=y, times=times, func=seir, parms=parameters))
      
    })
  })
  
  #prepare the simulated data for user download
  observeEvent(input$go, {
    datasetinput<-reactive({
      parameters <- c(beta=input$beta, gamma=input$gamma, sigma=input$sigma)
      y<-c(S=input$S, E=input$E, I=input$I, R=input$R)
      times=seq(0,input$tfinal,by=1)
      results <- data.frame(ode(y=y, times=times, func=seir, parms=parameters)) 
    })
    
    output$downloadData<-downloadHandler(filename = function(){paste("simulated_data-", Sys.Date(), ".csv", sep = "")},
                                         content = function(file){write.csv(datasetinput(), file, row.names=FALSE)})
      
  })
  
  
  #Summarize the results of the simulation to include min, max, mean, median using dplyr
  
  
  
  }


shinyApp(ui, server)


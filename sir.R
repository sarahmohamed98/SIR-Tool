library(reshape)
library(ggplot2)
library(plotly)
library(shiny)
library(shinydashboard)
library(plotly)
library(mathjaxr)
library(igraph)
library(epimdr)
body <- dashboardBody(
  fluidRow(
    tabBox(
      side = "right", height = "1000px", width = 15,
      selected = "Plot",
      tabPanel("About Tool", h4("
               This tool  solves and plots a deterministic,
               compartmental epidemic models (DCMs),
               The model simulations are driven by Rstudio.
               
               
               Models here are limited to basic one-group homo
               genous mixing models with a limited set of parameters,
               initial conditions, and control settings.
               
               This web application, built with Shiny may and shinydashboard")),
      tabPanel("Data", "Tabular Data",dataTableOutput("dataT")),
      tabPanel("Model Equations", "Equations",uiOutput("entered")),
      tabPanel("Plot","",
               
               tabBox( side = "left", width = 1000,
                       selected = "Interactive SIR Model",
                       tabPanel("Interactive SIR Model",
                                plotlyOutput("graph3"),
                       ),
                       tabPanel("Seperated SIR Model",
                                plotlyOutput("graph2"),
                       ),
                       tabPanel("Pie Chart",
                                plotlyOutput("graph4")
                       ),
                       tabPanel("Intervention",
                                
                                fluidRow(
                                  box(plotlyOutput("graph5")),
                                  box(plotlyOutput("graph6")),
                                ),
                                
                       ),
                       
                       tabPanel("SEIR Model",
                                plotlyOutput("graph7")
                       ),
                       tabPanel("R0",
                                plotOutput ("graph8")
                       ),
                       tabPanel("Parameter Estimation",
                                plotlyOutput ("graph10")
                       ),
                       
                       tabPanel("Testing",
                                plotlyOutput ("graph9")
                       ),
               ),
               
               hr(),
               
               div(
                 style = "display: inline-block",
                 sliderInput("infectivity", "Infectivity Rate", min = 0, max = 1, value = 0.1, step = 0.001, width = 150)
               ),
               div(style = "display: inline-block",
                   sliderInput("contactRate", "Contact Rate", min = 0, max = 5, value = 0.1, step = 0.001, width = 150)
               ),
               div(
                 style = "display: inline-block",
                 sliderInput("gamma", "Recovery Rate", min = 0, max = 1, value = 0.1, step = 0.001, width = 150)
               ),
               
               
               
               div(
                 style = "display: inline-block",
                 sliderInput("arrivalRate", "Arrival Rate (Birth,Immigration)", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               
               div(
                 style = "display: inline-block",
                 sliderInput("d1", "Death from S", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               div(
                 style = "display: inline-block",
                 sliderInput("d2", "Death from I", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               div(
                 style = "display: inline-block",
                 sliderInput("d3", "Death from R", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               div(
                 style = "display: inline-block",
                 sliderInput("d4", "Death from E", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               div(
                 style = "display: inline-block",
                 sliderInput("p", "Vaccine", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               
               div(
                 style = "display: inline-block",
                 sliderInput("q", "quarantine rate", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),
               
               div(
                 style = "display: inline-block",
                 sliderInput("delta", "delta", min = 0, max = 1, value = 0, step = 0.001, width = 150)
               ),      
      )
    )
  )
)

shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "SIR Model"),
    dashboardSidebar(
      #selectInput("select", h5("Model type"),
      #           choices = list("SIR Model" = 1, "SEIR Model" = 2, selected = 1)
      #),
      
      numericInput("totalPopulation",label=h5("Total Population N"),value=1,step=.1 ,min=0 , max=1),
      numericInput("initialS",label=h5("Initial Number of Susceptible"), value=0.9,step=.1 ,min=0 , max=1),
      numericInput("initialI",label=h5("Initial Number of Infected"), value=0.1,step=.1 ,min=0 , max=1),
      numericInput("initialR",label=h5("Initial Number of Recoverd"), value=0,step=.1, min=0 , max=1),
      numericInput("quarantine",label=h5("Quarantine"),value=.1,step=.1 ,min=0 , max=1),
      
      numericInput("timeDuration",label=h5("Time Duration"), value=100, min=0),
      numericInput("h",label=h5("Euler h"), value=0.01, max=1),
      
      #  selectInput("select", h5("Time Duration By"),
      #             choices = list("Days" = 1, "Months" = 2, "Years" = 3 ,selected = 1)),
      
      #  numericInput("stepSize",label=h5("Size Step"), value=0, min=0),
      hr(),
      numericInput("E0",label=h5("Initial Number of Exposed (SEIR)"), value=.1,step=.1 ,min=0 , max=1),
      numericInput("e",label=h5("Latent Period (SEIR)"), value=.1,step=.1 ,min=0 , max=1)
    ),
    
    
    
    body
  ),
  server = function(input, output) {
    
    
    ################ All in one SIR ########################
    SIR = function(t,h,N0,S0,I0,R0,infectivity,contactRate,gamma,b,d1,d2,d3,p,q,delta){
      lambda = infectivity * contactRate
      times = seq(0,t,h)
      
      S=rep(0,length(times))
      I=rep(0,length(times))
      R=rep(0,length(times))
      Q=rep(0,length(times))
      N=rep(0,length(times))
      
      S[1]=S0
      I[1]=I0
      R[1]=R0
      N[1]=N0
      Q[0]=0
      
      for(j in 2:length(times)){
        
        S[j] <- S[j-1]+(b*N[j-1]-(lambda*S[j-1]*I[j-1] / N[j-1])- d1 * S[j-1]- p * S[j-1])*h
        I[j] <- I[j-1]+((lambda*S[j-1]*I[j-1] / N[j-1]) - (gamma+d2+q)* I[j-1])*h
        R[j] <- R[j-1]+(gamma * I[j-1]-d3*R[j-1]+p * S[j-1]+delta*Q[j-1])*h
        Q[j] <- Q[j-1]+ (q*I[j-1]-Q[j-1]*delta)*h
        N[j] <- Q[j-1]+S[j-1]+I[j-1]+R[j-1]+(b*N[j-1])*h
        
        S[j] = min(N[j],S[j])
        I[j] = min(N[j],I[j])
        R[j] = min(N[j],R[j])
        Q[j] = min(N[j],Q[j])
        
        S[j] = max(0,S[j])
        I[j] = max(0,I[j])
        R[j] = max(0,R[j])
        Q[j] = max(0,Q[j])
        
        
      }
      out <-as.data.frame(cbind(times,N,S,I,R,Q))
      
      out = out[(0:t)/h+1,]
      return (out[2:nrow(out),])
      
    }
    
    
    
    #_________________________________________________________________________SIR1 graph4 pie chart
    
    #_________________________________________________________________vaccine
    
    
    #_____________________________________________________________________Quarntine
    
    #_____________________________________________________________________________________SEIR
    
    SEIR = function(t,N0,S0,I0,R0,infectivity,contactRate,gamma,b,d1,d2,d3,d4,e,E0){
      
      lambda = infectivity * contactRate
      times=seq(0,t,1)
      
      S=rep(0,length(times))
      E=rep(0,length(times))
      I=rep(0,length(times))
      R=rep(0,length(times))
      N=rep(0,length(times))
      
      
      S[1]=S0
      E[1]=E0
      I[1]=I0
      R[1]=R0
      N[1]=N0
      
      
      
      for(j in 2:length(times)){
        
        
        
        
        
        
        S[j] <- S[j-1]+b*N[j-1]-(lambda*S[j-1]*I[j-1] / N)- d1 * S[j-1]
        E[j] <- E[j-1]+(lambda*S[j-1]*I[j-1] / N) - (e+d2)* E[j-1]
        I[j] <- I[j-1]+ e*E[j-1] -(gamma+d3)* I[j-1]
        R[j] <- R[j-1]+gamma * I[j-1] - d4*R[j-1]
        N[j] <- N[j-1]+b*N[j-1]-d1*S[j-1]-d2*E[j-1]-d3*I[j-1]-d4*R[j-1]
        
        
        
        if( S[j]<0)
          
        {  S[j]=  0
        }
        
        
        if( S[j]==Inf)
          
        {S[j]=  3.428008e+162
        }
        
        if( E[j]<0)
        {
          E[j]=  0
        }
        
        
        if( E[j]==Inf)
        {
          E[j]=  3.428008e+162
        }
        
        
        if( I[j]<0)
        {
          I[j]=  0
        }
        
        
        if( I[j]==Inf)
        {
          I[j]=  3.428008e+162
        }
        
        
        
        if( R[j]<0)
        {
          R[j]=  0
        }
        
        
        if( R[j]==Inf)
        {
          R[j]=  3.428008e+162
        }
        
        
        if( N[j]<0)
        {
          N[j]=  0
        }
        
        
        if( N[j]==Inf)
        {
          N[j]=  3.428008e+162
        }
        
      }
      out <-as.data.frame(cbind(times,N,S,E,I,R))
      return (out)
    }
    
    #_____________________________________________________________________________________R0
    
    SIRr0 = function(infectivity,contactRate,gamma) {
      lambda = infectivity * contactRate
      r0 =  lambda /  gamma
      return(r0)
      
    }
    SIRtest = function (t,h,S0,I0,R0,infectivity,contactRate,gamma)
    {
      lambda = infectivity * contactRate
      times = seq(0,t,h)
      
      S = rep(0,length(times))
      I = rep(0,length(times))
      R = rep(0,length(times))
      N = rep(0,length(times))
      
      S[1]=S0
      I[1]=I0
      R[1]=R0
      N[1]=S0+I0+R0
      
      
      for(j in 2:length(times)){
        
        S[j] <- S[j-1]-((lambda*S[j-1]*I[j-1])/N[j-1])*h
        I[j] <- I[j-1]+((lambda*S[j-1]*I[j-1]/N[j-1])-gamma*I[j-1])*h                            
        R[j] <- R[j-1]+(gamma*I[j-1])*h
        N[j] <- S[j]+I[j]+R[j]
        
        
        
        S[j]=max(S[j],0)
        I[j]=max(I[j],0)
        R[j]=max(R[j],0)
        N[j]=max(N[j],0)
        
        
        S[j]=min(S[j],N[j])
        I[j]=min(I[j],N[j])
        R[j]=min(R[j],N[j])
        N[j]=min(N[j],N[j])
        
      }
      
      
      out <-as.data.frame(cbind(times,N,S,I,R))
      
      
      
      out = out[(0:t)/h+1,]
      return (out[2:nrow(out),])
      
      
    }
    
    #_____________________________________________________________________________________graphs
    
    
    
    output$graph2 <- renderPlotly({
      mydata= SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      SIR_data_full <- melt(as.data.frame(mydata), id="times")
      graph2 <- SIR_data_full %>%
        ggplot() + geom_line(aes(x=times, y=value, color=variable))+
        facet_wrap( ~ variable, scales="free")
    })
    
    output$graph3 <- renderPlotly({
      mydata= SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      SIR_data_full <- melt(as.data.frame(mydata), id="times")
      graph3 <- SIR_data_full %>%
        ggplot(aes(x=times , y=value , color=variable))+
        geom_point(size=0.5)+
        geom_line()+
        theme_bw()+
        labs(title= "SIR model")+
        xlab("Time by days")+
        ylab("propotion of S , I , R") +
        scale_y_continuous(labels=scales::comma)+
        ##geom_smooth(method="lm" )-> g
        geom_smooth(method="loess", span=0.25) -> g
      
      ggplotly(g)
      
      
    })
    
    #SIR1 = function(t,N0,S0,I0,R0,lambda,gamma,b,d1,d2,d3){
    output$graph4 <- renderPlotly({
      
      out = SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      mydata1 = data.frame()
      
      for(i in 1:nrow(out)){
        
        mydata1 = rbind(mydata1,rbind( cbind(out$times[i],"S",out$S[i]),cbind(out$times[i],"I",out$I[i]),
                                       cbind(out$times[i],"R",out$R[i]),cbind(out$times[i],"Q",out$Q[i])))
      }
      
      names(mydata1)=c("times","compartment","cases")
      
      
      
      fig <- plot_ly(mydata1, labels = ~compartment, frame=~times,values = ~cases,marker = list(colors = c('blue', 'red','green','orange')), type = 'pie')
      
      
      fig
      
      
    })
    #    SIRV = function(t,N0,S0,I0,R0,lambda,gamma,b,d1,d2,d3,p){
    
    output$graph5 <- renderPlotly({
      mydata= SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      SIR_data_full <- melt(as.data.frame(mydata), id="times")
      graph5 <- SIR_data_full %>%
        ggplot(aes(x=times , y=value , color=variable))+
        geom_point(size=0.5)+
        geom_line()+
        theme_bw()+
        labs(title= "SIR model (Vaccine)")+
        xlab("Time by days")+
        ylab("propotion of S , I , R") +
        scale_y_continuous(labels=scales::comma)+
        ##geom_smooth(method="lm" )-> g
        geom_smooth(method="loess", span=0.25) -> g
      
      ggplotly(g)
      
      
    })
    
    #SIRQ = function(t,S0,I0,R0,Q0,lambda,gamma,b,d1,d2,d3,q,delta){    Quarantine
    
    output$graph6 <- renderPlotly({
      
      mydata= SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate*input$h,gamma=input$gamma*input$h,b=input$arrivalRate*input$h,d1=input$d1*input$h,d2=input$d2*input$h,d3=input$d3*input$h,p=input$p*input$h,q=input$q*input$h,delta=input$delta*input$h)
      
      SIR_data_full <- melt(as.data.frame(mydata), id="times")
      graph6 <- SIR_data_full %>%
        ggplot(aes(x=times , y=value , color=variable))+
        geom_point(size=0.5)+
        geom_line()+
        theme_bw()+
        labs(title= "SIR model (Quarantine)")+
        xlab("Time by days")+
        ylab("propotion of S , I , R") +
        scale_y_continuous(labels=scales::comma)+
        ##geom_smooth(method="lm" )-> g
        geom_smooth(method="loess", span=0.25) -> g
      
      ggplotly(g)
      
      
    })
    #_________________________________________________________graph SEIR
    #SIRE = function(t,N0,S0,I0,R0,lambda,gamma,b,d1,d2,d3,d4,e,E0){
    output$graph7 <- renderPlotly({
      mydata = SEIR(input$timeDuration,input$totalPopulation,input$initialS,input$initialI,input$initialR,
                    input$infectivity,input$contactRate,input$gamma,input$arrivalRate,input$d1,input$d2,input$d3,input$d4,input$e,input$E0)
      
      SIR_data_full <- melt(as.data.frame(mydata), id="times")
      graph7 <- SIR_data_full %>%
        ggplot(aes(x=times , y=value , color=variable))+
        geom_point(size=0.5)+
        geom_line()+
        theme_bw()+
        labs(title= "SEIR model")+
        xlab("Time by days")+
        ylab("propotion of S , E , I , R") +
        scale_y_continuous(labels=scales::comma)+
        ##geom_smooth(method="lm" )-> g
        geom_smooth(method="loess", span=0.25) -> g
      
      ggplotly(g)
      
      
    })
    
    output$graph8 <- renderPlot({
      
      
      
      r = SIRr0(input$infectivity,input$contactRate,input$gamma)
      
      if((r-floor(r))>=0.5)
        r0 = floor(0.5+r)
      else
        r0 = floor(r)
      if (r0>=1){
        g = make_tree((r0+1)+r0*r0, children=r0,mode="undirected")
        
        
      }else {
        g = make_tree(1, children=1,mode="undirected")
      }
      plot(g)
    })
    
    output$graph9 <- renderPlotly({
      
      
      
      data(flu)
      dataset = as.data.frame(flu)
      #View(dataset)
      
      dataset = data.frame(dataset)
      
      
      
      
      
      
      mydata = SIR(t=nrow(dataset),h=1,N0=764,S0=761,I0=3,R0=0,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      
      fig <- plot_ly()
      fig <- fig %>%
        add_markers(data=dataset, name="Actual curve", x = ~day, y = ~cases,type = "scatter", mode = "markers+lines")
      fig <- fig %>%
        add_trace(data=mydata, name=" Predictive curve", x = ~times , y =~I,type = "scatter", mode = "markers+lines")
      # show figure
      fig
    })
    
    output$graph10 <- renderPlotly({
      dataset = read.csv("/cloud/project/owidmonkeypox.csv")
      #View(dataset)
      
      dataset = data.frame(dataset)
      dataset = dataset[dataset$location=="United States",]
      #View(dataset)
      
      dataset$S=dataset$total_cases
      
      for(i in 22: nrow(dataset)){
        dataset[i,"S"] = dataset[i,"total_cases"]-dataset[(i-21),"total_cases"]
        
      }
      dataset$day = 1:nrow(dataset)
      N= 350*1000000
      
      
      sirmod = function(t, y, params) {
        S = y[1]
        I = y[2]
        R = y[3]
        with(as.list(params), {
          dS = -beta * S * I/N
          dI = beta * S * I/N - gamma * I
          dR = gamma * I
          res = c(dS, dI, dR)
          list(res)
        })
      }
      
      lfn2 = function(p, I, N) {
        times = seq(1, 245, by = 1)
        start = c(S = N, I = 19, R = 0)
        paras = c(beta = p[1], gamma = p[2], N = N)
        out = as.data.frame(ode(start, times = times,
                                sirmod, paras))
        n = length(I)
        rss = sum((I - out$I)^2)
        
        return(log(rss) * (n/2) - n * (log(n) -
                                         log(2 * pi) - 1)/2)
      }
      
      
      
      paras0 = c(1.2, 1/2)
      flufit = optim(paras0, lfn2, I = dataset$S, N = N,
                     hessian = TRUE)
      
      
      times = seq(1, nrow(dataset), by=1)
      start = c(S=350*1000000, I=19, R = 0)
      paras=c(beta=flufit$par[1], gamma=flufit$par[2], N=350*1000000)
      out = as.data.frame(ode(start, times=times,
                              sirmod, paras))
      #plot(out$time, out$I, ylab="Prevalence",
      #    xlab="Day", type="l")
      #points(dataset$day, dataset$S)
      
      
      
      mydata = SIR(t=nrow(dataset),h=1,N0=350*1000000,S0=350*1000000-9,I0=9,R0=0,infectivity=input$infectivity,contactRate=input$contactRate,gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
      
      
      
      fig <- plot_ly()
      fig <- fig %>%
        add_markers(data=dataset, name="Actual curve", x = ~day, y = ~S,type = "scatter", mode = "markers+lines")
      fig <- fig %>%
        add_trace(data=mydata, name="Predictive curve", x = ~times , y =~I,type = "scatter", mode = "markers+lines")
      # show figure
      fig
      
    })
    
    # The currently selected tab from the first box
    output$tabset1Selected <- renderText({
      input$tabset1
    })
    
    output$dataT = renderDataTable({mydata = SIR(t=input$timeDuration,h=input$h,N0=input$totalPopulation,
                                                 S0=input$initialS,I0=input$initialI,R0=input$initialR,infectivity=input$infectivity,contactRate=input$contactRate,
                                                 gamma=input$gamma,b=input$arrivalRate,d1=input$d1,d2=input$d2,d3=input$d3,p=input$p,q=input$q,delta=input$delta)
    
    })
    
    output$entered = renderUI({
      withMathJax(
        helpText("Population $$\\frac{dN}{dt} =
\\ bN - d1S - d2I - d3R$$"),
        helpText("Susceptible $$\\frac{dS}{dt} =
\\ - ??IS$$"),
        helpText("Infecitous $$\\frac{dI}{dt} =
\\ ??IS  -   ??R$$"),
        helpText("Removed $$\\frac{dR}{dt} =
\\  ?I(t) ? d4R(t)$$"),
        hr(),
        
        helpText("Population $$\\frac{dN}{dt} =
\\ S(t) + E(t) + I(t) + R(t)$$"),
        helpText("Susceptible $$\\frac{dS}{dt} =
\\ bN(t) - ??I(t)S(t) - d1S(t)$$"),
        helpText("Exposed $$\\frac{dE}{dt} =
\\ ?I(t)S(t) ? (? + d2)E(t)$$"),
        helpText("Infectious $$\\frac{dI}{dt} =
\\ ?E(t) ? (? + d3)I(t)$$"),
        helpText("Removed $$\\frac{dR}{dt} =
\\   ?I(t) ? d4R(t)$$")
      )
    })
  }
)

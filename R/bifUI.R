bifUI <- function(state, parms, plotopts, numopts) {
  choices <- c("Time series" = 1, "1 parameter bifurcation" = 2, "2 parameter bifurcation" = 3)
  myTabs <- lapply(1:length(choices), function(i) {tabPanel(title = names(choices)[i], plotOutput(outputId=paste0("plot", choices[i]), height = "100%"), value = choices[i])})

  ui <- shinydashboardPlus::dashboardPage(
    shinyjs::useShinyjs(),
    ########## Header of window
    header = shinydashboardPlus::dashboardHeader(
      # Set height of dashboardHeader
      title = tagList(
        span(class = "logo-lg", "Bifurcation analysis"),
        icon("compass"), tags$style(".fa-compass {color:#E87722}")),
      leftUi = tagList(span(class = "help-button", icon("question-circle", verify_fa = FALSE),
                            tags$style(".fa-question-circle {font-size: 24px; color:#66CC66; left: 235px; top: 13px; position: fixed;}")),
                       tags$li(class = "dropdown", actionButton("showODEs", "Show ODEs", class = "show-odes"),
                               tags$style(".show-odes {font-size: 13px;
                                      border-width:2px;
                                      width: 85px; height: 30px;
                                      text-indent: -8px;
                                      left: 270px; top: 11px; position: fixed;}"))),
      titleWidth = 220,
      controlbarIcon = shiny::icon("cogs", verify_fa = FALSE)
    ),
    ########## Left side-bar
    sidebar = shinydashboardPlus::dashboardSidebar(
      width = 220,
      tags$style(type='text/css', "#computefwrd1 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computefwrd2 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computefwrd3 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computebwrd1 { font-size: 13px;}"),
      tags$style(type='text/css', "#computebwrd2 { font-size: 13px;}"),
      tags$style(type='text/css', "#computebwrd3 { font-size: 13px;}"),
      tags$head(
        tags$style(HTML("
          .selectize-dropdown-content .active {
            background: #2196f3 !important;
            color: white !important;
          }
  "))
      ),
      # conditionPanels inside sidebarMenu do not really work well, so the entire sidebarMenu should be wrapped inside a conditionaPanel
      conditionalPanel(
        condition = "input.plottab == 1",
        h4("Initial values", align = "center"),
        selectInput('selectpoint1', "", c("User specified" = 0), selected=0),
        sidebarMenu(
          id = "varsparscurvesmenu",
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_1"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tags$head(
              tags$style(
                HTML(
                  '
              .control-label {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              height: 16px !important;
              }
              .form-group {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .form-control {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 0px !important;
              border-top: 0 !important;
              border-bottom: 0px !important;
              height: 25px !important;
              }
              '))),
            tabName = "state1tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_1"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tags$head(
              tags$style(
                HTML(
                  '
              .control-label {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              height: 16px !important;
              }
              .form-group {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .form-control {
              padding-top: 0 !important;
              padding-bottom: 0px !important;
              margin-top: 0 !important;
              margin-bottom: 0px !important;
              border-top: 0 !important;
              border-bottom: 0px !important;
              height: 25px !important;
              }
              '))),
            tabName = "pars1tab"
          ),
          br(),
          splitLayout(cellWidths = c("52%", "48%"),
                      actionButton("computebwrd1", "Compute", icon("backward")),
                      actionButton("computefwrd1", "Compute", icon("forward"))),
          br(),
          h4("Curve management", align = "center"),
          menuItem(h4("Load curves"),
                   textInput("loadcurve1", "Give a valid R variable name"),
                   splitLayout(cellWidths = c("50%", "50%"),
                               actionButton("appendbtn1", "Append"),
                               actionButton("replacebtn1", "Replace")),
                   tabName = "loadtab1"
          ),
          menuItem(h4("Save curves"),
                   selectInput('savecurve1', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   textInput("curvename1", "Give a valid R variable name"),
                   actionButton("savebtn1", "Save curve"),
                   tabName = "savetab1"
          ),
          menuItem(h4("Delete curves"),
                   selectInput('deletecurve1', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   actionButton("deletebtn1", "Delete curve"),
                   tabName = "deletetab1"
          )
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 2",
        sidebarMenu(
          id = "varsparscurvesmenu",
          h4("Initial values", align = "center"),
          selectInput('selectpoint2', "", c("User specified" = 0), selected=0),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_2"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tabName = "state2tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_2"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tabName = "pars2tab"
          ),
          radioButtons("curvetype2", "Curve type to compute", choices = c("EQ", "LC"), inline = TRUE),
          br(),
          splitLayout(cellWidths = c("52%", "48%"),
                      actionButton("computebwrd2", "Compute", icon("backward")),
                      actionButton("computefwrd2", "Compute", icon("forward"))),
          br(),
          h4("Curve management", align = "center"),
          menuItem(h4("Load curves"),
                   textInput("loadcurve2", "Give a valid R variable name"),
                   splitLayout(cellWidths = c("50%", "50%"),
                               actionButton("appendbtn2", "Append"),
                               actionButton("replacebtn2", "Replace")),
                   tabName = "loadtab2"
          ),
          menuItem(h4("Save curves"),
                   selectInput('savecurve2', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   textInput("curvename2", "Give a valid R variable name"),
                   actionButton("savebtn2", "Save curve"),
                   tabName = "savetab2"
          ),
          menuItem(h4("Delete curves"),
                   selectInput('deletecurve2', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   actionButton("deletebtn2", "Delete curve"),
                   tabName = "deletetab2"
          )
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 3",
        sidebarMenu(
          id = "varsparscurvesmenu",
          h4("Initial values", align = "center"),
          selectInput('selectpoint3', "", c("User specified" = 0), selected=0),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_3"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tabName = "state3tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_3"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tabName = "pars3tab"
          ),
          radioButtons("curvetype3", "Curve type to compute", choices = c("BP", "HP", "LP"), inline = TRUE),
          br(),
          splitLayout(cellWidths = c("52%", "48%"),
                      actionButton("computebwrd3", "Compute", icon("backward")),
                      actionButton("computefwrd3", "Compute", icon("forward"))),
          br(),
          h4("Curve management", align = "center"),
          menuItem(h4("Load curves"),
                   textInput("loadcurve3", "Give a valid R variable name"),
                   splitLayout(cellWidths = c("50%", "50%"),
                               actionButton("appendbtn3", "Append"),
                               actionButton("replacebtn3", "Replace")),
                   tabName = "loadtab3"
          ),
          menuItem(h4("Save curves"),
                   selectInput('savecurve3', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   textInput("curvename3", "Give a valid R variable name"),
                   actionButton("savebtn3", "Save curve"),
                   tabName = "savetab3"
          ),
          menuItem(h4("Delete curves"),
                   selectInput('deletecurve3', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                   actionButton("deletebtn3", "Delete curve"),
                   tabName = "deletetab3"
          )
        )
      ),
      sidebarMenu(
        id = "restmenu",
        radioButtons("reportlevel", "Console report level", choiceNames = c("None", "Normal", "Full"), choiceValues = (0:2), inline = TRUE),
        div(style="line-height: 18px !important", br()),
        conditionalPanel(
          condition = "input.plottab == 2 | input.plottab == 3",
          shinyjs::hidden(actionButton("pausebtn", "Pause", icon("pause-circle", verify_fa = FALSE), style = "color: white;
                       background-color: green;
                       font-size: 16px;
                       position: relative;
                       left: 18px;
                       height: 45px;
                       width: 150px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 6px;
                       border-width: 2px")),
          div(style="line-height: 8px !important", br()),
          shinyjs::hidden(actionButton("stopbtn", "Stop", icon("stop-circle", verify_fa = FALSE), style = "color: white;
                       background-color: red;
                       font-size: 16px;
                       position: relative;
                       left: 18px;
                       height: 45px;
                       width: 150px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 6px;
                       border-width: 2px"))
        )
      )
    ),
    ########## Main panels
    body = shinydashboard::dashboardBody(
      tags$script(
        "
        Shiny.addCustomMessageHandler('scrollCallback',
        function(color) {
        var objDiv = document.getElementById('console');
        objDiv.scrollTop = objDiv.scrollHeight;
        }
        );
        // Bind function to respond to the help button
        $('i.fa.fa-question-circle').on('click',function(){
          Shiny.onInputChange('helpClicked', Math.random());
        });
        // Bind function to the toggle sidebar button
        $('i.fa.fa-cogs').on('click',function(){
          setTimeout( function() { $(window).trigger('resize')}, 200); // Trigger resize event
        });
        var dimension = [0, 0];
        $(document).on('shiny:connected', function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange('dimension', dimension);
        });
        $(window).resize(function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange('dimension', dimension);
        });
        "
      ),
      # tags$head(
      #   tags$style("body { min-height: 611px; height: auto; min-width: 800px; margin: auto; }")
      # ),
      do.call(tabsetPanel, c(myTabs, id = "plottab")),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#saveplot { width: 105px; left: calc(100% - 125px); top: 107px; position: absolute;}"
      ))),
      downloadButton("saveplot", label = paste0("Save ", plotopts[[1]]$saveplotas)),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#console { font-size: 11px; width: calc(100%); left: calc(242px); height: 149px; overflow: auto; }"
      ))),
      tags$head(tags$style(HTML('.box {margin-bottom: 0px; margin-top: 0px;}'))),
      shinydashboard::box(width = NULL, height = "170px", verbatimTextOutput("console", placeholder = TRUE)),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#progress { font-size: 10px; width: calc(100%); left: calc(242px); height: 110px; overflow: auto;}"
      ))),
      tags$head(tags$style(HTML('.box {margin-bottom: 0px; margin-top: 0px;}'))),
      shinydashboard::box(width = NULL, height = "130px", verbatimTextOutput("progress", placeholder = TRUE))
    ),
    ########## Right side-bar
    controlbar = shinydashboardPlus::dashboardControlbar(
      controlbarMenu(
        id = "controlbartabs",
        controlbarItem(
          NULL,
          div(style="font-size: 20px; line-height: 22px; margin-top: 0px; margin-bottom: 5px;",
              ("Plot options")),
          selectInput('xcol', 'Variable(s) on X-axis',
                      c("Time" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$xcol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="xmin", label="Minimum", value=plotopts[[1]]$xmin),
                      numericInput(inputId="xmax", label="Maximum", value=plotopts[[1]]$xmax)),
          selectInput('logx', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
          div(style="line-height: 6px !important", br()),
          selectInput('ycol', 'Variable(s) on Y-axis',
                      c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$ycol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="ymin", label="Minimum", value=plotopts[[1]]$ymin),
                      numericInput(inputId="ymax", label="Maximum", value=plotopts[[1]]$ymax)),
          selectInput('logy', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
          div(style="line-height: 6px !important", br()),
          conditionalPanel(
            condition = "input.plottab != 3 && input.ycol > 1",
            selectInput('y2col', 'Variable on 2nd Y-axis',
                        c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$y2col),
            conditionalPanel(
              condition = "input.y2col > 1",
              splitLayout(cellWidths = c("50%", "50%"),
                          numericInput(inputId="y2min", label="Minimum", value=plotopts[[1]]$y2min),
                          numericInput(inputId="y2max", label="Maximum", value=plotopts[[1]]$y2max)),
              selectInput('logy2', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy2),
              splitLayout(cellWidths = c("50%", "50%"),
                          strong("Plot type"),
                          radioButtons("plot3d", NULL, choices = c("2D" = 0, "3D" = 1), selected = plotopts[[1]]$plot3d,
                                       inline = TRUE)),
              conditionalPanel(condition = "input.plot3d == 1",
                               sliderInput(inputId="theta", label="Viewing angle",
                                           min=-90, max=90, value=plotopts[[1]]$theta, step=1,
                                           ticks = FALSE, round=TRUE))
            )
          ),
          div(style="line-height: 12px !important", br()),
          actionButton("plotoptsapply", "Apply", icon("sync", verify_fa = FALSE)),
          tags$head(
            tags$style(
              HTML(
                '
              .form-group {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .shiny-split-layout {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              label {font-size: 13px;}
              #xmin{height: 30px}
              #xmax{height: 30px}
              #ymin{height: 30px}
              #ymax{height: 30px}
              #y2min{height: 30px}
              #y2max{height: 30px}
              '))),
          value = "control-sidebar-plotopttab-tab",
          icon = shiny::icon("chart-line")),
        controlbarItem(
          NULL,
          div(style="font-size: 20px; line-height: 22px; margin-top: 0px; margin-bottom: 5px;",
              ("Numerical options")),
          conditionalPanel(
            condition = "input.plottab == 1",
            h4("Time integration"),
            splitLayout(cellWidths = c("50%", "50%"),
                        numericInput(inputId="tmax", label="Maximum time", value=numopts$tmax),
                        numericInput(inputId="tstep", label="Time step", value=numopts$tstep, min=0, step=0.1)),
            selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4"),selected=numopts$odemethod)
          ),
          conditionalPanel(
            condition = "input.plottab > 1",
            h4("Curve continuation"),
            div(style="font-size: 16px; line-height: 0px; margin-top: 20px; margin-bottom: 12px; !important", ("Tolerances")),
            splitLayout(cellWidths = c("45%", "55%"),
                        textInput(inputId="rhstol", label="Function", value=sprintf("%.1E", numopts$rhstol)),
                        textInput(inputId="dytol", label="Variable", value=sprintf("%.1E", numopts$dytol))),
            splitLayout(cellWidths = c("45%", "55%"),
                        textInput(inputId="iszero", label="Zero identity", value=sprintf("%.1E", numopts$iszero)),
                        textInput(inputId="neartol", label="Neighbourhood", value=sprintf("%.1E", numopts$neartol))),
            textInput(inputId="jacdif", label="Jacobian perturbation", value=sprintf("%.1E", numopts$jacdif)),
            div(style="font-size: 16px; line-height: 0px; margin-top: 24px; margin-bottom: 12px; !important", ("Step size")),
            splitLayout(cellWidths = c("33%", "33%", "34%"),
                        textInput(inputId="initstepsize", label="Initial", value=sprintf("%.1E", numopts$initstepsize)),
                        textInput(inputId="minstepsize", label="Minimum", value=sprintf("%.1E", numopts$minstepsize)),
                        textInput(inputId="maxstepsize", label="Maximum", value=sprintf("%.1E", numopts$maxstepsize))),
            div(style="font-size: 16px; line-height: 0px; margin-top: 24px; margin-bottom: 12px; !important", ("Iterations")),
            numericInput(inputId="maxiter", label="Maximum iterations", value=numopts$maxiter),
            numericInput(inputId="maxpoints", label="Maximum number of points", value=numopts$maxpoints),
            numericInput(inputId="replotfreq", label="Points between plot updates", min=1, max=10000,
                         value=numopts$replotfreq, step=1),
            div(style="font-size: 16px; line-height: 0px; margin-top: 24px; margin-bottom: 12px; !important", ("Limit cycle continuation")),
            splitLayout(cellWidths = c("45%", "55%"),
                        numericInput(inputId="ninterval", label="Intervals", value=numopts$ninterval, min = 1, max = 40, step = 1),
                        numericInput(inputId="glorder", label="Order", value=numopts$glorder), min = 2, max = 7, step = 1),
            textInput(inputId="lcampl", label="Initial amplitude", value=sprintf("%.1E", numopts$lcampl))
          ),
          div(style="line-height: 12px !important", br()),
          actionButton("numoptsapply", "Apply", icon("sync", verify_fa = FALSE)),
          tags$head(
            tags$style(
              HTML(
                '
              .form-group {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .shiny-split-layout {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .shiny-bound-input {
              font-size: 12px !important;
              padding-left: 5px !important;
              padding-right: 2px !important;
              }
              label {font-size: 13px;}
              #tmax{height: 30px}
              #tstep{height: 30px}
              #rhstol{height: 30px}
              #dytol{height: 30px}
              #neartol{height: 30px}
              #iszero{height: 30px}
              #jacdif{height: 30px}
              #initstepsize{height: 30px}
              #minstepsize{height: 30px}
              #maxstepsize{height: 30px}
              #maxiter{height: 30px}
              #maxpoints{height: 30px}
              #replotfreq{height: 30px}
              #ninterval{height: 30px}
              #glorder{height: 30px}
              #lcampl{height: 30px}
              '))),
          value = "control-sidebar-numopttab-tab",
          icon = shiny::icon("tachometer-alt", verify_fa = FALSE)),
        selected = "control-sidebar-plotopttab-tab"),
      id = "controlbar",
      skin = "dark"),
    footer = NULL
  )

  return(ui)
}

#' Phaseplane analysis of a system of ODEs
#'
#' \code{phaseplane}
#'
#'
#'   phaseplane(model, state, parms, resume = TRUE, ...)
#'
#'
#' @param   model  (function, required)
#' \preformatted{}
#'               An R-function that computes the values of the derivatives
#'               in the ODE system (the model definition) at time t.
#'               The model must be defined as: model <- function(t, state, parms),
#'               where t is the current time point in the integration, state is
#'               the current value of the variables in the ODE #' system and
#'               parms is a vector or list of parameters.
#'               The return value of func should be a list, whose first and single
#'               element is a vector containing the derivatives of y with respect
#'               to time. The derivatives must be specified in the same order as
#'               the state variables state. The vector state and parms should both
#'               have name attributes for all their elements
#'
#' @param   state  (numeric vector, required)
#' \preformatted{}
#'               The initial (state) values for the ODE system. This vector should
#'               have name attributes for all its elements
#'
#' @param   parms  (numeric vector, required)
#' \preformatted{}
#'               The values of the parameters in the ODE system. This vector should
#'               have name attributes for all its elements
#'
#' @param   resume  (boolean, optional)
#' \preformatted{}
#'               If TRUE the program will try to load the curves computed during
#'               the last session from the global variable '<model>PhaseCurves' and try
#'               to restore the numerical and plot settings by importing them from
#'               the global variable '<model>PhaseSettings', where the substring
#'               '<model>' is the name of the function describing the dynamics, which
#'               is passed as first argument to 'bifurcation()'.
#'               The program saves the curves computed during a session and the
#'               numerical and plot settings of this last session in these global
#'               variables '<model>PhaseCurves' and '<model>PhaseSettings'.
#'
#' @param   ...  (optional arguments)
#' \preformatted{}
#'               Additional arguments that can be included at the command line to tweak
#'               graphical default values used by the application.
#'               Valid arguments are:
#' \preformatted{}
#'               \code{lwd}: Line width (default 2)
#' \preformatted{}
#'               \code{cex}: Base font size (default 1.2)
#' \preformatted{}
#'               \code{tcl.len}:     Length of axes ticks (default 0.03)
#'
#' \preformatted{}
#'               \code{saveplotas}:  Possible values: "pdf" or "png" (default).
#'                                   Save plot to PDF or PNG file.
#'
#' @return None.
#'
#' @examples
#' if(interactive()){
#' # The initial state of the system has to be specified as a named vector of state values.
#' state <- c(R=1, N=0.01)
#'
#' # Parameters has to be specified as a named vector of parameters.
#' parms <- c(r=1, K=1, a=1, c=1, delta=0.5)
#'
#' # The model has to be specified as a function that returns
#' # the derivatives as a list.
#' model <- function(t, state, parms) {
#'   with(as.list(c(state,parms)), {
#'
#'     dR <- r*R*(1 - R/K) - a*R*N
#'     dN <- c*a*R*N - delta*N
#'
#'    # The order of the derivatives in the returned list has to be
#'    # identical to the order of the state variables contained in
#'    # the argument "state"
#'     return(list(c(dR, dN)))
#'   })
#' }
#'
#' phaseplane(model, state, parms)
#' }
#' @import deSolve rootSolve
#' @import shiny
#' @importFrom shinydashboard dashboardBody box menuItem sidebarMenu
#' @importFrom shinydashboardPlus dashboardPage dashboardHeader dashboardSidebar dashboardControlbar controlbarItem controlbarMenu
#' @importFrom graphics contour legend lines par plot points text title axis mtext persp axTicks segments
#' @importFrom grDevices trans3d dev.off png pdf
#' @importFrom shinyjs useShinyjs click removeClass html
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom utils browseURL capture.output unzip
#' @export
phaseplane <- function(model, state, parms, resume = TRUE, ...) {

  modelname <- as.list(match.call())[[2]]
  savedSettingsName <- paste0(modelname, "PhaseSettings")
  savedCurvesName <- paste0(modelname, "PhaseCurves")
  if (interactive()) {
    if (length(unlist(model(0, state, parms))) != length(state))
      stop("The number of derivatives returned by the model function must equal the length of the state vector")

    # Get the names of the state variables and the parameters
    statenames <- names(state)
    parmsnames <- names(parms)

    # Initialize numerical options
    initnopts <- list(odemethod = "lsoda", tmax = 100, tstep = 0.1, eps = -.001, pgrid=5, ssgrid=10, npixels=200)

    # Initialize options for plotting etc.
    initpopts <- vector(mode = "list", 2)
    for (i in 1:length(initpopts)) {
      initpopts[[i]] <- list(xcol = 1, xmin = 0, xmax = 1, logx = 0, xlab = "",
                             ycol = 1, ymin = 0, ymax = 1, logy = 0, ylab = "",
                             y2col = 1, y2min = 0, y2max = 1, logy2 = 0, y2lab = "None", plot3d = 0,
                             plotmar = c(5,5,4,4), lwd = 2, pch = 20, tcl.len = 0.03, theta = -35,
                             cex = 1.2, cex.lab = 1.25, cex.axis = 1, cex.legend = 1,
                             colors = c("red","blue","darkgreen","darkorange","darkmagenta", "gold","darkorchid",
                                        "aquamarine","deeppink","gray",seq(2,991)),
                             saveplotas = "png")
    }
    initpopts[[2]]$ycol <- 2

    # Read options from the environment
    if (resume && exists(savedSettingsName, envir = .GlobalEnv)) {
      inlist    <- get(savedSettingsName, envir = .GlobalEnv)
      initnopts <- phaseCheckNumSettings(initnopts, inlist)
      initpopts <- phaseCheckPlotSettings(initpopts, inlist, state, parms)
    }

    # Read options from the command line
    adjustableopts <- c("lwd", "cex", "tcl.len")
    dots <- list(...)
    if (!is.null(dots)) {
      useropts <- dots[names(dots) %in% adjustableopts]
      if (!is.null(useropts)) {
        for (j in 1:length(useropts)) {
          for (i in 1:length(initpopts)) initpopts[[i]][names(useropts)[j]] <- useropts[j]
        }
      }
      if (("saveplotas" %in% names(dots)) && ((dots["saveplotas"] == "png") || (dots["saveplotas"] == "pdf"))) {
        for (i in 1:length(initpopts)) initpopts[[i]]["saveplotas"] <- dots["saveplotas"]
      }
    }
    initpopts[[1]]$xlab  <- "Time"
    initpopts[[1]]$xmax  <- initnopts$tmax
    initpopts[[1]]$ylab  <- ifelse(initpopts[[1]]$ycol  == 1, "State variables", statenames[initpopts[[1]]$ycol -1])
    initpopts[[1]]$y2lab <- ifelse(initpopts[[1]]$y2col == 1, "None",            statenames[initpopts[[1]]$y2col-1])
    initpopts[[2]]$xlab  <- statenames[initpopts[[2]]$xcol]
    initpopts[[2]]$ylab  <- statenames[initpopts[[2]]$ycol]
    initpopts[[2]]$y2lab <- ifelse(initpopts[[2]]$y2col == 1, "None",            statenames[initpopts[[2]]$y2col-1])

    # Read the curves from the environment
    initCurves <- list(Orbits = list(), TotalCurves = 0)
    if (resume && exists(savedCurvesName, envir = .GlobalEnv)) {
      inlist     <- get(savedCurvesName, envir = .GlobalEnv)
      initCurves <- phaseCheckInputCurves(NULL, inlist, statenames, parmsnames)
    }

    ui <- phaseUI(state, parms, initpopts, initnopts)

    # server = function(input, output) { }
    server <- function(input, output, session) {

      # Hide the button to collapse the left sidebar
      shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")

      curstate   <- state
      curparms   <- parms

      # Create the variable curveList as a reactive value, such that the plots will be updated when curveList changes
      curveList <- reactiveValues()
      curveList$Orbits <- initCurves$Orbits
      curveList$TotalCurves <- initCurves$TotalCurves
      curveListNames <- c('Orbits', 'TotalCurves')

      # Create the variable plotopts as a reactive value, such that the plots will be updated when plotopts changes
      plotopts <- reactiveValues()
      plotopts$Orbits <- initpopts[[1]]
      plotopts$PhasePlane <- initpopts[[2]]

      # Create the variable numopts as a reactive value, even though computations will only start after a button push
      numopts <- reactiveValues()
      lapply((1:length(initnopts)), function(i) {numopts[[(names(initnopts)[[i]])]] <- initnopts[[i]]})

      # Create the variable Derivatives
      Derivatives <- NULL

      # Create the variable consoleLog as a reactive value, such that the console will be updated when text is added to it
      consoleLog <- reactiveVal()

      # Create curveDirection, updatePlot and changeCurveMenu as reactive values
      curveDirection <- reactiveVal(0)
      updatePlot <- reactiveVal(0)
      changeCurveMenu <- reactiveVal(1)

      # React to changes in updatePlot by replotting the current plot
      observe({
        if (as.numeric(updatePlot()) == 0) return(NULL)

        isolate({
          curtab <- as.numeric(input$plottab)
          clist <- reactiveValuesToList(curveList)
          popts <- reactiveValuesToList(plotopts)
          nopts <- reactiveValuesToList(numopts)
          output[[paste0("plot", curtab)]] <- renderPlot({
            # Save plot options to restore on exit
            oldpar <- par(no.readonly = TRUE)
            on.exit(par(oldpar))

            if (curtab == 1) bifOrbitplot(session, clist$Orbits, popts$Orbits)
            else {
              if (length(curstate) == 1) {
                phasePlot1D(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts)
              }
              else {
                Derivatives <<- phasePlot2D(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts, zlst = Derivatives)
              }
              if (curtab >= 3) {
                msg <- allequi(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts)
                consoleLog(msg)
              }
              if ((curtab == 5) && (length(clist$Orbits) > 0)) {
                lapply((1:length(clist$Orbits)), function(i) {
                  if ((popts$PhasePlane$xlab %in% statenames) && (popts$PhasePlane$ylab %in% statenames) &&
                      (all(curparms == clist$Orbits[[i]]$parameters))) {
                    xcol <- 1 + match(popts$PhasePlane$xlab, statenames)
                    ycol <- 1 + match(popts$PhasePlane$ylab, statenames)
                    lines(clist$Orbits[[i]]$points[,xcol], clist$Orbits[[i]]$points[,ycol],
                          col="black", lwd=popts$PhasePlane$lwd)
                  }
                })
              }
            }
          },
          height = function() {setPlotHeight(session, input)},
          width = function() {setPlotWidth(session, input)})

          updatePlot(0)

          # Update the console log
          consoleLog(session$userData$alltext)
        })
      })

      # Initialise a computation. Triggered by a change in the reactive value curveDirection()
      observe({
        if (as.numeric(curveDirection()) == 0) return(NULL)

        isolate({
          # Close the rightSidebar
          shinyjs::removeClass(id = "controlbar", class = "control-sidebar-open")
          # Collapse the State variables and Parameters stacks
          # shinyjs::removeClass(selector = "li.treeview", class = "active")
          # shinyjs::hide(selector = "ul.menu-open");
          updateSelectInput(session, "deletecurve1", selected = 0)

          curtab <- as.numeric(input$plottab)
          clist <- reactiveValuesToList(curveList)
          nopts <- reactiveValuesToList(numopts)

          pointid <- as.numeric(input[['selectpoint1']])

          for (i in statenames) curstate[i] <<- input[[paste0(i, "_1")]]
          for (i in parmsnames) curparms[i] <<- input[[paste0(i, "_1")]]

          nopts$tstep <- as.numeric(curveDirection())*abs(nopts$tstep)

          newlist <- computeTimeseries(session, model, curstate, curparms, clist, pointid, nopts)
          if (!is.null(newlist)) lapply((1:length(newlist)),
                                        function(i) {curveList[[(curveListNames[[i]])]] <- newlist[[(curveListNames[[i]])]]})
          changeCurveMenu(1)
          updatePlot(1)
          curveDirection(0)
        })
      })


      # React to changes in curveList by updating the special point selection menu
      observe({
        if (as.numeric(changeCurveMenu()) == 0) return(NULL)
        else if (as.numeric(changeCurveMenu()) == 1)
          updateSpecialPointsList(session, reactiveValuesToList(curveList), as.numeric(isolate(input$selectpoint1)))
        else
          updateSpecialPointsList(session, reactiveValuesToList(curveList), 0)

        # Updating the save and delete curve menu
        updateCurveMenu(session, curveList$Orbits, 1)
        changeCurveMenu(0)
      })

      # React to changes in consoleLog, by rendering the new text
      observe({
        shinyjs::html(id = "console", html = HTML(gsub("\n", "<br>", consoleLog())))
        session$sendCustomMessage(type = "scrollCallback", 1)
      })

      # Handle requests to save the plot as an image file
      output$saveplot <- downloadHandler(
        filename = function() {
          tempfile(pattern = "Rplot", tmpdir = getwd(), fileext = paste0(".", initpopts[[1]]$saveplotas))
        },
        content = function(file) {
          curtab <- as.numeric(input$plottab)
          clist <- reactiveValuesToList(curveList)
          popts <- reactiveValuesToList(plotopts)
          nopts <- reactiveValuesToList(numopts)
          filetype <- tools::file_ext(file)
          destfile <- paste0(tools::file_path_sans_ext(file), ".", filetype)
          if (filetype == "png") {
            png(destfile,
                height = setPlotHeight(session, input),
                width = setPlotWidth(session, input))
          } else {
            pdf(file = destfile, onefile = T,
                height = 7 * setPlotHeight(session, input) / setPlotWidth(session, input),
                width = 7)
          }
          # Save plot options to restore on exit
          oldpar <- par(no.readonly = TRUE)
          on.exit(par(oldpar))

          if (curtab == 1) bifOrbitplot(session, clist$Orbits, popts$Orbits)
          else {
            if (length(curstate) == 1) {
              phasePlot1D(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts)
            }
            else {
              Derivatives <<- phasePlot2D(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts, zlst = Derivatives)
            }
            if (curtab >= 3) msg <- allequi(curtab, odes=model, state=curstate, parms=curparms, plotopts=popts$PhasePlane, numopts = nopts)
            if (curtab == 3) output[["console"]] <- renderText({msg})

            if ((curtab == 5) && (length(clist$Orbits) > 0)) {
              lapply((1:length(clist$Orbits)), function(i) {
                if ((popts$PhasePlane$xlab %in% statenames) && (popts$PhasePlane$ylab %in% statenames) &&
                    (all(curparms == clist$Orbits[[i]]$parameters))) {
                  xcol <- 1 + match(popts$PhasePlane$xlab, statenames)
                  ycol <- 1 + match(popts$PhasePlane$ylab, statenames)
                  lines(clist$Orbits[[i]]$points[,xcol], clist$Orbits[[i]]$points[,ycol],
                        col="black", lwd=popts$PhasePlane$lwd)
                }
              })
            }
          }
          dev.off()
        })

      # Respond to a change in plot tabs
      observeEvent(input$plottab, {
        curtab <- as.numeric(input$plottab)
        popts <- reactiveValuesToList(plotopts)
        # Update the plot options
        if (curtab == 1) {
          popts <- popts$Orbits
          updateSelectInput(session, "xcol",  label='Variable(s) on X-axis', choices=c("Time" = 1), selected=popts$xcol)
          if (length(curstate) == 1)
            updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis', choices=c(setNames(1, statenames)), selected=popts$ycol)
          else {
            updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis', choices=c("All" = 1, setNames((2:(length(statenames)+1)), statenames)), selected=popts$ycol)
            updateSelectInput(session, "y2col", label='Variable on 2nd Y-axis',  choices=c("None" = 1, setNames((2:(length(statenames)+1)), statenames)), selected=popts$y2col)
            updateSelectInput(session,  "logy2", selected=popts$logy2)
            updateNumericInput(session, "y2min", value=popts$y2min)
            updateNumericInput(session, "y2max", value=popts$y2max)
            popts$y2lab <- (c("None", statenames))[popts$y2col]
          }
        } else {
          popts <- popts$PhasePlane
          updateSelectInput(session, "xcol",  label='Variable(s) on X-axis', choices=c(setNames((1:(length(statenames))), statenames)), selected=popts$xcol)
          if (length(curstate) == 1)
            updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis', choices=c(setNames(1, paste0("d", statenames[1], "/dt"))), selected=popts$ycol)
          else
            updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis', choices=c(setNames((1:(length(statenames))), statenames)), selected=popts$ycol)
        }
        updateSelectInput(session,  "logx", selected=popts$logx)
        updateNumericInput(session, "xmin", value=popts$xmin)
        updateNumericInput(session, "xmax", value=popts$xmax)
        popts$xlab <- ifelse(curtab == 1, (c("Time", statenames))[popts$xcol], statenames[popts$xcol])

        updateSelectInput(session,  "logy", selected=popts$logy)
        updateNumericInput(session, "ymin", value=popts$ymin)
        updateNumericInput(session, "ymax", value=popts$ymax)
        if (length(curstate) == 1)
          popts$ylab <- paste0("d", statenames[1], "/dt")
        else
          popts$ylab <- ifelse(curtab == 1, (c("State variables", statenames))[popts$ycol], statenames[popts$ycol])
        updatePlot(1)
      })

      # React to selection of an initial point
      observeEvent(input$selectpoint1, {
        # In contrast to bifurcation() there is only a single button 'selectpoint1' and
        # all curves are stored under the 'Orbits' element of curveList, hence curtab
        # is always 1
        curtab <- 1
        clist <- reactiveValuesToList(curveList)
        updateSelectedPoint(session, curtab, clist, as.numeric(input[[paste0('selectpoint', curtab)]]),
                            statenames, parmsnames)
      })

      # Apply newly entered state and parameter values
      observeEvent(input$lapply, {
        isolate({
          curtab <- as.numeric(input$plottab)
          oldpars <- curparms
          for (i in statenames) curstate[i] <<- input[[paste0(i, "_1")]]
          for (i in parmsnames) curparms[i] <<- input[[paste0(i, "_1")]]

          if ((length(curstate) > 1) && (any(oldpars != curparms))) Derivatives <<- NULL

          updatePlot(1)
        })
      })

      # Apply newly entered plot settings
      observeEvent(input$plotoptsapply, {
        curtab <- as.numeric(input$plottab)
        curtabname <- ifelse(curtab == 1, "Orbits", "PhasePlane")
        popts <- reactiveValuesToList(plotopts)
        if (curtab == 1) {
          xlabs  <- "Time"
          ylabs  <- c("State variables", statenames)
        } else {
          xlabs  <- statenames
          ylabs  <- statenames
        }
        y2labs <- c("None", statenames)

        plotopts[[curtabname]] <- processPlotOptionsApply(session, input, curtab, popts[[curtabname]])

        plotopts[[curtabname]]$xlab <- ifelse(curtab == 1, "Time", statenames[plotopts[[curtabname]]$xcol])
        plotopts[[curtabname]]$ylab <- ifelse(curtab == 1, (c("State variables", statenames))[plotopts[[curtabname]]$ycol],
                                              statenames[plotopts[[curtabname]]$ycol])
        plotopts[[curtabname]]$y2lab <- (c("None", statenames))[plotopts[[curtabname]]$y2col]

        if ((curtabname == "PhasePlane") && (length(curstate) > 1)) Derivatives <<- NULL

        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Apply newly entered numerical settings
      observeEvent(input$numoptsapply, {
        processNumOptionsApply(session, input, 1, numopts)

        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Initiate a forward computation
      observeEvent(input$computefwrd1, {
        curveDirection(1)
      }, ignoreInit = TRUE)

      # Initiate a backward computation
      observeEvent(input$computebwrd1, {
        curveDirection(-1)
      }, ignoreInit = TRUE)

      # Delete one or more curves
      observeEvent(input$deletebtn1, {
        totalcurves <- as.numeric(length((curveList$Orbits)))
        deletenr <- as.numeric(input$deletecurve1)

        if ((totalcurves > 0) && ((deletenr > 0) || (deletenr == -1)) && (deletenr < (totalcurves + 1))) {
          if (deletenr == -1) {
            msg <- paste0("All curves deleted\n")
            curveList$Orbits <- list()
          } else {
            msg <- paste0("Curve '", curveList$Orbits[[deletenr]]$label, "' deleted\n")
            curveList$Orbits[[deletenr]] <- NULL
          }
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
        }
        changeCurveMenu(-1)
        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Save one or more curves
      observeEvent(input$savebtn1, {
        totalcurves <- as.numeric(length((curveList$Orbits)))
        savenr <- as.numeric(input$savecurve1)
        varname <- make.names(input$curvename1, unique = TRUE)

        if ((totalcurves > 0) && ((savenr > 0) || (savenr == -1)) && (savenr < (totalcurves + 1))) {
          if (exists(varname, envir = .GlobalEnv)) {
            rm(list = varname, envir = .GlobalEnv)
          }
          if (savenr == -1) {
            msg <- paste0("All curves saved to '", varname, "'\n")
            (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(varname, curveList$Orbits, 1L)
          } else {
            msg <- paste0("Curve '", curveList$Orbits[[savenr]]$label, "' saved to '", varname, "'\n")
            (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(varname, curveList$Orbits[[savenr]], 1L)
          }
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
        }

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Append one or more curves to the current curve list
      observeEvent(input$appendbtn1, {
        curvarname <- input$loadcurve1
        # Check whether the name is a valid R variable name
        if (make.names(curvarname) != curvarname) return(NULL)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, curvarname, statenames, parmsnames, replace = FALSE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}
        changeCurveMenu(-1)
        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Replace the current curve list with a stored list
      observeEvent(input$replacebtn1, {
        curvarname <- input$loadcurve1
        # Check whether the name is a valid R variable name
        if (make.names(curvarname) != curvarname) return(NULL)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, curvarname, statenames, parmsnames, replace = TRUE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}
        changeCurveMenu(-1)
        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Show the manual
      observeEvent(input$helpClicked, {
        # browseURL(paste0(system.file("manual", package = "deBif"), "/index.html"))
        oldwd <- getwd()
        on.exit(setwd(oldwd))
        tempDir <- tempdir()
        unlink(paste0(tempDir, "/manual"), recursive = TRUE)
        dir.create(paste0(tempDir, "/manual"))
        setwd(paste0(tempDir, "/manual"))
        unzip(paste0(system.file("manual", package = "deBif"), "/deBif-manual.zip"))
        setwd(oldwd)
        htmlFile <- file.path(tempDir, "manual/index.html")
        browseURL(paste0("file:/", htmlFile))
      })

      # Show the model
      observeEvent(input$showODEs, {
        systxt <- capture.output({print(model)})
        systxt <- paste0(systxt[1:(length(systxt))], collapse = "<br>")
        odestxt <- paste0("<b>State variables:</b><br>",
                          paste(names(state), ":", state, "<br>", collapse = ""),
                          "<br><b>Parameters:</b><br>",
                          paste(names(parms), ":", parms, "<br>", collapse = ""),
                          "<br><b>System of ODEs:</b><br>",
                          systxt, collapse = "")
        showModal(modalDialog(title = "System of ODEs", HTML(odestxt)))
      })

      # Actions to be carried out when the app is stopped
      onStop(fun = function() {
        isolate({
          cat("Saving curves and programs settings")
          # Save the current curve list
          if (exists(savedCurvesName, envir = .GlobalEnv)) {
            rm(list = savedCurvesName, envir = .GlobalEnv)
          }

          # global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
          (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(savedCurvesName,
                                                                               list(Orbits = curveList$Orbits, TotalCurves = curveList$TotalCurves), 1L)
          # Save the plot and numerical settings in the global environment
          if (exists(savedSettingsName, envir = .GlobalEnv)) {
            rm(list = savedSettingsName, envir = .GlobalEnv)
          }

          # global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
          (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(savedSettingsName,
                                                                               list(plotopts = reactiveValuesToList(plotopts), numopts = reactiveValuesToList(numopts)), 1L)
        })
      })
    }

    shinyApp(ui = ui, server = server, options = list(width = 900, height = 900))
  }
}


#' Phaseplane analysis of a system of ODEs
#'
#' \code{bifurcation}
#'
#'
#'   bifurcation(model, state, parms, resume = TRUE, ...)
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
#'               the last session from the global variable '<model>BifCurves' and try
#'               to restore the numerical and plot settings by importing them from
#'               the global variable '<model>BifSettings', where the substring
#'               '<model>' is the name of the function describing the dynamics, which
#'               is passed as first argument to 'bifurcation()'.
#'               The program saves the curves computed during a session and the
#'               numerical and plot settings of this last session in these global
#'               variables '<model>BifCurves' and '<model>BifSettings'.
#'
#' @param   ...  (optional arguments)
#' \preformatted{}
#'               Additional arguments that can be included at the command line to tweak
#'               graphical default values used by the application.
#'               Valid arguments are:
#' \preformatted{}
#'               \code{lwd}:         Line width (default 3)
#' \preformatted{}
#'               \code{cex}:         Base font size (default 1.2)
#' \preformatted{}
#'               \code{tcl.len}:     Length of axes ticks (default 0.03)
#' \preformatted{}
#'               \code{bifsym}:      Symbol used to mark a bifurcation point
#'                                     in an equilibrium curve (default: 8)
#' \preformatted{}
#'               \code{biflblpos}:   Location of label of a bifurcation point. Values
#'                                   of 1, 2, 3 and 4, respectively, indicate positions
#'                                   below, to the left of, above and to the right of
#'                                   the symbol marking the bifurcation point (default: 3)
#' \preformatted{}
#'               \code{unstablelty}: Line style of curve section representing unstable
#'                                   equilibrium points (default: 3 (refers to dotted lines))
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
#' bifurcation(model, state, parms)
#' }
#' @useDynLib deBif
#' @import deSolve rootSolve shiny
#' @importFrom shinydashboard dashboardBody box menuItem sidebarMenu
#' @importFrom shinydashboardPlus dashboardPage dashboardHeader dashboardSidebar dashboardControlbar controlbarItem controlbarMenu
#' @importFrom graphics contour legend lines par plot points text title axis mtext persp axTicks segments
#' @importFrom grDevices trans3d dev.off png pdf
#' @importFrom shinyjs useShinyjs click removeClass html
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom utils browseURL capture.output unzip
#' @export
bifurcation <- function(model, state, parms, resume = TRUE, ...) {

  modelname <- as.list(match.call())[[2]]
  savedSettingsName <- paste0(modelname, "BifSettings")
  savedCurvesName <- paste0(modelname, "BifCurves")
  if (interactive()) {
    if (length(unlist(model(0, state, parms))) != length(state))
      stop("The number of derivatives returned by the model function must equal the length of the state vector")

    # Get the names of the state variables and the parameters
    statenames <- names(state)
    parmsnames <- names(parms)

    # Initialize numerical options
    initnopts <- list(odemethod = "lsoda", tmax = 1000, tstep = 0.1,
                      args_run = unique(names(c(formals(deSolve::ode), formals(deSolve::lsoda)))),
                      methods_run = as.character(formals(deSolve::ode)$method),
                      rhstol = 1e-7, dytol = 1e-7, neartol = 0.05, jacdif = 1.0E-4, maxiter = 20,
                      maxpoints = 1000, iszero = 1.0E-5, minstepsize = 1.0E-3, maxstepsize = 0.05,
                      replotfreq = 10, ninterval = 10, glorder = 4, lcampl = 1.0E-3
    )

    # Initialize options for plotting etc.
    initpopts <- vector(mode = "list", 3)
    for (i in 1:length(initpopts)) {
      initpopts[[i]] <- list(xcol = 1, xmin = 0, xmax = 1, logx = 0, xlab = "",
                             ycol = 1, ymin = 0, ymax = 1, logy = 0, ylab = "",
                             y2col = 1, y2min = 0, y2max = 1, logy2 = 0, y2lab = "None", plot3d = 0,
                             lwd = 3, bifsym = 8, unstablelty = 3, tcl.len = 0.03, theta = -35,
                             cex = 1.2, cex.lab = 1.25, cex.axis = 1, cex.legend = 1, cex.sym = 1, biflblpos = 3,
                             colors = c("red","blue","darkgreen","darkorange","darkmagenta",
                                        "gold","darkorchid","aquamarine","deeppink","gray",seq(2,991)),
                             saveplotas = "png")
    }
    initpopts[[1]]$xmax <- initnopts$tmax
    initpopts[[3]]$ycol <- 2

    # Read options from the environment
    if (resume && exists(savedSettingsName, envir = .GlobalEnv)) {
      inlist    <- get(savedSettingsName, envir = .GlobalEnv)
      initnopts <- bifCheckNumSettings(initnopts, inlist)
      initpopts <- bifCheckPlotSettings(initpopts, inlist, state, parms)
    }

    # Read options from the command line
    adjustableopts <- c("lwd", "cex", "tcl.len", "bifsym", "biflblpos", "unstablelty")
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

    initpopts[[1]]$xlab  <- (c("Time", statenames))[initpopts[[1]]$xcol]
    initpopts[[1]]$ylab  <- ifelse(initpopts[[1]]$ycol  == 1, "State variables", statenames[initpopts[[1]]$ycol -1])
    initpopts[[1]]$y2lab <- ifelse(initpopts[[1]]$y2col == 1, "None",            statenames[initpopts[[1]]$y2col-1])
    initpopts[[2]]$xlab  <- parmsnames[initpopts[[2]]$xcol]
    initpopts[[2]]$ylab  <- ifelse(initpopts[[2]]$ycol  == 1, "State variables", statenames[initpopts[[2]]$ycol -1])
    initpopts[[2]]$y2lab <- ifelse(initpopts[[2]]$y2col == 1, "None",            statenames[initpopts[[2]]$y2col-1])
    initpopts[[3]]$xlab  <- parmsnames[initpopts[[3]]$xcol]
    initpopts[[3]]$ylab  <- parmsnames[initpopts[[3]]$ycol]

    # Read the curves from the environment
    initCurves <- list(Orbits = list(), BifurcationCurves = list(), BifurcationBounds = list(), TotalCurves = 0)
    if (resume && exists(savedCurvesName, envir = .GlobalEnv)) {
      inlist     <- get(savedCurvesName, envir = .GlobalEnv)
      initCurves <- bifCheckInputCurves(NULL, inlist, statenames, parmsnames)
    }

    ui <- bifUI(state, parms, initpopts, initnopts)

    # server <- function(input, output, session) {bifServer(input, output, session)}

    ############################################# BEGIN SERVER FUNCTION ####################################################
    server <- function(input, output, session) {

      # Hide the button to collapse the left sidebar
      shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")

      # Create the variable curveList as a reactive value, such that the plots will be updated when curveList changes
      curveList <- reactiveValues()
      curveList$Orbits <- initCurves$Orbits
      curveList$BifurcationCurves <- initCurves$BifurcationCurves
      curveList$BifurcationBounds <- initCurves$BifurcationBounds
      curveList$TotalCurves <- initCurves$TotalCurves
      curveListNames <- c('Orbits', 'BifurcationCurves', 'BifurcationBounds', 'TotalCurves')

      updateSpecialPointsList(session, initCurves, c(0, 0, 0))

      # Create the variable plotopts as a reactive value, such that the plots will be updated when plotopts changes
      plotopts <- reactiveValues()
      plotopts$Orbits <- initpopts[[1]]
      plotopts$BifurcationCurves <- initpopts[[2]]
      plotopts$BifurcationBounds <- initpopts[[3]]

      # Create the variable numopts as a reactive value, even though computations will only start after a button push
      numopts <- reactiveValues()
      lapply((1:length(initnopts)), function(i) {numopts[[(names(initnopts)[[i]])]] <- initnopts[[i]]})

      # Create the variable consoleLog as a reactive value, such that the console will be updated when text is added to it
      consoleLog <- reactiveVal()

      busyComputing <- reactiveVal(0)
      updatePlot <- reactiveVal(0)
      curveDirection <- reactiveVal(0)
      changeCurveMenu <- reactiveVal(0)

      # React to changes in curveList or plotopts, by replotting the current plot
      observe({
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        if (as.numeric(updatePlot()) == 0) return(NULL)

        output[[paste0("plot", curtab)]] <- renderPlot({
          if (curtab == 1) bifOrbitplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else if (curtab == 2) bif1parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else bif2parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
        },
        height = function() {setPlotHeight(session, input)},
        width = function() {setPlotWidth(session, input)})
        updatePlot(0)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Handle requests to save the plot as an image file
      output$saveplot <- downloadHandler(
        filename = function() {
          tempfile(pattern = "Rplot", tmpdir = getwd(), fileext = paste0(".", initpopts[[1]]$saveplotas))
        },
        content = function(file) {
          curtab <- as.numeric(isolate(input$plottab))
          curtabname <- curveListNames[curtab]
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
          if (curtab == 1) bifOrbitplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else if (curtab == 2) bif1parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else bif2parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          dev.off()
        })

      # React to changes in curveList by updating the special point selection menu
      observe({
        if (as.numeric(isolate(busyComputing())) == 1) return(NULL)

        if (as.numeric(changeCurveMenu()) == 0) return(NULL)
        else if (as.numeric(changeCurveMenu()) == 1)
          updateSpecialPointsList(session, reactiveValuesToList(curveList),
                                  c(as.numeric(isolate(input$selectpoint1)), as.numeric(isolate(input$selectpoint2)), as.numeric(isolate(input$selectpoint3))))
        else
          updateSpecialPointsList(session, reactiveValuesToList(curveList), c(0, 0, 0))

        # Updating the save and delete curve menu
        curtabname <- curveListNames[as.numeric(isolate(input$plottab))]
        updateCurveMenu(session, curveList[[curtabname]], 3)
        changeCurveMenu(0)
      })

      # React to changes in consoleLog, by rendering the new text
      observe({
        shinyjs::html(id = "console", html = HTML(gsub("\n", "<br>", consoleLog())))
        session$sendCustomMessage(type = "scrollCallback", 1)
      })

      # React to selection of an initial point
      observeEvent(c(input$selectpoint1, input$selectpoint2, input$selectpoint3), {
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        clist <- reactiveValuesToList(curveList)
        updateSelectedPoint(session, curtab, clist, as.numeric(input[[paste0('selectpoint', curtab)]]), statenames, parmsnames)
      })

      # Initialise a computation. Triggered by a change in the reactive value curveDirection()
      observe({
        if (as.numeric(curveDirection()) == 0) return(NULL)

        isolate({
          # Close the rightSidebar
          shinyjs::removeClass(id = "controlbar", class = "control-sidebar-open")
          # Collapse the State variables and Parameters stacks
          shinyjs::removeClass(selector = "li.treeview", class = "active")
          shinyjs::hide(selector = "ul.menu-open");
          updateSelectInput(session, "deletecurve1", selected = 0)
          updateSelectInput(session, "deletecurve2", selected = 0)
          updateSelectInput(session, "deletecurve3", selected = 0)

          curtab <- as.numeric(input$plottab)
          curtabname <- curveListNames[curtab]
          clist <- reactiveValuesToList(curveList)
          popts <- reactiveValuesToList(plotopts)
          pointid <- as.numeric(input[[paste0('selectpoint', curtab)]])

          if (curtab == 1) {
            initstate <- state
            initparms <- parms
            for (i in statenames) initstate[i] <- input[[paste0(i, "_1")]]
            for (i in parmsnames) initparms[i] <- input[[paste0(i, "_1")]]
            numopts$tstep <- as.numeric(curveDirection())*abs(numopts$tstep)

            newlist <- computeTimeseries(session, model, initstate, initparms, clist, pointid, numopts)
            if (!is.null(newlist)) lapply((1:length(newlist)),
                                          function(i) {curveList[[(curveListNames[[i]])]] <- newlist[[(curveListNames[[i]])]]})
            changeCurveMenu(1)
            curveDirection(0)
          } else {
            if (curtabname == 'BifurcationCurves') curvetype <- input$curvetype2
            else if (curtabname == 'BifurcationBounds') curvetype <- input$curvetype3

            # Get the starting point
            if (pointid > 0) {
              ind1 <- round(pointid/1000000)
              ind2 <- round((pointid-ind1*1000000)/1000)
              ind3 <- round(pointid-ind1*1000000-ind2*1000)
              cln1 <- curveListNames[ind1]
              ii <- ifelse((ind1 == 3), 2, 1) # 2 parameter bifurcation points have 2 columns before state, otherwise only 1

              initstate <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (ii + (1:length(state)))])
              initparms <- as.numeric(clist[[cln1]][[ind2]]$parameters)
              inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
              if (inittype != "TS")
                initparms[as.numeric(clist[[cln1]][[ind2]]$bifpars)] <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (1:ii)])
              inittanvec <- clist[[cln1]][[ind2]]$tanvec[ind3,]
              } else {
              initstate <- state
              initparms <- parms
              for (i in statenames) initstate[i] <- input[[paste0(i, "_", curtab)]]
              for (i in parmsnames) initparms[i] <- input[[paste0(i, "_", curtab)]]
              inittype <- "US"
              inittanvec <- NULL
            }
            names(initstate) <- names(state)
            names(initparms) <- names(parms)

            newlist <- initCurveContinuation(session, model, initstate, initparms, inittanvec, curtabname, clist,
                                             curvetype, inittype, popts[[curtabname]], numopts, as.numeric(input$reportlevel),
                                             as.numeric(curveDirection()))
            if (!is.null(newlist)) {
              lapply((1:length(newlist)), function(i) {curveList[[(curveListNames[[i]])]] <- newlist[[(curveListNames[[i]])]]})
              updatePlot(1)
              if (numopts$maxpoints > 1) {
                busyComputing(1)
                updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
                shinyjs::show("pausebtn")
                shinyjs::show("stopbtn")
              } else {
                changeCurveMenu(1)
                curveDirection(0)
              }
            } else {
              curveDirection(0)
            }
          }

          # Update the console log
          consoleLog(session$userData$alltext)
        })
      })

      # Compute the next batch batch of solution points along a 1- or 2-parameter bifurcaiton curve
      observe({
        if ((as.numeric(busyComputing()) != 1) || is.null(session$userData$curveData)) return(NULL)

        isolate({
          curveData <- session$userData$curveData
          nsol <- tryCatch(nextCurvePoints(isolate(round(as.numeric(numopts$replotfreq))), session$userData$curveData,
                                           plotopts[[curveData$tabname]], numopts, session = session),
                           warning = function(e) {
                             msg <- gsub(".*:", "Warning in nextCurvePoints:", e)
                             if (!is.null(session)) updateConsoleLog(session, msg)
                             else cat(msg)
                             session$userData$curveData <- NULL
                             return(NULL)
                           },
                           error = function(e) {
                             msg <- gsub(".*:", "Error in nextCurvePoints:", e)
                             if (!is.null(session)) updateConsoleLog(session, msg)
                             else cat(msg)
                             session$userData$curveData <- NULL
                             return(NULL)
                           })

          if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {

            newcurve <- curveList[[curveData$tabname]][[curveData$newcurvenr]]
            newcurve$points <- rbind(newcurve$points, nsol$points)
            newcurve$eigvals <- rbind(newcurve$eigvals, nsol$eigvals)
            newcurve$tanvec <- rbind(newcurve$tanvec, nsol$tanvec)
            newcurve$special.points <- rbind(newcurve$special.points, nsol$special.points)
            newcurve$special.tags <- rbind(newcurve$special.tags, nsol$special.tags)
            newcurve$special.eigvals <- rbind(newcurve$special.eigvals, nsol$special.eigvals)
            newcurve$special.tanvec <- rbind(newcurve$special.tanvec, nsol$special.tanvec)

            curveList[[curveData$tabname]][[curveData$newcurvenr]] <- newcurve
          }
        })

        # Update the console log
        consoleLog(session$userData$alltext)

        # Scheduling the computation of the subsequent batch of points: Invalidate this for later if computation has not ended
        if (!is.null(session$userData$curveData)) {
          updatePlot(1)
          invalidateLater(10, session)
        } else {
          changeCurveMenu(1)
          updatePlot(1)
          curveDirection(0)
          busyComputing(0)
          shinyjs::hide("pausebtn")
          shinyjs::hide("stopbtn")
        }
      })

      # Respond to the Pause button: observeEvent is non-reactive, it only reacts to invalidation of the specified event
      observeEvent(input$pausebtn, {
        busycomp <- as.numeric(isolate(busyComputing()))
        if (busycomp == 0) return(NULL)
        else if (busycomp == 1) {
          updateActionButton(session, "pausebtn", label = "Continue", icon = icon("forward"))
          busyComputing(-1)
          changeCurveMenu(1)
          updatePlot(1)
        }
        else {
          updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
          busyComputing(1)
        }
      })

      # Respond to the Stop button: observeEvent is non-reactive, it only reacts to invalidation of the specified event
      observeEvent(input$stopbtn, {
        if (as.numeric(isolate(busyComputing())) == 0) return(NULL)
        updateConsoleLog(session, "Computation interrupted by the user\n")

        curveData <- session$userData$curveData

        # Add the final points as special point
        newcurve <- curveList[[curveData$tabname]][[curveData$newcurvenr]]
        if (!is.null(newcurve) && !is.null(newcurve$points)) {
          if (curveData$curvetype == "LC") {
            statedim <- length(newcurve$initstate)
            freeparsdim <- length(newcurve$bifpars)

            vals <- lapply((1:statedim), function(i) {
              indxrange <- statedim*(1:(numopts$ninterval*numopts$glorder))
              yname <- names(newcurve$points[1, freeparsdim+i])
              y <- newcurve$points[nrow(newcurve$points), freeparsdim+i+indxrange]
              return(paste0("Min.", yname, "=", round(min(y), 3), ", Max.", yname, "=", round(max(y), 3)))
            })
            endPnt <- c("Type" = curveData$curvetype,
                        "Description" = paste0(names(newcurve$points[1, 1]), "=",
                                               round(newcurve$points[nrow(newcurve$points), 1], 3), " ",
                                               names(newcurve$points[1, ncol(newcurve$points)]), "=",
                                               round(newcurve$points[nrow(newcurve$points), ncol(newcurve$points)], 3), " ",
                                               paste(unlist(vals), collapse = ' ')))
          } else {
            endPnt <- c("Type" = curveData$curvetype,
                        "Description" = paste(unlist(lapply(1:length(newcurve$points[1,]),
                                                            function(i) {
                                                              paste0(names(newcurve$points[1,i]), "=",
                                                                     round(newcurve$points[nrow(newcurve$points), i], 3))
                                                            })),
                                              collapse=', '))
          }
          updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
          endPnt["Description"] <- paste0(sprintf("%04d: ", (curveData$pntnr-1)), endPnt["Description"])

          newcurve$special.points <- rbind(newcurve$special.points, newcurve$points[nrow(newcurve$points),])
          newcurve$special.eigvals <- rbind(newcurve$special.eigvals, newcurve$special.eigvals[nrow(newcurve$special.eigvals),])
          newcurve$special.tanvec <- rbind(newcurve$special.tanvec, newcurve$special.tanvec[nrow(newcurve$special.tanvec),])
          newcurve$special.tags <- rbind(newcurve$special.tags, c(endPnt))

          curveList[[curveData$tabname]][[curveData$newcurvenr]] <- newcurve
        }

        session$userData$curveData <- NULL
        changeCurveMenu(1)
        updatePlot(1)
        curveDirection(0)
        busyComputing(0)
        updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
        shinyjs::hide("pausebtn")
        shinyjs::hide("stopbtn")

        # Update the console log
        consoleLog(session$userData$alltext)
      }, ignoreInit = TRUE)

      # Respond to a change in plot tabs
      observeEvent(input$plottab, {
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        clist <- reactiveValuesToList(curveList)
        popts <- reactiveValuesToList(plotopts)
        plotopts[[curtabname]] <- updatePlotOptionEntries(session, curtab, popts[[curtabname]], statenames, parmsnames)
        updatePlot(1)
        updateCurveMenu(session, curveList[[curtabname]], 3)
      })

      # Apply newly entered plot or numerical settings
      observeEvent(c(input$plotoptsapply, input$numoptsapply), {
        # Close the rightSidebar
        # shinyjs::removeClass(id = "controlbar", class = "control-sidebar-open")

        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        popts <- reactiveValuesToList(plotopts)
        plotopts[[curtabname]] <- processPlotOptionsApply(session, input, curtab, popts[[curtabname]])

        plotopts[[curtabname]]$xlab <- ifelse(curtab == 1, (c("Time", statenames))[plotopts[[curtabname]]$xcol], parmsnames[plotopts[[curtabname]]$xcol])
        plotopts[[curtabname]]$ylab <- ifelse(curtab  < 3, (c("State variables", statenames))[plotopts[[curtabname]]$ycol], parmsnames[plotopts[[curtabname]]$ycol])
        plotopts[[curtabname]]$y2lab <- (c("None", statenames))[plotopts[[curtabname]]$y2col]

        processNumOptionsApply(session, input, curtab, numopts)
        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Initiate a forward computation
      observeEvent(c(input$computefwrd1,input$computefwrd2,input$computefwrd3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curveDirection(1)
      }, ignoreInit = TRUE)

      # Initiate a backward computation
      observeEvent(c(input$computebwrd1,input$computebwrd2,input$computebwrd3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curveDirection(-1)
      }, ignoreInit = TRUE)

      # Delete one or more curves
      observeEvent(c(input$deletebtn1, input$deletebtn2, input$deletebtn3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        curveList[[curtabname]] <- processDeleteCurve(session, curtab, curveList[[curtabname]], as.numeric(input[[paste0('deletecurve', curtab)]]))
        changeCurveMenu(-1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Save one or more curves
      observeEvent(c(input$savebtn1, input$savebtn2, input$savebtn3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        processSaveCurve(session, curtab, curveList[[curtabname]], as.numeric(input[[paste0('savecurve', curtab)]]),
                         make.names(input[[paste0('curvename', curtab)]], unique = TRUE))

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Append one or more curves to the current curve list
      observeEvent(c(input$appendbtn1, input$appendbtn2, input$appendbtn3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curvarname <- input[[paste0('loadcurve', curtab)]]
        # Check whether the name is a valid R variable name
        if (make.names(curvarname) != curvarname) return(NULL)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, curvarname, statenames, parmsnames, replace = FALSE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}

        updateCurveMenu(session, curveList[[curveListNames[curtab]]], 3)
        updateSpecialPointsList(session, reactiveValuesToList(curveList), c(0, 0, 0))

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Replace the current curve list with a stored list
      observeEvent(c(input$replacebtn1, input$replacebtn2, input$replacebtn3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curvarname <- input[[paste0('loadcurve', curtab)]]
        # Check whether the name is a valid R variable name
        if (make.names(curvarname) != curvarname) return(NULL)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, curvarname, statenames, parmsnames, replace = TRUE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}

        updateCurveMenu(session, curveList[[curveListNames[curtab]]], 3)
        updateSpecialPointsList(session, reactiveValuesToList(curveList), c(0, 0, 0))

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      # Show the manual
      observeEvent(input$helpClicked, {
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
                                                                               list(Orbits = curveList$Orbits, BifurcationCurves = curveList$BifurcationCurves,
                                                                                    BifurcationBounds = curveList$BifurcationBounds, TotalCurves = curveList$TotalCurves), 1L)
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
    ############################################## END SERVER FUNCTION #####################################################

    shinyApp(ui = ui, server = server, options = list(width = 910, height = 950))
  }
}

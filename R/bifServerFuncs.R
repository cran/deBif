updateConsoleLog <- function(session, addtext) {
  if ("alltext" %in% names(session$userData)) session$userData$alltext <- paste(session$userData$alltext, addtext, sep = "")
  else session$userData$alltext <- addtext
}

updateCurveMenu <- function(session, clist, ntabs) {
  # Update the save and delete curve menu
  lbls <- do.call("rbind", lapply(clist, "[[", "label"))
  ids <- c((0:length(clist)), -1)
  names(ids) <- c("None", lbls[,1], "All")[1:length(ids)]

  for (i in (1:ntabs)) {
    updateSelectInput(session, paste0("deletecurve", i), choices=ids, selected=0)
    updateSelectInput(session, paste0("savecurve", i), choices=ids, selected=0)
  }
}

updatePlotOptionEntries <- function(session, curtab, popts, snames, pnames) {
  # Update the plot options
  if (curtab == 1){
    updateSelectInput(session, "xcol",  label='Variable(s) on X-axis',
                      choices=c("Time" = 1, setNames((2:(length(snames)+1)), snames)), selected=popts$xcol)
    updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis',
                      choices=c("All" = 1, setNames((2:(length(snames)+1)), snames)), selected=popts$ycol)
    updateSelectInput(session, "y2col", label='Variable on 2nd Y-axis',
                      choices=c("None" = 1, setNames((2:(length(snames)+1)), snames)), selected=popts$y2col)
  } else if (curtab == 2){
    updateSelectInput(session, "xcol",  label='Bifurcation parameter',
                      choices=c(setNames((1:(length(pnames))), pnames)), selected=popts$xcol)
    updateSelectInput(session, "ycol",  label='Variable(s) on Y-axis',
                      choices=c("All" = 1, setNames((2:(length(snames)+1)), snames)), selected=popts$ycol)
    updateSelectInput(session, "y2col", label='Variable on 2nd Y-axis',
                      choices=c("None" = 1, setNames((2:(length(snames)+1)), snames)), selected=popts$y2col)
  } else {
    updateSelectInput(session, "xcol",  label='1st bifurcation parameter',
                      choices=c(setNames((1:(length(pnames))), pnames)), selected=popts$xcol)
    updateSelectInput(session, "ycol",  label='2nd bifurcation parameter',
                      choices=c(setNames((1:(length(pnames))), pnames)), selected=popts$ycol)
  }
  updateSelectInput(session,  "logx", selected=popts$logx)
  updateNumericInput(session, "xmin", value=popts$xmin)
  updateNumericInput(session, "xmax", value=popts$xmax)
  popts$xlab <- ifelse(curtab == 1, (c("Time", snames))[popts$xcol], pnames[popts$xcol])

  updateSelectInput(session,  "logy", selected=popts$logy)
  updateNumericInput(session, "ymin", value=popts$ymin)
  updateNumericInput(session, "ymax", value=popts$ymax)
  popts$ylab <- ifelse(curtab < 3, (c("State variables", snames))[popts$ycol], pnames[popts$ycol])

  updateSelectInput(session,  "logy2", selected=popts$logy2)
  updateNumericInput(session, "y2min", value=popts$y2min)
  updateNumericInput(session, "y2max", value=popts$y2max)
  popts$y2lab <- (c("None", snames))[popts$y2col]
  return(popts)
}

updateSelectedPoint <- function(session, curtab, clist, pointid, snames, pnames) {
  if (pointid > 0) {
    ind1 <- round(pointid/1000000)
    ind2 <- round((pointid-ind1*1000000)/1000)
    ind3 <- round(pointid-ind1*1000000-ind2*1000)
    cln1 <- (c('Orbits', 'BifurcationCurves', 'BifurcationBounds'))[ind1]

    lapply(snames,
           function(x) {
             updateNumericInput(session, paste0(x, "_", curtab), value=as.numeric(clist[[cln1]][[ind2]]$special.points[[ind3,x]]))
             })
    lapply(pnames,
           function(x) {
             cnames <- colnames(clist[[cln1]][[ind2]]$special.points)
             if (x %in% cnames)
               updateNumericInput(session, paste0(x, "_", curtab), value=as.numeric(clist[[cln1]][[ind2]]$special.points[[ind3,x]]))
             else
               updateNumericInput(session, paste0(x, "_", curtab), value=as.numeric(clist[[cln1]][[ind2]]$parameters[[x]]))
           })

    inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
    if ((curtab == 3) && (inittype %in% c("BP", "HP", "LP"))) {
      updateSelectInput(session, "curvetype3", selected = inittype)
    } else if (curtab == 2) {
      if (inittype == "HP") updateSelectInput(session, "curvetype2", selected = "LC")
      else updateSelectInput(session, "curvetype2", selected = "EQ")
    }
  }
}

updateSpecialPointsList <- function(session, clist, selected) {
  splist <- list()
  listlbls <- NULL
  ntabs <- length(selected)
  for (i in (1:ntabs)) {
    cln <- (c('Orbits', 'BifurcationCurves', 'BifurcationBounds'))[i]
    if (length(clist[[cln]]) > 0) {
      listlbls <- c(listlbls, unlist(lapply((1:length(clist[[cln]])), function(j){return(clist[[cln]][[j]]$label)})))
      splist <- c(splist, lapply((1:length(clist[[cln]])),
                                 function(j){
                                   lbls <- unlist(clist[[cln]][[j]]$special.tags[,"Description"], use.names=FALSE);
                                   ids <- ((i*1000000)+j*1000)+(1:length(lbls)); names(ids) <- lbls;
                                   return(ids)}))
    }
  }
  if (length(splist) > 0) {
    names(splist) <- listlbls
    for (i in (1:ntabs)) {
      updateSelectInput(session, paste0("selectpoint", i), choices=c(list("User specified" = 0), splist), selected=selected[i])
    }
  } else {
    for (i in (1:ntabs)) {
      updateSelectInput(session, paste0("selectpoint", i), choices=c(list("User specified" = 0), splist), selected=0)
    }
  }
}

processPlotOptionsApply <- function(session, input, curtab, popts) {
  popts$xcol <- as.numeric(input[["xcol"]])
  popts$logx <- as.numeric(input[["logx"]])
  popts$xmin <- ifelse(input[["logx"]] == 1, max(as.numeric(input[["xmin"]]), 1.0E-10), as.numeric(input[["xmin"]]))
  popts$xmax <- as.numeric(input[["xmax"]])
  popts$ycol <- as.numeric(input[["ycol"]])
  popts$logy <- as.numeric(input[["logy"]])
  popts$ymin <- ifelse(input[["logy"]] == 1, max(as.numeric(input[["ymin"]]), 1.0E-10), as.numeric(input[["ymin"]]))
  popts$ymax <- as.numeric(input[["ymax"]])

  if (popts$xmax < 1.0001*popts$xmin) {
    updateConsoleLog(session, "Maximum of x-axis not significantly different from its minimum\n")
    popts$xmax <- 1.0001*popts$xmin
  }
  if (popts$ymax < 1.0001*popts$ymin) {
    updateConsoleLog(session, "Maximum of y-axis not significantly different from its minimum\n")
    popts$ymax <- 1.0001*popts$ymin
  }

  popts$plot3d <- as.numeric(0)
  if ((curtab < 3) && (popts$ycol > 1)) {
    popts$y2col <- as.numeric(input[["y2col"]])
    popts$logy2 <- as.numeric(input[["logy2"]])
    popts$y2min <- ifelse(input[["logy2"]] == 1, max(as.numeric(input[["y2min"]]), 1.0E-10), as.numeric(input[["y2min"]]))
    popts$y2max <- as.numeric(input[["y2max"]])
    popts$plot3d <- ifelse((popts$y2col > 1), as.numeric(input[["plot3d"]]), 0)
    popts$theta <- as.numeric(input[["theta"]])
    popts$theta <- max(popts$theta, -90)
    popts$theta <- min(popts$theta, 90)

    if (popts$y2max < 1.0001*popts$y2min) {
      updateConsoleLog(session, "Maximum of y-axis not significantly different from its minimum\n")
      popts$y2max <- 1.0001*popts$y2min
    }
  }
  return(popts)
}

processNumOptionsApply <- function(session, input, curtab, nopts) {
  text2numeric <- function(oldval, newval){
    newval2 <- gsub("[^0-9.E+-]*", "", newval)
    if (!is.na(suppressWarnings(as.numeric(newval2)))) return (as.numeric(newval2))
    else return(as.numeric(oldval))
  }

  if (curtab == 1) {
    nopts$tmax <- as.numeric(input[["tmax"]])
    nopts$tstep <- abs(as.numeric(input[["tstep"]]))
    nopts$odemethod <- input[["method"]]
    if ("ssgrid" %in% names(input)) nopts$ssgrid <- min(max(input[["ssgrid"]], 1), 50)
    if ("pgrid" %in% names(input)) nopts$pgrid <- min(max(input[["pgrid"]], 1), 10)
  }
  else {
    nopts$rhstol <- max(text2numeric(nopts$rhstol, input[["rhstol"]]), 1.0E-10)
    updateTextInput(session, "rhstol", value=sprintf("%.1E", nopts$rhstol))
    nopts$dytol <- max(text2numeric(nopts$dytol, input[["dytol"]]), 1.0E-10)
    updateTextInput(session, "dytol", value=sprintf("%.1E", nopts$dytol))
    nopts$iszero <- max(text2numeric(nopts$iszero, input[["iszero"]]), 1.0E-10)
    updateTextInput(session, "iszero", value=sprintf("%.1E", nopts$iszero))
    nopts$neartol <- max(text2numeric(nopts$neartol, input[["neartol"]]), 1.0E-10)
    updateTextInput(session, "neartol", value=sprintf("%.1E", nopts$neartol))
    nopts$jacdif <- max(text2numeric(nopts$jacdif, input[["jacdif"]]), 1.0E-10)
    updateTextInput(session, "jacdif", value=sprintf("%.1E", nopts$jacdif))

    nopts$minstepsize <- max(as.numeric(input[["minstepsize"]]), 1.0E-10)
    nopts$maxstepsize <- max(as.numeric(input[["maxstepsize"]]), nopts$minstepsize)
    nopts$maxiter <- max(as.numeric(input[["maxiter"]]), 1)
    nopts$maxpoints <- max(as.numeric(input[["maxpoints"]]), 1)
    nopts$replotfreq <- max(round(as.numeric(input[["replotfreq"]])), 1)

    nopts$ninterval <- min(max(round(as.numeric(input[["ninterval"]])), 1), 40)
    updateTextInput(session, "ninterval", value=nopts$ninterval)
    nopts$glorder <- min(max(round(as.numeric(input[["glorder"]])), 2), 7)
    updateTextInput(session, "glorder", value=nopts$glorder)
    nopts$lcampl <- max(text2numeric(nopts$lcampl, input[["lcampl"]]), 1.0E-7)
    updateTextInput(session, "lcampl", value=sprintf("%.1E", nopts$lcampl))
  }
  return(nopts)
}

processDeleteCurve <- function(session, curtab, clist, deletenr) {
  totalcurves <- as.numeric(length((clist)))

  if ((totalcurves > 0) && ((deletenr > 0) || (deletenr == -1)) && (deletenr < (totalcurves + 1))) {

    if (deletenr == -1) {
      msg <- paste0("All curves deleted\n")
      clist <- list()
    }
    else {
      msg <- paste0("Curve '", clist[[deletenr]]$label, "' deleted\n")
      clist[[deletenr]] <- NULL
    }
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
  }

  return(clist)
}

processSaveCurve <- function(session = NULL, curtab, clist, savenr, varname) {
  totalcurves <- as.numeric(length((clist)))

  if ((totalcurves > 0) && ((savenr > 0) || (savenr == -1)) && (savenr < (totalcurves + 1))) {
    if (exists(varname, envir = .GlobalEnv)) {
      rm(list = varname, envir = .GlobalEnv)
    }
    # if (savenr == -1) assign(varname, clist, envir = .GlobalEnv)
    # else assign(varname, clist[[savenr]], envir = .GlobalEnv)
    # global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
    if (savenr == -1) {
      msg <- paste0("All curves saved to '", varname, "'\n")
      (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(varname, clist, 1L)
    }
    else {
      msg <- paste0("Curve '", clist[[savenr]]$label, "' saved to '", varname, "'\n")
      (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(varname, clist[[savenr]], 1L)
    }
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
  }
}

processLoadCurve <- function(session = NULL, clist, varname, snames, pnames, replace = FALSE) {
  inlist <- NULL
  if (exists(varname, envir = .GlobalEnv)) {
    inlist <- get(varname, envir = .GlobalEnv)
  }

  if (!any(c("Orbits", "BifurcationCurves", "BifurcationBounds") %in% names(inlist))) {
    if ("type" %in% names(inlist)) {
      ctype <- inlist$type
      inlist <- list(inlist)
    }
    else if ("type" %in% names(inlist[[1]])) ctype <- inlist[[1]]$type
    else ctype <- NULL
    if (!is.null(ctype)) {
      if (ctype == "TS") lbl <- "Orbits"
      else if ((ctype == "EQ") || (ctype == "LC")) lbl <- "BifurcationCurves"
      else lbl <- "BifurcationBounds"

      inlist <- list(inlist)
      names(inlist) <- lbl
    }
  }

  if (is.null(inlist) || (!is.list(inlist)) || (length(inlist) == 0)) {
    msg <- paste0("Curves not loaded: variable '", varname, "' not a valid list in the global environment\n")
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
    return(NULL)
  }

  nlist <- bifCheckInputCurves(NULL, inlist, snames, pnames)

  if ((length(nlist$Orbits) == 0) && (length(nlist$BifurcationCurves) == 0) && (length(nlist$BifurcationBounds) == 0)) {
    msg <- paste0("Curves not loaded: variable '", varname, "' does not contain any valid curve data\n")
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
    return(NULL)
  }
  outlist <- clist
  if (length(nlist$Orbits) > 0) {
    if (replace) outlist$Orbits <- nlist$Orbits
    else outlist$Orbits <- c(clist$Orbits, nlist$Orbits)
  }
  if (length(nlist$BifurcationCurves) > 0) {
    if (replace) outlist$BifurcationCurves <- nlist$BifurcationCurves
    else outlist$BifurcationCurves <- c(clist$BifurcationCurves, nlist$BifurcationCurves)
  }
  if (length(nlist$BifurcationBounds) > 0) {
    if (replace) outlist$BifurcationBounds <- nlist$BifurcationBounds
    else outlist$BifurcationBounds <- c(clist$BifurcationBounds, nlist$BifurcationBounds)
  }
  if (replace) outlist$TotalCurves <- nlist$TotalCurves
  else outlist$TotalCurves <- clist$TotalCurves + nlist$TotalCurves

  msg <- paste0(nlist$TotalCurves, " curves loaded from variable '", varname, "'\n")
  if (!is.null(session)) updateConsoleLog(session, msg)
  else cat(msg)

  return(outlist)
}

setPlotHeight <- function(session, input) {
  curtab <- as.numeric(input$plottab)
  if (is.null(input$dimension))
    return(0.75*session$clientData[[paste0("output_plot", curtab, "_width")]])
  else return(min((as.numeric(input$dimension[2]) - 110), 0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]))
}

setPlotWidth <- function(session, input) {
  curtab <- as.numeric(input$plottab)
  return(0.99*session$clientData[[paste0("output_plot", curtab, "_width")]])
}

bifCheckNumSettings <- function(oldopts, inlist) {
  initnopts <- oldopts
  if (is.list(inlist) && ("numopts" %in% names(inlist))) {
    nlist <- inlist$numopts
    if (("odemethod" %in% names(nlist)) && (nlist$odemethod %in% c("lsoda", "ode23", "ode45", "rk4")))
      initnopts$odemethod <- nlist$odemethod
    if (("tmax" %in% names(nlist)) && is.numeric(nlist$tmax) && (nlist$tmax > 0))
      initnopts$tmax <- nlist$tmax
    if (("tstep" %in% names(nlist)) && is.numeric(nlist$tstep) && (nlist$tstep > 0))
      initnopts$tstep <- nlist$tstep
    if (("rhstol" %in% names(nlist)) && is.numeric(nlist$rhstol) && (nlist$rhstol > 0) && (nlist$rhstol < 1) )
      initnopts$rhstol <- nlist$rhstol
    if (("dytol" %in% names(nlist)) && is.numeric(nlist$dytol) && (nlist$dytol > 0) && (nlist$dytol < 1) )
      initnopts$dytol <- nlist$dytol
    if (("neartol" %in% names(nlist)) && is.numeric(nlist$neartol) && (nlist$neartol > 0) && (nlist$neartol < 1) )
      initnopts$neartol <- nlist$neartol
    if (("iszero" %in% names(nlist)) && is.numeric(nlist$iszero) && (nlist$iszero > 0) && (nlist$iszero < 1) )
      initnopts$iszero <- nlist$iszero
    if (("jacdif" %in% names(nlist)) && is.numeric(nlist$jacdif) && (nlist$jacdif > 0) && (nlist$jacdif < 1) )
      initnopts$jacdif <- nlist$jacdif
    if (("minstepsize" %in% names(nlist)) && is.numeric(nlist$minstepsize)) {
      initnopts$minstepsize <- max(nlist$minstepsize, 1.0E-10)
    }
    if (("maxstepsize" %in% names(nlist)) && is.numeric(nlist$maxstepsize))
      initnopts$maxstepsize <- max(nlist$maxstepsize, initnopts$minstepsize)
    else
      initnopts$maxstepsize <- max(initnopts$maxstepsize, initnopts$minstepsize)

    if (("initstepsize" %in% names(nlist)) && is.numeric(nlist$initstepsize)) {
      initnopts$initstepsize <- max(nlist$initstepsize, initnopts$minstepsize)
      initnopts$initstepsize <- min(initnopts$initstepsize, initnopts$maxstepsize)
    }
    else {
      initnopts$initstepsize <- max(initnopts$minstepsize, initnopts$maxstepsize / 10.0)
    }

    if (("maxiter" %in% names(nlist)) && is.numeric(nlist$maxiter) && (nlist$maxiter > 0))
      initnopts$maxiter <- as.integer(nlist$maxiter)
    if (("maxpoints" %in% names(nlist)) && is.numeric(nlist$maxpoints) && (nlist$maxpoints > 0))
      initnopts$maxpoints <- as.integer(nlist$maxpoints)
    if (("replotfreq" %in% names(nlist)) && is.numeric(nlist$replotfreq) && (nlist$replotfreq >= 1) && (nlist$replotfreq <= 10000) )
      initnopts$replotfreq <- as.integer(nlist$replotfreq)
    if (("ninterval" %in% names(nlist)) && is.numeric(nlist$ninterval) && (nlist$ninterval >= 1) && (nlist$ninterval <= 40) )
      initnopts$ninterval <- as.integer(nlist$ninterval)
    if (("glorder" %in% names(nlist)) && is.numeric(nlist$glorder) && (nlist$glorder >= 2) && (nlist$glorder <= 7) )
      initnopts$glorder <- as.integer(nlist$glorder)
    if (("lcampl" %in% names(nlist)) && is.numeric(nlist$lcampl) && (nlist$lcampl > 0))
      initnopts$lcampl <- nlist$lcampl
  }
  return(initnopts)
}


bifCheckPlotSettings <- function(oldopts, inlist, state, parms) {
  initpopts <- oldopts
  if (is.list(inlist) && ("plotopts" %in% names(inlist)) &&
      all(names(inlist$plotopts) %in% c("Orbits", "BifurcationCurves", "BifurcationBounds"))) {
    plist <- inlist$plotopts
    for (i in 1:length(initpopts)) {
      j <- (1:3)[names(plist) == (c("Orbits", "BifurcationCurves", "BifurcationBounds"))[i]]
      if ((i == 1) && ("xcol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$xcol) && (plist[[j]]$xcol %in% (1:(1+length(state)))))
        initpopts[[i]]$xcol <- plist[[j]]$xcol
      else if (("xcol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$xcol) && (plist[[j]]$xcol %in% (1:length(parms))))
        initpopts[[i]]$xcol <- plist[[j]]$xcol
      if ((i == 3) && ("ycol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$ycol) && (plist[[j]]$ycol %in% (1:length(parms))))
        initpopts[[i]]$ycol <- plist[[j]]$ycol
      else if (("ycol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$ycol) && (plist[[j]]$ycol %in% (1:(1+length(state)))))
        initpopts[[i]]$ycol <- plist[[j]]$ycol
      if ((i == 3) && ("y2col" %in% names(plist[[j]])) && is.numeric(plist[[j]]$y2col) && (plist[[j]]$y2col %in% (1:length(parms))))
        initpopts[[i]]$y2col <- plist[[j]]$y2col
      else if (("y2col" %in% names(plist[[j]])) && is.numeric(plist[[j]]$y2col) && (plist[[j]]$y2col %in% (1:(1+length(state)))))
        initpopts[[i]]$y2col <- plist[[j]]$y2col

      for (ax in c("x", "y", "y2")) {
        lbl <- paste0("log", ax)
        if ((lbl %in% names(plist[[j]])) && is.numeric(plist[[j]][[lbl]]) && (plist[[j]][[lbl]] %in% (0:1)))
          initpopts[[i]][[lbl]] <- plist[[j]][[lbl]]
        lbl <- paste0(ax, "min")
        lbl2 <- paste0(ax, "max")
        if ((lbl %in% names(plist[[j]])) && (lbl2 %in% names(plist[[j]]))
            && is.numeric(plist[[j]][[lbl]]) && is.numeric(plist[[j]][[lbl2]]) && (plist[[j]][[lbl]] < plist[[j]][[lbl2]])) {
          initpopts[[i]][[lbl]]  <- plist[[j]][[lbl]]
          initpopts[[i]][[lbl2]] <- plist[[j]][[lbl2]]
        }
      }
      if (("lwd" %in% names(plist[[j]])) && is.numeric(plist[[j]]$lwd) && (plist[[j]]$lwd >= 0))
        initpopts[[i]]$lwd <- plist[[j]]$lwd
      if (("cex" %in% names(plist[[j]])) && is.numeric(plist[[j]]$cex) && (plist[[j]]$cex >  0))
        initpopts[[i]]$cex <- plist[[j]]$cex
      if (("tcl.len" %in% names(plist[[j]])) && is.numeric(plist[[j]]$tcl.len) && (plist[[j]]$tcl.len >  0))
        initpopts[[i]]$tcl.len <- plist[[j]]$tcl.len
      if (("theta" %in% names(plist[[j]])) && is.numeric(plist[[j]]$theta) && (plist[[j]]$theta >=  -90) && (plist[[j]]$theta <=  90))
        initpopts[[i]]$theta <- plist[[j]]$theta
      if (("plot3d" %in% names(plist[[j]])) && is.numeric(plist[[j]]$plot3d) && (plist[[j]]$plot3d %in% (0:1)))
        initpopts[[i]]$plot3d <- plist[[j]]$plot3d
      if (("bifsym" %in% names(plist[[j]])) && is.numeric(plist[[j]]$bifsym) && (plist[[j]]$bifsym %in% (0:25)))
        initpopts[[i]]$bifsym <- plist[[j]]$bifsym
      if (("biflblpos" %in% names(plist[[j]])) && is.numeric(plist[[j]]$biflblpos) && (plist[[j]]$biflblpos %in% (1:4)))
        initpopts[[i]]$biflblpos <- plist[[j]]$biflblpos
      if (("unstablelty" %in% names(plist[[j]])) && is.numeric(plist[[j]]$unstablelty) && (plist[[j]]$unstablelty %in% (1:6)))
        initpopts[[i]]$unstablelty <- plist[[j]]$unstablelty
    }
  }
  return(initpopts)
}


bifCheckInputCurves <- function(oldcurves, inlist, snames, pnames) {
  if (is.null(oldcurves)) clist <- list(Orbits = list(), BifurcationCurves = list(), BifurcationBounds = list(), TotalCurves = 0)
  else clist <- oldcurves

  if (!is.null(inlist) && is.list(inlist)) {
    if (("Orbits" %in% names(inlist)) && (length(inlist$Orbits) > 0)) {
      nlist <- inlist$Orbits
      for (i in (1:length(nlist))) {
        if (all(c("label", "type", "initstate", "parameters", "points", "special.points", "special.tags") %in% names(nlist[[i]]))
            && (ncol(nlist[[i]]$points) == (length(snames)+1))
            && all(colnames(nlist[[i]]$points) == c("Time", snames))
            && (length(names(nlist[[i]]$parameters)) == length(pnames))
            && all(names(nlist[[i]]$parameters) == pnames)) {
          clist$Orbits[[length((clist$Orbits))+1]] <- nlist[[i]]
          clist$TotalCurves <- clist$TotalCurves + 1
        }
      }
    }
    if (("BifurcationCurves" %in% names(inlist))  && (length(inlist$BifurcationCurves) > 0)) {
      nlist <- inlist$BifurcationCurves
      for (i in (1:length(nlist))) {
        if (nlist[[i]]$type == "EQ") {
          if (all(c("label", "type", "initstate", "parameters", "bifpars", "points", "eigvals",
                    "tanvec", "special.points", "special.eigvals", "special.tanvec", "special.tags") %in% names(nlist[[i]]))
              && (ncol(nlist[[i]]$points) == (length(snames)+1))
              && ((colnames(nlist[[i]]$points))[1] %in% pnames)
              && all(colnames(nlist[[i]]$points)[2:(length(snames)+1)] == snames)
              && (length(names(nlist[[i]]$parameters)) == length(pnames))
              && all(names(nlist[[i]]$parameters) == pnames)) {
            clist$BifurcationCurves[[length((clist$BifurcationCurves))+1]] <- nlist[[i]]
            clist$TotalCurves <- clist$TotalCurves + 1
          }
        }
        if (nlist[[i]]$type == "LC") {
          if (all(c("label", "type", "initstate", "parameters", "bifpars", "points",
                    "tanvec", "special.points", "special.tanvec", "special.tags") %in% names(nlist[[i]]))
              && ((colnames(nlist[[i]]$points))[1] %in% pnames)
              && all(colnames(nlist[[i]]$points)[2:(length(snames)+1)] == snames)
              && (length(names(nlist[[i]]$parameters)) == length(pnames))
              && all(names(nlist[[i]]$parameters) == pnames)) {
            clist$BifurcationCurves[[length((clist$BifurcationCurves))+1]] <- nlist[[i]]
            clist$TotalCurves <- clist$TotalCurves + 1
          }
        }
      }
    }
    if (("BifurcationBounds" %in% names(inlist))  && (length(inlist$BifurcationBounds) > 0)) {
      nlist <- inlist$BifurcationBounds
      for (i in (1:length(nlist))) {
        if (all(c("label", "type", "initstate", "parameters", "bifpars", "points", "eigvals",
                  "tanvec", "special.points", "special.eigvals", "special.tanvec", "special.tags") %in% names(nlist[[i]]))
            && (ncol(nlist[[i]]$points) == (length(snames)+2))
            && all((colnames(nlist[[i]]$points))[1:2] %in% pnames)
            && all(colnames(nlist[[i]]$points)[3:(length(snames)+2)] == snames)
            && (length(names(nlist[[i]]$parameters)) == length(pnames))
            && all(names(nlist[[i]]$parameters) == pnames))
        {
          clist$BifurcationBounds[[length((clist$BifurcationBounds))+1]] <- nlist[[i]]
          clist$TotalCurves <- clist$TotalCurves + 1
        }
      }
    }
  }
  return(clist)
}


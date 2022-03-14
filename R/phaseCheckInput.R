phaseCheckNumSettings <- function(oldopts, inlist) {
  initnopts <- oldopts
  if (is.list(inlist) && ("numopts" %in% names(inlist))) {
    nlist <- inlist$numopts
    if (("odemethod" %in% names(nlist)) && (nlist$odemethod %in% c("lsoda", "ode23", "ode45", "rk4")))
      initnopts$odemethod <- nlist$odemethod
    if (("tmax" %in% names(nlist)) && is.numeric(nlist$tmax) && (nlist$tmax > 0))
      initnopts$tmax <- nlist$tmax
    if (("tstep" %in% names(nlist)) && is.numeric(nlist$tstep) && (nlist$tstep > 0))
      initnopts$tstep <- nlist$tstep
    if (("ssgrid" %in% names(nlist)) && is.numeric(nlist$ssgrid))
      initnopts$ssgrid <- min(max(nlist$ssgrid, 1), 50)
    if (("pgrid" %in% names(nlist)) && is.numeric(nlist$pgrid))
      initnopts$pgrid <- min(max(nlist$pgrid, 1), 10)
  }
  return(initnopts)
}

phaseCheckPlotSettings <- function(oldopts, inlist, state, parms) {
  initpopts <- oldopts
  if (is.list(inlist) && ("plotopts" %in% names(inlist)) &&
      all(names(inlist$plotopts) %in% c("Orbits", "PhasePlane"))) {
    plist <- inlist$plotopts
    for (i in 1:length(initpopts)) {
      j <- (1:2)[names(plist) == (c("Orbits", "PhasePlane"))[i]]
      if ((i > 1) && ("xcol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$xcol) && (plist[[j]]$xcol %in% (1:length(state))))
        initpopts[[i]]$xcol <- plist[[j]]$xcol

      if ((i == 1) && ("ycol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$ycol) && (plist[[j]]$ycol %in% (1:(1+length(state)))))
        initpopts[[i]]$ycol <- plist[[j]]$ycol
      else if (("ycol" %in% names(plist[[j]])) && is.numeric(plist[[j]]$ycol) && (plist[[j]]$ycol %in% (1:length(state))))
        initpopts[[i]]$ycol <- plist[[j]]$ycol

      if ((i == 1) && ("y2col" %in% names(plist[[j]])) && is.numeric(plist[[j]]$y2col) && (plist[[j]]$y2col %in% (1:length(state))))
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
    }
  }
  return(initpopts)
}

phaseCheckInputCurves <- function(oldcurves, inlist, snames, pnames) {
  if (is.null(oldcurves)) clist <- list(Orbits = list(), TotalCurves = 0)
  else clist <- oldcurves

  if (!is.null(inlist) && is.list(inlist)) {
    if (("Orbits" %in% names(inlist)) && (length(inlist$Orbits) > 0)) {
      nlist <- inlist$Orbits
      for (i in (1:length(nlist))) {
        if (all(c("label", "type", "initstate", "parameters", "points", "special.points", "special.tags") %in% names(nlist[[i]]))
            && (ncol(nlist[[i]]$points) == (length(snames)+1))
            && all(colnames(nlist[[i]]$points) == c("Time", snames))) {
          clist$Orbits[[length((clist$Orbits))+1]] <- nlist[[i]]
          clist$TotalCurves <- clist$TotalCurves + 1
        }
      }
    }
  }
  return(clist)
}

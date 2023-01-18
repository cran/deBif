phasePlot2D <- function(curtab, odes, state, parms, plotopts, numopts, zlst = NULL) {

  if (is.null(zlst)) {
    derivs <- nullclines(odes, state, parms, plotopts, numopts)
  } else {
    derivs <- zlst
  }

  xcol <- as.numeric(plotopts$xcol)
  ycol <- as.numeric(plotopts$ycol)
  xmin <- as.numeric(plotopts$xmin)
  xmax <- as.numeric(plotopts$xmax)
  ymin <- as.numeric(plotopts$ymin)
  ymax <- as.numeric(plotopts$ymax)
  logx <- (as.numeric(plotopts$logx) == 1)
  logy <- (as.numeric(plotopts$logy) == 1)
  xlab <- plotopts$xlab
  ylab <- plotopts$ylab
  grid <- numopts$pgrid
  npixels <- numopts$npixels

  # Make a phase plane with nullclines and/or phase portrait
  ishows <- c(xcol, ycol)
  lvec <- 50                         # length of vector

  logxy <- ""
  if (logx) {
    logxy <- paste0(logxy, "x")
    xmin <- max(xmin, 1.0E-10)
  }
  if (logy) {
    logxy <- paste0(logxy, "y")
    ymin <- max(ymin, 1.0E-10)
  }

  par(cex = as.numeric(plotopts["cex"]), mar = plotopts$plotmar)
  plot(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, log=logxy,
       cex.lab=as.numeric(plotopts["cex.lab"]), cex.axis=as.numeric(plotopts["cex.axis"]))

  for (i in (1:length(ishows)))
    {
      zc <- derivs[[i + 2]]                            # [[1]] = xc, [[2]] = yc, [[3]] = dxc, [[4]] = dyc
      goodrows <- !(rowSums(is.nan(zc)) == npixels)    # Ignore rows with only NaN
      goodcols <- !(colSums(is.nan(zc)) == npixels)    # Ignore columns with only NaN
      if (any(!is.nan(zc)))
        contour(derivs$xc[goodrows], derivs$yc[goodcols], zc[goodrows,goodcols],
                levels=0, drawlabels=FALSE, add=TRUE, col=plotopts$colors[ishows[i]], lwd=plotopts["lwd"])
    }

  if ((curtab == 4) || (curtab == 6)) {
    if (logx) {dx <- (log10(xmax)-log10(xmin))/grid; vx <- 1+3.32*grid*dx/lvec}
    else {dx <- (xmax-xmin)/grid; vx = grid*dx/lvec}
    if (logy) {dy <- (log10(ymax)-log10(ymin))/grid; vy <- 1+3.32*grid*dy/lvec}
    else {dy <- (ymax-ymin)/grid; vy = grid*dy/lvec}

    times <- seq(0, numopts$tmax, by=numopts$tstep)
    for (i in seq(1,grid)) {
      if (logx) state[xcol] <- 10^((i-1)*dx + dx/2 + log10(xmin))
      else state[xcol] <- (i-1)*dx + dx/2 + xmin
      for (j in seq(1,grid,1)) {
        if (logy) state[ycol] <- 10^((j-1)*dy + dy/2 + log10(ymin))
        else state[ycol] <- (j-1)*dy + dy/2 + ymin
        # Draw vectorfield
        if (curtab == 4) {
          dt <- sign(unlist(odes(0, state, parms)))
          if (logx) lines(c(state[xcol], state[xcol]*vx^dt[xcol]), c(state[ycol], state[ycol]))
          else lines(c(state[xcol], state[xcol]+vx*dt[xcol]), c(state[ycol], state[ycol]))
          if (logy) lines(c(state[xcol], state[xcol]), c(state[ycol], state[ycol]*vy^dt[ycol]))
          else lines(c(state[xcol], state[xcol]), c(state[ycol], state[ycol]+vy*dt[ycol]))
        }
        # Draw phase portrait
        if (curtab == 6) {
          points(state[xcol], state[ycol], pch=as.numeric(plotopts["pch"]))
          nsol <- as.data.frame(do.call('ode', c(list(times=times, func=odes, y=state, parms=parms), method=numopts$method)))
          lines(cbind(nsol[xcol+1], nsol[ycol+1]), col="black")
        }
      }
    }
  }
  legend("topright", legend=paste0('d', names(state)[ishows], '/dt'), col=plotopts$colors[ishows], lwd=plotopts["lwd"], cex=as.numeric(plotopts["cex.legend"]), bg = "white")

  return(derivs)
}

phasePlot2D <- function(curtab, odes, state, parms, plotopts, numopts) {
  # Save plot options to restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

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
  eps  <- numopts$eps
  npixels <- numopts$npixels

  # Make a phase plane with nullclines and/or phase portrait
  ishows <- c(xcol, ycol)
  nvar <- length(state)
  lvec <- 50                         # length of vector

  logxy <- ""
  if (logx) {
    xc <- 10^seq(log10(xmin), log10(xmax), length.out=npixels)
    logxy <- paste0(logxy, "x")
    xmin <- max(xmin, 1.0E-10)
  } else xc <- seq(xmin+eps, xmax, length.out=npixels)
  if (logy) {
    yc <- 10^seq(log10(ymin), log10(ymax), length.out=npixels)
    logxy <- paste0(logxy, "y")
    ymin <- max(ymin, 1.0E-10)
  } else yc <- seq(ymin+eps, ymax, length.out=npixels)

  par(cex = as.numeric(plotopts["cex"]), mar = plotopts$plotmar)
  plot(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, log=logxy,
       cex.lab=as.numeric(plotopts["cex.lab"]), cex.axis=as.numeric(plotopts["cex.axis"]))

  vstate <- as.list(state)
  npixels2 <- npixels^2
  #vstate<-lapply(vstate,rep,vstate,npixels2);vstate[[xcol]]<-0;vstate[[ycol]]<-0
  for (j in seq(1, nvar)) if (j!=xcol & j!=ycol) vstate[[j]]<-rep(vstate[[j]],npixels2);
  FUN <- function(xc, yc, i){  # wrapper around model()
    vstate[[xcol]] <- xc; vstate[[ycol]] <- yc
    odes(0, vstate, parms)[[1]][seq((i-1)*npixels2+1, i*npixels2)]
  }
  for (i in ishows)
    contour(xc, yc, outer(xc, yc, FUN, i), levels=0, drawlabels=FALSE, add=TRUE, col=plotopts$colors[i], lwd=plotopts["lwd"])

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
}

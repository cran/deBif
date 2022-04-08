phasePlot1D <- function(curtab, odes, state, parms, plotopts, numopts) {

  x <- 1
  xmin <- as.numeric(plotopts$xmin)
  xmax <- as.numeric(plotopts$xmax)
  ymin <- as.numeric(plotopts$ymin)
  ymax <- as.numeric(plotopts$ymax)
  logx <- (as.numeric(plotopts$logx) == 1)
  logy <- (as.numeric(plotopts$logy) == 1)
  eps  <- numopts$eps
  npixels <- numopts$npixels

  logxy <- ""
  if (logx) {
    xc <- 10^seq(log10(xmin), log10(xmax), length.out=npixels)
    logxy <- paste0(logxy, "x")
    xmin <- max(xmin, 1.0E-10)
  } else xc <- seq(xmin, xmax, length.out=npixels)
  if (logy) {
    logxy <- paste0(logxy, "y")
    ymin <- max(ymin, 1.0E-10)
  }
  xlab <- names(state)[x];
  ylab <- paste0("d", xlab, "/dt")
  par(cex = as.numeric(plotopts["cex"]), mar = plotopts$plotmar)
  do.call('plot', c(list(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, log=logxy,
                         cex.lab=as.numeric(plotopts["cex.lab"]), cex.axis=as.numeric(plotopts["cex.axis"]))))
  lines(c(xmin,xmax), c(0,0), col="black", lwd=1, lty=2)
  # legend("topright",legend=ylab, col=plotopts$colors[1], lty=1, lwd=plotopts["lwd"], cex=as.numeric(plotopts["cex.legend"]), bg = "white")

  dxdt <- as.numeric(lapply(xc, function(i) {state[1] <- i; odes(0, state, parms)[[1]][[1]]}))

  lines(xc, dxdt, lwd=plotopts["lwd"], col=plotopts$colors[1])
}


bif2parplot <- function(session = NULL, curvelist = NULL, popts) {
  # Save plot options to restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Plot the bifurcation curves

  par(cex = popts$cex, mar = as.numeric(c(1.6*(popts$cex.lab+popts$cex.axis), 1.6*(popts$cex.lab+popts$cex.axis), 2, 2)))
  logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
  plot(NULL, type='n', xlab="", ylab="", xaxs = "i", yaxs = "i",
       xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy, cex.axis=popts$cex.axis)
  title(xlab=popts$xlab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))
  title(ylab=popts$ylab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))

  if (!is.null(curvelist) && (length(curvelist) > 0)) {
    lapply((1:length(curvelist)), function(i) {
      cnames <- colnames(curvelist[[i]]$points)
      if (cnames[1] == popts$xlab) clx <- 1
      else {
        if (cnames[2] == popts$xlab) clx <- 2
        else {
          msg <- paste0("Curve plotting skipped: parameter '", popts$xlab, "' not one of the curve variables\n")
          if (!is.null(session)) updateConsoleLog(session, msg)
          return(NA)
        }
      }
      if (cnames[1] == popts$ylab) cly <- 1
      else {
        if (cnames[2] == popts$ylab) cly <- 2
        else {
          msg <- paste0("Curve plotting skipped: parameter '", popts$ylab, "' not one of the curve variables\n")
          if (!is.null(session)) updateConsoleLog(session, msg)
          return(NA)
        }
      }
      colindx <- match(curvelist[[i]]$type, c("BP", "HP", "LP"))
      lines(curvelist[[i]]$points[,clx], curvelist[[i]]$points[,cly], col=popts$colors[colindx], lwd=popts$lwd)

      if (!is.null(curvelist[[i]]$special.points)) {
        lbls <- c(curvelist[[i]]$special.tags[,1])
        bps <- (lbls %in% c("BT", "CP"))
        if (any(bps)) {
          x <- curvelist[[i]]$special.points[bps,clx]
          y <- curvelist[[i]]$special.points[bps,cly]
          points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
          text(x, y, labels=lbls[bps], pos = popts$biflblpos)
        }
      }
    })
    legend("topright", legend=c("BP", "HP", "LP"), col=popts$colors[c(1, 2, 3)], lty=1, lwd=popts$lwd, cex=popts$cex.legend, bg = "white")
  }
}

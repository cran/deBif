bifOrbitplot <- function(session = NULL, curvelist = NULL, popts) {
  # Save plot options to restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Plot the computed time series
  if (popts$plot3d == 1) {
    if (popts$theta < 0) par(oma=c(0,0,0,0), cex = popts$cex, mar = c(2.5,5,1,2))
    else par(oma=c(0,0,0,0), cex = popts$cex, mar = c(2.5,2,1,5))

    get_axp <- function(x) 10^c(ceiling(x[1]), floor(x[2]))
    pmat <- persp( c(0,1), c(0,1),  matrix(0, nrow=2, ncol=2), scale = FALSE,
                   zlim=c(0,1),  xlab='', ylab='', zlab='', ticktype='detailed', box=F, axes=F, border = NA,
                   expand = 0.8, col=NULL, shade=NA, theta=popts$theta, phi=5, d = 2)

    # Plot the x-axis
    alims <- c(popts$xmin,popts$xmax)
    if (popts$logx == 1) {
      usr.i <- log10(alims)
      aT.x <- axTicks(side = 1, usr = usr.i, axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = 5)
    } else {
      aT.x <- pretty(alims)
    }
    aT.xval <- converty2y(aT.x, 0, 1, 0, popts$xmin, popts$xmax, popts$logx)
    aT.xlab <- as.character(pretty(aT.x))
    # Axis
    lines(trans3d(c(0:1), 0, 0, pmat) , col="black")
    tick.start <- trans3d(aT.xval, 0, 0, pmat)
    tick.end <- trans3d(aT.xval, popts$tcl.len, 0, pmat)
    segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)
    label.pos <- trans3d(aT.xval, -0.02, -0.02, pmat)
    if (popts$theta < 0){
      text(label.pos$x, label.pos$y, labels=aT.xlab, adj=c(0.2,1), cex=popts$cex.axis, xpd = T)
    } else {
      text(label.pos$x, label.pos$y, labels=aT.xlab, adj=c(1,1), cex=popts$cex.axis, xpd = T)
    }
    title.pos <- trans3d(0.5, -0.08, -0.08, pmat)
    text(title.pos$x, title.pos$y, labels=popts$xlab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd = T)
    # Dotted backplane
    tick.start <- trans3d(aT.xval, 0, 0, pmat)
    tick.end <- trans3d(aT.xval, 1, 0, pmat)
    segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')
    tick.start <- trans3d(aT.xval, 1, 0, pmat)
    tick.end <- trans3d(aT.xval, 1, 1, pmat)
    segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')

    # Plot the vertical (z-)axis
    alims <- c(popts$ymin,popts$ymax)
    if (popts$logy == 1) {
      usr.i <- log10(alims)
      aT.z <- axTicks(side = 2, usr = usr.i, axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = 5)
    } else {
      aT.z <- pretty(alims)
    }
    aT.zval <- converty2y(aT.z, 0, 1, 0, popts$ymin, popts$ymax, popts$logy)
    aT.zlab <- as.character(pretty(aT.z))
    # Axis
    if (popts$theta < 0){
      lines(trans3d(0, 1, c(0:1), pmat) , col="black")
      tick.start <- trans3d(0, 1, aT.zval, pmat)
      tick.end <- trans3d(popts$tcl.len, 1, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

      # Dotted back plane
      tick.start <- trans3d(1, 1, aT.zval, pmat)
      tick.end <- trans3d(1, 0, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')
      tick.start <- trans3d(0, 1, aT.zval, pmat)
      tick.end <- trans3d(1, 1, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')

      label.pos <- trans3d(-0.02, 1.02, aT.zval, pmat)
      text(label.pos$x, label.pos$y, labels=aT.zlab, adj=c(1, 0), cex=popts$cex.axis, xpd = T)
      title.pos <- trans3d(-0.13, 1.13, 0.5, pmat)
      text(title.pos$x, title.pos$y, labels=popts$ylab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd=T)
    } else {
      lines(trans3d(1, 1, c(0:1), pmat) , col="black")
      tick.start <- trans3d(1, 1, aT.zval, pmat)
      tick.end <- trans3d(1-popts$tcl.len, 1, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

      # Dotted back plane
      tick.start <- trans3d(0, 1, aT.zval, pmat)
      tick.end <- trans3d(0, 0, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')
      tick.start <- trans3d(0, 1, aT.zval, pmat)
      tick.end <- trans3d(1, 1, aT.zval, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')

      label.pos <- trans3d(1.02, 1.02, aT.zval, pmat)
      text(label.pos$x, label.pos$y, labels=aT.zlab, adj=c(0, 0), cex=popts$cex.axis, xpd = T)
      title.pos <- trans3d(1.13, 1.13, 0.5, pmat)
      text(title.pos$x, title.pos$y, labels=popts$ylab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd=T)
    }

    # Plot the depth (y-)axis
    alims <- c(popts$y2min,popts$y2max)
    if (popts$logy2 == 1) {
      usr.i <- log10(alims)
      aT.y <- axTicks(side = 2, usr = usr.i, axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = 5)
    } else {
      aT.y <- pretty(alims)
    }
    aT.yval <- converty2y(aT.y, 0, 1, 0, popts$y2min, popts$y2max, popts$logy2)
    aT.ylab <- as.character(pretty(aT.y))
    if (popts$theta < 0) {
      # Axis
      lines(trans3d(0, c(0:1), 0, pmat) , col="black")
      tick.start <- trans3d(0, aT.yval, 0, pmat)
      tick.end <- trans3d(popts$tcl.len, aT.yval, 0, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

      # Dotted back plane
      tick.start <- trans3d(0, aT.yval, 0, pmat)
      tick.end <- trans3d(1, aT.yval, 0, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')
      tick.start <- trans3d(1, aT.yval, 0, pmat)
      tick.end <- trans3d(1, aT.yval, 1, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')

      label.pos <- trans3d(-0.02, aT.yval, -0.02, pmat)
      text(label.pos$x, label.pos$y, labels=aT.ylab, adj=c(1,0.8), cex=popts$cex.axis, xpd=T)
      title.pos <- trans3d(-0.1, 0.5, -0.1, pmat)
      text(title.pos$x, title.pos$y, labels=popts$y2lab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd=T)
    } else {
      lines(trans3d(1, c(0:1), 0, pmat) , col="black")
      tick.start <- trans3d(1, aT.yval, 0, pmat)
      tick.end <- trans3d(1-popts$tcl.len, aT.yval, 0, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

      # Dotted back plane
      tick.start <- trans3d(0, aT.yval, 0, pmat)
      tick.end <- trans3d(1, aT.yval, 0, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')
      tick.start <- trans3d(0, aT.yval, 0, pmat)
      tick.end <- trans3d(0, aT.yval, 1, pmat)
      segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y, lty='dotted', lwd=1, col='grey')

      label.pos <- trans3d(1.02, aT.yval, -0.02, pmat)
      text(label.pos$x, label.pos$y, labels=aT.ylab, adj=c(0,1), cex=popts$cex.axis, xpd=T)
      title.pos <- trans3d(1.1, 0.5, -0.1, pmat)
      text(title.pos$x, title.pos$y, labels=popts$y2lab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd=T)
    }

    if (!is.null(curvelist) && (length(curvelist) > 0)) {
      lapply((1:length(curvelist)), function(i) {
        cnames <- colnames(curvelist[[i]]$points)
        if (popts$xlab %in% cnames) {
          xcol <- match(popts$xlab, cnames)
        } else {
          msg <- paste0("Curve plotting skipped: parameter '", popts$xlab, "' not one of the curve variables\n")
          if (!is.null(session)) updateConsoleLog(session, msg)
          return(NA)
        }
        if (cnames[as.numeric(popts$ycol)] != popts$ylab) {
          msg <- paste0("Curve plotting skipped: variable '", popts$ylab, "' not one of the curve variables\n")
          if (!is.null(session)) updateConsoleLog(session, msg)
          return(NA)
        }
        if (cnames[as.numeric(popts$y2col)] != popts$y2lab) {
          msg <- paste0("Curve plotting skipped: variable '", popts$y2lab, "' not one of the curve variables\n")
          if (!is.null(session)) updateConsoleLog(session, msg)
          return(NA)
        }
        x <- converty2y(curvelist[[i]]$points[,xcol], 0, 1, 0, popts$xmin, popts$xmax, popts$logx)
        y <- converty2y(curvelist[[i]]$points[,popts$y2col], 0, 1, 0, popts$y2min, popts$y2max, popts$logy2)
        z <- converty2y(curvelist[[i]]$points[,popts$ycol], 0, 1, 0, popts$ymin, popts$ymax, popts$logy)
        lines(trans3d(x, y, z, pmat), col=popts$colors[1], lwd=popts$lwd)
      })
    }
  } else {
    if ((popts$ycol > 1) && (popts$y2col > 1)) {
      par(cex = popts$cex, mar = as.numeric(c(1.6*(popts$cex.lab+popts$cex.axis), 1.6*(popts$cex.lab+popts$cex.axis), 2, 1.6*(popts$cex.lab+popts$cex.axis))))
      logxy <- ifelse(popts$logx == 1, ifelse(popts$logy2 == 1, "xy", "x"), ifelse(popts$logy2 == 1, "y", ""))
      plot(NULL, type='n', xlab="", ylab="", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
           xlim=c(popts$xmin,popts$xmax), ylim=c(popts$y2min,popts$y2max), log=logxy)
      axis(4, cex.axis=popts$cex.axis)
      mtext(popts$y2lab, side = 4, line = (popts$cex.lab+popts$cex.axis), cex=(popts$cex.lab*popts$cex))
      par(new = TRUE)
    } else {
      par(cex = popts$cex, mar = as.numeric(c(1.6*(popts$cex.lab+popts$cex.axis), 1.6*(popts$cex.lab+popts$cex.axis), 2, 2)))
    }

    logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
    plot(NULL, type='n', xlab="", ylab="", xaxs = "i", yaxs = "i",
         xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy,
         cex.axis=popts$cex.axis, mgp = c((popts$cex.lab+popts$cex.axis), popts$cex.axis, 0))
    title(xlab=popts$xlab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))
    title(ylab=popts$ylab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))

    if (!is.null(curvelist) && (length(curvelist) > 0)) {
      if (popts$ycol == 1) {
        lapply((1:length(curvelist)), function(i) {
          cnames <- colnames(curvelist[[i]]$points)
          if (popts$xlab %in% cnames) {
            xcol <- match(popts$xlab, cnames)
            lapply(2:ncol(curvelist[[i]]$points), function(j) {
              lines(curvelist[[i]]$points[,xcol], curvelist[[i]]$points[,j], col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
            })
          } else {
            msg <- paste0("Curve plotting skipped: '", popts$xlab, "' not one of the curve variables\n")
            if (!is.null(session)) updateConsoleLog(session, msg)
            return(NA)
          }
        })
        legend("topright", legend=colnames(curvelist[[1]]$points)[2:ncol(curvelist[[1]]$points)], col=popts$colors[1:(ncol(curvelist[[1]]$points)-1)], lty=1, lwd=popts$lwd, cex=popts$cex.legend, bg = "white")
      } else {
        lapply((1:length(curvelist)), function(i) {
          cnames <- colnames(curvelist[[i]]$points)
          if (popts$xlab %in% cnames) {
            xcol <- match(popts$xlab, cnames)
            if (cnames[as.numeric(popts$ycol)] != popts$ylab) {
              msg <- paste0("Curve plotting skipped: variable '", popts$ylab, "' not one of the curve variables\n")
              if (!is.null(session)) updateConsoleLog(session, msg)
              return(NA)
            }
            lines(curvelist[[i]]$points[,xcol], curvelist[[i]]$points[,popts$ycol], col=popts$colors[1], lwd=popts$lwd)
          } else {
            msg <- paste0("Curve plotting skipped: '", popts$xlab, "' not one of the curve variables\n")
            if (!is.null(session)) updateConsoleLog(session, msg)
            return(NA)
          }
        })
        if (popts$y2col > 1) {
          lapply((1:length(curvelist)), function(i) {
            cnames <- colnames(curvelist[[i]]$points)
            if (popts$xlab %in% cnames) {
              xcol <- match(popts$xlab, cnames)
              if (cnames[as.numeric(popts$y2col)] != popts$y2lab) {
                msg <- paste0("Curve plotting skipped: variable '", popts$y2lab, "' not one of the curve variables\n")
                if (!is.null(session)) updateConsoleLog(session, msg)
                return(NA)
              }
              lines(curvelist[[i]]$points[,xcol],
                    converty2y(curvelist[[i]]$points[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2),
                    col=popts$colors[2], lwd=popts$lwd)
            } else {
              msg <- paste0("Curve plotting skipped: '", popts$xlab, "' not one of the curve variables\n")
              if (!is.null(session)) updateConsoleLog(session, msg)
              return(NA)
            }
          })
          legend("topright", legend=colnames(curvelist[[1]]$points)[c(popts$ycol, popts$y2col)], col=popts$colors[c(1, 2)], lty=1, lwd=popts$lwd, cex=popts$cex.legend, bg = "white")
        }
      }
    }
  }
}

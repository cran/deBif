cliptrans3d <- function(x0, y0, z0, ...) {
  x <- as.numeric(c(x0))
  xin <- ((x >= 0) & (x <= 1))
  y <- as.numeric(c(y0))
  yin <- ((y >= 0) & (y <= 1))
  z <- as.numeric(c(z0))
  zin <- ((z >= 0) & (z <= 1))
  allin <- xin & yin & zin
  if (any(allin))
    return(trans3d(x[allin], y[allin], z[allin], ...))
  else return(NULL)
}

bif1parplot <- function(session = NULL, curvelist = NULL, popts) {
  # Save plot options to restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Plot the bifurcation curves
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
      text(label.pos$x, label.pos$y, labels=aT.ylab, adj=c(1,0.8), cex=popts$cex.axis, xpd = T)
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
      text(label.pos$x, label.pos$y, labels=aT.ylab, adj=c(0,1), cex=popts$cex.axis, xpd = T)
      title.pos <- trans3d(1.1, 0.5, -0.1, pmat)
      text(title.pos$x, title.pos$y, labels=popts$y2lab, adj=c(0.5,0.5), cex=popts$cex.lab, xpd=T)
    }

    if (!is.null(curvelist) && (length(curvelist) > 0)) {
      lapply((1:length(curvelist)), function(i) {
        cnames <- colnames(curvelist[[i]]$points)
        if (cnames[1] != popts$xlab) {
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
        if (curvelist[[i]]$type == "LC") {
          sdim <- length(curvelist[[i]]$initstate)
          mdim <- (length(curvelist[[i]]$points[1,]) - 2)/sdim
          mrange <- sdim*(0:(mdim-1))
          for (j in (1:nrow(curvelist[[i]]$points))){
            x <- converty2y(curvelist[[i]]$points[j,1], 0, 1, 0, popts$xmin, popts$xmax, popts$logx)
            y <- converty2y(c(curvelist[[i]]$points[j,popts$y2col + mrange]), 0, 1, 0, popts$y2min, popts$y2max, popts$logy2)
            z <- converty2y(c(curvelist[[i]]$points[j,popts$ycol + mrange]), 0, 1, 0, popts$ymin, popts$ymax, popts$logy)
            x <- rep(x, each=length(mrange))
            tr3d <- cliptrans3d(x, y, z, pmat)
            if (!is.null(tr3d)) lines(tr3d, colvar = NULL, col=popts$colors[2], lwd=1)
          }
        } else {
          evmax <- Re(curvelist[[i]]$eigvals[,1])
          sp <- (evmax < 0)
          x <- converty2y(curvelist[[i]]$points[,1], 0, 1, 0, popts$xmin, popts$xmax, popts$logx)
          y <- converty2y(curvelist[[i]]$points[,popts$y2col], 0, 1, 0, popts$y2min, popts$y2max, popts$logy2)
          z <- converty2y(curvelist[[i]]$points[,popts$ycol], 0, 1, 0, popts$ymin, popts$ymax, popts$logy)
          runlen <- rle(sp)$lengths
          runval <- rle(sp)$values
          runend <- cumsum(runlen)
          runstart <- (c(1, (1 + runend)))[1:length(runend)]
          for (j in (1:length(runstart))) {
            tr3d <- cliptrans3d(x[runstart[j]:runend[j]], y[runstart[j]:runend[j]], z[runstart[j]:runend[j]], pmat)
            if (!is.null(tr3d)) {
              if (runval[j])
                lines(tr3d, colvar = NULL, col=popts$colors[1], lwd=popts$lwd)
              else
                lines(tr3d, colvar = NULL, col=popts$colors[1], lty=popts$unstablelty, lwd=popts$lwd)
            }
          }
          if (!is.null(curvelist[[i]]$special.points)) {
            lbls <- c(curvelist[[i]]$special.tags[,1])
            bps <- (lbls %in% c("BP", "HP", "LP"))
            if (any(bps)) {
              x <- converty2y(curvelist[[i]]$special.points[bps,1], 0, 1, 0, popts$xmin, popts$xmax, popts$logx)
              y <- converty2y(curvelist[[i]]$special.points[bps,popts$y2col], 0, 1, 0, popts$y2min, popts$y2max, popts$logy2)
              z <- converty2y(curvelist[[i]]$special.points[bps,popts$ycol], 0, 1, 0, popts$ymin, popts$ymax, popts$logy)
              tr3d <- cliptrans3d(x, y, z, pmat)
              if (!is.null(tr3d)) {
                points(tr3d, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
                text(tr3d, labels=lbls[bps], pos = popts$biflblpos)
              }
            }
          }
        }
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
         cex.axis=popts$cex.axis, , mgp = c((popts$cex.lab+popts$cex.axis), popts$cex.axis, 0))
    title(xlab=popts$xlab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))
    title(ylab=popts$ylab, cex.lab=popts$cex.lab, line = (popts$cex.lab+popts$cex.axis))

    if (!is.null(curvelist) && (length(curvelist) > 0)) {
      if (popts$ycol == 1) {
        lapply((1:length(curvelist)), function(i) {
          cnames <- colnames(curvelist[[i]]$points)
          if (cnames[1] != popts$xlab) {
            msg <- paste0("Curve plotting skipped: parameter '", popts$xlab, "' not one of the curve variables\n")
            if (!is.null(session)) updateConsoleLog(session, msg)
            return(NA)
          }
          if (curvelist[[i]]$type == "LC") {
            sdim <- length(curvelist[[i]]$initstate)
            mdim <- (length(curvelist[[i]]$points[1,]) - 2)/sdim
            mrange <- sdim*(0:(mdim-1))
            lapply(2:(sdim+1), function(j) {
              x <- curvelist[[i]]$points[,1]
              if (nrow(curvelist[[i]]$points) > 1)
                y <- apply(curvelist[[i]]$points[, j+mrange], 1, min)
              else
                y <- min(curvelist[[i]]$points[, j+mrange])
              lines(x, y, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
              if (nrow(curvelist[[i]]$points) > 1)
                y <- apply(curvelist[[i]]$points[, j+mrange], 1, max)
              else
                y <- max(curvelist[[i]]$points[, j+mrange])
              lines(x, y, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
            })
          } else {
            lapply(2:ncol(curvelist[[i]]$points), function(j) {
              evmax <- Re(curvelist[[i]]$eigvals[,1])
              sp <- (evmax < 0)
              x <- curvelist[[i]]$points[,1]
              y <- curvelist[[i]]$points[,j]
              x[!sp] <- NA
              y[!sp] <- NA
              lines(x, y, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
              x <- curvelist[[i]]$points[,1]
              y <- curvelist[[i]]$points[,j]
              x[sp] <- NA
              y[sp] <- NA
              lines(x, y, lty=popts$unstablelty, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
              if (!is.null(curvelist[[i]]$special.points)) {
                lbls <- c(curvelist[[i]]$special.tags[,1])
                bps <- (lbls %in% c("BP", "HP", "LP"))
                if (any(bps)) {
                  x <- curvelist[[i]]$special.points[bps,1]
                  y <- curvelist[[i]]$special.points[bps,j]
                  points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
                  text(x, y, labels=lbls[bps], pos = popts$biflblpos)
                }
              }
            })
          }
        })
        legend("topright", legend=colnames(curvelist[[1]]$points)[2:ncol(curvelist[[1]]$points)], col=popts$colors[1:(ncol(curvelist[[1]]$points)-1)], lty=1, lwd=popts$lwd, cex=popts$cex.legend, bg = "white")
      } else {
        lapply((1:length(curvelist)), function(i) {
          cnames <- colnames(curvelist[[i]]$points)
          if (cnames[1] != popts$xlab) {
            msg <- paste0("Curve plotting skipped: parameter '", popts$xlab, "' not one of the curve variables\n")
            if (!is.null(session)) updateConsoleLog(session, msg)
            return(NA)
          }
          if (cnames[as.numeric(popts$ycol)] != popts$ylab) {
            msg <- paste0("Curve plotting skipped: variable '", popts$ylab, "' not one of the curve variables\n")
            if (!is.null(session)) updateConsoleLog(session, msg)
            return(NA)
          }
          if (curvelist[[i]]$type == "LC") {
            sdim <- length(curvelist[[i]]$initstate)
            mdim <- (length(curvelist[[i]]$points[1,]) - 2)/sdim
            mrange <- sdim*(0:(mdim-1))
            x <- curvelist[[i]]$points[,1]
            if (nrow(curvelist[[i]]$points) > 1)
              y <- apply(curvelist[[i]]$points[, popts$ycol+mrange], 1, min)
            else
              y <- min(curvelist[[i]]$points[, popts$ycol+mrange])
            lines(x, y, col=popts$colors[1], lwd=popts$lwd)
            if (nrow(curvelist[[i]]$points) > 1)
              y <- apply(curvelist[[i]]$points[, popts$ycol+mrange], 1, max)
            else
              y <- max(curvelist[[i]]$points[, popts$ycol+mrange])
            lines(x, y, col=popts$colors[1], lwd=popts$lwd)
          } else {
            evmax <- Re(curvelist[[i]]$eigvals[,1])
            sp <- (evmax < 0)
            x <- curvelist[[i]]$points[,1]
            y <- curvelist[[i]]$points[,popts$ycol]
            x[!sp] <- NA
            y[!sp] <- NA
            lines(x, y, col=popts$colors[1], lwd=popts$lwd)
            x <- curvelist[[i]]$points[,1]
            y <- curvelist[[i]]$points[,popts$ycol]
            x[sp] <- NA
            y[sp] <- NA
            lines(x, y, lty=popts$unstablelty, col=popts$colors[1], lwd=popts$lwd)
            if (!is.null(curvelist[[i]]$special.points)) {
              lbls <- c(curvelist[[i]]$special.tags[,1])
              bps <- (lbls %in% c("BP", "HP", "LP"))
              if (any(bps)) {
                x <- curvelist[[i]]$special.points[bps,1]
                y <- curvelist[[i]]$special.points[bps,popts$ycol]
                points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
                text(x, y, labels=lbls[bps], pos = popts$biflblpos)
              }
            }
          }
        })
        if (popts$y2col > 1) {
          lapply((1:length(curvelist)), function(i) {
            cnames <- colnames(curvelist[[i]]$points)
            if (cnames[as.numeric(popts$y2col)] != popts$y2lab) {
              msg <- paste0("Curve plotting skipped: variable '", popts$y2lab, "' not one of the curve variables\n")
              if (!is.null(session)) updateConsoleLog(session, msg)
              return(NA)
            }
            if (curvelist[[i]]$type == "LC") {
              sdim <- length(curvelist[[i]]$initstate)
              mdim <- (length(curvelist[[i]]$points[1,]) - 2)/sdim
              mrange <- sdim*(0:(mdim-1))
              x <- curvelist[[i]]$points[,1]
              if (nrow(curvelist[[i]]$points) > 1)
                y <- apply(curvelist[[i]]$points[, popts$y2col+mrange], 1, min)
              else
                y <- min(curvelist[[i]]$points[, popts$y2col+mrange])
              y <- converty2y(y, popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
              lines(x, y, col=popts$colors[2], lwd=popts$lwd)
              if (nrow(curvelist[[i]]$points) > 1)
                y <- apply(curvelist[[i]]$points[, popts$y2col+mrange], 1, max)
              else
                y <- max(curvelist[[i]]$points[, popts$y2col+mrange])
              y <- converty2y(y, popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
              lines(x, y, col=popts$colors[2], lwd=popts$lwd)
            } else {
              evmax <- Re(curvelist[[i]]$eigvals[,1])
              sp <- (evmax < 0)
              x <- curvelist[[i]]$points[,1]
              y <- converty2y(curvelist[[i]]$points[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
              x[!sp] <- NA
              y[!sp] <- NA
              lines(x, y, col=popts$colors[2], lwd=popts$lwd)
              x <- curvelist[[i]]$points[,1]
              y <- converty2y(curvelist[[i]]$points[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
              x[sp] <- NA
              y[sp] <- NA
              lines(x, y, lty=popts$unstablelty, col=popts$colors[2], lwd=popts$lwd)
              if (!is.null(curvelist[[i]]$special.points)) {
                lbls <- c(curvelist[[i]]$special.tags[,1])
                bps <- (lbls %in% c("BP", "HP", "LP"))
                if (any(bps)) {
                  x <- curvelist[[i]]$special.points[bps,1]
                  y <- converty2y(curvelist[[i]]$special.points[bps,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
                  points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
                  text(x, y, labels=lbls[bps], pos = popts$biflblpos)
                }
              }
            }
          })
          legend("topright", legend=colnames(curvelist[[1]]$points)[c(popts$ycol, popts$y2col)], col=popts$colors[c(1, 2)], lty=1, lwd=popts$lwd, cex=popts$cex.legend, bg = "white")
        }
      }
    }
  }
}

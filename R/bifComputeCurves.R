# Integration of ODE system forward or backward
computeTimeseries <- function(session, model, state, parms, clist, pointid, nopts) {

  curvescomputed <- as.numeric(clist[['TotalCurves']])

  if (pointid > 0) {
    ind1 <- round(pointid/1000000)
    ind2 <- round((pointid-ind1*1000000)/1000)
    ind3 <- round(pointid-ind1*1000000-ind2*1000)
    cln1 <- (c('Orbits', 'BifurcationCurves', 'BifurcationBounds'))[ind1]
    ii <- ifelse((ind1 == 3), 2, 1)   # 2 parameter bifurcation points have 2 columns before the state, otherwise 1 only

    initstate <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (ii + (1:length(state)))])
    initparms <- as.numeric(clist[[cln1]][[ind2]]$parameters)
    inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
    if (inittype != "TS")
      initparms[as.numeric(clist[[cln1]][[ind2]]$bifpars)] <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (1:ii)])
    names(initstate) <- names(state)
    names(initparms) <- names(parms)
  } else {
    initstate <- state
    initparms <- parms
    inittype <- "US"
  }

  newcurvenr <- length((clist[['Orbits']]))+1

  times <- seq(0, nopts$tmax, by=abs(nopts$tstep))
  if (nopts$tstep < 0.0) times <- nopts$tmax - times
  nsol <- as.data.frame(do.call('ode', c(list(times=times, func=model, y=initstate, parms=initparms), method=nopts$odemethod)))

  names(nsol) <- c("Time", names(state))

  startPnt <- c("Type" = "TS",
                "Description" = paste0('T=', times[1], ' ',
                                       paste(unlist(lapply(1:length(state),
                                                           function(i) {paste0(names(state[i]), "=", round(nsol[1, (1+i)], 3))})),
                                             collapse=' ')))
  endPnt <- c("Type" = "TS",
              "Description" = paste0('T=', times[length(times)], ' ',
                                     paste(unlist(lapply(1:length(state),
                                                         function(i) {paste0(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 3))})),
                                           collapse=' ')))

  updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
  curvescomputed <- curvescomputed + 1

  lbl <- paste0("TS", sprintf("%02d", curvescomputed),": ", startPnt["Description"])
  startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])
  endPnt["Description"] <- paste0(sprintf("%04d: ", nrow(nsol)), endPnt["Description"])

  newcurve <- list(label = lbl, type = "TS", initstate = initstate, parameters = initparms, points = nsol,
                   special.points = rbind(c(nsol[1,]), c(nsol[nrow(nsol),])), special.tags = rbind(startPnt, endPnt))

  clist$Orbits[[newcurvenr]] <- newcurve
  clist$TotalCurves <- curvescomputed

  return(clist)
}


# Find first points on curve for continuation purposes
initCurveContinuation <- function(session, model, initstate, initparms, tanvec, curtabname, clist,
                                  curvetype, inittype, popts, nopts, reportlevel, direction) {

  if (exists("deBifverbose", envir = .GlobalEnv)) verbose <- get("deBifverbose", envir = .GlobalEnv)
  else verbose <- FALSE

  if (!(curvetype %in% c("EQ", "BP", "HP", "LP", "LC"))) {
    msg <- paste0("Computation aborted:\nContinuation for curve type ", curvetype, " not implemented\n")
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
    return(NULL)
  }

  if (curtabname == 'BifurcationCurves') {
    freepars = c(as.numeric(popts[["xcol"]]))
  } else {
    freepars = c(as.numeric(popts[["xcol"]]), as.numeric(popts[["ycol"]]))
  }

  freeparsdim <- length(freepars)
  statedim <- length(initstate)

  y <- c(initparms[freepars], initstate)
  fixedpars <- initparms[!(1:length(initparms)) %in% freepars]

  varnames <- names(y)
  eignames <- unlist(lapply((1:statedim), function(i){paste0("Eigenvalue", i)}))
  tvnames <- unlist(lapply((1:length(y)), function(i){paste0("d", varnames[i])}))

  if ((curvetype == "EQ") || (curvetype == "LC")) {
    condfun <- NULL
  } else {
    if (exists(paste0(curvetype, "continuation"), mode = "function"))
      condfun <- list(get(paste0(curvetype, "continuation"), mode = "function"))
    else {
      msg <- paste0("Computation aborted:\nAdditional condition function for curve type ", curvetype, " not found\n")
      if (!is.null(session)) updateConsoleLog(session, msg)
      else cat(msg)
      return(NULL)
    }
  }

  # Test functions will only be implemented for EQ curves
  if (exists(paste0("analyse", curvetype), mode = "function")) {
    analysefun <- get(paste0("analyse", curvetype), mode = "function")
    testvals <- NULL
  }
  else analysefun <- NULL

  cData <- list()
  cData$guess <- y
  cData$yold <- y
  cData$tanvec <- tanvec
  cData$testvals <- NULL
  cData$pntnr <- 1
  cData$stepsize <- as.numeric(direction)*abs(nopts$maxstepsize) / 10.0
  cData$direction <- as.numeric(direction)

  cData$tabname <- curtabname
  cData$model <- model
  cData$fixedpars <- fixedpars
  cData$condfun <- condfun
  cData$extSys <- switch(curvetype, "BP" = ExtSystemEQ, "EQ" = ExtSystemEQ, "HP" = ExtSystemEQ, "LP" = ExtSystemEQ, "LC" = ExtSystemLC)

  if (curvetype == "LC") cData$jacfun <- ExtSystemLCjac
  else cData$jacfun <- jacobian.full
  cData$analysefun <- analysefun
  cData$statedim <- statedim
  cData$freeparsdim <- freeparsdim
  cData$pointdim <- freeparsdim + statedim
  cData$inittype <- inittype
  cData$curvetype <- curvetype
  cData$varnames <- varnames
  cData$eignames <- eignames
  cData$tvnames <- tvnames
  cData$reportlevel <- reportlevel
  cData$newcurvenr <- length((clist[[curtabname]]))+1

  # If present run the initializer
  if (exists(paste0("init", curvetype), mode = "function")) {
    initfun <- get(paste0("init", curvetype), mode = "function")
    names(y) <- varnames
    initres <- initfun(y, fixedpars, cData, nopts, session)
    if (is.null(initres)) {
      return(NULL)
    } else {
      if (!is.null(initres$curveData)) cData <- initres$curveData
      if (!is.null(initres$y)) {
        y <- c(as.numeric(initres$y))
        names(y) <- cData$varnames
      }
      if (!is.null(initres$tanvec)) tanvec <- c(as.numeric(initres$tanvec))
    }
  }

  if (is.null(tanvec) || (length(tanvec) != length(y)))
    tanvec <- c(rep(1.0, length(freepars)), rep(0, (length(y) - length(freepars))))
  else if (abs(tanvec[1]) > as.numeric(nopts$iszero)) tanvec <- sign(tanvec[1])*tanvec

  cData$guess <- y
  cData$yold <- y
  cData$tanvec <- tanvec

  nsol <- tryCatch(nextCurvePoints(1, cData, popts, nopts, session = session),
                   warning = function(e) {
                     msg <- gsub(".*:", "Warning in nextCurvePoints:", e)
                     if (!is.null(session)) updateConsoleLog(session, msg)
                     else cat(msg)
                     return(NULL)
                   },
                   error = function(e) {
                     msg <- gsub(".*:", "Error in nextCurvePoints:", e)
                     if (!is.null(session)) updateConsoleLog(session, msg)
                     else cat(msg)
                     return(NULL)
                   })

  if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {
    if (curvetype == "LC") {
      if (sign(cData$stepsize * nsol$tanvec[1]) != direction) {
        session$userData$curveData$stepsize <- -1.0 * as.numeric(cData$stepsize)
      }

      vals <- lapply((1:statedim), function(i) {
        indxrange <- statedim*(1:(nopts$ninterval*nopts$glorder))
        y <- nsol$points[1, freeparsdim+i+indxrange]
        yname <- names(nsol$points[1, freeparsdim+i])
        paste0("Min.", yname, "=", round(min(y), 3), ", Max.", yname, "=", round(max(y), 3))
      })
      startPnt <- c("Type" = curvetype,
                    "Description" = paste0(names(nsol$points[1, 1]), "=", round(nsol$points[1, 1], 3), " ",
                                           names(nsol$points[1, ncol(nsol$points)]), "=",
                                           round(nsol$points[1, ncol(nsol$points)], 3), " ",
                                           paste(unlist(vals), collapse = ', ')))
    } else {
      startPnt <- c("Type" = curvetype,
                    "Description" = paste(unlist(lapply(1:length(nsol$points[1,]),
                                                        function(i) {paste0(names(nsol$points[1,i]), "=",
                                                                            round(nsol$points[1, i], 3))})),
                                          collapse=' '))
    }
    lbl <- paste0(curvetype, sprintf("%02d", (clist$TotalCurves + 1)),": ", startPnt["Description"])
    startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])

    newcurve <- list(label = lbl, type = curvetype, initstate = initstate, parameters = initparms, bifpars = freepars,
                     points = nsol$points, eigvals = nsol$eigvals, tanvec = nsol$tanvec,
                     special.points = nsol$points, special.eigvals = nsol$eigvals,
                     special.tanvec = nsol$tanvec, special.tags = rbind(NULL, c(startPnt)))

    clist[[cData$tabname]][[cData$newcurvenr]] <- newcurve
    clist$TotalCurves <- (clist$TotalCurves + 1)
    return(clist)
  }

  return(NULL)
}

# Extend curve during continuation purposes
nextCurvePoints <- function(maxpoints, curveData, popts, nopts, session = NULL) {

  if (exists("deBifverbose", envir = .GlobalEnv) && (get("deBifverbose", envir = .GlobalEnv)))
    verbose <- TRUE
  else verbose <- FALSE

  cData <- curveData
  curvetype <- cData$curvetype
  pntnr <- cData$pntnr
  statedim <- cData$statedim
  freeparsdim <- cData$freeparsdim

  if (curvetype == "LC") maxpoints <- 1

  allsols <- NULL
  alltvs <- NULL
  alleigs <- NULL
  specialsols <- NULL
  specialtvs <- NULL
  specialeigs <- NULL
  specialtags <- NULL

  corrections <- 1
  while ((pntnr - cData$pntnr) < as.numeric(maxpoints)) {
    result <- tryCatch(.Call("deBif", curvetype, cData$model, cData$guess, cData$fixedpars, cData$tanvec,
                             nopts$rhstol, nopts$dytol, nopts$jacdif, nopts$maxiter, nopts$glorder, nopts$ninterval,
                             cData[c("statedim", "finemeshdim", "upoldp", "wi", "wt", "wp", "wpvec")]),
                       warning = function(e) {
                         msg <- gsub(".*:", "Warning in deBif:", e)
                         if (!is.null(session)) updateConsoleLog(session, msg)
                         else cat(msg)
                         return(NULL)
                       },
                       error = function(e) {
                         msg <- gsub(".*:", "Error in deBif:", e)
                         if (!is.null(session)) updateConsoleLog(session, msg)
                         else cat(msg)
                         return(NULL)
                       })
    if (is.null(result)) {
      newsolution <- FALSE
    } else {
      newsolution <- TRUE
      if ((cData$pntnr > 1) && (curvetype != "LC")) {
        nonzeroyold <- (abs(cData$yold) > as.numeric(nopts$iszero))
        if (length(nonzeroyold) > cData$pointdim) nonzeroyold[((cData$pointdim+1):length(nonzeroyold))] <- FALSE
        if (any(nonzeroyold))
          newsolution <- (max((cData$yold[nonzeroyold] - result$y[nonzeroyold])^2/cData$yold[nonzeroyold]^2) < nopts$neartol)
      }
    }
    if (!is.null(result) && newsolution) {
      y <- result$y
      iternr <- result$niter

      if ("Jacobian" %in% names(result)) {
        jac <- result$Jacobian
      } else {
        # Compute the Jacobian w.r.t. to free parameters (the first cData$freeparsdim
        # columns) and the state variables
        jac <- cData$jacfun(y=y, func=cData$extSys, parms=cData$fixedpars, pert = nopts$jacdif, curveData = cData, nopts = nopts)
      }

      if (curvetype == "LC") {
        cData$upoldp <- updateRefSol(0, y, cData$fixedpars, cData, nopts)
      } else {
        # Compute the eigenvalues of the restricted Jacobian (exclude the first cData$freeparsdim columns
        # and takink only as many rows as there are state variables)
        eig <- eigen(jac[(1:cData$statedim),((cData$freeparsdim+1):cData$pointdim)])
        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <- cData$eignames
      }

      if ("tanvec" %in% names(result)) {
        cData$tanvec <- result$tanvec
      } else {
        # Append the current tangent vector as the last row to the jacobian to
        # preserve direction. See the matcont manual at
        # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
        # Notice that cData$jacfun returns a square matrix with NA values on the last row
        jac[nrow(jac),] <- cData$tanvec
        if (rcond(jac) > nopts$rhstol) {
          tvnew <- solve(jac, c(rep(0, (length(cData$tanvec)-1)), 1))
          tvnorm <- sqrt(sum(tvnew^2))
          tvnew <- tvnew/tvnorm
          names(tvnew) <- cData$tvnames
          cData$tanvec <- tvnew
        } else tvnew <- cData$tanvec
      }
      ############## Execute the test functions and store the results
      if (!is.null(cData$analysefun)) {
        testvals <- cData$analysefun(state = y, parms = cData$fixedpars, cData, nopts, session = session)
        testvals2report <- testvals
        if (!is.null(testvals) && ("y" %in% names(testvals) > 0) && (length(testvals$y) > 0)) {

          allsols <- rbind(allsols, c(testvals$y[1:cData$pointdim]))
          alltvs <- rbind(alltvs, c(testvals$tanvec[1:cData$pointdim]))
          if (!is.null(testvals$eigval)) {
            alleigs <- rbind(alleigs, c(testvals$eigval[1:length(eigval)]))
          }

          specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
          specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
          if (!is.null(testvals$eigval)) {
            specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))
          }

          dscp <- paste(unlist(lapply(1:cData$pointdim, function(i) {paste0(names(testvals$y[i]), "=", round(c(testvals$y)[i], 5))})),
                        collapse=', ')
          specialtags <- rbind(specialtags, c("Type" = testvals$biftype, "Description" = paste0(testvals$biftype, ": ", dscp)))

          msg <- switch(testvals$biftype,
                        BP = "Branching", HP = "Hopf bifurcation", LP = "Limit", BT = "Bogdanov-Takens", CP = "Cusp")
          msg <- paste0("Solution ", pntnr, ": ", msg, " point found:\n", dscp, "\n")

          if (curvetype != "LC") {
            msg <- paste0(msg, "Eigenvalues:\n",
                          paste(unlist(lapply(1:length(testvals$eigval),
                                              function(i) {rcprintf("%12.5E", testvals$eigval[i])})), collapse=' '), "\n")
          }
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)

          if (cData$reportlevel >= 1) {
            yy <- c(as.numeric(testvals$y[1:cData$pointdim]), testvals$eigval)
            if (cData$reportlevel == 2) {
              specvar <- c("bpval", "hpval", "lpval", "btval", "cpval")
              for (testname in specvar) {
                if (testname %in% names(testvals)) yy <- c(yy, as.numeric(testvals[[testname]]))
              }
            }
            names(yy) <- NULL
            cat(paste(unlist(lapply(yy, function(x) {rcprintf("%12.5E", x)})), collapse = " "),
                sprintf("**%s**\n", testvals$biftype))
          }
          # Reset the detected bifurcation test value to NULL
          specvar <- c("bpval", "hpval", "lpval", "btval", "cpval")
          speclbl <- c("BP", "HP", "LP", "BT", "CP")
          ii <- (1:length(speclbl))[testvals$biftype == speclbl]
          testvals[[specvar[ii]]] <- NULL

          testvals$y <- NULL
          testvals$tanvec <- NULL
          if (curvetype != "LC") testvals$eigval <- NULL
          testvals$biftype <- NULL

          pntnr <- pntnr + 1
        }
        cData$testvals <- testvals
      }

      ############## Report the solution point and the eigenvalues
      if (curvetype == "LC") {
        msg <- paste0("Solution ", pntnr, " found:\n",
                      "     ",
                      paste(unlist(lapply(c((1:(cData$freeparsdim+cData$statedim)), cData$pointdim),
                                          function(i) {sprintf("%12s", names(y[i]))})), collapse = " "),
                      "\nMin. ",
                      sprintf("%12.5E", y[cData$freeparsdim]), " ",
                      paste(unlist(lapply((1:cData$statedim),
                                          function(i) {sprintf("%12.5E", min(y[cData$freeparsdim+i+cData$statedim*(1:(nopts$ninterval*nopts$glorder))]))})),
                            collapse = " "), " ",
                      sprintf("%12.5E", y[cData$pointdim]),
                      "\nMax. ",
                      sprintf("%12.5E", y[cData$freeparsdim]), " ",
                      paste(unlist(lapply((1:cData$statedim),
                                          function(i) {sprintf("%12.5E", max(y[cData$freeparsdim+i+cData$statedim*(1:(nopts$ninterval*nopts$glorder))]))})),
                            collapse = " "), " ",
                      sprintf("%12.5E", y[cData$pointdim]),
                      "\n")
      } else {
        msg <- paste0("Solution ", sprintf("%-5d", pntnr), ":  ",
                      paste(unlist(lapply(1:cData$pointdim,
                                          function(i) {paste0(names(y[i]), "=", sprintf("%12.5E", y[i]))})), collapse=', '), "\n")
        msg <- paste0(msg, "Eigenvalues   :  ",
                      paste(unlist(lapply(1:length(eigval), function(i) {rcprintf("%12.5E", eigval[i])})), collapse=' '), "\n")
      }
      if (!is.null(session)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
      else cat(msg)

      if (cData$reportlevel >= 1) {
        if (curvetype == "LC") {
          yy <- c(as.numeric(y[1:cData$freeparsdim]),
                  unlist(lapply((1:cData$statedim),
                                function(i) {c(min(y[cData$freeparsdim+i+cData$statedim*(1:(nopts$ninterval*nopts$glorder))]),
                                               max(y[cData$freeparsdim+i+cData$statedim*(1:(nopts$ninterval*nopts$glorder))]))})),
                  as.numeric(y[cData$pointdim]))
          if (pntnr == 1) {
            namesyy <- c(names(y[1:cData$freeparsdim]),
                         unlist(lapply((1:cData$statedim),
                                       function(i){c(paste0("min.", cData$varnames[i+1]),
                                                     paste0("max.", cData$varnames[i+1]))})),
                         names(y[cData$pointdim]))
            cat("\n\nPoint: ", paste(unlist(lapply(namesyy, function(x) {sprintf("%12s", x)})), collapse = " "), "\n")
          }
        } else {
          specvar <- c("bpval", "hpval", "lpval", "btval", "cpval")
          speclbl <- c("BP", "HP", "LP", "BT", "CP")
          if (pntnr == 1) {
            msg <- paste(unlist(lapply(names(y[1:cData$pointdim]), function(x) {sprintf("%12s", x)})), collapse = " ")
            for (ii in (1:length(cData$eignames))) {
              if (abs(Im(eigval[ii])) < nopts$iszero) msg <- paste(msg, sprintf("%12s", cData$eignames[ii]))
              else msg <- paste(msg, sprintf("%25s", cData$eignames[ii]))
            }
            if (cData$reportlevel == 2) {
              for (ii in (1:length(specvar))) {
                if (specvar[ii] %in% names(testvals2report)) msg <- paste(msg, sprintf("%7s test", speclbl[ii]))
              }
            }
            cat("\n\nPoint  ", msg, "\n")
          }
          yy <- c(as.numeric(y[1:cData$pointdim]), eigval)
          if (cData$reportlevel == 2) {
            for (testname in specvar) {
              if (testname %in% names(testvals2report)) yy <- c(yy, as.numeric(testvals2report[[testname]]))
            }
          }
        }
        cat(sprintf("%5d: ", pntnr), paste(unlist(lapply(yy, function(x) {rcprintf("%12.5E", x)})), collapse = " "), "\n")
      }

      ############## Store the results
      allsols <- rbind(allsols, (c(y)[1:cData$pointdim]))
      alltvs <- rbind(alltvs, (c(cData$tanvec)[1:cData$pointdim]))
      if (curvetype != "LC") alleigs <- rbind(alleigs, c(eigval))

      ############## Determine the new step along the curve
      if ((corrections == 1) && (iternr < 4)) {
        cData$stepsize <- sign(cData$stepsize) * min(1.3 * abs(cData$stepsize), abs(nopts$maxstepsize));
      }
      corrections <- 1
      if (curvetype == "LC") cData$guess <- y + cData$stepsize * cData$tanvec
      else cData$guess <- y + setStepSize(y, cData, as.numeric(nopts$minstepsize), as.numeric(nopts$dytol))
      cData$yold <- y

      ############## Stop the curve if outside the visible plotting region
      pntnr <- pntnr + 1

      ############## Stop the curve if maximum points reached or curve outside the visible plotting region
      curvedone <- FALSE
      if (pntnr > as.numeric(nopts$maxpoints)) {
        msg <- "Computation halted:\nMaximum number of points along the curve reached\n"
        curvedone <- TRUE
      }

      if ((!curvedone ) && (pntnr > 10)) {
        bndtol <- as.numeric(nopts$iszero)
        if ((!curvedone ) && ((as.numeric(y[1]) < (as.numeric(popts$xmin) - bndtol)) ||
                              (as.numeric(y[1]) > (as.numeric(popts$xmax) + bndtol)))) {
          msg <- "Computation halted:\nMinimum or maximum of x-axis domain reached\n"
          curvedone <- TRUE
        }
        if ((!curvedone ) && (curvetype == "EQ")) {
          if (popts$ycol == 1) {
            if ((!curvedone ) && (any(as.numeric(y[(2:(cData$statedim+1))]) < (as.numeric(popts$ymin) - bndtol)))) {
              msg <- "Computation halted:\nMinimum of y-axis domain reached for one of the y-axis variables\n"
              curvedone <- TRUE
            }
            if ((!curvedone ) && (any(as.numeric(y[(2:(cData$statedim+1))]) > (as.numeric(popts$ymax) + bndtol)))) {
              msg <- "Computation halted:\nMaximum of y-axis domain reached for one of the y-axis variables\n"
              curvedone <- TRUE
            }
          } else {
            if ((!curvedone ) && ((as.numeric(y[popts$ycol]) < (as.numeric(popts$ymin) - bndtol)) ||
                                  (as.numeric(y[popts$ycol]) > (as.numeric(popts$ymax) + bndtol)))) {
              msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
              curvedone <- TRUE
            }
            if ((!curvedone ) && (popts$y2col > 1) &&
                ((as.numeric(y[popts$y2col]) < (as.numeric(popts$y2min) - bndtol)) ||
                 (as.numeric(y[popts$y2col]) > (as.numeric(popts$y2max) + bndtol)))) {
              msg <- "Computation halted:\nMaximum or maximum of y-axis domain reached for 2nd y-axis variable\n"
              curvedone <- TRUE
            }
          }
        } else {
          if ((!curvedone ) && (curvetype == "HP")) {
            if (!any(is.complex(eigval))) {
              msg <- "Computation halted:\nEigenvalues have become real valued\n"
              curvedone <- TRUE
            }
          } else if ((!curvedone ) && (curvetype != "LC")) {
            if ((as.numeric(y[2]) < (as.numeric(popts$ymin) - bndtol)) ||
                (as.numeric(y[2]) > (as.numeric(popts$ymax) + bndtol))) {
              msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
              curvedone <- TRUE
            }
          }
        }
      }

      if (curvedone) {
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cData <- NULL
        break
      }
    } else {                                      # Solution not found
      if (pntnr == 1) {
        msg <- "No convergence at initial point\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        return(NULL)
      }
      if (corrections > as.numeric(nopts$maxiter)) {
        msg <- "Unable to find next solution point with smallest step size\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cData <- NULL
        break
      }
      corrections <- corrections + 1
      cData$stepsize <- cData$stepsize * 0.5
      dyfailed <- cData$guess - cData$yold
      cData$guess <- cData$yold + dyfailed * 0.5
    }
  }
  if (is.null(cData)) {
    if (!is.null(allsols)) {
      if (curvetype == "LC") {
        vals <- lapply((1:statedim), function(i) {
          indxrange <- statedim*(1:(nopts$ninterval*nopts$glorder))
          yname <- names(allsols[1, freeparsdim+i])
          y <- c(allsols[1, freeparsdim+i+indxrange])
          paste0("Min.", yname, "=", round(min(y), 3), ", Max.", yname, "=", round(max(y), 3))
        })
        endPnt <- c("Type" = curvetype,
                    "Description" = paste0(names(allsols[1, 1]), "=", round(allsols[nrow(allsols), 1], 3), ' ',
                                           names(allsols[1, ncol(allsols)]), "=",
                                           round(allsols[nrow(allsols), ncol(allsols)], 3), ' ',
                                           paste(unlist(vals), collapse = ' ')))
      } else {
        endPnt <- c("Type" = curvetype,
                    "Description" = paste(unlist(lapply(1:length(allsols[1,]),
                                                        function(i) {paste0(names(allsols[1,i]), "=",
                                                                            round(allsols[nrow(allsols), i], 3))})),
                                          collapse=' '))
      }
      updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
      endPnt["Description"] <- paste0(sprintf("%04d: ", pntnr-1), endPnt["Description"])

      specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
      specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
      if (curvetype != "LC") specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))
      specialtags <- rbind(specialtags, c(endPnt))
    }
  } else {
    cData$pntnr <- pntnr
  }
  session$userData$curveData <- cData

  if (!is.null(allsols)) {
    return(list("points" = allsols, "eigvals" = alleigs, "tanvec" = alltvs,
                "special.points" = specialsols, "special.eigvals" = specialeigs,
                "special.tanvec" = specialtvs, "special.tags" = specialtags))
  } else {
    return(NULL)
  }
}

HPcontinuation <- function(state, parms, curveData, nopts, rhsval) {

  #  When curveData$freeparsdim == 2 the Jacobian equals the following square (n+2)x(n+2)
  #  matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  # Otherwise, when curveData$freeparsdim == 1 the Jacobian equals the following square
  # (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  A <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

  # And return the determinant of the bialternate product
  twoAI <-bialt2AI(A)

  return(c(unlist(rhsval), det(twoAI)))
}

analyseHP <- function(state, parms, curveData, nopts, session) {

  btval <- HP_BTtest(state, parms, curveData, nopts, NULL)
  names(btval) <- NULL

  lastvals <- curveData$testvals
  testvals <- list()
  testvals$btval <- unlist(btval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    if (!is.null(lastvals$btval) && ((lastvals$btval)*btval < -(nopts$rhstol*nopts$rhstol))) {
      biftype = "BT"
    }
    if (!is.null(biftype)) {
      cData <- curveData
      cData$guess <- NULL
      cData$tanvec <- NULL
      cData$condfun <- list(HPcontinuation, get(paste0("HP_", biftype, "test"), mode = "function"))
      res <- tryCatch(stode(state, time = 0, func = ExtSystemEQ, parms = parms,
                            atol = nopts$rhstol, ctol = nopts$dytol, rtol = nopts$rhstol,
                            maxiter = nopts$maxiter, verbose = FALSE, curveData = cData, nopts = nopts),
                      warning = function(e) {
                        msg <- gsub(".*:", "Warning in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      },
                      error = function(e) {
                        msg <- gsub(".*:", "Error in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      })
      if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
        y <- res$y[1:length(state)]
        names(y) <- names(state)

        # Compute the Jacobian w.r.t. to the free parameter and the state variables
        jac <- jacobian.full(y=y, func=curveData$model, parms=parms, pert = nopts$jacdif)

        # Discard all th rows and columns that pertain to additional variables
        jac <- jac[(1:curveData$pointdim), (1:curveData$pointdim)]

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:curveData$statedim), function(i){paste0("Eigenvalue", i)}))

        # Append the current tangent vector as the last row to the jacobian to
        # preserve direction. See the matcont manual at
        # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
        # Notice that jacobian.full returns a square matrix with NA values on the last row
        jac[nrow(jac),] <- curveData$tanvec[(1:curveData$pointdim)]
        if (rcond(jac) > nopts$rhstol) {
          tvnew <- solve(jac, c(rep(0, (curveData$pointdim-1)), 1))
          tvnorm <- sqrt(sum(tvnew^2))
          tvnew <- tvnew/tvnorm
          names(tvnew) <- unlist(lapply((1:length(state)), function(i){paste0("d", names(state)[i])}))
        } else {
          tvnew <- curveData$tanvec[(1:curveData$pointdim)]
        }

        if (biftype == "BT") testvals$btval <- HP_BTtest(y, parms, curveData, nopts, NULL)

        testvals[["y"]] <- y
        testvals[["tanvec"]] <- tvnew
        testvals[["eigval"]] <- eigval
        testvals[["biftype"]] <- biftype
      } else {
        if (biftype == "BT") msg <- "Locating Bogdanov-Takens point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}

HP_BTtest <- function(state, parms, curveData, nopts, rhsval) {

  res <- TangentVecEQ(state, parms, curveData, nopts)
  jac <- res$jac[(1:curveData$statedim),(2:(curveData$statedim+1))]

  return(c(unlist(rhsval), det(jac)))
}

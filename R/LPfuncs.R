initLP <- function(state, parms, curveData, nopts, session = NULL) {

  # To continue an LP curve we use the extended system of eqs. (10.76) on pg. 504 in
  # Kuznetsov (1998), as the matrix of this system has full rank at a Bogdanov-Takens
  # point. As additional values we need then an approximation to the null vector of
  # the Jacobian q, which we also use as additional vector q0.

  # Compute the Jacobian
  jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system.
  # Because the initial point is BP, curveData$freeparsdim is known to equal 2
  Fx <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

  # Initialize the eigenvalue estimate to 0
  eigval <- 0.0

  # Find the eigenvector of F_x(x, p) pertaining to the eigenvalue with
  # smallest absolute value
  q <- approxNullVec(Fx)

  # Find the eigenvector of (F_x(x, p))^T pertaining to the eigenvalue with
  # smallest absolute value
  p <- approxNullVec(t(Fx))

  return(list(y = c(state, eigval, as.numeric(q), as.numeric(p))))
}

LPcontinuation <- function(state, parms, curveData, nopts, rhsval) {

  # Continuation of a limitpoint is carried out using the defining system
  # (10.97) on page 515 of Kuznetsov (1998):
  #
  #  F(x, p)                 = 0
  #  F_x(x, p) q             = 0
  #  (F_x(x, p))^T p - eps*p = 0
  #  q^T q - 1               = 0
  #  p^T p - 1               = 0
  #
  #  with initial conditions b = 0 and v the eigenvector of the matrix
  #  (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
  #  The unknowns are:
  #
  #  p:   the bifurcation parameter
  #  x:   the solution point
  #  eps: an additional value
  #  q:   the eigenvector of F_x pertaining to the 0 eigenvalue
  #  p:   the eigenvector of (F_x)^T pertaining to the 0 eigenvalue
  #
  # The advantage of this defining system is that detection of Bogdanov-Takens
  # points and cusp point are straightforward

  # This function is only called when curveData$freeparsdim == 2, in which case the
  # Jacobian equals the following square (n+2)x(n+2) matrix of partial
  # derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  A <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

  # Extract the additional variables from the state
  eps = unlist(state[curveData$freeparsdim + curveData$statedim + 1]);
  q = unlist(state[(curveData$freeparsdim + curveData$statedim + 2):(curveData$freeparsdim + curveData$statedim + 1 + curveData$statedim)]);
  p = unlist(state[(curveData$freeparsdim + curveData$statedim + 2 + curveData$statedim):(curveData$freeparsdim + curveData$statedim + 1 + 2*curveData$statedim)]);

  # Add the additional values
  # A q = 0
  rhsval <- c(unlist(rhsval), c(A %*% q))

  # A^T p - eps*p = 0
  rhsval <- c(unlist(rhsval), c((t(A) %*% p) - eps*p))

  # <q, q) - 1 =
  rhsval <- c(unlist(rhsval), c((q %*% q) - 1))

  # <p, p) - 1 =
  rhsval <- c(unlist(rhsval), c((p %*% p) - 1))

  return(unlist(rhsval))
}

analyseLP <- function(state, parms, curveData, nopts, session) {

  btval <- LP_BTtest(state, parms, curveData, nopts, NULL)
  names(btval) <- NULL
  cpval <- LP_CPtest(state, parms, curveData, nopts, NULL)
  names(cpval) <- NULL

  lastvals <- curveData$testvals
  testvals <- list()
  testvals$btval <- unlist(btval)
  testvals$cpval <- unlist(cpval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    if (!is.null(lastvals$cpval) && ((lastvals$cpval)*cpval < -(nopts$rhstol*nopts$rhstol))) {
      biftype = "CP"
    }
    if (!is.null(lastvals$btval) && ((lastvals$btval)*btval < -(nopts$rhstol*nopts$rhstol))) {
      biftype = "BT"
    }
    if (!is.null(biftype)) {
      cData <- curveData
      cData$guess <- NULL
      cData$tanvec <- NULL
      cData$condfun <- list(LPcontinuation, get(paste0("LP_", biftype, "test"), mode = "function"))
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

        if (biftype == "BT") testvals$btval <- LP_BTtest(y, parms, curveData, nopts, NULL)
        else testvals$cpval <- LP_CPtest(y, parms, curveData, nopts, NULL)

        testvals[["y"]] <- y
        testvals[["tanvec"]] <- tvnew
        testvals[["eigval"]] <- eigval
        testvals[["biftype"]] <- biftype
      } else {
        if (biftype == "BT") msg <- "Locating Bogdanov-Takens point failed\n"
        else msg <- "Locating cusp point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}

LP_BTtest <- function(state, parms, curveData, nopts, rhsval) {
  # Extract the additional variables from the state
  q = unlist(state[(curveData$pointdim + 2):(curveData$pointdim + 1 + curveData$statedim)]);
  p = unlist(state[(curveData$pointdim + 2 + curveData$statedim):(curveData$pointdim + 1 + 2*curveData$statedim)]);

  # See pg 515 of Kuznetsov (1998), just below eq. (10.97)

  return(c(unlist(rhsval), as.numeric(p %*% q)))
}

LP_CPtest <- function(state, parms, curveData, nopts, rhsval) {

  # Extract the additional variables from the state
  qval = unlist(state[(curveData$pointdim + 2):(curveData$pointdim + 1 + curveData$statedim)]);
  pval = unlist(state[(curveData$pointdim + 2 + curveData$statedim):(curveData$pointdim + 1 + 2*curveData$statedim)]);

  # The quantity B(q, q) can be computed most easily using a directional derivative
  # see eq. (10.52) and its approximation some lines below (10.52) on page 490 of
  # Kuznetsov (1998)
  stmp <- state[1:curveData$pointdim]
  hv <- rep(0, length(stmp))
  hv[(curveData$freeparsdim+1):curveData$pointdim] <- sqrt(nopts$jacdif)*qval
  rhsP <- unlist(curveData$model(0, (stmp + hv), parms))
  rhsM <- unlist(curveData$model(0, (stmp - hv), parms))
  Bqq <- (rhsP + rhsM)/nopts$jacdif

  # See bottom of pg 514 of Kuznetsov (1998), just above eq. (10.97)

  return(c(unlist(rhsval), as.numeric(pval %*% Bqq)))
}

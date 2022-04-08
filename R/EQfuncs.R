initEQ <- function(state, parms, curveData, nopts, session = NULL) {
  if (curveData$inittype == "BP") {
    #  The Jacobian equals the following square (n+1)x(n+1) matrix of partial
    #  derivatives:
    #
    #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
    #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
    #           |   .       .    ...    .   |
    #      Df = |   .       .    ...    .   |
    #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
    #           |   NA      NA   ...    NA  |
    #
    # Here is state the vector of state variables (length n)

    jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)
    # Append the current tangent vector as the last row to the jacobian to
    # obtain the matrix D(0) (see eq.(10.61)-(10.63) on pg. 499 of Kuznetsov
    # (1996))
    jac[nrow(jac),] <- curveData$tanvec
    D <- jac

    # Approximate the vector q1 with the stored tangent vector
    q1 <- as.numeric(c(curveData$tanvec))
    names(q1) <- names(state)

    # Solve the vector q2 from D q2 = 0
    eig <- eigen(D)
    # minindx <- which.min(abs(Re(eig$values)))
    # q2 <- as.numeric(c(eig$vectors[,minindx]))
    allindx <- (1:length(eig$values))[Im(eig$values) == 0]
    minindx <- allindx[which.min(abs(Re(eig$values[Im(eig$values) == 0])))]
    q2 <- as.numeric(Re(c(eig$vectors[,minindx])))
    names(q2) <- names(state)

    # Extract the (n+1)xn matrix (J(0))^T (see eq.(10.61)-(10.63) on pg. 499 of
    # Kuznetsov (1996))
    JT <- t(jac[(1:(nrow(jac)-curveData$freeparsdim)), (1:ncol(jac))])

    # Solve the vector phi from (J(0))^T phi = 0. According to Kuznetsov (1996,
    # pg. 497, eq. (10.59)) there is a unique vector (up to scalar multiple) phi
    # of length n satisfying this linear equation. To handle the arbitrary
    # scaling vector, drop that row of the matrix which yields the least
    # singular restricted matrix, such that eigenvaues can be determined easily
    # Estimate the condition numbers
    if (nrow(JT) == 1) {
      phi <- 1.0
    } else {
      condnrs <- unlist(lapply((1:nrow(JT)),
                               function(i){
                                 incrows <- rep(TRUE, nrow(JT))
                                 incrows[i] <- FALSE
                                 jacc <- JT[incrows,]
                                 if (abs(det(jacc)) > 1.0E-5) return(rcond(jacc))
                                 else return(0.0)
                               }))
      maxind <- which.max(condnrs)
      incrows <- rep(TRUE, nrow(JT))
      incrows[maxind] <- FALSE

      # For the restricted matrix solve for the eigenvalues and select the
      # eigenvector belonging to the eigenvalue closest to 0
      eig <- eigen(JT[incrows,])
      # minindx <- which.min(abs(Re(eig$values)))
      # phi <- c(eig$vectors[,minindx])
      allindx <- (1:length(eig$values))[Im(eig$values) == 0]
      minindx <- allindx[which.min(abs(Re(eig$values[Im(eig$values) == 0])))]
      phi <- Re(c(eig$vectors[,minindx]))
    }

    # Now compute B(q1,q2) and B(q2, q2) as explained. The function B(x,y) is
    # defined in eq. (10.56) on pg. 496 of Kuznetsov (1996). The vectors
    # B(q1,q2) and B(q2, q2) of length n (n equals the number of state
    # variables) are needed to compute the factors b12 and b22 (see below)

    # The quantity B(q2, q2) can be computed most easily using a directional derivative
    # see eq. (10.52) and its approximation some lines below (10.52) on page 490 of
    # Kuznetsov (1996)
    rhsP1 <- unlist(curveData$model(0, state + sqrt(nopts$jacdif)*(q1 + q2), parms))
    rhsM1 <- unlist(curveData$model(0, state - sqrt(nopts$jacdif)*(q1 + q2), parms))
    rhsP2 <- unlist(curveData$model(0, state + sqrt(nopts$jacdif)*(q1 - q2), parms))
    rhsM2 <- unlist(curveData$model(0, state - sqrt(nopts$jacdif)*(q1 - q2), parms))

    Bq1q2 <- (rhsP1 + rhsM1 - rhsP2 - rhsM2)/(4*nopts$jacdif)

    rhsP <- unlist(curveData$model(0, state + sqrt(nopts$jacdif)*q2, parms))
    rhsM <- unlist(curveData$model(0, state - sqrt(nopts$jacdif)*q2, parms))
    Bq2q2 <- (rhsP + rhsM)/nopts$jacdif

    # The inner product of phi with the vectors B(q1,q2) and B(q2, q2) of length
    # n (n equals the number of state variables) yields the factors b12 and b22
    b12 <- c(phi %*% Bq1q2)
    b22 <- c(phi %*% Bq2q2)
    v2 <- (-b22/(2*b12))*q1 + q2
    v2 <- v2/sqrt(sum(v2^2))

    ############## Determine the new step along the curve
    guess <- state + curveData$stepsize * v2

    return(list(y = guess, tanvec = v2))
  } else return(list())
}

ExtSystemEQ <- function(t, state, parms, curveData, nopts = NULL) {
  # Evaluate the equations
  rhsval <- unlist(curveData$model(t, state, parms))
  names(rhsval) <- NULL

  # add possible conditions for curveData$freeparsdim > 1 continuation
  if (!is.null(curveData$condfun)) {
    for (i in (1:length(curveData$condfun))) {
      rhsval <- do.call(curveData$condfun[[i]], list(state, parms, curveData, nopts = nopts, rhsval))
    }
    names(rhsval) <- NULL
  }

  # add the dot product of the difference between the current state and
  # the initial guess with tangent vector for pseudo-arclength continuation
  if (!is.null(curveData$guess) && !is.null(curveData$tanvec)) {
    rhsval <- c(unlist(rhsval), c((state - curveData$guess) %*% curveData$tanvec))
  }

  return(list(rhsval))
}

TangentVecEQ <- function(state, parms, curveData, nopts) {

  # Routine calculates the tangent vector to the curve determined by the equation
  #
  #                           F(x, p1) = 0
  #
  # This is a system of n equations, in which n is the length of the state variable
  # vector (given by curveData$statedim)

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

  # If curveData$freeparsdim == 2 delete the last row and 2nd column of the Jacobian, the
  # latter being the derivative w.r.t. second free parameters
  jac <- jac[(1:(curveData$statedim+1)), c(1, ((curveData$freeparsdim+1):curveData$pointdim))]

  tvdone <- FALSE
  # Append the current tangent vector as the last row to the jacobian to
  # preserve direction. See the matcont manual at
  # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  # But only do this if the tangent vector passed as argument has the same
  # dimension as the equilibrium curve
  if (!is.null(curveData$tanvec) && (length(curveData$tanvec) == (curveData$statedim + 1))) {
    jac[nrow(jac),] <- curveData$tanvec
    if (rcond(jac) > nopts$rhstol) {
      # Solve for the tangent vector to the curve of the 1st free parameter and
      # the state variables.
      tvnew <- solve(jac, c(rep(0, (length(curveData$tanvec)-1)), 1))
      tvdone <- TRUE
    }
  }

  if (!tvdone) {
    # Set the last row with NA to 0
    jac[nrow(jac),] <- 0

    # Estimate the condition numbers
    condnrs <- unlist(lapply((1:ncol(jac)),
                             function(i){
                               jacc <- jac;
                               jacc[nrow(jac),i] <- 1;
                               if (abs(det(jacc)) > 1.0E-5) return(rcond(jacc))
                               else return(0.0)
                             }))
    maxind <- which.max(condnrs)
    jac[nrow(jac),maxind] <- 1

    # Solve for the tangent vector to the curve of the 1st free parameter and
    # the state variables.
    tvnew <- solve(jac, c(rep(0, (nrow(jac)-1)), 1))
  }
  tvnorm <- sqrt(sum(tvnew^2))
  tvnew <- tvnew/tvnorm
  names(tvnew) <- NULL
  jac[nrow(jac),] <- NA

  return (list(jacobian = jac, tanvec = tvnew))
}

analyseEQ <- function(state, parms, curveData, nopts, session) {

  bpval <- EQ_BPtest(state, parms, curveData, nopts, NULL)
  names(bpval) <- NULL
  if (length(state) > 2) {
    hpval <- EQ_HPtest(state, parms, curveData, nopts, NULL)
    names(hpval) <- NULL
  } else {
    hpval <- 1.0
  }
  lpval <- EQ_LPtest(state, parms, curveData, nopts, NULL)
  names(lpval) <- NULL

  lastvals <- curveData$testvals
  testvals <- list()
  testvals$bpval <- unlist(bpval)
  testvals$hpval <- unlist(hpval)
  testvals$lpval <- unlist(lpval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    # if (!is.null(lastvals$hpval) && ((lastvals$hpval)*hpval < 0.0)) {
    if (!is.null(lastvals$hpval) && ((lastvals$hpval)*hpval < -(nopts$rhstol*nopts$rhstol))) {
        biftype = "HP"
    }

    # Because the test function for a limit point involves the condition that the test function for
    # a branching point should be non-zero (See below eq. (35)-(37) in the MATCONT manual from 2012
    # (ManualSep2012.pdf), we first test for a branching point and only after that test for a
    # limit point

    # if (!is.null(lastvals$bpval) && ((lastvals$bpval)*bpval < 0.0)) {
    if (!is.null(lastvals$bpval) && ((lastvals$bpval)*bpval < -(nopts$rhstol*nopts$rhstol))) {
        biftype = "BP"
    } else {
      # if (!is.null(lastvals$lpval) && ((lastvals$lpval)*lpval < 0.0)) {
      if (!is.null(lastvals$lpval) && ((lastvals$lpval)*lpval < -(nopts$rhstol*nopts$rhstol))) {
          biftype = "LP"
      }
    }
    if (!is.null(biftype)) {
      cData <- curveData
      cData$guess <- NULL
      cData$tanvec <- NULL
      cData$condfun <- list(get(paste0("EQ_", biftype, "test"), mode = "function"))
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

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:curveData$statedim), function(i){paste0("Eigenvalue", i)}))

        NeutralSaddle <- (biftype == "HP") && (abs(Im(eigval[1])) < as.numeric(nopts$iszero))
        if (!NeutralSaddle) {
          # Append the current tangent vector as the last row to the jacobian to
          # preserve direction. See the matcont manual at
          # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
          # Notice that jacobian.full returns a square matrix with NA values on the last row
          jac[nrow(jac),] <- curveData$tanvec
          if (rcond(jac) > nopts$rhstol) {
            tvnew <- solve(jac, c(rep(0, (length(curveData$tanvec)-1)), 1))
            tvnorm <- sqrt(sum(tvnew^2))
            tvnew <- tvnew/tvnorm
            names(tvnew) <- unlist(lapply((1:length(state)), function(i){paste0("d", names(state)[i])}))
          } else {
            tvnew <- curveData$tanvec
          }

          if (biftype == "BP") testvals$bpval <- EQ_BPtest(y, parms, curveData, nopts, NULL)
          else if (biftype == "HP") testvals$hpval <- EQ_HPtest(y, parms, curveData, nopts, NULL)
          else testvals$lpval <- EQ_LPtest(y, parms, curveData, nopts, NULL)

          testvals[["y"]] <- y
          testvals[["tanvec"]] <- tvnew
          testvals[["eigval"]] <- eigval
          testvals[["biftype"]] <- biftype
        }
      } else {
        if (biftype == "BP") msg <- "Locating branching point failed\n"
        else if (biftype == "HP") msg <- "Locating Hopf bifurcation point failed\n"
        else msg <- "Locating limit point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}

EQ_BPtest <- function(state, parms, curveData, nopts, rhsval) {

  # The following defines as test function the determinant of the extended Jacobian of the
  # system:
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA   ...    NA  |
  #
  # with the last row replaced by the tangent vector to the curve. See eq. (35) in the MATCONT
  # manual from 2012 (ManualSep2012.pdf)

  res <- TangentVecEQ(state, parms, curveData, nopts)
  jac <- res$jacobian
  jac[nrow(jac),] <- res$tanvec

  return(c(unlist(rhsval), det(jac)))
}

EQ_HPtest <- function(state, parms, curveData, nopts, rhsval) {

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

EQ_LPtest <- function(state, parms, curveData, nopts, rhsval) {

  res <- TangentVecEQ(state, parms, curveData, nopts)

  # The following defines as test function the parameter component of the tangent vector to the
  # solution curve. See the remark below eq. (10.39) in Kuznetsov (1998) and eq. (37) in the MATCONT
  # manual from 2012 (ManualSep2012.pdf)
  return(c(unlist(rhsval), res$tanvec[1]))

  # The following defines as test function the determinant of the restricted Jacobian f_x(x, p) which
  # equals the produce of the eigenvalues of this Jacobian. See eq. (10.39) in Kuznetsov (1998)
  # return(c(unlist(rhsval), det(res$jac[(1:curveData$statedim),((curveData$freeparsdim+1):curveData$pointdim)])))
}

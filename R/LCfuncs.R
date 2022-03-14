glPoints <- function(n) {
  #  For 2 <= N <= 7
  #  N point Gauss-Legendre quadrature rule over the interval [-1,1].
  p <- switch (as.numeric(n)-1,
               {c <- 1/sqrt(3); c(-c, c)},
               {c <- sqrt(6/10); c(-c, 0, c)},
               {c <- sqrt(24/245);  c1 <- sqrt(3/7+c); c2 <- sqrt(3/7-c); c(-c1, -c2, c2, c1)},
               {c1 <- 0.90617984593866399280; c2 <- 0.53846931010568309104; c(-c1, -c2, 0, c2, c1)},
               {c1 <- 0.93246951420315202781;  c2 <- 0.66120938646626451366; c3 <- 0.23861918608319690863; c(-c1, -c2, -c3, c3, c2, c1)},
               {c1 <- 0.949107991234275852452; c2 <- 0.74153118559939443986; c3 <- 0.40584515137739716690; c(-c1, -c2, -c3, 0, c3, c2, c1)}
  )
  return(p)
}

glWeights <- function(n) {
  # For 2 <= N <= 7
  # N point Gauss-Legendre quadrature rule over the interval [-1,1].
  p <- switch (as.numeric(n)-1,
               c(1, 1),
               c(0.55555555555555555555, 0.88888888888888888888, 0.55555555555555555555),
               c(0.34785484513745385737, 0.65214515486254614263, 0.65214515486254614263, 0.34785484513745385737),
               c(0.23692688505618908751, 0.47862867049936646804, 0.56888888888888888888, 0.47862867049936646804,
                 0.23692688505618908751),
               c(0.17132449237917034504, 0.36076157304813860757, 0.46791393457269104739, 0.46791393457269104739,
                 0.36076157304813860757, 0.17132449237917034504),
               c(0.12948496616886969327, 0.27970539148927666790, 0.38183005050511894495, 0.41795918367346938775,
                 0.38183005050511894495, 0.27970539148927666790, 0.12948496616886969327)
  )
  return(p)
}

ncWeights <- function(n) {
  #  For 2 <= N <= 7, returns the weights W of an
  #  N+1 point Newton-Cotes quadrature rule.
  p <- switch (as.numeric(n)-1,
               c(1/6, 4/6, 1/6),
               c(1/8, 3/8, 3/8, 1/8),
               c(7/90, 32/90, 12/90, 32/90, 7/90),
               c(19/288, 25/96, 25/144, 25/144, 25/96, 19/288),
               c(41/840, 9/35, 9/280, 34/105, 9/280, 9/35, 41/840),
               c(751/17280, 3577/17280, 49/640, 2989/17280, 2989/17280, 49/640, 3577/17280, 751/17280)
  )
  return(p)
}

fastkron <- function(c,p,A,B) {
  t = p:((c+2)*p-1)
  return(A[rep(1,p), floor(t/p)]*B[, (t%%p)+1])
}

initLC <- function(state, parms, curveData, nopts, session = NULL) {
  if (curveData$inittype == "HP") {
    cData <- curveData

    # Compute the Jacobian
    jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)

    # Extract the restricted Jacobian of the system.
    Fx <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

    # Solve for the eigenvalues and eigenvectors
    eig <- eigen(Fx)

    # Get the index ordered such that the eigenvalues closest
    orderedindx <- order(abs(Re(eig$values)))
    eigvals <- eig$values[orderedindx]
    eigvecs <- eig$vectors[,orderedindx]

    if ((as.numeric(Im(eigvals[1])) == 0) && (as.numeric(Im(eigvals[2])) == 0)) {
      msg <- paste0("Computation aborted:\nStarting point is not a Hopf bifurcation point\n")
      if (!is.null(session)) updateConsoleLog(session, msg)
      else cat(msg)
    }

    # Get imaginary part and corresponding eigenvector
    omega <- abs(Im(eigvals[1]))
    Qv <- eigvecs[,1]

    d <- t(Re(Qv)) %*% Re(Qv)
    s <- t(Im(Qv)) %*% Im(Qv)
    r <- t(Re(Qv)) %*% Im(Qv)

    Qv <- Qv %*% exp((1i)*atan2(2*r,s-d)/2)
    Qv <- Qv/norm(Re(Qv), type = "2")

    cData$mesh <- (1.0/nopts$ninterval)*(0:nopts$ninterval)
    cData$finemesh <- (1.0/(nopts$ninterval*nopts$glorder))*(0:(nopts$ninterval*nopts$glorder))
    cData$finemeshdim <- nopts$ninterval*nopts$glorder+1
    t <- t(kronecker(exp(2*pi*(1i)*cData$finemesh), t(Qv)))

    cData$upoldp <- -Im(t)
    v0 <- c(0, Re(c(t)), 0)

    y <- c(state[1], rep(state[2:cData$pointdim], length(cData$finemesh)), 2*pi/omega) + nopts$lcampl*v0
    varnames <- c(names(state[1]), rep(names(state[2:cData$pointdim]), length(cData$finemesh)), "LCperiod")
    names(y) <- varnames
    v0 <- v0/norm(v0, type = "2")

    cData$freeparsdim <- 1
    cData$pointdim <- length(y)
    cData$varnames <- varnames
    cData$eignames <- NULL
    cData$tvnames <- unlist(lapply((1:length(y)), function(i){paste0("d", varnames[i])}))

    # Calculate weights
    zm <- glPoints(nopts$glorder)/2 + 0.5
    xm <- (0:nopts$glorder)/nopts$glorder

    xindx <- (1:(nopts$glorder+1))
    wghts <- matrix(0, nopts$glorder+1, nopts$glorder)
    wpvec <- wghts

    for (j in xindx) {
      xmi <- xm[setdiff(xindx,j)]
      denom <- prod(xm[j]-xmi)
      zxmat <- matrix(zm, nopts$glorder, length(zm), byrow = T) - matrix(xmi, nopts$glorder, length(xmi), byrow = F)
      wghts[j,] <- unlist(lapply((1:ncol(zxmat)), function(i){prod(zxmat[,i])}))/denom
      sval <- rep(0, nopts$glorder)
      for (l in setdiff(xindx, j)){
        xmil <- xm[setdiff(xindx,c(j, l))]
        zxmat <- matrix(zm, nopts$glorder-1, length(zm), byrow = T) - matrix(xmil, nopts$glorder-1, length(xmi), byrow = F)
        if (nopts$glorder > 2) sval <- sval + unlist(lapply((1:ncol(zxmat)), function(i){prod(zxmat[,i])}))
        else sval <- sval + zxmat
      }
      wpvec[j,] <- sval/denom
    }

    cData$wi  <- ncWeights(nopts$glorder)
    cData$pwi <- matrix(cData$wi, nrow = cData$statedim, ncol = length(cData$wi), byrow = T)
    cData$wt  <- wghts
    cData$wpvec <- wpvec
    cData$wp <- kronecker(t(wpvec), diag(cData$statedim))

    return(list(y = y, curveData = cData, tanvec = v0))
  } else {
    msg <- paste0("Computation aborted:\nLimit cycle continuation starting from ", curveData$inittype, " not implemented\n")
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
    return(NULL)
  }
}

ExtSystemLC <- function(t, state, parms, curveData, nopts = NULL) {

  # Extract a matrix with state variable values. The matrix has statedim rows and a number of colums equal
  # to the finemesh dimension
  state0 <- state[c(1:(curveData$freeparsdim+curveData$statedim),curveData$pointdim)]
  ups <- matrix(state[curveData$freeparsdim + (1:(curveData$statedim*curveData$finemeshdim))],
                curveData$statedim, curveData$finemeshdim, byrow = F)

  # Setup a copy for the derivatives, filled with 0
  rhsval <- matrix(0, curveData$statedim, curveData$finemeshdim)

  # Compute the values for the integral condition and setup a matrix for the integral condition, filled with 0
  ficd <- colSums(ups * curveData$upoldp)
  ficdmat <- matrix(0, nopts$glorder+1, nopts$ninterval)

  # Evaluate the derivatives at all the nodal points of the mesh
  dt <- 1.0/nopts$ninterval
  range1 <- (1:(nopts$glorder+1))
  for (i in (1:nopts$ninterval)) {
    # Value of polynomial on each collocation point
    xp <- ups[, range1] %*% curveData$wt

    # Derivative of polynomial on each collocation point
    t  <- (ups[, range1] %*% curveData$wpvec)/dt

    # Evaluate function value on each collocation point
    for (j in (1:nopts$glorder)) {
      sval <- state0
      sval[curveData$freeparsdim + (1:curveData$statedim)] <- c(xp[,j])
      rhsval[, j + (i-1)*nopts$glorder] <- c(t[,j] - state["LCperiod"]*unlist(curveData$model(t, sval, parms)))
    }

    # Put the appropriate values of the integral condition in the matrix
    ficdmat[,i] <- ficd[range1]

    range1 <- range1 + nopts$glorder
  }

  # Ciruclar boundary conditions
  rhsval[, curveData$finemeshdim] <- ups[, 1] - ups[, curveData$finemeshdim]

  # Integral constraint
  rhsval <- c(c(rhsval), sum(dt*(curveData$wi %*% ficdmat)))

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

ExtSystemLCblockjac <- function(xp, state, parms, curveData, nopts = NULL) {

  jaccol <- curveData$statedim*(nopts$glorder+1)+curveData$freeparsdim+1
  jacrow <- curveData$statedim*nopts$glorder
  state0 <- state[names(state) != "LCperiod"]
  blockjac <- matrix(0, jacrow, jaccol)

  dt <- 1.0/nopts$ninterval
  wploc = curveData$wp/dt;

  # xp:value of polynomial on each collocation point
  range1 <- (1:curveData$statedim)
  for (j in (1:nopts$glorder)) {
    sval <- state0
    sval[curveData$freeparsdim + (1:curveData$statedim)] <- c(xp[,j])
    jac <- jacobian.full(y=sval, func=curveData$model, parms=parms, pert = nopts$jacdif)
    sysjac <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):(curveData$freeparsdim+curveData$statedim))]

    blockjac[range1,1] <- -state["LCperiod"]*jac[(1:curveData$statedim), 1]
    blockjac[range1,1+(1:(jaccol-2))] <- wploc[range1,] - state["LCperiod"]*fastkron(nopts$glorder, curveData$statedim, t(curveData$wt[,j]), sysjac)
    blockjac[range1,jaccol] <- -t(unlist(curveData$model(0, sval, parms)))
    range1 <- range1 + curveData$statedim
  }

  return(blockjac)
}

ExtSystemLCjac <- function(y, func, parms, pert = NULL, curveData, nopts = NULL) {

  # Extract a matrix with state variable values. The matrix has statedim rows and a number of colums equal
  # to the finemesh dimension
  state0 <- y[c(1:(curveData$freeparsdim+curveData$statedim),curveData$pointdim)]
  ups <- matrix(y[curveData$freeparsdim + (1:(curveData$statedim*curveData$finemeshdim))],
                curveData$statedim, curveData$finemeshdim, byrow = F)

  # Setup a copy for the jacobian matrix, filled with 0
  jacdim <- curveData$statedim*curveData$finemeshdim+curveData$freeparsdim+1
  fulljac <- matrix(0, jacdim, jacdim)
  dt <- 1.0/nopts$ninterval

  # Evaluate the block jacobians
  blockrow <- curveData$statedim*nopts$glorder
  blockcol <- curveData$statedim*(nopts$glorder+1)
  rowrange <- (1:blockrow)
  colrange <- (1:blockcol)
  range1 <- (1:(nopts$glorder+1))
  range2 <- (1:blockcol)
  ic <- rep(0, jacdim);

  for (i in (1:nopts$ninterval)) {
    # Value of polynomial on each collocation point
    xp <- ups[, range1] %*% curveData$wt

    partjac <- ExtSystemLCblockjac(xp, state0, parms, curveData, nopts)
    fulljac[rowrange, 1] <- partjac[(1:blockrow), 1]
    fulljac[rowrange, 1 + colrange] <- partjac[(1:blockrow), 1 + (1:blockcol)]
    fulljac[rowrange, jacdim] <- partjac[(1:blockrow), blockcol+2]

    # Derivative of the integral constraint
    p <- dt*(curveData$upoldp[,range1]*curveData$pwi)
    ic[1 + colrange] <- ic[1 + colrange] + p[1:blockcol]

    range1 <- range1 + nopts$glorder
    rowrange <- rowrange + blockrow
    colrange <- colrange + blockrow
  }

  # Derivative of the boundary condition
  bc <- cbind(matrix(0, curveData$statedim, 1), diag(curveData$statedim), matrix(0, curveData$statedim, curveData$statedim*(nopts$ninterval*nopts$glorder -1)),
              -diag(curveData$statedim), matrix(0, curveData$statedim, 1))
  rowrange <- (jacdim - (curveData$statedim + 2)) + (1:curveData$statedim)
  fulljac[rowrange,] <- bc

  # Derivative of the integral constraint
  fulljac[jacdim-1,] <- ic

  return(fulljac)
}

updateRefSol <- function(t, state, parms, curveData, nopts = NULL) {

  ##############################################################
  # Extract a matrix with state variable values. The matrix has statedim rows and a number of colums equal
  # to the finemesh dimension
  state0 <- state[c(1:(curveData$freeparsdim+curveData$statedim),curveData$pointdim)]
  ups <- matrix(state[curveData$freeparsdim + (1:(curveData$statedim*curveData$finemeshdim))],
                curveData$statedim, curveData$finemeshdim, byrow = F)

  # Setup a copy for the derivatives, filled with 0
  rhsval <- matrix(0, curveData$statedim, curveData$finemeshdim)

  # Evaluate the derivatives at all the nodal points of the mesh
  for (j in (1:curveData$finemeshdim)) {
    # Evaluate function value on each collocation point
    sval <- state0
    sval[curveData$freeparsdim + (1:curveData$statedim)] <- c(ups[,j])
    rhsval[,j] <- state["LCperiod"]*unlist(curveData$model(t, sval, parms))
  }

  return(rhsval)
}


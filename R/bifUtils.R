bialt2AI <- function(A) {
  # Constructs the bialternate matrix product of 2A and I.
  #
  # See pages 485-487 in Kuznetsov, 1995; Elements of Applied Bifurcation
  # Analysis In particular, figure 10.8 and the equation at the top of page 487.
  n <- ncol(A)
  m <-n*(n-1)/2

  # p = 2, 3, 3, 4, 4, 4, ....
  p <- unlist(lapply((2:n), function(i) rep(i, (i-1))))
  # q = 1, 1, 2, 1, 2, 3, ....
  q <- unlist(lapply((2:n), function(i) (1:(i-1))))
  r <- p
  s <- q
  P <- matrix(p, m, m)
  Q <- matrix(q, m, m)
  R <- t(P)
  S <- t(Q)
  Aps <- matrix(c(A[p,]), length(p), n)[,s]
  Apr <- matrix(c(A[p,]), length(p), n)[,r]
  App <- matrix(c(A[p,]), length(p), n)[,p]
  Aqq <- matrix(c(A[q,]), length(q), n)[,q]
  Aqs <- matrix(c(A[q,]), length(q), n)[,s]
  Aqr <- matrix(c(A[q,]), length(q), n)[,r]

  twoAI <- (R == Q)*(-Aps) + ((R != P) & (S == Q))*Apr + ((R == P) & (S == Q))*(App + Aqq) + ((R == P) & (S != Q))*Aqs + (S == P)*(-Aqr)

  return(twoAI)
}

approxNullVec <- function(A) {
  # Find the eigenvector of the matrix A pertaining to the eigenvalue with
  # smallest absolute value
  eig <- eigen(A)
  minindx <- which.min(abs(Re(eig$values)))
  eigvec <- c(eig$vectors[,minindx])
  return (as.numeric(eigvec))
}

rcprintf <- function(fmt, x) {
  # Prints x to string using sprintf(fmt,...), handling both real and complex numbers
  if (is.complex(x) && (abs(Im(x)) > 1.0E-99)) {
    if (Im(x) > 0) return(sprintf("%12.5E+%11.5Ei", ifelse(abs(Re(x)) > 1.0E-99, Re(x), 0), abs(Im(x))))
    else return(sprintf("%12.5E-%11.5Ei", ifelse(abs(Re(x)) > 1.0E-99, Re(x), 0), abs(Im(x))))
  } else return(sprintf("%12.5E", ifelse(abs(x) > 1.0E-99, x, 0)))
}

converty2y <- function(y, ymin1, ymax1, logy1, ymin2, ymax2, logy2) {
  # Convert values on the second y axis to corresponding values on the first y-axis
  if (logy2 == 0) {
    ytmp <- pmax((y - ymin2)/(ymax2 - ymin2), -0.5)
  } else {
    ytmp <- pmax(log(pmax(y, 1.0E-15)/ymin2)/(log(ymax2/ymin2)), -0.5)
  }

  if (logy1 == 0) {
    yval <- ymin1 + ytmp*(ymax1 - ymin1)
  } else {
    yval <- exp(log(ymin1) + ytmp*(log(ymax1/ymin1)))
  }
  return(yval)
}

setStepSize <- function(y, cData, minstep, iszero) {
  minstepsize <- as.numeric(minstep)

  ############## Determine the default step along the curve
  dydef <- as.numeric(cData$stepsize) * cData$tanvec

  # Determine the relative change in the components and the index of the largest change
  yabs <- abs(y)
  indx0s <- (yabs < as.numeric(iszero))                     # Indices of zero elements of y. Ignore their relative change
  yabs[indx0s] <- 1.0
  dyrel <- abs(dydef) / yabs
  dyrel[indx0s] <- 0
  if (length(dyrel) > cData$pointdim) dyrel[((cData$pointdim+1):length(dyrel))] <- 0
  dyind <- which.max(dyrel)                                 # Index with maximum relative change

  # Scale the step in such a way that the change in the fastest changing component equals 1
  dy <- dydef / abs(dydef[dyind])

  # Compute the target step size as a relative change equal to as.numeric(cData$stepsize) in the fastest changing y component
  dyfinal <- max(abs(as.numeric(cData$stepsize) * y[dyind]), minstepsize) * dy

  # If there are invalid entries, fall back to the default step along the curve
  if (any(is.infinite(dyfinal)) || any(is.na(dyfinal))) dyfinal <- dydef

  # If all components change less than the absolute minimum step size times the tangent fall back to that absolute minimum step
  if (all(abs(dyfinal) < abs(minstepsize * cData$tanvec))) {
    dyfinal <- as.numeric(cData$direction) * minstepsize * cData$tanvec
  }

  return(dyfinal)
}

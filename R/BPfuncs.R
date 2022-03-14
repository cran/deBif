initBP <- function(state, parms, curveData, nopts, session = NULL) {

  # Compute the Jacobian
  jac <- jacobian.full(y=state, func=curveData$model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system.
  # Because the initial point is BP, curveData$freeparsdim is known to equal 2
  Fx <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

  # Initialize the eigenvalue estimate to 0
  eigval <- 0.0

  # Find the eigenvector of (F_x(x, p))^T pertaining to the eigenvalue with
  # smallest absolute value
  eigvec <- approxNullVec(t(Fx))

  return(list(y = c(state, eigval, eigvec), tanvec = curveData$tanvec))
}

BPcontinuation <- function(state, parms, curveData, nopts, rhsval) {

  # Continuation of a branching point is carried out using the defining system
  # eq. 41 on page 36 of the Matcont documentation (august 2011):
  #
  #  F(x, p) + b*v   = 0
  #  (F_x(x, p))^T v = 0
  #  v^T F_p(x, p)   = 0
  #  v^T v - 1       = 0
  #
  #  with initial conditions b = 0 and v the eigenvector of the matrix
  #  (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
  #  The unknowns are:
  #
  #  p:  the bifurcation parameter
  #  x:  the solution point
  #  b:  an additional value
  #  v:  the eigenvector (same dimension as x)
  #
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
  jac <- jacobian.full(y=state[1:curveData$pointdim], func=curveData$model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  Fx <- jac[(1:curveData$statedim), ((curveData$freeparsdim+1):curveData$pointdim)]

  # Extract F_p(x, p)
  Fp <- jac[(1:curveData$statedim), 1]

  # Extract the value of b and v from the state
  eigval = unlist(state[curveData$pointdim + 1]);
  eigvec = unlist(state[(curveData$pointdim + 2):length(state)]);

  #  F(x, p) + b*v   = 0
  rhsval <- unlist(rhsval)
  rhsval[1:length(eigvec)] <- rhsval[1:length(eigvec)] + eigval*eigvec

  # (F_x(x, p))^T v = 0
  rhsval <- c(rhsval, c(t(Fx) %*% eigvec))

  # v^T F_p(x, p)   = 0
  rhsval <- c(rhsval, c(eigvec %*% Fp))

  # v^T v - 1       = 0
  rhsval <- c(rhsval, (c(eigvec %*% eigvec) - 1))

  return(rhsval)
}

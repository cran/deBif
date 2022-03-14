# The juvenile-adult cannibalism model
# ------------------------------------

# Equations:
# ----------
#
# dJ
# -- = alpha A - delta J - beta J A
# dt
#
# dA                  mu A
# -- = delta J -  ------------
# dt              rho + beta J
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(J = 5.0, A = 100.0)

# Parameters have to be specified as a named vector of parameters.
parms <- c(alpha = 1.0, beta = 0.2, delta = 0.2, mu = 0.01, rho = 0.5)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
cannibalism <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    dJ = alpha*A - delta*J - beta*J*A
    dA = delta*J - mu*A/(rho + beta*J)

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dJ, dA)))
  })
}

bifurcation(cannibalism, state, parms)

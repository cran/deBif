# The Lotka-Volterra predator-prey model with logistic prey growth
# ----------------------------------------------------------------

# Equations:
# ----------
#
# dR            R
# -- = r R (1 - -) - a R N
# dt            K
#
# dN
# -- = c a R N - delta N
# dt
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(R=1, N=0.01)

# Parameters have to be specified as a named vector of parameters.
parms <- c(r=1, K=1, a=1, c=1, delta=0.5)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dR, dN)))
  })
}

phaseplane(model, state, parms)


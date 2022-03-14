# The logistic growth model
# -------------------------

# Equations:
# ----------
#
# dR            R
# -- = r R (1 - -)
# dt            K
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(R=0.01)

# Parameters have to be specified as a named vector of parameters.
parms <- c(r=1, K=1)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR <- r*R*(1 - R/K)

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dR)))
  })
}

phaseplane(model, state, parms)

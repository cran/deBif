# The spruce budworm model
# ------------------------

# Equations:
# ----------
#
# dN             N            N^2
# -- = r N (1 - ---) - E ------------- P
# dt            q A      (f A)^2 + N^2
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(N=1.0)

# Parameters have to be specified as a named vector of parameters.
parms <- c(A = 0.5, q = 20, E = 0.314, f = 0.474, P = 0.7, r = 0.1)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
library(deBif)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    K=q*A
    N0=f*A
    dN <- r*N*(1 - N/K) - E*P*N^2/(N0^2 + N^2)

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dN)))
  })
}

# phaseplane(model, state, parms)
bifurcation(model, state, parms)

# The model for vegetation catastrophes
# -------------------------------------

# Equations:
# ----------
#
# dW     P + k2 W0             W
# -- = R ---------  -  cmax ------ P  - rw W
# dt      P + k2            W + k1
#
# dP           W
# -- = gmax ------ P  - d P  -  b P
# dt        W + k1
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(W = 1.0, P = 1.0)

# Parameters have to be specified as a named vector of parameters.
parms <- c(R = 2.0, k1 = 3.0, k2 = 5.0, cmax = 0.05, gmax = 0.5, rw = 0.1, d = 0.1, W0 = 0.9, b = 0.2)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
vegetation <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    Win = R*(P + k2*W0)/(P + k2)
    gW  = W/(W + k1)
    dW  = Win - cmax*gW*P - rw*W
    dP  = gmax*gW*P - d*P -b*P

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dW, dP)))
  })
}

bifurcation(vegetation, state, parms)


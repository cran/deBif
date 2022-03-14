# The Rosenzweig-McArthur predator-prey model with density-dependent predator mortality
# -------------------------------------------------------------------------------------

# Equations:
# ----------
#
# dR            R       R
# -- = r R (1 - -) - ------- C
# dt            K    1 + ahR
#
# dC          R                 F C^2
# -- = eps -------C - mu C - ----------
# dt       1 + ahR           1 + Ad C^2
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(R = 0.05, C = 0.1)

# Parameters have to be specified as a named vector of parameters.
parms <- c(r = 2.0, K = 0.1, a = 9.1463, h = 0.667, eps = 0.5, mu = 0.1, Ad = 44.444, F = 6.667)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
homoclinic <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dR = r*R*(1 - R/K) - a*R*C/(1 + a*h*R)
    dC = eps*a*R*C/(1 + a*h*R) - mu*C - F*C^2/(1 + Ad*C^2)

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dR, dC)))
  })
}

bifurcation(homoclinic, state, parms)


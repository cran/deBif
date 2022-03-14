# Feedback-control model, Example 5.4 in Kuznetsov (1998; pg. 178)
# ----------------------------------------------------------------

# Equations:
# ----------
#
# dx
# -- = y
# dt
#
# dy
# -- = z
# dt
#
# dz
# -- = - alpha z - beta y - x + x^2
# dt
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(x = 0.0, y = 0.0, z = 0.0)

# Parameters have to be specified as a named vector of parameters.
parms <- c(alpha = 0.5, beta = 1.0)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
adaptx <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dx = y
    dy = z
    dz = -alpha*z - beta*y - x + x^2

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state` 
    return(list(c(dx, dy, dz)))
  })
}

bifurcation(adaptx, state, parms)


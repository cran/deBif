# The Lotka-Volterra competition model
# ------------------------------------

# Equations:
# ----------
#
# dN1              N1 + beta12 N2
# --- = r1 N1 (1 - --------------)
#  dt                     K1
#
# dN2              N2 + beta21 N1
# --- = r2 N2 (1 - --------------)
#  dt                     K2
#

# The initial state of the system has to be specified as a named vector of state values.
state <- c(N1=1, N2=1)

# Parameters have to be specified as a named vector of parameters.
parms <- c(r1=1, r2=1, K1=100, K2=100, beta12=4/3, beta21=3/2)

# The model has to be specified as a function that returns
# the derivatives as a list. You can adapt the body below
# to represent your model
competition <- function(t, state, parms) {
  with(as.list(c(state,parms)), {

    dN1 <- r1*N1*(1 - (N1 + beta12*N2)/K1)
    dN2 <- r2*N2*(1 - (N2 + beta21*N1)/K2)

    # The order of the derivatives in the returned list has to be
    # identical to the order of the state variables contained in 
    # the argument `state`
    return(list(c(dN1, dN2)))
  })
}

# phaseplane(competition, state, parms)

state2 <- c(N1=1, N2=1.98615)
parms2 <-  c(r1=0.1, r2=0.08, K1=100, K2=100, beta12=3/2, beta21=3/2)
phaseplane(competition, state2, parms2)

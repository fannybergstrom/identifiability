## Deterministic SIR model ##
# Input: times, initial states and parameters for the modified SIR model
# Outouts: States for each compartment at each time point


# Define the SIR model
sir_model <- function(time, state, params) {
  with(as.list(c(state, params)), {
    dS <- - (beta_r * I_r + beta_u * I_u) * S / N
    dI_r <- p * (beta_r * I_r + beta_u * I_u) * S / N - gamma * I_r
    dI_u <- (1-p) * (beta_r * I_r + beta_u * I_u) * S / N  - gamma * I_u
    dR_r <- gamma * I_r
    dR_u <- gamma * I_u
    return(list(c(dS, dI_r,dI_u, dR_r, dR_u)))
  })
}

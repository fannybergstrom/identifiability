# Simulation of the stochastic modified SIR model using the Gillespie algorithm

SIR_stoch <- function(params, initial, end_time) {
  # Simulation of the stochastic SIR model using Gillespie algorithm
  
  state <- initial  # Holds the state variables n, I_r, I_u, R_r, R_u
  
  state$S <- initial$pop*(1-params$pi)
  
  result <- list(time = numeric(0), 
                 S = numeric(0), 
                 I_r = numeric(0), 
                 I_u = numeric(0), 
                 R_r = numeric(0), 
                 R_u = numeric(0),
                 N_1 = numeric(0),
                 N_2 = numeric(0),
                 N_3 = numeric(0),
                 N_4 = numeric(0)) # Initialize result list
  
  time <- 0  # Start time
  
  while (time < end_time && (state$I_r + state$I_u) > 0) {  # Run until end time or no infected hosts
    
    # Calculate rates for current state
    rate <- rates(state = list(S = state$S, 
                                I_r = state$I_r, 
                                I_u = state$I_u, 
                                R_r = state$R_r, 
                                R_u = state$R_u), 
                   params)
    
    # When does the next process happen?
    tau <- rexp(1, sum(rate))  # Exponential random variable for the next event time
    
    # Update time
    time <- time + tau
    
    # Determine which process happens
    which <- process(rate)
    
    # Update state based on the selected process
    if (which == 1) {
      state$S <- state$S - 1
      state$I_r <- state$I_r + 1
      state$N_1 <- state$N_1 + 1
    } else if (which == 2) {
      state$S <- state$S - 1
      state$I_u <- state$I_u + 1
      state$N_3 <- state$N_3 + 1
    } else if (which == 3) {  
      state$I_r <- state$I_r - 1
      state$R_r <- state$R_r + 1
      state$N_2 <- state$N_2 + 1
    } else if (which == 4) {
      state$I_u <- state$I_u - 1
      state$R_u <- state$R_u + 1
      state$N_4 <- state$N_4 + 1
    }
    
    # Store results
    result$time <- c(result$time, time)
    result$S <- c(result$S, state$S)
    result$I_r <- c(result$I_r, state$I_r)
    result$I_u <- c(result$I_u, state$I_u)
    result$R_r <- c(result$R_r, state$R_r)
    result$R_u <- c(result$R_u, state$R_u)
    result$N_1 <- c(result$N_1, state$N_1)
    result$N_2 <- c(result$N_2, state$N_2)
    result$N_3 <- c(result$N_3, state$N_3)
    result$N_4 <- c(result$N_4, state$N_4)
  }
  
  return(result)
}

# Function to determine which process happens
process <- function(probs) {
  t <- runif(1) * sum(probs)  # Generate random number
  which <- 1
  s <- probs[1]
  
  # Loop through the probabilities to determine which event occurs
  while (t > s) {
    which <- which + 1
    s <- s + probs[which]
  }
  
  return(which)
}

# Function to calculate process rates for the given state and parameters
rates <- function(state, params) {
  a <- numeric(4)

    # Compute process rates according to the model
  a[1] <-  params$p  *  (params$b_r * state$I_r + params$b_u * state$I_u) * state$S  / pop # Reported infection 
  a[2] <- (1-params$p) *  (params$b_r * state$I_r + params$b_u * state$I_u) * state$S / pop # Unreported infection 
  a[3] <- r * state$I_r  # Reported recovery 
  a[4] <- r * state$I_u  # Unreported recovery
  
  
  
  return(a)
}

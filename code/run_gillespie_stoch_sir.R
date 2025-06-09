## Code for running the Gillespie algorithm

source("code/gillespie_stoch_sir.R")

SIR_stoch_run <- function(pop = 10000, r_f = 1, 
                          p_f = 0.4, pi_f = 0.3, 
                          b_r_f = 2.5, b_u_f = 1.5,
                          I_r_f = pop*0.001, I_u_f = pop*0.001*(1-p_f)/p_f) {
  
  # Simulation of the stochastic SIR model
  
  # Define parameters
  params <- list(pi = pi_f,  # initially immune
                 b_r = b_r_f,   # infection rate reported
                 b_u = b_u_f,   # infection rate unreported
                 r = r_f,   # recovery rate
                 p = p_f    # reporting probability
                  )
  
  
  
  # Initial conditions
  initial <- list(pop = pop,  # population
                  I_r = I_r_f,   # number of infected individuals
                  I_u = I_u_f,
                  R_r = pop*params$pi*params$p,
                  R_u = pop*params$pi*(1-params$p),
                  N_1 = 0,
                  N_2 = 0,
                  N_3 = 0,
                  N_4 = 0)  
  
  # Simulation time and number of runs
  end_time <- 1000  # end of simulation time span starting at 0
  run_count <- 1   # number of runs
  
  # Result containers
  result <- list(time = numeric(0), 
                 S = numeric(0), 
                 I_r = numeric(0),
                 I_u = numeric(0),
                 R_r = numeric(0),
                 R_u = numeric(0),
                 N_1 = numeric(0),
                 N_2 = numeric(0),
                 N_3 = numeric(0),
                 N_4 = numeric(0))
  
  # Simulate several stochastic SIR models and collect data
  for (r in 1:run_count) {
    out <- SIR_stoch(params, initial, end_time)  # Call to SIR_stoch
    
    # Collect results
    result$time <- c(result$time, out$time)
    result$S <- c(result$S, out$S)
    result$I_r <- c(result$I_r, out$I_r)
    result$I_u <- c(result$I_u, out$I_u)
    result$R_r <- c(result$R_r, out$R_r)
    result$R_u <- c(result$R_u, out$R_u)
    result$N_1 <- c(result$N_1, out$N_1)
    result$N_2 <- c(result$N_2, out$N_2)
    result$N_3 <- c(result$N_3, out$N_3)
    result$N_4 <- c(result$N_4, out$N_4)
  }
  
  result
  
}



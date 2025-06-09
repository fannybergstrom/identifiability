### This script load the simulations and estimates p, beta* and pi and saves the estimates.

### Load Required Libraries
pacman::p_load(
  tidyverse, xtable, tidyr, EpiEstim, incidence,
  openxlsx, lubridate, here, ggplot2, viridis,
  ISOweek, zoo, RcppRoll, patchwork, knitr, kableExtra,
  readxl, ggsci, scales, hrbrthemes, deSolve
)


### Initial Parameter Setup
pi_1 = 0.3  # Initial proportion immune
p_1 <- 0.4  # Proportion of population in group 1
b_r = 2.5   # Infection rate for group r
b_u = 1.5   # Infection rate for group u
beta_1 = 2.5*p_1 + 1.5*(1-p_1)  # Transmission rate
r = gamma = 1   # Recovery rate
pop = 10000  # Total population size


# Compute Derived Parameters
rho_1 <- beta_1 * (1 - pi_1) - r
R_E_1 = beta_1 * (1 - pi_1) / r
z_r_1 <- z_r(R_E_1, p = p_1, pi = pi_1)

# Initialize Results
res$N_1 <- 0 
n_sim <- 100  # Number of simulations
# Load External Script for SIR Model Simulation
source(here::here("code/run_gillespie_stoch_sir.R"))


### Function to Find Proportion p
find_p <- function(R_E, pi, z_r) {
  final_eq <- function(p) {
    1 - z_r / (p * (1 - pi)) - exp(-R_E * z_r / (p * (1 - pi)))
  }
  result1 <- uniroot(final_eq, lower = 1e-12, upper = 1 - 1e-12)$root
  return(result1)
}


### Function to Find Proportion pi
find_pi <- function(R_E, p, z_r) {
  final_eq <- function(pi) {
    1 - z_r / (p * (1 - pi)) - exp(-R_E * z_r / (p * (1 - pi)))
  }
  result1 <- uniroot(final_eq, lower = 1e-12, upper = 1 - 1e-12)$root
  return(result1)
}


# Estimate pi from sample first
est_params_1 <- matrix(NA, n_sim, 5)
for(i in 1:n_sim){
  
  n <- str_c("gil_sim_", i)
  # Read simulation file
  res <- read_rds(here::here(str_c("simulations/", n, ".rds")))
  
  t_min <- sum(res$N_1==0)+1
  t_max <- sum(res$N_1<pop*0.075)
  
  lm_res <- lm( log(res$N_1[t_min:t_max]) ~ res$time[t_min:t_max]  )
  rho_hat <- lm_res$coefficients[2]
  
  
  R_E_hat <- (rho_hat + gamma)/gamma
  z_r_hat <- res$N_1[length(res$N_1)]/pop
  
  
  pi_hat <- sample(rbernoulli(pop, 0.3), 1000) %>% sum / 1000
  beta_hat <- beta_tilde(pi_hat, rho = rho_hat)
  p_hat <- find_p(R_E_hat, pi_hat, z_r_hat)
  est_params_1[i,] <- c(beta = beta_hat, pi = pi_hat, p = p_hat, 
                        rho = rho_hat, R_e = R_E_hat )  
  
}

## Now we estimate p first (at the peak)
est_params_2 <- matrix(NA, n_sim, 5)
for(i in 1:n_sim){
  
  n <- str_c("gil_sim_", i)
  # Read simulation file
  res <- read_rds(here::here(str_c("simulations/", n, ".rds")))
  
  t_min <- sum(res$N_1==0)+1
  t_max <- sum(res$N_1<pop*0.075)
  
  lm_res <- lm( log(res$N_1[t_min:t_max]) ~ res$time[t_min:t_max]  )
  rho_hat <- lm_res$coefficients[2]
  
  
  R_E_hat <- (rho_hat + gamma)/gamma
  z_r_hat <- res$N_1[length(res$N_1)]/pop
  
  ind <- res$N_1 %>% which.max() # Find peak
  neg <- rep(0, pop- (res$N_1[ind]+ res$N_3[ind]))
  pos_r <- rep(1, res$N_1[ind])
  pos_u <- rep(2, res$N_3[ind])
  s <- sample(c(pos_r,pos_u,neg), 1000)
  p_hat <- sum(s==1)/(sum(s==1)+sum(s==2))
  z_r_hat <- res$N_1[length(res$N_1)]/pop
  pi_hat <- find_pi(R_E_hat, p_hat, z_r_hat)
  beta_hat <- beta_tilde(pi_hat, rho = rho_hat)
  est_params_2[i,] <- c(beta = beta_hat, pi = pi_hat, p = p_hat, rho = rho_hat, R_e = R_E_hat)
}

# Save results
data.frame(param = c("beta", "pi", "p", "rho", "R_E"),
           truth = c(beta_1, pi_1, p_1, rho_1, R_E_1),
           from_pi = est_params_1 %>% apply(2, mean),
           from_p = est_params_2 %>% apply(2, mean),
           sd_from_pi = est_params_1 %>% apply(2, sd),
           sd_from_p = est_params_2 %>% apply(2, sd)) %>% write_rds("results/summary_tbl.rds")


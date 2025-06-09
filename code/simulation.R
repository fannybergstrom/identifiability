## Code for simulations

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

# Initialize Results
res <- c()
res$N_1 <- 0 
n_sim <- 100  # Number of simulations

# Load External Script for SIR Model Simulation
source(here::here("code/run_gillespie_stoch_sir.R"))



for(i in 1:n_sim){
n <- str_c("gil_sim_", i)
while(sum(res$N_1) < 500) {
  res <- SIR_stoch_run(pop = 10000, 
                       r_f = 1, 
                       p_f = p_1, 
                       pi_f = pi_1,       
                       b_r_f = b_r,
                       b_u_f = b_u,
                       I_r_f = pop*0.001, 
                       I_u_f = pop*0.001*(1-p_1)/p_1)
  
  
}
name <- str_c("gil_sim_", i)
write_rds(res, here::here(str_c("simulations/", name, ".rds")))
res$N_1 <- 0
}

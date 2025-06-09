# Figure 3

# functions
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



find_p <- function(R_E, pi, z_r) {
  final_eq <- function(p) {
    ##########################
    # give the expression here
    ##########################
    1 - z_r/(p*(1-pi))  - exp(-R_E * z_r / (p*(1-pi)) )
  }
  "
  Note: This final size equation is always a solution tau = 0.
  When R_0 <=1, it gives only one solution tau = 0;
  if R_0 > 1, then there are two solutions.
  "
  result1 <- uniroot(final_eq, lower = 1e-12, upper = 1- 1e-12)$root
  #result2 <- uniroot(final_eq, lower = 1e-12, upper = 11e-12, extendInt = "yes")$root
  
  # fill in here (minimum or maximum of result1 and 2 or others?)
  #result <- max(result1, result2)
  return(result1)
}




z_r <- function(R_E, pi, p) {
  final_eq <- function(z_r) {
    1 - z_r/(p*(1-pi))  - exp(-R_E * z_r / (p*(1-pi)) )
  }
  "
  Note: This final size equation is always a solution z = 0.
  When R_0 <=1, it gives only one solution tau = 0;
  if R_0 > 1, then there are two solutions.
  "
  result1 <- uniroot(final_eq, lower = 0, upper = 1)$root
  result2 <- uniroot(final_eq, lower = 1e-12, upper = 1, extendInt = "yes")$root
  
  # fill in here (minimum or maximum of result1 and 2 or others?)
  result <- max(result1, result2)
  
  return(result)
}



# Inital set of parametes
pi_1 = 0.3  # initially immune
p_1 <- 0.4
b_r = 2.5   # infection rate
b_u = 1.5   # infection rate
beta_1 = 1.9
r = 1   # recovery rate
pop = 10000
rho_1 <- beta_1 * (1-pi_1)-r
R_E_1 = beta_1 * (1-pi_1)/r
z_r_1 <- z_r(R_E_1, p = p_1, pi = pi_1)


pi_2 <-  0
beta_2 <- (rho_1+r)/(1-pi_2)
R_E_2 = beta_2 * (1-pi_2)/r
p_2 <- beta_2 * p_1 / beta_1 #find_p(R_E = R_E_2, pi = pi_2, z_r = z_r_1)

z_r_2 <- z_r(R_E_2, p = p_2, pi = pi_2)


# Set the initial values and parameters

init_state <- c(S = pop*(1-pi_1), I_r = pop*0.001, I_u = pop*0.001*(1-p_1)/p_1, R_r = pop*pi_1*p_1, R_u =pop*pi_1*(1-p_1))

init_state2 <- c(S = pop*(1-pi_2), I_r = pop*0.001, I_u = pop*0.001*(1-p_2)/p_2, R_r = pop*pi_2*p_2, R_u = pop*pi_2*(1-p_2))

params <- c(beta_r = beta_1, 
            beta_u = beta_1,
            gamma = r, 
            N = pop,
            p = p_1)

params2 <- c(beta_r = beta_2, 
             beta_u = beta_2,
             gamma = r, 
             N = pop,
             p = p_2)

# Simulate the model

times <- seq(0.1, 35, by = 0.1)
sir_out <- ode(y = init_state, 
               times = times, 
               func = sir_model, 
               parms = params)

sir_df <- data.frame(time = sir_out[,1], 
                     S = sir_out[,2], 
                     I = sir_out[,3]+sir_out[,4], 
                     R = sir_out[,5]+sir_out[,6],
                     I_r = sir_out[,3], 
                     R_r = sir_out[,5])

sir_out2 <- ode(y = init_state2, 
                times = times, 
                func = sir_model, 
                parms = params2)

# Create a ggplot object

sir_df2 <- data.frame(time = sir_out2[,1], 
                      S = sir_out2[,2], 
                      I = sir_out2[,3]+sir_out2[,4], 
                      R = sir_out2[,5]+sir_out2[,6])


data.frame(time = sir_out[,1], 
           I_r = sir_out[,3],
           I_r2 = sir_out2[,3]) %>% 
  pivot_longer(I_r:I_r2) %>% 
  ggplot() + geom_line(aes(time, value/pop, col = name, linetype = name), linewidth= 1) +
  scale_color_lancet(label = expression("1", "2"),  name = expression("Parameter set")) +  
  scale_linetype_manual(values = c(1,2), label = expression("1", "2"),  name = expression("Parameter set")) + theme(legend.title = element_text()) +
  ylab("Pop. share reported infectious") + xlab("Time") + theme_minimal() + theme(legend.position = "bottom")

ggsave("figures/fig3.png", height = 4.5, width = 7)




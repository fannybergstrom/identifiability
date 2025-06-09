# Figure 2

pacman::p_load(
  tidyverse, xtable, tidyr, EpiEstim, incidence,
  openxlsx, lubridate, here, ggplot2, viridis,
  ISOweek, zoo, RcppRoll, patchwork, knitr, kableExtra,
  readxl, ggsci, scales, hrbrthemes, deSolve, ggsci
)


res_1 <- read_rds(here::here("simulations/gil_sim_1.rds"))
t_min <- sum(res_1$N_1==0)+1
t_max <- sum(res_1$N_1<pop*0.075)
lm_res <- lm(log(res_1$N_1[t_min:t_max]) ~ res_1$time[t_min:t_max]  )

p1 <- bind_cols(t = res_1$time-min(res_1$time), Cases = log(res_1$N_1), Model = lm_res$coefficients[1] + lm_res$coefficients[2]* res_1$time) %>% 
  filter(t < 15) %>% 
  pivot_longer(Cases:Model) %>%  
  ggplot(aes(t, value)) + geom_line(aes(col = name)) +
  geom_vline(xintercept=res_1$time[t_max], linetype="dashed") +
  ylab("Cumulative reported incidence (log scale)") + xlab("Time") +
  scale_color_npg(labels = expression(Cases, "Fitted exponential growth rate")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        legend.title = element_blank()) 


# Assuming $p$ known. The solution for equation \@ref(eq:rho).
rho_hat <- lm_res$coefficients[2]


pi_hat <- 0:99/100
beta_hat <- beta_tilde(pi_hat, rho = rho_hat)
p_hat <- c()
n_max <- 63
for(i in 1:n_max){
  p_hat[i] <- find_p(rho_hat+1, pi_hat[i], res$N_1[length(res_1$N_1)]/pop)
}
df_na <- bind_cols(p = p_hat[1:n_max], beta = beta_hat[1:n_max], pi = pi_hat[1:63])

p2<- ggplot(df_na, aes(p, pi)) +
  geom_line(aes(colour = beta), size = 2) +
  ylab(expression(pi))+
  scale_colour_gradient(low = "yellow", high = "red", na.value = NA,  guide = "colourbar", name = expression(beta)) + 
  theme_minimal() +
  theme(legend.title = element_text(), legend.position = "bottom")

p1 + p2 & plot_annotation(tag_levels = "A")

ggsave("figures/fig2.png", height = 4.5, width = 7)

# plot results from SIR model (figure 1)
# param: result list from SIR model

plot_sir <- function(result, pop){
  
  # Extract unique times and the corresponding data
  unique_times <- unique(result$time)
  m <- match(unique_times, result$time)  # Find indices of unique times
  S <- result$S[m]
  I_r <- result$I_r[m]
  I_u <- result$I_u[m]
  R_r <- result$R_r[m]
  R_u <- result$R_u[m]
  N_1 <- result$N_1[m]
  N_2 <- result$N_2[m]
  N_3 <- result$N_3[m]
  N_4 <- result$N_4[m]
  
  # Plot results
  
  p1 <- bind_cols(t = unique_times, S= S, I_r = I_r, I_u = I_u, 
                  R_r = R_r, R_u = R_u) %>% 
    pivot_longer(S:R_u) %>% 
    mutate(name = factor(name, levels = c("S", "I_r", "I_u", "R_r", "R_u")),
           comp = str_extract(name, "[A-Z]"),
           comp = factor(comp, levels = c("S", "I", "R"))) %>% 
    group_by(t, comp) %>% 
    summarise(value = sum(value)) %>% 
    ggplot(aes(x = t, y = value/pop, col = comp)) + geom_line() +
    scale_color_npg()
  
  p2 <- bind_cols(t = unique_times, N_1= N_1, N_2 = N_2, N_3 = N_3, N_4= N_4) %>% 
    pivot_longer(N_1:N_4) %>% 
    ggplot(aes(x = t, y = value/pop, col = name)) + geom_line() +
    scale_color_npg(labels = expression(N[1],N[2], N[3], N[4]))
  
  p <- p1 + p2  &
    theme_minimal()& 
    theme(legend.position = "bottom", 
          legend.title = element_blank()) &
    ylab("Pop. share") & xlab("Time") & plot_annotation(tag_levels = "A")

  p
}


res <- read_rds(here::here("simulations/gil_sim_1.rds"))
plot_sir(res, 10000)
ggsave("figures/fig1.png", height = 4, width = 7)

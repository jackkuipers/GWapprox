
library(tidyverse)

plot_df <- NULL

for (n in 2*(3:6)) {

load(paste0("../sim_results/g_approx",n,".rData"))

out_data %>% mutate(seed = factor(seed), y = -rel_err) %>% filter(gamma < 50, gamma > 10) %>% 
  mutate(nlo = (IG_rat_app_nlo-IG_rat_app)*gamma) -> local_data

results <- local_data %>%
  group_by(seed) %>%
  summarize(
    p = mean(p),
    nlo_app = mean(nlo),
    nlo_est = coef(lm(gamma ~ I(1/y)))[2],
    .groups = "drop"
  )

plot_df <- rbind(plot_df, results)

}



p <- ggplot(plot_df, aes(x = nlo_est, y = nlo_app, colour = as.factor(p), fill = as.factor(p))) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_point(shape = 21, alpha = 0.5, stroke = 1, size = 3) + 
  scale_x_log10() + scale_y_log10() +
  theme_bw(base_size = 16) + labs(colour = "n", fill = "n") +
  xlab("Next-order correction (Monte Carlo estimate)") + 
  ylab("Next-order correction (asymptotics)")

print(p)

ggsave("nlo.pdf", height = 5, width = 6)



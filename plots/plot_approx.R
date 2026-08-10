
library(tidyverse)

for (n in 2*(3:6)) {

if (n == 6) {
  shift_gamma <- 1
  purp_intercept <- -0.7
}

if (n == 8) {
  shift_gamma <- 2.4 # weird number
  purp_intercept <- -0.4
}

if (n == 10) {
  shift_gamma <- 3.4
  purp_intercept <- -0.2
}

if (n == 12) {
  shift_gamma <- 4.3
  purp_intercept <- -0.2
}


load(paste0("../sim_results/g_approx",n,".rData"))

out_data %>% mutate(seed = factor(seed), y = -rel_err) %>% filter(gamma < 50) -> plot_data

ggplot(plot_data, aes(x = gamma-3.5, y = y, colour = as.numeric(seed), group = seed)) +
  geom_line(alpha = 0.4) + scale_x_log10() + scale_y_log10() +
  theme_bw() + scale_color_gradient(low = "blue", high = "darkblue") +
  theme(legend.position = "none") + xlab(expression(gamma - 3.5)) + ylab("Relative error") +
  geom_smooth(data = subset(plot_data, gamma >= 8),
              method = "lm", se = FALSE, linetype = "dotted", linewidth = 0.4)



fit_shifted_log <- function(x, y, min_shift = 1e-6) {
  
  obj <- function(shift) {
    if (shift <= min_shift) return(Inf)
    fit <- lm(log(y) ~ log(x - shift))
    sum(residuals(fit)^2)
  }
  
  opt <- optim(
    par     = min_shift,
    fn      = obj,
    method  = "Brent",
    lower   = min_shift,
    upper   = 5.5
  )
  
  shift_opt <- opt$par
  fit <- lm(log(y) ~ log(x - shift_opt))
  
  list(
    shift = shift_opt,
    model = fit,
    r2    = summary(fit)$r.squared
  )
}


fits <- plot_data %>% filter(gamma > 5) %>% drop_na() %>% 
  #filter(as.numeric(seed) < 11) %>%
  group_by(seed) %>%
  summarise(
    fit = list(
      fit_shifted_log(gamma, y)
    ),
    .groups = "drop"
  )


fit_summary <- fits %>%
  mutate(
    shift = map_dbl(fit, "shift"),
    r2    = map_dbl(fit, "r2"),
    slope = map_dbl(fit, function(f) {coef(f$model)[2]})
  ) %>%
  select(seed, shift, r2, slope)

print(n)
print(mean(fit_summary$shift))
print(mean(fit_summary$slope))
print(sd(fit_summary$slope))

plot_data %>% mutate(log_y = log10(y)) %>% group_by(gamma) %>% summarise(mean_y = mean(log_y)) 

p <- ggplot(plot_data, aes(x = gamma-shift_gamma, y = -rel_err, colour = as.numeric(seed), group = seed)) +
  geom_line(alpha = 0.4) + scale_x_log10() + scale_y_log10() +
  theme_bw(base_size = 16) + scale_color_gradient(low = "blue", high = "darkblue") +
  theme(legend.position = "none") + xlab(bquote(gamma ~ "-" ~ .(shift_gamma))) + ylab("Error") +
  geom_smooth(data = subset(plot_data, gamma >= 8),
              method = "lm", se = FALSE, linetype = "dotted", linewidth = 0.4) + 
  geom_abline(slope = -1, intercept = purp_intercept, colour ="darkorchid4", linewidth=1.5) +
  annotate("text", x = Inf, y = Inf, label = as.expression(bquote(n == .(n))), 
           hjust = 1.5, vjust = 1.5, size = 6)

print(p)

ggsave(paste0("convergeapprox",n,".pdf"), p, height = 4, width = 6)

}


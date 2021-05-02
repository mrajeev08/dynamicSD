# param fig

library(ggplot2)
priors <- readRDS(fp("analysis/out/priors.rds"))

set.seed(125)
prior_dists <- data.table(R0 = priors$R0(1e5), 
                          iota = priors$iota(1e5),
                          k = priors$k(1e5))
priors <- melt(prior_dists)
ints <- data.frame(ints = c(exp(0.2), exp(0), exp(0.25)), 
                   variable = c("R0", "k", "iota"))

ggplot() +
  geom_density(data = priors, aes(x = value), fill = "blue", color = NA, 
               alpha = 0.75) +
  # geom_vline(data = ints, aes(xintercept = ints), color = "darkbl")  +
  facet_wrap(~variable, scales = "free", ncol = 1) + 
  cowplot::theme_half_open(font_size = 12) +
  labs(x = "Parameter value", y = "Prior density") -> priors_out

ggsave("analysis/figs/fig_priors.jpeg", priors_out, height = 8, width = 8)

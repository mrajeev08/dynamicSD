# Parameter estimation figures

library(data.table)
library(ggplot2)
library(foreach)
library(tibble)
library(readr)
library(dplyr)
library(magrittr)
library(glue)
library(ggtext)
library(tidytext)
source("R/utils.R")
source("R/figure_funs.R")

# params for plotting
ppars <- plotpars()

# get best models at each scale
mod_predicted <- read_csv("analysis/out/mod_summary.csv")
mod_predicted %>%
  filter(type == 'full') %>%
  group_by(scale) %>% 
  filter(votes.V1 == max(votes.V1)) -> chosen

# Combine all files -------
# path to files
fpath <- "analysis/out/par_ests/"
fp <- function(x) paste0(fpath, x)

# combine files
test <- list.files(fpath)
test <- test[grep(".rds", test)]
out_all <- lapply(test, function(x) {
  reslt <- readRDS(fp(x)) 
  reslt <-lapply(reslt, append_col, col_name = "type", val = x)
}
)

out_all <- do.call(c, out_all)
out_names <- unique(names(out_all))
out_all <-
  lapply(out_names, 
         function(x) {
           tt <- rbindlist(out_all[names(out_all) == x], fill = TRUE)
           if(nrow(tt) > 0 & "nsim" %in% names(tt)) {
             setnafill(tt, cols = "nsim", fill = 0)
           }
           tt
         }
  )
names(out_all) <- out_names

# Variable importance plots -----
var_imp <- out_all$var_imp
var_imp %<>%
  left_join(select(ppars$summ_stat_labs, 
                   type_stat = type, 
                   summstat)) %>%
  left_join(ppars$mod_labs) %>%
  mutate(type = case_when(nsim == 0 ~ "full", 
                          nsim > 0 ~ "se")) %>%
  filter(interaction(param, incs) != "iota.Fixed")
  
var_importance_grouped <- 
  ggplot(var_imp, aes(x = incs, y = importance, 
                      color = type_stat)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5) + 
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  facet_grid(param ~ scale, scales = "free") +
  ppars$theme_proj() 

var_imp$lab <- ppars$slabs[match(var_imp$summstat, names(ppars$slabs))]
  
var_importance_plot_best <- 
  ggplot(filter(var_imp, type == "full", modname %in% chosen$modname)) + 
  geom_point(aes(x = importance, y = reorder_within(lab, importance, param), 
                 color = type_stat, 
                 group = interaction(scale, nsim), 
                 shape = scale)) +
  scale_y_reordered() +
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  scale_shape_discrete(labels = ppars$scale_labs, name = "Model scale") +
  labs(y = "", tag = "A") +
  facet_wrap(~ param, scales = "free", ncol = 1) +
  ppars$theme_proj() +
  theme(axis.text.y = element_markdown(), 
        axis.text.x = element_text(angle = 0))

# for fpo
var_imp$param <- factor(var_imp$param, levels = c("R0", "iota", "k"))
var_imp$param <- forcats::fct_recode(var_imp$param,
                                   `R[0]` = "R0", 
                                   `iota` = "iota",
                                   `k` = "k")
var_importance_plot_fpo <- 
  ggplot(filter(var_imp, type == "full", modname %in% chosen$modname)) + 
  geom_point(aes(x = importance, y = reorder_within(lab, importance, param), 
                 color = type_stat, 
                 group = interaction(scale, nsim), 
                 shape = scale)) +
  scale_y_reordered() +
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  scale_shape_discrete(labels = ppars$scale_labs, name = "Model scale") +
  labs(y = "") +
  facet_wrap(~ param, scales = "free", ncol = 3, labeller = labeller(param = label_parsed)) +
  ppars$theme_proj() +
  theme(axis.text.y = element_markdown(), 
        axis.text.x = element_text(angle = 45))

ggsave("analysis/figs/fpo/var_imp_pars.jpeg", var_importance_plot_fpo, height = 7, width = 10)

# Prior error rates ----
err <- out_all$err
err %<>%
  left_join(ppars$mod_labs) %>%
  mutate(type = case_when(nsim == 0 ~ "full", 
                          nsim > 0 ~ "se")) %>%
  filter(interaction(param, incs) != "iota.Fixed")
err$param <- factor(err$param, levels = c("R0", "iota", "k"))
err$param <- forcats::fct_recode(err$param,
                                 `R[0]` = "R0", 
                                 `iota` = "iota",
                                 `k` = "k")
err$type[is.na(err$type)] <- "full"
err$nsim[is.na(err$nsim)] <- 0

err_plot_supp <-
  ggplot(err) + 
  geom_line(aes(x = ntree, y = oob_mse, 
                color = scale, 
                group = interaction(scale, type, nsim, move, incs, limits), 
                linetype = type)) +
  scale_linetype_discrete(labels = ppars$run_labs, name = "Model run") + 
  scale_color_brewer(palette = "Dark2", labels = ppars$scale_labs, name = "Model scale") +
  cowplot::theme_half_open(font_size = 12) +
  facet_grid(param ~ incs, scales = "free", labeller = labeller(param = label_parsed)) + 
  labs(x = "Number of trees", y = "Prior error rate")

# Priors vs. posteriors ---- 
posts <- out_all$preds
posts <- posts[data.table(ppars$mod_labs), on = "modname"]
posts <- posts[interaction(param, incs) != "iota.Fixed"]
posts[, type := ifelse(nsim == 0, "full", "se")]
posts$param <- factor(posts$param, levels = c("R0", "iota", "k"))
posts$param <- forcats::fct_recode(posts$param,
                                       `R[0]` = "R0", 
                                       `iota` = "iota",
                                       `k` = "k")

post_dens <- posts[, as.list(density(resp, weights = V1)[c("x", "y")]), 
                   by = c("param", "move", "incs", "scale",
                          "limits", "nsim", "type")]

prior_dens <- posts[, as.list(density(resp)[c("x", "y")]), 
                    by = c("param", "move", "incs", "scale",
                           "limits", "nsim", "type")]

full_posts <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
               fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = post_dens, 
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, type, incs, limits, scale), 
                   color = interaction(move, limits), 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_discrete(name = "Model movement") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, labeller = labeller(param = label_parsed)) 


post_joint <- posts[, .(check = all(V1 > 0)), 
                    by = c("move", "incs", "scale", "modindex",
                           "limits", "nsim", "type", "sim_id", "modname")][check == TRUE]
post_joint[, include := paste(sim_id, modname, nsim, sep = "_")]
posts[, include := paste(sim_id, modname, nsim, sep = "_")]
post_joint <- posts[include %in% post_joint$include] 
post_joint_dens <- post_joint[, as.list(density(resp, 
                                                weights = V1/sum(V1))[c("x", "y")]), 
                         by = c("param", "move", "incs", "scale",
                                "limits", "nsim", "type")]

# Table of parameter estimates ---
out_all$stats %>%
  left_join(ppars$mod_labs) %>%
  mutate(type = case_when(nsim == 0 ~ "full", 
                          nsim > 0 ~ "se")) %>%
  filter(interaction(param, incs) != "iota.Fixed", 
         type == "full") -> par_stats_ind
par_stats_ind$param <- factor(par_stats_ind$param, levels = c("R0", "iota", "k"))
par_stats_ind$param <- forcats::fct_recode(par_stats_ind$param,
                                           `R[0]` = "R0", 
                                           `iota` = "iota",
                                           `k` = "k")
posts_sample <- post_joint_dens[, sample(x, 1e4, prob = y, replace = TRUE), 
                                by = c("param", "move", "incs", "scale",
                                      "limits", "nsim", "type")]
par_stats_joint <- posts_sample[, .(mean = mean(V1), 
                                     median = median(V1), 
                                     quantile_0.975 = quantile(V1, 0.975), 
                                     quantile_0.025 = quantile(V1, 0.025), 
                                    par_type = "Joint"), 
                                 by = c("param", "move", "incs", "scale",
                                        "limits", "nsim", "type")]
prior_stats <- posts[, .(mean = mean(resp), 
                         median = median(resp), 
                         quantile_0.975 = quantile(resp, 0.975), 
                         quantile_0.025 = quantile(resp, 0.025), 
                         par_type = "Prior"), 
                     by = c("param", "move", "incs", "scale",
                            "limits", "nsim", "type")]
par_stats_ind <- par_stats_ind[, mean := expectation][, -c("variance_postmse", "modname", "variance_cdf", "expectation", "post_nmae_mean", "modindex")][, par_type := "Independent"]

par_all <- rbind(prior_stats, par_stats_joint, par_stats_ind)

write_csv(par_all, "analysis/out/par_stats_comp.csv")

joint_posts <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
                   fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = post_joint_dens, 
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, type, move, incs, limits, scale), 
                   color = interaction(move, limits), 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_discrete(name = "Model movement") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, 
             labeller = labeller(param = label_parsed)) 

## FPO the best one only ----
chosen %>%
  mutate(chosen_names = interaction(move, incs, scale, limits)) %$%
  chosen_names -> cnames

posts_best <- post_joint_dens[interaction(move, incs, scale, limits) %in% cnames]

post_best_fig <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
               fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = posts_best[nsim == 0],
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, incs, limits, scale), 
                   fill = scale),
               alpha = 0.5, color = NA) +
  facet_wrap(~ param, scales = "free", ncol = 1, 
             labeller = labeller(param = label_parsed)) +
  scale_fill_brewer(name = "Scale", palette = "Dark2") +
  cowplot::theme_half_open(font_size = 12) +
  labs(x = "Parameter estimate", y = "Density")
ggsave("analysis/figs/fpo/posts_best.jpeg", post_best_fig, height = 7, width = 7)


posts_all <-
  rbind(post_joint_dens[type == "full"][, post_type := "Joint"],
      post_dens[type == "full"][, post_type := "Independent"])

joint_vs_full <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
               fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = posts_all,
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, incs, limits, scale, 
                                       post_type), 
                   color = post_type, 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_brewer(palette = "Set1", name = "Posterior type") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, labeller = labeller(param = label_parsed)) +
  cowplot::theme_half_open(font_size = 12)

joint_vs_full <-
  ggplot() + 
  geom_density(data = prior_dens[], 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
               fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = posts_all,
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, incs, limits, scale, 
                                       post_type), 
                   color = post_type, 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_brewer(palette = "Set1", name = "Posterior type") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, labeller = labeller(param = label_parsed)) +
  cowplot::theme_half_open(font_size = 12)


# Ability to recover all params
test_preds <- out_all$preds_test

test_preds %<>%
  left_join(ppars$mod_labs) %>%
  filter(interaction(param, incs) != "iota.Fixed")
test_preds$param <- factor(test_preds$param, levels = c("R0", "iota", "k"))
test_preds$param <- forcats::fct_recode(test_preds$param,
                                       `R[0]` = "R0", 
                                       `iota` = "iota",
                                       `k` = "k")

test_preds %>%
  filter(interaction(move, incs, scale, limits) %in% cnames) -> test_preds_best

test_param_recovery <-
  ggplot(test_preds) + 
  geom_pointrange(aes(x = true_val, y = expectation,
                      ymin = quantile_0.025, ymax = quantile_0.975, 
                      color = move, shape = incs),
               alpha = 0.5) +
  scale_shape(name = "Introductions") +
  scale_color_discrete(name = "Model movement") +
  labs(x = "Simulated (true) value", y = "Estimated value") + 
  cowplot::theme_half_open(font_size = 12) +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, 
             labeller = labeller(param = label_parsed)) 

test_param_recovery_best <-
  ggplot(test_preds_best) + 
  geom_pointrange(aes(x = true_val, y = expectation,
                      ymin = quantile_0.025, ymax = quantile_0.975,
                      color = scale),
                  alpha = 0.5) +
  facet_wrap(~ param, scales = "free", ncol = 1, 
             labeller = labeller(param = label_parsed)) +
  scale_color_brewer(name = "Scale", palette = "Dark2") +
  cowplot::theme_half_open(font_size = 12) +
  labs(x = "Simulated (true) value", y = "Estimated value")
ggsave("analysis/figs/sfig_param_recov.jpeg", test_param_recovery_best, height = 8, width = 8)

# out all figs ----
ggsave("analysis/figs/mfig_posts_best.jpeg", post_best_fig, height = 8, width = 8)
ggsave("analysis/figs/mfig_posts_vimp.jpeg", var_importance_plot_best, height = 8, width = 8)

# supplementary
ggsave("analysis/figs/sfig_posts_vimpgr.jpeg", var_importance_grouped, height = 8, width = 8)
ggsave("analysis/figs/sfig_posts_err.jpeg", err_plot_supp, height = 8, width = 8)
ggsave("analysis/figs/sfig_posts_jntvsfull.jpeg", joint_vs_full,  height = 8, width = 8)


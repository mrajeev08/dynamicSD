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
           if(nrow(tt) > 0) {
             setnafill(tt, cols = "nsim", fill = 0)
             tt
           }
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

# Prior error rates ----
err <- out_all$err
err %<>%
  left_join(ppars$mod_labs) %>%
  mutate(type = case_when(nsim == 0 ~ "full", 
                          nsim > 0 ~ "se")) %>%
  filter(interaction(param, incs) != "iota.Fixed")

err_plot_supp <-
  ggplot(err) + 
  geom_line(aes(x = ntree, y = oob_mse, 
                color = scale, 
                group = interaction(scale, type, nsim, move, incs, limits), 
                linetype = type)) +
  scale_linetype_discrete(labels = ppars$run_labs, name = "Model run") + 
  scale_color_brewer(palette = "Dark2", labels = ppars$scale_labs, name = "Model scale") +
  cowplot::theme_half_open(font_size = 12) +
  facet_grid(param ~ incs, scales = "free", labeller = param_labeller) + 
  labs(x = "Number of trees", y = "Out-of-bag error rate")

# Table of parameter estimates ---
out_all$stats %<>%
  left_join(ppars$mod_labs) %>%
  mutate(type = case_when(nsim == 0 ~ "full", 
                          nsim > 0 ~ "se")) %>%
  filter(interaction(param, incs) != "iota.Fixed", 
         type == "full") -> par_stats
write_csv(par_stats, "analysis/out/par_stats.csv")

# Priors vs. posteriors ---- 
posts <- out_all$preds
posts <- posts[data.table(ppars$mod_labs), on = "modname"]
posts <- posts[interaction(param, incs) != "iota.Fixed"]
posts[, type := ifelse(nsim == 0, "full", "se")]
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
                   group = interaction(nsim, param, move, type, limits, scale), 
                   color = interaction(move, limits), 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_discrete(name = "Model movement") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2) 


post_joint <- posts[, .(check = all(V1 > 0)), 
                    by = c("move", "incs", "scale",
                           "limits", "nsim", "type", "sim_id")][check == TRUE]

post_joint <-  posts[sim_id %in% post_joint$sim_id]
post_joint_dens <- post_joint[, as.list(density(resp, 
                                                weights = V1/sum(V1))[c("x", "y")]), 
                         by = c("param", "move", "incs", "scale",
                                "limits", "nsim", "type")]

joint_posts <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
                   fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = post_joint_dens, 
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, limits, scale), 
                   color = interaction(move, limits), 
                   linetype = incs)) +
  scale_linetype(name = "Introductions") +
  scale_color_discrete(name = "Model movement") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2) 

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
  facet_wrap(param ~ scale, scales = "free", ncol = 2) +
  cowplot::theme_half_open(font_size = 12)

chosen %>%
  mutate(chosen_names = interaction(move, incs, scale, limits)) %$%
  chosen_names -> cnames
         
posts_best <- posts_all[ interaction(move, incs, scale, limits) %in% cnames]

posts_best$param <- factor(posts_best$param, levels = c("R0", "iota", "k"))
posts_best$param <- forcats::fct_recode(posts_best$param,
                                        `R[0]` = "R0", 
                                        `iota` = "iota",
                                        `k` = "k")
prior_dens$param <- factor(prior_dens$param, levels = c("R0", "iota", "k"))
prior_dens$param <- forcats::fct_recode(prior_dens$param,
                                        `R[0]` = "R0", 
                                        `iota` = "iota",
                                        `k` = "k")

post_best_fig <-
  ggplot() + 
  geom_density(data = prior_dens, 
               aes(x = x, weight = y,
                   group = interaction(nsim, param)), 
               fill = "grey", color = "NA", alpha = 0.5) +
  geom_density(data = posts_best,
               aes(x = x, weight = y, 
                   group = interaction(nsim, param, move, incs, limits, scale, 
                                       post_type), 
                   fill = post_type), alpha = 0.5, color = NA) +
  scale_fill_brewer(palette = "Set1", name = "Posterior type") +
  facet_wrap(param ~ scale, scales = "free", ncol = 2, 
             labeller = labeller(param = label_parsed)) +
  cowplot::theme_half_open(font_size = 12) +
  labs(x = "Parameter estimate", y = "Density")

# out all figs ----
ggsave("analysis/figs/mfig_posts_best.jpeg", post_best_fig)
ggsave("analysis/figs/mfig_posts_vimp.jpeg", var_importance_plot_best)

# supplementary

ggsave("analysis/figs/sfig_posts_vimpgr.jpeg", var_importance_grouped)
ggsave("analysis/figs/sfig_posts_err.jpeg", err_plot_supp)
ggsave("analysis/figs/sfig_posts_jntvsfull.jpeg", joint_vs_full)


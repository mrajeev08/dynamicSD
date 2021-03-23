# Main and supplementary figures and tables for model selection --------

library(data.table)
library(ggplot2)
library(foreach)
library(tibble)
library(readr)
library(dplyr)
library(magrittr)
library(glue)
library(ggtext)
source("R/utils.R")
source("R/figure_funs.R")

# params for plotting
ppars <- plotpars()

# Combine all files -------
# path to files
fpath <- "analysis/out/mod_comp/"
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
             tt[, c("type", "scale") := tstrsplit(type, "_", fixed = TRUE)]
             tt[, scale := gsub(".rds", "", scale)]
             tt
           }
        }
    )
names(out_all) <- out_names

# Variable importance plots ----
var_imp <- out_all$var_imp
var_imp %<>%
  left_join(select(ppars$summ_stat_labs, 
                   type_stat = type, 
                   summstat)) 

var_importance_plot_supp <- 
  ggplot(var_imp) + 
  geom_point(aes(x = importance, y = reorder(summstat, importance), 
                 color = type_stat, 
                 group = interaction(scale, nsim), 
                 shape = type)) +
  scale_y_discrete(labels = ppars$slabs) +
  scale_shape_manual(values = c(16, 1), 
                     labels = ppars$run_labs, name = "Model run") + 
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  facet_wrap(~scale, labeller = labeller(scale = ppars$scale_labs)) +
  labs(y = "") +
  ppars$theme_proj() +
  theme(axis.text.y = element_markdown())


var_importance_plot_main <- 
  ggplot(filter(var_imp, type == "full")) + 
  geom_point(aes(x = importance, y = reorder(summstat, importance), 
                 color = type_stat, 
                 group = interaction(scale, nsim), 
                 shape = scale)) +
  scale_y_discrete(labels = ppars$slabs) +
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  scale_shape_discrete(labels = ppars$scale_labs, name = "Model scale") +
  labs(y = "", tag = "A") +
  ppars$theme_proj() +
  theme(axis.text.y = element_markdown())


# oob errors ----
err <- out_all$err
err_plot_supp <-
  ggplot(err) + 
  geom_line(aes(x = ntree, y = error.rate, 
                color = scale, 
                group = interaction(scale, type, nsim), 
                linetype = type)) +
  scale_linetype_discrete(labels = ppars$run_labs, name = "Model run") + 
  scale_color_brewer(palette = "Dark2", labels = ppars$scale_labs, name = "Model scale") +
  ppars$theme_proj() +
  labs(x = "Number of trees", y = "Prior error rate")

# model ranks -----
out_all$mod_predicted %>%
  mutate(scale = case_when(scale == "grid" ~ "1x1 km",
                           scale == "vill" ~ "Village")) %>%
  left_join(ppars$mod_labs) ->  mod_predicted

model_ranks_supp <-
  ggplot(filter(mod_predicted, type == "full")) + 
  geom_col(aes(x = reorder(interaction(move, limits), votes.V1), 
               y = votes.V1/500*100, fill = scale), 
           position = "dodge") + 
  facet_wrap(~incs,  labeller = labeller(incs = ppars$fct_incs)) + 
  scale_fill_brewer(palette = "Dark2", 
                    labels = ppars$scale_labs, 
                    name = "Model scale") +
  ppars$theme_proj() +
  labs(x = "Movement", y = "Votes (% of trees that \n selected given model)")
write_csv(mod_predicted, "analysis/out/mod_comp/summary.csv")

# lda plots ----
mod_predicted %>%
  filter(type == 'full') %>%
  group_by(scale) %>% 
  filter(votes.V1 == max(votes.V1)) -> chosen

lda_preds <- out_all$lda_predicted
lda_obs <- out_all$lda_obs
lda_preds %>%
  mutate(scale = case_when(scale == "grid" ~ "1x1 km",
                           scale == "vill" ~ "Village"), 
         type_obs = "Predicted") %>%
  right_join(chosen) -> lda_chosen

lda_obs %>%
  mutate(scale = case_when(scale == "grid" ~ "1x1 km",
                           scale == "vill" ~ "Village"), 
         type_obs = "Observed") %>%
  right_join(chosen) %>%
  bind_rows(lda_chosen) -> lda_all

# pick the best village & grid cell mod and show the lda plot with the data
lda_main <- 
  ggplot(data = lda_all) + 
  geom_point(aes(x = LD1, y = LD2, color = scale, shape = type_obs, 
                 alpha = type_obs)) +
  scale_shape_manual(values = c(2, 16), name = "") + 
  scale_alpha_manual(values = c(1, 0.1), name = "") +
  scale_color_brewer(palette = "Dark2", 
                    labels = ppars$scale_labs, 
                    name = "Model scale") +
  ppars$theme_proj() +
  labs(x = "LD Axis 1", y = "LD Axis 2")

# supplementary plot with hex bins showing all comparisons
lda_preds %>%
  mutate(scale = case_when(scale == "grid" ~ "1x1 km",
                           scale == "vill" ~ "Village"), 
         type_obs = "Predicted") %>%
  right_join(mod_predicted) -> lda_all

lda_obs %>%
  mutate(scale = case_when(scale == "grid" ~ "1x1 km",
                           scale == "vill" ~ "Village"), 
         type_obs = "Observed") %>%
  left_join(mod_predicted) -> lda_obs_all

lda_supp <- 
  ggplot(filter(lda_all, type == "full")) +
  geom_hex(aes(x = LD1, y = LD2, colour = ..count.., fill = scale)) +
  geom_point(data = filter(lda_obs_all, type == "full"), 
             aes(x = LD1, y = LD2, fill = scale), shape = 2) +
  scale_color_distiller(palette = "Greys", direction = 1) +
  scale_fill_brewer(palette = "Dark2", 
                     labels = ppars$scale_labs, 
                     name = "Model scale") +
  facet_grid(limits ~ interaction(move, incs)) + 
  cowplot::theme_half_open(font_size = 12)


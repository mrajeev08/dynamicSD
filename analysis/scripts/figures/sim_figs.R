# Simulations from posterior figures

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
library(lubridate)
source("R/utils.R")
source("R/figure_funs.R")

# params for plotting
ppars <- plotpars()
date_brks <- seq(0, 228, by = 24)
date_labs <- format(seq(ymd('2002-01-01'), ymd('2020-12-31'), by='24 months'), "%Y")

# get best models at each scale
# get best models at each scale
mod_predicted <- read_csv("analysis/out/mod_summary.csv")
mod_predicted %>%
  filter(type == 'full') %>%
  group_by(scale) %>% 
  filter(votes.V1 == max(votes.V1)) -> chosen

# Combine all files -------
fpath <- "analysis/out/post_sims/"
fp <- function(x) paste0(fpath, x)

# combine files
test <- list.files(fpath)
fp_curves <- test[grep("curve_scores", test)]
fp_scores <- test[grep("ts_scores", test)]
fp_sims <- test[grep("sims", test)]

out_curves <- rbindlist(lapply(fp_curves, 
                               function(x) fread(fp(x))))
out_scores <-  rbindlist(lapply(fp_scores, 
                                function(x) fread(fp(x))))

# getting sims in vaccination setting ----

# all ----
# bring in observed data
case_ts <- fread("data/case_ts.csv")

# merge with correct labels
sims_vacc <- out_scores[vacc_type == 'vill_vacc']
sims_vacc <- sims_vacc[data.table(ppars$mod_labs), on = "modname"]

# melt to show top 5 sims to data and top 5 most central sims
best_sims <- melt(sims_vacc, 
                  id.vars = c("cal_month", "scale",
                              "incs", "limits", "move", 
                              "vacc_type", "post_type", "modname"))[grep("best_times|best_dat", variable)]

best_sims$score_type <- case_when(grepl("best_times", best_sims$variable) ~ "Most central", 
                                  grepl("best_dat", best_sims$variable) ~ "Best fit to data")
best_sims[, c("best", "class", "rank") := tstrsplit(variable, "_", fixed = TRUE)]
best_sims$rank <- as.numeric(best_sims$rank)

all_post_sims_joint <- 
  ggplot(sims_vacc[post_type == "joint"]) +
  geom_ribbon(aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims[post_type == "joint"],
            aes(x = cal_month, y = as.numeric(value), 
                group =interaction(variable, modname), 
                color = score_type)) +
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "black") +
  scale_color_manual(values = c("#d95f02", "#7570b3"), 
                     name = "") + 
  facet_wrap(scale ~ incs, scales = "free") +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_post_sims_full <-
  ggplot(sims_vacc[post_type == "full"]) +
  geom_ribbon(aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims[post_type == "joint"],
            aes(x = cal_month, y = as.numeric(value), 
                group =interaction(variable, modname), 
                color = score_type)) +
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "black") +
  scale_color_manual(values = c("#d95f02", "#7570b3"), 
                     name = "") + 
  facet_wrap(scale ~ incs, scales = "free") +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Cases")


# comparing scores (all supp) ----
score_crps <-
  ggplot(sims_vacc) +
  geom_line(aes(x = cal_month, y = crps_score, 
                group = interaction(modname, vacc_type, post_type), 
                color = scale)) +
  scale_color_brewer(palette = "Dark2", labels = ppars$scale_labs, name = "Model scale") +
  facet_grid(post_type ~ incs, scales = "free") +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "CRPS")
  

score_dss <-
  ggplot(sims_vacc) +
  geom_line(aes(x = cal_month, y = dss_score, 
                group = interaction(modname, vacc_type, post_type), 
                color = scale)) +
  scale_color_brewer(palette = "Dark2", labels = ppars$scale_labs, name = "Model scale") +
  facet_grid(post_type ~ incs, scales = "free") +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "DSS")

# getting endemic sims -----
out_endemic <- rbindlist(lapply(fp_sims, function(x) {
  dt <- fread(fp(x))
  dt[, grp := rleid(I_ts), by = c("modname", "vacc_type", "post_type", "sim")]
  dt[, .(cal_month = cal_month[which.max(grp)]), by = c("modname", "vacc_type",                                                       "post_type", "sim")]
}))

count_endemic <-
  out_endemic[, .(month = 1:228,
                  count = cumsum(tabulate(cal_month, nbins = 228))/.N),
              by = c("modname", "vacc_type", "post_type")]

count_endemic <- count_endemic[data.table(ppars$mod_labs), on = "modname"]

endemic_sims_all <-
  ggplot(count_endemic[month < 225 & vacc_type == "endemic"]) +
  geom_line(aes(x = month, y = count, color = scale, 
                group = modname)) +
  facet_grid(post_type ~ incs, scales = "free") +
  scale_color_brewer(palette = "Dark2", 
                     labels = ppars$scale_labs, 
                     name = "Model scale") +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Proportion of sims with \n pop decline > 15 %")

# best ----
chosen %>%
  mutate(chosen_names = interaction(move, incs, scale, limits)) %$%
  chosen_names -> cnames

sims_vacc_best <- sims_vacc[interaction(move, incs, scale, limits) %in% cnames]
best_sims_best <- best_sims[interaction(move, incs, scale, limits) %in% cnames]

best_posts_joint <- 
  ggplot() +
  geom_ribbon(data = sims_vacc_best[post_type == "joint"],
              aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims_best[post_type == "joint" & rank < 4],
            aes(x = cal_month, y = as.numeric(value), 
                group =interaction(variable, modname), 
                color = score_type)) +
  scale_color_manual(values = c("#d95f02", "#7570b3"), 
                     name = "") + 
  scale_x_discrete(breaks = seq(0, 228, 12), 
                   labels = ) + 
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "black") +
  facet_wrap(scale ~ incs, scales = "fixed", ncol = 1) +
  
  ppars$theme_proj() +
  labs(x = "", y = "Cases")

count_endemic <- count_endemic[interaction(move, incs, scale, limits) %in% cnames]

best_endemic_prop <-
  ggplot(count_endemic[month < 225 & vacc_type == "endemic" & post_type == "joint"]) +
  geom_line(aes(x = month, y = count, color = scale, 
                group = modname)) +
  scale_color_brewer(palette = "Dark2", 
                     labels = ppars$scale_labs, 
                     name = "Model scale") + 
  ppars$theme_proj() +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  labs(x = "", y = "Proportion of sims with \n pop decline > 15 %")

# out figs

ggsave("analysis/figs/mfig_sims_best.jpeg", best_posts_joint)
ggsave("analysis/figs/sfig_sims_endemic.jpeg", best_endemic_prop) 
ggsave("analysis/figs/sfig_sims_full.jpeg", all_post_sims_full)
ggsave("analysis/figs/sfig_sims_joint.jpeg", all_post_sims_joint)
ggsave("analysis/figs/sfig_sims_crps.jpeg", score_crps)
ggsave("analysis/figs/sfig_sims_dss.jpeg", score_dss)

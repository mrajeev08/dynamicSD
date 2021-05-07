# Sim figs counterfactuals 

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
library(tidyr)
source("R/utils.R")
source("R/figure_funs.R")

# params for plotting
ppars <- plotpars()
date_brks <- seq(0, 228, by = 24)
date_labs <- format(seq(ymd('2002-01-01'), ymd('2020-12-31'), by='24 months'), "%Y")

# Combine all files -------
fp <- here::here

# combine files
test <- list.files("analysis/out/post_sims", full.names = TRUE)
test_names <- list.files("analysis/out/post_sims")

# Filter to best mod at each scale
mod_comp <- readRDS(fp("analysis/out/mod_comp/full_all.rds"))
mod_comp$modlookup %>%
  mutate(modindex = paste0("g", modindex)) %>%
  separate(modname, into = c("scale", "limits", "intros", "move"), sep = "_", 
           remove = FALSE) %>%
  left_join(mod_comp$mod_predicted) %>%
  mutate(N = sum(votes.V1), prob = votes.V1/N) -> lookup

lookup %>%
  group_by(scale) %>%
  arrange(desc(votes.V1)) %>%
  slice(1) -> chosen

test <- test[grep("grid1x1_nolim_estincs_kern|vill_nolim_estincs_kern", test)]
test_names <- test_names[grep("grid1x1_nolim_estincs_kern|vill_nolim_estincs_kern", test_names)]
fp_curves <- test[grep("curve_scores", test)]
fp_scores <- test[grep("ts_scores", test)]
fp_sims <- test[grep("sims", test_names)]
fp_vax <- test[grep("vax", test_names)]
out_curves <- rbindlist(lapply(fp_curves, 
                               function(x) fread(fp(x))))
out_scores <-  rbindlist(lapply(fp_scores, 
                                function(x) fread(fp(x))))

out_sims <-  rbindlist(lapply(fp_sims, 
                                function(x) fread(fp(x))), fill = TRUE)

# bring in observed data
case_ts <- fread("data/case_ts.csv")

# merge with correct labels vill vacc original ----
sims_vacc <- out_scores[vacc_type == 'vill_vacc']
sims_vacc <- sims_vacc[data.table(ppars$mod_labs), on = "modname",  nomatch = 0]

# melt to show top 5 sims to data and top 5 most central sims
best_sims <- melt(sims_vacc, 
                  id.vars = c("cal_month", "scale",
                              "incs", "limits", "move", 
                              "vacc_type", "post_type", "modname"))[grep("best_times|best_dat", variable)]

best_sims$score_type <- case_when(grepl("best_times", best_sims$variable) ~ "Most central", 
                                  grepl("best_dat", best_sims$variable) ~ "Best fit to data")
best_sims[, c("best", "class", "rank") := tstrsplit(variable, "_", fixed = TRUE)]
best_sims$rank <- as.numeric(best_sims$rank)

sims_best_sd <- 
  ggplot(sims_vacc) +
  geom_ribbon(aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims,
            aes(x = cal_month, y = as.numeric(value), 
                group =interaction(variable, modname), 
                color = score_type)) +
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "black") +
  scale_color_manual(values = c("#d95f02", "#7570b3"), 
                     name = "") + 
  facet_wrap(~scale, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = date_brks, labels = date_labs) +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# merge with correct labels vill vacc smoothed ----
sims_vacc <- out_scores[vacc_type == 'vill_vacc_smoothed']
sims_vacc <- sims_vacc[data.table(ppars$mod_labs), on = "modname",  nomatch = 0]

# melt to show top 5 sims to data and top 5 most central sims
best_sims <- melt(sims_vacc, 
                  id.vars = c("cal_month", "scale",
                              "incs", "limits", "move", 
                              "vacc_type", "post_type", "modname"))[grep("best_times|best_dat", variable)]

best_sims$score_type <- case_when(grepl("best_times", best_sims$variable) ~ "Most central", 
                                  grepl("best_dat", best_sims$variable) ~ "Best fit to data")
best_sims[, c("best", "class", "rank") := tstrsplit(variable, "_", fixed = TRUE)]
best_sims$rank <- as.numeric(best_sims$rank)

ggplot(sims_vacc) +
  geom_ribbon(aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims,
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

# vill vacc all ----
sims_vacc <- out_scores[vacc_type == 'vill_vacc_all']
sims_vacc <- sims_vacc[data.table(ppars$mod_labs), on = "modname",  nomatch = 0]

# melt to show top 5 sims to data and top 5 most central sims
best_sims <- melt(sims_vacc, 
                  id.vars = c("cal_month", "scale",
                              "incs", "limits", "move", 
                              "vacc_type", "post_type", "modname"))[grep("best_times|best_dat", variable)]

best_sims$score_type <- case_when(grepl("best_times", best_sims$variable) ~ "Most central", 
                                  grepl("best_dat", best_sims$variable) ~ "Best fit to data")
best_sims[, c("best", "class", "rank") := tstrsplit(variable, "_", fixed = TRUE)]
best_sims$rank <- as.numeric(best_sims$rank)

ggplot(sims_vacc) +
  geom_ribbon(aes(x = cal_month, ymin = min, ymax = max,
                  group = modname), 
              fill = "grey", alpha = 0.5) + 
  geom_line(data= best_sims,
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

# merge with correct labels for endemic ----
sims_vacc_all <- out_sims[vacc_type == 'endemic']
sims_vacc <- # group and get the stats (min & max)
  sims_vacc_all[, .(i_loc = mean(I_ts - incs_ts), i_loc_max = max(I_ts - incs_ts),
                    i_loc_min = min(I_ts - incs_ts),
                    i_incs = mean(incs_ts), i_incs_max = max(incs_ts),
                    i_incs_min = min(incs_ts)), by = c("cal_month", "modname", "stopped")]
sims_vacc <- sims_vacc[data.table(ppars$mod_labs), on = "modname", nomatch = 0]


endemic_scenario <-
  ggplot() +
  geom_ribbon(data = sims_vacc,
              aes(x = cal_month, ymin = i_loc_min, ymax = i_loc_max, 
                  group = interaction(stopped, scale)),
              fill = "grey", alpha = 0.5) +
  geom_line(data = sims_vacc_all, 
            aes(x = cal_month, y = I_ts - incs_ts, 
                group = interaction(stopped, modname, sim), 
                color = "District transmission"), alpha = 0.1) +
  geom_line(data = sims_vacc_all, 
            aes(x = cal_month, y = incs_ts, 
                group = interaction(stopped, modname, sim),
                color = "Introductions"),
                alpha = 0.1) +
  facet_grid(stopped~scale) +
  scale_color_manual(values = c("navy", "darkred"), 
                     name = "Type of transmission") +
  cowplot::theme_half_open(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Simulation month", y = "Detected cases")
ggsave("analysis/figs/fpo/endemic_scenario.jpeg", endemic_scenario, 
       height = 6, width = 6)

# Last bit = chain stats ----
best_chain_stats <- list.files("analysis/out/trees/", full.names = TRUE)
out_sims <-  rbindlist(lapply(best_chain_stats, 
                              function(x) fread(fp(x))), fill = TRUE)

consens <- fread("analysis/out/ttrees/consensus_links_ln0.95.csv")
load("data/sd_case_data.rda")
sd_case_data %>%
  select(date_symptoms = symptoms_started, x_coord = utm_easting, 
         y_coord = utm_northing, id) %>%
  left_join(consens, by = c("id" = "id_case")) %>%
  mutate(date_symptoms = dmy(date_symptoms)) %>%
  as.data.table(.) -> consens

attrs <- consens[, .(id, t = date_symptoms - min(date_symptoms), 
                     x_coord, 
                     y_coord)]
library(igraph)
I_gr <- data.table(from = consens$id_progen[!is.na(consens$id_progen)], 
                   to = consens$id[!is.na(consens$id_progen)])
gr <- graph_from_data_frame(d = I_gr,
                            vertices = attrs,
                            directed = TRUE)
V(gr)$membership <- components(gr)$membership

summ_stats <- 
  data.table(membership = vertex_attr(gr, "membership"),
             id_case = vertex_attr(gr, "name"),
             t_days = vertex_attr(gr, "t"))
summ_stats[, c("start_date", "end_date") := .(min(t_days), max(t_days)), by = "membership"]
chain_dt <-
  summ_stats[, .(length_wks = (end_date[1] - start_date[1])/7, 
                                   size = .N), 
                               by = "membership"]
out_sims %>%
  select(membership, length_wks, size, Type = modname) %>%
  mutate(Type = case_when(Type == "grid1x1_nolim_estincs_kern" ~ "1x1 km", 
                          Type == "vill_nolim_estincs_kern" ~ "Village")) %>%
  bind_rows(mutate(chain_dt, Type = "Empirical")) -> chain_stats

sizes <-
  ggplot(chain_stats) +
  geom_histogram(aes(x = size, fill = Type, y = stat(density*width)), 
                 binwidth = 1, position = "dodge") +
  scale_x_continuous(trans = "sqrt", limits = c(0, 1000)) +
  scale_fill_manual(values = c('#1b9e77', '#7570b3', '#d95f02')) +
  labs(x = "Chain size", y = "Proportion") +
  cowplot::theme_minimal_hgrid()

pers <-
  ggplot(chain_stats) +
  geom_histogram(aes(x = length_wks, fill = Type, y = stat(density*width)), 
                 binwidth = 1, position = "dodge") +
  scale_x_continuous(trans = "sqrt", limits = c(0, 1000)) +
  scale_fill_manual(values = c('#1b9e77', '#7570b3', '#d95f02')) +
  labs(x = "Persistence of chain (weeks)", y = "Proportion") +
  cowplot::theme_minimal_hgrid()
library(patchwork)

out <- sizes / pers + plot_layout(guides = "collect")

ggsave("analysis/figs/fpo/chain_stats.jpeg", out, height = 5, width = 7)

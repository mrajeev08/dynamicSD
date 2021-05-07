# Figures & stats for FPO talk

# Simulate endemic scenarios across priors -------

# sub_cmd:=-t 8 -n 16 -jn camps -wt 1m -md \'gdal\' -ar \'1\' -cmd \'100\' -sn

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = FALSE)

# set up args and cluster if applicable ----
vacc_ind <- ifelse(!set_up$slurm, 1, as.numeric(arg[1]))
sim_vacc <- if(vacc_ind == 1) "fixed" else "random"
nsims <- ifelse(!set_up$slurm, 5, as.numeric(arg[2]))

cl <- make_cl(set_up$ncores)
register_cl(cl)
print(paste("Cluster size:", cl_size(cl)))

if(!set_up$slurm) fp <- here::here else fp <- cl_safe

# Dependencies
library(raster)
library(data.table)
library(sf)
library(tidyr)
library(dplyr)
library(magrittr)
library(simrabid) # remotes::install_github("mrajeev08/simrabid")
library(foreach)
library(iterators)
library(doRNG)
library(lubridate)
library(scoringRules)
library(igraph)
library(ggraph)
library(ggplot2)
library(cowplot)
library(patchwork)
library(glue)
library(ggtext)
library(readr)
library(graphlayouts)

# params for plotting
ppars <- plotpars()
# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
incursions <- fread(fp("analysis/out/incursions.csv"))
load("data/sd_case_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/get_observed_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")
source("R/run_mods.R")
source("R/score_preds.R")
source("R/conn_metrics.R")
source("R/figure_funs.R")

# viz reconstructed tree ---
consensus_tree <- fread("analysis/out/ttrees/consensus_links_ln0.95.csv")
I_gr <- setnames(consensus_tree[, c("id_progen", "id_case")][!is.na(id_progen)], c("from", "to"))
I_names <- data.table(names = unique(consensus_tree$id_case))

gr <- graph_from_data_frame(d = I_gr,
                            vertices = I_names,
                            directed = TRUE)
V(gr)$membership <- components(gr)$membership

intro_graph <- 
  ggraph(gr, layout= "stress")+
  geom_edge_link0(width = 0.2, colour = "grey")+
  geom_node_point(aes(col = factor(membership)), size = 0.3)+
  scale_color_discrete(guide = "none") + 
  theme_graph()

priors <- readRDS("analysis/out/priors.rds")
iota_prior <-
  ggplot(data =  data.frame(iota = priors$iota(1e5))) +
  geom_density(aes(x = iota), fill = "grey", color = NA) +
  cowplot::theme_minimal_hgrid() +
  labs(x = expression(iota), y = "Density")
r0_prior <- 
  ggplot(data =  data.frame(R0 = priors$R0(1e5))) +
  geom_density(aes(x = R0), fill = "grey", color = NA) +
  geom_vline(xintercept = 1.2, linetype = 2) +
  cowplot::theme_minimal_hgrid() +
  labs(x = expression(R[0]), y = "Density")
k_prior <- 
  ggplot(data =  data.frame(k = priors$R0(1e5))) +
  geom_density(aes(x = k), fill = "grey", color = NA) +
  cowplot::theme_minimal_hgrid() +
  labs(x = expression(k), y = "Density")

ggsave("analysis/figs/fpo/intro_graph.jpeg", intro_graph, height = 5, width = 7)
ggsave("analysis/figs/fpo/iota_prior.jpeg", iota_prior, height = 5, width = 7)
ggsave("analysis/figs/fpo/r0_prior.jpeg", r0_prior, height = 5, width = 7)
ggsave("analysis/figs/fpo/k_prior.jpeg", k_prior, height = 5, width = 7)

# monthly binomial reporting probability
report_mod <- data.frame(density = dbeta(seq(0.5, 1.01, by = 0.01), 
                                         shape2 = param_defaults$detect_beta, 
                                         shape1 = param_defaults$detect_alpha), 
                         detect = seq(0.5, 1.01, by = 0.01))
report_fig <-
  ggplot(data =  report_mod) +
  geom_density(aes(x = detect, weight = density), fill = "grey", color = NA) +
  cowplot::theme_minimal_hgrid() +
  labs(x = "Detection probability", y = "Density")
ggsave("analysis/figs/fpo/report_mod.jpeg", report_fig)

# viz model scale ---

out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

vill_scale <- 
  ggplot(sd_shapefile) + 
  geom_sf(aes(fill = pop_2002)) +
  scale_fill_distiller(trans = "log",
                       name = "Starting population", 
                       labels = scales::label_number(accuracy = 100)) +
  labs(title = "Village scale") +
  theme_map()

grids <- as.data.frame(out$rast, xy = TRUE)
grids$pop <- out$start_pop
grid_scale <- 
  ggplot() + 
  geom_raster(data = grids, aes(x = x, y = y, fill = pop)) +
  scale_fill_distiller(trans = "log", name = "Starting population", 
                       labels = scales::label_number(accuracy = 3), 
                       na.value = "white") +
  geom_sf(data = sd_shapefile, color = "grey", fill = NA) +
  labs(title = "1x1 km grid scale") +
  theme_map()

scale_fig <- grid_scale | vill_scale

ggsave("analysis/figs/fpo/scale_fig.jpeg")

# summarize model selection ----
mod_comp <- readRDS(fp("analysis/out/mod_comp/full_all.rds"))
mod_comp$modlookup %>%
  mutate(modindex = paste0("g", modindex)) %>%
  separate(modname, into = c("scale", "limits", "intros", "move"), sep = "_") %>%
  left_join(mod_comp$mod_predicted) %>%
  mutate(N = sum(votes.V1), prob = votes.V1/N) -> lookup

mod_ldas$lda_predicted %>%
  left_join(select(lookup, - N)) -> mod_ldas

ggplot(mod_ldas) +
  geom_point(aes(x = LD1, y = LD2, color = scale))

lookup %>%
  group_by(scale) %>%
  summarize(votes = sum(votes.V1), prob = votes/N[1]) %>%
  mutate(compare = "Scale", group = scale) -> scale_comp

lookup %>%
  group_by(intros) %>%
  summarize(votes = sum(votes.V1), prob = votes/N[1])  %>%
  mutate(compare = "Introductions", group = intros) -> intro_comp

lookup %>%
  group_by(move) %>%
  summarize(votes = sum(votes.V1), prob = votes/N[1]) %>%
  mutate(compare = "Movement", group = move) -> move_comp

lookup %>%
  group_by(limits) %>%
  summarize(votes = sum(votes.V1), prob = votes/N[1]) %>%
  mutate(compare = "Movement constraints", group = limits) -> limits_comp

mod_comp <- bind_rows(limits_comp, move_comp, intro_comp, scale_comp)
mod_comp %>% group_by(compare) %>% summarize(min = min(prob), max = max(prob)) -> mod_lines

xlabs <- c("Intros" = "Empirical vs.\n Estimated Introductions", 
           "Movement" = "Random Walk vs. \nKernel-based movement", 
           "Movement constraints" = "Constraints on \n movemements", 
           "Scale" = "Village vs. \n 1x1 km spatial scale")

mod_comp_fig <-
  ggplot(mod_comp) + 
  geom_point(aes(y = prob, x = compare, color = compare)) +
  geom_linerange(data = mod_lines, aes(x = compare, ymin = min, ymax = max, 
                                       color = compare)) +
  scale_x_discrete(labels = xlabs) +
  coord_flip() +
  scale_color_brewer(guide = "none", palette = "Dark2") +
  labs(x = "", y = "Posterior probability") +
  cowplot::theme_minimal_hgrid()
ggsave("analysis/figs/fpo/mod_comp.jpeg",mod_comp_fig, height = 6, width = 7)

lookup %>%
  group_by(scale) %>%
  arrange(desc(votes.V1)) %>%
  slice(1)

library(ggbeeswarm)
ggplot(lookup) + 
  geom_boxplot(aes(x = scale, y = prob)) + 
  labs(x = "Scale of mixing", y = "Posterior probability") 

# var importance plots ----
mod_comp <- readRDS(fp("analysis/out/mod_comp/full_all.rds"))
var_imp <- mod_comp$var_imp
var_imp %<>%
  left_join(select(ppars$summ_stat_labs, 
                   type_stat = type, 
                   summstat)) 


var_importance_plot_main <- 
  ggplot(filter(var_imp, !is.na(type_stat))) + 
  geom_point(aes(x = importance, y = reorder(summstat, importance), 
                 color = type_stat)) +
  scale_y_discrete(labels = ppars$slabs) +
  scale_color_manual(values = ppars$scols, name = "Type of \n summary stat") +
  labs(y = "") +
  ppars$theme_proj() +
  theme(axis.text.y = element_markdown()) 
ggsave("analysis/figs/fpo/var_imp_mod.jpeg", var_importance_plot_main, height = 6, width = 7)

# best sims


mod_comp <- readRDS("analysis/out/mod_comp/full_nostops.rds")
mod_comp_se <- readRDS("analysis/out/mod_comp/se_nostops.rds")
ggplot(mod_comp$err) + geom_line(aes(x = ntree, y = error.rate))

ggplot(mod_comp$lda_predicted) + 
  geom_hex(aes(x = LD1, y = LD2)) + 
  geom_point(data = mod_comp$lda_obs, aes(x = LD1, y = LD2), color = "red") + 
  facet_wrap(~modindex) 
ggplot(mod_comp$lda_predicted) + 
  geom_point(aes(x = LD1, y = LD2, color = factor(modindex)), alpha = 0.1) +
  geom_point(data = mod_comp$lda_obs, aes(x = LD1, y = LD2), color = "red")

test <- list.files("analysis/out/par_ests/")
out_all <- lapply(test, function(x) readRDS(paste0("analysis/out/par_ests/", x)))

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY = FALSE, fill = TRUE)
}

tt <-
  foreach(i = seq_len(length(out_all)), .combine = comb, 
        .multicombine = TRUE) %do% {
          
          out_all[[i]]
  
        }

preds <- tt$preds
preds[, c("scale", "movement", "incursions", "dispersal") := tstrsplit(modname, "_", fixed = TRUE)]
library(ggplot2)
setnafill(preds, fill = 0, cols = "nsim")

var_imp <- tt$var_imp
var_imp[, c("scale", "movement", "incursions", "dispersal") := tstrsplit(modname, "_", fixed = TRUE)]
ggplot(preds[param == "R0" & modname == "grid1x1_flim_estincs_kern"]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))

ggplot(preds[param == "k" & modname == "grid1x1_flim_estincs_kern"]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))

ggplot(preds[param == "iota"]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))


inds_joint <- preds[, .(check = all(V1 > 0)), by = c("sim_id", "nsim", "modname")][check == TRUE]

post_joint <-  preds[sim_id %in% inds_joint$sim_id]

ggplot(post_joint[param == "iota" & nsim == 0]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))

ggplot(post_joint[param == "k" & nsim == 0]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))


ggplot(post_joint[param == "R0" & nsim == 0]) + 
  geom_density(aes(x = resp), color = "NA", fill = "grey", alpha = 0.75) +
  geom_density(aes(x = resp, weight = V1, 
                   color = interaction(dispersal, movement), 
                   group = interaction(nsim, modname))) +
  facet_wrap(~interaction(scale, incursions))

library(tidyr)
preds %>% 
  pivot_wider(id_cols = c("modname", "nsim", "sim_id"), 
              names_from = "param", 
              values_from = c("V1", "resp")) -> preds_wider

preds_wider %>%
  filter(across(starts_with("V1")) > 0) -> preds_joint

ggplot(preds_joint) + 
  geom_density(aes(x = resp_k, weight = V1_k, color = factor(nsim))) +
  geom_density(data = preds[param == "k"], aes(x = resp), color = "black") +
  facet_wrap(~modname)


ggplot(preds_joint) + 
  geom_density(aes(x = resp_R0, weight = V1_R0, color = factor(nsim))) +
  geom_density(data = preds[param == "R0"], aes(x = resp), color = "black") +
  facet_wrap(~modname)


ggplot(preds_joint) + 
  geom_density(aes(x = resp_iota, weight = V1_iota, color = factor(nsim))) +
  geom_density(data = preds[param == "iota"], aes(x = resp), color = "black") +
  facet_wrap(~modname)

ggplot() + 
  geom_point(data = preds_wider, 
             aes(x = resp_iota, y = resp_k), color = "grey", alpha = 0.25) +
  geom_point(data = preds_joint, 
             aes(x = resp_iota, y = resp_k), color = "black") +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  facet_wrap(~modname) 

ggplot() + 
  geom_point(data = preds_wider, 
             aes(x = resp_iota, y = resp_R0), color = "grey", alpha = 0.25) +
  geom_point(data = preds_joint, 
             aes(x = resp_iota, y = resp_R0), color = "black") +
  facet_wrap(~modname)

ggplot() + 
  geom_point(data = preds_wider, 
             aes(x = resp_R0, y = resp_k), color = "grey", alpha = 0.25) +
  geom_point(data = preds_joint, aes(x = resp_R0, y = resp_k), color = "black") +
  facet_wrap(~modname)

# do the joint posteriors look any better?
mod_comp <- readRDS("analysis/out/archive/mod_comp/full_grid.rds")
lda_vals <- mod_comp$lda_predicted

ggplot() + 
  geom_point(data = lda_vals,
             aes(x = LD1, y = LD2), alpha = 0.1) +
  geom_point(data = mod_comp$lda_obs,
             aes(x = LD1, y = LD2), alpha = 1, color = "red") +
  facet_wrap(~modindex)

# example sims (w / 1e4)
test_scores <- fread("analysis/out/archive/post_sims/grid1x1_plim_estincs_kern_ts_scores.csv")
case_ts <- fread("data/case_ts.csv")
test_sims <- fread("analysis/out/archive/post_sims/grid1x1_plim_estincs_kern_sims.csv")

ggplot(data = test_scores, aes(x = cal_month)) + 
  geom_line(aes(y = crps_score, color = vacc_type)) +
  facet_wrap(~ post_type)

ggplot(data = test_scores, aes(x = cal_month)) + 
  geom_ribbon(aes(x = cal_month, ymin = min_times_0.5, ymax = max_times_0.5), 
              fill = "grey") +
  geom_line(aes(y = best_dat_1), color = "blue") +
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "purple") +
  facet_grid(vacc_type ~ post_type, scales = "free")

ggplot(data = test_sims) + 
  geom_hex(aes(x = cal_month, y = I_ts, color=..density.., fill=..density..)) + 
  geom_line(data = case_ts, aes(x = cal_month, y = cases), col = "purple") +
  facet_wrap(post_type ~ vacc_type, scales = "free")

ggplot(data = test_sims) + 
  geom_hex(aes(x = cal_month, y = cov, color=..density.., fill=..density..)) + 
  facet_wrap(post_type ~ vacc_type, scales = "free")

ggplot(data = test_sims) + 
  geom_line(aes(x = cal_month, y = 1- cov, group = sim)) + 
  facet_wrap(post_type ~ vacc_type, scales = "free")

ggplot(data = test_sims, aes(y = factor(cal_month), x = I_ts)) + 
  geom_density_ridges(color = NA, fill = "blue") + 
  geom_point(data = case_ts, aes(y = cal_month, x = cases), color = "grey", alpha = 0.2) +
  facet_wrap(post_type ~ vacc_type, scales = "free")

out <- test_sims[, .(x = density(I_ts)$x, density = density(I_ts)$y), 
                 by = c("cal_month", "post_type", "vacc_type", "modname")]

ggplot(out) + geom_tile(aes(x = cal_month, y = x, fill = density)) + facet_wrap(post_type ~ vacc_type)

test_sims[, .(med = median(incs_ts)), by = c("cal_month", "post_type", "vacc_type", "modname")]

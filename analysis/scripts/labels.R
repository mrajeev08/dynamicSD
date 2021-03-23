library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)

# summ stat names & tables
# Tribble with summary and classification of summary stats
summ_stat_labs <- 
  tribble(~summstat, ~lab, ~descr, ~type, 
          "max_I", "Max I", "Maximum monthly number of cases", "Temporal",
          "median_I", "Median I", "Maximum monthly number of cases", "Temporal",
          "mean_I", "Mean I", "Mean monthly number of cases", "Temporal",
          "kb_temp", "Temporal KB", "Kullback-Leibler Divergence between observed and simulated distribution of cases temporally", "Temporal", 
          "kb_spat", "Spatial KB", "Kullback-Leibler Divergence between observed and simulated distribution of cases spatially", "Spatial", 
          "mean_dist_4wks", "Mean distance 4 weeks", "Mean distance between cases within an 4 week window", "Spatiotemporal",
          "mean_dist_8wks", "Mean distance 8 weeks","Mean distance betweeen cases within an 8 week window", "Spatiotemporal",
          "mean_dist_4wks_norm", "Mean distance 4 weeks <br> normalized", "Mean distances of cases within an 4 week window, normalized by the average distance between all cases", "Spatiotemporal",
          "mean_dist_8wks_norm", "Mean distance 8 weeks <br> normalized", "Mean distances of cases within an 8 week window, normalized by the average distance between all cases", "Spatiotemporal", 
          "spat_rmse", "Spatial RMSE", "Root mean squared error between the observed vs. simulated density of case locations", "Spatial",
          "spat_loss", "Spatial loss" , "Loss value of the absolute difference between observed and simulated cases locations", "Spatial",
          "temp_rmse", "Temporal RMSE", "Root mean squared error between the observed vs. simulated numbers of monthly cases", "Spatial",
          "temp_loss", "Temporal loss", "Loss value of the absolute difference between observed and simulated monthly cases", "Temporal",
          "mean_dist_all", "Mean distance all cases", "Average distance between all cases", "Spatial",
          NA, "ACF Monthly 1-6", "Autocorrelation values for monthly timeseries at 1 - 6 month lag", "Temporal",
          NA, "ACF Weekly 1-10", "Autocorrelation values for weekly timeseries at 1 - 10 week lag", "Temporal",
          "acf_month1", "ACF monthly-lag 1", NA, "Temporal",
          "acf_month2", "ACF monthly-lag 2", NA, "Temporal",
          "acf_month3", "ACF monthly-lag 3", NA, "Temporal",
          "acf_month4", "ACF monthly-lag 4", NA, "Temporal",
          "acf_month5", "ACF monthly-lag 5", NA, "Temporal",
          "acf_month6", "ACF monthly-lag 6", NA, "Temporal",
          "acf_week1", "ACF weekly-lag 1", NA, "Temporal",
          "acf_week2", "ACF weekly-lag 2", NA, "Temporal",
          "acf_week3", "ACF weekly-lag 3", NA, "Temporal",
          "acf_week4", "ACF weekly-lag 4", NA, "Temporal",
          "acf_week5", "ACF weekly-lag 5", NA, "Temporal",
          "acf_week6", "ACF weekly-lag 6", NA, "Temporal",
          "acf_week7", "ACF weekly-lag 7", NA, "Temporal",
          "acf_week8", "ACF weekly-lag 8", NA, "Temporal",
          "acf_week9", "ACF weekly-lag 9", NA, "Temporal",
          "acf_week10", "ACF weekly-lag 10",  NA, "Temporal",
          "LD1", NA, NA, "Model based",
          "LD2", NA, NA, "Model based",
          "LD3", NA, NA, "Model based",
          "LD4", NA, NA, "Model based",
          "LD5", NA, NA, "Model based",
          "LD6", NA, NA, "Model based",
          "LD7", NA, NA, "Model based",
          "LD8", NA, NA, "Model based",
          "LD9", NA, NA, "Model based",
          "LD10", NA, NA, "Model based",
          "LD11", NA, NA, "Model based")
write_csv(summ_stat_labs, "analysis/out/summ_stats_labs.csv")

# model names 
fpath <- "analysis/out/mod_comp/"
fp <- function(x) paste0(fpath, x)

# combine files
test <- list.files(fpath)
out_names <- lapply(test, function(x) {
  reslt <- readRDS(fp(x))$modlookup
  reslt[, type := x]
  reslt
})

# distinct out names
out_names <- rbindlist(out_names, fill = TRUE)
out_names <- distinct(out_names, modname, modindex) # 24 candidate models
out_names %<>%
  mutate(scale = case_when(grepl("grid", modname) ~ "1x1 km",
                           grepl("vill", modname) ~ "Village"), 
         move = case_when(grepl("kern", modname) ~ "Kernel",
                          grepl("seq", modname) ~ "Random Walk"), 
         incs = case_when(grepl("fixedincs", modname) ~ "Fixed",
                          grepl("estincs", modname) ~ "Estimated"), 
         limits = case_when(grepl("flim", modname) ~ "Full",
                            grepl("plim", modname) ~ "Patch", 
                            grepl("nolim", modname) ~ "None"), 
         modindex = paste0("g", modindex))
write_csv(out_names, "analysis/out/model_labs.csv")

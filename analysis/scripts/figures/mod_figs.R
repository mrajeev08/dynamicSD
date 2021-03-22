# Main and supplementary figures and tables for model selection

library(data.table)
library(ggplot2)
library(foreach)
source("R/utils.R")

# path to files
fpath <- "analysis/out/archive/mod_comp/"
fp <- function(x) paste0(fpath, x)

# combine files
test <- list.files(fpath)
out_all <- lapply(test, function(x) {
    reslt <- readRDS(fp(x)) 
    reslt$name <- x
    reslt
  }
)
modcomp_list <-
  foreach(i = seq_len(length(out_all)), .combine = comb, 
          .multicombine = TRUE) %do% {
            browser()
            tnow <- out_all[[i]]
            tnow <- lapply(tnow, append_col, col_name = "type", val = tnow$name)
  }

# oob errors
ggplot(mod_comp_grid$err) + 
  geom_line(aes(x = ntree, y = error.rate)) +
  geom_line(data = mod_comp_se$err, aes(x = ntree, y = error.rate, 
                                        group = nsim, color = factor(nsim)))

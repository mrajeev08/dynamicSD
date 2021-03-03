## sim.ibm testing code
# ## Unhash for testing! -------------------------------------------------------------------------
grid = as.matrix(SD_raster); data = grid_data; vacc = vacc_mat;
I_seed = 2; start_vacc = 0.2;
R_0 = 1.05; k = 0.2;
inc_rate = 1;
# weekly number of incursions and scaling factor
sigma = get.prob(rate = 7/22.3, step = 1); # weekly rate to prob exp -> inf
births = get.prob(rate = 0.44 + (grid_data$growth - 1), step = 52); # annual birth rate to prob
mu = get.prob(rate = 0.44, step = 52); # annual death rate to prob
ntimes = 100; # maximum time step
nu = get.prob(rate = 0.33, step = 52);
p_revacc = 0.5; p_vacc = 1; nlocs = nrow(grid_data); cell_id = grid_data$cell_id;
# annual waning rate to prob + probability of revaccination
dispersalShape = 0.3484; dispersalScale = 41.28/100;
x_topl = bbox(SD_raster)[1, "min"];
y_topl = bbox(SD_raster)[2, "max"]; res_m = res(SD_raster)[1]; return_coords = FALSE;
# No vacc scenario
# vacc_mat[vacc_mat > 0 ] <- 0
# start_vacc = 0

## Testing sim.bites
secondaries = I_coords_now$secondaries; x_coord = I_coords_now$x_coord;
y_coord = I_coords_now$y_coord; tstep = t; ids = I_coords_now$ID; counter; grid; res_m; 
cells_pop; dispersalShape; dispersalScale; x_topl; y_topl

## Testing sim vacc
locid = data$villcode; Sdogs = S[, 1]; Ndogs = N[, 1]; Vdogs = V[, 1];
vill_vacc = vacc[, 77]; p_revacc = p_revacc; p_allocate = runif(nrow(S), 0, 1);

library(microbenchmark)
check <- as.matrix(grid)
microbenchmark({cell <- grid[ceiling(-(y_new - y_topl)/res_m), ceiling((x_new - x_topl)/res_m)]},
               {cell <- check[ceiling(-(y_new - y_topl)/res_m), ceiling((x_new - x_topl)/res_m)]}, times = 100)

check <- 1:100
microbenchmark({
  within <- ifelse(cell %in% cells_pop, 1, 0)
  }, {
    cells_pop <- cell_id[which(N[, 1] > 0)] ## populated cells at beginning
    inds <- t((sapply(cells_pop, function(x) which(grid == x, arr.ind = TRUE))))
    row_inds <- inds[, 1]
    col_inds <- inds[, 2]
    
}, times = 100)

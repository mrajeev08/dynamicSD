## Getting adjacency graph for shapefile
get.adjgraph.shape <- function(shapefile){
  adj_vill <- poly2nb(shapefile)
  adj_mat <- nb2mat(adj_vill, style="B", zero.policy = TRUE)
  adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  return(adj_graph)
}

## Getting adjacency graph for raster
get.adjgraph.raster <- function(raster, cell_id = grid_data$cell_id, 
                                row_id = 1:nrow(grid_data)){
  adj_pairs <- as.data.table(adjacent(raster, cells = 1:ncell(raster),
                                      directions = 8, pairs = TRUE))
  adj_mat <- sparseMatrix(i = adj_pairs$from, j = adj_pairs$to,  x = 1, use.last.ij = TRUE)
  adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  V(adj_graph)$row_id <- row_id[match(V(adj_graph), cell_id)]
  adj_graph <- delete.vertices(adj_graph, which(is.na(V(adj_graph)$row_id)))
  return(adj_graph)
}

## apply this function to each column in cov mat (needs to be merged with all grid cells!)
get.connmetric <- function(x, adj_graph = adj_graph_rast, threshold = 0.8){
  V(adj_graph)$sus <- x
  adj_graph <- delete.vertices(adj_graph, which(V(adj_graph)$sus <= threshold))
  npatches <- components(adj_graph)$csize
  conn <- sum(npatches^2)
  return(conn)
}
# metric of connected patches with susceptibility greater than 0.8 (or threshold)
# or if cov_mat then, metric of connected patches with coverage greater than 0.8 (or threshold)
adj_vill <- poly2nb(SD_shape)
adj_mat <- nb2mat(adj_vill, style="B", zero.policy = TRUE)
adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
vill_vacc <- check[[7]]
vill_vacc <- vill_vacc[order(match(rownames(vill_vacc), SD_shape$VILLCODE)), ]
SD_shape$VILLCODE == as.factor(rownames(vill_vacc)) ## to check order is right
quartz()
plot(apply(vill_vacc, 2, get.connmetric, adj_graph = adj_graph, threshold = 0.2), type = "l")

## connectivity raster
adj_pairs <- as.data.table(adjacent(SD_raster, cells = 1:ncell(SD_raster),
                                    directions = 8, pairs = TRUE))
adj_mat <- sparseMatrix(i = adj_pairs$from, j = adj_pairs$to,  x = 1, use.last.ij = TRUE)
adj_graph_rast <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
cell_id <- grid_data$cell_id
row_id <- 1:nrow(grid_data)
V(adj_graph_rast)$row_id <- row_id[match(V(adj_graph_rast), cell_id)]
adj_graph_rast <- delete.vertices(adj_graph_rast, which(is.na(V(adj_graph_rast)$row_id)))
cov_mat <- 1 - S/N
quartz()
plot(apply(cov_mat, 2, get.connmetric, adj_graph = adj_graph_rast, threshold = 0.20)/nrow(cov_mat)^2, type = "l")

## percolation metrics!
## Getting connectivity metrics
## @ vill level
library(igraph)
library(spdep)
library(Matrix)
adj_vill <- poly2nb(SD_shape)
adj_mat <- nb2mat(adj_vill, style="B", zero.policy = TRUE)
adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
V(adj_graph)$sus <- 1 - vill_vacc[, 2] ## need to double check this!
adj_graph <- delete.vertices(adj_graph, which(V(adj_graph)$sus >= 0.8))
npatches <- components(adj_graph)$csize
conn <- sum(npatches^2)

## gridded version--apply to cov_mat
adj_pairs <- as.data.table(adjacent(SD_raster, cells = 1:ncell(SD_raster),
                                    directions = 8, pairs = TRUE))
adj_mat <- sparseMatrix(i = adj_pairs$from, j = adj_pairs$to,  x = 1, use.last.ij = TRUE)
adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
cell_id <- grid_data$cell_id
row_id <- 1:nrow(grid_data)
V(adj_graph)$row_id <- row_id[match(V(adj_graph), cell_id)]
adj_graph <- delete.vertices(adj_graph, which(is.na(V(adj_graph)$row_id)))
## need to double check this!
## Way to match between cell ids and row ids do this outside of cov mat?
V(adj_graph)$sus <- S[row_id[match(V(adj_graph), cell_id)], 2]/
  N[row_id[match(V(adj_graph), cell_id)], 2]  ## need to double check this!
adj_graph <- delete.vertices(adj_graph, which(is.na(V(adj_graph)$sus), 
                                              which(V(adj_graph)$sus <= 0.8)))
npatches <- components(adj_graph)$csize
conn <- sum(npatches^2)
max_conn <- nrow(grid_data)^2

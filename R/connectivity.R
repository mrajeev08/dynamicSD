## Getting adjacency graph for shapefile
library(rgdal)
library(spdep)
library(igraph)
SD_shape <- readOGR("~/Documents/Projects/dynamics_SD/data/SD_shape/Serengeti_villages UTM_region.shp")
get.adjgraph.shape <- function(shapefile){
  adj_vill <- poly2nb(shapefile)
  adj_mat <- nb2mat(adj_vill, style="B", zero.policy = TRUE)
  adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  return(adj_graph)
}

adj_graph <- get.adjgraph.shape(SD_shape)
saveRDS(adj_graph, "~/Documents/Projects/rabies_sfunk/data/vill_graph.rds")

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
get.connmetric <- function(x, adj_graph = adj_graph_rast, threshold = 0.2){
  V(adj_graph)$cov <- x
  adj_graph <- delete.vertices(adj_graph, which(V(adj_graph)$cov <= threshold))
  npatches <- components(adj_graph)$csize
  conn <- sum(npatches^2)
  return(conn)
}

## Getting connectivity metrics
## @ vill level
library(igraph)
library(spdep)
library(Matrix)
# metric of connected patches with susceptibility greater than 0.8 (or threshold)
# or if cov_mat then, metric of connected patches with coverage greater than 0.8 (or threshold)
library(rgdal)
SD_shape <- readOGR("~/Documents/Projects/dynamics_SD/data/SD_shape/Serengeti_villages UTM_region.shp")

adj_vill <- poly2nb(SD_shape, row.names = SD_shape$VILLCODE)
adj_mat <- nb2mat(adj_vill, style="B", zero.policy = TRUE)
adj_graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE, add.rownames = "villcodes")
saveRDS(adj_graph, "~/Documents/Projects/rabies_sfunk/data/vill_graph.rds")

vill_vacc <- check[[7]]

vill_vacc <- vill_vacc[order(match(rownames(vill_vacc), V(adj_graph)$villcodes)), ]
V(adj_graph)$villcodes == rownames(vill_vacc) ## to check order is right

## For writing to seb's folder
# vill_vacc[, 1] <- 0.2
# vill_vacc <- vill_vacc[, 1:676]
# vill_vacc <- vill_vacc[, seq(1, 676, by = 4)]
# vill_vacc <- vill_vacc[, 1:157]
# write.csv(vill_vacc, "~/Documents/Projects/rabies_sfunk/data/vill_vaccmat.csv", row.names = TRUE)

# write.csv(1- colSums(S)/colSums(N), "cov_ts.csv", row.names = FALSE)

quartz()
plot(apply(vill_vacc, 2, get.connmetric, adj_graph = adj_graph, threshold = 0.2)/75^2, 
     type = "l", 
     bty = "n", xlab = "Time (weeks)", ylab = "Connectivity metric", col = "red")
points( 1- colSums(S)/colSums(N), type = "l", col = "blue")
     

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
points( 1- colSums(S)/colSums(N), type = "l", col = "blue")

## percolation metrics!


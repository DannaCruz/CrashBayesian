suppressPackageStartupMessages({
  library(igraph)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
repo <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
data_dir <- file.path(repo, "data", "processed")
out_dir <- file.path(repo, "model-output", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

covs <- readRDS(file.path(data_dir, "model_covariates_bogota.rds"))
g_seg <- readRDS(file.path(data_dir, "graph_segments_car.rds"))
g_int <- readRDS(file.path(data_dir, "graph_intersections_car.rds"))

graph_data <- function(graph, n_expected) {
  stopifnot(vcount(graph) == n_expected)
  edges <- as_edgelist(graph, names = FALSE)
  membership <- components(graph)$membership
  ids <- seq_len(max(membership))
  nodes <- unlist(lapply(ids, function(z) which(membership == z)), use.names = FALSE)
  sizes <- tabulate(membership, nbins = length(ids))
  list(
    N_edges = nrow(edges), node1 = as.integer(edges[, 1]),
    node2 = as.integer(edges[, 2]), N_components = length(sizes),
    component_start = as.integer(c(1, 1 + cumsum(sizes))),
    component_nodes = as.integer(nodes)
  )
}

base_data <- function(y, x, offset) list(
  N = length(y), K = ncol(x), Y = as.integer(y),
  X = unname(split(unclass(x), row(unclass(x)))),
  log_offset = as.numeric(offset)
)

seg <- base_data(covs$Y_segment, covs$X_segment, covs$log_offset_segment)
int <- base_data(covs$Y_intersection, covs$X_intersection, covs$log_offset_intersection)
stopifnot(seg$N == 7856L, int$N == 7992L,
          sum(seg$Y) == 18158L, sum(int$Y) == 47953L)

write_json(seg, file.path(out_dir, "segment_nonspatial.json"), auto_unbox = TRUE, digits = NA)
write_json(int, file.path(out_dir, "intersection_nonspatial.json"), auto_unbox = TRUE, digits = NA)
write_json(c(seg, graph_data(g_seg, seg$N)), file.path(out_dir, "segment_car.json"), auto_unbox = TRUE, digits = NA)
write_json(c(int, graph_data(g_int, int$N)), file.path(out_dir, "intersection_car.json"), auto_unbox = TRUE, digits = NA)
message("Stan data written to ", out_dir)


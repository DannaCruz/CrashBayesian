suppressPackageStartupMessages({
  library(sf); library(sfnetworks); library(igraph); library(Matrix); library(RSpectra)
})
root <- normalizePath('outputs/modelos_espectrales', mustWork = TRUE)
src <- Sys.getenv('CRASHBAYESIAN_SOURCE_DIR')
if (!nzchar(src)) stop('Set CRASHBAYESIAN_SOURCE_DIR to the original RData directory')
set.seed(20261120)
audits <- list(); cov_audits <- list()
for (city in c('Barcelona', 'Bogota')) {
  env <- new.env()
  filename <- if (city == 'Bogota') 'Bogota_network_ready_for_models.RData' else 'Eixample_network_ready_for_models.RData'
  source_file <- file.path(src, filename)
  load(source_file, envir = env)
  roads <- if (city == 'Bogota') env$roads_ls else env$roads_sf
  net_edges <- st_as_sf(env$net, 'edges')
  stopifnot(nrow(roads) == ecount(env$g),
    identical(st_as_binary(st_geometry(roads)), st_as_binary(st_geometry(net_edges))))
  ep <- ends(env$g, E(env$g), names = FALSE)
  stopifnot(!any(ep[, 1] == ep[, 2]))
  comp <- components(env$g)
  edge_comp <- comp$membership[ep[, 1]]
  chosen <- which.max(tabulate(edge_comp))
  keep <- which(edge_comp == chosen)
  # Preserve original road-segment order, including legitimate parallel roads.
  r <- roads[keep, ]; ep <- ep[keep, , drop = FALSE]
  nodes <- sort(unique(as.vector(ep)))
  ep <- matrix(match(ep, nodes), ncol = 2)
  incidence <- sparseMatrix(i = as.vector(ep), j = rep(seq_len(nrow(r)), 2),
                           x = 1, dims = c(length(nodes), nrow(r)))
  A <- crossprod(incidence); diag(A) <- 0
  A@x[] <- 1; A <- drop0(A)
  deg <- as.numeric(rowSums(A)); stopifnot(all(deg > 0), isSymmetric(A))
  S <- Diagonal(x = 1 / sqrt(deg)) %*% A %*% Diagonal(x = 1 / sqrt(deg))
  eig <- eigs_sym(S, k = 61, which = 'LA', opts = list(tol = 1e-10, maxitr = 20000))
  ord <- order(eig$values, decreasing = TRUE)
  eig$values <- eig$values[ord]; eig$vectors <- eig$vectors[, ord, drop = FALSE]
  # Fix arbitrary eigenvector signs deterministically for reproducible maps.
  for (j in seq_len(ncol(eig$vectors))) {
    anchor <- which.max(abs(eig$vectors[, j]))
    if (eig$vectors[anchor, j] < 0) eig$vectors[, j] <- -eig$vectors[, j]
  }
  residual <- vapply(seq_len(61), function(j)
    sqrt(sum((S %*% eig$vectors[, j] - eig$values[j] * eig$vectors[, j])^2)), numeric(1))
  stopifnot(max(residual) < 1e-6, abs(eig$values[1] - 1) < 1e-6)
  candidate <- eig$vectors[, 2:61, drop = FALSE]
  B60 <- scale(candidate, center = TRUE, scale = FALSE)
  B60 <- sweep(B60, 2, sqrt(colMeans(B60^2)), '/')
  B <- B60[, 1:20, drop = FALSE]
  lambda <- eig$values[2:21]
  counts <- expm1(r$Y)
  stopifnot(all(is.finite(counts)), max(abs(counts - round(counts))) < 1e-7,
            all(counts >= 0))
  y <- as.integer(round(counts))
  length_km <- as.numeric(st_length(r)) / 1000
  stopifnot(all(is.finite(length_km)), all(length_km > 0))
  highway <- factor(r$highway, levels = c('residential','tertiary','secondary','primary'))
  stopifnot(!anyNA(highway))
  lanes_num <- suppressWarnings(as.numeric(r$lanes))
  invalid_lanes <- !is.na(r$lanes) & (is.na(lanes_num) | lanes_num < 1 |
                          lanes_num > 10 | lanes_num != round(lanes_num))
  lanes_num[invalid_lanes] <- NA_real_
  lanes_group <- ifelse(is.na(lanes_num), 'unknown',
                        ifelse(lanes_num == 1, 'one', ifelse(lanes_num == 2, 'two', 'three_plus')))
  lanes <- factor(lanes_group, levels = c('two','one','three_plus','unknown'))
  surface_group <- ifelse(is.na(r$surface) | r$surface == '', 'unknown',
                          ifelse(r$surface == 'asphalt', 'asphalt','other_recorded'))
  surface <- factor(surface_group, levels = c('asphalt','other_recorded','unknown'))
  cov_df <- data.frame(highway, lanes, surface)
  X <- model.matrix(~ highway + lanes + surface, data = cov_df)[, -1, drop = FALSE]
  stopifnot(nrow(X) == length(y), qr(cbind(1, X))$rank == ncol(X) + 1)
  r$source_segment_row <- keep; r$analysis_segment_id <- seq_len(nrow(r))
  r$count_recovered <- y; r$length_km_verified <- length_km
  r$lanes_group <- lanes; r$surface_group <- surface
  xy <- st_coordinates(st_centroid(st_geometry(r)))
  # Fixed geographic folds, built without looking at the response. Saved for later blocked validation.
  set.seed(20261120)
  blocks <- kmeans(scale(xy), centers = 5, nstart = 30)$cluster
  input <- list(city = city, y = y, offset = log(length_km), X = X, B = B,
    lambda = lambda, B60 = B60, lambda60 = eig$values[2:61], U_raw = candidate,
    eigenvalues = eig$values, eigen_residual = residual, adjacency = A,
    segment_endpoints = ep, segments = r, geographic_fold = blocks,
    source_file = source_file, source_md5 = unname(tools::md5sum(source_file)),
    response_note = 'Recovered integer counts from expm1(roads$Y); event-level provenance not validated.',
    spectral_note = 'Raw normalized-adjacency eigenvectors; fitted columns are centered and RMS-scaled. Perron mode omitted, as in manuscript.',
    offset_note = 'Road length in km, over each dataset observation period; no traffic-volume exposure.',
    covariate_note = 'OSM fields; dates not harmonized. Missing categories are nuisance labels, not absence of infrastructure.')
  saveRDS(input, file.path(root, 'datos', paste0(city, '_preparado.rds')))
  flat <- data.frame(analysis_segment_id = r$analysis_segment_id, source_segment_row = keep,
    osm_id = r$osm_id, count_recovered = y, length_km = length_km,
    cov_df, geographic_fold = blocks, x = xy[,1], y_coord = xy[,2])
  write.csv(flat, file.path(root, 'datos', paste0(city, '_segmentos.csv')), row.names = FALSE)
  audits[[city]] <- data.frame(city, original_segments = nrow(roads), retained_segments = nrow(r),
    excluded_segments = nrow(roads) - nrow(r), retained_intersections = length(nodes),
    stored_Y_edges_length = length(env$Y_edges), response_used = 'roads$Y, inverted with expm1',
    retained_count_total = sum(y), zero_segments = sum(y == 0), length_km = sum(length_km),
    invalid_lanes = sum(invalid_lanes), maximum_eigen_residual = max(residual),
    source_md5 = input$source_md5)
  cov_audits[[city]] <- do.call(rbind, lapply(names(cov_df), function(v) {
    tab <- table(cov_df[[v]])
    data.frame(city, variable = v, level = names(tab), n = as.integer(tab))
  }))
  write.csv(data.frame(index = 1:61, eigenvalue = eig$values, residual = residual),
            file.path(root,'resultados',paste0(city,'_espectro.csv')), row.names = FALSE)
  cat(city, 'ready:',length(y),'segments,',sum(y),'counts\n')
  rm(env); gc()
}
write.csv(do.call(rbind, audits), file.path(root,'resultados','auditoria_datos.csv'), row.names = FALSE)
write.csv(do.call(rbind, cov_audits), file.path(root,'resultados','covariables_niveles.csv'), row.names = FALSE)

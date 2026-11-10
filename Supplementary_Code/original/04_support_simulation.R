#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MASS))

set.seed(20260824)

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

make_grid <- function(side = 10L) {
  nodes <- expand.grid(x = seq_len(side), y = seq_len(side))
  id <- matrix(seq_len(side * side), nrow = side, ncol = side)
  horizontal <- cbind(as.vector(id[, -side]), as.vector(id[, -1]))
  vertical <- cbind(as.vector(id[-side, ]), as.vector(id[-1, ]))
  edges <- rbind(horizontal, vertical)
  list(nodes = nodes, edges = edges)
}

safe_nb <- function(formula, data) {
  tryCatch(glm.nb(formula, data = data, control = glm.control(maxit = 100)),
           error = function(e) glm(formula, data = data, family = quasipoisson()))
}

one_rep <- function(rep_id, side = 10L, size = 2) {
  g <- make_grid(side)
  nn <- nrow(g$nodes)
  ne <- nrow(g$edges)
  node_x <- scale(rnorm(nn))[, 1]
  edge_x <- scale((node_x[g$edges[, 1]] + node_x[g$edges[, 2]]) / 2 + rnorm(ne, sd = 0.5))[, 1]
  edge_mid_x <- (g$nodes$x[g$edges[, 1]] + g$nodes$x[g$edges[, 2]]) / 2
  edge_mid_y <- (g$nodes$y[g$edges[, 1]] + g$nodes$y[g$edges[, 2]]) / 2
  edge_spatial <- scale(sin(pi * edge_mid_x / side) + cos(pi * edge_mid_y / side))[, 1]
  beta <- 0.8
  eta_edge <- -0.2 + beta * edge_x + 0.7 * edge_spatial
  y_edge <- rnbinom(ne, mu = exp(eta_edge), size = size)

  native <- safe_nb(y ~ x + sx + sy,
                    data.frame(y = y_edge, x = edge_x,
                               sx = scale(edge_mid_x)[, 1], sy = scale(edge_mid_y)[, 1]))
  pred_native <- predict(native, type = "response")

  # Preserve the total event count while relocating every edge event to one of
  # its endpoints, which mimics changing from edge to intersection support.
  y_node <- integer(nn)
  for (e in seq_len(ne)) {
    left <- rbinom(1, y_edge[e], 0.5)
    y_node[g$edges[e, 1]] <- y_node[g$edges[e, 1]] + left
    y_node[g$edges[e, 2]] <- y_node[g$edges[e, 2]] + y_edge[e] - left
  }
  node_fit <- safe_nb(y ~ x + sx + sy,
                      data.frame(y = y_node, x = node_x,
                                 sx = scale(g$nodes$x)[, 1], sy = scale(g$nodes$y)[, 1]))
  node_mu <- predict(node_fit, type = "response")
  mapped_edge <- (node_mu[g$edges[, 1]] + node_mu[g$edges[, 2]]) / 2
  mapped_edge <- mapped_edge * sum(y_edge) / sum(mapped_edge)

  data.frame(
    replicate = rep_id,
    native_RMSE = sqrt(mean((y_edge - pred_native)^2)),
    changed_support_RMSE = sqrt(mean((y_edge - mapped_edge)^2)),
    native_Spearman = suppressWarnings(cor(y_edge, pred_native, method = "spearman")),
    changed_support_Spearman = suppressWarnings(cor(y_edge, mapped_edge, method = "spearman")),
    native_beta_error = abs(coef(native)["x"] - beta)
  )
}

res <- do.call(rbind, lapply(seq_len(250), one_rep))
write.csv(res, file.path(root, "results/support_simulation_replicates.csv"), row.names = FALSE)

summary <- data.frame(
  metric = c("RMSE", "Spearman correlation"),
  native_median = c(median(res$native_RMSE), median(res$native_Spearman)),
  changed_support_median = c(median(res$changed_support_RMSE), median(res$changed_support_Spearman)),
  native_better_percent = c(
    mean(res$native_RMSE < res$changed_support_RMSE) * 100,
    mean(res$native_Spearman > res$changed_support_Spearman) * 100
  )
)
write.csv(summary, file.path(root, "results/support_simulation_summary.csv"), row.names = FALSE)

cat("Support simulation complete.\n")
print(summary)


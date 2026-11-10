#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(MASS)
})

set.seed(20260824)

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

acc <- readRDS(file.path(root, "data/processed/accidents_classified.rds"))
roads <- readRDS(file.path(root, "data/processed/roads_domain_counts.rds"))
nodes <- readRDS(file.path(root, "data/processed/nodes_domain_sf.rds"))
covs <- readRDS(file.path(root, "data/processed/model_covariates_bogota.rds"))

train_years <- 2007:2019
test_years <- 2020:2024

acc$unit_id <- NA_integer_
is_int <- acc$support == "intersection"
acc$unit_id[is_int] <- acc$nearest_node[is_int]
is_seg <- acc$support == "segment"
acc$unit_id[is_seg] <- st_nearest_feature(acc[is_seg, ], roads)

count_period <- function(support_name, years, n_units) {
  z <- acc$support == support_name & acc$ANO_OCURRENCIA_ACC %in% years
  tabulate(acc$unit_id[z], nbins = n_units)
}

score_predictions <- function(y, mu, model, support) {
  mu <- pmax(as.numeric(mu), 1e-10)
  data.frame(
    support = support,
    model = model,
    MAE = mean(abs(y - mu)),
    RMSE = sqrt(mean((y - mu)^2)),
    Spearman = suppressWarnings(cor(y, mu, method = "spearman")),
    top_decile_capture = sum(y[mu >= quantile(mu, 0.9)]) / sum(y),
    mean_observed = mean(y),
    mean_predicted = mean(mu),
    stringsAsFactors = FALSE
  )
}

fit_support <- function(support) {
  if (support == "segment") {
    n <- nrow(roads)
    X <- as.data.frame(covs$X_segment)
    base_offset <- covs$log_offset_segment
  } else {
    n <- nrow(nodes)
    X <- as.data.frame(covs$X_intersection)
    base_offset <- covs$log_offset_intersection
  }
  names(X) <- paste0("x", seq_len(ncol(X)))
  y_train <- count_period(support, train_years, n)
  y_test <- count_period(support, test_years, n)
  d <- cbind(data.frame(y = y_train, off = base_offset + log(length(train_years))), X)
  full <- glm.nb(y ~ . - off + offset(off), data = d, control = glm.control(maxit = 100))
  exposure <- glm.nb(y ~ 1 + offset(off), data = d, control = glm.control(maxit = 100))
  nd <- cbind(data.frame(off = base_offset + log(length(test_years))), X)
  mu_full <- predict(full, newdata = nd, type = "response")
  mu_exposure <- predict(exposure, newdata = nd, type = "response")
  list(
    metrics = rbind(
      score_predictions(y_test, mu_exposure, "Exposure-only NB", support),
      score_predictions(y_test, mu_full, "Covariate NB", support)
    ),
    coefficients = data.frame(
      support = support,
      term = names(coef(full)),
      estimate = unname(coef(full)),
      stringsAsFactors = FALSE
    ),
    predictions = data.frame(
      support = support,
      unit_id = seq_len(n),
      observed = y_test,
      predicted_exposure = mu_exposure,
      predicted_covariate = mu_full
    )
  )
}

segment <- fit_support("segment")
intersection <- fit_support("intersection")

write.csv(rbind(segment$metrics, intersection$metrics),
          file.path(root, "results/temporal_validation_metrics.csv"), row.names = FALSE)
write.csv(rbind(segment$coefficients, intersection$coefficients),
          file.path(root, "results/temporal_validation_coefficients.csv"), row.names = FALSE)
saveRDS(rbind(segment$predictions, intersection$predictions),
        file.path(root, "results/temporal_validation_predictions.rds"))

cat("Temporal validation complete.\n")
print(rbind(segment$metrics, intersection$metrics))

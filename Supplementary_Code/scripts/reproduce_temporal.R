suppressPackageStartupMessages(library(MASS))
set.seed(20260824)
a <- commandArgs(FALSE); script <- sub("^--file=", "", a[grep("^--file=",a)])
package_root <- normalizePath(file.path(dirname(script),".."))
args <- commandArgs(TRUE)
if(length(args)!=1) stop("Supply a NEW output directory")
root <- args[1]
if(file.exists(root)) stop("Refusing to overwrite existing output")
dir.create(file.path(root,"results"),recursive=TRUE)
inputs <- readRDS(file.path(package_root,"inputs/temporal_aggregate_inputs.rds"))
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
  z <- inputs[[support]]
  X <- as.data.frame(z$X)
  base_offset <- z$base_offset
  n <- nrow(X)
  names(X) <- paste0("x", seq_len(ncol(X)))
  y_train <- z$y_train
  y_test <- z$y_test
  d <- cbind(data.frame(y = y_train, off = base_offset + log(13)), X)
  full <- glm.nb(y ~ . - off + offset(off), data = d, control = glm.control(maxit = 100))
  exposure <- glm.nb(y ~ 1 + offset(off), data = d, control = glm.control(maxit = 100))
  nd <- cbind(data.frame(off = base_offset + log(5)), X)
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

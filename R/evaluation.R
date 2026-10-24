# Métricas, intervalos de confiança e testes de significância.

rmse <- function(y, p) sqrt(mean((p - y)^2))
mae <- function(y, p) mean(abs(p - y))

#' Concordância de Lin: correlação e calibração ao mesmo tempo.
ccc <- function(y, p) {
  2 * mean((y - mean(y)) * (p - mean(p))) /
    (mean((y - mean(y))^2) + mean((p - mean(p))^2) + (mean(y) - mean(p))^2)
}

#' Proporção de pares de misturas ordenados corretamente (empate vale meio).
c_index <- function(y, p) {
  pairs <- utils::combn(length(y), 2)
  dy <- sign(y[pairs[1, ]] - y[pairs[2, ]])
  dp <- sign(p[pairs[1, ]] - p[pairs[2, ]])
  ok <- dy != 0
  mean((dy[ok] == dp[ok]) + 0.5 * (dp[ok] == 0))
}

safe_cor <- function(y, p) if (stats::sd(p) > 1e-12) stats::cor(y, p) else NA_real_

metric_table <- function(pred, by = "model") {
  groups <- split(pred, pred[by], drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    out <- g[1, by, drop = FALSE]
    cbind(out, RMSE = rmse(g$true, g$pred), MAE = mae(g$true, g$pred),
          Pearson = safe_cor(g$true, g$pred), CCC = ccc(g$true, g$pred),
          C_index = c_index(g$true, g$pred))
  }))
  rownames(out) <- NULL
  out
}

#' IC 95% do RMSE por bootstrap percentil, reamostrando dentro de cada doador.
bootstrap_rmse <- function(pred, n_boot = CONFIG$n_bootstrap, seed = CONFIG$seed) {
  set.seed(seed)
  out <- do.call(rbind, lapply(split(pred, pred$model), function(g) {
    strata <- split(seq_len(nrow(g)), g$donor)
    vals <- replicate(n_boot, {
      idx <- unlist(lapply(strata, function(s) s[sample.int(length(s), replace = TRUE)]))
      rmse(g$true[idx], g$pred[idx])
    })
    data.frame(model = g$model[1], estimate = rmse(g$true, g$pred),
               ci_low = unname(stats::quantile(vals, 0.025)),
               ci_high = unname(stats::quantile(vals, 0.975)))
  }))
  rownames(out) <- NULL
  out
}

#' Wilcoxon pareado do erro absoluto (melhor modelo contra cada outro) + Holm.
#'
#' O p-valor por mistura vale para estes três doadores; a coluna de doadores
#' em que o melhor modelo venceu é o que fala de generalização.
paired_tests <- function(pred, reference) {
  pred$abs_err <- abs(pred$pred - pred$true)
  ref <- pred[pred$model == reference, ]
  others <- setdiff(unique(pred$model), reference)
  out <- do.call(rbind, lapply(others, function(m) {
    o <- pred[pred$model == m, ]
    stopifnot(identical(o$mixture_id, ref$mixture_id))
    w <- stats::wilcox.test(ref$abs_err, o$abs_err, paired = TRUE, alternative = "less")
    by_donor <- tapply(ref$abs_err - o$abs_err, ref$donor, mean)
    data.frame(comparison = paste(reference, "vs", m),
               median_diff = stats::median(ref$abs_err - o$abs_err),
               p_value = w$p.value, donors_better = sum(by_donor < 0),
               n_donors = length(by_donor))
  }))
  out$p_holm <- stats::p.adjust(out$p_value, method = "holm")
  out
}

#' AUC com IC de DeLong (pROC).
auc_ci <- function(score, positive) {
  roc <- pROC::roc(positive, score, levels = c(FALSE, TRUE), direction = "<", quiet = TRUE)
  ci <- as.numeric(suppressWarnings(pROC::ci.auc(roc, method = "delong")))
  c(auc = ci[2], ci_low = ci[1], ci_high = min(ci[3], 1))
}

#' AUC rápida pela estatística de Mann-Whitney.
auc_mw <- function(score, positive) {
  r <- rank(score)
  n1 <- sum(positive)
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * sum(!positive))
}

empirical_p <- function(observed, null) (1 + sum(null >= observed)) / (1 + length(null))

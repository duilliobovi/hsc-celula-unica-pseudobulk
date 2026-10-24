# Modelos que estimam a fração de HSC num pseudo-bulk.
# Interface comum: fit_<modelo>(treino, referência) devolve uma lista com
# $predict(contagens), que recebe contagens (genes x misturas).

#' Perfis de referência: contagens médias por célula de cada grupo.
#'
#' Usar contagens por célula (e não CPM) é o que permite estimar frações de
#' CÉLULAS: uma mistura é a soma das células, então mistura = S x n, em que n
#' é o número de células de cada tipo. Com CPM, o NNLS estimaria frações de
#' RNA, e tipos com mais RNA por célula pareceriam mais abundantes.
reference_profiles <- function(compartment, donors) {
  cells <- compartment$cells
  groups <- c(CONFIG$hsc_cluster, CONFIG$progenitor_clusters)
  per_donor <- lapply(donors, function(d) {
    sapply(groups, function(g) {
      sel <- cells$donor == d & cells$BioClassification == g
      Matrix::rowMeans(compartment$counts[, sel, drop = FALSE])
    })
  })
  Reduce(`+`, per_donor) / length(per_donor)
}

#' Genes mais específicos de cada grupo (log2 da razão contra o maior vizinho).
specific_genes <- function(profiles, n = CONFIG$n_specific_genes) {
  cpm <- t(t(profiles) / colSums(profiles)) * 1e6
  logp <- log2(cpm + 1)
  expressed <- apply(cpm, 1, max) >= CONFIG$min_cpm
  genes <- character(0)
  for (k in seq_len(ncol(logp))) {
    spec <- logp[, k] - apply(logp[, -k, drop = FALSE], 1, max)
    spec[!expressed] <- -Inf
    genes <- union(genes, names(sort(spec, decreasing = TRUE))[seq_len(n)])
  }
  genes
}

#' Escore de assinatura numa matriz log2(CPM + 1).
signature_score <- function(log_expr, signature) {
  up <- intersect(signature$up, rownames(log_expr))
  down <- intersect(signature$down, rownames(log_expr))
  s <- Matrix::colMeans(log_expr[up, , drop = FALSE])
  if (length(down) > 0) s <- s - Matrix::colMeans(log_expr[down, , drop = FALSE])
  s
}

fit_baseline <- function(train, ref) {
  m <- mean(train$fraction)
  list(name = "Referência (média)", predict = function(counts) rep(m, ncol(counts)))
}

fit_signature_lm <- function(train, ref) {
  x <- signature_score(log_cpm(train$counts), ref$signature)
  lm_fit <- stats::lm(y ~ x, data = data.frame(y = train$fraction, x = x))
  list(name = "Escore da assinatura", predict = function(counts) {
    x <- signature_score(log_cpm(counts), ref$signature)
    pmin(pmax(stats::predict(lm_fit, data.frame(x = x)), 0), 1)
  })
}

fit_nnls <- function(train, ref) {
  S <- ref$profiles[ref$genes, ]
  list(name = "NNLS (referência single-cell)", predict = function(counts) {
    apply(counts[ref$genes, , drop = FALSE], 2, function(m) {
      w <- nnls::nnls(S, m)$x
      if (sum(w) > 0) w[1] / sum(w) else 0
    })
  })
}

fit_glmnet <- function(train, ref) {
  x <- t(log_cpm(train$counts)[ref$genes, ])
  # Com dois doadores de treino não há como fazer validação interna agrupada
  # por doador (o glmnet exige ao menos 3 partições). A penalização é escolhida
  # por 5 partições aleatórias do treino; isso pode deixar o lambda um pouco
  # otimista, mas não contamina a avaliação, feita no doador de fora.
  set.seed(CONFIG$seed)
  model <- glmnet::cv.glmnet(x, train$fraction, family = "gaussian", alpha = 0.5,
                             nfolds = 5)
  list(name = "Elastic Net (glmnet)", predict = function(counts) {
    p <- stats::predict(model, newx = t(log_cpm(counts)[ref$genes, , drop = FALSE]),
                        s = "lambda.min")
    pmin(pmax(as.numeric(p), 0), 1)
  })
}

MODELS <- list(baseline = fit_baseline, signature = fit_signature_lm,
               nnls = fit_nnls, glmnet = fit_glmnet)

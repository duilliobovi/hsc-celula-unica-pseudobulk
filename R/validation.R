# Validação deixando um doador de fora (leave-one-donor-out), em dois níveis:
#   1. célula: a assinatura separa HSC de progenitores no doador novo?
#   2. pseudo-bulk: os modelos estimam a fração de HSC no doador novo?
# Em cada rodada, assinatura, perfis de referência, genes e modelos usam só
# os dois doadores de treino.

#' Escore por célula e AUC no doador de teste, com nula de genes aleatórios.
cell_level_lodo <- function(compartment, n_sets = CONFIG$n_random_sets,
                            seed = CONFIG$seed) {
  cells <- compartment$cells
  results <- list()
  scores <- list()
  nulls <- list()
  for (d in CONFIG$donors) {
    train <- setdiff(CONFIG$donors, d)
    sig <- consistent_signature(compartment, train)
    test <- cells$donor == d
    x <- log_normalize(compartment$counts[, test])
    score <- signature_score(x, sig)
    markers <- list(up = intersect(KNOWN_HSC_MARKERS, rownames(x)), down = character(0))
    marker_score <- signature_score(x, markers)
    hsc <- cells$is_hsc[test]
    a_sig <- auc_ci(score, hsc)
    a_mark <- auc_ci(marker_score, hsc)
    # Nula: conjuntos aleatórios do mesmo tamanho, sorteados entre genes
    # expressos no doador de teste (pelo menos 1% das células).
    expressed <- rownames(x)[Matrix::rowMeans(x > 0) >= 0.01]
    set.seed(seed)
    null <- vapply(seq_len(n_sets), function(i) {
      g <- sample(expressed, length(sig$up) + length(sig$down))
      s <- signature_score(x, list(up = g[seq_along(sig$up)],
                                   down = g[length(sig$up) + seq_along(sig$down)]))
      auc_mw(s, hsc)
    }, numeric(1))
    results[[d]] <- data.frame(
      donor = d, n_hsc = sum(hsc), n_progenitor = sum(!hsc),
      auc_signature = a_sig[["auc"]], ci_low = a_sig[["ci_low"]], ci_high = a_sig[["ci_high"]],
      auc_markers = a_mark[["auc"]], markers_ci_low = a_mark[["ci_low"]],
      markers_ci_high = a_mark[["ci_high"]],
      null_q95 = unname(stats::quantile(null, 0.95)),
      p_random = empirical_p(a_sig[["auc"]], null))
    scores[[d]] <- data.frame(donor = d, cluster = cells$BioClassification[test],
                              is_hsc = hsc, score = score)
    nulls[[d]] <- data.frame(donor = d, auc = null)
  }
  summary <- do.call(rbind, results)
  rownames(summary) <- NULL
  list(summary = summary, scores = do.call(rbind, scores), null = do.call(rbind, nulls))
}

#' Rodadas de pseudo-bulk: treina em dois doadores e testa no terceiro.
pseudobulk_lodo <- function(compartment, seed = CONFIG$seed) {
  preds <- list()
  for (i in seq_along(CONFIG$donors)) {
    d <- CONFIG$donors[i]
    train_donors <- setdiff(CONFIG$donors, d)
    train <- simulate_pseudobulks(compartment, train_donors,
                                  CONFIG$n_train_per_donor * length(train_donors),
                                  seed + i)
    test <- simulate_pseudobulks(compartment, d, CONFIG$n_test, seed + 100 + i)
    profiles <- reference_profiles(compartment, train_donors)
    signature <- consistent_signature(compartment, train_donors)
    ref <- list(profiles = profiles, signature = signature,
                genes = union(specific_genes(profiles),
                              intersect(c(signature$up, signature$down), rownames(profiles))))
    for (m in names(MODELS)) {
      model <- MODELS[[m]](train, ref)
      preds[[length(preds) + 1]] <- data.frame(
        donor = d, mixture_id = paste(d, seq_len(CONFIG$n_test)), model = m,
        true = test$fraction, pred = model$predict(test$counts))
    }
  }
  do.call(rbind, preds)
}

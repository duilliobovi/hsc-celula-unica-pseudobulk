# Ponte single-cell -> bulk: a assinatura aprendida nas células funciona num
# RNA-seq bulk de outro laboratório (GSE74246, Corces 2016)?

#' log2(CPM + 1) das populações normais purificadas do GSE74246.
read_corces <- function(path) {
  counts <- data.table::fread(path, data.table = FALSE)
  genes <- counts[[1]]
  counts <- rowsum(as.matrix(counts[, -1]), genes)
  donor <- sub("-.*$", "", colnames(counts))
  normal <- !startsWith(donor, "SU")
  counts <- counts[, normal]
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  list(log_cpm = log2(cpm[rowSums(cpm > 1) >= 3, ] + 1),
       samples = data.frame(sample = colnames(counts), donor = donor[normal],
                            cell_type = sub("^[^-]+-", "", colnames(counts))))
}

#' Escore em postos (imune à escala da plataforma), AUC e nula aleatória.
external_validation <- function(corces, signature, n_sets = CONFIG$n_random_sets,
                                seed = CONFIG$seed) {
  x <- corces$log_cpm
  ranks <- apply(x, 2, function(v) (rank(v) - 0.5) / length(v))
  rownames(ranks) <- rownames(x)
  up <- intersect(signature$up, rownames(ranks))
  down <- intersect(signature$down, rownames(ranks))
  score <- colMeans(ranks[up, ]) - colMeans(ranks[down, ])
  hsc <- corces$samples$cell_type == "HSC"
  set.seed(seed)
  null <- vapply(seq_len(n_sets), function(i) {
    g <- sample(rownames(ranks), length(up) + length(down))
    s <- colMeans(ranks[g[seq_along(up)], ]) - colMeans(ranks[g[length(up) + seq_along(down)], ])
    auc_mw(s, hsc)
  }, numeric(1))
  auc <- auc_mw(score, hsc)
  mpp <- corces$samples$cell_type == "MPP"
  list(scores = data.frame(corces$samples, score = score),
       genes_used = c(up = length(up), down = length(down)),
       auc = auc, p_random = empirical_p(auc, null),
       hsc_vs_mpp_p = stats::wilcox.test(score[hsc], score[mpp], alternative = "greater",
                                         exact = TRUE)$p.value,
       hsc_rank = rank(-score)[hsc])
}

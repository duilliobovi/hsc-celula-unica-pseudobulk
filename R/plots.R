# Gráficos do estudo (ggplot2).

MODEL_LABELS <- c(baseline = "Referência (média)", signature = "Escore da assinatura",
                  nnls = "NNLS (referência single-cell)", glmnet = "Elastic Net (glmnet)")

theme_study <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

save_plot <- function(plot, name, width = 7, height = 4.5) {
  dir.create(CONFIG$figures_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(CONFIG$figures_dir, name), plot, width = width,
                  height = height, dpi = 130, bg = "white")
  invisible(plot)
}

plot_umap <- function(cells) {
  cells$group <- ifelse(cells$BioClassification == CONFIG$hsc_cluster, "HSC",
                        ifelse(cells$BioClassification %in% CONFIG$progenitor_clusters,
                               "Progenitores CD34+", "Outras células"))
  cells <- cells[order(cells$group == "HSC"), ]
  ggplot2::ggplot(cells, ggplot2::aes(UMAP1, UMAP2, colour = group)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.6) +
    ggplot2::scale_colour_manual(values = c(HSC = "#d95f02",
                                            "Progenitores CD34+" = "#1b9e77",
                                            "Outras células" = "grey80"), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::labs(title = "35.582 células de medula, sangue e CD34+ (UMAP dos autores)") +
    theme_study()
}

plot_composition <- function(cells) {
  tab <- as.data.frame(table(donor = cells$donor, cluster = cells$BioClassification))
  ggplot2::ggplot(tab, ggplot2::aes(donor, Freq, fill = cluster)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(title = "Composição do compartimento CD34+ por doador",
                  x = NULL, y = "Células", fill = "Grupo") +
    theme_study()
}

plot_pseudoreplication <- function(pb, cell) {
  df <- merge(pb[, c("gene", "p_value")], cell[, c("gene", "p_value")], by = "gene",
              suffixes = c("_pb", "_cell"))
  ggplot2::ggplot(df, ggplot2::aes(-log10(p_value_pb), -log10(pmax(p_value_cell, 1e-300)))) +
    ggplot2::geom_point(size = 0.4, alpha = 0.3) +
    ggplot2::geom_abline(linetype = 2, colour = "#b2182b") +
    ggplot2::labs(title = "Pseudorreplicação: célula como réplica infla a significância",
                  subtitle = "Cada ponto é um gene; acima da linha, o teste por célula é mais 'confiante'",
                  x = "-log10(p) pseudo-bulk (unidade = doador)",
                  y = "-log10(p) por célula (unidade = célula)") +
    theme_study()
}

plot_volcano <- function(table, markers = KNOWN_HSC_MARKERS) {
  table$class <- ifelse(table$fdr < CONFIG$fdr & abs(table$logFC) > CONFIG$min_log2_fc,
                        ifelse(table$logFC > 0, "maior em HSC", "menor em HSC"), "n.s.")
  lab <- table[table$gene %in% markers, ]
  ggplot2::ggplot(table, ggplot2::aes(logFC, -log10(p_value), colour = class)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), colour = "black",
                             size = 3, max.overlaps = 30, min.segment.length = 0) +
    ggplot2::scale_colour_manual(values = c("maior em HSC" = "#d95f02",
                                            "menor em HSC" = "#1b9e77", "n.s." = "grey75"),
                                 name = NULL) +
    ggplot2::labs(title = "HSC contra progenitores CD34+ (pseudo-bulk pareado por doador)",
                  subtitle = "Rótulos: marcadores clássicos de HSC", x = "log2 fold change",
                  y = "-log10(p)") +
    theme_study()
}

plot_cell_scores <- function(scores, summary) {
  scores$group <- ifelse(scores$is_hsc, "HSC", "Progenitores")
  labels <- stats::setNames(sprintf("%s (AUC %.3f)", summary$donor, summary$auc_signature),
                            summary$donor)
  ggplot2::ggplot(scores, ggplot2::aes(score, fill = group)) +
    ggplot2::geom_density(alpha = 0.55, colour = NA) +
    ggplot2::facet_wrap(~donor, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_fill_manual(values = c(HSC = "#d95f02", Progenitores = "#1b9e77"),
                               name = NULL) +
    ggplot2::labs(title = "Escore por célula no doador deixado de fora",
                  x = "Escore da assinatura", y = "Densidade") +
    theme_study() + ggplot2::theme(legend.position = "bottom")
}

plot_random_null <- function(null, summary) {
  ggplot2::ggplot(null, ggplot2::aes(auc)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey80", colour = "grey50") +
    ggplot2::geom_vline(data = summary, ggplot2::aes(xintercept = auc_signature),
                        colour = "#b2182b", linewidth = 1) +
    ggplot2::facet_wrap(~donor) +
    ggplot2::labs(title = "Assinatura contra 1.000 conjuntos aleatórios de genes",
                  x = "AUC por célula", y = "Conjuntos") +
    theme_study()
}

plot_predicted <- function(pred) {
  pred <- pred[pred$model != "baseline", ]
  pred$model <- factor(MODEL_LABELS[pred$model], levels = MODEL_LABELS[-1])
  ggplot2::ggplot(pred, ggplot2::aes(true, pred, colour = donor)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.6) +
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::facet_wrap(~model) +
    ggplot2::coord_equal(xlim = c(0, 0.35), ylim = c(0, 0.6)) +
    ggplot2::labs(title = "Fração de HSC em pseudo-bulks do doador deixado de fora",
                  x = "Fração real", y = "Fração prevista", colour = "Doador") +
    theme_study()
}

plot_rmse_ci <- function(ci) {
  ci$label <- factor(MODEL_LABELS[ci$model], levels = MODEL_LABELS[ci$model[order(-ci$estimate)]])
  ggplot2::ggplot(ci, ggplot2::aes(estimate, label)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0.2,
                           orientation = "y") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(title = "RMSE com IC 95% (bootstrap dentro de cada doador)",
                  x = "RMSE da fração de HSC", y = NULL) +
    theme_study()
}

plot_external <- function(scores) {
  scores$cell_type <- stats::reorder(scores$cell_type, -scores$score, FUN = stats::median)
  ggplot2::ggplot(scores, ggplot2::aes(cell_type, score)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92") +
    ggplot2::geom_jitter(ggplot2::aes(colour = cell_type == "HSC"), width = 0.15, size = 2,
                         show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey40")) +
    ggplot2::labs(title = "Assinatura single-cell aplicada ao RNA-seq bulk (GSE74246)",
                  x = NULL, y = "Escore da assinatura (postos)") +
    theme_study() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

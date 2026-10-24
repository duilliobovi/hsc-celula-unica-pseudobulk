# Executa o estudo completo e grava tabelas, figuras e métricas em reports/.
#
# Uso (na raiz do repositório):
#   Rscript main.R
# Tempo: cerca de 5 minutos; memória: cerca de 4 GB (o objeto tem 35.582 células).

suppressPackageStartupMessages({
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
})

run_study <- function() {
  set.seed(CONFIG$seed)
  dir.create(CONFIG$tables_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Dados
  data <- load_granja(download_if_missing(CONFIG$granja_url))
  all_cells <- data$cells
  compartment <- cd34_compartment(data)
  rm(data)
  invisible(gc())
  message(sprintf("Células no compartimento CD34+: %d (HSC: %d)",
                  nrow(compartment$cells), sum(compartment$cells$is_hsc)))

  # 2. Expressão diferencial: pseudo-bulk (correto) e por célula (comparação)
  de <- pseudobulk_de(compartment)
  signature <- build_signature(de$table)
  cell_de <- cell_level_welch(compartment, de$table$gene)
  markers_in <- intersect(KNOWN_HSC_MARKERS, de$table$gene)
  markers_hit <- intersect(markers_in, signature$up)
  marker_p <- stats::phyper(length(markers_hit) - 1, length(markers_in),
                            nrow(de$table) - length(markers_in), length(signature$up),
                            lower.tail = FALSE)

  # 3. Validação deixando um doador de fora
  cells_lodo <- cell_level_lodo(compartment)
  pred <- pseudobulk_lodo(compartment)
  overall <- metric_table(pred)
  per_donor <- metric_table(pred, by = c("model", "donor"))
  ci <- bootstrap_rmse(pred)
  best <- ci$model[ci$model != "baseline"][which.min(ci$estimate[ci$model != "baseline"])]
  tests <- paired_tests(pred, best)

  # 4. Ponte para o RNA-seq bulk
  corces <- read_corces(download_if_missing(CONFIG$corces_counts_url))
  external <- external_validation(corces, signature)

  # 5. Tabelas
  write <- function(x, name) utils::write.csv(x, file.path(CONFIG$tables_dir, name),
                                              row.names = FALSE)
  write(de$table[order(de$table$p_value), ], "pseudobulk_de.csv")
  write(data.frame(direction = rep(c("up", "down"), c(length(signature$up),
                                                      length(signature$down))),
                   gene = c(signature$up, signature$down)), "hsc_signature.csv")
  write(cells_lodo$summary, "cell_level_auc.csv")
  write(pred, "pseudobulk_predictions.csv")
  write(per_donor, "metrics_per_donor.csv")
  write(tests, "paired_tests.csv")
  write(external$scores, "external_scores.csv")

  # 6. Figuras
  save_plot(plot_umap(all_cells), "01_umap.png", width = 7.5, height = 5.5)
  save_plot(plot_composition(compartment$cells), "02_composition.png")
  save_plot(plot_pseudoreplication(de$table, cell_de), "03_pseudoreplication.png")
  save_plot(plot_volcano(de$table), "04_volcano.png")
  save_plot(plot_cell_scores(cells_lodo$scores, cells_lodo$summary), "05_cell_scores.png",
            width = 9, height = 3.8)
  save_plot(plot_random_null(cells_lodo$null, cells_lodo$summary), "06_random_null.png",
            width = 9, height = 3.5)
  save_plot(plot_predicted(pred), "07_predicted_vs_true.png", width = 9, height = 4)
  save_plot(plot_rmse_ci(ci), "08_rmse_ci.png", width = 6.5, height = 3)
  save_plot(plot_external(external$scores), "09_external_validation.png")

  # 7. Métricas
  metrics <- list(
    n_cells_compartment = nrow(compartment$cells),
    cells_by_donor = as.data.frame.matrix(table(compartment$cells$donor,
                                                ifelse(compartment$cells$is_hsc, "HSC",
                                                       "progenitor"))),
    genes_tested = nrow(de$table),
    signature = list(n_up = signature$n_up, n_down = signature$n_down,
                     top_up = head(signature$up, 20), top_down = head(signature$down, 20)),
    known_markers = list(measured = markers_in, in_top50_up = markers_hit,
                         p_value = signif(marker_p, 3)),
    pseudoreplication = list(
      significant_pseudobulk = sum(de$table$fdr < CONFIG$fdr),
      significant_cell_level = sum(cell_de$fdr < CONFIG$fdr),
      median_neglog10p_pseudobulk = stats::median(-log10(de$table$p_value)),
      median_neglog10p_cell_level = stats::median(-log10(pmax(cell_de$p_value, 1e-300)))),
    cell_level_lodo = cells_lodo$summary,
    pseudobulk_overall = overall,
    pseudobulk_per_donor = per_donor,
    rmse_ci = ci, best_model = best, paired_tests = tests,
    external = list(genes_used = as.list(external$genes_used), auc = external$auc,
                    p_random = external$p_random, hsc_vs_mpp_p = external$hsc_vs_mpp_p,
                    hsc_ranks = as.list(external$hsc_rank),
                    median_score = as.list(tapply(external$scores$score,
                                                  external$scores$cell_type, stats::median))),
    versions = list(R = R.version.string,
                    limma = as.character(utils::packageVersion("limma")),
                    glmnet = as.character(utils::packageVersion("glmnet")),
                    nnls = as.character(utils::packageVersion("nnls")))
  )
  jsonlite::write_json(metrics, file.path(CONFIG$reports_dir, "metrics.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 6)
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(CONFIG$reports_dir, "session_info.txt"))
  message("Métricas gravadas em reports/metrics.json")
  invisible(metrics)
}

if (sys.nframe() == 0L) run_study()

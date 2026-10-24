# Expressão diferencial: pseudo-bulk por doador (correto) e por célula
# (pseudorreplicação, mostrado só para comparar).

#' limma-voom pareado por doador: HSC contra progenitores CD34+.
#'
#' Cada doador contribui com dois pseudo-bulks (HSC e progenitores). O modelo
#' ~ doador + tipo compara HSC e progenitores dentro do mesmo doador, o que
#' remove as diferenças de doador, tecido (medula total ou CD34+) e
#' profundidade de sequenciamento. Com três doadores há só dois graus de
#' liberdade residuais; o encolhimento de variâncias do eBayes é o que torna o
#' teste viável.
pseudobulk_de <- function(compartment, min_cpm = CONFIG$min_cpm) {
  cells <- compartment$cells
  group <- paste(cells$donor, ifelse(cells$is_hsc, "HSC", "Prog"), sep = "_")
  pb <- pseudobulk(compartment$counts, group)
  info <- data.frame(id = colnames(pb), donor = sub("_.*$", "", colnames(pb)),
                     type = factor(sub("^.*_", "", colnames(pb)), levels = c("Prog", "HSC")))
  cpm <- t(t(pb) / colSums(pb)) * 1e6
  keep <- rowSums(cpm >= min_cpm) >= 2
  design <- stats::model.matrix(~ donor + type, data = info)
  v <- limma::voom(pb[keep, ], design, lib.size = colSums(pb))
  fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
  tab <- limma::topTable(fit, coef = "typeHSC", number = Inf, sort.by = "none")
  list(table = data.frame(gene = rownames(tab), logFC = tab$logFC, ave_expr = tab$AveExpr,
                          t = tab$t, p_value = tab$P.Value, fdr = tab$adj.P.Val,
                          stringsAsFactors = FALSE),
       pseudobulk = pb, info = info)
}

#' Teste t de Welch célula a célula (pseudorreplicação, só para comparação).
#'
#' Trata cada célula como réplica independente. Como células do mesmo doador
#' são correlacionadas, o número efetivo de réplicas é muito menor que o de
#' células, e os p-valores ficam artificialmente pequenos (Squair et al.,
#' 2021, Nature Communications 12:5692).
cell_level_welch <- function(compartment, genes) {
  x <- log_normalize(compartment$counts[genes, ])
  hsc <- compartment$cells$is_hsc
  stats_for <- function(m) {
    n <- ncol(m)
    mu <- Matrix::rowMeans(m)
    v <- (Matrix::rowSums(m^2) - n * mu^2) / (n - 1)
    list(n = n, mu = mu, v = v)
  }
  a <- stats_for(x[, hsc])
  b <- stats_for(x[, !hsc])
  se2 <- a$v / a$n + b$v / b$n
  t <- (a$mu - b$mu) / sqrt(se2)
  df <- se2^2 / ((a$v / a$n)^2 / (a$n - 1) + (b$v / b$n)^2 / (b$n - 1))
  p <- 2 * stats::pt(-abs(t), df)
  data.frame(gene = genes, diff = a$mu - b$mu, t = t, p_value = p,
             fdr = stats::p.adjust(p, method = "BH"), stringsAsFactors = FALSE)
}

#' Assinatura: genes com FDR e tamanho de efeito acima dos limiares.
build_signature <- function(table, fdr = CONFIG$fdr, min_lfc = CONFIG$min_log2_fc,
                            size = CONFIG$signature_size) {
  up <- table[table$fdr < fdr & table$logFC > min_lfc, ]
  down <- table[table$fdr < fdr & table$logFC < -min_lfc, ]
  list(up = head(up$gene[order(-up$t)], size), down = head(down$gene[order(down$t)], size),
       n_up = nrow(up), n_down = nrow(down))
}

#' Assinatura rápida para as rodadas de validação, com dois doadores de treino.
#'
#' Com só dois doadores, o limma teria um grau de liberdade; em vez disso,
#' exige que o gene mude na mesma direção nos dois doadores e ordena pelo
#' menor log2 fold change entre eles (o gene precisa ser forte em ambos).
consistent_signature <- function(compartment, train_donors, size = CONFIG$signature_size,
                                 min_cpm = CONFIG$min_cpm) {
  cells <- compartment$cells
  lfc <- sapply(train_donors, function(d) {
    sel <- cells$donor == d
    pb <- pseudobulk(compartment$counts[, sel], ifelse(cells$is_hsc[sel], "HSC", "Prog"))
    cpm <- t(t(pb) / colSums(pb)) * 1e6
    out <- log2(cpm[, "HSC"] + 1) - log2(cpm[, "Prog"] + 1)
    out[pmax(cpm[, "HSC"], cpm[, "Prog"]) < min_cpm] <- NA
    out
  })
  up_strength <- apply(lfc, 1, min)
  down_strength <- apply(lfc, 1, max)
  up <- names(sort(up_strength[!is.na(up_strength) & up_strength > 0], decreasing = TRUE))
  down <- names(sort(down_strength[!is.na(down_strength) & down_strength < 0]))
  list(up = head(up, size), down = head(down, size))
}

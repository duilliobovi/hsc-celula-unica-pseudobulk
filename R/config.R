# Configuração central: caminhos, URLs, sementes e decisões de análise.

CONFIG <- list(
  seed = 2026L,
  data_dir = "data",
  reports_dir = "reports",
  figures_dir = file.path("reports", "figures"),
  tables_dir = file.path("reports", "tables"),

  # Granja et al. (2019), Nature Biotechnology 37:1458 - objeto já anotado,
  # publicado pelos autores (link no README do repositório GreenleafLab).
  granja_url = paste0("https://jeffgranja.s3.amazonaws.com/MPAL-10x/",
                      "Supplementary_Data/Healthy-Data/",
                      "scRNA-Healthy-Hematopoiesis-191120.rds"),
  # Corces et al. (2016): RNA-seq bulk de populações purificadas (validação externa).
  corces_counts_url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/",
                             "GSE74246/suppl/GSE74246_RNAseq_All_Counts.txt.gz"),

  # Grupos celulares do compartimento CD34+ (anotação dos autores).
  hsc_cluster = "01_HSC",
  progenitor_clusters = c("02_Early.Eryth", "04_Early.Baso", "05_CMP.LMPP",
                          "06_CLP.1", "07_GMP", "15_CLP.2"),
  # Doadores com HSC: D1 (medula óssea total, duas corridas técnicas) e
  # D2, D3 (células CD34+ purificadas). D4 (sangue periférico) quase não tem HSC.
  donors = c("D1", "D2", "D3"),

  # Expressão diferencial em pseudo-bulk.
  min_cpm = 10,             # CPM mínimo no pseudo-bulk para testar o gene
  fdr = 0.05,
  min_log2_fc = 1,
  signature_size = 50L,

  # Misturas pseudo-bulk com células reais.
  cells_per_mixture = 500L,
  hsc_fraction_range = c(0, 0.3),
  n_train_per_donor = 300L,
  n_test = 200L,
  n_specific_genes = 40L,   # genes específicos por grupo na matriz de referência

  # Inferência.
  n_random_sets = 1000L,
  n_bootstrap = 2000L
)

# Marcadores clássicos de HSC humanas, fixados antes da análise.
KNOWN_HSC_MARKERS <- c("HLF", "CRHBP", "AVP", "MECOM", "PROM1", "HOPX",
                       "MLLT3", "PRDM16", "SPINK2", "MEIS1", "CD34", "HES1",
                       "PROCR", "THY1", "EMCN", "GATA2", "HOXA9")

# Baixa o objeto dos autores (≈ 270 MB) e o RNA-seq bulk do GEO para data/.
source("R/config.R"); source("R/data.R")
for (url in c(CONFIG$granja_url, CONFIG$corces_counts_url)) {
  message("Baixando ", basename(url)); download_if_missing(url)
}

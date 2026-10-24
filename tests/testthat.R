# Executa todos os testes: Rscript tests/testthat.R (na raiz do repositório).
library(testthat)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)

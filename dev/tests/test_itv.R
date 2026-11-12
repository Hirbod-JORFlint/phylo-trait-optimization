# Run from the repository root: Rscript dev/tests/test_itv.R
dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
source("dev/build/interactive_tree_viz_new.R")
library(ape)
set.seed(42)
n <- 20
tr <- rtree(n, tip.label = paste0("Tip", 1:n))
Td <- matrix(runif(n * 6), nrow = n, dimnames = list(tr$tip.label, paste0("Trait", 1:6)))
res <- interactive_tree_viz(tr, Td, reorder_olo = TRUE,
  output_file = "dev/tmp/test_default.html",
  save_newick = "dev/tmp/test_default.tre",
  save_heatmap_csv = "dev/tmp/test_default.csv",
  save_config = "dev/tmp/test_default.json")
cat("file:", file.exists("dev/tmp/test_default.html"), "\n")
cat("order ok:", setequal(res$order, tr$tip.label), "\n")
cat("config:", res$config$layout, res$config$theme, "\n")

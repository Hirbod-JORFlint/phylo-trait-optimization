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
  config = tree_viz_config(
    heatmap_view_width = 1100, heatmap_view_height = 720,
    trait_view_width = 900, trait_view_height = 620,
    auto_fit_subviews = TRUE, trait_subview_alignment = "center"
  ),
  save_newick = "dev/tmp/test_default.tre",
  save_heatmap_csv = "dev/tmp/test_default.csv",
  save_config = "dev/tmp/test_default.json")
cat("file:", file.exists("dev/tmp/test_default.html"), "\n")
cat("order ok:", setequal(res$order, tr$tip.label), "\n")
cat("config:", res$config$layout, res$config$theme, "\n")
stopifnot(res$config$heatmap_view_width == 1100,
          res$config$heatmap_view_height == 720,
          res$config$trait_view_width == 900,
          res$config$trait_view_height == 620,
          isTRUE(res$config$auto_fit_subviews),
          res$config$trait_subview_alignment == "center")
stopifnot(res$widget$width == 1100, res$widget$height == 720)

# Extract the interactive-tree-viz-function chunk from the Rmd and smoke-test it.
# Run from the repository root: Rscript dev/tests/test_embedded.R
dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
rmd <- readLines("diagnostics_on_tree_reorder.Rmd", warn = FALSE, encoding = "UTF-8")
start <- grep("^```\\{r interactive-tree-viz-function", rmd)
stopifnot(length(start) == 1L)
ends <- which(seq_along(rmd) > start & grepl("^```\\s*$", rmd))
end <- ends[1]
code <- rmd[(start + 1L):(end - 1L)]
writeLines(code, "dev/tmp/extracted_function.R", useBytes = TRUE)

source("dev/tmp/extracted_function.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(ape))
set.seed(42)
tree <- rtree(24, tip.label = paste0("Tip", 1:24))
traits <- matrix(runif(24 * 6), 24, 6,
                 dimnames = list(tree$tip.label, paste0("Trait", 1:6)))

res <- interactive_tree_viz(
  tree = tree, traits = traits, reorder_olo = TRUE, scale_traits = TRUE,
  output_file = "dev/tmp/embedded_test.html",
  config = tree_viz_config(
    layout = "circular", title = "Embedded raw-string smoke test",
    theme = "dark", height = 700, aspect_ratio = 1.6,
    highlight = list(trait = "Trait1", max = 0.5),
    heatmap_columns = colnames(traits)[1:4]
  ),
  save_newick = "dev/tmp/embedded_test.tre",
  save_heatmap_csv = "dev/tmp/embedded_test.csv",
  save_config = "dev/tmp/embedded_test.json"
)
stopifnot(file.exists(res$file), file.exists("dev/tmp/embedded_test.tre"),
          file.exists("dev/tmp/embedded_test.csv"),
          file.exists("dev/tmp/embedded_test.json"))
stopifnot(identical(res$config$layout, "circular"),
          identical(res$config$theme, "dark"))
stopifnot(isTRUE(all.equal(sort(res$order), sort(tree$tip.label))))

html <- paste(readLines(res$file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(!grepl("__PAYLOAD__", html, fixed = TRUE))
stopifnot(grepl("var TD =", html, fixed = TRUE),
          grepl("circular", html, fixed = TRUE))

exports <- export_tree_viz_data(res, dir = "dev/tmp/embedded_exports", prefix = "smoke")
stopifnot(all(file.exists(exports)))
cfg2 <- load_tree_viz_config(exports[["config_json"]])
stopifnot(identical(cfg2$layout, "circular"))

cat("Embedded function chunk OK; widget:", res$file, "\n")
cat("HTML size:", file.info(res$file)$size, "bytes\n")

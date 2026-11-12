# Run from the repository root: Rscript dev/tests/test_helpers.R
dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
source("dev/build/interactive_tree_viz_new.R")
library(ape)
set.seed(11)
tr <- rtree(12, tip.label = paste0("T", 1:12))
Td <- matrix(runif(12 * 4), 12, 4,
             dimnames = list(tr$tip.label, paste0("Q", 1:4)))

cfg <- tree_viz_config(layout = "fan", title = "Helper test",
                       tip_filter = list(trait = "Q1", min = 0.2, max = 0.8),
                       highlight = list(trait = "Q2", min = 0.5, color = "#ff0000"),
                       heatmap_columns = c("Q1", "Q3"),
                       clade = 13:14,
                       heatmap_view_width = 1100, heatmap_view_height = 720,
                       trait_view_width = 900, trait_view_height = 620)
stopifnot(cfg$layout == "fan", cfg$title == "Helper test")

p <- tree_viz_preset("dark", title = "Dark preset")
stopifnot(p$theme == "dark", p$title == "Dark preset", p$layout == "rectangular")

save_tree_viz_config(p, "dev/tmp/helper_cfg.json")
p2 <- load_tree_viz_config("dev/tmp/helper_cfg.json")
stopifnot(p2$theme == "dark", p2$title == "Dark preset")

p3 <- tree_viz_config(heatmap_view_width = 1200, trait_view_height = 600)
save_tree_viz_config(p3, "dev/tmp/helper_cfg2.json")
p4 <- load_tree_viz_config("dev/tmp/helper_cfg2.json")
stopifnot(p4$heatmap_view_width == 1200, p4$trait_view_height == 600)

res <- interactive_tree_viz(
  tr, Td, reorder_olo = FALSE, output_file = "dev/tmp/helper_widget.html",
  config = cfg, save_newick = "dev/tmp/helper_tree.tre",
  save_heatmap_csv = "dev/tmp/helper_heatmap.csv",
  save_config = "dev/tmp/helper_saved_cfg.json"
)
stopifnot(file.exists("dev/tmp/helper_widget.html"),
          file.exists("dev/tmp/helper_tree.tre"),
          file.exists("dev/tmp/helper_heatmap.csv"),
          file.exists("dev/tmp/helper_saved_cfg.json"))
stopifnot(res$config$layout == "fan", res$config$theme == "light")
stopifnot(identical(res$config$tip_filter$trait, "Q1"))
stopifnot(identical(res$config$heatmap_columns, c("Q1", "Q3")))
stopifnot(res$config$heatmap_view_width == 1100,
          res$config$heatmap_view_height == 720,
          res$config$trait_view_width == 900,
          res$config$trait_view_height == 620)

out <- export_tree_viz_data(res, dir = "dev/tmp/helper_out", prefix = "demo")
stopifnot(all(file.exists(out)))
hm <- read.csv(out[["heatmap_csv"]], check.names = FALSE)
stopifnot(ncol(hm) == 5L, nrow(hm) == 12L)
cfg_back <- jsonlite::read_json(out[["config_json"]], simplifyVector = TRUE)
stopifnot(cfg_back$layout == "fan")

# A config list may also be supplied through `preset`.
res2 <- interactive_tree_viz(
  tr, Td, reorder_olo = FALSE, output_file = "dev/tmp/helper_widget2.html",
  preset = tree_viz_preset("publication")
)
stopifnot(res2$config$heatmap_colorscale == "Cividis",
          res2$config$responsive == FALSE)

cat("helper tests OK\n")

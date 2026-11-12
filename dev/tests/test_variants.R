dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
source("dev/build/interactive_tree_viz_new.R")
library(ape)
set.seed(7)
n <- 30
tr <- rtree(n, tip.label = paste0("Sp", sprintf("%02d", 1:n)))
Td <- matrix(runif(n * 6), nrow = n,
             dimnames = list(tr$tip.label, paste0("Trait", 1:6)))
rownames(Td) <- tr$tip.label

base_args <- list(
  tree = tr, traits = Td, reorder_olo = FALSE,
  output_file = "dev/tmp/test_slanted.html",
  config = list(layout = "slanted", responsive = FALSE, width = 1200, height = 800)
)
do.call(interactive_tree_viz, base_args)

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_circular.html",
  config = list(layout = "circular", responsive = FALSE, width = 1000, height = 900)
)))

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_fan.html",
  config = list(layout = "fan", fan_angle = 240, responsive = FALSE,
                width = 1000, height = 900)
)))

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_dark.html",
  config = list(theme = "dark", responsive = FALSE, width = 1200, height = 800,
                heatmap_colorscale = "Plasma")
)))

clade_node <- n + 2
do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_filter.html",
  config = list(
    responsive = FALSE, width = 1200, height = 800,
    tip_filter = list(trait = "Trait1", min = 0.25, max = 0.75),
    highlight = list(trait = "Trait2", min = 0.6, color = "#e41a1c"),
    heatmap_columns = c("Trait1", "Trait3", "Trait5"),
    clade = clade_node,
    subtitle = "filtered + highlighted"
  )
)))

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_ids.html",
  config = list(
    responsive = FALSE, width = 1200, height = 800,
    title = "Custom title", subtitle = "custom subtitle",
    tip_labels = FALSE, show_internal_node_ids = TRUE,
    tip_node_color = "#d62728", internal_node_color = "#2ca02c",
    branch_color = "#333333", branch_width = 2.5,
    heatmap_colorscale = "Cividis", heatmap_range = c(0, 1),
    heatmap_title = "Traits", heatmap_colorbar_title = "scaled value",
    heatmap_panel_title = "Heatmap panel", tree_panel_title = "Tree panel"
  )
)))

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_preset.html", preset = "publication"
)))

do.call(interactive_tree_viz, modifyList(base_args, list(
  output_file = "dev/tmp/test_aspect.html",
  config = list(responsive = TRUE, aspect_ratio = 1.5)
)))

cat("variants done\n")

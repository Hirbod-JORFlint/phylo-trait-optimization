# Extract the demo chunk and run it against synthetic upstream objects.
# Run from the repository root: Rscript dev/tests/test_demo_chunk.R
dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
rmd <- readLines("diagnostics_on_tree_reorder.Rmd", warn = FALSE, encoding = "UTF-8")
start <- grep("^```\\{r interactive-tree-viz-demo", rmd)
stopifnot(length(start) == 1L)
ends <- which(seq_along(rmd) > start & grepl("^```\\s*$", rmd))
end <- ends[1]
demo <- rmd[(start + 1L):(end - 1L)]
writeLines(demo, "dev/tmp/extracted_demo.R", useBytes = TRUE)

source("dev/tmp/extracted_function.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(ape))
set.seed(7)
my_tree <- rtree(20, tip.label = paste0("Tip", 1:20))
traits <- matrix(runif(20 * 6), 20, 6,
                 dimnames = list(my_tree$tip.label, paste0("Trait", 1:6)))

# Run the demo inside dev/tmp so its generated files do not clutter the root.
old_wd <- setwd("dev/tmp")
on.exit(setwd(old_wd), add = TRUE)
invisible(suppressWarnings(suppressMessages(
  eval(parse(text = demo, keep.source = TRUE), envir = globalenv())
)))

stopifnot(exists("res_viz"), exists("res_dark"), exists("exported"))
stopifnot(file.exists("interactive_tree_heatmap.html"),
          file.exists("interactive_tree_heatmap_dark.html"))
stopifnot(length(exported) == 3L)
cat("Demo chunk OK\n")

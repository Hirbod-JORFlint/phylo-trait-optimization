# Build a minimal Rmd that contains the real function and demo chunks plus
# synthetic upstream objects, then render it.
# Run from the repository root: Rscript dev/tests/make_min_rmd.R
dir.create("dev/tmp", showWarnings = FALSE, recursive = TRUE)
rmd <- readLines("diagnostics_on_tree_reorder.Rmd", warn = FALSE, encoding = "UTF-8")

chunk_lines <- function(label) {
  start <- grep(paste0("^```\\{r ", label), rmd)
  stopifnot(length(start) == 1L)
  ends <- which(seq_along(rmd) > start & grepl("^```\\s*$", rmd))
  end <- ends[1]
  rmd[start:end]
}

min_rmd <- c(
  "---",
  'title: "interactive_tree_viz minimal render test"',
  "output: html_document",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = TRUE)",
  "```",
  "",
  "```{r synthetic-data}",
  "suppressPackageStartupMessages(library(ape))",
  "set.seed(11)",
  "my_tree <- rtree(18, tip.label = paste0('Tip', 1:18))",
  "traits <- matrix(runif(18 * 6), 18, 6,",
  "  dimnames = list(my_tree$tip.label, paste0('Trait', 1:6)))",
  "```",
  "",
  chunk_lines("interactive-tree-viz-function"),
  "",
  chunk_lines("interactive-tree-viz-demo")
)
writeLines(min_rmd, "dev/tmp/test_rmd_min.Rmd", useBytes = TRUE)

if (!rmarkdown::pandoc_available())
  stop("pandoc is not available")
old_wd <- setwd("dev/tmp")
on.exit(setwd(old_wd), add = TRUE)
rmarkdown::render("test_rmd_min.Rmd", output_file = "test_rmd_min.html",
                  quiet = TRUE, envir = new.env())
stopifnot(file.exists("test_rmd_min.html"))
cat("Rendered dev/tmp/test_rmd_min.html:",
    file.info("test_rmd_min.html")$size, "bytes\n")

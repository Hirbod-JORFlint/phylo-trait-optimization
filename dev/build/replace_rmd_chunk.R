# Replace the contents of a knitr code chunk in the Rmd with a file's contents.
# Run from the repository root: Rscript dev/build/replace_rmd_chunk.R
rmd_path <- "diagnostics_on_tree_reorder.Rmd"
rmd_raw  <- readBin(rmd_path, "raw", file.info(rmd_path)$size)
eol      <- if (any(rmd_raw == as.raw(13L))) "\r\n" else "\n"
rmd  <- sub("\r$", "", readLines(rmd_path, warn = FALSE, encoding = "UTF-8"))
body <- sub("\r$", "", readLines("dev/build/interactive_tree_viz_final.R",
                                 warn = FALSE, encoding = "UTF-8"))

label <- "interactive-tree-viz-function"
start <- grep(paste0("^```\\{r ", label, "(,|\\})"), rmd)
stopifnot(length(start) == 1L)
ends <- which(seq_along(rmd) > start & grepl("^```\\s*$", rmd))
end <- ends[1]
stopifnot(length(end) == 1L, end > start)

out <- c(rmd[seq_len(start)], body, rmd[end:length(rmd)])
con <- file(rmd_path, open = "wb")
writeLines(out, con, sep = eol, useBytes = TRUE)
close(con)
cat("Replaced chunk", label, "lines", start, "-", end,
    "with", length(body), "lines; file now", length(out), "lines\n")

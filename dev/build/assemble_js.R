# Concatenate the widget JavaScript parts into tree_viz_widget.js.
# Run from the repository root: Rscript dev/build/assemble_js.R
parts <- c("dev/build/js_part1.js",
           "dev/build/js_part2.js",
           "dev/build/js_part3.js")
js <- unlist(lapply(seq_along(parts), function(i) {
  lines <- sub("\r$", "", readLines(parts[i], warn = FALSE, encoding = "UTF-8"))
  if (i < length(parts)) c(lines, "") else lines
}), use.names = FALSE)

out_path <- "dev/build/tree_viz_widget.js"
con <- file(out_path, open = "wb")
writeLines(js, con, sep = "\n", useBytes = TRUE)
close(con)
cat("Wrote", out_path, "with", length(js), "lines\n")

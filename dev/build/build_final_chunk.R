# Build a self-contained version of interactive_tree_viz_new.R by embedding the
# widget JavaScript as an R raw string.
# Run from the repository root: Rscript dev/build/build_final_chunk.R
src_path <- "dev/build/interactive_tree_viz_new.R"
src_raw  <- readBin(src_path, "raw", file.info(src_path)$size)
eol      <- if (any(src_raw == as.raw(13L))) "\r\n" else "\n"
new <- sub("\r$", "", readLines(src_path, warn = FALSE, encoding = "UTF-8"))
js  <- sub("\r$", "", readLines("dev/build/tree_viz_widget.js",
                                warn = FALSE, encoding = "UTF-8"))

start <- grep("^  # ----- build the plotly htmlwidget -----$", new)
end   <- grep("^  js_full <- gsub\\(\"__PAYLOAD__\"", new)
stopifnot(length(start) == 1L, length(end) == 1L, start < end)

if (any(grepl(')---"', js, fixed = TRUE)))
  stop("JS contains the raw-string terminator")

block <- c(
  "  # ----- build the plotly htmlwidget -----",
  "  # The widget JavaScript is embedded as a raw R string so this chunk is fully",
  "  # self-contained. __PAYLOAD__ is replaced by the per-dataset JSON payload.",
  '  js_body <- r"---(',
  js,
  ')---"',
  '  js_full <- gsub("__PAYLOAD__", payload_json, js_body, fixed = TRUE)'
)

out <- c(new[seq_len(start - 1L)], block, new[(end + 1L):length(new)])
out_path <- "dev/build/interactive_tree_viz_final.R"
con <- file(out_path, open = "wb")
writeLines(out, con, sep = eol, useBytes = TRUE)
close(con)
cat("Wrote dev/build/interactive_tree_viz_final.R with", length(out), "lines\n")

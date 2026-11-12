# Development workspace

This folder keeps the interactive-tree-viz build pipeline, tests, and debug
artifacts out of the repository root. The rendered R Markdown file and its
tracked demo outputs remain at the root.

## Layout

```text
dev/
  build/                 Widget JavaScript parts and R build scripts
  tests/                 R/Python/JS smoke and UI tests
  debug/                 Ignored debug HTML dumps, screenshots and one-off scripts
  outputs/demo/          Ignored archived demo exports
  archive/               Ignored local backups
  tmp/                   Ignored output directory for test runs
```

## Rebuild the embedded widget chunk

Run these from the repository root, in order:

```powershell
Rscript dev/build/assemble_js.R
node --check dev/build/tree_viz_widget.js
Rscript dev/build/build_final_chunk.R
Rscript dev/build/replace_rmd_chunk.R
```

`assemble_js.R` concatenates `js_part1.js`, `js_part2.js` and `js_part3.js`
into `tree_viz_widget.js`. `build_final_chunk.R` embeds it into
`dev/build/interactive_tree_viz_final.R`, and `replace_rmd_chunk.R` copies that
file into the `interactive-tree-viz-function` chunk of
`diagnostics_on_tree_reorder.Rmd`.

## Run the tests

All tests assume the working directory is the repository root. Generated files
are written to `dev/tmp/`, which is ignored by git.

```powershell
Rscript dev/tests/test_embedded.R
Rscript dev/tests/test_demo_chunk.R
Rscript dev/tests/test_helpers.R
Rscript dev/tests/test_itv.R
Rscript dev/tests/test_variants.R
Rscript dev/tests/make_min_rmd.R
```

`test_demo_chunk.R` expects `dev/tmp/extracted_function.R` from
`test_embedded.R`, so run `test_embedded.R` first.

For the DOM/UI test:

```powershell
Rscript dev/tests/test_itv.R
python dev/tests/make_ui_test.py
# then open dev/tmp/ui_test2.html in a headless browser
```

`probe_html.py` injects a DOM probe into any widget HTML:

```powershell
python dev/tests/probe_html.py dev/tmp/ui_test2.html dev/tmp/ui_test2_probe.html
```

## Ignored folders

`dev/debug/`, `dev/tmp/`, `dev/archive/` and `dev/outputs/demo/` are ignored in
`.gitignore`; they are for local diagnostics and generated files only.

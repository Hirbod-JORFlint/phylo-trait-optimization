# =============================================================================
# interactive_tree_viz()
# -----------------------------------------------------------------------------
# Plotly-based interactive visualization of an ape phylo tree + trait heatmap.
#   - Rectangular, slanted, circular and fan layouts
#   - Hover tips / internal nodes / heatmap cells for details
#   - Click an internal-node diamond to flip its child clades (topology safe)
#   - Optional OLO (Bar-Joseph et al. 2001) optimal leaf ordering
#   - Second view: tree + species names + values of a selected trait
#   - Appearance, titles, per-view sizing, filtering, highlighting, dark/light mode
#   - Export PNG / SVG / PDF (print), Newick, heatmap CSV and config JSON
#   - Built-in and user-saved configuration presets
# =============================================================================

# -----------------------------------------------------------------------------
# Self-contained OLO helpers (suffixed _2 to avoid clobbering the project's
# versions in 02_reorder.Rmd -- identical algorithm)
# -----------------------------------------------------------------------------
get_children_list_2 <- function(phy) {
  n_tip  <- length(phy$tip.label)
  n_node <- phy$Nnode
  children <- vector("list", n_tip + n_node)
  for (i in seq_len(nrow(phy$edge))) {
    p  <- phy$edge[i, 1]
    ch <- phy$edge[i, 2]
    children[[p]] <- c(children[[p]], ch)
  }
  children
}

olo_phylo_2 <- function(phy, d) {
  n_tip <- length(phy$tip.label)
  if (!is.matrix(d)) d <- as.matrix(d)
  stopifnot(nrow(d) == n_tip, ncol(d) == n_tip)
  children <- get_children_list_2(phy)
  memo <- vector("list", n_tip + phy$Nnode)
  dp <- function(node) {
    if (!is.null(memo[[node]])) return(memo[[node]])
    kids <- children[[node]]
    if (is.null(kids)) {
      res <- list(leaves = node,
                  table = matrix(0, 1, 1, dimnames = list(node, node)))
      memo[[node]] <<- res
      return(res)
    }
    if (length(kids) != 2) {
      stop("Node ", node, " has ", length(kids),
           " children; tree must be fully bifurcating (resolve with ape::multi2di).")
    }
    left  <- dp(kids[1])
    right <- dp(kids[2])
    Lleaves <- left$leaves;  Rleaves <- right$leaves
    LT <- left$table;        RT <- right$table
    L_full <- matrix(Inf, length(Lleaves), length(Lleaves),
                     dimnames = list(Lleaves, Lleaves))
    if (length(Lleaves) == 1) L_full[1, 1] <- 0 else {
      L_full[rownames(LT), colnames(LT)] <- LT
      L_full[colnames(LT), rownames(LT)] <- t(LT)
    }
    R_full <- matrix(Inf, length(Rleaves), length(Rleaves),
                     dimnames = list(Rleaves, Rleaves))
    if (length(Rleaves) == 1) R_full[1, 1] <- 0 else {
      R_full[rownames(RT), colnames(RT)] <- RT
      R_full[colnames(RT), rownames(RT)] <- t(RT)
    }
    tab  <- matrix(Inf, length(Lleaves), length(Rleaves),
                   dimnames = list(Lleaves, Rleaves))
    ptrB <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    ptrC <- matrix(NA_integer_, length(Lleaves), length(Rleaves))
    for (ui in seq_along(Lleaves)) {
      u <- Lleaves[ui]
      Lrow <- L_full[as.character(u), ]
      for (wi in seq_along(Rleaves)) {
        w <- Rleaves[wi]
        Rcol <- R_full[, as.character(w)]
        best <- Inf; bestB <- NA_integer_; bestC <- NA_integer_
        for (bi in seq_along(Lleaves)) {
          Lcost <- Lrow[bi]
          if (!is.finite(Lcost)) next
          b <- Lleaves[bi]
          total_vec <- Lcost + Rcol + d[as.character(b), as.character(Rleaves)]
          ci <- which.min(total_vec)
          if (total_vec[ci] < best) { best <- total_vec[ci]; bestB <- b; bestC <- Rleaves[ci] }
        }
        tab[ui, wi] <- best; ptrB[ui, wi] <- bestB; ptrC[ui, wi] <- bestC
      }
    }
    res <- list(leaves = c(Lleaves, Rleaves), table = tab,
                Lleaves = Lleaves, Rleaves = Rleaves,
                ptrB = ptrB, ptrC = ptrC, left = left, right = right)
    memo[[node]] <<- res
    res
  }
  root <- unique(setdiff(phy$edge[, 1], phy$edge[, 2]))
  stopifnot(length(root) == 1)
  top <- dp(root)
  get_order <- function(node_res, start, end) {
    if (length(node_res$leaves) == 1) return(node_res$leaves)
    Lleaves <- node_res$Lleaves; Rleaves <- node_res$Rleaves
    if (start %in% Lleaves && end %in% Rleaves) {
      ui <- match(start, Lleaves); wi <- match(end, Rleaves)
      b <- node_res$ptrB[ui, wi]; c_ <- node_res$ptrC[ui, wi]
      c(get_order(node_res$left, start, b), get_order(node_res$right, c_, end))
    } else if (start %in% Rleaves && end %in% Lleaves) {
      rev(get_order(node_res, end, start))
    } else stop("invalid boundary pair")
  }
  min_cost <- min(top$table)
  idx <- which(top$table == min_cost, arr.ind = TRUE)[1, ]
  best_order <- get_order(top, top$Lleaves[idx[1]], top$Rleaves[idx[2]])
  list(order = as.integer(best_order), cost = min_cost)
}

make_trait_dist_2 <- function(traits, phy, method = "euclidean", scale_traits = TRUE) {
  stopifnot(all(phy$tip.label %in% rownames(traits)))
  traits <- traits[phy$tip.label, , drop = FALSE]
  if (scale_traits) traits <- scale(traits)
  d <- as.matrix(dist(traits, method = method))
  rownames(d) <- colnames(d) <- as.character(seq_len(nrow(d)))
  d
}

# Local layout helpers (self-contained, rely only on ape edge matrix)
leaf_order_2 <- function(phy) {
  n <- length(phy$tip.label)
  phy$tip.label[phy$edge[phy$edge[, 2] <= n, 2]]
}

layout_cost_2 <- function(phy, traits) {
  ord <- leaf_order_2(phy)
  X <- traits[ord, , drop = FALSE]
  sum(sqrt(rowSums(diff(X)^2)))
}

# -----------------------------------------------------------------------------
# Configuration helpers
# -----------------------------------------------------------------------------
.itv_or <- function(x, y) if (is.null(x)) y else x

tree_viz_config <- function(
    layout = "rectangular",
    title = "Interactive phylogenetic tree + trait heatmap",
    subtitle = "",
    theme = "light",
    tip_labels = TRUE,
    tip_label_size = 10,
    show_internal_nodes = TRUE,
    internal_node_size = 9,
    show_internal_node_ids = FALSE,
    internal_node_id_size = 9,
    tip_node_color = "#1f77b4",
    tip_node_size = 7,
    internal_node_color = "#e08a1e",
    branch_color = "#8a8a8a",
    branch_width = 1.4,
    heatmap_colorscale = "Viridis",
    heatmap_range = NULL,
    heatmap_title = "Trait",
    heatmap_colorbar_title = "value",
    heatmap_panel_title = "Trait heatmap",
    tree_panel_title = "Phylogenetic tree",
    species_panel_title = "Species",
    trait_panel_title = "",
    responsive = TRUE,
    width = NULL,
    height = NULL,
    aspect_ratio = NULL,
    heatmap_view_width = NULL,
    heatmap_view_height = NULL,
    trait_view_width = NULL,
    trait_view_height = NULL,
    tip_filter = NULL,
    clade = NULL,
    heatmap_columns = NULL,
    highlight = NULL,
    fan_angle = 300,
    show_heatmap = TRUE
) {
  list(
    layout = layout, title = title, subtitle = subtitle, theme = theme,
    tip_labels = tip_labels, tip_label_size = tip_label_size,
    show_internal_nodes = show_internal_nodes,
    internal_node_size = internal_node_size,
    show_internal_node_ids = show_internal_node_ids,
    internal_node_id_size = internal_node_id_size,
    tip_node_color = tip_node_color, tip_node_size = tip_node_size,
    internal_node_color = internal_node_color,
    branch_color = branch_color, branch_width = branch_width,
    heatmap_colorscale = heatmap_colorscale,
    heatmap_range = heatmap_range,
    heatmap_title = heatmap_title,
    heatmap_colorbar_title = heatmap_colorbar_title,
    heatmap_panel_title = heatmap_panel_title,
    tree_panel_title = tree_panel_title,
    species_panel_title = species_panel_title,
    trait_panel_title = trait_panel_title,
    responsive = responsive, width = width, height = height,
    aspect_ratio = aspect_ratio,
    heatmap_view_width = heatmap_view_width,
    heatmap_view_height = heatmap_view_height,
    trait_view_width = trait_view_width,
    trait_view_height = trait_view_height,
    tip_filter = tip_filter, clade = clade,
    heatmap_columns = heatmap_columns, highlight = highlight,
    fan_angle = fan_angle, show_heatmap = show_heatmap
  )
}

merge_tree_viz_config <- function(base, override) {
  if (is.null(override)) return(base)
  if (!is.list(override)) stop("A configuration must be a list.", call. = FALSE)
  for (nm in names(override)) {
    if (is.null(nm) || !nzchar(nm)) next
    base[[nm]] <- override[[nm]]
  }
  base
}

normalize_tree_viz_config <- function(cfg, traits, tree) {
  defaults <- tree_viz_config()
  cfg <- merge_tree_viz_config(defaults, cfg)
  trait_names <- colnames(traits)
  if (is.null(trait_names)) trait_names <- paste0("Trait", seq_len(ncol(traits)))
  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  chr1 <- function(x, d) {
    if (is.null(x) || length(x) != 1L || is.na(x)) d else as.character(x)
  }
  num1 <- function(x, d) {
    if (is.null(x) || length(x) != 1L || !is.finite(suppressWarnings(as.numeric(x)))) d else as.numeric(x)
  }
  lgl1 <- function(x, d) if (is.null(x) || length(x) != 1L || is.na(x)) d else isTRUE(as.logical(x))

  cfg$layout <- match.arg(chr1(cfg$layout, "rectangular"),
                          c("rectangular", "slanted", "circular", "fan"))
  cfg$theme <- match.arg(chr1(cfg$theme, "light"), c("light", "dark"))
  cfg$title <- chr1(cfg$title, defaults$title)
  cfg$subtitle <- chr1(cfg$subtitle, "")
  cfg$tip_labels <- lgl1(cfg$tip_labels, TRUE)
  cfg$show_internal_nodes <- lgl1(cfg$show_internal_nodes, TRUE)
  cfg$show_internal_node_ids <- lgl1(cfg$show_internal_node_ids, FALSE)
  cfg$show_heatmap <- lgl1(cfg$show_heatmap, TRUE)
  cfg$tip_label_size <- num1(cfg$tip_label_size, 10)
  cfg$internal_node_size <- num1(cfg$internal_node_size, 9)
  cfg$internal_node_id_size <- num1(cfg$internal_node_id_size, 9)
  cfg$tip_node_size <- num1(cfg$tip_node_size, 7)
  cfg$branch_width <- num1(cfg$branch_width, 1.4)
  cfg$fan_angle <- max(10, min(360, num1(cfg$fan_angle, 300)))
  cfg$tip_node_color <- chr1(cfg$tip_node_color, "#1f77b4")
  cfg$internal_node_color <- chr1(cfg$internal_node_color, "#e08a1e")
  cfg$branch_color <- chr1(cfg$branch_color, "#8a8a8a")
  cfg$heatmap_colorscale <- chr1(cfg$heatmap_colorscale, "Viridis")
  cfg$heatmap_title <- chr1(cfg$heatmap_title, "Trait")
  cfg$heatmap_colorbar_title <- chr1(cfg$heatmap_colorbar_title, "value")
  cfg$heatmap_panel_title <- chr1(cfg$heatmap_panel_title, "Trait heatmap")
  cfg$tree_panel_title <- chr1(cfg$tree_panel_title, "Phylogenetic tree")
  cfg$species_panel_title <- chr1(cfg$species_panel_title, "Species")
  cfg$trait_panel_title <- chr1(cfg$trait_panel_title, "")

  hr <- cfg$heatmap_range
  if (!is.null(hr) && length(hr) == 2L && all(is.finite(suppressWarnings(as.numeric(hr))))) {
    cfg$heatmap_range <- sort(as.numeric(hr))
  } else cfg$heatmap_range <- NULL

  tf <- cfg$tip_filter
  if (is.character(tf)) tf <- list(tips = tf)
  if (is.logical(tf) && !is.null(names(tf))) tf <- list(tips = names(tf)[tf])
  if (is.list(tf)) {
    out <- list()
    if (!is.null(tf$trait) && length(tf$trait) == 1L && tf$trait %in% trait_names) {
      out$trait <- as.character(tf$trait)
      mn <- suppressWarnings(as.numeric(tf$min)); mx <- suppressWarnings(as.numeric(tf$max))
      out$min <- if (length(mn) == 1L && is.finite(mn)) mn else NULL
      out$max <- if (length(mx) == 1L && is.finite(mx)) mx else NULL
    } else if (!is.null(tf$trait)) {
      warning("tip_filter$trait not found in trait columns; ignoring trait filter.", call. = FALSE)
    }
    if (!is.null(tf$tips)) {
      tips_keep <- intersect(as.character(tf$tips), tree$tip.label)
      if (length(tips_keep)) out$tips <- tips_keep
    }
    cfg$tip_filter <- if (length(out)) out else NULL
  } else cfg$tip_filter <- NULL

  cl <- cfg$clade
  if (!is.null(cl)) {
    cl <- suppressWarnings(as.integer(cl))
    cl <- cl[is.finite(cl) & cl >= n_tip + 1L & cl <= n_tip + n_node]
    cfg$clade <- if (length(cl)) cl else NULL
  }

  hc <- cfg$heatmap_columns
  if (!is.null(hc)) {
    if (is.numeric(hc)) {
      hc <- trait_names[hc[hc >= 1 & hc <= length(trait_names)]]
    } else {
      hc <- intersect(as.character(hc), trait_names)
    }
    cfg$heatmap_columns <- hc
  }

  hl <- cfg$highlight
  if (is.character(hl)) hl <- list(tips = hl)
  if (is.logical(hl) && !is.null(names(hl))) hl <- list(tips = names(hl)[hl])
  if (is.list(hl)) {
    out <- list(color = chr1(hl$color, "#e41a1c"))
    if (!is.null(hl$trait) && length(hl$trait) == 1L && hl$trait %in% trait_names) {
      out$trait <- as.character(hl$trait)
      mn <- suppressWarnings(as.numeric(hl$min)); mx <- suppressWarnings(as.numeric(hl$max))
      out$min <- if (length(mn) == 1L && is.finite(mn)) mn else NULL
      out$max <- if (length(mx) == 1L && is.finite(mx)) mx else NULL
    }
    if (!is.null(hl$tips)) {
      tips_hl <- intersect(as.character(hl$tips), tree$tip.label)
      if (length(tips_hl)) out$tips <- tips_hl
    }
    cfg$highlight <- if (length(out) > 1L) out else NULL
  } else cfg$highlight <- NULL

  cfg$responsive <- lgl1(cfg$responsive, TRUE)
  pos1 <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) == 1L && is.finite(v) && v > 0) v else NULL
  }
  cfg$width <- pos1(cfg$width)
  cfg$height <- pos1(cfg$height)
  cfg$aspect_ratio <- pos1(cfg$aspect_ratio)
  cfg$heatmap_view_width <- pos1(cfg$heatmap_view_width)
  cfg$heatmap_view_height <- pos1(cfg$heatmap_view_height)
  cfg$trait_view_width <- pos1(cfg$trait_view_width)
  cfg$trait_view_height <- pos1(cfg$trait_view_height)
  if (!cfg$responsive && is.null(cfg$width) && is.null(cfg$height)) {
    cfg$width <- 1000; cfg$height <- 700
  }
  if (!is.null(cfg$aspect_ratio)) {
    if (!is.null(cfg$width) && is.null(cfg$height)) cfg$height <- cfg$width / cfg$aspect_ratio
    if (is.null(cfg$width) && !is.null(cfg$height)) cfg$width <- cfg$height * cfg$aspect_ratio
  }
  cfg
}

# Built-in visualization presets -------------------------------------------------
tree_viz_presets <- local({
  base <- tree_viz_config()
  list(
    default = base,
    publication = merge_tree_viz_config(base, list(
      title = "Phylogenetic tree and trait heatmap",
      tip_label_size = 8, tip_node_size = 4, internal_node_size = 6,
      branch_width = 1, heatmap_colorscale = "Cividis",
      width = 1200, height = 900, responsive = FALSE
    )),
    dark = merge_tree_viz_config(base, list(
      theme = "dark", heatmap_colorscale = "Plasma"
    )),
    compact = merge_tree_viz_config(base, list(
      tip_labels = FALSE, show_internal_nodes = FALSE,
      height = 500, heatmap_panel_title = ""
    )),
    colorblind = merge_tree_viz_config(base, list(
      heatmap_colorscale = "Cividis",
      tip_node_color = "#0072B2", internal_node_color = "#D55E00",
      branch_color = "#555555", highlight = list(color = "#009E73")
    ))
  )
})

tree_viz_preset <- function(name = "default", ...) {
  if (!is.character(name) || length(name) != 1L || !name %in% names(tree_viz_presets)) {
    stop("Unknown preset '", name, "'. Available: ",
         paste(names(tree_viz_presets), collapse = ", "), call. = FALSE)
  }
  merge_tree_viz_config(tree_viz_presets[[name]], list(...))
}

save_tree_viz_config <- function(config, file) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required.", call. = FALSE)
  jsonlite::write_json(config, file, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", digits = NA)
  invisible(file)
}

load_tree_viz_config <- function(file) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required.", call. = FALSE)
  x <- jsonlite::read_json(file, simplifyVector = TRUE)
  merge_tree_viz_config(tree_viz_config(), x)
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================
interactive_tree_viz <- function(
    tree,                       # ape::phylo object
    traits,                     # numeric matrix, rownames == tree$tip.label
    reorder_olo = TRUE,         # reorder with OLO (Bar-Joseph) before plotting
    scale_traits = TRUE,        # scale traits for the OLO distance / cost readout
    output_file = "interactive_tree.html",
    title = NULL,
    save_newick = NULL,         # optional path to write the (possibly reordered) tree
    selfcontained = TRUE,       # inline plotly.js bundle in the HTML
    config = NULL,              # list, or name of a built-in preset
    preset = NULL,              # name of a built-in preset (or a config list)
    save_heatmap_csv = NULL,    # optional path for the heatmap data
    save_config = NULL,         # optional path for the configuration JSON
    ...                         # any tree_viz_config() field, e.g. layout = "circular"
) {
  for (p in c("ape", "plotly", "htmlwidgets", "jsonlite")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required.", call. = FALSE)
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be an ape 'phylo' object.")
  if (is.data.frame(traits)) {
    if (!all(vapply(traits, is.numeric, logical(1))))
      stop("'traits' must be numeric.")
    traits <- as.matrix(traits)
  }
  if (!is.matrix(traits)) stop("'traits' must be a numeric matrix or data.frame.")
  if (is.null(rownames(traits))) stop("'traits' must have rownames matching tree$tip.label.")
  missing_tips <- setdiff(tree$tip.label, rownames(traits))
  if (length(missing_tips))
    stop("Traits missing for tip(s): ", paste(missing_tips, collapse = ", "))
  n <- length(tree$tip.label)
  if (n < 2) stop("Tree must have at least 2 tips.")
  traits <- traits[tree$tip.label, , drop = FALSE]
  if (is.null(colnames(traits))) colnames(traits) <- paste0("Trait", seq_len(ncol(traits)))

  # ----- resolve the initial configuration -----
  defaults <- tree_viz_config()
  cfg <- defaults
  if (!is.null(preset)) {
    if (is.character(preset) && length(preset) == 1L) {
      if (!preset %in% names(tree_viz_presets))
        stop("Unknown preset '", preset, "'. Available: ",
             paste(names(tree_viz_presets), collapse = ", "), call. = FALSE)
      cfg <- merge_tree_viz_config(cfg, tree_viz_presets[[preset]])
    } else if (is.list(preset)) {
      cfg <- merge_tree_viz_config(cfg, preset)
    } else stop("'preset' must be a preset name or a configuration list.", call. = FALSE)
  }
  if (!is.null(config)) {
    if (is.character(config) && length(config) == 1L) {
      if (!config %in% names(tree_viz_presets))
        stop("Unknown preset '", config, "'. Available: ",
             paste(names(tree_viz_presets), collapse = ", "), call. = FALSE)
      cfg <- merge_tree_viz_config(cfg, tree_viz_presets[[config]])
    } else if (is.list(config)) {
      cfg <- merge_tree_viz_config(cfg, config)
    } else stop("'config' must be a preset name or a configuration list.", call. = FALSE)
  }
  extra <- list(...)
  if (length(extra)) {
    unknown <- setdiff(names(extra), names(defaults))
    if (length(unknown))
      warning("Unknown configuration field(s) ignored: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    extra <- extra[intersect(names(extra), names(defaults))]
    cfg <- merge_tree_viz_config(cfg, extra)
  }
  if (!is.null(title)) cfg$title <- title
  cfg <- normalize_tree_viz_config(cfg, traits, tree)

  # scale guard: constant columns produce NaN after scale()
  if (scale_traits) {
    sds <- apply(traits, 2, stats::sd)
    if (any(!is.finite(sds) | sds == 0)) {
      warning("Constant trait column(s) detected; disabling scaling for OLO.", call. = FALSE)
      scale_traits <- FALSE
    }
  }

  olo_order <- NULL
  cost_mat  <- traits
  out_cost  <- layout_cost_2(tree, traits)

  if (reorder_olo) {
    node_counts <- table(tree$edge[, 1])
    if (any(node_counts > 2))
      stop("OLO requires a fully bifurcating tree; resolve polytomies with ape::multi2di().")
    d <- make_trait_dist_2(traits, tree, method = "euclidean", scale_traits = scale_traits)
    res <- olo_phylo_2(tree, d)
    ord_labels <- tree$tip.label[res$order]
    tree <- ape::rotateConstr(tree, ord_labels)
    olo_order <- ord_labels
    cost_mat  <- if (scale_traits) scale(traits) else traits
    out_cost  <- res$cost
    cat("OLO minimal adjacent trait-dissimilarity:", out_cost, "\n")
  }

  if (!is.null(save_newick)) {
    ape::write.tree(tree, file = save_newick)
    cat("Step tree written to:", save_newick, "\n")
  }
  if (!is.null(save_heatmap_csv)) {
    ord <- leaf_order_2(tree)
    hm <- data.frame(tip = ord, traits[ord, , drop = FALSE],
                     check.names = FALSE, stringsAsFactors = FALSE)
    utils::write.csv(hm, save_heatmap_csv, row.names = FALSE)
    cat("Heatmap data written to:", save_heatmap_csv, "\n")
  }
  if (!is.null(save_config)) {
    save_tree_viz_config(cfg, save_config)
    cat("Configuration written to:", save_config, "\n")
  }

  # ----- serialize tree structure + data for the browser -----
  children <- get_children_list_2(tree)
  n_node <- tree$Nnode
  node_count <- n + n_node
  edge_len <- rep(1, node_count)
  if (!is.null(tree$edge.length)) edge_len[tree$edge[, 2]] <- tree$edge.length
  root <- unique(setdiff(tree$edge[, 1], tree$edge[, 2]))
  root <- root[1]

  child_map <- list()
  for (nd in (n + 1):node_count) child_map[[as.character(nd)]] <- children[[nd]]

  payload <- list(
    tips       = tree$tip.label,
    children   = child_map,
    edgeLen    = as.list(stats::setNames(edge_len, as.character(seq_len(node_count)))),
    root       = root,
    nodeCount  = node_count,
    traits     = lapply(seq_len(n), function(i) as.numeric(traits[i, ])),
    traitNames = colnames(traits),
    costMat    = lapply(seq_len(n), function(i) as.numeric(cost_mat[i, ])),
    oloOrder   = olo_order,
    title      = cfg$title,
    config     = cfg,
    defaults   = defaults,
    presets    = tree_viz_presets,
    costScaled = isTRUE(scale_traits),
    outputFile = basename(output_file)
  )
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = NA,
                                   null = "null")

  # ----- build the plotly htmlwidget -----
  # For development the JS lives in dev/build/tree_viz_widget.js. The generated
  # Rmd chunk embeds the same JS as a raw string (see build_final_chunk.R).
  js_path <- "dev/build/tree_viz_widget.js"
  if (!file.exists(js_path)) js_path <- "tree_viz_widget.js"
  js_body <- if (file.exists(js_path)) {
    paste(readLines(js_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    stop("tree_viz_widget.js not found (expected in dev/build/).", call. = FALSE)
  }
  js_full <- gsub("__PAYLOAD__", payload_json, js_body, fixed = TRUE)

  fig <- suppressWarnings(suppressMessages(plotly::plot_ly()))
  widget <- suppressWarnings(suppressMessages(
    htmlwidgets::onRender(fig, paste0("function(el, x) {\n", js_full, "\n}"))
  ))

  # ----- responsive / fixed / aspect-ratio sizing -----
  # The widget JS sizes the plot div inside a flex column (control panel on
  # top, plot filling the rest). Each view (heatmap vs species+trait) can have
  # its own width/height; the initial view is the heatmap view.
  initial_width  <- .itv_or(cfg$heatmap_view_width, cfg$width)
  initial_height <- .itv_or(cfg$heatmap_view_height, cfg$height)
  if (!is.null(cfg$aspect_ratio)) {
    if (!is.null(initial_width) && is.null(initial_height))
      initial_height <- initial_width / cfg$aspect_ratio
    if (is.null(initial_width) && !is.null(initial_height))
      initial_width <- initial_height * cfg$aspect_ratio
  }
  widget$width <- .itv_or(initial_width, "100%")
  widget$height <- .itv_or(initial_height, "100%")

  suppressWarnings(suppressMessages(
    htmlwidgets::saveWidget(widget, output_file, selfcontained = selfcontained)
  ))

  cat("Saved interactive tree visualization to:", output_file, "\n")

  invisible(list(
    tree         = tree,
    order        = leaf_order_2(tree),
    cost         = out_cost,
    olo          = olo_order,
    file         = output_file,
    widget       = widget,
    config       = cfg,
    heatmap      = traits,
    heatmap_order = traits[leaf_order_2(tree), , drop = FALSE],
    newick       = save_newick,
    heatmap_csv  = save_heatmap_csv,
    config_json  = save_config
  ))
}

# -----------------------------------------------------------------------------
# Export helpers
# -----------------------------------------------------------------------------
export_tree_viz_data <- function(result, dir = ".", prefix = "tree_viz") {
  if (!is.list(result) || is.null(result$tree))
    stop("'result' must be the object returned by interactive_tree_viz().", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  tre <- file.path(dir, paste0(prefix, ".tre"))
  csv <- file.path(dir, paste0(prefix, "_heatmap.csv"))
  cfg <- file.path(dir, paste0(prefix, "_config.json"))
  ape::write.tree(result$tree, file = tre)
  hm <- result$heatmap_order
  hm <- data.frame(tip = rownames(hm), hm, check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(hm, csv, row.names = FALSE)
  save_tree_viz_config(result$config, cfg)
  invisible(c(tree = tre, heatmap_csv = csv, config_json = cfg))
}

export_tree_viz_image <- function(html_file, output_file,
                                  format = c("png", "pdf", "svg"),
                                  width = 1200, height = 800, delay = 2) {
  format <- match.arg(format)
  if (!file.exists(html_file)) stop("HTML file not found: ", html_file, call. = FALSE)
  url <- paste0("file:///", normalizePath(html_file, winslash = "/", mustWork = TRUE))
  if (requireNamespace("webshot2", quietly = TRUE)) {
    if (format == "svg")
      stop("webshot2 cannot export SVG; use format = 'png' or 'pdf'.", call. = FALSE)
    webshot2::webshot(url, file = output_file, vwidth = width, vheight = height,
                      delay = delay, zoom = 1, selector = ".js-plotly-plot")
    return(invisible(output_file))
  }
  if (!requireNamespace("chromote", quietly = TRUE))
    stop("Install 'webshot2' (recommended) or 'chromote' to export images.", call. = FALSE)
  b <- chromote::ChromoteSession$new()
  on.exit(b$close(), add = TRUE)
  b$set_viewport_size(width, height)
  b$Page$navigate(url)
  b$Page$loadEventFired()
  Sys.sleep(delay)
  if (format == "pdf") {
    b$Page$printToPDF(landscape = width > height, printBackground = TRUE,
                      paperWidth = width / 96, paperHeight = height / 96,
                      marginTop = 0, marginBottom = 0, marginLeft = 0, marginRight = 0,
                      out = output_file)
  } else if (format == "png") {
    b$screenshot(output_file, cliprect = list(x = 0, y = 0, width = width, height = height, scale = 1))
  } else {
    js <- sprintf("(function(){var el=document.querySelector('.js-plotly-plot');return Plotly.toImage(el,{format:'svg',width:%d,height:%d});})()",
                  width, height)
    res <- b$Runtime$evaluate(expression = js, awaitPromise = TRUE, returnByValue = TRUE)
    data_url <- res$result$value
    svg <- rawToChar(jsonlite::base64_dec(sub("^data:image/svg\\+xml;base64,", "", data_url)))
    writeLines(svg, output_file, useBytes = TRUE)
  }
  invisible(output_file)
}

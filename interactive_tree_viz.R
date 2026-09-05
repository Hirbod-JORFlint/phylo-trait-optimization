# =============================================================================
# interactive_tree_viz()
# -----------------------------------------------------------------------------
# Plotly-based interactive visualization of an ape phylo tree + trait heatmap.
#   - Tree (rectangular phylogram) on the left, aligned trait heatmap on right
#   - Hover tips / internal nodes / heatmap cells for details
#   - Click an internal-node diamond to flip (rotate) its two child clades
#     without changing the topology
#   - Optional OLO (Bar-Joseph et al. 2001) optimal leaf ordering before plot,
#     with a button to restore the OLO-optimal order after manual edits
#   - Output: a standalone interactive HTML file (plotly.js htmlwidget)
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
    selfcontained = TRUE        # inline plotly.js bundle in the HTML
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
  missing_tips <- setdiff(tree$tip.label, rownames(traits))
  if (length(missing_tips))
    stop("Traits missing for tip(s): ", paste(missing_tips, collapse = ", "))
  n <- length(tree$tip.label)
  if (n < 2) stop("Tree must have at least 2 tips.")
  traits <- traits[tree$tip.label, , drop = FALSE]
  if (is.null(title)) title <- "Interactive phylogenetic tree + trait heatmap"
  if (is.null(colnames(traits))) colnames(traits) <- paste0("Trait", seq_len(ncol(traits)))

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
    if (!requireNamespace("ape", quietly = TRUE)) stop("ape required for OLO.")
    node_counts <- table(tree$edge[, 1])
    if (any(node_counts > 2))
      stop("OLO requires a fully bifurcating tree; resolve polytomies with ape::multi2di().")
    d <- make_trait_dist_2(traits, tree, method = "euclidean", scale_traits = scale_traits)
    res <- olo_phylo_2(tree, d)
    ord_labels <- tree$tip.label[res$order]
    # topology-preserving rotation into the OLO-optimal layout
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

  # ----- serialize tree structure + data for the browser -----
  children <- get_children_list_2(tree)
  n_node <- tree$Nnode
  node_count <- n + n_node
  edge_len <- rep(1, node_count)
  if (!is.null(tree$edge.length)) edge_len[tree$edge[, 2]] <- tree$edge.length
  root <- unique(setdiff(tree$edge[, 1], tree$edge[, 2]))

  # only internal nodes are rotatable -> send those (tips have no children)
  child_map <- list()
  for (nd in (n + 1):node_count) child_map[[as.character(nd)]] <- children[[nd]]
  root <- root[1]

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
    title      = title
  )
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = NA)

  # ----- build the plotly htmlwidget -----
  .js_body <- r"-(
var TD = __PAYLOAD__;
var tips = TD.tips;
var children = TD.children;
var edgeLen = TD.edgeLen;
var root = TD.root;
var traitsHeat = TD.traits;
var traitsCost = TD.costMat || TD.traits;
var tNames = TD.traitNames;
var n = tips.length;
var oloOrder = TD.oloOrder || null;
var ttl = TD.title || "Phylogenetic tree";

var rev = {};
function initRev(){ rev = {}; for (var k in children) rev[k] = false; }
initRev();

function leafOrder(r){
  var order = [];
  (function dfs(node){
    var ch = children[node];
    if(!ch){ order.push(node*1); return; }
    var kk = r[node] ? ch.slice().reverse() : ch;
    for(var i=0;i<kk.length;i++) dfs(kk[i]);
  })(root);
  return order;
}

function leafCount(node){
  var ch = children[node];
  if(!ch) return 1;
  var s = 0;
  for(var i=0;i<ch.length;i++) s += leafCount(ch[i]);
  return s;
}

function layoutCalc(r){
  var order = leafOrder(r);
  var y = {}, slot = {};
  for(var i=0;i<order.length;i++){ slot[order[i]] = i+1; y[order[i]] = i+1; }
  var mn = {}, mx = {};
  (function ranges(node){
    var ch = children[node];
    if(!ch){ mn[node]=slot[node]; mx[node]=slot[node]; return; }
    var lo = Infinity, hi = -Infinity;
    for(var i=0;i<ch.length;i++){
      ranges(ch[i]);
      if(mn[ch[i]] < lo) lo = mn[ch[i]];
      if(mx[ch[i]] > hi) hi = mx[ch[i]];
    }
    mn[node] = lo; mx[node] = hi;
    y[node] = (lo + hi) / 2;
  })(root);
  var x = {}; x[root] = 0;
  (function xs(node){
    var ch = children[node] || [];
    for(var i=0;i<ch.length;i++){
      x[ch[i]] = x[node] + (edgeLen[ch[i]] || 1);
      xs(ch[i]);
    }
  })(root);
  return {x:x, y:y, order:order};
}

function edgeTrace(L){
  var ex = [], ey = [];
  (function go(node){
    var ch = children[node] || [];
    for(var i=0;i<ch.length;i++){
      var c = ch[i];
      ex.push(L.x[node], L.x[c], null); ey.push(L.y[node], L.y[node], null);
      ex.push(L.x[c], L.x[c], null);    ey.push(L.y[node], L.y[c], null);
      go(c);
    }
  })(root);
  return {type:"scatter", mode:"lines", x:ex, y:ey,
          line:{color:"#8a8a8a", width:1.4}, hoverinfo:"skip", showlegend:false};
}

function tipTrace(L){
  var x=[], y=[], txt=[];
  L.order.forEach(function(t){ x.push(L.x[t]); y.push(L.y[t]); txt.push(tips[t-1]); });
  var tr = {type:"scatter", mode:"markers", x:x, y:y, text:txt,
            hovertemplate:"<b>%{text}</b><extra>tip</extra>",
            marker:{size:7, color:"#1f77b4", line:{width:1, color:"#ffffff"}},
            showlegend:false};
  if(n <= 40){ tr.mode = "markers+text"; tr.textposition = "middle right"; tr.textfont = {size:10}; }
  return tr;
}

function intTrace(L){
  var x=[], y=[], cd=[], tt=[];
  var keys = [];
  for(var k in children) keys.push(+k);
  keys.sort(function(a,b){ return a-b; });
  keys.forEach(function(nd){
    x.push(L.x[nd]); y.push(L.y[nd]); cd.push(nd);
    tt.push(leafCount(nd) + " tips in this clade - click to rotate");
  });
  return {type:"scatter", mode:"markers", x:x, y:y, customdata:cd, text:tt,
          hovertemplate:"<b>Internal node %{customdata}</b><br>%{text}<extra></extra>",
          marker:{symbol:"diamond", size:9, color:"#e08a1e", line:{width:1, color:"#9a5a00"}},
          showlegend:false};
}

function heatTrace(L){
  var Z = [];
  L.order.forEach(function(t){ Z.push(traitsHeat[t-1]); });
  return {type:"heatmap", z:Z,
          y:L.order.map(function(t){ return tips[t-1]; }),
          x:tNames, xaxis:"x2", yaxis:"y2",
          colorscale:"Viridis",
          colorbar:{thickness:10, len:0.6, y:0.5, title:{text:"value"}},
          hovertemplate:"<b>%{y}</b><br>%{x}: %{z:.3f}<extra></extra>",
          showlegend:false};
}

var lay = {
  title:{text: ttl},
  xaxis:{domain:[0,0.40], showgrid:false, zeroline:false, showticklabels:false},
  yaxis:{domain:[0.05,0.96], range:[0.5,n+0.5], showgrid:false, zeroline:false, showticklabels:false},
  xaxis2:{domain:[0.50,0.90], showgrid:false, title:{text:"Trait"}, tickangle:-30, automargin:true},
  yaxis2:{domain:[0.05,0.96], range:[0.5,n+0.5], showgrid:false, showticklabels:false},
  margin:{l:30,r:20,t:70,b:60},
  showlegend:false, hovermode:"closest"
};

var costEl = null;
function updateCost(L){
  var s = 0;
  for(var i=1;i<L.order.length;i++){
    var a = traitsCost[L.order[i-1]-1];
    var b = traitsCost[L.order[i]-1];
    var d2 = 0;
    for(var j=0;j<a.length;j++){ var dd = a[j]-b[j]; d2 += dd*dd; }
    s += Math.sqrt(d2);
  }
  if(costEl) costEl.textContent = s.toFixed(4);
}

function revForTarget(targetLabels){
  var tIdx = {};
  for(var i=0;i<targetLabels.length;i++) tIdx[targetLabels[i]] = i;
  var base = leafOrder({});
  var pos = {};
  base.forEach(function(t){ pos[t] = tIdx[tips[t-1]]; });
  var mn = {};
  (function calc(node){
    var ch = children[node];
    if(!ch){ mn[node] = pos[node]; return; }
    var lo = Infinity;
    for(var i=0;i<ch.length;i++){ calc(ch[i]); if(mn[ch[i]] < lo) lo = mn[ch[i]]; }
    mn[node] = lo;
  })(root);
  var r = {};
  (function orient(node){
    var ch = children[node];
    if(!ch) return;
    if(ch.length === 2) r[node] = mn[ch[1]] < mn[ch[0]];
    for(var i=0;i<ch.length;i++) orient(ch[i]);
  })(root);
  return r;
}

function redraw(){
  var L = layoutCalc(rev);
  Plotly.react(el, [edgeTrace(L), tipTrace(L), intTrace(L), heatTrace(L)], lay);
  updateCost(L);
}

var ui = document.createElement("div");
ui.style.cssText = "font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:13px;padding:6px 2px;line-height:1.7;";
var html = "<b>Adjacent trait distance:</b> <span style=\"color:#1f77b4;font-weight:700\" id=\"olocostLbl\">0</span>";
html += "<button id=\"oloresetBtn\" style=\"margin-left:16px\">Reset</button>";
if(oloOrder){ html += " <button id=\"olobestBtn\" title=\"Apply the OLO-optimal tip order\">Apply OLO-optimal order</button>"; }
html += "&nbsp; <span style=\"color:#777\">Click a <b style=\"color:#e08a1e\">diamond</b> (internal node) to flip its two child clades. Hover tips, internal nodes, or heatmap cells for details.</span>";
ui.innerHTML = html;
el.parentNode.insertBefore(ui, el);
costEl = document.getElementById("olocostLbl");
document.getElementById("oloresetBtn").addEventListener("click", function(){ initRev(); redraw(); });
if(oloOrder){
  document.getElementById("olobestBtn").addEventListener("click", function(){ rev = revForTarget(oloOrder); redraw(); });
}

el.on("plotly_click", function(data){
  var pts = data.points || [];
  for(var i=0;i<pts.length;i++){
    var p = pts[i];
    if(p.curveNumber === 2 && p.customdata !== undefined){
      var node = p.customdata;
      rev[node] = !rev[node];
      redraw();
      return;
    }
  }
});

redraw();
)-"

js_full <- gsub("__PAYLOAD__", payload_json, .js_body, fixed = TRUE)

  fig <- suppressWarnings(suppressMessages(plotly::plot_ly()))
  widget <- suppressWarnings(suppressMessages(
    htmlwidgets::onRender(fig, paste0("function(el, x) {\n", js_full, "\n}"))
  ))
  suppressWarnings(suppressMessages(
    htmlwidgets::saveWidget(widget, output_file, selfcontained = selfcontained)
  ))

  cat("Saved interactive tree visualization to:", output_file, "\n")

  invisible(list(
    tree  = tree,
    order = leaf_order_2(tree),
    cost  = out_cost,
    olo   = olo_order,
    file  = output_file,
    widget = widget
  ))
}

# -----------------------------------------------------------------------------
# quick smoke test if run directly
# -----------------------------------------------------------------------------
if (sys.nframe() == 0L && interactive()) {
  library(ape)
  set.seed(42)
  n <- 20
  tr <- rtree(n, tip.label = paste0("Tip", 1:n))
  Td <- matrix(runif(n * 6), nrow = n,
               dimnames = list(tr$tip.label, paste0("Trait", 1:6)))
  out <- interactive_tree_viz(tr, Td, reorder_olo = TRUE,
                              output_file = "test_interactive.html")
}
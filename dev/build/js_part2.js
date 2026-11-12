function edgeTrace(L) {
  var ex = [], ey = [];
  var radial = (cfg.layout === "circular" || cfg.layout === "fan");
  if (radial) {
    var r = L.r, ang = L.angle;
    (function go(node) {
      var ch = L.children[node] || [];
      if (!ch.length) return;
      var pr = r[node];
      for (var i = 0; i < ch.length; i++) {
        var c = ch[i], ac = ang[c], rc = r[c];
        ex.push(pr * Math.cos(ac), rc * Math.cos(ac), null);
        ey.push(pr * Math.sin(ac), rc * Math.sin(ac), null);
      }
      var a0 = ang[ch[0]], a1 = ang[ch[ch.length - 1]];
      if (a1 < a0) a1 += 2 * Math.PI;
      var steps = Math.max(6, Math.ceil((a1 - a0) / 0.06));
      for (var k = 0; k <= steps; k++) {
        var a = a0 + (a1 - a0) * k / steps;
        ex.push(pr * Math.cos(a)); ey.push(pr * Math.sin(a));
      }
      ex.push(null); ey.push(null);
      for (var j = 0; j < ch.length; j++) go(ch[j]);
    })(L.root);
  } else {
    var rect = (cfg.layout === "rectangular");
    (function go2(node) {
      var ch = L.children[node] || [];
      for (var i = 0; i < ch.length; i++) {
        var c = ch[i];
        if (rect) {
          ex.push(L.x[node], L.x[c], null); ey.push(L.y[node], L.y[node], null);
          ex.push(L.x[c], L.x[c], null); ey.push(L.y[node], L.y[c], null);
        } else {
          ex.push(L.x[node], L.x[c], null); ey.push(L.y[node], L.y[c], null);
        }
        go2(c);
      }
    })(L.root);
  }
  return { type: "scatter", mode: "lines", x: ex, y: ey,
    line: { color: cfg.branch_color || "#8a8a8a", width: +cfg.branch_width || 1.4 },
    hoverinfo: "skip", showlegend: false, meta: "edge" };
}

function tipMarker() {
  return { size: +cfg.tip_node_size || 7, color: cfg.tip_node_color || "#1f77b4",
    line: { width: 1, color: themeColors().nodeEdge } };
}

function tipTraces(L) {
  // The species-name subview already labels every tip; dropping duplicate tree
  // labels there saves width and prevents overlap in the auto-fit layout.
  var labels = !!cfg.tip_labels && !(cfg.auto_fit_subviews && state.viewMode === "trait");
  var mk = tipMarker();
  if (cfg.layout !== "circular" && cfg.layout !== "fan") {
    var x = [], y = [], txt = [], cd = [];
    L.order.forEach(function (t) { x.push(L.x[t]); y.push(L.y[t]); txt.push(tips[t - 1]); cd.push(t); });
    var tr = { type: "scatter", mode: labels ? "markers+text" : "markers",
      x: x, y: y, text: txt, customdata: cd,
      hovertemplate: "<b>%{text}</b><extra>tip</extra>",
      marker: mk, showlegend: false, meta: "tip" };
    if (labels) {
      tr.textposition = "middle right";
      tr.textfont = { size: +cfg.tip_label_size || 10, color: themeColors().font };
    }
    return [tr];
  }
  var groups = { right: [], left: [] };
  L.order.forEach(function (t) {
    var item = { x: L.x[t], y: L.y[t], text: tips[t - 1], cd: t };
    (Math.cos(L.angle[t]) >= 0 ? groups.right : groups.left).push(item);
  });
  function make(items, pos) {
    if (!items.length) return null;
    var tr = { type: "scatter", mode: labels ? "markers+text" : "markers",
      x: items.map(function (d) { return d.x; }), y: items.map(function (d) { return d.y; }),
      text: items.map(function (d) { return d.text; }),
      customdata: items.map(function (d) { return d.cd; }),
      hovertemplate: "<b>%{text}</b><extra>tip</extra>",
      marker: mk, showlegend: false, meta: "tip" };
    if (labels) {
      tr.textposition = pos;
      tr.textfont = { size: +cfg.tip_label_size || 10, color: themeColors().font };
    }
    return tr;
  }
  return [make(groups.right, "middle right"), make(groups.left, "middle left")];
}

function highlightTipTrace(L) {
  if (state.highlight.traitIdx === null && !(state.highlight.tips && state.highlight.tips.length)) return null;
  var x = [], y = [], txt = [];
  L.order.forEach(function (t) {
    if (tipMatches(t)) { x.push(L.x[t]); y.push(L.y[t]); txt.push(tips[t - 1]); }
  });
  if (!x.length) return null;
  return { type: "scatter", mode: "markers", x: x, y: y, text: txt,
    marker: { symbol: "circle-open", size: (+cfg.tip_node_size || 7) + 6,
      color: state.highlight.color, line: { width: 2.5, color: state.highlight.color } },
    hovertemplate: "<b>%{text}</b><br>highlighted<extra></extra>",
    showlegend: false, meta: "highlight" };
}

function intTrace(L) {
  if (!cfg.show_internal_nodes) return null;
  var x = [], y = [], cd = [], hover = [], labels = [], colors = [], sizes = [];
  var keys = [];
  for (var k in L.children) keys.push(+k);
  keys.sort(function (a, b) { return a - b; });
  keys.forEach(function (nd) {
    x.push(L.x[nd]); y.push(L.y[nd]); cd.push(nd);
    var cnt = leafCountFiltered(L, nd);
    hover.push(cnt + " visible tips in this clade - click to rotate");
    labels.push(cfg.show_internal_node_ids ? String(nd) : "");
    var hl = ((state.highlight.traitIdx !== null || (state.highlight.tips && state.highlight.tips.length)) && cladeMatches(nd));
    colors.push(hl ? state.highlight.color : (cfg.internal_node_color || "#e08a1e"));
    sizes.push((+cfg.internal_node_size || 9) + (hl ? 4 : 0));
  });
  return { type: "scatter", mode: cfg.show_internal_node_ids ? "markers+text" : "markers",
    x: x, y: y, customdata: cd, text: labels, hovertext: hover,
    hovertemplate: "<b>Internal node %{customdata}</b><br>%{hovertext}<extra></extra>",
    marker: { symbol: "diamond", size: sizes, color: colors,
      line: { width: 1, color: themeColors().nodeEdge } },
    textposition: (cfg.layout === "circular" || cfg.layout === "fan") ? "middle center" : "top center",
    textfont: { size: +cfg.internal_node_id_size || 9, color: themeColors().font },
    showlegend: false, meta: "internal" };
}

function heatTrace(L) {
  if (!cfg.show_heatmap) return null;
  var cols = state.visibleCols;
  if (!cols.length) return null;
  var Z = [], Y = [], CD = [];
  L.order.forEach(function (t) {
    var row = [];
    for (var j = 0; j < cols.length; j++) row.push(traitsHeat[t - 1][cols[j]]);
    Z.push(row); Y.push(L.slot[t]);
    var line = [];
    for (var j2 = 0; j2 < cols.length; j2++) line.push(tips[t - 1]);
    CD.push(line);
  });
  var tr = { type: "heatmap", z: Z,
    x: cols.map(function (i) { return tNames[i]; }), y: Y,
    xaxis: "x2", yaxis: "y2",
    colorscale: cfg.heatmap_colorscale || "Viridis",
    colorbar: { thickness: 10, len: 0.6, y: 0.5,
      title: { text: cfg.heatmap_colorbar_title || "value" } },
    hovertemplate: "<b>%{customdata}</b><br>%{x}: %{z:.3f}<extra></extra>",
    customdata: CD, showlegend: false, meta: "heatmap" };
  if (cfg.heatmap_range && cfg.heatmap_range.length === 2) {
    if (cfg.heatmap_range[0] !== null) tr.zmin = +cfg.heatmap_range[0];
    if (cfg.heatmap_range[1] !== null) tr.zmax = +cfg.heatmap_range[1];
  }
  return tr;
}

function visibleTraitValues() {
  var vals = [];
  for (var t = 1; t <= n; t++) if (state.visibleTips[t]) vals.push(traitsHeat[t - 1][state.traitIdx]);
  return vals;
}

function traitRange() {
  var v = visibleTraitValues();
  var lo = Infinity, hi = -Infinity;
  for (var i = 0; i < v.length; i++) {
    if (isFinite(v[i])) { if (v[i] < lo) lo = v[i]; if (v[i] > hi) hi = v[i]; }
  }
  if (!isFinite(lo) || !isFinite(hi)) return [0, 1];
  var pad = (hi - lo) * 0.12;
  if (!(pad > 0)) pad = 0.1;
  return [lo - pad, hi + pad];
}

function labelStride() {
  var m = state.lastLayout ? state.lastLayout.m : n;
  return Math.max(1, Math.ceil(m / 45));
}

function namesTrace(L) {
  var stride = labelStride();
  var x = [], y = [], txt = [];
  L.order.forEach(function (t, i) {
    if (i % stride !== 0) return;
    x.push(0.02); y.push(L.slot[t]); txt.push("<i>" + tips[t - 1] + "</i>");
  });
  return { type: "scatter", mode: "text", x: x, y: y, text: txt,
    xaxis: "x2", yaxis: "y2", textposition: "middle right",
    textfont: { size: 10, family: "Times New Roman, Georgia, serif", color: themeColors().font },
    hovertemplate: "<b>%{text}</b><extra>species</extra>",
    showlegend: false, meta: "names" };
}

function guideTrace(L, xr) {
  var x = [], y = [];
  L.order.forEach(function (t) {
    x.push(xr[0], xr[1], null);
    y.push(L.slot[t], L.slot[t], null);
  });
  return { type: "scatter", mode: "lines", x: x, y: y,
    xaxis: "x3", yaxis: "y2",
    line: { color: themeColors().grid, width: 0.8, dash: "dot" },
    hoverinfo: "skip", showlegend: false, meta: "guide" };
}

function valueTrace(L) {
  var x = [], y = [], txt = [], cd = [];
  L.order.forEach(function (t) {
    var v = traitsHeat[t - 1][state.traitIdx];
    x.push(v); y.push(L.slot[t]); txt.push("<i>" + tips[t - 1] + "</i>"); cd.push(v);
  });
  return { type: "scatter", mode: "markers", x: x, y: y, text: txt, customdata: cd,
    xaxis: "x3", yaxis: "y2",
    hovertemplate: "<b>%{text}</b><br>" + tNames[state.traitIdx] + ": %{customdata:.3f}<extra>trait</extra>",
    marker: { size: 8, color: "#228b22", line: { width: 1, color: themeColors().nodeEdge } },
    showlegend: false, meta: "value" };
}

function viewDimsCfg() {
  var traitView = (state.viewMode === "trait");
  var vw = numOrNull(traitView ? cfg.trait_view_width : cfg.heatmap_view_width);
  var vh = numOrNull(traitView ? cfg.trait_view_height : cfg.heatmap_view_height);
  return {
    width: (vw !== null) ? vw : numOrNull(cfg.width),
    height: (vh !== null) ? vh : numOrNull(cfg.height),
    aspect_ratio: numOrNull(cfg.aspect_ratio)
  };
}

function containerDims() {
  var d = viewDimsCfg();
  var w = d.width;
  var h = d.height;
  var ar = d.aspect_ratio;
  if (w === null) {
    var hostWidth = (wrap && wrap.clientWidth) || el.clientWidth || 900;
    var minAutoWidth = (state.viewMode === "trait") ? 700 : 640;
    w = cfg.auto_fit_subviews ? Math.max(hostWidth, minAutoWidth) : hostWidth;
  }
  if (h === null) {
    // Ignore the CSS min-height placeholder when no real height is available.
    h = (el.clientHeight && el.clientHeight > 240) ? el.clientHeight : 600;
  }
  if (ar !== null && ar > 0) {
    if (d.width !== null && d.height === null) h = w / ar;
    else if (d.height !== null && d.width === null) w = h * ar;
    else if (d.width === null && d.height === null) h = w / ar;
  }
  return { width: Math.max(320, Math.round(w)), height: Math.max(240, Math.round(h)) };
}

var tipWidthCache = {};
function naturalTipsPanelWidth(L) {
  var key = tips.join("\u0001");
  if (tipWidthCache[key] !== undefined) return tipWidthCache[key];
  var ctx = document.createElement("canvas").getContext("2d");
  ctx.font = "10px Times New Roman, Georgia, serif";
  var max = 0;
  for (var i = 0; i < L.order.length; i++) {
    var w = ctx.measureText(tips[L.order[i] - 1]).width;
    if (w > max) max = w;
  }
  // Small insets keep the label panel snug without clipping long names.
  tipWidthCache[key] = Math.max(54, Math.min(340, Math.ceil(max + 22)));
  return tipWidthCache[key];
}

function panelDomains(dims, L) {
  var radial = (cfg.layout === "circular" || cfg.layout === "fan");
  var traitView = (state.viewMode === "trait");
  var plotW = Math.max(1, dims.width - 50);
  var xMin = 0.01, xMax = 0.97;
  var available = (xMax - xMin) * plotW;
  var widths = [], gaps = [];
  var align = "left";

  if (traitView) {
    align = cfg.trait_subview_alignment || "left";
    if (cfg.auto_fit_subviews) {
      var treeW = Math.round(available * (radial ? 0.32 : 0.32));
      var tipsW = naturalTipsPanelWidth(L);
      var gap = 12;
      var valuesW = Math.round(available * 0.40);
      var deficit = treeW + tipsW + valuesW + 2 * gap - available;
      if (deficit > 0) {
        var treeShrink = Math.min(deficit, Math.max(0, treeW - 150));
        treeW -= treeShrink;
        deficit -= treeShrink;
        var valuesShrink = Math.min(deficit, Math.max(0, valuesW - 150));
        valuesW -= valuesShrink;
        deficit -= valuesShrink;
        tipsW -= Math.min(deficit, Math.max(0, tipsW - 54));
      }
      widths = [treeW, tipsW, valuesW];
      gaps = [gap, gap];
    } else {
      widths = [
        numOrNull(cfg.trait_tree_width) || Math.round(plotW * (radial ? 0.34 : 0.33)),
        numOrNull(cfg.tips_panel_width) || Math.round(plotW * (radial ? 0.23 : 0.24)),
        numOrNull(cfg.trait_values_width) || Math.round(plotW * 0.29)
      ];
      gaps = [plotW * (radial ? 0.02 : 0.02), plotW * (radial ? 0.09 : 0.08)];
    }
  } else {
    align = cfg.heatmap_subview_alignment || "left";
    if (cfg.auto_fit_subviews) {
      var treeHeatW = Math.round(available * (radial ? 0.50 : 0.46));
      var naturalHeatW = Math.max(180, Math.ceil(state.visibleCols.length * 48 + 90));
      var heatW = Math.min(Math.round(available * 0.45), naturalHeatW);
      var heatGap = 12;
      var excess = treeHeatW + heatW + heatGap - available;
      if (excess > 0) treeHeatW = Math.max(150, treeHeatW - excess);
      widths = [treeHeatW, heatW];
      gaps = [heatGap];
    } else {
      widths = [
        numOrNull(cfg.heatmap_tree_width) || Math.round(plotW * (radial ? 0.55 : 0.39)),
        numOrNull(cfg.heatmap_panel_width) || Math.round(plotW * (radial ? 0.29 : 0.40))
      ];
      gaps = [plotW * (radial ? 0.07 : 0.10)];
    }
  }

  var used = widths.reduce(function (a, b) { return a + b; }, 0) +
    gaps.reduce(function (a, b) { return a + b; }, 0);
  if (used > available) {
    var scale = Math.max(0.05, (available - gaps.reduce(function (a, b) { return a + b; }, 0)) /
      widths.reduce(function (a, b) { return a + b; }, 0));
    widths = widths.map(function (w) { return w * scale; });
    used = widths.reduce(function (a, b) { return a + b; }, 0) +
      gaps.reduce(function (a, b) { return a + b; }, 0);
  }
  var spare = Math.max(0, available - used);
  var offset = align === "center" ? spare / 2 : (align === "right" ? spare : 0);
  var cursor = xMin + offset / plotW;
  function take(width, trailingLabelWidth) {
    var start = cursor;
    var domainWidth = Math.max(1, width - (trailingLabelWidth || 0));
    cursor += width / plotW;
    return [start, start + domainWidth / plotW, 0.05, 0.92];
  }
  function advanceGap(px) { cursor += px / plotW; }

  var treeLabelGutter = (!traitView && cfg.auto_fit_subviews && cfg.tip_labels)
    ? naturalTipsPanelWidth(L) : 0;
  var tree = take(widths[0], treeLabelGutter);
  if (traitView) {
    advanceGap(gaps[0]);
    var names = take(widths[1]);
    advanceGap(gaps[1]);
    var values = take(widths[2]);
    return { tree: tree, names: names, values: values };
  }
  advanceGap(gaps[0]);
  return { tree: tree, heat: take(widths[1]) };
}

function xdom(rect) { return [rect[0], rect[1]]; }
function ydom(rect) { return [rect[2], rect[3]]; }

// Shrink a panel rectangle to the largest square (in pixels) centred inside
// it, so radial trees stay circular instead of stretching with the panel.
function squareTreeDomain(rect, dims) {
  var ml = 30, mr = 20, mt = 90, mb = 60;
  var plotW = Math.max(1, dims.width - ml - mr);
  var plotH = Math.max(1, dims.height - mt - mb);
  var panelW = (rect[1] - rect[0]) * plotW;
  var panelH = (rect[3] - rect[2]) * plotH;
  var side = Math.min(panelW, panelH);
  var cx = (rect[0] + rect[1]) / 2, cy = (rect[2] + rect[3]) / 2;
  var hw = side / plotW / 2, hh = side / plotH / 2;
  return [cx - hw, cx + hw, cy - hh, cy + hh];
}

function buildLayout(L) {
  var t = themeColors();
  var radial = (cfg.layout === "circular" || cfg.layout === "fan");
  var m = L.m;
  var dims = containerDims();
  var dom = panelDomains(dims, L);
  if (radial) dom.tree = squareTreeDomain(dom.tree, dims);
  var lay = {
    title: { text: cfg.title || "", x: 0.5, xanchor: "center", font: { size: 18, color: t.font } },
    paper_bgcolor: t.paper, plot_bgcolor: t.bg,
    font: { color: t.font, family: "-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif" },
    margin: { l: 30, r: 20, t: 90, b: 60 },
    showlegend: false, hovermode: "closest",
    autosize: !!cfg.responsive
  };
  if (cfg.subtitle) {
    lay.title.text = (cfg.title || "") + "<br><span style='font-size:12px;color:" + t.muted + "'>" + cfg.subtitle + "</span>";
  }
  if (!cfg.responsive) {
    lay.autosize = false;
    lay.width = dims.width; lay.height = dims.height;
  } else if (numOrNull(cfg.aspect_ratio) !== null) {
    lay.height = dims.height;
  }
  lay.xaxis = { domain: xdom(dom.tree), showgrid: false, zeroline: false, showticklabels: false, visible: false, fixedrange: true };
  lay.yaxis = { domain: ydom(dom.tree), showgrid: false, zeroline: false, showticklabels: false, visible: false, fixedrange: true };
  if (radial) {
    lay.xaxis.range = [-1.18, 1.18]; lay.yaxis.range = [-1.18, 1.18];
  } else {
    lay.yaxis.range = [0.5, m + 0.5];
  }
  lay.annotations = [];
  if (state.viewMode === "trait") {
    lay.xaxis2 = { domain: xdom(dom.names), anchor: "y2", range: [0, 1], type: "linear", showgrid: false, zeroline: false, showticklabels: false, visible: false, fixedrange: true };
    lay.xaxis3 = { domain: xdom(dom.values), anchor: "y2", type: "linear", showgrid: true, gridcolor: t.grid, zeroline: false, automargin: true,
      title: { text: cfg.trait_panel_title || tNames[state.traitIdx] || "Trait value", font: { size: 12, color: t.font } } };
    lay.yaxis2 = { domain: ydom(dom.names), anchor: "x2", range: [0.5, m + 0.5], visible: false, showgrid: false, zeroline: false, showticklabels: false, fixedrange: true };
    lay.annotations.push({ text: cfg.tree_panel_title || "Phylogenetic tree", xref: "paper", yref: "paper",
      x: (dom.tree[0] + dom.tree[1]) / 2, y: 1.01, showarrow: false, font: { size: 13, color: t.font } });
    lay.annotations.push({ text: (cfg.species_panel_title || "Species") + (labelStride() > 1 ? " (1 in " + labelStride() + ")" : ""),
      xref: "paper", yref: "paper", x: (dom.names[0] + dom.names[1]) / 2, y: 1.01, showarrow: false, font: { size: 13, color: t.font } });
    lay.annotations.push({ text: cfg.trait_panel_title || (tNames[state.traitIdx] + " value"),
      xref: "paper", yref: "paper", x: (dom.values[0] + dom.values[1]) / 2, y: 1.01, showarrow: false, font: { size: 13, color: t.font } });
  } else {
    lay.xaxis2 = { domain: xdom(dom.heat), anchor: "y2", type: "category", showgrid: false, zeroline: false, showticklabels: true, visible: true,
      title: { text: cfg.heatmap_title || "Trait", font: { size: 12, color: t.font } }, tickangle: -30, automargin: true, tickfont: { color: t.font } };
    lay.yaxis2 = { domain: ydom(dom.heat), anchor: "x2", range: [0.5, m + 0.5], visible: false, showgrid: false, zeroline: false, showticklabels: false, fixedrange: true };
    lay.xaxis3 = { visible: false };
    if (cfg.tree_panel_title) {
      lay.annotations.push({ text: cfg.tree_panel_title, xref: "paper", yref: "paper",
        x: (dom.tree[0] + dom.tree[1]) / 2, y: 1.01, showarrow: false, font: { size: 13, color: t.font } });
    }
    if (cfg.heatmap_panel_title) {
      lay.annotations.push({ text: cfg.heatmap_panel_title, xref: "paper", yref: "paper",
        x: (dom.heat[0] + dom.heat[1]) / 2, y: 1.01, showarrow: false, font: { size: 13, color: t.font } });
    }
    if (!cfg.show_heatmap || !state.visibleCols.length) {
      lay.annotations.push({ text: "Heatmap hidden or no trait columns selected.", xref: "paper", yref: "paper",
        x: (dom.heat[0] + dom.heat[1]) / 2, y: 0.5, showarrow: false, font: { size: 12, color: t.muted } });
    }
  }
  return lay;
}

function plotConfig() {
  var dims = containerDims();
  return { responsive: !!cfg.responsive, displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d"],
    toImageButtonOptions: { format: "png", filename: "tree_viz", width: dims.width, height: dims.height, scale: 1 } };
}

function applyContainerSizing() {
  var dims = containerDims();
  var d = viewDimsCfg();
  var arOnly = d.aspect_ratio !== null &&
    d.width === null && d.height === null;
  if (cfg.responsive && arOnly) {
    var hostWidth0 = wrap && wrap.clientWidth;
    el.style.width = (hostWidth0 && dims.width > hostWidth0) ? dims.width + "px" : "100%";
    el.style.height = dims.height + "px";
    el.style.flex = "0 0 auto";
    if (wrap) wrap.style.height = "auto";
  } else if (cfg.responsive) {
    // Fill a parent with a definite height (standalone page or a sized knitr
    // container). In R Markdown the container is usually auto-height, so give
    // the plot an explicit pixel height instead of collapsing to min-height.
    var parent = wrap && wrap.parentNode;
    var parentH = parent ? parent.clientHeight : 0;
    var parentSized = !!(parent && parent.style && parent.style.height &&
                         parent.style.height !== "auto");
    var wc = document.getElementById("htmlwidget_container");
    var standalone = !!(wc && wc.parentNode === document.body);
    var heightIsDerived = d.aspect_ratio !== null && d.width !== null;
    var fill = (standalone || parentSized) && parentH >= 240 &&
               d.height === null && !heightIsDerived;
    var hostWidth = wrap && wrap.clientWidth;
    var needsScroll = !!(hostWidth && dims.width > hostWidth);
    if (wrap) wrap.style.height = fill ? "100%" : "auto";
    if (wrap) wrap.style.overflowX = needsScroll ? "auto" : "";
    el.style.width = (d.width !== null || needsScroll) ? (dims.width + "px") : "100%";
    if (fill) {
      el.style.height = "";
      el.style.flex = "1 1 auto";
    } else {
      el.style.height = dims.height + "px";
      el.style.flex = "0 0 auto";
    }
  } else {
    if (wrap) wrap.style.height = "auto";
    if (wrap) wrap.style.overflowX = (wrap.clientWidth && dims.width > wrap.clientWidth) ? "auto" : "";
    el.style.width = dims.width + "px";
    el.style.height = dims.height + "px";
    el.style.flex = "0 0 auto";
  }
  el.style.minHeight = "240px";
}

function updateCost(L) {
  var s = 0;
  for (var i = 1; i < L.order.length; i++) {
    var a = traitsCost[L.order[i - 1] - 1];
    var b = traitsCost[L.order[i] - 1];
    var d2 = 0;
    for (var j = 0; j < a.length; j++) { var dd = a[j] - b[j]; d2 += dd * dd; }
    s += Math.sqrt(d2);
  }
  var e = q("cost");
  if (e) e.textContent = s.toFixed(4);
}

function updateStatus() {
  var vis = visibleTipCount();
  var msg = "Showing " + vis + " of " + n + " tips";
  if (state.visibleCols.length !== tNames.length) msg += "; " + state.visibleCols.length + " of " + tNames.length + " trait columns";
  var filters = [];
  if (state.filter.traitIdx !== null) {
    filters.push(tNames[state.filter.traitIdx] + " in [" +
      (state.filter.min === null ? "-Inf" : state.filter.min) + ", " +
      (state.filter.max === null ? "Inf" : state.filter.max) + "]");
  }
  if (state.filter.tips && state.filter.tips.length) filters.push(state.filter.tips.length + " named tips");
  if (state.clades.length) filters.push(state.clades.length + " clade(s)");
  if (filters.length) msg += " | filters: " + filters.join(", ");
  if (state.highlight.traitIdx !== null) msg += " | highlight: " + tNames[state.highlight.traitIdx];
  var e = q("status");
  if (e) e.textContent = msg;
}

function q(role) { return ui ? ui.querySelector('[data-role="' + role + '"]') : null; }
function qa(role) { return ui ? ui.querySelectorAll('[data-role="' + role + '"]') : []; }

function revForTarget(targetLabels) {
  var tIdx = {};
  for (var i = 0; i < targetLabels.length; i++) tIdx[targetLabels[i]] = i;
  var base = [];
  (function dfs(node) {
    var ch = children[node];
    if (!ch) { base.push(node * 1); return; }
    for (var j = 0; j < ch.length; j++) dfs(ch[j]);
  })(root);
  var pos = {};
  base.forEach(function (t) { pos[t] = tIdx[tips[t - 1]]; });
  var mn = {};
  (function calc(node) {
    var ch = children[node];
    if (!ch) { mn[node] = pos[node]; return; }
    var lo = Infinity;
    for (var j = 0; j < ch.length; j++) { calc(ch[j]); if (mn[ch[j]] < lo) lo = mn[ch[j]]; }
    mn[node] = lo;
  })(root);
  var r = {};
  (function orient(node) {
    var ch = children[node];
    if (!ch) return;
    if (ch.length === 2) r[node] = mn[ch[1]] < mn[ch[0]];
    for (var j = 0; j < ch.length; j++) orient(ch[j]);
  })(root);
  return r;
}

function redraw() {
  applyContainerSizing();
  var F = buildFiltered(state.rev);
  if (!F || leafOrderFiltered(F).length < 1) {
    state.lastLayout = null;
    var lay0 = buildLayout({ m: 1 });
    lay0.annotations.push({ text: "No tips match the current filters.", xref: "paper", yref: "paper",
      x: 0.5, y: 0.5, showarrow: false, font: { size: 16, color: themeColors().muted } });
    Plotly.react(el, [], lay0, plotConfig());
    var e0 = q("cost"); if (e0) e0.textContent = "0.0000";
    updateStatus();
    return;
  }
  var L = layoutCalcFiltered(F);
  if (cfg.layout === "circular" || cfg.layout === "fan") L = circularLayout(L);
  state.lastLayout = L;
  var traces = [];
  if (Object.keys(L.children).length) traces.push(edgeTrace(L));
  var tt = tipTraces(L);
  for (var i = 0; i < tt.length; i++) if (tt[i]) traces.push(tt[i]);
  var ht = highlightTipTrace(L); if (ht) traces.push(ht);
  var it = intTrace(L); if (it) traces.push(it);
  if (state.viewMode === "trait") {
    traces.push(namesTrace(L));
    traces.push(guideTrace(L, traitRange()));
    traces.push(valueTrace(L));
  } else {
    var hm = heatTrace(L); if (hm) traces.push(hm);
  }
  Plotly.react(el, traces, buildLayout(L), plotConfig());
  updateCost(L);
  updateStatus();
}

function syncThemeClass() {
  if (!ui) return;
  if (cfg.theme === "dark") ui.classList.add("itv-dark");
  else ui.classList.remove("itv-dark");
  var btn = q("theme");
  if (btn) btn.textContent = (cfg.theme === "dark") ? "Light mode" : "Dark mode";
  var t = themeColors();
  if (ui.parentNode) ui.parentNode.style.backgroundColor = t.paper;
}

function syncViewButton() {
  var b = q("view");
  if (b) b.textContent = (state.viewMode === "trait") ? "Back to heatmap view" : "Species + trait panel";
}

function markModified() {
  state.presetValue = "";
  var sel = q("preset_select");
  if (sel) sel.value = "";
}

function setVal(role, val) {
  var e = q(role);
  if (!e) return;
  if (e.type === "checkbox") e.checked = !!val;
  else if (val === null || val === undefined) e.value = "";
  else e.value = val;
}

function syncUI() {
  setVal("title", cfg.title || "");
  setVal("subtitle", cfg.subtitle || "");
  setVal("layout", cfg.layout || "rectangular");
  setVal("tip_labels", cfg.tip_labels !== false);
  setVal("tip_label_size", cfg.tip_label_size);
  setVal("show_internal_nodes", cfg.show_internal_nodes !== false);
  setVal("internal_node_size", cfg.internal_node_size);
  setVal("show_internal_node_ids", !!cfg.show_internal_node_ids);
  setVal("internal_node_id_size", cfg.internal_node_id_size);
  setVal("tip_node_color", cfg.tip_node_color || "#1f77b4");
  setVal("tip_node_size", cfg.tip_node_size);
  setVal("internal_node_color", cfg.internal_node_color || "#e08a1e");
  setVal("branch_color", cfg.branch_color || "#8a8a8a");
  setVal("branch_width", cfg.branch_width);
  setVal("heatmap_colorscale", cfg.heatmap_colorscale || "Viridis");
  setVal("heatmap_title", cfg.heatmap_title || "Trait");
  setVal("heatmap_colorbar_title", cfg.heatmap_colorbar_title || "value");
  setVal("heatmap_panel_title", cfg.heatmap_panel_title || "");
  setVal("tree_panel_title", cfg.tree_panel_title || "");
  setVal("species_panel_title", cfg.species_panel_title || "");
  setVal("trait_panel_title", cfg.trait_panel_title || "");
  setVal("fan_angle", cfg.fan_angle);
  setVal("show_heatmap", cfg.show_heatmap !== false);
  setVal("responsive", cfg.responsive !== false);
  setVal("width", cfg.width);
  setVal("height", cfg.height);
  setVal("aspect_ratio", cfg.aspect_ratio);
  setVal("heatmap_view_width", cfg.heatmap_view_width);
  setVal("heatmap_view_height", cfg.heatmap_view_height);
  setVal("trait_view_width", cfg.trait_view_width);
  setVal("trait_view_height", cfg.trait_view_height);
  setVal("auto_fit_subviews", !!cfg.auto_fit_subviews);
  setVal("heatmap_tree_width", cfg.heatmap_tree_width);
  setVal("heatmap_panel_width", cfg.heatmap_panel_width);
  setVal("trait_tree_width", cfg.trait_tree_width);
  setVal("tips_panel_width", cfg.tips_panel_width);
  setVal("trait_values_width", cfg.trait_values_width);
  setVal("heatmap_subview_alignment", cfg.heatmap_subview_alignment || "left");
  setVal("trait_subview_alignment", cfg.trait_subview_alignment || "left");
  ["heatmap_tree_width", "heatmap_panel_width", "trait_tree_width",
    "tips_panel_width", "trait_values_width"].forEach(function (role) {
      var input = q(role);
      if (input) input.disabled = !!cfg.auto_fit_subviews;
    });

  var auto = !(cfg.heatmap_range && cfg.heatmap_range.length === 2);
  setVal("heatmap_auto", auto);
  if (!auto) { setVal("heatmap_min", cfg.heatmap_range[0]); setVal("heatmap_max", cfg.heatmap_range[1]); }
  else { setVal("heatmap_min", ""); setVal("heatmap_max", ""); }
  var hmin = q("heatmap_min"), hmax = q("heatmap_max");
  if (hmin) hmin.disabled = auto;
  if (hmax) hmax.disabled = auto;

  setVal("filter_trait", state.filter.traitIdx === null ? -1 : state.filter.traitIdx);
  setVal("filter_min", state.filter.min);
  setVal("filter_max", state.filter.max);
  setVal("filter_tips", state.filter.tips ? state.filter.tips.join(", ") : "");
  var ft = q("filter_trait"), fmin = q("filter_min"), fmax = q("filter_max");
  if (fmin) fmin.disabled = (state.filter.traitIdx === null);
  if (fmax) fmax.disabled = (state.filter.traitIdx === null);
  var cl = q("clade");
  if (cl) for (var i = 0; i < cl.options.length; i++) cl.options[i].selected = state.clades.indexOf(+cl.options[i].value) >= 0;

  setVal("highlight_trait", state.highlight.traitIdx === null ? -1 : state.highlight.traitIdx);
  setVal("highlight_min", state.highlight.min);
  setVal("highlight_max", state.highlight.max);
  setVal("highlight_color", state.highlight.color || "#e41a1c");
  var hlt = q("highlight_trait"), hlmin = q("highlight_min"), hlmax = q("highlight_max");
  if (hlmin) hlmin.disabled = (state.highlight.traitIdx === null);
  if (hlmax) hlmax.disabled = (state.highlight.traitIdx === null);

  var boxes = qa("col_box");
  for (var b = 0; b < boxes.length; b++) boxes[b].checked = state.visibleCols.indexOf(+boxes[b].value) >= 0;

  syncThemeClass();
  syncViewButton();
  var sel = q("preset_select");
  if (sel) sel.value = state.presetValue || "";
  updateStatus();
}

function refreshCladeOptions() {
  var sel = q("clade");
  if (!sel) return;
  var keys = [];
  for (var k in children) keys.push(+k);
  keys.sort(function (a, b) { return a - b; });
  var h = "";
  for (var i = 0; i < keys.length; i++) {
    h += '<option value="' + keys[i] + '">Node ' + keys[i] + " (" + leafCountOriginal(keys[i]) + " tips)</option>";
  }
  sel.innerHTML = h;
}

function refreshColumnFilter() {
  var box = q("col_filter");
  if (!box) return;
  box.innerHTML = "";
  for (var i = 0; i < tNames.length; i++) {
    var lab = document.createElement("label");
    lab.className = "itv-check";
    var cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = i;
    cb.setAttribute("data-role", "col_box");
    cb.checked = state.visibleCols.indexOf(i) >= 0;
    cb.addEventListener("change", function () {
      var idx = +this.value;
      var pos = state.visibleCols.indexOf(idx);
      if (this.checked && pos < 0) state.visibleCols.push(idx);
      if (!this.checked && pos >= 0) state.visibleCols.splice(pos, 1);
      state.visibleCols.sort(function (a, b) { return a - b; });
      markModified(); redraw();
    });
    lab.appendChild(cb);
    lab.appendChild(document.createTextNode(" " + tNames[i]));
    box.appendChild(lab);
  }
}

var LS_KEY = "itv_presets_v1";
var memoryStore = {};
function readStore() {
  try {
    var raw = window.localStorage.getItem(LS_KEY);
    return raw ? JSON.parse(raw) : {};
  } catch (e) { return memoryStore; }
}
function writeStore(obj) {
  try { window.localStorage.setItem(LS_KEY, JSON.stringify(obj)); }
  catch (e) { memoryStore = obj; }
}

function refreshPresetOptions() {
  var sel = q("preset_select");
  if (!sel) return;
  var h = '<option value="">Custom / modified</option>';
  var builtins = TD.presets || {};
  for (var name in builtins) {
    h += '<option value="builtin:' + esc(name) + '">' + esc(name) + " (built-in)</option>";
  }
  var saved = readStore();
  for (var sname in saved) {
    h += '<option value="saved:' + esc(sname) + '">' + esc(sname) + " (saved)</option>";
  }
  sel.innerHTML = h;
  sel.value = state.presetValue || "";
}

function downloadBlob(content, filename, mime) {
  var blob = new Blob([content], { type: mime });
  var url = URL.createObjectURL(blob);
  var a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click();
  setTimeout(function () { URL.revokeObjectURL(url); if (a.parentNode) a.parentNode.removeChild(a); }, 1200);
}

function nwkLabel(s) {
  s = String(s);
  return /[(),:;'\s]/.test(s) ? "'" + s.replace(/'/g, "''") + "'" : s;
}
function newickNode(node, revMap) {
  var ch = children[node];
  var len = edgeLen[node] ? (":" + edgeLen[node]) : "";
  if (!ch) return nwkLabel(tips[node - 1]) + len;
  var ordered = revMap[node] ? ch.slice().reverse() : ch;
  var kids = ordered.map(function (c) { return newickNode(c, revMap); }).join(",");
  return "(" + kids + ")" + (node === root ? "" : len);
}
function newickString() {
  return newickNode(root, state.rev) + ";";
}

function snapshotConfig() {
  var c = {};
  for (var k in cfg) c[k] = cfg[k];
  c.tip_filter = null;
  if (state.filter.traitIdx !== null || (state.filter.tips && state.filter.tips.length)) {
    c.tip_filter = {};
    if (state.filter.traitIdx !== null) {
      c.tip_filter.trait = tNames[state.filter.traitIdx];
      c.tip_filter.min = state.filter.min;
      c.tip_filter.max = state.filter.max;
    }
    if (state.filter.tips && state.filter.tips.length) c.tip_filter.tips = state.filter.tips.slice();
  }
  c.clade = state.clades.length ? state.clades.slice() : null;
  c.heatmap_columns = state.visibleCols.map(function (i) { return tNames[i]; });
  c.highlight = null;
  if (state.highlight.traitIdx !== null || (state.highlight.tips && state.highlight.tips.length)) {
    c.highlight = { color: state.highlight.color };
    if (state.highlight.traitIdx !== null) {
      c.highlight.trait = tNames[state.highlight.traitIdx];
      c.highlight.min = state.highlight.min;
      c.highlight.max = state.highlight.max;
    }
    if (state.highlight.tips && state.highlight.tips.length) c.highlight.tips = state.highlight.tips.slice();
  }
  return c;
}

function applyConfig(newCfg, resetRotations) {
  cfg = mergeCfg(DEFAULT_CFG, newCfg || {});
  if (resetRotations) initRev();
  initStateFromConfig();
  syncUI();
  redraw();
}

function exportImage(format) {
  var dims = containerDims();
  Plotly.downloadImage(el, { format: format, width: dims.width, height: dims.height, filename: "tree_viz" });
  var e = q("export_status"); if (e) e.textContent = format.toUpperCase() + " export started.";
}

function exportPDF() {
  var dims = containerDims();
  Plotly.toImage(el, { format: "svg", width: dims.width, height: dims.height }).then(function (url) {
    var w = window.open("", "_blank");
    if (!w) { var e = q("export_status"); if (e) e.textContent = "Allow pop-ups to print / save as PDF."; return; }
    w.document.write('<html><head><title>tree_viz</title></head><body style="margin:0">' +
      '<img src="' + url + '" style="width:100%;height:auto">' +
      '<script>window.onload=function(){setTimeout(function(){window.print();},300);};<\/script>' +
      '</body></html>');
    w.document.close();
    var e2 = q("export_status"); if (e2) e2.textContent = "Print dialog opened; choose Save as PDF.";
  }).catch(function () {
    var e3 = q("export_status"); if (e3) e3.textContent = "PDF export failed in this browser.";
  });
}

function exportCSV() {
  var L = state.lastLayout;
  if (!L) { var e = q("export_status"); if (e) e.textContent = "Nothing to export."; return; }
  var cols = state.visibleCols;
  var lines = ["tip," + cols.map(function (i) { return tNames[i]; }).join(",")];
  L.order.forEach(function (t) {
    var row = [tips[t - 1]];
    for (var j = 0; j < cols.length; j++) row.push(traitsHeat[t - 1][cols[j]]);
    lines.push(row.join(","));
  });
  downloadBlob(lines.join("\n"), "tree_heatmap.csv", "text/csv");
  var e2 = q("export_status"); if (e2) e2.textContent = "CSV exported.";
}

function exportNewick() {
  downloadBlob(newickString(), "tree_viz.tre", "text/plain");
  var e = q("export_status"); if (e) e.textContent = "Newick tree exported.";
}

function exportConfig() {
  downloadBlob(JSON.stringify(snapshotConfig(), null, 2), "tree_viz_config.json", "application/json");
  var e = q("export_status"); if (e) e.textContent = "Configuration JSON exported.";
}

function importConfigFile(file) {
  var reader = new FileReader();
  reader.onload = function (ev) {
    try {
      var obj = JSON.parse(ev.target.result);
      applyConfig(obj, false);
      var e = q("export_status"); if (e) e.textContent = "Configuration loaded.";
    } catch (err) {
      var e2 = q("export_status"); if (e2) e2.textContent = "Invalid JSON: " + err.message;
    }
  };
  reader.readAsText(file);
}

function applyPreset() {
  var sel = q("preset_select");
  if (!sel || !sel.value) return;
  var v = sel.value;
  if (v.indexOf("builtin:") === 0) {
    var name = v.slice(8);
    if (TD.presets && TD.presets[name]) {
      applyConfig(TD.presets[name], false);
      state.presetValue = v; syncUI();
      var e = q("preset_status"); if (e) e.textContent = "Applied built-in preset: " + name;
    }
  } else if (v.indexOf("saved:") === 0) {
    var sname = v.slice(6);
    var saved = readStore();
    if (saved[sname]) {
      applyConfig(saved[sname], false);
      state.presetValue = v; syncUI();
      var e2 = q("preset_status"); if (e2) e2.textContent = "Applied saved preset: " + sname;
    }
  }
}

function savePreset() {
  var nameEl = q("preset_name");
  var name = nameEl ? nameEl.value.trim() : "";
  if (!name) name = "My preset";
  if (TD.presets && TD.presets[name]) name = name + " (custom)";
  var saved = readStore();
  saved[name] = snapshotConfig();
  writeStore(saved);
  state.presetValue = "saved:" + name;
  refreshPresetOptions();
  syncUI();
  var e = q("preset_status"); if (e) e.textContent = "Saved preset: " + name;
}

function deletePreset() {
  var sel = q("preset_select");
  if (!sel || sel.value.indexOf("saved:") !== 0) {
    var e = q("preset_status"); if (e) e.textContent = "Select a saved preset to delete.";
    return;
  }
  var name = sel.value.slice(6);
  var saved = readStore();
  delete saved[name];
  writeStore(saved);
  state.presetValue = "";
  refreshPresetOptions();
  syncUI();
  var e2 = q("preset_status"); if (e2) e2.textContent = "Deleted saved preset: " + name;
}

function restoreDefaults() {
  cfg = mergeCfg(DEFAULT_CFG, {});
  state.presetValue = "builtin:default";
  initRev();
  initStateFromConfig();
  refreshPresetOptions();
  syncUI();
  redraw();
  var e = q("preset_status"); if (e) e.textContent = "Restored default settings.";
}

function bindEvents() {
  function on(role, evt, fn) {
    var e = q(role);
    if (e) e.addEventListener(evt, fn);
  }
  on("reset", "click", function () { initRev(); redraw(); });
  if (oloOrder) on("olo", "click", function () { state.rev = revForTarget(oloOrder); redraw(); });
  on("view", "click", function () {
    state.viewMode = (state.viewMode === "trait") ? "heat" : "trait";
    syncViewButton(); redraw();
  });
  on("theme", "click", function () {
    cfg.theme = (cfg.theme === "dark") ? "light" : "dark";
    markModified(); syncThemeClass(); redraw();
  });

  on("title", "change", function () { cfg.title = this.value; markModified(); redraw(); });
  on("subtitle", "change", function () { cfg.subtitle = this.value; markModified(); redraw(); });
  on("layout", "change", function () { cfg.layout = this.value; markModified(); redraw(); });
  on("tip_labels", "change", function () { cfg.tip_labels = this.checked; markModified(); redraw(); });
  on("tip_label_size", "change", function () { cfg.tip_label_size = numOrNull(this.value); markModified(); redraw(); });
  on("show_internal_nodes", "change", function () { cfg.show_internal_nodes = this.checked; markModified(); redraw(); });
  on("internal_node_size", "change", function () { cfg.internal_node_size = numOrNull(this.value); markModified(); redraw(); });
  on("show_internal_node_ids", "change", function () { cfg.show_internal_node_ids = this.checked; markModified(); redraw(); });
  on("internal_node_id_size", "change", function () { cfg.internal_node_id_size = numOrNull(this.value); markModified(); redraw(); });
  on("tip_node_color", "change", function () { cfg.tip_node_color = this.value; markModified(); redraw(); });
  on("tip_node_size", "change", function () { cfg.tip_node_size = numOrNull(this.value); markModified(); redraw(); });
  on("internal_node_color", "change", function () { cfg.internal_node_color = this.value; markModified(); redraw(); });
  on("branch_color", "change", function () { cfg.branch_color = this.value; markModified(); redraw(); });
  on("branch_width", "change", function () { cfg.branch_width = numOrNull(this.value); markModified(); redraw(); });
  on("heatmap_colorscale", "change", function () { cfg.heatmap_colorscale = this.value; markModified(); redraw(); });
  on("heatmap_title", "change", function () { cfg.heatmap_title = this.value; markModified(); redraw(); });
  on("heatmap_colorbar_title", "change", function () { cfg.heatmap_colorbar_title = this.value; markModified(); redraw(); });
  on("heatmap_panel_title", "change", function () { cfg.heatmap_panel_title = this.value; markModified(); redraw(); });
  on("tree_panel_title", "change", function () { cfg.tree_panel_title = this.value; markModified(); redraw(); });
  on("species_panel_title", "change", function () { cfg.species_panel_title = this.value; markModified(); redraw(); });
  on("trait_panel_title", "change", function () { cfg.trait_panel_title = this.value; markModified(); redraw(); });
  on("fan_angle", "change", function () { cfg.fan_angle = numOrNull(this.value); markModified(); redraw(); });
  on("show_heatmap", "change", function () { cfg.show_heatmap = this.checked; markModified(); redraw(); });
  on("responsive", "change", function () { cfg.responsive = this.checked; markModified(); applyContainerSizing(); redraw(); });
  on("width", "change", function () { cfg.width = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("height", "change", function () { cfg.height = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("aspect_ratio", "change", function () { cfg.aspect_ratio = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("heatmap_view_width", "change", function () { cfg.heatmap_view_width = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("heatmap_view_height", "change", function () { cfg.heatmap_view_height = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("trait_view_width", "change", function () { cfg.trait_view_width = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("trait_view_height", "change", function () { cfg.trait_view_height = numOrNull(this.value); markModified(); applyContainerSizing(); redraw(); });
  on("auto_fit_subviews", "change", function () {
    cfg.auto_fit_subviews = this.checked;
    markModified(); syncUI(); redraw();
  });
  on("heatmap_tree_width", "change", function () { cfg.heatmap_tree_width = numOrNull(this.value); markModified(); redraw(); });
  on("heatmap_panel_width", "change", function () { cfg.heatmap_panel_width = numOrNull(this.value); markModified(); redraw(); });
  on("trait_tree_width", "change", function () { cfg.trait_tree_width = numOrNull(this.value); markModified(); redraw(); });
  on("tips_panel_width", "change", function () { cfg.tips_panel_width = numOrNull(this.value); markModified(); redraw(); });
  on("trait_values_width", "change", function () { cfg.trait_values_width = numOrNull(this.value); markModified(); redraw(); });
  on("heatmap_subview_alignment", "change", function () { cfg.heatmap_subview_alignment = this.value; markModified(); redraw(); });
  on("trait_subview_alignment", "change", function () { cfg.trait_subview_alignment = this.value; markModified(); redraw(); });
  on("heatmap_auto", "change", function () {
    if (this.checked) cfg.heatmap_range = null;
    else {
      var lo = numOrNull(q("heatmap_min") ? q("heatmap_min").value : null);
      var hi = numOrNull(q("heatmap_max") ? q("heatmap_max").value : null);
      cfg.heatmap_range = (lo !== null && hi !== null) ? [lo, hi] : null;
    }
    markModified(); syncUI(); redraw();
  });
  on("heatmap_min", "change", function () {
    if (q("heatmap_auto") && q("heatmap_auto").checked) return;
    var lo = numOrNull(this.value), hi = numOrNull(q("heatmap_max") ? q("heatmap_max").value : null);
    cfg.heatmap_range = (lo !== null && hi !== null) ? [lo, hi] : null;
    markModified(); redraw();
  });
  on("heatmap_max", "change", function () {
    if (q("heatmap_auto") && q("heatmap_auto").checked) return;
    var lo = numOrNull(q("heatmap_min") ? q("heatmap_min").value : null), hi = numOrNull(this.value);
    cfg.heatmap_range = (lo !== null && hi !== null) ? [lo, hi] : null;
    markModified(); redraw();
  });

  on("filter_trait", "change", function () {
    var v = +this.value;
    state.filter.traitIdx = (v >= 0) ? v : null;
    markModified(); syncUI(); recomputeVisibility(); redraw();
  });
  on("filter_min", "change", function () {
    state.filter.min = numOrNull(this.value); markModified(); recomputeVisibility(); redraw();
  });
  on("filter_max", "change", function () {
    state.filter.max = numOrNull(this.value); markModified(); recomputeVisibility(); redraw();
  });
  on("filter_tips", "change", function () {
    var arr = this.value.split(",").map(function (s) { return s.trim(); }).filter(function (s) { return s.length > 0; });
    state.filter.tips = arr.length ? arr : null;
    markModified(); recomputeVisibility(); redraw();
  });
  on("clade", "change", function () {
    var sel = this, ids = [];
    for (var i = 0; i < sel.options.length; i++) if (sel.options[i].selected) ids.push(+sel.options[i].value);
    state.clades = ids;
    markModified(); recomputeVisibility(); redraw();
  });
  on("clear_filters", "click", function () {
    state.filter = { traitIdx: null, min: null, max: null, tips: null };
    state.clades = [];
    markModified(); syncUI(); recomputeVisibility(); redraw();
  });
  on("highlight_trait", "change", function () {
    var v = +this.value;
    state.highlight.traitIdx = (v >= 0) ? v : null;
    markModified(); syncUI(); matchCache = {}; redraw();
  });
  on("highlight_min", "change", function () {
    state.highlight.min = numOrNull(this.value); markModified(); matchCache = {}; redraw();
  });
  on("highlight_max", "change", function () {
    state.highlight.max = numOrNull(this.value); markModified(); matchCache = {}; redraw();
  });
  on("highlight_color", "change", function () {
    state.highlight.color = this.value; markModified(); redraw();
  });
  on("clear_highlight", "click", function () {
    state.highlight = { traitIdx: null, min: null, max: null, tips: null, color: "#e41a1c" };
    markModified(); syncUI(); matchCache = {}; redraw();
  });

  on("export_png", "click", function () { exportImage("png"); });
  on("export_svg", "click", function () { exportImage("svg"); });
  on("export_pdf", "click", exportPDF);
  on("export_newick", "click", exportNewick);
  on("export_csv", "click", exportCSV);
  on("export_config", "click", exportConfig);
  on("import_config", "change", function () { if (this.files && this.files[0]) importConfigFile(this.files[0]); });
  on("preset_apply", "click", applyPreset);
  on("preset_save", "click", savePreset);
  on("preset_delete", "click", deletePreset);
  on("preset_restore", "click", restoreDefaults);
}

// ---------------------------------------------------------------------------
// Build the control panel and start the widget
// ---------------------------------------------------------------------------
initRev();
initStateFromConfig();

var styleEl = document.createElement("style");
styleEl.textContent = [
  ".itv-wrap{display:flex;flex-direction:column;width:100%;height:100%;min-height:0;box-sizing:border-box;overflow-x:auto;}",
  ".itv-ui{--itv-bg:#f7f7f7;--itv-fg:#222;--itv-border:#ddd;--itv-input:#fff;--itv-accent:#1f77b4;--itv-muted:#666;",
  "font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:13px;color:var(--itv-fg);",
  "background:var(--itv-bg);border:1px solid var(--itv-border);border-radius:6px;padding:6px 8px;margin:0 0 6px 0;box-sizing:border-box;flex:0 0 auto;}",
  ".itv-ui.itv-dark{--itv-bg:#171b21;--itv-fg:#e8e8e8;--itv-border:#2f363f;--itv-input:#0f1115;--itv-accent:#4da3ff;--itv-muted:#a8b0b8;}",
  ".itv-ui button{font:inherit;padding:2px 8px;margin:0 2px 0 0;border:1px solid var(--itv-border);border-radius:4px;background:var(--itv-input);color:var(--itv-fg);cursor:pointer;}",
  ".itv-ui button:hover{border-color:var(--itv-accent);}",
  ".itv-ui input,.itv-ui select{font:inherit;color:var(--itv-fg);background:var(--itv-input);border:1px solid var(--itv-border);border-radius:4px;padding:1px 4px;box-sizing:border-box;}",
  ".itv-ui input[type=text],.itv-ui input[type=number]{width:100%;min-width:90px;}",
  ".itv-ui input[type=color]{width:42px;height:24px;padding:0;}",
  ".itv-ui details{margin:4px 0;}",
  ".itv-ui summary{cursor:pointer;font-weight:600;padding:2px 0;}",
  ".itv-ui .itv-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:5px 12px;padding:4px 0;}",
  ".itv-ui .itv-grid label{display:flex;align-items:center;gap:5px;min-width:0;}",
  ".itv-ui .itv-grid label input[type=text],.itv-ui .itv-grid label input[type=number]{flex:1 1 auto;}",
  ".itv-ui .itv-wide{grid-column:1/-1;}",
  ".itv-ui .itv-toolbar{display:flex;flex-wrap:wrap;align-items:center;gap:4px 6px;}",
  ".itv-ui .itv-status{color:var(--itv-muted);font-size:12px;margin-left:4px;}",
  ".itv-ui .itv-hint{color:var(--itv-muted);font-size:12px;margin-top:4px;}",
  ".itv-ui .itv-checks{display:flex;flex-wrap:wrap;gap:3px 10px;}",
  ".itv-ui .itv-check{display:inline-flex;align-items:center;gap:3px;font-size:12px;}",
  ".itv-ui .itv-subhead{font-weight:600;margin:2px 0;}",
  ".itv-ui .itv-actions{display:flex;flex-wrap:wrap;align-items:center;gap:4px 6px;}",
  ".itv-ui .itv-file{display:inline-flex;align-items:center;gap:4px;font-size:12px;}"
].join("");

var h = [];
h.push('<div class="itv-toolbar">');
h.push('<span><b>Adjacent trait distance' + (costScaled ? " (scaled)" : "") + ':</b> <span data-role="cost">0</span></span>');
h.push('<button type="button" data-role="reset">Reset rotations</button>');
if (oloOrder) h.push('<button type="button" data-role="olo">Apply OLO-optimal order</button>');
h.push('<button type="button" data-role="view">Species + trait panel</button>');
h.push('<button type="button" data-role="theme">Dark mode</button>');
h.push('<span class="itv-status" data-role="status"></span>');
h.push('</div>');

var traitOpts = "";
for (var ti = 0; ti < tNames.length; ti++) traitOpts += '<option value="' + ti + '">' + esc(tNames[ti]) + "</option>";
var scaleNames = ["Viridis", "Plasma", "Inferno", "Magma", "Cividis", "Turbo", "Bluered", "RdBu", "RdYlBu", "RdYlGn", "Earth", "Portland", "Jet", "Hot", "Greys", "YlGnBu", "YlOrRd", "Picnic", "Rainbow", "Blackbody", "Electric", "Twilight", "Ice"];
var scaleOpts = "";
for (var si = 0; si < scaleNames.length; si++) scaleOpts += '<option value="' + scaleNames[si] + '">' + scaleNames[si] + "</option>";

h.push('<details><summary>Appearance &amp; titles</summary><div class="itv-grid">');
h.push('<label>Title<input type="text" data-role="title"></label>');
h.push('<label>Subtitle<input type="text" data-role="subtitle"></label>');
h.push('<label>Layout<select data-role="layout"><option value="rectangular">Rectangular</option><option value="slanted">Slanted</option><option value="circular">Circular</option><option value="fan">Fan</option></select></label>');
h.push('<label>Tree tip labels<input type="checkbox" data-role="tip_labels"></label>');
h.push('<label>Tip label size<input type="number" min="1" step="0.5" data-role="tip_label_size"></label>');
h.push('<label>Show internal nodes<input type="checkbox" data-role="show_internal_nodes"></label>');
h.push('<label>Internal node size<input type="number" min="1" step="0.5" data-role="internal_node_size"></label>');
h.push('<label>Show internal node IDs<input type="checkbox" data-role="show_internal_node_ids"></label>');
h.push('<label>Node ID font size<input type="number" min="1" step="0.5" data-role="internal_node_id_size"></label>');
h.push('<label>Tip node color<input type="color" data-role="tip_node_color"></label>');
h.push('<label>Tip node size<input type="number" min="1" step="0.5" data-role="tip_node_size"></label>');
h.push('<label>Internal node color<input type="color" data-role="internal_node_color"></label>');
h.push('<label>Branch color<input type="color" data-role="branch_color"></label>');
h.push('<label>Branch width<input type="number" min="0.1" step="0.1" data-role="branch_width"></label>');
h.push('<label>Heatmap colorscale<select data-role="heatmap_colorscale">' + scaleOpts + '</select></label>');
h.push('<label>Heatmap auto range<input type="checkbox" data-role="heatmap_auto"></label>');
h.push('<label>Heatmap min<input type="number" step="any" data-role="heatmap_min"></label>');
h.push('<label>Heatmap max<input type="number" step="any" data-role="heatmap_max"></label>');
h.push('<label>Heatmap axis title<input type="text" data-role="heatmap_title"></label>');
h.push('<label>Colorbar title<input type="text" data-role="heatmap_colorbar_title"></label>');
h.push('<label>Heatmap panel title<input type="text" data-role="heatmap_panel_title"></label>');
h.push('<label>Tree panel title<input type="text" data-role="tree_panel_title"></label>');
h.push('<label>Species panel title<input type="text" data-role="species_panel_title"></label>');
h.push('<label>Trait panel title<input type="text" data-role="trait_panel_title" placeholder="auto: trait name + value"></label>');
h.push('<label>Fan angle (deg)<input type="number" min="10" max="360" step="1" data-role="fan_angle"></label>');
h.push('<label>Show heatmap<input type="checkbox" data-role="show_heatmap"></label>');
h.push('</div></details>');
h.push('<details><summary>View sizing &amp; subviews</summary><div class="itv-grid">');
h.push('<div class="itv-wide itv-subhead">View sizing (blank = default)</div>');
h.push('<label>Heatmap view width (px)<input type="number" min="320" step="10" data-role="heatmap_view_width"></label>');
h.push('<label>Heatmap view height (px)<input type="number" min="240" step="10" data-role="heatmap_view_height"></label>');
h.push('<label>Species + trait view width (px)<input type="number" min="320" step="10" data-role="trait_view_width"></label>');
h.push('<label>Species + trait view height (px)<input type="number" min="240" step="10" data-role="trait_view_height"></label>');
h.push('<label>Default width (px)<input type="number" min="320" step="10" data-role="width"></label>');
h.push('<label>Default height (px)<input type="number" min="240" step="10" data-role="height"></label>');
h.push('<label>Aspect ratio (w/h)<input type="number" min="0.1" step="0.01" data-role="aspect_ratio"></label>');
h.push('<label>Responsive fill<input type="checkbox" data-role="responsive"></label>');
h.push('<div class="itv-wide itv-subhead">Subviews</div>');
h.push('<label>Auto-fit subview widths<input type="checkbox" data-role="auto_fit_subviews"></label>');
h.push('<label>Heatmap panel group alignment<select data-role="heatmap_subview_alignment"><option value="left">Left</option><option value="center">Center</option><option value="right">Right</option></select></label>');
h.push('<label>Tree panel width (px)<input type="number" min="100" step="10" data-role="heatmap_tree_width"></label>');
h.push('<label>Heatmap panel width (px)<input type="number" min="100" step="10" data-role="heatmap_panel_width"></label>');
h.push('<label>Species + trait panel group alignment<select data-role="trait_subview_alignment"><option value="left">Left</option><option value="center">Center</option><option value="right">Right</option></select></label>');
h.push('<label>Species + trait tree width (px)<input type="number" min="100" step="10" data-role="trait_tree_width"></label>');
h.push('<label>Tips panel width (px)<input type="number" min="60" step="10" data-role="tips_panel_width"></label>');
h.push('<label>Trait values panel width (px)<input type="number" min="100" step="10" data-role="trait_values_width"></label>');
h.push('<div class="itv-wide itv-hint">Auto-fit measures the species labels, reserves room for tree labels, and fits the main chart panel. The Species + trait view uses its dedicated tips panel instead of duplicate tree labels. Panel heights follow the view height to keep tip rows aligned.</div>');
h.push('</div></details>');

h.push('<details><summary>Filtering &amp; highlighting</summary><div class="itv-grid">');
h.push('<label>Trait filter<select data-role="filter_trait"><option value="-1">(none)</option>' + traitOpts + '</select></label>');
h.push('<label>Min<input type="number" step="any" data-role="filter_min"></label>');
h.push('<label>Max<input type="number" step="any" data-role="filter_max"></label>');
h.push('<label>Tip names (comma separated)<input type="text" data-role="filter_tips" placeholder="blank = all tips"></label>');
h.push('<label class="itv-wide">Clade filter (Ctrl/Cmd-click for multiple; none = all)<select multiple size="5" data-role="clade"></select></label>');
h.push('<div class="itv-wide"><div class="itv-subhead">Heatmap columns</div><div class="itv-checks" data-role="col_filter"></div></div>');
h.push('<label>Highlight trait<select data-role="highlight_trait"><option value="-1">(none)</option>' + traitOpts + '</select></label>');
h.push('<label>Highlight min<input type="number" step="any" data-role="highlight_min"></label>');
h.push('<label>Highlight max<input type="number" step="any" data-role="highlight_max"></label>');
h.push('<label>Highlight color<input type="color" data-role="highlight_color"></label>');
h.push('<div class="itv-wide itv-actions"><button type="button" data-role="clear_filters">Clear filters</button><button type="button" data-role="clear_highlight">Clear highlight</button></div>');
h.push('</div></details>');

h.push('<details><summary>Export</summary><div class="itv-actions">');
h.push('<button type="button" data-role="export_png">PNG</button>');
h.push('<button type="button" data-role="export_svg">SVG</button>');
h.push('<button type="button" data-role="export_pdf">PDF (print)</button>');
h.push('<button type="button" data-role="export_newick">Tree (Newick)</button>');
h.push('<button type="button" data-role="export_csv">Heatmap CSV</button>');
h.push('<button type="button" data-role="export_config">Configuration JSON</button>');
h.push('<label class="itv-file">Load configuration JSON<input type="file" accept=".json,application/json" data-role="import_config"></label>');
h.push('<span class="itv-status" data-role="export_status"></span>');
h.push('</div></details>');

h.push('<details><summary>Presets &amp; settings</summary><div class="itv-grid">');
h.push('<label>Preset<select data-role="preset_select"></select></label>');
h.push('<label>Preset name<input type="text" data-role="preset_name" placeholder="My preset"></label>');
h.push('<div class="itv-wide itv-actions"><button type="button" data-role="preset_apply">Apply preset</button><button type="button" data-role="preset_save">Save current as preset</button><button type="button" data-role="preset_delete">Delete saved preset</button><button type="button" data-role="preset_restore">Restore defaults</button><span class="itv-status" data-role="preset_status"></span></div>');
h.push('</div></details>');
h.push('<div class="itv-hint">Click an internal-node diamond to flip its child clades. Hover tips, nodes, heatmap cells, or trait dots for details.</div>');

ui = document.createElement("div");
ui.className = "itv-ui";
ui.innerHTML = h.join("");

var wrap = document.createElement("div");
wrap.className = "itv-wrap";
el.parentNode.insertBefore(wrap, el);
wrap.appendChild(ui);
wrap.appendChild(el);
document.head.appendChild(styleEl);

// In a standalone htmlwidgets page the container has no explicit height, so a
// "100%" responsive fill collapses. Give it the viewport height there only
// (never inside an R Markdown document, where knitr controls the size).
var widgetContainer = document.getElementById("htmlwidget_container");
if (widgetContainer && widgetContainer.parentNode === document.body) {
  document.documentElement.style.height = "100%";
  document.body.style.height = "100%";
  document.body.style.margin = "0";
  widgetContainer.style.height = "100%";
}

refreshCladeOptions();
refreshColumnFilter();
refreshPresetOptions();
syncUI();
bindEvents();
applyContainerSizing();
redraw();

window.addEventListener("resize", function () {
  if (cfg.responsive) {
    applyContainerSizing();
    Plotly.Plots.resize(el);
  }
});

el.on("plotly_click", function (data) {
  var pts = data.points || [];
  for (var i = 0; i < pts.length; i++) {
    var p = pts[i];
    if (p.data && p.data.meta === "internal" && p.customdata !== undefined) {
      var node = p.customdata;
      state.rev[node] = !state.rev[node];
      redraw();
      return;
    }
  }
});

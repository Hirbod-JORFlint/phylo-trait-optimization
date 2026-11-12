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
var costScaled = !!TD.costScaled;

function mergeCfg(base, over) {
  var out = {};
  for (var k in base) out[k] = base[k];
  if (over) for (var k2 in over) out[k2] = over[k2];
  return out;
}
function numOrNull(v) {
  if (v === undefined || v === null || v === "") return null;
  var x = +v;
  return isFinite(x) ? x : null;
}
function asArray(v) {
  if (v === undefined || v === null) return null;
  return Array.isArray(v) ? v : [v];
}
function inRange(v, lo, hi) {
  if (!isFinite(v)) return false;
  if (lo !== null && v < lo) return false;
  if (hi !== null && v > hi) return false;
  return true;
}
function esc(s) {
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

var DEFAULT_CFG = mergeCfg({
  layout: "rectangular",
  title: "Interactive phylogenetic tree + trait heatmap",
  subtitle: "",
  theme: "light",
  tip_labels: true,
  tip_label_size: 10,
  show_internal_nodes: true,
  internal_node_size: 9,
  show_internal_node_ids: false,
  internal_node_id_size: 9,
  tip_node_color: "#1f77b4",
  tip_node_size: 7,
  internal_node_color: "#e08a1e",
  branch_color: "#8a8a8a",
  branch_width: 1.4,
  heatmap_colorscale: "Viridis",
  heatmap_range: null,
  heatmap_title: "Trait",
  heatmap_colorbar_title: "value",
  heatmap_panel_title: "Trait heatmap",
  tree_panel_title: "Phylogenetic tree",
  species_panel_title: "Species",
  trait_panel_title: "",
  responsive: true,
  width: null,
  height: null,
  aspect_ratio: null,
  tip_filter: null,
  clade: null,
  heatmap_columns: null,
  highlight: null,
  fan_angle: 300,
  show_heatmap: true
}, TD.defaults || {});
var cfg = mergeCfg(DEFAULT_CFG, TD.config || {});

var state = {
  rev: {},
  viewMode: "heat",
  traitIdx: 0,
  visibleTips: {},
  visibleCols: [],
  filter: { traitIdx: null, min: null, max: null, tips: null },
  clades: [],
  highlight: { traitIdx: null, min: null, max: null, tips: null, color: "#e41a1c" },
  presetValue: "",
  lastLayout: null
};
var matchCache = {};
var ui = null;

function initRev() {
  state.rev = {};
  for (var k in children) state.rev[k] = false;
}

function recomputeVisibility() {
  for (var t = 1; t <= n; t++) state.visibleTips[t] = true;
  if (state.filter.tips && state.filter.tips.length) {
    var keep = {};
    for (var i = 0; i < state.filter.tips.length; i++) {
      var idx = tips.indexOf(state.filter.tips[i]);
      if (idx >= 0) keep[idx + 1] = true;
    }
    for (var t2 = 1; t2 <= n; t2++) if (!keep[t2]) state.visibleTips[t2] = false;
  }
  if (state.filter.traitIdx !== null) {
    for (var t3 = 1; t3 <= n; t3++) {
      var v = traitsHeat[t3 - 1][state.filter.traitIdx];
      if (!inRange(v, state.filter.min, state.filter.max)) state.visibleTips[t3] = false;
    }
  }
  if (state.clades.length) {
    var inClade = {};
    for (var c = 0; c < state.clades.length; c++) {
      (function mark(node) {
        var ch = children[node];
        if (!ch) { inClade[node] = true; return; }
        for (var j = 0; j < ch.length; j++) mark(ch[j]);
      })(state.clades[c]);
    }
    for (var t4 = 1; t4 <= n; t4++) if (!inClade[t4]) state.visibleTips[t4] = false;
  }
  matchCache = {};
}

function initStateFromConfig() {
  for (var t = 1; t <= n; t++) state.visibleTips[t] = true;
  state.visibleCols = [];
  for (var j = 0; j < tNames.length; j++) state.visibleCols.push(j);
  state.filter = { traitIdx: null, min: null, max: null, tips: null };
  state.clades = [];
  state.highlight = { traitIdx: null, min: null, max: null, tips: null, color: "#e41a1c" };

  var tf = cfg.tip_filter;
  if (tf) {
    tf = asArray(tf);
    if (Array.isArray(tf) && tf.length && typeof tf[0] !== "object") tf = { tips: tf };
    else if (Array.isArray(tf)) tf = tf[0] || {};
    if (tf.trait !== undefined && tf.trait !== null) {
      var ti = tNames.indexOf(tf.trait);
      if (ti >= 0) {
        state.filter.traitIdx = ti;
        state.filter.min = numOrNull(tf.min);
        state.filter.max = numOrNull(tf.max);
      }
    }
    if (tf.tips) state.filter.tips = asArray(tf.tips).slice();
  }
  var cl = asArray(cfg.clade);
  if (cl) {
    state.clades = cl.map(Number).filter(function (x) { return isFinite(x); });
  }
  var hc = asArray(cfg.heatmap_columns);
  if (hc) {
    var cols = [];
    for (var i = 0; i < hc.length; i++) {
      var idx = (typeof hc[i] === "number") ? hc[i] : tNames.indexOf(hc[i]);
      if (idx >= 0 && idx < tNames.length && cols.indexOf(idx) < 0) cols.push(idx);
    }
    state.visibleCols = cols;
  }
  var hl = cfg.highlight;
  if (hl) {
    hl = asArray(hl);
    if (Array.isArray(hl) && hl.length && typeof hl[0] !== "object") hl = { tips: hl };
    else if (Array.isArray(hl)) hl = hl[0] || {};
    if (hl.trait !== undefined && hl.trait !== null) {
      var hi = tNames.indexOf(hl.trait);
      if (hi >= 0) {
        state.highlight.traitIdx = hi;
        state.highlight.min = numOrNull(hl.min);
        state.highlight.max = numOrNull(hl.max);
      }
    }
    if (hl.tips) state.highlight.tips = asArray(hl.tips).slice();
    if (hl.color) state.highlight.color = hl.color;
  }
  recomputeVisibility();
}

function visibleTipCount() {
  var c = 0;
  for (var t = 1; t <= n; t++) if (state.visibleTips[t]) c++;
  return c;
}

function tipMatches(t) {
  if (state.highlight.tips && state.highlight.tips.length) {
    return state.highlight.tips.indexOf(tips[t - 1]) >= 0;
  }
  if (state.highlight.traitIdx === null) return false;
  var v = traitsHeat[t - 1][state.highlight.traitIdx];
  return inRange(v, state.highlight.min, state.highlight.max);
}

function cladeMatches(node) {
  if (matchCache[node] !== undefined) return matchCache[node];
  var ch = children[node];
  if (!ch) return matchCache[node] = (!!state.visibleTips[node] && tipMatches(node));
  var any = false;
  for (var i = 0; i < ch.length; i++) {
    if (cladeMatches(ch[i])) { any = true; break; }
  }
  return matchCache[node] = any;
}

function buildFiltered(revMap) {
  var fchildren = {};
  var flen = {};
  function rec(node) {
    var ch = children[node];
    if (!ch) {
      if (!state.visibleTips[node]) return null;
      return { id: node, len: (edgeLen[node] || 1), children: null };
    }
    var ordered = revMap[node] ? ch.slice().reverse() : ch.slice();
    var kept = [];
    for (var i = 0; i < ordered.length; i++) {
      var r = rec(ordered[i]);
      if (r) kept.push(r);
    }
    if (kept.length === 0) return null;
    if (kept.length === 1) {
      return { id: kept[0].id, len: (edgeLen[node] || 1) + kept[0].len, children: kept[0].children };
    }
    var ids = [];
    for (var j = 0; j < kept.length; j++) {
      ids.push(kept[j].id);
      flen[kept[j].id] = kept[j].len;
    }
    fchildren[node] = ids;
    return { id: node, len: (edgeLen[node] || 1), children: ids };
  }
  var top = rec(root);
  if (!top) return null;
  return { root: top.id, children: fchildren, len: flen };
}

function leafOrderFiltered(F) {
  var order = [];
  (function dfs(node) {
    var ch = F.children[node];
    if (!ch) { order.push(node * 1); return; }
    for (var i = 0; i < ch.length; i++) dfs(ch[i]);
  })(F.root);
  return order;
}

function layoutCalcFiltered(F) {
  var order = leafOrderFiltered(F);
  var m = order.length;
  var y = {}, slot = {};
  for (var i = 0; i < m; i++) { slot[order[i]] = i + 1; y[order[i]] = i + 1; }
  var mn = {}, mx = {};
  (function ranges(node) {
    var ch = F.children[node];
    if (!ch) { mn[node] = slot[node]; mx[node] = slot[node]; return; }
    var lo = Infinity, hi = -Infinity;
    for (var j = 0; j < ch.length; j++) {
      ranges(ch[j]);
      if (mn[ch[j]] < lo) lo = mn[ch[j]];
      if (mx[ch[j]] > hi) hi = mx[ch[j]];
    }
    mn[node] = lo; mx[node] = hi;
    y[node] = (lo + hi) / 2;
  })(F.root);
  var x = {};
  x[F.root] = 0;
  (function xs(node) {
    var ch = F.children[node] || [];
    for (var j = 0; j < ch.length; j++) {
      x[ch[j]] = x[node] + (F.len[ch[j]] || 1);
      xs(ch[j]);
    }
  })(F.root);
  return { x: x, y: y, order: order, root: F.root, children: F.children, len: F.len, m: m, slot: slot, angle: null };
}

function circularLayout(L) {
  var m = L.m;
  var fan = cfg.layout === "fan";
  var span = fan ? (Math.min(360, Math.max(10, +cfg.fan_angle || 300)) * Math.PI / 180) : 2 * Math.PI;
  var start = -Math.PI / 2;
  var angle = {}, count = {};
  for (var i = 0; i < m; i++) angle[L.order[i]] = start + ((i + 0.5) / m) * span;
  function calc(node) {
    var ch = L.children[node];
    if (!ch) { count[node] = 1; return angle[node]; }
    var sx = 0, sy = 0, c = 0;
    for (var j = 0; j < ch.length; j++) {
      var a = calc(ch[j]);
      var w = count[ch[j]];
      sx += w * Math.cos(a); sy += w * Math.sin(a); c += w;
    }
    count[node] = c;
    var a2 = Math.atan2(sy, sx);
    angle[node] = a2;
    return a2;
  }
  calc(L.root);
  var maxD = 0;
  for (var k in L.x) if (L.x[k] > maxD) maxD = L.x[k];
  if (!(maxD > 0)) maxD = 1;
  var cx = {}, cy = {}, r = {};
  for (var k2 in L.x) {
    r[k2] = L.x[k2] / maxD;
    cx[k2] = r[k2] * Math.cos(angle[k2]);
    cy[k2] = r[k2] * Math.sin(angle[k2]);
  }
  return { x: cx, y: cy, angle: angle, r: r, order: L.order, root: L.root, children: L.children, m: m, slot: L.slot, radial: true };
}

function leafCountOriginal(node) {
  var ch = children[node];
  if (!ch) return 1;
  var s = 0;
  for (var i = 0; i < ch.length; i++) s += leafCountOriginal(ch[i]);
  return s;
}

function leafCountFiltered(L, node) {
  var ch = L.children[node];
  if (!ch) return 1;
  var s = 0;
  for (var i = 0; i < ch.length; i++) s += leafCountFiltered(L, ch[i]);
  return s;
}

function themeColors() {
  if (cfg.theme === "dark") {
    return { bg: "#0f1115", paper: "#0f1115", font: "#e8e8e8", muted: "#a8b0b8",
      grid: "#2a3038", axis: "#444c55", nodeEdge: "#f0f0f0", uiBg: "#171b21",
      uiFg: "#e8e8e8", uiBorder: "#2f363f", uiInput: "#0f1115", uiAccent: "#4da3ff" };
  }
  return { bg: "#ffffff", paper: "#ffffff", font: "#222222", muted: "#666666",
    grid: "#eeeeee", axis: "#cccccc", nodeEdge: "#ffffff", uiBg: "#f7f7f7",
    uiFg: "#222222", uiBorder: "#dddddd", uiInput: "#ffffff", uiAccent: "#1f77b4" };
}

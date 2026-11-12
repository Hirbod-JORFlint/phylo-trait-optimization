(function () {
  var out = { steps: [] };
  function q(role) { return document.querySelector('[data-role="' + role + '"]'); }
  function setVal(role, val) {
    var e = q(role);
    if (!e) { out.steps.push("missing:" + role); return; }
    if (e.type === "checkbox") e.checked = !!val; else e.value = val;
    e.dispatchEvent(new Event("change", { bubbles: true }));
  }
  function plotEl() { return document.querySelector(".js-plotly-plot"); }
  function full() { return plotEl()._fullLayout || {}; }
  function metas() { return (plotEl()._fullData || []).map(function (d) { return d.meta || d.type; }); }
  function status() { var e = q("status"); return e ? e.textContent : ""; }

  var steps = [];
  function step(name, fn) { steps.push({ name: name, fn: fn }); }

  step("initial", function () {
    var p = plotEl();
    out.steps.push({ name: "initial", metas: metas(), paper: full().paper_bgcolor,
      status: status(), width: p.style.width, height: p.style.height });
  });
  step("dark", function () {
    q("theme").click();
    out.steps.push({ name: "dark", cls: document.querySelector(".itv-ui").className,
      paper: full().paper_bgcolor, btn: q("theme").textContent });
  });
  step("circular", function () {
    setVal("layout", "circular");
    out.steps.push({ name: "circular", xrange: full().xaxis.range, yrange: full().yaxis.range,
      metas: metas(), heatmap: !!document.querySelector(".hm") });
  });
  step("fan", function () {
    setVal("layout", "fan");
    out.steps.push({ name: "fan", xrange: full().xaxis.range, heatmap: !!document.querySelector(".hm") });
  });
  step("slanted", function () {
    setVal("layout", "slanted");
    out.steps.push({ name: "slanted", metas: metas() });
  });
  step("rectangular", function () {
    setVal("layout", "rectangular");
    out.steps.push({ name: "rectangular", metas: metas() });
  });
  step("filter", function () {
    setVal("filter_trait", "0");
    setVal("filter_min", "0.3");
    setVal("filter_max", "0.7");
    out.steps.push({ name: "filter", status: status(),
      y2range: full().yaxis2 && full().yaxis2.range,
      hmHeight: document.querySelector(".hm") ? Math.round(document.querySelector(".hm").getBoundingClientRect().height) : 0 });
  });
  step("highlight", function () {
    setVal("highlight_trait", "1");
    setVal("highlight_min", "0.5");
    setVal("highlight_max", "1");
    out.steps.push({ name: "highlight", status: status(), metas: metas() });
  });
  step("traitView", function () {
    q("view").click();
    var p = plotEl();
    out.steps.push({ name: "traitView", metas: metas(),
      width: p.style.width, height: p.style.height,
      x3: full().xaxis3 && full().xaxis3.domain,
      y2: full().yaxis2 && full().yaxis2.domain,
      y2range: full().yaxis2 && full().yaxis2.range });
  });
  step("heatView", function () {
    q("view").click();
    var p = plotEl();
    out.steps.push({ name: "heatView", metas: metas(),
      width: p.style.width, height: p.style.height });
  });
  step("presetDark", function () {
    setVal("preset_select", "builtin:dark");
    q("preset_apply").click();
    out.steps.push({ name: "presetDark", paper: full().paper_bgcolor,
      cls: document.querySelector(".itv-ui").className });
  });
  step("presetDefault", function () {
    setVal("preset_select", "builtin:default");
    q("preset_apply").click();
    out.steps.push({ name: "presetDefault", paper: full().paper_bgcolor,
      cls: document.querySelector(".itv-ui").className });
  });
  step("savePreset", function () {
    setVal("preset_name", "ui test preset");
    q("preset_save").click();
    var sel = q("preset_select");
    var names = [];
    for (var i = 0; i < sel.options.length; i++) names.push(sel.options[i].value);
    out.steps.push({ name: "savePreset", selected: sel.value,
      hasSaved: names.indexOf("saved:ui test preset") >= 0,
      status: q("preset_status").textContent });
  });
  step("exports", function () {
    var downloads = [];
    var origClick = HTMLAnchorElement.prototype.click;
    HTMLAnchorElement.prototype.click = function () { downloads.push(this.download); };
    q("export_newick").click();
    q("export_csv").click();
    q("export_config").click();
    HTMLAnchorElement.prototype.click = origClick;
    var images = [];
    var origDI = Plotly.downloadImage;
    Plotly.downloadImage = function (el, opts) { images.push(opts.format); };
    q("export_png").click();
    q("export_svg").click();
    Plotly.downloadImage = origDI;
    out.steps.push({ name: "exports", downloads: downloads, images: images,
      status: q("export_status").textContent });
  });
  step("clearFilters", function () {
    q("clear_filters").click();
    q("clear_highlight").click();
    out.steps.push({ name: "clearFilters", status: status(),
      y2range: full().yaxis2 && full().yaxis2.range });
  });

  var i = 0;
  function run() {
    if (i >= steps.length) {
      document.body.setAttribute("data-test", JSON.stringify(out));
      return;
    }
    var s = steps[i++];
    try { s.fn(); } catch (e) { out.steps.push({ name: s.name, error: String(e) }); }
    setTimeout(run, 350);
  }
  setTimeout(run, 2000);
})();

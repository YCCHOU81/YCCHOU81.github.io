"use strict";

const countNames = ["a", "b", "c"];
const order = counts => 1 + counts.reduce((sum, value) => sum + value, 0);
const validCounts = counts => counts.length === 3 &&
  counts.every(value => Number.isSafeInteger(value) && value >= 0) &&
  Number.isSafeInteger(order(counts));
const query = new URLSearchParams(location.search);
let selectedCounts = countNames.map((name, i) => query.has(name) ? Number(query.get(name)) : [1, 0, 0][i]);
if (!validCounts(selectedCounts)) selectedCounts = [1, 0, 0];

function refreshSharedInformation() {
  const params = new URLSearchParams(location.search);
  countNames.forEach((name, i) => params.set(name, selectedCounts[i]));
  const current = new URLSearchParams(location.search);
  for (const key of ["seed", "decoration", "qt", "qu", "ql", "qm", "modeA", "modeB", "countA", "countB", "qAB", "qAC", "qAD", "qBC", "qBD", "qCD", "bouquetCore", "bq0", "bq1", "bq2", "bq3", "bq4", "bq5"]) if (current.has(key)) params.set(key, current.get(key));
  document.querySelectorAll("a[data-page]").forEach(link => {
    const next = new URLSearchParams(params);
    if (link.dataset.page === "explore.html" && current.has("r") && current.has("s")) {
      next.set("r", current.get("r")); next.set("s", current.get("s"));
    }
    link.href = window.AtlasUI?.pageHref(link.dataset.page,next) ?? `${link.dataset.page}?${next}`;
  });
  const context = document.getElementById("selection-context");
  if (context) context.textContent = `Selected attachment: (${selectedCounts.join(", ")}). Order n = ${order(selectedCounts)}.`;
}

// TV has one explorer implementation, remounted in place for each seed.
// Keep display preferences separate from each family's bounded result cache.
refreshSharedInformation();
if (document.getElementById("monodromy-source")) {
  const settingIds = {view:"explore-view",quality:"quality",cuspT:"cusp-position",cuspDepth:"cusp-depth"};
  const preferred = Object.fromEntries(Object.entries(settingIds).map(([key,id]) => [key,query.get(key) ?? document.getElementById(id).value]));
  if (!query.has("quality")) preferred.quality = "8";
  let mounted = null, factory = null, applied = null, hiddenSource = null;
  let cacheReady = !window.AtlasResults?.ready?.then;
  const defaultClosed = node => node.getAttribute?.("data-default-closed") != null;
  const helpKey = node => node.id ? "id:" + node.id : "summary:" + (node.querySelector("summary")?.textContent.trim() ?? "");
  const readControls = () => Object.fromEntries(Object.entries(settingIds).map(([key,id]) => [key,document.getElementById(id).value]));

  function rememberSettings() {
    for (const [key,value] of Object.entries(readControls())) {
      if (!applied || value !== applied[key]) preferred[key] = value;
    }
  }
  function savePageDisplay() {
    rememberSettings();
    const detailsByKey = Object.fromEntries([...document.querySelectorAll("details")]
      .filter(node => !defaultClosed(node)).map(node => [helpKey(node),node.open]));
    window.AtlasResults?.set("explore-display-v1","current",{settings:{...preferred},detailsByKey});
  }
  function restorePageDisplay() {
    const saved = window.AtlasResults?.get("explore-display-v1","current");
    if (saved?.settings) Object.assign(preferred,saved.settings);
    return saved;
  }
  function mount(detail = {}, initial = false) {
    if (!cacheReady || !factory) return;
    if (validCounts(detail.counts ?? [])) selectedCounts = [...detail.counts];
    const details = [...document.querySelectorAll("details")].filter(node => !defaultClosed(node))
      .map(node => ({node,open:node.open,id:node.id,parent:node.parentElement?.id}));
    const recommendation = document.getElementById("preview-recommendation");
    if (recommendation) recommendation.textContent = "Apply & compute chooses a recommended preview for this graph. Saved results are restored automatically; the controls below refine either image.";
    mounted = factory({source:detail.source,settings:{...preferred},initial});
    applied = readControls();
    for (const state of details) {
      const node = state.node.isConnected ? state.node : state.id ? document.getElementById(state.id)
        : state.parent ? document.getElementById(state.parent)?.querySelector("details") : null;
      if (node) node.open = state.open;
    }
  }
  window.AtlasExplorer = {
    register(kind,nextFactory) {
      if (kind !== "tv") return;
      factory = nextFactory;
      if (!mounted) mount({},true);
    },
    get source() {return mounted?.source?.() ?? null;}
  };
  window.addEventListener("atlas:beforefamilychange",event => {
    rememberSettings();
    event.detail.source = mounted?.source?.() ?? event.detail.source;
    mounted?.dispose();
    mounted = null;
  });
  window.addEventListener("atlas:familychange",event => {
    // Direct change events have the same lifecycle as a before/change pair.
    if (mounted) {rememberSettings();mounted.dispose();mounted = null;}
    mount(event.detail ?? {});
  });
  window.addEventListener("pagehide",() => {
    savePageDisplay();
    hiddenSource = mounted?.source?.() ?? hiddenSource;
    mounted?.dispose();
    mounted = null;
  });
  window.addEventListener("atlas:cachedpage",() => {
    if (!mounted) mount({source:hiddenSource,counts:[...selectedCounts]});
  });
  window.addEventListener("atlas:beforepagechange",() => {savePageDisplay();mounted?.remember?.();});
  function applySavedDetails(saved) {
    for (const node of document.querySelectorAll("details")) {
      if (defaultClosed(node)) node.open = false;
      else {
        const open = saved?.detailsByKey?.[helpKey(node)];
        if (typeof open === "boolean") node.open = open;
      }
    }
  }
  if (cacheReady) applySavedDetails(restorePageDisplay());
  else Promise.resolve(window.AtlasResults.ready).catch(() => null).then(() => {
    const saved = restorePageDisplay();
    cacheReady = true;
    if (!mounted) mount({},true);
    applySavedDetails(saved);
  });
}

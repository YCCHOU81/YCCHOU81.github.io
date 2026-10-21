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
  const params = new URLSearchParams();
  countNames.forEach((name, i) => params.set(name, selectedCounts[i]));
  const current = new URLSearchParams(location.search);
  for (const key of ["seed", "decoration", "qt", "qu", "ql", "qm", "modeA", "modeB", "countA", "countB", "qAB", "qAC", "qAD", "qBC", "qBD", "qCD"]) if (current.has(key)) params.set(key, current.get(key));
  document.querySelectorAll("a[data-page]").forEach(link => {
    const next = new URLSearchParams(params);
    if (link.dataset.page === "explore.html" && current.has("r") && current.has("s")) {
      next.set("r", current.get("r")); next.set("s", current.get("s"));
    }
    link.href = `${link.dataset.page}?${next}`;
  });
  const context = document.getElementById("selection-context");
  if (context) context.textContent = `Selected attachment: (${selectedCounts.join(", ")}). Order n = ${order(selectedCounts)}.`;
}

// Nearest-point projection in (r,s) onto K_delta, delta = .04.
function projectCompact(r, s) {
  const low = 0.04, high = 0.48, diagonal = 0.52;
  if (![r, s].every(Number.isFinite)) return {r: 0.31, s: 0.29};
  if (r <= high && s <= high && r + s >= diagonal) return {r, s};
  const clamp = x => Math.max(low, Math.min(high, x));
  const x = clamp((r - s + diagonal) / 2);
  const candidates = [{r: high, s: clamp(s)}, {r: clamp(r), s: high}, {r: x, s: diagonal - x}];
  return candidates.reduce((best, p) =>
    (p.r - r) ** 2 + (p.s - s) ** 2 < (best.r - r) ** 2 + (best.s - s) ** 2 ? p : best);
}

const sourcePosition = (r, s) => window.AtlasView.sourceAtlasXY(r, s);

refreshSharedInformation();
if (document.getElementById("attachment-form") && !window.SeedCases?.isNew) connectExplorer();

function connectExplorer() {
  const get = id => document.getElementById(id);
  const source = get("monodromy-source");
  const attachedSource = get("monodromy-attached");
  get("explore-view").value = query.get("view") === "overlay" ? "overlay" : "separate";
  const initialPair = Number(query.get("cuspT"));
  get("cusp-position").value = initialPair >= .02 && initialPair <= .48 ? String(initialPair) : "0.2";
  get("cusp-depth").value = query.get("cuspDepth") === "10" ? "10" : "8";
  const meshCache = new Map();
  const cuspCache = new Map();
  let cuspWorkers = [], cuspVersion = 0, cuspTraces = Array(6).fill(null);
  let meshWorker = null, meshVersion = 0, currentMesh = null;
  let pointWorker = null, pointVersion = 0, pointBusy = false, queuedPoint = null;
  let activePoint = null, pointTimer = null, pointResult = null, pointerId = null, pointerSource = null;
  let point = projectCompact(query.has("r") ? Number(query.get("r")) : 0.31,
    query.has("s") ? Number(query.get("s")) : 0.29);
  countNames.forEach((name, i) => { get(name).value = selectedCounts[i]; });

  function message(id, text, error = false) {
    get(id).textContent = text;
    get(id).classList.toggle("error", error);
  }
  function redraw() {
    const overlap = get("explore-view").value === "overlay";
    get("separate-plots").hidden = overlap; get("overlap-plot").hidden = !overlap;
    const resolved = branch => branch?.status === "ok" && window.AtlasView.finiteTau(branch.tau) ? branch.tau : null;
    const base = resolved(pointResult?.base), attached = resolved(pointResult?.result);
    const baseLayer = {mesh: currentMesh, field: "base", color: "#285f65"};
    const attachedLayer = {mesh: currentMesh, field: "result", color: "#b55e3c"};
    const tracePaths = field => cuspTraces.filter(Boolean).map(trace => ({color: field === "base" ? "#285f65" : "#b55e3c", width: 1.5, radius: 1.6,
      points: trace.nodes.map(node => cuspTau(node, field))}));
    const traceBands = field => [0, 1, 2].flatMap(edge => {
      const left = cuspTraces[2 * edge], right = cuspTraces[2 * edge + 1];
      if (!left || !right) return [];
      const paired = new Map(right.nodes.map(node => [node.epsilon, node]));
      return [{color: field === "base" ? "#285f65" : "#b55e3c", left: left.nodes.map(node => cuspTau(node, field)),
        right: left.nodes.map(node => cuspTau(paired.get(node.epsilon), field))}];
    });
    const settings = {height: 520};
    if (overlap) {
      const chart = window.AtlasView.drawDomain(get("overlay-image"), {...settings,
        meshLayers: [baseLayer, attachedLayer], points: [{tau: base, color: baseLayer.color}, {tau: attached, color: attachedLayer.color}].filter(p => p.tau),
        paths: [...tracePaths("base"), ...tracePaths("result")],
        bands: [...traceBands("base"), ...traceBands("result")],
        arrow: base && attached ? {from: base, to: attached} : null});
      const distance = chart && base && attached ? Math.hypot(...chart.map(base).map((v, i) => v - chart.map(attached)[i])) : 0;
      get("overlap-note").textContent = !base || !attached ? "The attachment arrow appears when both selected period ratios resolve." :
        distance <= 12 ? "The selected points coincide or are too close to separate at this scale." : "The solid arrow compares the two selected period ratios; it is not a continuation path.";
    } else {
      window.AtlasView.drawDomain(get("base-image"), {...settings, meshLayers: [baseLayer], points: [{tau: base, color: baseLayer.color}], paths: tracePaths("base"), bands: traceBands("base")});
      window.AtlasView.drawDomain(get("result-image"), {...settings, meshLayers: [attachedLayer], points: [{tau: attached, color: attachedLayer.color}], paths: tracePaths("result"), bands: traceBands("result")});
    }
  }

  function cuspTau(node, branch) {
    const result = node?.[branch];
    if (result?.status !== "ok" || !result.tauExact) return null;
    const tau = {re: Number(result.tauExact[0]), im: Number(result.tauExact[1])};
    return window.AtlasView.inMarkedDomain(tau) ? tau : null;
  }
  function updateCuspEndpoints() {
    const fmt = number => Number(number.toPrecision(10)).toString();
    get("cusp-endpoints").textContent = cuspTraces.map((trace, index) => {
      const edge = Math.floor(index / 2);
      if (!trace) return `Cusp ${window.AtlasView.cuspForEdge[edge]}, series ${index % 2 + 1}: waiting for computation.`;
      const describe = branch => {
        const node = trace.nodes.findLast(node => cuspTau(node, branch));
        const tau = cuspTau(node, branch);
        return tau ? `${branch === "base" ? "base" : "attachment"}: Re τ = ${fmt(tau.re)}, Im τ = ${fmt(tau.im)} at ε = ${node.epsilon}` : `${branch}: unresolved`;
      };
      return `Cusp ${window.AtlasView.cuspForEdge[edge]}, t = ${trace.position} · ${trace.status}\n${describe("base")}\n${describe("result")}`;
    }).join("\n\n");
  }
  function stopCuspApproaches(text = "Stopped. Only computed samples are retained.") {
    ++cuspVersion; cuspWorkers.forEach(worker => worker.terminate()); cuspWorkers = [];
    cuspTraces.forEach(trace => {if (trace?.status === "computing") trace.status = "stopped";});
    updateCuspEndpoints();
    get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
    message("cusp-status", text);
  }
  function startCuspApproaches(event) {
    event?.preventDefault();
    if (!get("cusp-controls").reportValidity()) return;
    const position = String(Number(get("cusp-position").value));
    if (!(Number(position) >= .02 && Number(position) <= .48)) return;
    const positions = [position, String(1 - Number(position))];
    const target = get("cusp-depth").value === "10" ? "1e-10" : "1e-8";
    stopCuspApproaches(); const id = cuspVersion;
    const counts = [...selectedCounts], sum = counts.reduce((a, b) => a + b, 0);
    const tailRefinement = sum >= 8 ? 2 : sum >= 4 ? 1 : 0;
    const key = `computed-cusp-pairs-v3/${counts.join(",")}/${position}/${target}/${tailRefinement}`;
    let stored = event ? null : cuspCache.get(key);
    if (!event && !stored) {try {stored = JSON.parse(window.sessionStorage.getItem(key));} catch { /* Cache is optional. */ }}
    if (Array.isArray(stored) && stored.length === 6 && stored.every(trace => typeof trace?.status === "string" && Array.isArray(trace.nodes))) {
      cuspTraces = stored; updateCuspEndpoints(); redraw();
      message("cusp-status", `${stored.every(trace => trace.status === "ok") ? "Six computed cusp approaches" : "Partial computed approaches"} reused for (${counts.join(", ")}), t = ${positions.join(" and ")}. Recompute to run again.`); return;
    }
    cuspTraces = Array(6).fill(null); updateCuspEndpoints(); redraw();
    get("cusp-run").disabled = true; get("cusp-cancel").disabled = false;
    const url = new URL(location.href); url.searchParams.set("cuspT", position); url.searchParams.set("cuspDepth", get("cusp-depth").value); history.replaceState(null, "", url);
    const cores = Number(window.navigator?.hardwareConcurrency);
    const poolSize = cores > 0 && cores <= 2 ? 1 : 2;
    const slots = Array.from({length: poolSize}, () => ({worker: null, job: null}));
    let nextIndex = 0, completed = 0;
    function closeSlot(slot) {
      const worker = slot.worker;
      worker?.terminate(); cuspWorkers = cuspWorkers.filter(active => active !== worker);
      slot.worker = null; slot.job = null;
    }
    function finishAll() {
      if (id !== cuspVersion || completed !== 6) return;
      slots.forEach(closeSlot);
      get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
      const complete = cuspTraces.every(trace => trace?.status === "ok");
      if (cuspTraces.some(trace => trace.nodes.some(node => cuspTau(node, "base") || cuspTau(node, "result")))) {
        cuspCache.set(key, structuredClone(cuspTraces));
        if (cuspCache.size > 6) cuspCache.delete(cuspCache.keys().next().value);
        try {window.sessionStorage.setItem(key, JSON.stringify(cuspTraces));} catch { /* Quota does not affect computation. */ }
      }
      message("cusp-status", complete ? `Six computed approaches complete for (${counts.join(", ")}), t = ${positions.join(" and ")}. ${tailRefinement ? `${tailRefinement} extra tail samples per series for this attachment. ` : ""}No extrapolated endpoints.` : "Paired cusp computations finished with partial traces. Ribbons stop at their last shared resolved sample.");
      updateCuspEndpoints();
    }
    function finishJob(slot) {
      slot.job = null; completed++;
      updateCuspEndpoints(); redraw();
      if (completed === 6) finishAll(); else next(slot);
    }
    function next(slot) {
      if (id !== cuspVersion || slot.job || nextIndex === 6) return;
      const index = nextIndex++;
      const edge = Math.floor(index / 2), series = index % 2, pathPosition = positions[series], jobId = id * 8 + index;
      cuspTraces[index] = {nodes: [], status: "computing", edge, position: pathPosition};
      message("cusp-status", `Computing cusp ${window.AtlasView.cuspForEdge[edge]}, series ${series + 1}/2 (${index + 1}/6), t = ${pathPosition}…`);
      try {
        slot.job = {id: jobId, index, edge, position: pathPosition};
        if (!slot.worker) {
          const worker = new Worker("limits-worker.js?v=20261021-dual1");
          slot.worker = worker; cuspWorkers.push(worker);
          worker.onmessage = ({data}) => {
            const job = slot.job;
            if (id !== cuspVersion || slot.worker !== worker || !cuspWorkers.includes(worker) || !job || data.id !== job.id) return;
            const index = job.index, edge = job.edge, series = index % 2;
            if (data.type === "progress") {
              cuspTraces[index].nodes.push(data.node);
              message("cusp-status", `Cusp ${window.AtlasView.cuspForEdge[edge]}, series ${series + 1}/2 (${index + 1}/6): ${data.completed}/${data.total}; ε ≈ ${Number(data.node.epsilon).toExponential(2)}.`);
              updateCuspEndpoints(); redraw();
            }
            if (data.type === "done") {cuspTraces[index] = {nodes: data.result.nodes, status: data.result.status, edge, position: job.position}; finishJob(slot);}
            if (data.type === "error") {cuspTraces[index].status = "unresolved"; finishJob(slot);}
          };
          worker.onerror = () => {
            if (id !== cuspVersion || slot.worker !== worker || !slot.job) return;
            cuspTraces[slot.job.index].status = "unavailable";
            closeSlot(slot); finishJob(slot);
          };
        }
        slot.worker.postMessage({type: "trace", id: jobId, counts, path: {kind: "edge", index: edge, position: pathPosition}, startEpsilon: ".04", target, samplesPerDecade: 2, tailRefinement});
      } catch {closeSlot(slot); cuspTraces[index].status = "unavailable"; finishJob(slot);}
    }
    slots.forEach(next);
  }
  function showSource() {
    const [x, y] = sourcePosition(point.r, point.s);
    for (const id of ["source-dot", "source-halo"]) {
      get(id).setAttribute("cx", x); get(id).setAttribute("cy", y);
    }
    get("source-coordinates").textContent = `r = ${point.r.toFixed(6)}, s = ${point.s.toFixed(6)}`;
    get("source-r").value = String(point.r); get("source-s").value = String(point.s);
    source.setAttribute("aria-label", `Monodromy point r ${point.r.toFixed(6)}, s ${point.s.toFixed(6)}. Drag inside the shaded compact region or use arrow keys.`);
    const api = window.AtlasView, mapped = api.monodromyMap(point.r, point.s, selectedCounts);
    attachedSource.setAttribute("aria-label", `Attached monodromy point R ${mapped.r.toFixed(6)}, S ${mapped.s.toFixed(6)} in triangle Δ${mapped.index}. Drag inside its shaded compact region or use arrow keys to move both points.`);
    const compact = [[.48, .04], [.48, .48], [.04, .48]];
    get("source-full").setAttribute("points", api.monodromyVertices([0, 0, 0]).map(p => api.sourceAtlasXY(...p).join(",")).join(" "));
    get("source-compact").setAttribute("points", compact.map(p => api.sourceAtlasXY(...p).join(",")).join(" "));
    get("attached-source-full").setAttribute("points", api.monodromyVertices(selectedCounts).map(p => api.sourceAtlasXY(...p).join(",")).join(" "));
    get("attached-source-compact").setAttribute("points", compact.map(([r, s]) => {const p = api.monodromyMap(r, s, selectedCounts); return api.sourceAtlasXY(p.r, p.s).join(",");}).join(" "));
    const [ax, ay] = api.sourceAtlasXY(mapped.r, mapped.s);
    for (const name of ["attached-source-dot", "attached-source-halo"]) {get(name).setAttribute("cx", ax); get(name).setAttribute("cy", ay);}
    get("attached-source-label").textContent = `Attachment (${selectedCounts.join(", ")}) · Δ${mapped.index}`;
    get("attached-source-coordinates").textContent = `R = ${mapped.r.toFixed(6)}, S = ${mapped.s.toFixed(6)} · half shifts (${mapped.shiftR}, ${mapped.shiftS})`;
  }
  function saveSelection() {
    const url = new URL(location.href);
    countNames.forEach((name, i) => url.searchParams.set(name, selectedCounts[i]));
    url.searchParams.set("r", point.r.toPrecision(15)); url.searchParams.set("s", point.s.toPrecision(15));
    url.searchParams.set("view", get("explore-view").value);
    history.replaceState(null, "", url);
    refreshSharedInformation();
  }
  function clearPointReadouts(text = "Computing…") {
    pointResult = null;
    for (const name of ["base", "result"]) {
      get(name + "-tau").textContent = text;
      get(name + "-tau").setAttribute("aria-busy", text === "Computing…" ? "true" : "false");
      get(name + "-tau").removeAttribute("title");
      get(name + "-diagnostic").textContent = "";
    }
    redraw();
  }
  function showPointResult(result) {
    pointResult = result;
    for (const name of ["base", "result"]) {
      const branch = result[name], output = get(name + "-tau");
      output.setAttribute("aria-busy", "false");
      if (branch?.status === "ok" && Number.isFinite(branch.tau?.re) && Number.isFinite(branch.tau?.im)) {
        const z = branch.tau, fmt = x => Number(x.toPrecision(10)).toString();
        output.textContent = `τ = ${fmt(z.re)} + ${fmt(z.im)}i`;
        output.title = `τ = ${z.re.toPrecision(16)} + ${z.im.toPrecision(16)}i (numerical)`;
        const estimate = branch.correctionEstimate;
        get(name + "-diagnostic").textContent = Number.isFinite(estimate) ?
          `Estimated Newton correction: ${estimate.toExponential(1)}; not a certificate.` : "Numerically resolved on the selected branch.";
      } else {
        output.textContent = "Unresolved";
        get(name + "-diagnostic").textContent = "This branch needs more precision here; no replacement root is shown.";
      }
    }
    message("point-status", result.status === "ok" ? "Both period ratios resolved at this source point." :
      "One or both branches are unresolved at this point.", result.status !== "ok");
    redraw();
  }
  function dispatchPoint() {
    clearTimeout(pointTimer); pointTimer = null;
    if (pointBusy || !queuedPoint) return;
    const job = queuedPoint; queuedPoint = null;
    try {
      if (!pointWorker) {
        pointWorker = new Worker("numerics-worker.js?v=20261021-dual1");
        pointWorker.onmessage = ({data}) => {
          if (!["point", "pointError"].includes(data.type)) return;
          if (data.id !== activePoint?.id) return;
          const sent = activePoint;
          pointBusy = false; activePoint = null;
          if (data.id === pointVersion && sent?.id === data.id) {
            const result = data.result;
            const matches = result && Array.isArray(result.counts) &&
              result.counts.every((v, i) => v === selectedCounts[i]) &&
              Math.abs(result.r - point.r) < 1e-10 && Math.abs(result.s - point.s) < 1e-10;
            if (data.type === "point" && matches) showPointResult(result);
            else {
              clearPointReadouts("Unresolved");
              message("point-status", data.message ?? "The point result did not match the current selection.", true);
            }
          }
          dispatchPoint();
        };
        pointWorker.onerror = () => {
          const failed = activePoint;
          pointWorker?.terminate(); pointWorker = null; pointBusy = false; activePoint = null;
          if (failed?.id === pointVersion) {
            clearPointReadouts("Unavailable");
            message("point-status", "The point solver could not run. Move the point to retry, or reload.", true);
          }
          dispatchPoint();
        };
      }
      pointBusy = true; activePoint = job;
      pointWorker.postMessage(job);
    } catch (error) {
      pointBusy = false; activePoint = null;
      clearPointReadouts("Unavailable"); message("point-status", error.message, true);
    }
  }
  function requestPoint() {
    const id = ++pointVersion;
    clearPointReadouts();
    queuedPoint = {type: "solve", id, counts: [...selectedCounts], r: point.r, s: point.s, inset: 0.04};
    message("point-status", "Solving at the selected point…");
    if (pointTimer === null) pointTimer = setTimeout(dispatchPoint, 35);
  }
  function movePoint(r, s, save = false) {
    const next = projectCompact(r, s);
    if (Math.abs(next.r - point.r) + Math.abs(next.s - point.s) > 1e-14) {
      point = next; showSource(); requestPoint();
    }
    if (save) saveSelection();
  }
  function pointerPoint(event, svg) {
    const rect = svg.getBoundingClientRect();
    const [, , width, height] = svg.getAttribute("viewBox").trim().split(/\s+/).map(Number);
    const scale = Math.min(rect.width / width, rect.height / height);
    if (!(scale > 0)) return;
    const x = (event.clientX - rect.left - (rect.width - scale * width) / 2) / scale;
    const y = (event.clientY - rect.top - (rect.height - scale * height) / 2) / scale;
    const api = window.AtlasView;
    let coordinates = api.sourceAtlasRS(x, y);
    if (svg === attachedSource) coordinates = api.monodromyInverse(coordinates.r, coordinates.s, selectedCounts);
    movePoint(coordinates.r, coordinates.s);
  }
  function stopMesh() {
    meshWorker?.terminate(); meshWorker = null;
    get("compute").disabled = false; get("cancel").disabled = true;
  }
  function acceptMesh(mesh, elapsed, cached = false) {
    currentMesh = mesh;
    const count = name => mesh.samples.filter(sample => sample?.[name]?.status === "ok").length;
    message("attachment-status", `Compact mesh: base ${count("base")}/${mesh.points.length}; result ${count("result")}/${mesh.points.length} samples. ` +
      (cached ? "Reused from this tab." : `${elapsed.toFixed(2)} seconds.`) + (mesh.status === "ok" ? "" : " Unresolved faces omitted."));
    redraw();
  }
  function compute(event) {
    event?.preventDefault();
    if (!get("attachment-form").reportValidity()) return;
    const counts = countNames.map(name => get(name).valueAsNumber);
    if (!validCounts(counts) || order(counts) > 12) {
      message("attachment-status", "Enter nonnegative whole numbers with a + b + c ≤ 11.", true); return;
    }
    stopMesh(); const id = ++meshVersion;
    selectedCounts = counts;
    showSource();
    saveSelection();
    get("result-heading").textContent = `Attachment (${counts.join(", ")}): n = ${order(counts)}`;
    currentMesh = null; requestPoint();
    startCuspApproaches();
    const resolution = Number(get("quality").value), key = `${counts.join(",")}/${resolution}/0.04`;
    if (meshCache.has(key)) { acceptMesh(meshCache.get(key), 0, true); return; }
    get("compute").disabled = true; get("cancel").disabled = false;
    message("attachment-status", "Sampling the compact interior image…");
    const started = performance.now();
    try {
      meshWorker = new Worker("numerics-worker.js?v=20261021-dual1");
      meshWorker.onmessage = ({data}) => {
        if (data.id !== meshVersion) return;
        if (data.type === "progress") message("attachment-status", `Compact mesh samples: ${data.completed}/${data.total}…`);
        if (data.type === "done") {
          stopMesh(); meshCache.set(key, data.mesh);
          if (meshCache.size > 8) meshCache.delete(meshCache.keys().next().value);
          acceptMesh(data.mesh, (performance.now() - started) / 1000);
        }
        if (data.type === "error") { stopMesh(); message("attachment-status", data.message, true); }
      };
      meshWorker.onerror = () => {
        if (id !== meshVersion) return;
        stopMesh(); message("attachment-status", "The interior mesh could not run. Point solving is independent.", true);
      };
      meshWorker.postMessage({type: "compute", id, counts, resolution, inset: 0.04});
    } catch (error) { stopMesh(); message("attachment-status", error.message, true); }
  }

  get("attachment-form").addEventListener("submit", compute);
  get("cusp-controls").addEventListener("submit", startCuspApproaches);
  get("cusp-cancel").addEventListener("click", () => stopCuspApproaches());
  get("explore-view").addEventListener("change", () => {saveSelection(); redraw();});
  get("cancel").addEventListener("click", () => {
    ++meshVersion; stopMesh(); message("attachment-status", "Mesh stopped. The draggable point remains active.");
  });
  get("point-form").addEventListener("submit", event => {
    event.preventDefault();
    if (!get("point-form").reportValidity()) return;
    movePoint(get("source-r").valueAsNumber, get("source-s").valueAsNumber, true);
    dispatchPoint();
  });
  function changeNumber(id, direction, large = false) {
    const input = get(id), current = input.valueAsNumber;
    if (!Number.isFinite(current)) return;
    const step = (id === "cusp-position" ? .01 : .002) * (large ? 5 : 1);
    const next = Math.max(Number(input.min), Math.min(Number(input.max),
      Number((current + direction * step).toPrecision(15))));
    if (id === "cusp-position") {
      input.value = String(next);
      stopCuspApproaches("Paired positions changed. Recompute cusp series to apply them; existing computed samples are retained.");
    } else {
      const r = id === "source-r" ? next : get("source-r").valueAsNumber;
      const s = id === "source-s" ? next : get("source-s").valueAsNumber;
      if (![r, s].every(Number.isFinite)) return;
      movePoint(r, s, true); dispatchPoint();
    }
  }
  for (const id of ["source-r", "source-s", "cusp-position"]) {
    for (const [suffix, direction] of [["up", 1], ["down", -1]]) {
      get(id + "-" + suffix).addEventListener("click", event => changeNumber(id, direction, event.shiftKey));
    }
    get(id).addEventListener("keydown", event => {
      if (!["ArrowUp", "ArrowDown"].includes(event.key)) return;
      event.preventDefault(); changeNumber(id, event.key === "ArrowUp" ? 1 : -1, event.shiftKey);
    });
  }
  get("cusp-position").addEventListener("input", () => {
    stopCuspApproaches("Paired positions changed. Recompute cusp series to apply them; existing computed samples are retained.");
  });
  for (const svg of [source, attachedSource]) {
    svg.addEventListener("pointerdown", event => {
      if (event.button != null && event.button !== 0) return;
      pointerId = event.pointerId; pointerSource = svg;
      svg.setPointerCapture(pointerId); pointerPoint(event, svg);
    });
    svg.addEventListener("pointermove", event => {
      if (pointerSource === svg && pointerId === event.pointerId) pointerPoint(event, svg);
    });
    function finishPointer(event) {
      if (pointerSource !== svg || pointerId !== event.pointerId) return;
      pointerId = null; pointerSource = null; dispatchPoint(); saveSelection();
    }
    svg.addEventListener("pointerup", finishPointer);
    svg.addEventListener("pointercancel", finishPointer);
    svg.addEventListener("lostpointercapture", () => {
      if (pointerSource === svg) { pointerId = null; pointerSource = null; }
    });
    svg.addEventListener("keydown", event => {
      const step = event.shiftKey ? 0.01 : 0.002;
      const delta = {ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, step], ArrowDown: [0, -step]}[event.key];
      if (!delta) return;
      event.preventDefault();
      if (svg === attachedSource) {
        const api = window.AtlasView, mapped = api.monodromyMap(point.r, point.s, selectedCounts);
        const next = api.monodromyInverse(mapped.r + delta[0], mapped.s + delta[1], selectedCounts);
        movePoint(next.r, next.s, true);
      } else movePoint(point.r + delta[0], point.s + delta[1], true);
    });
  }
  window.addEventListener("resize", redraw);
  window.addEventListener("pagehide", () => {
    clearTimeout(pointTimer); pointWorker?.terminate(); meshWorker?.terminate(); ++cuspVersion; cuspWorkers.forEach(worker => worker.terminate()); cuspWorkers = [];
  });
  showSource(); compute();
}


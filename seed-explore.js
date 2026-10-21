"use strict";

// New seed families have their own marked character domains. Keep the Lamé
// controller separate: its source triangle and Γ(2) parking are different.
(() => {
  const api = window.SeedCases;
  if (!api?.isNew || !document.getElementById("monodromy-source")) return;
  const active = api.active, get = id => document.getElementById(id);
  const params = new URLSearchParams(location.search), inset = .04;
  const colors = {base: "#285f65", result: "#b55e3c"};
  const meshCache = new Map(), traceCache = new Map();
  let point = api.project(params.has("r") ? Number(params.get("r")) : active.start[0],
    params.has("s") ? Number(params.get("s")) : active.start[1], inset);
  let pointResult = null, mesh = null, traces = Array(6).fill(null);
  let pointWorker = null, pointBusy = false, pointQueued = null, pointVersion = 0, pointTimer = null;
  let meshWorker = null, meshVersion = 0, traceVersion = 0, traceWorkers = [];
  let pointer = null, redrawFrame = null;
  const fmt = value => Number.isFinite(Number(value)) ? Number(Number(value).toPrecision(9)).toString() : "unresolved";
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
  const tau = branch => branch?.status === "ok" && window.AtlasView.finiteTau(branch.tau) ? branch.tau : null;
  const nodeTau = (node, field) => tau(node?.[field]);
  const makeWorker = () => new Worker("repeated-worker.mjs?v=20261021-repeats4", {type: "module"});
  const caseRequest = () => ({caseId: active.id, caseConfig: active.caseConfig});

  function message(id, text, error = false) {
    get(id).textContent = text;
    get(id).classList.toggle("error", error);
  }
  function applyLabels() {
    for (const id of ["a", "b", "c"]) get(id).closest("div").hidden = true;
    get("compute").textContent = "Compute compact image";
    get("quality").innerHTML = active.isRepeated ? '<option value="4">Quick — 15 points</option><option value="8">Standard — 45 points</option><option value="12">Fine — 91 points</option>' : '<option value="20">Precomputed — 231 points</option>';
    get("quality").value = active.isRepeated ? (["4", "8", "12"].includes(params.get("quality")) ? params.get("quality") : "8") : "20";
    get("attachment-caption").textContent = active.isRepeated ? "Apply the labeled counts above, then drag either source point. Interior and cusp computations run independently; increasing the counts or mesh size can take longer. Only resolved samples are drawn." : "The compact seed and single-decoration images are precomputed. Dragging and cusp series use the selected geometric branch; the cusp solver increases precision near the boundary.";
    get("base-source-label").textContent = `Base: ${active.baseTitle}`;
    get("attached-source-label").textContent = active.operation === "seed" ? "Same seed · no decoration" : `After decoration: ${active.title}`;
    get("base-heading").textContent = `Base: ${active.baseTitle}`;
    get("result-heading").textContent = active.operation === "seed" ? "Same seed · no decoration" : active.title;
    get("tau-workspace-heading").textContent = "Images in the marked τ chart";
    get("explore-chart-caption").textContent = "This is a marked upper-half-plane chart with Γ(2) tile guides; the points are not moved into a single fundamental domain. Height is compressed as Im τ / (1 + Im τ). The original τ values appear above. Unresolved samples are omitted; numerical checks are not interval certificates.";
    get("tau-image-caption").textContent = "Shading shows only the computed compact interior. Cusp curves contain numerical samples, stop at the last resolved point, and use no forced endpoint. The raw cusp at ½, when present, is retained in this marking even though it is Γ(2)-equivalent to ∞.";
    get("explore-sampling-note").textContent = "Dark shading samples the compact barycentric interior Kδ with δ = 0.04; its dashed contour is an artificial cutoff. Two numerical series approach each open source edge. Their light ribbon interpolates between computed paths and is not a certified image boundary. Source vertices have direction-dependent finite boundary limits.";
    get("source-domain-note").textContent = "Geometric characters: χ(1) = exp(−2πiR), χ(τ) = exp(2πiS). The translated triangle is an unwrapped lift of ordinary monodromy. Drag either point: the half-period shift is inverted automatically. Every decoration preserves projective monodromy.";
    document.querySelector('label[for="source-r"]').textContent = "Base R";
    document.querySelector('label[for="source-s"]').textContent = "Base S";
    get("source-r").min = "0"; get("source-r").max = "1";
    get("source-s").min = "0"; get("source-s").max = String(active.top);
    get("source-r").step = "any"; get("source-s").step = "any";
    const details = get("cusp-controls").closest("details");
    details.querySelector("p.caption").textContent = "Two numerical approaches per open edge, at paired boundary positions t and 1 − t. Samples use the selected branch and extra precision near the cusp. Results are cached in this tab; point dragging stays independent.";
    get("source-region-details").innerHTML = `<summary>Source region and cusp correspondence</summary><p class="caption">The dashed full triangle is excluded. Dragging stays inside its compact barycentric inset, with every barycentric coordinate at least ${inset}. Coordinates use the undeformed seed triangle.</p><div class="edge-cusp-key"><strong>Open source edges → raw marked cusps</strong>${active.edgeLabels.map(label => `<span>${escape(label)}</span>`).join("")}</div><p class="caption">The limit along an open edge is a cusp. Approaching a source vertex can instead give a finite limit depending on direction.</p>`;
    get("explore-view").value = params.get("view") === "overlay" ? "overlay" : "separate";
    const p = Number(params.get("cuspT"));
    get("cusp-position").value = p >= .02 && p <= .48 ? String(p) : "0.2";
    get("cusp-depth").innerHTML = '<option value="8">10⁻⁸</option><option value="10">10⁻¹⁰</option><option value="20">10⁻²⁰</option><option value="40">10⁻⁴⁰</option><option value="80">10⁻⁸⁰</option>';
    get("cusp-depth").value = ["8", "10", "20", "40", "80"].includes(params.get("cuspDepth")) ? params.get("cuspDepth") : "8";
  }
  function boundaryPoint(edge, position) {
    if (edge === 0) return [1, active.top * position];
    if (edge === 1) return [position, active.top];
    return [position, active.top * (1 - position)];
  }
  function save() {
    const url = new URL(location.href);
    url.searchParams.set("r", point.r.toPrecision(15));
    url.searchParams.set("s", point.s.toPrecision(15));
    url.searchParams.set("view", get("explore-view").value);
    url.searchParams.set("cuspT", get("cusp-position").value);
    url.searchParams.set("cuspDepth", get("cusp-depth").value);
    url.searchParams.set("quality", get("quality").value);
    history.replaceState(null, "", url);
    if (typeof refreshSharedInformation === "function") refreshSharedInformation();
  }
  function showSource() {
    api.renderMonodromy(get("monodromy-source"), {point: [point.r, point.s], attached: false, title: active.baseTitle});
    api.renderMonodromy(get("monodromy-attached"), {point: [point.r, point.s], attached: true, title: active.title});
    get("source-r").value = String(point.r); get("source-s").value = String(point.s);
    get("source-coordinates").textContent = `Base R = ${fmt(point.r)}, S = ${fmt(point.s)}`;
    const mapped = api.map(point.r, point.s, true);
    get("attached-source-coordinates").textContent = `After decoration: R′ = ${fmt(mapped.r)}, S′ = ${fmt(mapped.s)} · ordinary shift (${active.shift.join(", ")})`;
    get("monodromy-source").setAttribute("aria-label", `${active.baseTitle}. R ${fmt(point.r)}, S ${fmt(point.s)}. Drag the point or use arrow keys.`);
    get("monodromy-attached").setAttribute("aria-label", `${active.title}. R ${fmt(mapped.r)}, S ${fmt(mapped.s)}. Drag the point or use arrow keys.`);
  }
  function scheduleDraw() {
    if (redrawFrame === null) redrawFrame = requestAnimationFrame(() => {redrawFrame = null; redraw();});
  }
  function boundaryPaths(field) {
    const boundaries = mesh?.boundaries;
    if (!boundaries) return [];
    const list = Array.isArray(boundaries) ? boundaries : (boundaries[field] ?? []);
    return list.map(arc => ({color: colors[field], width: 1.3, radius: 0,
      points: (arc.nodes ?? arc.points ?? arc).map(node => {
        const z = Array.isArray(node) ? {re: Number(node[0]), im: Number(node[1])} : node;
        return nodeTau(node, field) ?? (window.AtlasView.finiteTau(z) ? z : null);
      })}));
  }
  function redraw() {
    const overlap = get("explore-view").value === "overlay";
    get("separate-plots").hidden = overlap; get("overlap-plot").hidden = !overlap;
    const layers = field => ({mesh, field, color: colors[field]});
    const paths = field => [...boundaryPaths(field), ...traces.filter(Boolean).map(trace => ({color: colors[field], width: 1.35, radius: 1.5,
      points: trace.nodes.map(node => nodeTau(node, field))}))];
    const bands = field => [0, 1, 2].flatMap(edge => {
      const left = traces[2 * edge], right = traces[2 * edge + 1];
      if (!left || !right) return [];
      const paired = new Map(right.nodes.map(node => [String(node.epsilon), node]));
      return [{color: colors[field], left: left.nodes.map(node => nodeTau(node, field)), right: left.nodes.map(node => nodeTau(paired.get(String(node.epsilon)), field))}];
    });
    const base = tau(pointResult?.base), result = tau(pointResult?.result);
    const settings = {height: 520, markedStrip: true};
    if (overlap) {
      const chart = window.AtlasView.drawDomain(get("overlay-image"), {...settings, meshLayers: [layers("base"), layers("result")],
        paths: [...paths("base"), ...paths("result")], bands: [...bands("base"), ...bands("result")],
        points: [{tau: base, color: colors.base}, {tau: result, color: colors.result}].filter(p => p.tau),
        arrow: base && result ? {from: base, to: result} : null});
      const distance = chart && base && result ? Math.hypot(...chart.map(base).map((value, index) => value - chart.map(result)[index])) : 0;
      get("overlap-note").textContent = !base || !result ? "The comparison arrow appears after both points resolve." : distance < 12 ? "The two points coincide or are too close to distinguish at this scale." : "The solid arrow compares base τ and decorated τ at the same source character. It is not a continuation path.";
    } else {
      for (const field of ["base", "result"]) window.AtlasView.drawDomain(get(field + "-image"), {...settings, meshLayers: [layers(field)],
        paths: paths(field), bands: bands(field), points: [{tau: field === "base" ? base : result, color: colors[field]}].filter(p => p.tau)});
    }
  }
  function showPoint(result) {
    pointResult = result;
    for (const field of ["base", "result"]) {
      const z = tau(result?.[field]), output = get(field + "-tau");
      output.setAttribute("aria-busy", "false");
      output.textContent = z ? `τ = ${fmt(z.re)} + ${fmt(z.im)}i` : "Unresolved";
      if (z) output.title = `τ = ${z.re.toPrecision(16)} + ${z.im.toPrecision(16)}i (numerical)`;
      else output.removeAttribute("title");
      get(field + "-diagnostic").textContent = z ? "Resolved on the selected geometric branch. Numerical result; not an interval certificate." : "No replacement root is displayed. Move the point to retry.";
    }
    const both = tau(result?.base) && tau(result?.result);
    message("point-status", both ? "Both marked period ratios resolved." : "One or both branches are unresolved here.", !both);
    scheduleDraw();
  }
  function clearPoint(text = "Computing…") {
    pointResult = null;
    for (const field of ["base", "result"]) {get(field + "-tau").textContent = text; get(field + "-tau").setAttribute("aria-busy", String(text === "Computing…")); get(field + "-diagnostic").textContent = "";}
    scheduleDraw();
  }
  function dispatchPoint() {
    clearTimeout(pointTimer); pointTimer = null;
    if (pointBusy || !pointQueued) return;
    const job = pointQueued; pointQueued = null;
    if (!pointWorker) {
      pointWorker = makeWorker();
      pointWorker.onmessage = ({data}) => {
        if (data.type === "progress") return;
        pointBusy = false;
        if (data.id === pointVersion) {
          if (data.type === "point" && data.result) showPoint(data.result);
          else {clearPoint("Unresolved"); message("point-status", data.message ?? data.error ?? "The selected branch could not be resolved.", true);}
        }
        dispatchPoint();
      };
      pointWorker.onerror = () => {
        pointWorker?.terminate(); pointWorker = null; pointBusy = false;
        clearPoint("Unavailable"); message("point-status", "The point solver could not run. Drag to retry.", true); dispatchPoint();
      };
    }
    pointBusy = true; pointWorker.postMessage(job);
  }
  function requestPoint() {
    clearPoint();
    pointQueued = {type: "point", id: ++pointVersion, ...caseRequest(), R: point.r, S: point.s};
    message("point-status", "Resolving the selected point…");
    if (pointTimer === null) pointTimer = setTimeout(dispatchPoint, 30);
  }
  function move(r, s, persist = false) {
    const next = api.project(r, s, inset);
    if (Math.abs(next.r - point.r) + Math.abs(next.s - point.s) > 1e-14) {point = next; showSource(); requestPoint();}
    if (persist) save();
  }
  function stopMesh(text) {
    meshWorker?.terminate(); meshWorker = null; ++meshVersion;
    get("compute").disabled = false; get("cancel").disabled = true;
    if (text) message("attachment-status", text);
  }
  function acceptMesh(result, elapsed = 0, cached = false) {
    mesh = result;
    const total = result.points.length, count = field => result.samples.filter(sample => tau(sample?.[field])).length;
    message("attachment-status", `Compact image: base ${count("base")}/${total}, decorated ${count("result")}/${total} samples. ${cached ? "Reused from this tab." : result.cached ? "Loaded the verified precomputed mesh." : `${elapsed.toFixed(2)} seconds.`} Dashed contours are sampling cutoffs.`);
    scheduleDraw();
  }
  function computeMesh(event) {
    event?.preventDefault(); stopMesh(); const id = meshVersion;
    const resolution = Number(get("quality").value), key = `${active.id}/${resolution}/${inset}`;
    save();
    if (meshCache.has(key)) {acceptMesh(meshCache.get(key), 0, true); return;}
    get("compute").disabled = true; get("cancel").disabled = false;
    message("attachment-status", `Computing ${active.title} compact image…`);
    const started = performance.now();
    meshWorker = makeWorker();
    meshWorker.onmessage = ({data}) => {
      if (data.id !== meshVersion) return;
      if (data.type === "progress") {const p = data.result ?? data; message("attachment-status", `Compact image: ${p.completed}/${p.total} samples…`);}
      else if (data.type === "mesh" && data.result) {const result = data.result; stopMesh(); meshCache.set(key, result); acceptMesh(result, (performance.now() - started) / 1000);}
      else {stopMesh(); message("attachment-status", data.message ?? data.error ?? "Compact sampling was interrupted. Point dragging is independent.", true);}
    };
    meshWorker.onerror = () => {if (id === meshVersion) {stopMesh(); message("attachment-status", "Compact sampling could not run. Point dragging is independent.", true);}};
    meshWorker.postMessage({type: "mesh", id, ...caseRequest(), R: point.r, S: point.s, resolution, inset});
  }
  function updateEndpoints() {
    get("cusp-endpoints").textContent = traces.map((trace, index) => {
      const edge = Math.floor(index / 2);
      if (!trace) return `${active.edgeLabels[edge]}, series ${index % 2 + 1}: not computed.`;
      const describe = field => {
        const node = trace.nodes.findLast(node => nodeTau(node, field)), z = nodeTau(node, field);
        return z ? `${field === "base" ? "Base" : "Decorated"}: Re τ = ${fmt(z.re)}, Im τ = ${fmt(z.im)} at path ε = ${node.epsilon}` : `${field}: unresolved`;
      };
      return `${active.edgeLabels[edge]}, t = ${trace.position} · ${trace.status}\n${describe("base")}\n${describe("result")}`;
    }).join("\n\n");
  }
  function showAsymptotics() {
    const element = get("explore-case-asymptotics"), position = Number(get("cusp-position").value);
    element.hidden = false;
    const describe = (edge, base) => {
      const a = api.asymptotic(edge, ...boundaryPoint(edge, position), base);
      return a && Number.isFinite(a.alpha) ? `cusp ${a.cusp}; Re w → ${a.xLabel ?? fmt(a.x0)}; Im w ∼ (${fmt(a.alpha)}/π) log(${fmt(a.C)}/δ)` : "Computed branch samples; no leading formula is asserted for this labeled decoration.";
    };
    element.innerHTML = `<details><summary>Cusp asymptotics for this case</summary><p class="caption">At boundary position t = ${escape(position)}, away from the vertices. The trace parameter ε need not equal the edge distance δ.</p><div class="table-scroll"><table><thead><tr><th>Source edge and coordinate</th><th>Base seed</th><th>Selected decoration</th></tr></thead><tbody>${active.edgeLabels.map((label, edge) => {const a = api.asymptotic(edge, ...boundaryPoint(edge, position), false) ?? api.asymptotic(edge, ...boundaryPoint(edge, position), true); return `<tr><th>${escape(label)}${a ? `<br><small>${escape(a.chartLabel)}; δ = ${escape(a.deltaLabel)}</small>` : ""}</th><td>${escape(describe(edge, true))}</td><td>${escape(describe(edge, false))}</td></tr>`;}).join("")}</tbody></table></div><p class="caption">Limits &amp; Zoom shows the resolved horizontal and vertical behavior. Formulas are displayed only for calibrated families; other labeled graphs are followed numerically.</p></details>`;
  }
  function stopTraces(text) {
    ++traceVersion; traceWorkers.forEach(worker => worker.terminate()); traceWorkers = [];
    traces.forEach(trace => {if (trace?.status === "computing") trace.status = "stopped";});
    get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
    updateEndpoints(); if (text) message("cusp-status", text);
  }
  function computeTraces(event) {
    event?.preventDefault();
    if (!get("cusp-controls").reportValidity()) return;
    const position = Number(get("cusp-position").value);
    if (!(position >= .02 && position <= .48)) return;
    stopTraces(); const version = traceVersion;
    const target = `1e-${get("cusp-depth").value}`;
    const key = `repeated-tails-v1/${active.id}/${point.r.toPrecision(15)}/${point.s.toPrecision(15)}/${position}/${target}`;
    save(); showAsymptotics();
    let stored = event ? null : traceCache.get(key);
    if (!stored && !event) {try {stored = JSON.parse(sessionStorage.getItem(key));} catch { /* Optional cache. */ }}
    if (Array.isArray(stored) && stored.length === 6 && stored.every(trace => Array.isArray(trace?.nodes))) {
      traces = stored; updateEndpoints(); scheduleDraw();
      message("cusp-status", "Six computed cusp series reused for the current case. No extrapolated endpoint."); return;
    }
    traces = Array(6).fill(null); updateEndpoints(); scheduleDraw();
    get("cusp-run").disabled = true; get("cusp-cancel").disabled = false;
    let nextIndex = 0, completed = 0;
    const poolSize = navigator.hardwareConcurrency > 2 ? 2 : 1;
    function finish() {
      if (completed !== 6 || version !== traceVersion) return;
      traceWorkers.forEach(worker => worker.terminate()); traceWorkers = [];
      get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
      const complete = traces.every(trace => trace?.status === "ok");
      if (traces.some(trace => trace.nodes.some(node => nodeTau(node, "base") || nodeTau(node, "result")))) {
        traceCache.set(key, structuredClone(traces));
        try {sessionStorage.setItem(key, JSON.stringify(traces));} catch { /* Optional cache. */ }
      }
      message("cusp-status", complete ? "Six cusp series computed. Curves stop at their last numerical sample." : "Cusp computations finished with partial series. Unresolved samples are omitted.");
      updateEndpoints(); scheduleDraw();
    }
    function next(worker) {
      if (version !== traceVersion || nextIndex >= 6) {finish(); return;}
      const index = nextIndex++, edge = Math.floor(index / 2), t = index % 2 ? 1 - position : position;
      const id = version * 8 + index;
      traces[index] = {nodes: [], edge, position: t, status: "computing"};
      worker.onmessage = ({data}) => {
        if (version !== traceVersion || data.id !== id) return;
        if (data.type === "progress") {
          const p = data.result ?? data;
          if (p.node) traces[index].nodes.push(p.node);
          message("cusp-status", `${active.edgeLabels[edge]}, t = ${fmt(t)}: ${p.completed}/${p.total} samples (${completed}/6 series finished).`);
          updateEndpoints(); scheduleDraw();
        } else if (data.type === "trace" && data.result) {
          traces[index] = {...traces[index], ...data.result, status: data.result.status ?? "ok"};
          completed++; updateEndpoints(); scheduleDraw(); next(worker); finish();
        } else {traces[index].status = "unresolved"; completed++; next(worker); finish();}
      };
      worker.onerror = () => {
        if (version !== traceVersion) return;
        traces[index].status = "unavailable"; completed++; worker.terminate();
        traceWorkers = traceWorkers.filter(w => w !== worker);
        if (nextIndex < 6) {const replacement = makeWorker(); traceWorkers.push(replacement); next(replacement);} finish();
      };
      worker.postMessage({type: "trace", id, ...caseRequest(), R: point.r, S: point.s, path: {kind: "edge", index: edge, position: t}, startEpsilon: .04, target, samplesPerDecade: 2});
    }
    for (let i = 0; i < poolSize; i++) {const worker = makeWorker(); traceWorkers.push(worker); next(worker);}
  }
  function changeNumber(id, direction, large = false) {
    const input = get(id), current = input.valueAsNumber;
    if (!Number.isFinite(current)) return;
    const delta = (id === "cusp-position" ? .01 : .002) * (large ? 5 : 1);
    const value = Math.max(Number(input.min), Math.min(Number(input.max), Number((current + direction * delta).toPrecision(15))));
    if (id === "cusp-position") {input.value = String(value); stopTraces("Paired positions changed. Recompute to apply them; the existing numerical curves remain visible."); showAsymptotics();}
    else {move(id === "source-r" ? value : point.r, id === "source-s" ? value : point.s, true); dispatchPoint();}
  }
  function pointerPoint(event, svg, attached) {
    const rect = svg.getBoundingClientRect(), box = svg.viewBox.baseVal;
    const scale = Math.min(rect.width / box.width, rect.height / box.height);
    if (!(scale > 0)) return;
    const x = box.x + (event.clientX - rect.left - (rect.width - scale * box.width) / 2) / scale;
    const y = box.y + (event.clientY - rect.top - (rect.height - scale * box.height) / 2) / scale;
    const next = api.svgInverse(x, y, attached); move(next.r, next.s);
  }

  applyLabels(); showSource(); showAsymptotics();
  get("attachment-form").addEventListener("submit", computeMesh);
  get("cancel").addEventListener("click", () => stopMesh("Compact image computation stopped. Point dragging stays active."));
  get("point-form").addEventListener("submit", event => {event.preventDefault(); if (get("point-form").reportValidity()) {move(get("source-r").valueAsNumber, get("source-s").valueAsNumber, true); dispatchPoint();}});
  get("explore-view").addEventListener("change", () => {save(); scheduleDraw();});
  get("cusp-controls").addEventListener("submit", computeTraces);
  get("cusp-cancel").addEventListener("click", () => stopTraces("Stopped. Computed samples remain visible."));
  get("cusp-position").addEventListener("input", () => {stopTraces("Paired positions changed. Recompute to apply them."); showAsymptotics();});
  for (const id of ["source-r", "source-s", "cusp-position"]) {
    for (const [suffix, direction] of [["up", 1], ["down", -1]]) get(`${id}-${suffix}`).addEventListener("click", event => changeNumber(id, direction, event.shiftKey));
    get(id).addEventListener("keydown", event => {if (["ArrowUp", "ArrowDown"].includes(event.key)) {event.preventDefault(); changeNumber(id, event.key === "ArrowUp" ? 1 : -1, event.shiftKey);}});
  }
  for (const [id, attached] of [["monodromy-source", false], ["monodromy-attached", true]]) {
    const svg = get(id);
    svg.addEventListener("pointerdown", event => {if (event.button != null && event.button !== 0) return; pointer = {svg, id: event.pointerId}; svg.setPointerCapture(event.pointerId); pointerPoint(event, svg, attached);});
    svg.addEventListener("pointermove", event => {if (pointer?.svg === svg && pointer.id === event.pointerId) pointerPoint(event, svg, attached);});
    const finish = event => {if (pointer?.svg === svg && pointer.id === event.pointerId) {pointer = null; dispatchPoint(); save();}};
    svg.addEventListener("pointerup", finish); svg.addEventListener("pointercancel", finish);
    svg.addEventListener("lostpointercapture", () => {if (pointer?.svg === svg) pointer = null;});
    svg.addEventListener("keydown", event => {
      const step = event.shiftKey ? .01 : .002, d = {ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, step], ArrowDown: [0, -step]}[event.key];
      if (d) {event.preventDefault(); move(point.r + d[0], point.s + d[1], true);}
    });
  }
  window.addEventListener("resize", scheduleDraw);
  window.addEventListener("pagehide", () => {clearTimeout(pointTimer); pointWorker?.terminate(); stopMesh(); stopTraces(); if (redrawFrame !== null) cancelAnimationFrame(redrawFrame);});
  requestPoint(); computeMesh(); computeTraces();
})();

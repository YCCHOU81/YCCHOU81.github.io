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

const sourcePosition = (r, s) => {
  const weights = [1 - 2 * s, 2 * (r + s) - 1, 1 - 2 * r];
  return [48 * weights[0] + 200 * weights[1] + 352 * weights[2], 272 - 238 * weights[1]];
};

refreshSharedInformation();
if (document.getElementById("attachment-form")) connectExplorer();

function connectExplorer() {
  const get = id => document.getElementById(id);
  const source = get("monodromy-source");
  const meshCache = new Map();
  let meshWorker = null, meshVersion = 0, currentMesh = null;
  let pointWorker = null, pointVersion = 0, pointBusy = false, queuedPoint = null;
  let activePoint = null, pointTimer = null, pointResult = null, pointerId = null;
  let point = projectCompact(query.has("r") ? Number(query.get("r")) : 0.31,
    query.has("s") ? Number(query.get("s")) : 0.29);
  countNames.forEach((name, i) => { get(name).value = selectedCounts[i]; });

  function message(id, text, error = false) {
    get(id).textContent = text;
    get(id).classList.toggle("error", error);
  }
  function redraw() {
    drawImage(get("base-image"), currentMesh, "base", "#285f65", pointResult?.base);
    drawImage(get("result-image"), currentMesh, "result", "#b55e3c", pointResult?.result);
  }
  function showSource() {
    const [x, y] = sourcePosition(point.r, point.s);
    for (const id of ["source-dot", "source-halo"]) {
      get(id).setAttribute("cx", x); get(id).setAttribute("cy", y);
    }
    get("source-coordinates").textContent = `r = ${point.r.toFixed(6)}, s = ${point.s.toFixed(6)}`;
    get("source-r").value = point.r.toFixed(6); get("source-s").value = point.s.toFixed(6);
    source.setAttribute("aria-label", `Monodromy point r ${point.r.toFixed(6)}, s ${point.s.toFixed(6)}. Drag inside the shaded compact region or use arrow keys.`);
  }
  function saveSelection() {
    const url = new URL(location.href);
    countNames.forEach((name, i) => url.searchParams.set(name, selectedCounts[i]));
    url.searchParams.set("r", point.r.toPrecision(15)); url.searchParams.set("s", point.s.toPrecision(15));
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
        pointWorker = new Worker("numerics-worker.js");
        pointWorker.onmessage = ({data}) => {
          if (!["point", "pointError"].includes(data.type)) return;
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
  function pointerPoint(event) {
    const rect = source.getBoundingClientRect();
    const scale = Math.min(rect.width / 400, rect.height / 320);
    if (!(scale > 0)) return;
    const x = (event.clientX - rect.left - (rect.width - scale * 400) / 2) / scale;
    const y = (event.clientY - rect.top - (rect.height - scale * 320) / 2) / scale;
    const b1 = (272 - y) / 238, b2 = (x - 48 - 152 * b1) / 304, b0 = 1 - b1 - b2;
    movePoint((1 - b2) / 2, (1 - b0) / 2);
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
    saveSelection();
    get("result-heading").textContent = `Attachment (${counts.join(", ")}): n = ${order(counts)}`;
    currentMesh = null; requestPoint();
    const resolution = Number(get("quality").value), key = `${counts.join(",")}/${resolution}/0.04`;
    if (meshCache.has(key)) { acceptMesh(meshCache.get(key), 0, true); return; }
    get("compute").disabled = true; get("cancel").disabled = false;
    message("attachment-status", "Sampling the compact interior image…");
    const started = performance.now();
    try {
      meshWorker = new Worker("numerics-worker.js");
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
  get("cancel").addEventListener("click", () => {
    ++meshVersion; stopMesh(); message("attachment-status", "Mesh stopped. The draggable point remains active.");
  });
  get("point-form").addEventListener("submit", event => {
    event.preventDefault();
    if (!get("point-form").reportValidity()) return;
    movePoint(get("source-r").valueAsNumber, get("source-s").valueAsNumber, true);
    dispatchPoint();
  });
  source.addEventListener("pointerdown", event => {
    if (event.button != null && event.button !== 0) return;
    pointerId = event.pointerId; source.setPointerCapture(pointerId); pointerPoint(event);
  });
  source.addEventListener("pointermove", event => { if (pointerId === event.pointerId) pointerPoint(event); });
  function finishPointer(event) {
    if (pointerId !== event.pointerId) return;
    pointerId = null; dispatchPoint(); saveSelection();
  }
  source.addEventListener("pointerup", finishPointer);
  source.addEventListener("pointercancel", finishPointer);
  source.addEventListener("lostpointercapture", () => { pointerId = null; });
  source.addEventListener("keydown", event => {
    const step = event.shiftKey ? 0.01 : 0.002;
    const delta = {ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, step], ArrowDown: [0, -step]}[event.key];
    if (!delta) return;
    event.preventDefault(); movePoint(point.r + delta[0], point.s + delta[1], true);
  });
  window.addEventListener("resize", redraw);
  window.addEventListener("pagehide", () => {
    clearTimeout(pointTimer); pointWorker?.terminate(); meshWorker?.terminate();
  });
  showSource(); compute();
}

// The neutral background is the coordinate chart; only mesh faces are data.
function drawImage(canvas, mesh, field, color, selected) {
  const width = canvas.getBoundingClientRect().width, height = 340;
  if (width < 1) return;
  const ratio = Math.min(devicePixelRatio || 1, 2);
  canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  const ctx = canvas.getContext("2d"); ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  const left = 42, top = 26, w = Math.max(20, width - 58), h = height - 72;
  const px = x => left + x * w, py = u => top + (1 - u) * h;
  const point = z => [px(z.re), py(z.im / (1 + z.im))];
  const line = (a, b) => { ctx.moveTo(...a); ctx.lineTo(...b); };
  const chart = new Path2D(); chart.moveTo(px(0), py(1));
  for (let i = 0; i <= 120; i++) {
    const x = i / 120; chart.lineTo(...point({re: x, im: Math.sqrt(x * (1 - x))}));
  }
  chart.lineTo(px(1), py(1)); chart.closePath();
  ctx.fillStyle = "#f3f3ee"; ctx.fill(chart);
  ctx.save(); ctx.clip(chart); ctx.strokeStyle = "#e0e3dc"; ctx.lineWidth = 0.7; ctx.beginPath();
  for (const y of [0.5, 1, 2, 4]) line(point({re: 0, im: y}), point({re: 1, im: y}));
  for (const x of [0.25, 0.5, 0.75]) line([px(x), py(0)], [px(x), py(1)]);
  ctx.stroke();
  if (mesh) {
    const pts = mesh.samples.map(sample => {
      const branch = sample?.[field];
      return branch?.status === "ok" && Number.isFinite(branch.tau?.re) && Number.isFinite(branch.tau?.im) ? point(branch.tau) : null;
    });
    const edges = new Map(); ctx.fillStyle = color + "18"; ctx.strokeStyle = color + "36"; ctx.lineWidth = 0.4;
    for (const face of mesh.triangles) {
      for (let i = 0; i < 3; i++) {
        const a = face[i], b = face[(i + 1) % 3], key = [Math.min(a, b), Math.max(a, b)].join(",");
        edges.set(key, (edges.get(key) || 0) + 1);
      }
      const p = face.map(i => pts[i]); if (p.some(v => !v)) continue;
      ctx.beginPath(); ctx.moveTo(...p[0]); ctx.lineTo(...p[1]); ctx.lineTo(...p[2]); ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    // All outer sampled edges are artificial cutoffs, drawn dashed.
    ctx.strokeStyle = color; ctx.lineWidth = 1.4; ctx.setLineDash([5, 4]); ctx.beginPath();
    for (const [key, count] of edges) {
      if (count !== 1) continue;
      const [a, b] = key.split(",").map(Number); if (pts[a] && pts[b]) line(pts[a], pts[b]);
    }
    ctx.stroke(); ctx.setLineDash([]);
  }
  if (selected?.status === "ok" && Number.isFinite(selected.tau?.re) && Number.isFinite(selected.tau?.im)) {
    const p = point(selected.tau);
    ctx.fillStyle = color + "28"; ctx.beginPath(); ctx.arc(...p, 10, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = color; ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(...p, 4.7, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
  }
  ctx.restore(); ctx.strokeStyle = "#afb6aa"; ctx.lineWidth = 0.8; ctx.stroke(chart);
  ctx.fillStyle = "#667261"; ctx.font = "10px system-ui"; ctx.textAlign = "center";
  ctx.fillText("INTERIOR SAMPLE · CUSPS OMITTED", left + w / 2, 13);
  for (const x of [0, 0.5, 1]) ctx.fillText(x, px(x), height - 24);
  ctx.fillText("Re τ", left + w / 2, height - 5); ctx.textAlign = "right";
  for (const y of [0.5, 1, 2, 4]) ctx.fillText(y, left - 6, point({re: 0, im: y})[1] + 4);
  ctx.fillText("∞", left - 6, top + 4);
  ctx.save(); ctx.translate(11, top + h / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText("Im τ", 0, 0); ctx.restore();
}

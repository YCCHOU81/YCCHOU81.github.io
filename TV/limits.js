"use strict";

const initializeLimits = () => {
  const $ = id => document.getElementById(id);

  if (!$("limits-controls")) return;

  const cases = window.SeedCases;
  const activeCase = cases?.active;
  let baseLabel = activeCase.baseTitle;
  let resultLabel = activeCase.title;
  // Persistent controls keep their listeners and display choices across families.
  const initialText = Object.fromEntries(["limit-source-controls-note","limit-depth-advice"]
    .map(id => [id,$(id)?.textContent ?? ""]));
  const familyDepths = new Map();
  let depthRevision = 0, expectedDepth = null;
  function trackDepthChoice() {
    const depth = $("limit-depth").value;
    if (expectedDepth !== null && depth !== expectedDepth) depthRevision++;
    expectedDepth = depth;
  }
  let selectedFamily = activeCase?.seed ?? window.AtlasFamily.defaultSeed;

  const color = "#b55e3c";
  const baseColor = "#285f65";
  const LENS_HEIGHT = 260;
  const LENS_PLOT_HEIGHT = 155;
  const LENS_LEFT = 100;
  const LENS_RIGHT = 16;

  const s = {
    job: 0,
    workers: {},
    pending: new Set(),
    nodes: [],
    mesh: null,
    meshMeta: null,
    step: -1,
    chart: "raw",
    fit: "tail",
    view: null,
    result: null,
    request: null,
    traceJob: null,
    progressTotal: 0,
    zoomCorners: [],
    playing: false,
    playTimer: null,
    detachedRestore: false,
    extension: null,
    draftDirty: false,
    start: ({r:activeCase.start[0],s:activeCase.start[1]})
  };

  // Every numerical cache shares one memory budget. In particular, repeated TV
  // computations cannot accumulate in an unbounded local Map or web storage.
  const cacheNamespace = kind => "limits-" + kind + "-v1";
  const DISPLAY_RESULTS = "limits-display-v1";
  function graphDescriptor(request = null) {
    return request
      ? request.caseId
        ? {family: request.seed, caseId: request.caseId, caseConfig: request.caseConfig ?? null}
        : {family: "lame", counts: request.counts}
      : ({family: activeCase.seed, caseId: activeCase.id, caseConfig: activeCase.caseConfig ?? null});
  }
  function saveDisplaySettings() {
    window.AtlasResults?.set(DISPLAY_RESULTS, "current", {
      chart:s.chart, fit:s.fit, zoom:$("limit-zoom").value,
      metric:$("limit-metric").value, compare:$("limit-compare").checked,
      resolveHorizontal:$("limit-resolve-horizontal").checked,
      quantitiesOpen:$("limits-quantities").open, customOpen:$("limit-custom-settings").open
    });
  }
  function restoreDisplaySettings() {
    const saved=window.AtlasResults?.get(DISPLAY_RESULTS,"current");
    if(!saved)return;
    // Explicit links may choose a chart, but page navigation keeps the other
    // viewing choices, including magnification and the comparison overlays.
    if(!query.has("chart") && ["raw","zero","one","half"].includes(saved.chart))s.chart=saved.chart;
    if(["trace","tail","image"].includes(saved.fit))s.fit=saved.fit;
    if(Number.isFinite(Number(saved.zoom)))$("limit-zoom").value=String(saved.zoom);
    if([...$("limit-metric").options].some(option=>option.value===saved.metric))$("limit-metric").value=saved.metric;
    if(typeof saved.compare==="boolean")$("limit-compare").checked=saved.compare;
    if(typeof saved.resolveHorizontal==="boolean")$("limit-resolve-horizontal").checked=saved.resolveHorizontal;
    if(typeof saved.quantitiesOpen==="boolean")$("limits-quantities").open=saved.quantitiesOpen;
    if(typeof saved.customOpen==="boolean")$("limit-custom-settings").open=saved.customOpen;
    $("limit-chart").value=s.chart;
  }

  function publishSharedResults(resultOverride = null) {
    const shared=window.AtlasNumerics;
    if(!shared || !s.request)return;
    const graph=graphDescriptor(s.request);
    if(s.mesh && s.meshMeta)shared.putMesh(graph,{...s.meshMeta,result:s.mesh});
    if(!s.nodes.length)return;
    const result=resultOverride ?? (s.pending.has("trace") || !s.result ? {
      ...(s.result??{}),status:"partial",nodes:s.nodes,
      completed:s.nodes.filter(good).length,total:s.progressTotal||s.nodes.length,
      message:"Checked samples saved before leaving the page; this approach is incomplete.",
      stop:{kind:"page_change",resumeAvailable:false}
    } : {...s.result,nodes:s.nodes});
    shared.putTrace(graph,s.request,result);
    // ε = 1 is exactly the chosen segment start. Reuse its checked pair as
    // Explore's clicked point; a near-boundary overview node is not that point.
    const start=s.request.path?.start;
    const first=Array.isArray(start) ? s.nodes.find(node=>good(node) && /^1(?:\.0*)?(?:e[+-]?0+)?$/i.test(String(node.epsilon))) : null;
    if(first){
      const displayBranch=branch=>{
        const re=Number(branch.tau?.re ?? branch.tauExact?.[0]),im=Number(branch.tau?.im ?? branch.tauExact?.[1]);
        return Number.isFinite(re)&&Number.isFinite(im)&&im>0?{...branch,tau:{re,im}}:null;
      };
      const base=displayBranch(first.base),attached=displayBranch(first.result);
      if(base&&attached)shared.putPoint(graph,{source:start,result:{r:Number(start[0]),s:Number(start[1]),status:"ok",base,result:attached}});
    }
  }
  function restoreSharedMesh() {
    const saved=window.AtlasNumerics?.getMesh(graphDescriptor(),{inset:.04});
    if(!saved?.result?.samples)return false;
    s.mesh=saved.result;s.meshMeta={resolution:saved.resolution,inset:saved.inset};
    status("image",`${s.mesh.status==="partial"?"Partial interior":"Saved interior"}: ${s.mesh.samples.length} points.`);
    return true;
  }

  function requestForControls() {
    const path=controlPath(),target="1e-"+$("limit-depth").value;
    const sampling={path,target,startEpsilon:importedStartEpsilon,samplesPerDecade:Number($("limit-density").value)};
    return {...graphDescriptor(),seed:activeCase.seed,decoration:activeCase.operation,weights:activeCase.weights,...sampling};
  }

  function restoreExactSharedTrace() {
    const request=requestForControls();
    const saved=window.AtlasNumerics?.getTrace(graphDescriptor(),request);
    if(!saved?.result?.nodes?.length)return false;
    const own=window.AtlasResults?.get(FAMILY_RESULTS,graphResultKey(request))?.trace;
    const sameSavedPath=own?.request && window.AtlasNumerics.traceKey(own.request)===window.AtlasNumerics.traceKey(request);
    const selectedEpsilon=sameSavedPath && good(own.nodes?.[own.step]) ? own.nodes[own.step].epsilon : null;
    s.extension=null;
    s.request=request;s.draftDirty=false;s.detachedRestore=true;s.traceJob=null;
    s.view=null;s.zoomCorners=[];
    const result={...saved.result,stop:saved.result.stop?{...saved.result.stop,resumeAvailable:false}:undefined};
    delete result.resumeToken;delete result.resumeAvailable;
    finishTrace(result);
    if(selectedEpsilon!==null){
      const index=s.nodes.findIndex(node=>node.epsilon===selectedEpsilon && good(node));
      if(index>=0){s.step=index;s.view=null;render();if(Number($("limit-zoom").value)>0)centerSelected();saveFamilyResults();}
    }
    status("trace",`${saved.result.status==="ok"?"Saved trace":"Partial trace"}: ${s.nodes.filter(good).length} checked samples.`);
    $("limits-applied").textContent=`Restored this graph’s exact source path and sampling settings from Explore / Limits & Zoom. ε from ${request.startEpsilon} to ${request.target}; ${request.samplesPerDecade} samples per decade. Display settings are retained.`;
    return true;
  }
  const epsilonKey=value=>window.AtlasNumerics?.decimalKey(value)??String(Number(value));
  function extensionView(request) {
    const compatible=s.request && window.AtlasNumerics &&
      window.AtlasNumerics.graphKey(graphDescriptor(s.request))===window.AtlasNumerics.graphKey(graphDescriptor(request)) &&
      window.AtlasNumerics.traceKey({...s.request,target:null})===window.AtlasNumerics.traceKey({...request,target:null});
    return {selectedEpsilon:compatible?s.nodes[s.step]?.epsilon:null,view:compatible?s.view:null};
  }
  function mergeExtensionNodes(incoming) {
    const merged=new Map((incoming??[]).filter(Boolean).map(node=>[epsilonKey(node.epsilon),node]));
    for(const node of s.extension?.prefix??[])merged.set(epsilonKey(node.epsilon),node);
    return [...merged.values()].sort((a,b)=>Number(b.epsilon)-Number(a.epsilon));
  }
  function keepExtensionSelection() {
    const selected=s.extension?.selectedEpsilon;
    const index=selected==null?-1:s.nodes.findIndex(node=>good(node)&&epsilonKey(node.epsilon)===epsilonKey(selected));
    s.step=index>=0?index:lastGood();
  }
  function prepareExtension(candidate,view) {
    s.extension=candidate?{prefix:candidate.result.nodes,selectedEpsilon:view.selectedEpsilon}:null;
    if(!candidate)return;
    s.nodes=[...candidate.result.nodes];s.view=view.view;s.draftDirty=false;
    keepExtensionSelection();
    status("trace",`Keeping ${s.nodes.length} checked samples. Continuing from ε = ${s.nodes.at(-1).epsilon} toward ${s.request.target}; only new points are computed.`);
    $("limits-applied").textContent+=` ${s.nodes.length} existing checked samples are retained; the last point is rechecked at the required precision before continuing.`;
  }
  function acceptTraceProgress(progress) {
    if(progress.solverVersion)s.result={...(s.result??{}),solverVersion:progress.solverVersion,branchValidation:progress.branchValidation,geometryChecks:progress.geometryChecks,
      ...(progress.bouquetContinuation?{bouquetContinuation:progress.bouquetContinuation,tailBranchValidationStatus:progress.tailBranchValidationStatus??"continued"}:{})};
    if(s.extension){
      const incoming=Array.isArray(progress.nodes)?progress.nodes:[...s.nodes,...(progress.node?[progress.node]:[])];
      s.nodes=mergeExtensionNodes(incoming);keepExtensionSelection();
      status("trace",`Retained ${s.extension.prefix.length}; added ${Math.max(0,s.nodes.filter(good).length-s.extension.prefix.length)} checked samples. Last computed ε = ${progress.node?.epsilon??s.nodes.findLast(good)?.epsilon??"—"}.`);
    }else{
      if(Array.isArray(progress.nodes))s.nodes=progress.nodes;else if(progress.node)s.nodes.push(progress.node);
      s.step=lastGood();
      status("trace",`Trace samples: ${progress.completed}/${progress.total}; ε = ${progress.node?.epsilon??"—"}.`);
    }
  }
  function readResultCache(kind, key) {
    return window.AtlasResults?.get(cacheNamespace(kind), key) ?? null;
  }

  function storeResultCache(kind, key, value) {
    window.AtlasResults?.set(cacheNamespace(kind), key, value);
  }

  // A family snapshot contains numerical data only. It never owns a worker or
  // restores a continuation token. The common bounded cache serializes a copy,
  // so later solver progress cannot mutate a saved family view.
  const FAMILY_RESULTS = "limits-family-v1";
  function graphResultKey(request = null) {
    const identity = request
      ? request.caseId
        ? {family: request.seed, caseId: request.caseId, caseConfig: request.caseConfig}
        : {family: "lame", counts: request.counts}
      : ({family: activeCase.seed, caseId: activeCase.id, caseConfig: activeCase.caseConfig});
    return window.AtlasResults?.key(identity) ?? JSON.stringify(identity);
  }

  function saveFamilyResults() {
    const cache = window.AtlasResults;
    if (!cache || !s.request || (!s.mesh && !s.nodes.length)) return;
    const key = graphResultKey(s.request), previous = cache.get(FAMILY_RESULTS, key);
    let trace = previous?.trace ?? null;
    if (s.nodes.length) {
      const running = s.pending.has("trace");
      const result = running || !s.result ? {
        ...(s.result ?? {}), status: "partial", nodes: s.nodes,
        completed: s.nodes.filter(good).length,
        total: s.progressTotal || s.result?.total || s.nodes.length,
        message: running
          ? "The calculation was stopped when leaving this family. Checked samples are retained."
          : "Checked samples are retained; the full approach was not completed.",
        stop: {kind: "family_change", resumeAvailable: false}
      } : s.result;
      const {nodes: _storedSeparately, ...metadata} = result;
      trace = {
        request: s.request, nodes: s.nodes, result: metadata, step: s.step,
        precisionNote: $("limit-precision-note").textContent
      };
    }
    cache.set(FAMILY_RESULTS, key, {mesh: s.mesh ?? previous?.mesh ?? null, meshMeta:s.meshMeta ?? previous?.meshMeta ?? null, trace});
    publishSharedResults();
  }

  function restoreFamilyResults() {
    const saved = window.AtlasResults?.get(FAMILY_RESULTS, graphResultKey());
    if (!saved || (!saved.mesh && !saved.trace?.nodes?.length)) return false;
    s.mesh = saved.mesh ?? null;
    s.meshMeta=saved.meshMeta ?? null;
    s.view = null; s.zoomCorners = [];
    if (s.mesh) status("image", `${s.mesh.status==="partial"?"Partial interior":"Saved interior"}: ${s.mesh.samples.length} points.`);
    const trace = saved.trace, request = trace?.request;
    const path = request?.path;
    const pathName = path ? path.kind + path.index : "";
    const samplingMatches = request?.target === "1e-" + $("limit-depth").value &&
      Number(request.samplesPerDecade) === Number($("limit-density").value);
    const validStart=Array.isArray(path?.start) && path.start.length===2 && path.start.every(x=>Number.isFinite(Number(x)));
    const validPath = path && [...$("limit-path").options].some(option => option.value === pathName) && validStart;
    if (!trace?.nodes?.length || !samplingMatches || !validPath || graphResultKey(request) !== graphResultKey()) {
      if (!s.mesh) return false;
      status("trace", trace?.nodes?.length
        ? "The compact image is restored. Its saved trace uses different sampling settings; Compute approach to use the current settings."
        : "The compact image is restored. Compute approach to evaluate this source path.");
      $("limits-applied").textContent = "Computed compact images restored for this labelled graph. The source path and sampling settings above are retained.";
      $("limit-result-summary").textContent = "Compact images restored · trace ready to compute.";
      return true;
    }

    // The path belongs to the computation, not the display. Restore its actual
    // start and target together; a mapped path from another family must never
    // silently relabel coordinates that were computed on this one.
    $("limit-path").value = pathName;
    $("limit-position").value = String(path.position);
    importedStartEpsilon=request.startEpsilon??"1";
    const actualStart=validStart?path.start:legacySourcePoint(path,importedStartEpsilon);
    s.start = {r: Number(actualStart[0]), s: Number(actualStart[1])};
    s.request = request;
    s.nodes = trace.nodes;
    s.result = {...trace.result, nodes: s.nodes, resumeAvailable: false,
      stop: trace.result?.stop ? {...trace.result.stop, resumeAvailable: false} : undefined};
    delete s.result.resumeToken;
    s.traceJob = null; s.detachedRestore = true; s.draftDirty = false;
    s.step = Number.isInteger(trace.step) && good(s.nodes[trace.step]) ? trace.step : lastGood();
    s.progressTotal = Number(s.result.total) || s.nodes.length;
    if (trace.precisionNote) $("limit-precision-note").textContent = trace.precisionNote;
    updateProgress(s.result.completed ?? s.nodes.filter(good).length, s.progressTotal, s.nodes.findLast(good)?.epsilon);
    const complete = s.result.status === "ok";
    status("trace", `${complete ? "Computed trace" : "Partial trace"} restored from this browser: ${s.nodes.filter(good).length} checked samples. No new computation was started.` +
      (complete ? "" : " A worker checkpoint is not retained after switching families; Compute approach to retry."));
    $("limits-applied").textContent = `Saved ${complete ? "computed" : "partial"} family view restored. Path ${pathName}, t = ${path.position}; start (${actualStart.join(", ")}); ε from ${request.startEpsilon ?? "1"} to ${request.target}; ${request.samplesPerDecade} samples per decade. Display settings are retained.`;
    return true;
  }

  function syncComputationURL() {
    const url = new URL(location.href);
    for (const [key, value] of Object.entries({
      path: $("limit-path").value, t: $("limit-position").value,
      fromR: s.start.r, fromS: s.start.s, depth: $("limit-depth").value,
      density: $("limit-density").value, chart: s.chart, autocompute: 0
    })) url.searchParams.set(key, String(value));
    url.searchParams.set("startEpsilon",String(importedStartEpsilon));
    url.searchParams.delete("sourceMode");
    url.searchParams.delete("fromOverview");
    query.delete("fromOverview");
    history.replaceState(null, "", url);
    window.AtlasUI?.refreshNavigation(url.searchParams);
  }
  const query = new URLSearchParams(location.search);
  let importedStartEpsilon=Number(query.get("startEpsilon"))>0 && Number(query.get("startEpsilon"))<=1
    ? query.get("startEpsilon") : "1";

  ["a", "b", "c"].forEach((name, i) => {
    $("limit-" + name).value = query.get(name) ?? [1, 0, 0][i];
  });

  if (/^(edge|vertex)[012]$/.test(query.get("path"))) {
    $("limit-path").value = query.get("path");
  }

  function ensureCoordinateOptions() {
    if (![...$("limit-chart").options].some(option => option.value === "half")) {
      const option=document.createElement("option");
      option.value="half";option.textContent="w = τ/(1 − 2τ) (marked cusp 1/2 → ∞)";
      $("limit-chart").append(option);
    }
    for (const [value, label] of [["8", "10⁻⁸ (Explore overview)"], ["20", "10⁻²⁰"], ["40", "10⁻⁴⁰"], ["80", "10⁻⁸⁰"]]) {
      if (![...$("limit-depth").options].some(option => option.value === value)) {
        const option=document.createElement("option");option.value=value;option.textContent=label;
        $("limit-depth").append(option);
      }
    }
    for (const option of $("limit-depth").options) {
      option.disabled=false;
      if (Number(option.value)>10) option.textContent=`10⁻${option.value.replace(/\d/g,n=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(n)])}${""}`;
    }
  }
  ensureCoordinateOptions();

  if (["3", "5", "7", "8", "10", "20", "40", "80"].includes(query.get("depth"))) {
    $("limit-depth").value = query.get("depth");
  }

  if (query.has("t")) {
    $("limit-position").value = query.get("t");
  }

  if (["1", "2", "3"].includes(query.get("density"))) $("limit-density").value = query.get("density");

  if (query.has("fromR") && query.has("fromS")) {
    s.start = importedStart(Number(query.get("fromR")), Number(query.get("fromS")));
  }

  // Keep bookmarked choices. A fresh costly mixed graph starts with a short
  // simulation preview; extra working precision stays inside the solver.
  let mixedCase=(activeCase.operation === "repeated") && activeCase.seed !== "bouquet" && !activeCase.caseConfig?.connector;
  const costlyBouquet=()=>(activeCase.seed==="bouquet") && (Math.max(...activeCase.weights)>12 || activeCase.weights.reduce((a,b)=>a+b,0)>24);
  if((mixedCase||costlyBouquet()) && !query.has("depth")) $("limit-depth").value="5";
  // Display sampling is separate from adaptive branch continuation and precision.
  let presets = {preview: [5, 1], compare: [10, 1], deeper: [20, 2]};
  function updateComputationSummary() {
    const depth=Number($("limit-depth").value), density=Number($("limit-density").value);
    const sum=0;
    const extra=0;
    const match=Object.entries(presets).find(([,v])=>v[0]===depth&&v[1]===density);
    $("limit-preset").value=match?.[0] ?? "custom";
    const exponent=String(depth).replace(/\d/g,n=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(n)]);
    $("limit-computation-summary").textContent=`ε: ${importedStartEpsilon} → 10⁻${exponent} · ${Math.max(1,Math.ceil((depth+Math.log10(Number(importedStartEpsilon)))*density)+1)+extra} samples per branch${extra?` (${extra} extra tail checks)`:""}.`+
      (depth>20 ? " Research depth: partial traces are possible." : "");
    const desired=requestForControls(),exact=window.AtlasNumerics?.getTrace(graphDescriptor(),desired);
    const continuation=exact?.result?.status==="ok"?null:window.AtlasNumerics?.getExtension(graphDescriptor(),desired);
    $("limits-run").textContent=continuation?(Number(continuation.request.target)===Number(desired.target)?"Continue saved trace":"Continue to deeper ε"):"Compute approach";
    if(continuation)$("limit-computation-summary").textContent+=` ${continuation.result.nodes.length} checked samples retained.`;
    if(s.result?.extended&&s.request&&window.AtlasNumerics?.graphKey(graphDescriptor(s.request))===window.AtlasNumerics.graphKey(graphDescriptor())&&window.AtlasNumerics.traceKey(s.request)===window.AtlasNumerics.traceKey(desired))
      $("limit-computation-summary").textContent=`ε: ${importedStartEpsilon} → 10⁻${exponent} · ${s.nodes.filter(good).length} checked samples (${s.result.reusedCount} retained; ${s.result.newCount} new).`;
  }
  function refreshDepthAdvice() {
    let advice=initialText["limit-depth-advice"]+" Optional deeper checks up to 80 decades are available for TV families.";
    if((activeCase.seed==="bouquet")){
      const path=controlPath(),f=path.kind==="edge"?window.TVCuspFormulas?.lookup(activeCase.caseConfig,path.index):null;
      advice+=` For this three-loop case${f?`, M = ${f.M}`:""}, start with ${costlyBouquet()?"5":"10"} decades and one displayed sample per decade. The solver selects working precision adaptively; display digits are not a precision request.`;
      if(f?.M>8)advice+=" This large resonance makes the leading asymptotic slow to settle. A tiny ε alone does not imply that the proposed K is already accurate; compare successive computed values before requesting a longer run.";
    }
    $("limit-depth-advice").textContent=advice;
  }
  function refreshSamplingLabels() {
    mixedCase=(activeCase.operation === "repeated") && activeCase.seed !== "bouquet" && !activeCase.caseConfig?.connector;
    presets={preview:[5,1],compare:[10,1],deeper:[20,2]};
    $("limit-preset").options[0].textContent=mixedCase || costlyBouquet() ? "Quick preview · 5 decades (recommended for larger mixed graphs)" : "Quick preview · 5 decades, 1 sample/decade";
    $("limit-preset").options[1].textContent="Compare limits · 10 decades, 1 sample/decade";
    $("limit-preset").options[2].textContent="Deeper check · 20 decades, 2 samples/decade";
    $("limit-research-advice").hidden=!true;
    refreshDepthAdvice();
  }
  refreshSamplingLabels();
  $("limit-preset").onchange=()=>{
    const chosen=presets[$("limit-preset").value];
    if(chosen){$("limit-depth").value=String(chosen[0]);$("limit-density").value=String(chosen[1]);changeDraft();}
    else $("limit-custom-settings").open=true;
  };
  updateComputationSummary();
  if($("limit-preset").value === "custom" || Number($("limit-depth").value)>20) $("limit-custom-settings").open=true;

  const good = node =>
    node?.status === "ok" && node.result?.status === "ok";

  const values = (node, branch = "result") =>
    node?.[branch]?.views?.[s.chart];

  const lastGood = () =>
    s.nodes.findLastIndex(good);

  const fmt = (x, digits = 6) =>
    Number.isFinite(Number(x))
      ? Number(Number(x).toPrecision(digits)).toString()
      : "—";

  const escape = x =>
    String(x)
      .replaceAll("&", "&amp;")
      .replaceAll("<", "&lt;")
      .replaceAll('"', "&quot;");

  /*
    Subtract decimal strings BEFORE converting the small difference
    to an ordinary JavaScript number.
  */
  function difference(a, b) {
    function parts(text) {
      const [mantissa, exponent = "0"] =
        String(text).toLowerCase().split("e");

      const [whole, fraction = ""] = mantissa.split(".");

      return {
        n: BigInt(whole + fraction),
        e: Number(exponent) - fraction.length
      };
    }

    const x = parts(a);
    const y = parts(b);
    const e = Math.min(x.e, y.e);

    const integer = x.n * 10n ** BigInt(x.e - e) -
      y.n * 10n ** BigInt(y.e - e);
    if (integer === 0n) return 0;
    // Convert only the leading significant digits. Directly converting a
    // 900-digit integer overflows even when the final small difference is finite.
    const negative = integer < 0n, digits = String(negative ? -integer : integer);
    const head = digits.slice(0, 17);
    return (negative ? -1 : 1) * Number(head) * 10 ** (e + digits.length - head.length);
  }

  /* Axis coordinates are origin + small offset. Add as integers before
     rounding the label to the visible scale; Number(origin) would erase
     every tick difference when the lens is below double precision. */
  function absoluteTick(origin, displacement, span) {
    function parts(text) {
      const [mantissa, exponent = "0"] = String(text).toLowerCase().split("e");
      const [whole, fraction = ""] = mantissa.split(".");
      return {n: BigInt(whole + fraction), e: Number(exponent) - fraction.length};
    }
    const a = parts(origin), b = parts(displacement);
    const exponent = Math.min(a.e, b.e);
    let integer = a.n * 10n ** BigInt(a.e - exponent) + b.n * 10n ** BigInt(b.e - exponent);
    const digits = Math.max(0, Math.min(290, Math.ceil(-Math.log10(Math.abs(span) || 1)) + 2));
    const shift = exponent + digits;
    if (shift >= 0) integer *= 10n ** BigInt(shift);
    else {
      const divisor = 10n ** BigInt(-shift), negative = integer < 0n;
      const magnitude = negative ? -integer : integer;
      integer = (magnitude + divisor / 2n) / divisor * (negative ? -1n : 1n);
    }
    const sign = integer < 0n ? "−" : "";
    let text = String(integer < 0n ? -integer : integer);
    if (digits) {
      text = text.padStart(digits + 1, "0");
      text = text.slice(0, -digits) + "." + text.slice(-digits);
      text = text.replace(/0+$/, "").replace(/\.$/, "");
    }
    return sign + text;
  }

  const offset = (z, origin) => [
    difference(z[0], origin[0]),
    difference(z[1], origin[1])
  ];

  function meshValue(sample, field = "result") {
    const branch = sample?.[field];
    if (branch?.views?.[s.chart]) return branch.views[s.chart];

    if (branch?.status !== "ok") return null;

    const {re: x, im: y} = branch.tau;

    if (s.chart === "raw") {
      return [String(x), String(y)];
    }

    if (s.chart === "half") {
      const den=(1-2*x)**2+4*y*y;
      return den>0 ? [String((x-2*x*x-2*y*y)/den),String(y/den)] : null;
    }
    const u = s.chart === "zero" ? x : 1 - x;
    const den = u * u + y * y;

    return [
      String(s.chart === "zero" ? 1 - x / den : u / den),
      String(y / den)
    ];
  }

  function importedStart(r,t) {return compactPoint(r,t);}

  function compactPoint(r,t) {return cases.project(r,t,.04);}

  function pathNote() {
    const vertex = $("limit-path").value.startsWith("vertex");
    $("limit-position").disabled = vertex;
    $("limit-path-note").textContent = vertex
      ? "The diamond is fixed at the selected vertex. Drag the hollow start circle to choose its direction; the edge-position control does not apply to a vertex."
      : "Drag the diamond along a boundary edge, or choose the boundary and its position t above. Drag the hollow circle anywhere inside the compact starting region.";
    if ($("limit-start-note")) $("limit-start-note").textContent =
      (("ε = 1 is the hollow starting circle; ε → 0 approaches the boundary diamond along the displayed segment."+(Number(importedStartEpsilon)!==1?` The imported overview trace samples ε from ${importedStartEpsilon}; the original segment parameter is retained.`:"")))+
      " Drag either handle in either diagram. No new computation starts until Compute approach. Dragging chooses a new compact start and resets the sampling start to ε = 1.";
  }

  function controlPath() {
    const name = $("limit-path").value;
    return {
      kind: name.startsWith("edge") ? "edge" : "vertex",
      index: Number(name.at(-1)),
      position: $("limit-position").value,
      ...(({start: [String(s.start.r), String(s.start.s)]}))
    };
  }

  // Kept only to migrate links made before freely chosen starts were added.
  function legacySourcePoint(path, epsilon) {
    const e = Number(epsilon), t = Number(path.position);
    const weights = [], other = [0,1,2].filter(i => i !== path.index);
    weights[path.index] = path.kind === "edge" ? e : 1-e;
    const remaining = path.kind === "edge" ? 1-e : e;
    weights[other[0]] = remaining*t;
    weights[other[1]] = remaining*(1-t);
    return [(1-weights[2])/2, (1-weights[0])/2];
  }

  function boundaryPoint(path) {
    
    const h=activeCase.top, x=activeCase.maxR ?? 1, t=Number(path.position);
    return path.kind === "vertex" ? activeCase.sourceVertices[path.index]
      : [[x,h*t],[x*t,h],[x*t,h*(1-t)]][path.index];
  }

  function sourcePoint(path, epsilon) {
    if (!path.start) return legacySourcePoint(path, epsilon);
    const e = Number(epsilon), end = boundaryPoint(path);
    return end.map((x,i) => (1-e)*x + e*Number(path.start[i]));
  }

  function nodeSource(node, path) {
    return Number.isFinite(Number(node?.r)) && Number.isFinite(Number(node?.s))
      ? [Number(node.r), Number(node.s)]
      : node ? sourcePoint(path, node.epsilon) : null;
  }

  function changeDraft(event) {
    if(event?.target && ["limit-path","limit-position"].includes(event.target.id))resetImportedPath();
    trackDepthChoice();
    updateComputationSummary();
    cancel();
    s.draftDirty=true;
    pathNote();
    renderSource();
    renderEvidence();
    renderInterruption();
    if(s.request) $("limits-applied").textContent="Settings changed. Result plots below retain the previous computation until Compute approach is pressed.";
  }
  function resetImportedPath() {
    importedStartEpsilon="1";
    query.delete("sourceMode");query.delete("fromOverview");
  }

  function updateSourceHandle(kind,p) {
    resetImportedPath();
    if(kind==="start") s.start=compactPoint(p.r,p.s);
    else snapBoundary(p.r,p.s);
    changeDraft();
  }

  function connectSourceEditor(svg,attachment) {
    let active=null;
    const eventPoint=event=>{
      const box=svg.getBoundingClientRect(), scale=Math.min(box.width/400,box.height/400);
      const x=(event.clientX-box.left-(box.width-400*scale)/2)/scale;
      const y=(event.clientY-box.top-(box.height-400*scale)/2)/scale;
      return cases.svgInverse(x,y,attachment);
    };
    svg.addEventListener("pointerdown",event=>{
      const handle=event.target.closest?.("[data-handle]");
      if(!handle || event.button>0) return;
      event.preventDefault();
      const point=eventPoint(event), target=boundaryPoint(controlPath());
      const startDistance=(point.r-s.start.r)**2+(point.s-s.start.s)**2;
      const targetDistance=(point.r-target[0])**2+(point.s-target[1])**2;
      const kind=Math.abs(startDistance-targetDistance)<1e-14
        ? handle.dataset.handle : startDistance<targetDistance ? "start" : "target";
      active={kind,id:event.pointerId};
      svg.setPointerCapture(event.pointerId);
      updateSourceHandle(active.kind,point);
    });
    svg.addEventListener("pointermove",event=>{
      if(active?.id===event.pointerId) updateSourceHandle(active.kind,eventPoint(event));
    });
    const finish=event=>{
      if(active?.id!==event.pointerId) return;
      const kind=active.kind;
      active=null;
      svg.querySelector(`[data-handle="${kind}"]`)?.focus();
    };
    svg.addEventListener("pointerup",finish);
    svg.addEventListener("pointercancel",finish);
    svg.addEventListener("lostpointercapture",()=>{active=null;});
    svg.addEventListener("keydown",event=>{
      const handle=event.target.closest?.("[data-handle]");
      const delta={ArrowLeft:[-1,0],ArrowRight:[1,0],ArrowUp:[0,1],ArrowDown:[0,-1]}[event.key];
      if(!handle||!delta) return;
      event.preventDefault();
      const kind=handle.dataset.handle, p=kind==="start"?[s.start.r,s.start.s]:boundaryPoint(controlPath());
      const mapped=cases.map(...p,attachment);
      const step = event.shiftKey ? .01 : .002;
      updateSourceHandle(kind,inverse({r:mapped.r+delta[0]*step,s:mapped.s+delta[1]*step}));
      svg.querySelector(`[data-handle="${kind}"]`)?.focus();
    });
  }

  function infinityReference(base=false) {
const path = s.request?.path ?? controlPath();
if (path.kind !== "edge") return null;

      const end=boundaryPoint(path), info=cases.asymptotic(path.index,end[0],end[1],base);
      if(window.TVCuspPredictions){
        return window.TVCuspPredictions.reference(info,s.nodes.map(node=>window.TVCuspPredictions.chartValues(node,base?"base":"result",s.chart)),{weights:base?[1,0,0,0]:activeCase.weights,chart:s.chart});
      }
      if (!info || info.known === false || info.chart !== s.chart || !Number.isFinite(info.x0)) return null;
      let denominator=info.denominator ?? 1, numerator=info.numerator;
      if (!Number.isSafeInteger(numerator) || !Number.isSafeInteger(denominator) || denominator <= 0) {
        denominator=1;
        while (denominator<4096 && Math.abs(info.x0*denominator-Math.round(info.x0*denominator))>1e-10) denominator++;
        if (Math.abs(info.x0*denominator-Math.round(info.x0*denominator))>1e-10) return null;
        numerator=Math.round(info.x0*denominator);
      }
      return {numerator,denominator,value:info.x0,label:info.xLabel??(denominator===1?String(numerator):`${numerator}/${denominator}`),status:"proposed",kind:"Proposed graph limit",info};
    
  }

  function fractionDecimal(numerator, denominator, digits = 950) {
    const n = BigInt(numerator), d = BigInt(denominator), negative=n<0n;
    const scaled = (negative?-n:n) * 10n ** BigInt(digits) / d;
    const text = String(scaled).padStart(digits + 1, "0");
    return (negative?"-":"") + text.slice(0, -digits) + "." + text.slice(-digits);
  }

  function rawTau(node, branch = "result") {
    if (node?.[branch]?.status !== "ok" || (branch === "result" && !good(node))) return null;
    const z = node[branch].tauExact ?? node[branch].views?.raw;
    const tau = z ? {re: Number(z[0]), im: Number(z[1])} : node[branch].tau;
    return tau && Number.isFinite(tau.re) && Number.isFinite(tau.im) && tau.im > 0 ? tau : null;
  }

  function inverseChart(x, y) {
    if (s.chart === "raw") return {re:x, im:y};
    if (s.chart === "half") {
      const den=(1+2*x)**2+4*y*y;
      return den>0 ? {re:(x+2*x*x+2*y*y)/den,im:y/den} : null;
    }
    const u = s.chart === "zero" ? 1 - x : x;
    const den = u*u + y*y;
    if (!(den > 0)) return null;
    return s.chart === "zero"
      ? {re:u/den, im:y/den}
      : {re:1-x/den, im:y/den};
  }

  function renderOverview() {
    const api = window.AtlasView;
    if (!api?.drawDomain) return;
    const result = s.nodes.map(node => rawTau(node));
    const base = s.nodes.map(node => rawTau(node, "base"));
    const paths = [];
    if ((s.mesh?.boundaries)) {
      for (const field of [...($("limit-compare").checked?["base"]:[]),"result"]) {
        for (const arc of s.mesh.boundaries[field]??[]) paths.push({points:arc.points.map(z=>({re:Number(z[0]),im:Number(z[1])})),color:field==="base"?baseColor+"88":color+"aa",width:1.1,radius:0});
      }
    }
    if ($("limit-compare").checked) paths.push({points:base, color:baseColor, selected:s.step});
    paths.push({points:result, color, selected:s.step});
    const ref = s.chart === "raw" ? infinityReference() : null;
    api.drawDomain($("limit-overview"), {
      height:560,
      markedStrip:true,
      meshLayers:s.mesh ? [
        ...($("limit-compare").checked ? [{mesh:s.mesh, field:"base", color:baseColor}] : []),
        {mesh:s.mesh, field:"result", color}
      ] : [],
      paths,
      references:ref ? [
        ...($("limit-compare").checked && infinityReference(true) ? [{re:infinityReference(true).value, label:`${baseLabel}: Re τ → ${infinityReference(true).label}`, color:baseColor}] : []),
        {re:ref.value, label:`Re τ → ${ref.label}`, color}
      ] : [],
      points:result[s.step] ? [{tau:result[s.step], color}] : [],
      zoomWindow:s.zoomCorners.length ? {corners:s.zoomCorners} : null
    });
    $("limit-overview-note").textContent = (((s.mesh?.boundaries?"Thin curves are computed finite boundary arcs, with their ends truncated. ":"")+"Shading is a compact interior image. "))+"Solid traces join consecutive resolved samples; unresolved gaps remain open and no segment is extrapolated to an ideal cusp. Dashed image contours bound the compact interior sample. The dark dashed box locates the local view; very small footprints are enlarged for visibility." + (s.chart === "raw" ? "" : " In transformed coordinates, the box approximately projects the lens corners; the actual curved viewport boundary is not drawn.");
  }

  function pausePlayback() {
    s.playing = false;
    clearTimeout(s.playTimer);
    s.playTimer = null;
    $("limit-play").textContent = "▶ Play approach";
    $("limit-play").setAttribute("aria-pressed", "false");
  }

  function startPlayback() {
    if (s.playing) { pausePlayback(); return; }
    const available = s.nodes.map((node, i) => good(node) ? i : -1).filter(i => i >= 0);
    if (available.length < 2) return;
    s.playing = true;
    $("limit-play").textContent = "Ⅱ Pause";
    $("limit-play").setAttribute("aria-pressed", "true");
    if (s.step >= available.at(-1)) selectStep(available[0]);
    const tick = () => {
      if (!s.playing) return;
      const next = available.find(i => i > s.step);
      if (next == null) { pausePlayback(); return; }
      selectStep(next);
      if (next === available.at(-1)) pausePlayback();
      else s.playTimer = setTimeout(tick, 650);
    };
    s.playTimer = setTimeout(tick, 650);
  }

  function status(kind, text, error = false) {
    const el = $("limits-" + kind + "-status");

    el.textContent = text;
    el.classList.toggle("error", error);
  }

  function settle(kind,keepWorker=false) {
    if(!keepWorker){s.workers[kind]?.terminate();delete s.workers[kind];}
    s.pending.delete(kind);

    $("limits-run").disabled = s.pending.size > 0;
    $("limits-cancel").disabled = !s.pending.size;
    $("limit-play").disabled = s.pending.has("trace") || s.nodes.filter(good).length < 2;
    renderInterruption();
  }

  function cancel(showNotice=false) {
    const wasRunning=s.pending.has("trace");
    pausePlayback();
    ++s.job;

    for (const kind of [...s.pending]) {
      settle(kind);

      status(
        kind === "mesh" ? "image" : "trace",
        "Cancelled. Computed samples are retained."
      );
    }
    for(const worker of Object.values(s.workers))worker.terminate();
    s.workers={};s.traceJob=null;
    if(showNotice && wasRunning){s.result={status:"partial",nodes:s.nodes,message:"Stopped by the user",stop:{kind:"cancelled",resumeAvailable:false}};renderInterruption();}
  }

  function startWorker(kind, file, data, handle) {
    s.pending.add(kind);
    if(kind === "trace"){
      s.traceJob={file,data,handle};s.progressTotal=0;
      $("limit-job-progress").hidden=false;$("limit-progress").removeAttribute("value");
      $("limit-progress-label").textContent="Preparing the marked branch…";
      $("limit-interruption").hidden=true;
    }

    $("limits-run").disabled = true;
    $("limits-cancel").disabled = false;

    try {
      // Module detection must use the path, even when the caller supplied
      // a version query. Keep a single current version on the worker URL.
      const workerPath = file.split("?")[0];
      const worker = new Worker(workerPath + "?v=20261025-tv-maintenance1", workerPath.endsWith(".mjs") ? {type:"module"} : undefined);

      s.workers[kind] = worker;

      worker.onmessage = ({data: reply}) => {
        if (reply.id === s.job) {
          if(kind === "trace" && reply.type === "progress"){
            const progress=reply.result??reply;
            updateProgress(progress.completed,progress.total,progress.node?.epsilon);
          }
          handle(reply);
        }
      };

      worker.onerror = () => {
        if (s.workers[kind] !== worker) return;

        settle(kind);

        status(
          kind === "mesh" ? "image" : "trace",
          "Could not run " + file +
          ". Check the file is beside limits.html and open " +
          "the site through GitHub Pages or a local server.",
          true
        );
        if(kind === "trace")traceError("The numerical worker could not run. Check the local server and reload the page.");
      };

      worker.postMessage(data);
    } catch (error) {
      settle(kind);

      status(
        kind === "mesh" ? "image" : "trace",
        error.message,
        true
      );
      if(kind === "trace")traceError(error.message);
    }
  }

  function finishTrace(result) {
    if(s.extension){
      const nodes=mergeExtensionNodes(result.nodes),reused=s.extension.prefix.length;
      result={...result,nodes,extended:true,reusedCount:reused,newCount:Math.max(0,nodes.filter(good).length-reused)};
    }
    s.result = result;
    updateProgress(result.completed??result.nodes.filter(good).length,result.total??result.nodes.length,result.nodes.findLast(good)?.epsilon);
    {
      const precisions = result.nodes.flatMap(node => [node.base?.precisionDigits, node.result?.precisionDigits])
        .map(Number).filter(value => Number.isFinite(value) && value > 0);
      const range = precisions.length ? `${Math.min(...precisions)}–${Math.max(...precisions)} working digits` : "adaptive working precision";
      $("limit-precision-note").textContent=(result.precisionDescription ??
        `Stored coordinates retain the solver's decimal strings (${range}). Compact images are locators; cusp points use stable evaluation and adaptive precision where available. Error figures are numerical diagnostics. Unresolved continuation stops without adding a boundary endpoint.`)+(result.tailBranchValidation?` ${result.tailBranchValidation}.`:"")+(result.message?` ${result.message}`:"");
    }
    s.nodes = result.nodes;
    if(s.extension)keepExtensionSelection();else s.step = lastGood();

    if (
      s.fit !== "image" &&
      Number($("limit-zoom").value) === 0
    ) {
      s.view = null;
    }

    const completed = s.nodes.filter(good).length;

    status(
      "trace",
      (
        result.status === "ok"
          ? "Trace complete: "
          : "Trace stopped early: "
      ) +
      completed + (result.cached ? " stored samples. " : " checked samples. ") +
      (result.extended ? `Kept ${result.reusedCount} previous samples; ${result.newCount} new checked samples. ` : "") +
      (result.status!=="ok"?(result.message ?? ""):"") +
      (result.status!=="ok"&&result.tailBranchValidation ? ` ${result.tailBranchValidation}.` : "") +
      (
        result.elapsedMs != null
          ? ` ${(result.elapsedMs / 1000).toFixed(1)} seconds.`
          : ""
      )
    );

    render();

    if (Number($("limit-zoom").value) > 0) {
      centerSelected();
    }
    saveFamilyResults();
    updateComputationSummary();
  }

  function traceError(message){
    s.result={status:s.nodes.some(good)?"partial":"unresolved",nodes:s.nodes,
      completed:s.nodes.filter(good).length,total:s.progressTotal||s.nodes.length,
      message:String(message??"The numerical worker stopped."),stop:{kind:"error",resumeAvailable:false}};
    renderInterruption();
    saveFamilyResults();
  }

  function refreshCaseLabels() {
    
    baseLabel=activeCase.baseTitle;resultLabel=activeCase.title;
    $("limit-base-source-title").textContent=`Base · ${baseLabel}`;
    $("limit-attached-source-title").textContent=resultLabel;
    $("limit-compare-label").textContent=`Include ${baseLabel}: image and trace`;
    $("limit-base-legend").textContent=`● ${baseLabel}`;
    $("limit-result-legend").textContent=`● ${activeCase.operation === "seed" ? "Selected seed" : "Decorated trace"}`;
    $("limit-overview-caption").textContent="Marked τ representatives in the upper half-plane. The Γ(2) guides are a reference; the horizontal range expands to retain the selected branch. Height is compactified as Im τ / (1 + Im τ). Shading is the computed compact interior.";
    $("limit-precision-note").textContent="Cusp coordinates retain decimal strings from stable evaluation and adaptive precision. Compact image meshes are locators. The trace reports the precision and checks actually used and stops if its branch cannot be resolved. Diagnostics are not certified bounds.";
    $("limit-source-controls-note").textContent=activeCase.seed === "bouquet"
      ? "Choose an open source edge or drag its target diamond. Position t runs from 0 to 1 along that edge, with physical endpoint coordinate t/2; the cusp factor is sin(πt). Drag the starting circle to choose the approach."
      : initialText["limit-source-controls-note"];
    const selectedPath=$("limit-path").value;
    $("limit-path").innerHTML=activeCase.edgeLabels.map((label,i)=>`<option value="edge${i}">${escape(label)}</option>`).join("") + (activeCase.seed === "bouquet" ? "" : activeCase.sourceVertices.map((p,i)=>`<option value="vertex${i}">Vertex ${["A","B","C"][i]}: (${p.join(", ")}) · direction-dependent finite boundary</option>`).join(""));
    $("limit-path").value=(activeCase.seed === "bouquet" ? /^edge[012]$/ : /^(edge|vertex)[012]$/).test(selectedPath)?selectedPath:activeCase.seed === "bouquet" && /^vertex[012]$/.test(selectedPath)?selectedPath.replace("vertex","edge"):"edge0";
  }

  function refreshFamilyControls() {
    const selectedPath=$("limit-path").value;
    $("limits-legacy-attachment").hidden=true;
    ["limit-a","limit-b","limit-c"].forEach(id=>{$(id).disabled=true;});
    ensureCoordinateOptions();refreshSamplingLabels();
    refreshCaseLabels();
    $("limit-chart").value=s.chart;
  }

  function initializeCase() {
    
    $("limits-legacy-attachment").hidden=true;
    ["limit-a","limit-b","limit-c"].forEach(id=>{$(id).disabled=true;});
    refreshCaseLabels();
    if ((activeCase.seed === "bouquet" ? /^edge[012]$/ : /^(edge|vertex)[012]$/).test(query.get("path"))) $("limit-path").value=query.get("path");
    else $("limit-path").value="edge0";
    if (query.has("chart") && ["raw","zero","one","half"].includes(query.get("chart"))) s.chart=query.get("chart");
    else s.chart=controlPath().kind === "edge" ? cases.asymptotic(controlPath().index,...boundaryPoint(controlPath()))?.chart ?? "raw" : "raw";
    $("limit-chart").value=s.chart;
  }

  function snapBoundary(r,t) {
    const h=activeCase.top, x=activeCase.maxR ?? 1, pairs=[[[x,0],[x,h]],[[0,h],[x,h]],[[0,h],[x,0]]];
    let best=null;
    pairs.forEach(([a,b],index)=>{
      const dx=(b[0]-a[0])/x,dy=(b[1]-a[1])/h;
      const position=Math.max(.02,Math.min(.98,((r-a[0])/x*dx+(t-a[1])/h*dy)/(dx*dx+dy*dy)));
      const point=[a[0]+position*(b[0]-a[0]),a[1]+position*(b[1]-a[1])];
      const distance=((point[0]-r)/x)**2+((point[1]-t)/h)**2;
      if (!best||distance<best.distance) best={index,position,distance};
    });
    $("limit-path").value="edge"+best.index;
    $("limit-position").value=String(Number(best.position.toPrecision(12)));
  }

  function renderSource() {
    const path=controlPath(), start=[s.start.r,s.start.s], end=boundaryPoint(path);
    refreshDepthAdvice();
    const formulaPath=path,formulaConfig=activeCase.caseConfig;
    const formulaHTML=formulaPath.kind==="edge"?window.TVCuspFormulas?.html(formulaConfig,formulaPath.index):"";
    if($("limit-count-formula")){$("limit-count-formula").hidden=!formulaHTML;$("limit-count-formula").innerHTML=formulaHTML??"";}
    if($("limit-count-formula-help"))$("limit-count-formula-help").hidden=!formulaHTML;
    const node=!s.draftDirty?s.nodes[s.step]:null, current=nodeSource(node,path);
    const text=p=>`(${fmt(p[0],7)}, ${fmt(p[1],7)})`;
    function draw(svg,attached,tint) {
      cases.renderMonodromy(svg,{point:null,attached,title:attached?resultLabel:baseLabel});
      const map=p=>cases.svgPoint(...p,attached), a=map(start),b=map(end);
      const selected=current?map(current):null;
      const handle=(kind,q,shape,label)=>`<g data-handle="${kind}" class="source-handle source-handle-${kind}" tabindex="0" role="slider" aria-label="${label}" aria-valuemin="0" aria-valuemax="1" aria-valuenow="${kind === "start" ? start[0] : end[0]}" aria-valuetext="${text(kind === "start"?start:end)}" style="cursor:grab;touch-action:none"><circle cx="${q[0]}" cy="${q[1]}" r="16" fill="transparent"></circle>${shape}</g>`;
      const samples=!s.draftDirty?s.nodes.map(n=>{const p=nodeSource(n,path);if(!p)return "";const q=map(p);return `<circle cx="${q[0]}" cy="${q[1]}" r="2.4" fill="${good(n)?tint:"#929b92"}"></circle>`;}).join(""):"";
      svg.insertAdjacentHTML("beforeend",`<path d="M ${a.join(" ")} L ${b.join(" ")}" stroke="${tint}" stroke-width="1.8" stroke-dasharray="6 4" fill="none"></path>${samples}${selected?`<circle cx="${selected[0]}" cy="${selected[1]}" r="4.5" fill="${tint}" stroke="white" stroke-width="1.3"></circle>`:""}${handle("start",a,`<circle cx="${a[0]}" cy="${a[1]}" r="7" fill="white" stroke="${tint}" stroke-width="2.5"></circle>`,"Starting point. Arrow keys move R and S; Shift moves faster.")}${handle("target",b,`<path d="M ${b[0]} ${b[1]-8} l 8 8 l -8 8 l -8 -8 Z" fill="${tint}" stroke="white" stroke-width="1.8"></path>`,"Boundary target. Arrow keys move along the nearest edge; Shift moves faster.")}`);
    }
    draw($("limit-source"),false,baseColor);draw($("limit-source-attached"),true,color);
    const after=cases.map(...start,true), endAfter=cases.map(...end,true);
    $("limit-start-readout").textContent=`Start ○ · base (R,S) = ${text(start)} · selected character = ${text([after.r,after.s])}`;
    $("limit-target-readout").textContent=`Target ◆ · base (R,S) = ${text(end)} · selected character = ${text([endAfter.r,endAfter.s])}`;
    $("limit-source-readout").textContent=current?`Selected source (R,S) = ${text(current)} · ε = ${node.epsilon}`:"The diagrams share one path. Drag either starting circle or target diamond.";
    $("limit-source-caption").textContent=`${activeCase.title}. Dashed outline: full source; inner triangle: compact starting region. ○ starts the path, ◆ is the limiting target. Both diagrams edit the same path.`;
    if (path.kind === "edge") {
      const info=cases.asymptotic(path.index,...end), baseInfo=cases.asymptotic(path.index,...end,true);
      const deltas=activeCase.seed === "bouquet" ? ["1/2 − R","1/2 − S","R + S − 1/2"] : ["1 − R",`(${activeCase.top === .5 ? "1/2" : "1"}) − S`,activeCase.top === .5?"R + 2S − 1":"R + S − 1"];
      const chartName={raw:"τ",zero:"w = 1 − 1/τ",half:"w = τ/(1 − 2τ)",one:"w = 1/(1 − τ)"};
      $("limit-cusp-target").textContent=info?.known === false
        ? `Selected chart targets cusp ${info.cusp}; follow the resolved numerical trace for this mixed decoration. ${baseLabel}: cusp ${baseInfo?.cusp ?? "undetermined"}.`
        : `Selected case: cusp ${info?.cusp ?? "undetermined"} · ${baseLabel}: cusp ${baseInfo?.cusp ?? "undetermined"}.`;
      const coefficient=info?.logCoefficientLabel ?? (info ? ({[.25]:"1/(4π)",[.5]:"1/(2π)",[1]:"1/π",[2]:"2/π",[1/3]:"1/(3π)",[2/3]:"2/(3π)"}[info.alpha]??`${fmt(info.alpha)}/π`) : "");
      const factor=info?.KLabel ?? (info ? Math.abs(info.constantFactor-2*Math.sqrt(6))<1e-10?"2√6":Math.abs(info.constantFactor-2*Math.sqrt(3))<1e-10?"2√3":Math.abs(info.constantFactor-4*Math.sqrt(3))<1e-10?"4√3":fmt(info.constantFactor) : "");
      const known=info && info.known !== false && Number.isFinite(info.x0) && Number.isFinite(info.alpha) && (Number.isFinite(info.C)||Number.isFinite(info.logC));
      $("limit-cusp-detail").textContent=known
        ? `Proposed asymptotic: cusp coordinate ${chartName[info.chart]}; Re ${info.chart === "raw" ? "τ" : "w"} → ${info.xLabel ?? fmt(info.x0)}. Im ${info.chart === "raw" ? "τ" : "w"} = ${coefficient} log(C/δ) + o(1), where C = ${factor} ${activeCase.seed === "bouquet" ? "sin(πt)" : info.sineLabel} ≈ ${Number.isFinite(info.C)?fmt(info.C):`exp(${fmt(info.logC)})`}, and δ = ${info.deltaLabel ?? deltas[path.index]}. ${activeCase.seed === "bouquet" ? `The graph proposal uses M = ${info.M} and full marked phase p = ${info.p}. The constant ${factor} is the proposed K. Fix the decoration counts before taking the limit. ` : ""}This simulation reference is compared with the numerical data and applies away from vertices. On the displayed segment, δ is proportional to ε.`
        : (info?.detail ?? `The selected decorated branch is followed numerically. Its full edge asymptotic has not been established for this count combination; no exact horizontal limit is drawn. δ = ${deltas[path.index]} is proportional to ε along this path.`);
    } else {
      $("limit-cusp-target").textContent="Vertex approach · direction-dependent finite boundary arc.";
      $("limit-cusp-detail").textContent="The direction is determined by the starting circle. Edge cusp formulas do not apply to this vertex path; the selected finite endpoint is approached through numerical continuation.";
    }
    $("limit-draft-status").textContent=s.draftDirty?"Settings changed. Compute approach to update the plots below.":s.request?"":"Choose a start and boundary target, then compute the approach.";
    $("limit-draft-status").classList.toggle("settings-dirty",s.draftDirty);
    $("limit-results-state").hidden=!s.draftDirty||!s.request;
    $("limit-results-state").textContent="Previous computed path — the source settings have changed. Compute approach to update these plots.";
  }

  function run(event) {
    event?.preventDefault();
    if (!$("limits-controls").reportValidity()) return;
    const path=controlPath(), target="1e-"+$("limit-depth").value,density=Number($("limit-density").value);
    const priorView=extensionView({caseId:activeCase.id,caseConfig:activeCase.caseConfig,seed:activeCase.seed,path,target,startEpsilon:importedStartEpsilon,samplesPerDecade:density});
    saveFamilyResults();cancel();const id=++s.job;
    s.nodes=[];s.mesh=null;s.result=null;s.step=-1;s.view=null;s.draftDirty=false;
    s.detachedRestore=false;
    s.request={caseId:activeCase.id,caseConfig:activeCase.caseConfig,seed:activeCase.seed,decoration:activeCase.operation,weights:activeCase.weights,path,target,startEpsilon:importedStartEpsilon,samplesPerDecade:density};
    // Coordinate choice is a display setting; retain it on recomputation.
    $("limit-chart").value=s.chart;
    $("limits-applied").textContent=`${activeCase.title}; weights (${activeCase.weights.join(", ")}). Path ${path.kind}${path.index}, t = ${fmt(path.position)}; start (${path.start.join(", ")}); ε from ${importedStartEpsilon} to ${target}; ${density} samples per decade. The selected source path is evaluated; precision details are reported below.`;
    const url=new URL(location.href);
    for (const key of ["a","b","c","start","r","s"]) url.searchParams.delete(key);
    Object.entries({seed:activeCase.seed,decoration:activeCase.operation,path:path.kind+path.index,t:path.position,fromR:path.start[0],fromS:path.start[1],depth:$("limit-depth").value,density,chart:s.chart,startEpsilon:importedStartEpsilon}).forEach(([k,v])=>url.searchParams.set(k,v));
    url.searchParams.delete("sourceMode");
    url.searchParams.delete("autocompute");
    history.replaceState(null,"",url);
    window.AtlasUI?.refreshNavigation(url.searchParams);
    const resolution=activeCase.isRepeated ? 6 : 10;
    const meshKey=JSON.stringify([activeCase.id,activeCase.caseConfig,resolution,.04]),traceKey=JSON.stringify([activeCase.id,activeCase.caseConfig,path,target,density,importedStartEpsilon]);
    const sharedMesh=window.AtlasNumerics?.getMesh(graphDescriptor(s.request),{inset:.04});
    const savedMesh=sharedMesh?.result ?? readResultCache("mesh",meshKey),savedTrace=window.AtlasNumerics?.getTrace(graphDescriptor(s.request),s.request)?.result ?? readResultCache("trace",traceKey);
    const continuedResolution=['partial','unresolved'].includes(savedMesh?.status)&&Number.isInteger(savedMesh?.resolution)&&savedMesh.resolution>=3&&savedMesh.resolution<=24?savedMesh.resolution:resolution;
    s.meshMeta={resolution:sharedMesh?.resolution ?? resolution,inset:sharedMesh?.inset ?? .04};
    const extension=savedTrace?.status==="ok"?null:window.AtlasNumerics?.getExtension(graphDescriptor(s.request),s.request);
    s.extension=null;
    status("image","Computing the base and selected compact τ-images…");status("trace","Continuing both geometric branches toward the chosen boundary…");
    prepareExtension(extension,priorView);
    render();
    const baseMesh=window.AtlasNumerics?.getBaseMesh(graphDescriptor(s.request),{inset:.04});
    const reuseMesh=savedMesh?.samples?{...graphDescriptor(s.request),result:savedMesh}:null;
    if (savedMesh?.samples&&!['partial','unresolved'].includes(savedMesh.status)) {s.mesh=savedMesh;status("image",`Computed compact images reused (${savedMesh.samples.length} samples); no image computation needed.`);}
    else startWorker("mesh","repeated-worker.mjs?v=20261025-tv-maintenance1",{type:"mesh",id,caseId:activeCase.id,caseConfig:activeCase.caseConfig,R:s.start.r,S:s.start.s,resolution:continuedResolution,inset:.04,...(baseMesh?{baseMesh}:{}),...(reuseMesh?{reuseMesh}:{})},reply=>{
      if (reply.type === "progress") {const q=reply.result??reply;status("image",`Compact image samples: ${q.completed??"…"}/${q.total??"…"}.`);}
      else if (reply.type === "error") {settle("mesh");status("image",reply.message??reply.result?.message??"Image computation failed.",true);}
      else if (reply.type === "mesh" || reply.type === "done") {
        settle("mesh");s.mesh=reply.result??reply.mesh;s.meshMeta={resolution:s.mesh?.resolution??resolution,inset:s.mesh?.inset??.04};
        if(s.mesh?.samples){storeResultCache("mesh",meshKey,s.mesh);const count=s.mesh.samples.filter(x=>x.result?.status==="ok").length,total=s.mesh.total??s.mesh.samples.length;status("image",`Compact image: ${count}/${total} resolved points.${s.mesh.metrics?.reusedSamples?` Kept ${s.mesh.metrics.reusedSamples} previously checked samples.`:""}${s.mesh.status!=="ok"?" Partial image; only checked faces are drawn. Compute approach continues its missing samples.":""}${s.mesh.message?" "+s.mesh.message:""}`);}
        else status("image","The compact image could not be resolved.",true);
        if (s.fit === "image")s.view=null;render();saveFamilyResults();
      }
    });
    if (savedTrace?.status === "ok" && Array.isArray(savedTrace.nodes)){finishTrace(savedTrace);status("trace",`Exact source path reused from the shared cache: ${s.nodes.filter(good).length} samples.`);}
    else startWorker("trace","repeated-worker.mjs?v=20261025-tv-maintenance1",{type:"trace",id,caseId:activeCase.id,caseConfig:activeCase.caseConfig,R:s.start.r,S:s.start.s,path,target,startEpsilon:importedStartEpsilon,samplesPerDecade:density,budgetMs:45000,...(extension?{extendFrom:extension}:{})},reply=>{
      if (reply.type === "progress") {const q=reply.result??reply;if(!q.node)return;acceptTraceProgress(q);if(s.fit!=="image"&&Number($("limit-zoom").value)===0)s.view=null;render();}
      else if (reply.type === "error") {settle("trace");const message=reply.message??reply.result?.message??"Trace computation failed.";status("trace",message,true);traceError(message);}
      else if (reply.type === "trace" || reply.type === "done") {const result=reply.result;settle("trace",!!result?.resumeToken && !!result?.stop?.resumeAvailable);if(result?.nodes){if(result.status==="ok")storeResultCache("trace",traceKey,result);finishTrace(result);}else status("trace","No resolved trace was returned.",true);}
    });
    render();
  }

  function updateProgress(completed,total,epsilon){
    $("limit-job-progress").hidden=false;
    if(total>0){s.progressTotal=total;$("limit-progress").max=total;$("limit-progress").value=completed??0;}
    $("limit-progress-label").textContent=total>0?`${completed??0}/${total} checked samples${epsilon?`; last ε = ${epsilon}`:""}`:"Preparing the marked branch…";
    $("limit-result-summary").textContent=total>0?`${completed??0}/${total} checked samples${s.pending.has("trace")?" · computing…":""}`:"Preparing the marked branch…";
  }
  function renderInterruption(){
    const result=s.result,box=$("limit-interruption");if(!box)return;
    box.hidden=!result||result.status === "ok"||s.pending.has("trace");
    $("limit-recovery-link").hidden=box.hidden;
    if(result&&!s.pending.has("trace"))$("limit-result-summary").textContent=result.status === "ok"?`Trace complete · ${s.nodes.filter(good).length} checked samples.`:`Approach incomplete · ${s.nodes.filter(good).length}/${result.total||s.progressTotal||s.nodes.length} checked samples.`;
    if(box.hidden){$("limit-resume").hidden=true;return;}
    const info=window.AtlasUI?.explainStop(result)??{title:"Trace incomplete",explanation:"Checked samples have been retained.",advice:"Try a shorter approach."};
    $("limit-interruption-heading").textContent=info.title;
    $("limit-interruption-explanation").textContent=info.explanation;
    const portablePrefix=!s.draftDirty&&s.detachedRestore&&s.request&&window.AtlasNumerics?.getExtension(graphDescriptor(s.request),s.request);
    $("limit-interruption-advice").textContent=s.draftDirty?"Settings have changed. Compute a new approach to use them.":portablePrefix?"Saved checked points are restored. Use the main computation button to continue the remaining trace without recomputing that prefix.":s.detachedRestore?"Saved checked points are restored. The previous worker checkpoint ended when this page or family was left; Compute approach to retry this path.":info.advice;
    const node=s.nodes.findLast(good),target=s.request?.target??"—";
    $("limit-interruption-range").textContent=`Resolved ${s.nodes.filter(good).length}/${result.total||s.progressTotal||s.nodes.length} samples · last checked ε = ${node?.epsilon??"none"} · requested ε = ${target}${result.stop?.failedEpsilon?` · next unresolved ε = ${result.stop.failedEpsilon}`:""}.`;
    $("limit-interruption-detail").textContent=result.message??result.stopReason??"No additional solver message.";
    const resumable=!!result.resumeToken && !!result.stop?.resumeAvailable && !!s.workers.trace && !s.draftDirty;
    $("limit-resume").hidden=!resumable;$("limit-resume").disabled=s.pending.size>0;
    $("limit-retry").disabled=s.pending.size>0;
    $("limit-retry").textContent=Number($("limit-depth").value)<=5 && $("limit-density").value==="1"?"Retry these settings":"Try shorter preview";
  }
  function resumeTrace(){
    const job=s.traceJob,worker=s.workers.trace,token=s.result?.resumeToken;
    if(!job||!worker||!token||s.draftDirty||s.pending.size)return;
    pausePlayback();const id=++s.job;s.pending.add("trace");
    $("limits-run").disabled=true;$("limits-cancel").disabled=false;
    $("limit-interruption").hidden=true;$("limit-progress-label").textContent="Continuing from the saved marked-branch checkpoint…";
    status("trace","Continuing this exact source path; previous checked points are kept.");
    worker.postMessage({...job.data,id,resumeToken:token,resume:true});
  }
  $("limit-resume").onclick=resumeTrace;
  $("limit-retry").onclick=()=>{$("limit-depth").value=String(Math.min(5,Number($("limit-depth").value)));$("limit-density").value="1";changeDraft();run();};

  function applyCaseChange() {
    
    // A changed graph cannot reuse a continuation checkpoint, old samples,
    // image bounds or numerical references from its predecessor.
    saveFamilyResults();cancel();
    s.extension=null;
    s.nodes=[];s.mesh=null;s.meshMeta=null;s.result=null;s.request=null;s.traceJob=null;
    s.step=-1;s.view=null;s.zoomCorners=[];s.progressTotal=0;s.draftDirty=true;
    s.detachedRestore=false;
    refreshCaseLabels();refreshSamplingLabels();updateComputationSummary();
    $("limit-job-progress").hidden=true;
    $("limit-progress").removeAttribute("value");
    $("limit-progress-label").textContent="";
    $("limits-applied").textContent="Decoration applied. Compute approach to evaluate this graph along the source path shown above.";
    $("limit-result-summary").textContent="Decoration changed · ready to compute.";
    status("image","Decoration changed. Compute approach to rebuild the compact image.");
    status("trace","Source path and sampling settings retained. Press Compute approach for the selected decoration.");
    const restored=restoreFamilyResults();pathNote();render();
    if(restored)syncComputationURL();
  }

  function applyFamilyChange(event) {
    const detail=event.detail ?? {};
    // Invalidate every callback and continuation checkpoint before exposing
    // the new graph. Editing controls and canvases themselves stay mounted.
    familyDepths.set(selectedFamily,{depth:$("limit-depth").value,revision:depthRevision});
    cancel();
    resetImportedPath();
    selectedFamily=detail.seed ?? activeCase?.seed ?? window.AtlasFamily.defaultSeed;
    const previousDepth=$("limit-depth").value;
    const previousPath=$("limit-path").value;
    refreshFamilyControls();
    const rememberedDepth=familyDepths.get(selectedFamily);
    // A newer explicit depth choice takes precedence over an older family
    // snapshot; otherwise restore that family's previous setting.
    const preferredDepth=rememberedDepth?.revision===depthRevision ? rememberedDepth.depth : previousDepth;
    $("limit-depth").value=preferredDepth;
    expectedDepth=$("limit-depth").value;
    if (Array.isArray(detail.counts)) detail.counts.forEach((value,index)=>{$("limit-"+["a","b","c"][index]).value=String(value);});
    const source=detail.source;
    s.start=source && [source.r,source.s].every(Number.isFinite)
      ? {r:source.r,s:source.s} : compactPoint(...(activeCase?.start ?? [.31,.29]));
    s.extension=null;
    s.nodes=[];s.mesh=null;s.meshMeta=null;s.result=null;s.request=null;s.traceJob=null;
    s.step=-1;s.view=null;s.zoomCorners=[];s.progressTotal=0;s.draftDirty=true;
    s.detachedRestore=false;
    $("limit-job-progress").hidden=true;
    $("limit-progress").removeAttribute("value");$("limit-progress-label").textContent="";
    const adjustments=(previousPath!==$("limit-path").value ? " This family uses open edges; the former vertex target is now on the corresponding edge." : "")+
      "";
    $("limits-applied").textContent="Seed family changed. The source path is mapped to this triangle."+adjustments+" Press Compute approach when ready.";
    $("limit-result-summary").textContent="Family changed · ready to compute.";
    status("image","Family changed. Compute approach to build this compact image.");
    status("trace","Display settings retained."+adjustments+" Press Compute approach for this family.");
    restoreFamilyResults();updateComputationSummary();pathNote();render();
    syncComputationURL();
  }

  /* Shared canvas axes */

  function frame(canvas, bounds, xLabel, yLabel, tickFormat = null) {
    const width = Math.max(
      200,
      canvas.getBoundingClientRect().width
    );

    const height = canvas.id === "limit-lens" || canvas.id === "limit-lens-scale" ? LENS_HEIGHT : 340;
    const ratio = Math.min(devicePixelRatio || 1, 2);

    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    canvas.style.height = height + "px";

    const c = canvas.getContext("2d");

    c.setTransform(ratio, 0, 0, ratio, 0, 0);

    const isLens = canvas.id === "limit-lens" || canvas.id === "limit-lens-scale";
    const box = {
      x: isLens ? LENS_LEFT : 68,
      y: 18,
      w: width - (isLens ? LENS_LEFT + LENS_RIGHT : 88),
      h: isLens ? LENS_PLOT_HEIGHT : height - 68
    };

    const map = p => [
      box.x +
        (p[0] - bounds[0]) /
        (bounds[1] - bounds[0]) * box.w,

      box.y +
        (bounds[3] - p[1]) /
        (bounds[3] - bounds[2]) * box.h
    ];

    c.font = "11px system-ui";
    c.fillStyle = "#596e62";
    c.strokeStyle = "#d6ded2";

    const tickCount = isLens ? 2 : 4;
    for (let i = 0; i <= tickCount; i++) {
      const x =
        bounds[0] + (bounds[1] - bounds[0]) * i / tickCount;

      const y =
        bounds[2] + (bounds[3] - bounds[2]) * i / tickCount;

      const [px] = map([x, bounds[2]]);
      const [, py] = map([bounds[0], y]);

      c.beginPath();
      c.moveTo(px, box.y);
      c.lineTo(px, box.y + box.h);
      c.moveTo(box.x, py);
      c.lineTo(box.x + box.w, py);
      c.stroke();

      if (isLens && tickFormat) {
        c.save();
        c.translate(px, box.y + box.h + 13);
        c.rotate(-Math.PI / 4);
        c.textAlign = "right";
        c.font = "10px ui-monospace, monospace";
        c.fillText(tickFormat.x(x), 0, 0, 82);
        c.restore();
        c.textAlign = "right";
        c.font = "10px ui-monospace, monospace";
        c.fillText(tickFormat.y(y), box.x - 7, py + 4, box.x - 30);
        c.font = "11px system-ui";
      } else {
        c.textAlign = "center";
        c.fillText(fmt(x), px, box.y + box.h + 18);
        c.textAlign = "right";
        c.fillText(fmt(y), box.x - 7, py + 4);
      }
    }

    c.textAlign = "center";
    c.fillText(
      xLabel,
      box.x + box.w / 2,
      height - 6
    );

    c.save();
    c.translate(12, box.y + box.h / 2);
    c.rotate(-Math.PI / 2);
    c.fillText(yLabel, 0, 0);
    c.restore();

    c.save();
    c.beginPath();
    c.rect(box.x, box.y, box.w, box.h);
    c.clip();

    return {
      c,
      map,
      box,
      finish: () => c.restore()
    };
  }

  function drawSeries(plot, points, tint, selected = -1) {
    const {c, map} = plot;

    c.strokeStyle = tint;
    c.lineWidth = 1.8;
    c.beginPath();

    points.forEach((p, i) => {
      if (!p) return;

      if (i && points[i - 1]) {
        c.lineTo(...map(p));
      } else {
        c.moveTo(...map(p));
      }
    });

    c.stroke();
    c.fillStyle = tint;

    points.forEach((p, i) => {
      if (!p) return;

      c.beginPath();

      c.arc(
        ...map(p),
        i === selected ? 4.5 : 2.5,
        0,
        2 * Math.PI
      );

      c.fill();
    });
  }

  /* Magnified image and trajectory */

  function renderLens() {
    const chosen = s.nodes[s.step];
    const trace = s.nodes.map(node => good(node) ? values(node) : null);
    const mesh = s.mesh?.samples.map(sample => meshValue(sample)) ?? [];
    const baseMesh = $("limit-compare").checked
      ? s.mesh?.samples.map(sample => meshValue(sample, "base")) ?? [] : [];
    const available = [...mesh, ...trace].filter(Boolean);
    const width = Math.max(200, $("limit-lens-scale").getBoundingClientRect().width);
    const aspect = Math.max(1, width - LENS_LEFT - LENS_RIGHT) / LENS_PLOT_HEIGHT;
    const ref = infinityReference();
    const resolveHorizontal = $("limit-resolve-horizontal").checked;
    $("limit-resolve-wrap").hidden = false;
    $("limit-horizontal-view").hidden = !resolveHorizontal;

    if (!s.view) {
      const origin = values(chosen) ?? trace.find(Boolean) ?? available[0] ?? ["0", "0"];
      let fitted = s.fit === "image" ? mesh.filter(Boolean) : trace.filter(Boolean);
      if (s.fit === "tail") {
        fitted = s.nodes.slice(Math.max(0, s.step - 3), s.step + 1)
          .filter(good).map(node => values(node)).filter(Boolean);
      }
      if (!fitted.length) fitted = available;
      const offsets = fitted.map(z => offset(z, origin));
      const xs = offsets.map(p => p[0]), ys = offsets.map(p => p[1]);
      // A rational reference is independent of the computed points. Include
      // it in the horizontal fit without inventing a trajectory endpoint.
      if (ref && s.fit !== "image" && offsets.length) {
        xs.push(difference(fractionDecimal(ref.numerator, ref.denominator), origin[0]));
      }
      const x0 = xs.length ? Math.min(...xs) : -.5;
      const x1 = xs.length ? Math.max(...xs) : .5;
      const y0 = ys.length ? Math.min(...ys) : -.5;
      const y1 = ys.length ? Math.max(...ys) : .5;
      const span = Math.max(y1-y0, (x1-x0)/aspect)*1.2 || .001;
      s.view = {
        origin, x:(x0+x1)/2, y:(y0+y1)/2, span,
        spanX:offsets.length ? Math.min(span*aspect,Math.max((x1-x0)*1.3, 1e-280)) : span*aspect
      };
    }

    const v = s.view;
    const zoom = 10 ** Number($("limit-zoom").value);
    const span = v.span / zoom;
    const normalSpanX = span*aspect;
    const detailSpanX = (v.spanX ?? v.span*aspect) / zoom;
    const variable = s.chart === "raw" ? "τ" : "w";
    const boundsFor = spanX => [v.x-spanX/2, v.x+spanX/2, v.y-span/2, v.y+span/2];
    const normalBounds = boundsFor(normalSpanX), detailBounds = boundsFor(detailSpanX);
    const activeBounds = resolveHorizontal ? detailBounds : normalBounds;

    // The overview is only an ordinary-precision locator. Both local plots
    // retain the decimal origin and subtract it before converting offsets.
    s.zoomCorners = [
      [activeBounds[0],activeBounds[2]], [activeBounds[1],activeBounds[2]],
      [activeBounds[1],activeBounds[3]], [activeBounds[0],activeBounds[3]]
    ].map(([x,y]) => inverseChart(Number(v.origin[0])+x, Number(v.origin[1])+y))
      .filter(z => z && Number.isFinite(z.re) && Number.isFinite(z.im));

    function drawLocal(canvas, bounds, spanX, equalUnits) {
      canvas.dataset.coordinateSystem = variable;
      canvas.dataset.horizontalSpan = String(spanX);
      canvas.dataset.verticalSpan = String(span);
      canvas.dataset.unitScale = equalUnits ? "equal" : "magnified";
      const plot = frame(canvas, bounds, "Re "+variable, "Im "+variable, {
        x:x => absoluteTick(v.origin[0],x,spanX),
        y:y => absoluteTick(v.origin[1],y,span)
      });
      const showMesh = Math.min(span,spanX) >= 1e-6;
      const references = ref ? [
        ...($("limit-compare").checked && infinityReference(true) ? [{decimal:fractionDecimal(infinityReference(true).numerator,infinityReference(true).denominator),tint:baseColor,label:`${infinityReference(true).label} (${baseLabel})`}] : []),
        {decimal:fractionDecimal(ref.numerator,ref.denominator),tint:color,label:ref.label,kind:(ref.kind??"Proposed limit")}
      ] : [];
      references.forEach((reference,i) => {
        const x = difference(reference.decimal,v.origin[0]);
        const px = plot.map([x,bounds[3]])[0];
        plot.c.save();
        plot.c.strokeStyle=reference.tint; plot.c.lineWidth=1.3;
        plot.c.setLineDash([5,4]); plot.c.beginPath();
        plot.c.moveTo(...plot.map([x,bounds[2]]));
        plot.c.lineTo(...plot.map([x,bounds[3]])); plot.c.stroke();
        if (px >= plot.box.x && px <= plot.box.x+plot.box.w) {
          const label=`${reference.kind ? reference.kind+": " : ""}Re ${variable} → ${reference.label}`;
          plot.c.font="bold 11px system-ui";
          const labelWidth=plot.c.measureText(label).width;
          const labelX=Math.max(plot.box.x+5,Math.min(plot.box.x+plot.box.w-labelWidth-5,px+6));
          const labelY=plot.box.y+17+i*20;
          plot.c.fillStyle="#fcfcf7ee";plot.c.fillRect(labelX-3,labelY-13,labelWidth+6,18);
          plot.c.fillStyle=reference.tint;plot.c.textAlign="left";
          plot.c.fillText(label,labelX,labelY);
        }
        plot.c.restore();
      });
      for (const layer of [{data:baseMesh,tint:baseColor},{data:mesh,tint:color}]) {
        if (!s.mesh || !showMesh || !layer.data.length) continue;
        plot.c.fillStyle=layer.tint+"28";plot.c.strokeStyle=layer.tint+"55";plot.c.lineWidth=.5;
        for (const face of s.mesh.triangles) {
          const points=face.map(i => layer.data[i] ? plot.map(offset(layer.data[i],v.origin)) : null);
          if (points.some(p => !p)) continue;
          plot.c.beginPath();plot.c.moveTo(...points[0]);
          plot.c.lineTo(...points[1]);plot.c.lineTo(...points[2]);
          plot.c.closePath();plot.c.fill();plot.c.stroke();
        }
      }
      if ($("limit-compare").checked) drawSeries(plot,s.nodes.map(node => {
        const z=node.base?.status === "ok" ? values(node,"base") : null;
        return z ? offset(z,v.origin) : null;
      }),baseColor,s.step);
      drawSeries(plot,trace.map(z => z ? offset(z,v.origin) : null),color,s.step);
      plot.finish();
      return showMesh;
    }

    const normalMesh = drawLocal($("limit-lens-scale"),normalBounds,normalSpanX,true);
    const detailMesh = resolveHorizontal
      ? drawLocal($("limit-lens"),detailBounds,detailSpanX,false) : normalMesh;
    $("limit-zoom-value").textContent=fmt(zoom)+"×";
    const selected=good(chosen) ? values(chosen) : null;
    const error=chosen?.result?.viewErrors?.[s.chart];
    $("limit-rational-target").hidden=!ref && !(((s.request?.path??controlPath()).kind==="edge"));
    if (ref) {
      const delta=selected ? Math.abs(difference(selected[0],fractionDecimal(ref.numerator,ref.denominator))) : null;
      const referenceLabel=(ref.kind??"Proposed graph limit");
      const baseRef=$("limit-compare").checked?infinityReference(true):null;
      const methodOpen=$("limit-rational-target").querySelector("details")?.open??false;
      $("limit-rational-target").dataset.predictionStatus=ref.status??"predicted";
      $("limit-rational-target").innerHTML=`<span class="tv-prediction-kind">${escape(referenceLabel)}</span> <span class="tv-rational-number">Re ${variable} → ${window.TVCuspPredictions?.fractionHTML(ref)??escape(ref.label)}</span>`+
        `<span class="tv-prediction-detail">${delta!=null?`Distance of selected sample from the reference line: ${escape(fmt(delta,5))}. `:""}${baseRef?`${escape(baseLabel)}: ${escape(baseRef.label)}. `:""}</span>`+
        ((`<details class="tv-reference-method"${methodOpen?" open":""}><summary>About this proposal</summary><p>${escape(ref.detail??"Simulation reference; no exact-limit proof is claimed.")}${Number.isFinite(ref.uncertainty)?` Tail-fit variation: ${escape(fmt(ref.uncertainty,4))} (a diagnostic, not an error bound).`:""}</p></details>`));
    } else {
      const p=s.request?.path??controlPath(),end=boundaryPoint(p),info=cases.asymptotic(p.index,...end);
      $("limit-rational-target").dataset.predictionStatus="unresolved";
      $("limit-rational-target").textContent=info?.chart!==s.chart
        ? `To inspect the horizontal rational limit, select the cusp-opening coordinate ${info?.chartLabel??"for this edge"}. The raw finite cusp is a different quantity.`
        : "Rational candidate not resolved yet. Continue the checked trace deeper; a candidate is displayed only when consecutive samples stabilize near the same simple rational.";
    }
    $("limit-lens-coordinate").textContent=selected
      ? `Selected ${variable} = ${fmt(selected[0],14)} + ${fmt(selected[1],14)} i`
      : `Select a resolved sample to display ${variable}.`;
    $("limit-lens-coordinate").title=selected
      ? `Re ${variable} = ${selected[0]}; Im ${variable} = ${selected[1]}` : "";
    if ($("limit-lens-scale-caption")) $("limit-lens-scale-caption").textContent=
      `Equal units on both axes. Width ${fmt(normalSpanX)}; height ${fmt(span)}. This view preserves the direction of approach.`;
    if ($("limit-horizontal-note")) $("limit-horizontal-note").textContent=
      `Same vertical range; horizontal scale ${fmt(normalSpanX/detailSpanX,4)}× the equal-units view. Width ${fmt(detailSpanX)}. Drag this graph to pan both views.`;
    const activeSpanX=resolveHorizontal ? detailSpanX : normalSpanX;
    $("limit-lens-note").textContent=
      `Axes show absolute ${variable} coordinates; labels are rounded to the viewing scale. `+
      (!detailMesh ? "Only the checked trajectory is shown below the ordinary-precision mesh scale. " : "")+
      (error && Math.min(span,activeSpanX)<10*Number(error)
        ? "This view is below the numerical precision diagnostic; further zoom adds no verified detail." : "");
    renderOverview();
  }

  /* Behaviour against decreasing epsilon */

  function renderTrend() {
    const chosen = s.nodes[lastGood()];
    const metric = $("limit-metric").value;
    const origin = values(chosen) ?? ["0", "0"];

    function series(branch) {
      return s.nodes.map((node, i) => {
        const z = values(node, branch);

        if (!good(node) || !z) return null;

        let y;

        if (metric === "change") {
          const prior = i
            ? values(s.nodes[i - 1], branch)
            : null;

          if (!prior || !good(s.nodes[i - 1])) {
            return null;
          }

          const distance = Math.hypot(...offset(z, prior));

          if (!(distance > 0)) return null;

          y = Math.log10(distance);
        } else {
          const index = metric === "real" ? 0 : 1;
          y = difference(z[index], origin[index]);
        }

        return [
          -Math.log10(Number(node.epsilon)),
          y
        ];
      });
    }

    const result = series("result");

    const base = $("limit-compare").checked
      ? series("base")
      : [];

    const points = [...result, ...base].filter(Boolean);
    const ys = points.map(p => p[1]);

    let low = ys.length ? Math.min(...ys) : 0;
    let high = ys.length ? Math.max(...ys) : 1;

    const pad =
      (high - low) * 0.12 ||
      Math.max(Math.abs(low) * 0.1, 1e-12);

    low -= pad;
    high += pad;

    const right = Math.max(
      3,
      ...points.map(p => p[0])
    );

    const variable = s.chart === "raw" ? "τ" : "w";

    const label = metric === "change"
      ? `log₁₀ |Δ${variable}|`
      : (metric === "real" ? "Δ Re " : "Δ Im ") + variable;

    const plot = frame(
      $("limit-trend"),
      [1, right, low, high],
      "−log₁₀ ε",
      label
    );

    drawSeries(plot, base, baseColor, s.step);
    drawSeries(plot, result, color, s.step);

    plot.finish();

    $("limit-trend-note").textContent =
      metric === "change"
        ? "Consecutive computed changes, on a logarithmic vertical " +
          "scale. These are not errors from an assumed limit."
        : "Vertical origin is the last checked result: " +
          `${metric === "real" ? "Re" : "Im"} ${variable}* = ` +
          `${origin[metric === "real" ? 0 : 1]}. ` +
          `Orange: result; teal: optional ${baseLabel} comparison.`;
  }

  function coordinateCheck(node) {
    const value=node.result?.viewErrors?.[s.chart] ?? (s.chart === "raw" ? node.result?.crossPrecisionTauError : null);
    if(value == null || !Number.isFinite(Number(value))) return null;
    return Math.abs(Number(value));
  }

  function leadingConstantComparison() {
    if((activeCase.seed!=="bouquet") || s.draftDirty || s.request?.caseId!==activeCase.id)return "";
    const path=s.request?.path,config=s.request?.caseConfig;
    if(path?.kind!=="edge"||!config)return "";
    const formula=window.TVCuspFormulas?.lookup(config,path.index);
    if(!formula||!Number.isFinite(formula.logConstantFactor))return "";
    const node=s.nodes.findLast(good),z=node?.result?.views?.[formula.chart]??(formula.chart==="raw"?node?.result?.tauExact:null);
    if(!z)return "";
    const start=(path.start??[]).map(Number),initial=[.5-start[0],.5-start[1],start[0]+start[1]-.5][path.index];
    const epsilon=Number(node.epsilon),sine=Math.sin(Math.PI*Number(path.position)),height=Number(z[1]);
    if(!(initial>0&&epsilon>0&&sine>0&&height>0))return "";
    // This is a post-computation data check. The straight path has exactly
    // delta=initial*epsilon; subtracting Number(R) from1/2 would erase a cusp.
    const logDelta=Math.log(initial)+Math.log(epsilon),logScale=Math.PI*formula.M*height;
    const logRatio=logDelta+logScale-formula.logConstantFactor-Math.log(sine);
    const repeat=Number(node.result?.viewErrors?.[formula.chart]??(formula.chart==="raw"?node.result?.crossPrecisionTauError:0))||0;
    const resolution=Math.max(1e-10,Math.PI*formula.M*Math.abs(repeat),32*Number.EPSILON*(Math.abs(logDelta)+Math.abs(logScale)+Math.abs(formula.logConstantFactor)+1));
    const relative=Math.expm1(logRatio),ratio=Math.exp(logRatio);
    const agreement=Math.abs(logRatio)<=resolution
      ? "Khat/K agrees with 1 within practical plot precision"
      : Number.isFinite(relative)&&Number.isFinite(ratio)
        ? `Khat/K ≈ ${fmt(ratio,6)}; relative discrepancy ${fmt(Math.abs(relative),4)}`
        : "the leading-constant estimate is outside the practical plotting range";
    return ` Cusp-scale check at ε = ${node.epsilon}: ${agreement}. Khat = δ exp(π M Im ${formula.chart==="raw"?"τ":"w"})/sin(πt), evaluated in this edge's cusp-opening coordinate. This is leading-term agreement, not numerical root error.`;
  }

  function renderEvidence() {
    const variable=s.chart === "raw" ? "τ" : "w", ref=infinityReference();
    const all=s.nodes.map((node,index)=>({node,index,z:good(node)?values(node):null})).filter(x=>x.z);
    const rows=all.slice(-3), reference=ref?fractionDecimal(ref.numerator,ref.denominator):null;
    $("limit-evidence-coordinate").textContent=`Re ${variable}`;
    $("limit-evidence-distance").textContent=ref ? "Distance to reference" : "Change from prior sample";
    $("limit-evidence-reference").textContent=ref
      ? `${(ref.kind??"Proposed graph limit")}: Re ${variable} → ${ref.label}`
      : "No reference limit asserted in this coordinate";
    $("limit-evidence-rows").innerHTML=rows.map(({node,index,z})=>{
      const prior=all[all.findIndex(x=>x.index===index)-1], check=coordinateCheck(node);
      const distance=ref ? Math.abs(difference(z[0],reference)) : prior ? Math.abs(difference(z[0],prior.z[0])) : null;
      const checkText=check==null ? "Unavailable" : check===0 ? "Below comparison scale" : fmt(check,3);
      const exact=String(z[0]).replace(/^(-?)\./,(_match,sign)=>sign+"0.");
      const shown=exact.length>30 ? exact.slice(0,29)+"…" : exact;
      return `<tr class="${index===s.step ? "selected" : ""}"><td><button type="button" data-step="${index}">${escape(node.epsilon)}</button></td><td class="stored-real" title="${escape(z[0])}">${escape(shown)}</td><td>${distance==null ? "—" : distance===0 ? "Below stored scale" : fmt(distance,5)}</td><td>${checkText}</td></tr>`;
    }).join("");
    if(!rows.length){$("limit-evidence-note").textContent="Compute a trace to compare the last three resolved samples. The two local graphs keep absolute coordinates; this table makes their small horizontal motion readable.";return;}
    let note=`Last ${rows.length} resolved samples of the selected graph, in the coordinate used by both local plots. `;
    if(ref && rows.length>1){
      const gaps=rows.map(x=>Math.abs(difference(x.z[0],reference))), last=gaps.at(-1), check=coordinateCheck(rows.at(-1).node);
      const decreasing=gaps.every((x,i)=>i===0||x<=gaps[i-1]);
      note+=check!=null && (check===0 ? last===0 : last<=10*check)
        ? "The final horizontal distance is at or below the comparison scale; do not infer more digits from this run. "
        : decreasing ? "These samples move closer to the reference line. This is numerical evidence; extend the depth if you need a later regime. "
        : "The last distances are not steadily decreasing. Check a deeper approach or a different path before inferring a limit. ";
    } else note+="Consecutive changes are empirical observations, not an asserted limit formula. ";
    if(s.result?.status === "partial") note+="The full approach is incomplete; some intermediate samples or branch checks remain unresolved. ";
    if(s.draftDirty) note+="These values belong to the previous computed path.";
    $("limit-evidence-note").textContent=note+leadingConstantComparison();
  }

  /* Sample selection and table */

  function render() {
    renderSource();
    renderLens();
    renderTrend();
    renderEvidence();
    renderInterruption();

    const node = s.nodes[s.step];

    $("limit-step").max = Math.max(0, s.nodes.length - 1);
    $("limit-step").value = Math.max(0, s.step);
    $("limit-step").disabled = !s.nodes.length;
    $("limit-play").disabled = s.pending.has("trace") || s.nodes.filter(good).length < 2;

    const checkedCount=s.nodes.filter(good).length,hasGaps=checkedCount<s.nodes.length;
    $("limit-step-label").textContent = node
      ? `${hasGaps?"Sample slot":"Sample"} ${s.step + 1}/${s.nodes.length}${hasGaps?` · ${checkedCount} checked`:""}; ε = ${node.epsilon}${good(node) ? "" : " · unresolved"}`
      : "No checked sample selected";

    const selectedCoordinates = good(node) ? values(node) : null;
    const selectedVariable = s.chart === "raw" ? "τ" : "w";
    $("limit-selected").textContent = selectedCoordinates
      ? `${selectedVariable} = (${selectedCoordinates[0]}, ${selectedCoordinates[1]}). `+
        "Numerical check in this coordinate: "+
        `${coordinateCheck(node)==null ? "no independent comparison available" : coordinateCheck(node)===0 ? "below the stored comparison scale" : fmt(coordinateCheck(node),5)}. Not a certified bound.`
      : node ? `This sample is unresolved; no ${selectedVariable} coordinate is plotted.` : "";

    $("limit-table").innerHTML = s.nodes.map((n, i) => {
      const resolved = good(n);

      return `
        <tr class="${i === s.step ? "selected" : ""}">
          <td>
            <button type="button"
                    data-step="${i}"
                    ${resolved ? "" : "disabled"}>
              ${escape(n.epsilon)}
            </button>
          </td>
          <td>
            ${resolved ? fmt(n.result.tauExact[0], 14) : "Unresolved"}
          </td>
          <td>
            ${resolved ? fmt(n.result.tauExact[1], 14) : "—"}
          </td>
          <td>
            ${
              resolved
                ? escape(fmt(n.result.viewErrors?.raw ?? n.result.crossPrecisionTauError))
                : "—"
            }
          </td>
        </tr>
      `;
    }).join("");

    $("limit-export").disabled = !s.nodes.length;
  }

  /* Zoom, pan, and selection */

  function pan(dx, dy, equalScale = false) {
    if (!s.view) return;
    const canvas = equalScale || !$("limit-resolve-horizontal").checked
      ? $("limit-lens-scale") : $("limit-lens");
    const zoom=10 ** Number($("limit-zoom").value);
    const span=s.view.span/zoom;
    const spanX=Number(canvas.dataset.horizontalSpan);
    if (!Number.isFinite(spanX) || !(spanX>0)) return;
    s.view.x += dx*spanX;
    s.view.y += dy*span;
    renderLens();
  }

  function centerSelected() {
    if (!s.view) renderLens();

    const z = values(s.nodes[s.step]);

    if (z && s.view) {
      [s.view.x, s.view.y] = offset(z, s.view.origin);
    }

    renderLens();
  }

  function selectStep(index) {
    if(s.extension)s.extension.selectedEpsilon=s.nodes[index]?.epsilon;
    s.step = index;

    if (s.fit === "tail") s.view = null;

    render();

    if (Number($("limit-zoom").value) > 0) {
      centerSelected();
    }
  }

  $("limits-controls").addEventListener("submit", run);
  $("limits-cancel").onclick = ()=>cancel(true);
  ["limit-a", "limit-b", "limit-c", "limit-position", "limit-depth", "limit-density"].forEach(id => {
    $(id).addEventListener("input", changeDraft);
    $(id).addEventListener("change", changeDraft);
  });
  $("limit-path").onchange = changeDraft;
  connectSourceEditor($("limit-source"), false);
  connectSourceEditor($("limit-source-attached"), true);
  $("limit-play").onclick = startPlayback;

  $("limit-chart").onchange = () => {
    s.chart = $("limit-chart").value;
    s.view = null;
    $("limit-zoom").value = 0;
    const url=new URL(location.href);
    url.searchParams.set("chart",s.chart);
    history.replaceState(null,"",url);
    window.AtlasUI?.refreshNavigation(url.searchParams);
    render();
  };

  $("limit-zoom").oninput = centerSelected;
  $("limit-center").onclick = centerSelected;
  $("limit-metric").onchange = renderTrend;
  $("limit-compare").onchange = () => {renderLens();renderTrend();};
  $("limit-resolve-horizontal").onchange = renderLens;

  document.querySelectorAll("[data-fit]").forEach(button => {
    button.onclick = () => {
      s.fit = button.dataset.fit;
      s.view = null;
      $("limit-zoom").value = 0;
      render();
    };
  });

  document.querySelectorAll("[data-pan]").forEach(button => {
    button.onclick = () => {
      pan(
        ...button.dataset.pan
          .split(",")
          .map(x => Number(x) * 0.15)
      );
    };
  });

  $("limit-step").oninput = () => {
    pausePlayback();
    selectStep(Number($("limit-step").value));
  };

  $("limit-evidence-rows").onclick = event => {
    const button=event.target.closest("button[data-step]");
    if(button){pausePlayback();selectStep(Number(button.dataset.step));}
  };

  $("limit-table").onclick = event => {
    const button = event.target.closest("button[data-step]");

    if (button) {
      pausePlayback();
      selectStep(Number(button.dataset.step));
    }
  };

  for (const id of ["limit-lens-scale","limit-lens"]) {
    const canvas=$(id);
    let drag=null;
    const equalScale=id === "limit-lens-scale";
    canvas.onpointerdown=event => {
      canvas.setPointerCapture(event.pointerId);
      drag={x:event.clientX,y:event.clientY,id:event.pointerId};
    };
    canvas.onpointermove=event => {
      if (!drag || drag.id!==event.pointerId) return;
      const width=Math.max(200,canvas.getBoundingClientRect().width);
      pan((drag.x-event.clientX)/(width-LENS_LEFT-LENS_RIGHT),
        (event.clientY-drag.y)/LENS_PLOT_HEIGHT,equalScale);
      drag={x:event.clientX,y:event.clientY,id:event.pointerId};
    };
    canvas.onpointerup=canvas.onpointercancel=()=>{drag=null;};
    canvas.onlostpointercapture=()=>{drag=null;};
    canvas.onkeydown=event => {
      const delta={ArrowLeft:[-.15,0],ArrowRight:[.15,0],ArrowUp:[0,.15],ArrowDown:[0,-.15]}[event.key];
      if (delta) {event.preventDefault();pan(...delta,equalScale);}
    };
  }

  /* Export the stored decimal coordinates and diagnostics */

  $("limit-export").onclick = () => {
    const data = {
      ...s.request,
      ...(s.result ?? {}),
      status:
        s.result?.status ??
        (s.pending.has("trace") ? "in_progress" : "stopped"),
      nodes: s.nodes,
      displayChart: s.chart,
      note:
        ($("limit-precision-note").textContent || "Adaptive numerical estimates; no interval certificate.")
    };

    const blob = new Blob(
      [JSON.stringify(data, null, 2)],
      {type: "application/json"}
    );

    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");

    link.href = url;
    link.download = "hemisphere-limit-trace.json";
    link.click();

    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };

  window.addEventListener("resize", () => {
    s.view = null;
    render();
  });
  window.addEventListener("atlas:casechange",applyCaseChange);
  window.addEventListener("atlas:beforefamilychange",event=>{trackDepthChoice();saveFamilyResults();event.detail.source={r:s.start.r,s:s.start.s};});
  window.addEventListener("atlas:familychange",applyFamilyChange);
  window.addEventListener("atlas:cachedpage", () => {
    const selectedEpsilon=s.nodes[s.step]?.epsilon,view=s.view;
    restoreSharedMesh();
    const restored=restoreExactSharedTrace();
    if(restored && selectedEpsilon!==undefined){
      const index=s.nodes.findIndex(node=>node.epsilon===selectedEpsilon && good(node));
      if(index>=0){s.step=index;s.view=view;}
    }
    pathNote();render();
  });

  $("limits-quantities").ontoggle = () => { if ($("limits-quantities").open) renderTrend(); };
  window.addEventListener("atlas:beforepagechange", () => {
    saveDisplaySettings();saveFamilyResults();
  });
  window.addEventListener("pagehide", () => {
    saveDisplaySettings();saveFamilyResults();
    const wasRunning=s.pending.has("trace");
    cancel();
    if(wasRunning){
      s.result={status:"partial",nodes:s.nodes,completed:s.nodes.filter(good).length,
        total:s.progressTotal||s.nodes.length,message:"The page was left while computing. Checked samples were retained.",stop:{kind:"page_change",resumeAvailable:false}};
      s.detachedRestore=true;
    }
    render();
  });

  initializeCase();
  
  refreshFamilyControls();
  restoreDisplaySettings();
  expectedDepth=$("limit-depth").value;
  pathNote();
  render();
  const savedPath = (query.has("fromR") && query.has("fromS")) ||
    (query.has("start") && query.has("path"));
  // Navigation links are authoritative about source paths. A remembered family
  // view must not replace an explicitly imported start or boundary selection.
  const explicitPath=savedPath || query.has("path") || query.has("t");
  const restoredFamily=!explicitPath && restoreFamilyResults();
  const restoredMesh=s.mesh ? true : restoreSharedMesh();
  const restoredTrace=s.nodes.length ? true : restoreExactSharedTrace();
  if(restoredFamily || restoredMesh || restoredTrace){pathNote();render();}
  if (restoredTrace) return;
  if (savedPath && query.get("autocompute")!=="0" && Number($("limit-depth").value)<=20) {
    const retainedZoom=$("limit-zoom").value;
    run();
    $("limit-zoom").value=retainedZoom;
    render();
  }
  else if(savedPath) $("limit-draft-status").textContent=Number($("limit-depth").value)>20 ? "Research settings loaded. Press Compute approach to start this deeper calculation." : "Saved source path loaded. Press Compute approach when ready.";
};
const bootLimits=()=>{
  try{initializeLimits();}
  catch(error){
    console.error("Limits controls could not initialize",error);
    const status=document.getElementById("limit-draft-status");
    if(status){status.textContent="The page controls could not load. Reload this page to try again.";status.dataset.initializationError=String(error?.message??error);}
  }
};
if(window.AtlasResults?.ready?.then)window.AtlasResults.ready.then(bootLimits,bootLimits);
else bootLimits();


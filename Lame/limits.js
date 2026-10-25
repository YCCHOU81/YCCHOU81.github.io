"use strict";

const initializeLimits = () => {
  const $ = id => document.getElementById(id);

  if (!$("limits-controls")) return;

  const baseLabel = "n = 1";
  const depthAdvice = $("limit-depth-advice").textContent;

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
    start: ({r: .31, s: .29})
  };

  // Every numerical cache shares one memory budget. In particular, repeated TV
  // computations cannot accumulate in an unbounded local Map or web storage.
  const cacheNamespace = kind => "limits-" + kind + "-v1";
  const DISPLAY_RESULTS = "limits-display-v1";
  function graphDescriptor(request = null) {
    return {family:"lame", counts:request ? request.counts : draftCounts()};
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

  const tailRefinementFor=total=>importedTailRefinement??(total>=8?2:total>=4?1:0);
  function requestForControls() {
    const path=controlPath(),target="1e-"+$("limit-depth").value;
    const sampling={path,target,startEpsilon:importedStartEpsilon,samplesPerDecade:Number($("limit-density").value)};
    
    const counts=draftCounts(),total=counts.reduce((a,b)=>a+b,0);
    return {counts,...sampling,tailRefinement:tailRefinementFor(total)};
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
      if(index>=0){s.step=index;s.view=null;render();if(Number($("limit-zoom").value)>0)centerSelected();saveCaseResults();}
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
    const identity = graphDescriptor(request);
    return window.AtlasResults?.key(identity) ?? JSON.stringify(identity);
  }

  function saveCaseResults() {
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

  function restoreCaseResults() {
    const saved = window.AtlasResults?.get(FAMILY_RESULTS, graphResultKey());
    if (!saved || (!saved.mesh && !saved.trace?.nodes?.length)) return false;
    s.mesh = saved.mesh ?? null;
    s.meshMeta=saved.meshMeta ?? null;
    s.view = null; s.zoomCorners = [];
    if (s.mesh) status("image", `${s.mesh.status==="partial"?"Partial interior":"Saved interior"}: ${s.mesh.samples?.length??0} points.`);
    const trace = saved.trace, request = trace?.request;
    const path = request?.path;
    const pathName = path ? path.kind + path.index : "";
    const samplingMatches = request?.target === "1e-" + $("limit-depth").value &&
      Number(request.samplesPerDecade) === Number($("limit-density").value);
    const validStart=Array.isArray(path?.start) && path.start.length===2 && path.start.every(x=>Number.isFinite(Number(x)));
    const validRay=path && !path.start && Number(request?.startEpsilon)>0;
    const validPath = path && [...$("limit-path").options].some(option => option.value === pathName) && (validStart||validRay);
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
    importedLegacyRay=validRay;importedStartEpsilon=request.startEpsilon??"1";
    importedTailRefinement=[0,1,2].includes(request.tailRefinement)?request.tailRefinement:null;
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
    status("trace", `${complete ? "Saved trace" : "Partial trace"}: ${s.nodes.filter(good).length} checked samples.` +
      (complete ? "" : " A worker checkpoint is not retained after switching families; Compute approach to retry."));
    $("limits-applied").textContent = `Saved ${complete ? "computed" : "partial"} family view restored. Path ${pathName}, t = ${path.position}; start (${actualStart.join(", ")}); ε from ${request.startEpsilon ?? "1"} to ${request.target}; ${request.samplesPerDecade} samples per decade. Display settings are retained.`;
    return true;
  }

  const query = new URLSearchParams(location.search);
  let importedTailRefinement=["0","1","2"].includes(query.get("tailRefinement"))?Number(query.get("tailRefinement")):null;
  let importedLegacyRay=(query.get("sourceMode")==="ray");
  let importedStartEpsilon=Number(query.get("startEpsilon"))>0 && Number(query.get("startEpsilon"))<=1
    ? query.get("startEpsilon") : importedLegacyRay ? query.get("start")||".04" : "1";

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
      option.disabled=(Number(option.value)>10);
      if (Number(option.value)>10) option.textContent=`10⁻${option.value.replace(/\d/g,n=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(n)])}${" (TV families only)"}`;
    }
  }
  ensureCoordinateOptions();

  if (["3", "5", "7", "8", "10"].includes(query.get("depth"))) {
    $("limit-depth").value = query.get("depth");
  }

  if (query.has("t")) {
    $("limit-position").value = query.get("t");
  }

  if (["1", "2", "3"].includes(query.get("density"))) $("limit-density").value = query.get("density");

  if(importedLegacyRay) {
    const legacy=legacySourcePoint(controlPath(),importedStartEpsilon);
    s.start={r:legacy[0],s:legacy[1]};
  } else if (query.has("fromR") && query.has("fromS")) {
    s.start = importedStart(Number(query.get("fromR")), Number(query.get("fromS")));
  } else if (query.has("start") && Number(query.get("start")) > 0 && Number(query.get("start")) <= .3) {
    const legacy = legacySourcePoint(controlPath(), query.get("start"));
    s.start = importedStart(...legacy);
  }

  // Mixed non-connector graphs need an expensive branch continuation even at
  // small weights. Start with a shorter path unless a bookmarked depth says otherwise.
  
  // Display sampling is separate from adaptive branch continuation and precision.
  let presets = {preview: [5, 1], compare: [10, 1], deeper: [10, 3]};
  function updateComputationSummary() {
    const depth=Number($("limit-depth").value), density=Number($("limit-density").value);
    const sum=["a","b","c"].reduce((sum,key)=>sum+Number($("limit-"+key).value||0),0);
    const extra=tailRefinementFor(sum);
    const match=Object.entries(presets).find(([,v])=>v[0]===depth&&v[1]===density);
    $("limit-preset").value=match?.[0] ?? "custom";
    const exponent=String(depth).replace(/\d/g,n=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(n)]);
    $("limit-computation-summary").textContent=`ε: ${importedStartEpsilon} → 10⁻${exponent} · ${Math.max(1,Math.ceil((depth+Math.log10(Number(importedStartEpsilon)))*density)+1)+extra} samples per branch${extra?` (${extra} extra tail checks)`:""}.`+
      (depth>20 ? " Research depth: partial traces are possible." : "");
    const desired=requestForControls(),exact=window.AtlasNumerics?.getTrace(graphDescriptor(),desired);
    const continuation=exact?.result?.status==="ok"?null:window.AtlasNumerics?.getExtension(graphDescriptor(),desired);
    $("limits-run").textContent=continuation?"Continue to deeper ε":"Compute approach";
    if(continuation)$("limit-computation-summary").textContent+=` ${continuation.result.nodes.length} checked samples retained.`;
    if(s.result?.extended&&s.request&&window.AtlasNumerics?.graphKey(graphDescriptor(s.request))===window.AtlasNumerics.graphKey(graphDescriptor())&&window.AtlasNumerics.traceKey(s.request)===window.AtlasNumerics.traceKey(desired))
      $("limit-computation-summary").textContent=`ε: ${importedStartEpsilon} → 10⁻${exponent} · ${s.nodes.filter(good).length} checked samples (${s.result.reusedCount} retained; ${s.result.newCount} new).`;
  }
  function refreshSamplingLabels() {
    presets={preview:[5,1],compare:[10,1],deeper:[10,3]};
    $("limit-preset").options[0].textContent="Quick preview · 5 decades, 1 sample/decade";
    $("limit-preset").options[1].textContent="Compare limits · 10 decades, 1 sample/decade";
    $("limit-preset").options[2].textContent="Denser check · 10 decades, 3 samples/decade";
    $("limit-research-advice").hidden=!false;
    $("limit-depth-advice").textContent=depthAdvice+" Lamé supports up to 10 decades; deeper choices are available for TV families. A previous TV depth is restored when you return to that family.";
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

  function importedStart(r, t) {
    // Overview rays can begin outside the compact drag region. Preserve their
    // exact interior start rather than silently choosing a different ray.
    if ((query.get("fromOverview")==="1") && [r,t].every(Number.isFinite) && r>0 && t>0 && r<.5 && t<.5 && r+t>.5) return {r,s:t};
    return compactPoint(r,t);
  }

  function compactPoint(r, t) {
    
    const low = .04, high = .48, diagonal = .52;
    if (![r, t].every(Number.isFinite)) return {r:.31, s:.29};
    if (r <= high && t <= high && r + t >= diagonal) return {r, s:t};
    const clamp = x => Math.max(low, Math.min(high, x));
    const x = clamp((r - t + diagonal) / 2);
    return [{r:high,s:clamp(t)}, {r:clamp(r),s:high}, {r:x,s:diagonal-x}]
      .reduce((a,b) => (a.r-r)**2+(a.s-t)**2 <= (b.r-r)**2+(b.s-t)**2 ? a : b);
  }

  function draftCounts() {
    return ["a", "b", "c"].map(name => Number($("limit-" + name).value));
  }

  function pathNote() {
    const vertex = $("limit-path").value.startsWith("vertex");
    $("limit-position").disabled = vertex;
    $("limit-path-note").textContent = vertex
      ? "The diamond is fixed at the selected vertex. Drag the hollow start circle to choose its direction; the edge-position control does not apply to a vertex."
      : "Drag the diamond along a boundary edge, or choose the boundary and its position t above. Drag the hollow circle anywhere inside the compact starting region.";
    if ($("limit-start-note")) $("limit-start-note").textContent =
      (importedLegacyRay?`This imported overview ray begins at ε = ${importedStartEpsilon}, shown by the hollow circle. Its original ε parameter is retained so the computed samples can be shared between pages.`:
      "ε = 1 is the hollow starting circle; ε → 0 approaches the boundary diamond along the displayed segment."+(Number(importedStartEpsilon)!==1?` The imported overview trace samples ε from ${importedStartEpsilon}; the original segment parameter is retained.`:""))+
      " Drag either handle in either diagram. No new computation starts until Compute approach. Dragging chooses a new compact start and resets the sampling start to ε = 1.";
  }

  function controlPath() {
    const name = $("limit-path").value;
    return {
      kind: name.startsWith("edge") ? "edge" : "vertex",
      index: Number(name.at(-1)),
      position: $("limit-position").value,
      ...(importedLegacyRay ? {} : {start: [String(s.start.r), String(s.start.s)]})
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
    return legacySourcePoint(path,0);
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

  function renderSource() {
    
    const api = window.AtlasView;
    if (!api?.sourceAtlasXY || !api?.monodromyMap || !api?.monodromyVertices) return;
    const path = controlPath(), counts = draftCounts();
    window.SeedCases?.updateLameGraph?.(counts,{draft:s.draftDirty});
    if (!counts.every(x => Number.isSafeInteger(x) && x >= 0)) return;
    const node = !s.draftDirty ? s.nodes[s.step] : null;
    const current = nodeSource(node, path);
    const pointText = p => `(${fmt(p[0], 7)}, ${fmt(p[1], 7)})`;
    const draftStart = [s.start.r, s.start.s], draftEnd = boundaryPoint(path);
    const draw = (svg, attachment, tint) => {
      const map = p => {const z=api.monodromyMap(...p,attachment); return api.sourceAtlasXY(z.r,z.s);};
      const vertices = api.monodromyVertices(attachment).map(p => api.sourceAtlasXY(...p));
      const origin=api.sourceAtlasXY(0,0), right=api.sourceAtlasXY(1,0), top=api.sourceAtlasXY(0,1);
      const gridPaths=[[[0,0],[1,0]],[[1,0],[0,1]],[[0,1],[0,0]],[[0,.5],[.5,0]],[[.5,0],[.5,.5]],[[.5,.5],[0,.5]]]
        .map(([a,b]) => `M ${api.sourceAtlasXY(...a).join(" ")} L ${api.sourceAtlasXY(...b).join(" ")}`).join(" ");
      const endpoint=map(draftEnd), start=map(draftStart), other=[0,1,2].filter(i=>i!==path.index);
      const compact=[[.48,.04],[.48,.48],[.04,.48]].map(map);
      const boundary=path.kind==="edge"
        ? `<path d="M ${vertices[other[0]].join(" ")} L ${vertices[other[1]].join(" ")}" stroke="${tint}" stroke-width="4" fill="none"></path>` : "";
      const samples=!s.draftDirty ? s.nodes.map(n=>{
        const q=nodeSource(n,path); if(!q) return ""; const p=map(q);
        return `<circle cx="${p[0]}" cy="${p[1]}" r="2.3" fill="${good(n)?tint:"#929b92"}"></circle>`;
      }).join("") : "";
      const selected=current?map(current):null;
      const handle=(kind,p,shape,label)=>{
        const character=api.monodromyMap(...(kind==="start"?draftStart:draftEnd),attachment);
        return `<g id="${svg.id}-${kind}" data-handle="${kind}" class="source-handle source-handle-${kind}" tabindex="0" role="slider" aria-label="${label}" aria-valuemin="0" aria-valuemax="1" aria-valuenow="${character.r}" aria-valuetext="r = ${fmt(character.r,7)}, s = ${fmt(character.s,7)}" style="cursor:grab;touch-action:none"><circle cx="${p[0]}" cy="${p[1]}" r="16" fill="transparent"></circle>${shape}</g>`;
      };
      svg.innerHTML=`<path d="${gridPaths}" fill="none" stroke="#ccd4ca" stroke-width="1"></path>
        <polygon points="${vertices.map(p=>p.join(",")).join(" ")}" fill="${tint}12" stroke="${tint}" stroke-width="1.5"></polygon>
        <polygon points="${compact.map(p=>p.join(",")).join(" ")}" fill="${tint}10" stroke="${tint}70" stroke-width="1" stroke-dasharray="3 3"></polygon>
        <g fill="#829084" font-size="12" font-family="system-ui">${[[1/3,1/3,"Δ₀"],[1/6,1/6,"Δ₁"],[2/3,1/6,"Δ₂"],[1/6,2/3,"Δ₃"]].map(([r,t,label])=>{const p=api.sourceAtlasXY(r,t);return `<text x="${p[0]}" y="${p[1]}" text-anchor="middle">${label}</text>`;}).join("")}</g>
        ${boundary}<path d="M ${start.join(" ")} L ${endpoint.join(" ")}" stroke="${tint}" fill="none" stroke-width="1.8" stroke-dasharray="6 4"></path>
        ${samples}${selected?`<circle cx="${selected[0]}" cy="${selected[1]}" r="4.6" fill="${tint}" stroke="white" stroke-width="1.3"></circle>`:""}
        ${handle("start",start,`<circle cx="${start[0]}" cy="${start[1]}" r="7" fill="white" stroke="${tint}" stroke-width="2.5"></circle>`,"Starting point; left and right change r, up and down change s; Shift moves faster")}
        ${handle("target",endpoint,`<path d="M ${endpoint[0]} ${endpoint[1]-8} l 8 8 l -8 8 l -8 -8 Z" fill="${tint}" stroke="white" stroke-width="1.8"></path>`,"Boundary target; arrow keys move r and s and snap to the nearest edge; Shift moves faster")}
        <g fill="#596e62" font-size="12" font-family="system-ui"><text x="${origin[0]-13}" y="${origin[1]+19}">0</text>
        <text x="${api.sourceAtlasXY(.5,0)[0]-9}" y="${origin[1]+19}">1/2</text><text x="${right[0]-4}" y="${origin[1]+19}">1</text>
        <text x="${origin[0]-29}" y="${api.sourceAtlasXY(0,.5)[1]+4}">1/2</text><text x="${origin[0]-15}" y="${top[1]+4}">1</text><text x="${right[0]+17}" y="${origin[1]+4}">r</text><text x="${origin[0]-5}" y="${top[1]-12}">s</text></g>`;
    };
    draw($("limit-source"),[0,0,0],baseColor);
    draw($("limit-source-attached"),counts,color);
    const attached=api.monodromyMap(...draftStart,counts), endAttached=api.monodromyMap(...draftEnd,counts);
    const cell=api.monodromyMap(.35,.3,counts).index;
    $("limit-attached-source-title").textContent=`After attachment · n = ${1+counts.reduce((a,b)=>a+b,0)} · Δ${["₀","₁","₂","₃"][cell]}`;
    $("limit-start-readout").textContent=`Start ○: n = 1 ${pointText(draftStart)} · attached ${pointText([attached.r,attached.s])}`;
    $("limit-target-readout").textContent=`Target ◆: n = 1 ${pointText(draftEnd)} · attached ${pointText([endAttached.r,endAttached.s])}`;
    $("limit-source-readout").textContent=current
      ? `Selected sample: n = 1 ${pointText(current)} · ε = ${node.epsilon}`
      : "The two diagrams edit the same start and target. The small filled dot appears during computed playback.";
    $("limit-source-caption").textContent=`${path.kind==="edge"?"Edge":"Vertex"} ${path.index}${path.kind==="edge"?`: t = ${fmt(path.position,6)}`:""}. Hollow circle: start; diamond: limiting target. Both handles are linked between the diagrams.`;
    const info=api.cuspInfo?.(path,counts);
    $("limit-cusp-target").textContent=typeof info==="string"?info:info?.label??"Boundary behaviour depends on the selected path and attachment.";
    $("limit-cusp-detail").textContent=info?.detail??"Computed trajectories do not by themselves prove an exact limiting value.";
    $("limit-draft-status").textContent=s.draftDirty
      ? "Settings changed. Compute approach to update the result plots below."
      : s.request ? "" : "Choose a start and boundary target, then compute the approach.";
    $("limit-draft-status").classList.toggle("settings-dirty",s.draftDirty);
    $("limit-results-state").hidden = !s.draftDirty || !s.request;
    $("limit-results-state").textContent = "Previous computed path — the source settings above have changed. Compute approach to update these plots.";
  }

  function changeDraft(event) {
    if(event?.target && ["limit-path","limit-position"].includes(event.target.id))resetImportedPath();
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
    importedLegacyRay=false;importedStartEpsilon="1";
    query.delete("sourceMode");query.delete("fromOverview");
  }

  function snapBoundary(r,t) {
    
    const vertices=[[.5,0],[.5,.5],[0,.5]];
    let best=null;
    for(let i=0;i<3;i++) {
      const other=[0,1,2].filter(j=>j!==i), a=vertices[other[0]], b=vertices[other[1]];
      const dx=a[0]-b[0],dy=a[1]-b[1];
      const position=Math.max(.02,Math.min(.98,((r-b[0])*dx+(t-b[1])*dy)/(dx*dx+dy*dy)));
      const point=[b[0]+position*dx,b[1]+position*dy];
      const distance=(r-point[0])**2+(t-point[1])**2;
      if(!best||distance<best.distance) best={index:i,position,distance};
    }
    $("limit-path").value="edge"+best.index;
    $("limit-position").value=String(Number(best.position.toPrecision(12)));
  }

  function updateSourceHandle(kind,p) {
    resetImportedPath();
    if(kind==="start") s.start=compactPoint(p.r,p.s);
    else snapBoundary(p.r,p.s);
    changeDraft();
  }

  function connectSourceEditor(svg,attachment) {
    let active=null;
    const counts=()=>attachment?draftCounts():[0,0,0];
    const inverse=p=>window.AtlasView.monodromyInverse(p.r,p.s,counts());
    const eventPoint=event=>{
      const box=svg.getBoundingClientRect(), scale=Math.min(box.width/400,box.height/400);
      const x=(event.clientX-box.left-(box.width-400*scale)/2)/scale;
      const y=(event.clientY-box.top-(box.height-400*scale)/2)/scale;
      return inverse(window.AtlasView.sourceAtlasRS(x,y));
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
      const mapped=window.AtlasView.monodromyMap(...p,counts());
      const step = event.shiftKey ? .01 : .002;
      updateSourceHandle(kind,inverse({r:mapped.r+delta[0]*step,s:mapped.s+delta[1]*step}));
      svg.querySelector(`[data-handle="${kind}"]`)?.focus();
    });
  }

  function infinityReference(base=false) {
    const path = s.request?.path ?? controlPath();
    if (path.kind !== "edge") return null;
    
    const counts = s.request?.counts ?? draftCounts();
    const [a,b,c] = counts;
    // These cyclic permutations are fixed by this site's edge and chart
    // conventions, and were checked against independent computed traces.
    const permuted = s.chart === "raw" && path.index === 0 ? [a,b,c]
      : s.chart === "zero" && path.index === 2 ? [b,c,a]
      : s.chart === "one" && path.index === 1 ? [c,a,b] : null;
    if (!permuted) return null;
    if (base) return {numerator:1,denominator:2,value:.5,label:"1/2",status:"verified"};
    const ref = window.AtlasView?.cuspRational?.(permuted);
    if (!ref || !(ref.denominator > 0) || !Number.isFinite(ref.value)) return null;
    return ref;
  }

  function fractionDecimal(numerator, denominator, digits = 95) {
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
    
    if ($("limit-compare").checked) paths.push({points:base, color:baseColor, selected:s.step});
    paths.push({points:result, color, selected:s.step});
    const ref = s.chart === "raw" ? infinityReference() : null;
    api.drawDomain($("limit-overview"), {
      height:560,
      markedStrip:false,
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
    $("limit-overview-note").textContent = ""+"Solid traces join consecutive resolved samples; unresolved gaps remain open and no segment is extrapolated to an ideal cusp. Dashed image contours bound the compact interior sample. The dark dashed box locates the local view; very small footprints are enlarged for visibility." + (s.chart === "raw" ? "" : " In transformed coordinates, the box approximately projects the lens corners; the actual curved viewport boundary is not drawn.");
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
      const worker = new Worker(file + "?v=20261025-maintenance1");

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
    if($("limit-run-note"))$("limit-run-note").textContent=[result.message,result.tailBranchValidation].filter(Boolean).join(" ");
    updateProgress(result.completed??result.nodes.filter(good).length,result.total??result.nodes.length,result.nodes.findLast(good)?.epsilon);
    
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
      completed + " checked samples. " +
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
    saveCaseResults();
    updateComputationSummary();
  }

  function traceError(message){
    s.result={status:s.nodes.some(good)?"partial":"unresolved",nodes:s.nodes,
      completed:s.nodes.filter(good).length,total:s.progressTotal||s.nodes.length,
      message:String(message??"The numerical worker stopped."),stop:{kind:"error",resumeAvailable:false}};
    renderInterruption();
    saveCaseResults();
  }

  function run(event) {
    event?.preventDefault();

    if (!$("limits-controls").reportValidity()) return;

    const counts = ["a", "b", "c"].map(name =>
      $("limit-" + name).valueAsNumber
    );

    if (
      !counts.every(x =>
        Number.isSafeInteger(x) && x >= 0 && x <= 20
      ) ||
      counts.reduce((a, b) => a + b) > 20
    ) {
      status(
        "trace",
        "Enter nonnegative whole numbers with a + b + c ≤ 20.",
        true
      );
      return;
    }

    const pathValue = $("limit-path").value;
    const path = controlPath();

    const target = "1e-" + $("limit-depth").value;
    const density = Number($("limit-density").value);
    const totalAttachments = counts.reduce((a,b) => a+b,0);
    const tailRefinement = tailRefinementFor(totalAttachments);
    const startEpsilon = importedStartEpsilon;
    const priorView=extensionView({counts,path,target,startEpsilon,samplesPerDecade:density,tailRefinement});

    saveCaseResults();
    cancel();

    const id = ++s.job;

    s.nodes = [];
    s.mesh = null;
    s.result = null;
    s.step = -1;
    s.view = null;
    s.draftDirty = false;
    s.detachedRestore = false;

    s.request = {
      counts,
      path,
      target,
      startEpsilon,
      samplesPerDecade: density,
      tailRefinement
    };

    // Zoom, chart and overlays are viewing choices, retained on deepening.

    $("limits-applied").textContent =
      `Applied attachment: (${counts.join(", ")}), ` +
      `n = ${1 + counts.reduce((a, b) => a + b)}. ` +
      `Path: ${pathValue}, t = ${path.position}; ` +
      `start = (${(path.start ?? [s.start.r,s.start.s]).join(", ")}); ` +
      `ε from ${startEpsilon} to ${target}; ${density} samples per decade.` +
      (tailRefinement ? ` ${tailRefinement} additional tail ${tailRefinement === 1 ? "sample" : "samples"} automatically for this attachment.` : "");

    const url = new URL(location.href);

    ["a", "b", "c"].forEach((name, i) => {
      url.searchParams.set(name, counts[i]);
    });

    url.searchParams.set("path", pathValue);
    url.searchParams.set("t", path.position);
    url.searchParams.set("depth", $("limit-depth").value);
    url.searchParams.set("startEpsilon",String(startEpsilon));
    if(importedLegacyRay){url.searchParams.set("sourceMode","ray");url.searchParams.set("start",String(startEpsilon));url.searchParams.delete("fromR");url.searchParams.delete("fromS");}
    else {url.searchParams.delete("sourceMode");url.searchParams.delete("start");url.searchParams.set("fromR",path.start[0]);url.searchParams.set("fromS",path.start[1]);}
    url.searchParams.set("density", density);

    url.searchParams.delete("autocompute");
    history.replaceState(null, "", url);

    window.AtlasUI?.refreshNavigation(url.searchParams);

    const meshResolution=totalAttachments>11?4:8;
    const meshKey = counts.join(",")+(meshResolution===8?"":"/quick4");
    const traceKey = JSON.stringify([
      counts,
      path,
      target,
      density,
      startEpsilon,
      tailRefinement
    ]);
    const sharedMesh=window.AtlasNumerics?.getMesh(graphDescriptor(s.request),{inset:.04});
    const savedMesh = sharedMesh?.result ?? readResultCache("mesh", meshKey);
    const savedTrace = window.AtlasNumerics?.getTrace(graphDescriptor(s.request),s.request)?.result ?? readResultCache("trace", traceKey);
    s.meshMeta={resolution:sharedMesh?.resolution ?? meshResolution,inset:sharedMesh?.inset ?? .04};
    const extension=savedTrace?.status==="ok"?null:window.AtlasNumerics?.getExtension(graphDescriptor(s.request),s.request);
    s.extension=null;

    status("image", `Computing a ${meshResolution===4?15:45}-point interior image…`);

    status(
      "trace",
      "Following the geometric branch and checking each sample " +
      "at two precisions…"
    );

    prepareExtension(extension,priorView);
    render();

    if (savedMesh?.samples && savedMesh?.triangles) {
      s.mesh = savedMesh;

      status(
        "image",
        `Computed interior image reused (${savedMesh.samples.length} samples); no image computation needed.`
      );
    } else {
      startWorker(
        "mesh",
        "numerics-worker.js",
        {
          type: "compute",
          id,
          counts,
          resolution: meshResolution,
          inset: 0.04
        },
        reply => {
          if (reply.type === "progress") {
            status(
              "image",
              `Image samples: ${reply.completed}/${reply.total}.`
            );
          }

          if (reply.type === "error") {
            settle("mesh");
            status("image", reply.message, true);
          }

          if (reply.type === "done") {
            settle("mesh");

            s.mesh = reply.mesh;s.meshMeta={resolution:meshResolution,inset:.04};
            storeResultCache("mesh", meshKey, s.mesh);

            const count = s.mesh.samples.filter(x =>
              x.result?.status === "ok"
            ).length;

            status(
              "image",
              `Interior image: ${count}/${s.mesh.points.length} ` +
              "resolved samples. Unresolved faces are omitted."
            );

            if (
              s.fit === "image" &&
              Number($("limit-zoom").value) === 0
            ) {
              s.view = null;
            }

            render();
            saveCaseResults();
          }
        }
      );
    }

    if (savedTrace?.status === "ok" && Array.isArray(savedTrace.nodes)) {
      finishTrace(savedTrace);
      status("trace", `Exact source path reused from the shared cache: ${s.nodes.filter(good).length} checked samples.`);
    } else {
      startWorker(
        "trace",
        "limits-worker.js",
        {
          type: "trace",
          id,
          counts,
          path,
          target,
          startEpsilon,
          samplesPerDecade: density,
          tailRefinement,
          ...(extension?{extendFrom:extension}:{})
        },
        reply => {
          if (reply.type === "progress") {
            acceptTraceProgress(reply);

            if (
              (s.fit !== "image" || !s.view) &&
              Number($("limit-zoom").value) === 0
            ) {
              s.view = null;
            }

            render();

            if (Number($("limit-zoom").value) > 0) {
              centerSelected();
            }
          }

          if (reply.type === "error") {
            settle("trace");
            status("trace", reply.message, true);
            traceError(reply.message);
          }

          if (reply.type === "done") {
            settle("trace",!!reply.result?.resumeToken && !!reply.result?.stop?.resumeAvailable);

            if (reply.result.status === "ok") {
              storeResultCache("trace", traceKey, reply.result);

            }

            finishTrace(reply.result);
          }
        }
      );
    }

    render();
  }

  function initializeControls() {
    $("limits-legacy-attachment").hidden = false;
    ["limit-a","limit-b","limit-c"].forEach(id => {$(id).disabled = false;});
    ensureCoordinateOptions();
    refreshSamplingLabels();
    $("limit-chart").value = s.chart;
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
    $("limit-interruption-advice").textContent=s.draftDirty?"Settings have changed. Compute a new approach to use them.":s.detachedRestore?"Saved checked points are restored. The previous worker checkpoint ended when this page or family was left; Compute approach to retry this path.":info.advice;
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
        spanX:offsets.length ? Math.min(span*aspect,Math.max((x1-x0)*1.3, 1e-88)) : span*aspect
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
        {decimal:fractionDecimal(ref.numerator,ref.denominator),tint:color,label:ref.label}
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
          const label=`Re ${variable} → ${reference.label}`;
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
    $("limit-rational-target").hidden=!ref;
    if (ref) {
      const delta=selected ? Math.abs(difference(selected[0],fractionDecimal(ref.numerator,ref.denominator))) : null;
      const theoretical=ref.status === "verified" || ref.verified === true;
      $("limit-rational-target").textContent=
        `${theoretical ? "Horizontal limit" : "Predicted horizontal limit"}: Re ${variable} → ${ref.label}`+
        (delta != null ? ` · distance from this line = ${fmt(delta,5)}` : "")+
        ($("limit-compare").checked && infinityReference(true) ? ` · ${baseLabel}: ${infinityReference(true).label}` : "");
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

  function renderEvidence() {
    const variable=s.chart === "raw" ? "τ" : "w", ref=infinityReference();
    const all=s.nodes.map((node,index)=>({node,index,z:good(node)?values(node):null})).filter(x=>x.z);
    const rows=all.slice(-3), reference=ref?fractionDecimal(ref.numerator,ref.denominator):null;
    $("limit-evidence-coordinate").textContent=`Re ${variable}`;
    $("limit-evidence-distance").textContent=ref ? "Distance to reference" : "Change from prior sample";
    $("limit-evidence-reference").textContent=ref
      ? `${ref.status === "verified" ? "Theoretical limit" : "Calibrated prediction"}: Re ${variable} → ${ref.label}`
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
    $("limit-evidence-note").textContent=note;
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
        "Numerical estimates checked at two working precisions; no interval certificate."
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
    saveDisplaySettings();saveCaseResults();
  });
  window.addEventListener("pagehide", () => {
    saveDisplaySettings();saveCaseResults();
    const wasRunning=s.pending.has("trace");
    cancel();
    if(wasRunning){
      s.result={status:"partial",nodes:s.nodes,completed:s.nodes.filter(good).length,
        total:s.progressTotal||s.nodes.length,message:"The page was left while computing. Checked samples were retained.",stop:{kind:"page_change",resumeAvailable:false}};
      s.detachedRestore=true;
    }
    render();
  });

  if(["raw","zero","one","half"].includes(query.get("chart"))) {s.chart=query.get("chart");$("limit-chart").value=s.chart;}
  initializeControls();
  restoreDisplaySettings();
  pathNote();
  render();
  const savedPath = (query.has("fromR") && query.has("fromS")) ||
    (query.has("start") && query.has("path"));
  // Navigation links are authoritative about source paths. A remembered family
  // view must not replace an explicitly imported start or boundary selection.
  const explicitPath=savedPath || query.has("path") || query.has("t");
  const restoredFamily=!explicitPath && restoreCaseResults();
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
if(window.AtlasResults?.ready?.then)window.AtlasResults.ready.then(initializeLimits);
else initializeLimits();


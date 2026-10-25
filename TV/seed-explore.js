"use strict";

// TV Explore owns point/image/trace scheduling for the active labelled graph.
// main.js owns mounting and disposal when switching between the three families.
(() => {
  if(!window.AtlasExplorer || !document.getElementById("monodromy-source"))return;
  window.AtlasExplorer.register("tv",mountSeedExplorer);
  function mountSeedExplorer({source:nextSource,settings={},initial=true}={}) {
  const api = window.SeedCases;
  const events=new AbortController();
  const listen=(target,type,callback)=>target.addEventListener(type,callback,{signal:events.signal});
  let disposed=false;
  const active = api.active, get = id => document.getElementById(id);
  const params = new URLSearchParams(location.search), inset = .04;
  for(const [key,value] of Object.entries(settings))if(value!=null)params.set(key,value);
  const traceDensity=()=>active.isRepeated && (active.weights.reduce((a,b)=>a+b,0)>12 || (active.seed!=="bouquet"&&!active.caseConfig?.connector)) ? 1 : 2;
  let overviewDensity=traceDensity(), qualityRepeated=null;
  let previewEpoch=0,autoPreview=null,meshCompletion=null;
  const colors = {base: "#285f65", result: "#b55e3c"};
  const graphKey = () => window.AtlasResults?.key({seed:active.seed,id:active.id,caseConfig:active.caseConfig??null})??JSON.stringify({seed:active.seed,id:active.id,caseConfig:active.caseConfig??null});
  const graphDescriptor=()=>({family:active.seed,caseId:active.id,caseConfig:active.caseConfig==null?null:JSON.parse(JSON.stringify(active.caseConfig))});
  // active is changed in place by the selector; retain the identity of data
  // already displayed so a new decoration cannot inherit the old graph's data.
  let appliedGraphKey=graphKey(),pointArtifact=null,meshArtifact=null,traceSettings=null;
  let appliedGraph=graphDescriptor();
  const savedResults=()=>window.AtlasResults?.get("explore-tv",appliedGraphKey);
  // Drag previews marked with 15-digit precision are transient. Use one policy
  // for local snapshots, cross-page sharing, and restoration; keep checked
  // released-point results and explicit unresolved branches as before.
  const isTransientPoint=result=>[result?.base,result?.result].some(branch=>branch?.precisionDigits===15);
  function rememberResults(){
    if(!window.AtlasResults)return;
    const snapshot={...(window.AtlasResults.get("explore-tv",appliedGraphKey)??{})};
    if(pointArtifact&&!isTransientPoint(pointArtifact.result))snapshot.point=pointArtifact;
    if(meshArtifact)snapshot.mesh=meshArtifact;
    if(traceSettings&&traces.some(Boolean))snapshot.traces={...traceSettings,traces:traces.map(trace=>trace?{...trace,status:trace.status==="computing"?"stopped":trace.status}:null)};
    if(snapshot.point||snapshot.mesh||snapshot.traces)window.AtlasResults.set("explore-tv",appliedGraphKey,snapshot);
    publishSharedResults();
  }
  function traceRequest(trace){
    if(trace.sharedRequest)return trace.sharedRequest;
    if(!traceSettings)return null;
    return {path:trace.path??{kind:"edge",index:trace.edge,position:trace.position,start:trace.start??traceSettings.start},target:traceSettings.target,startEpsilon:.04,samplesPerDecade:traceSettings.density,tailRefinement:0};
  }
  function publishSharedResults(){
    const shared=window.AtlasNumerics;if(!shared)return;
    if(pointArtifact&&!isTransientPoint(pointArtifact.result))shared.putPoint(appliedGraph,pointArtifact);
    if(meshArtifact)shared.putMesh(appliedGraph,{resolution:meshArtifact.quality,inset:meshArtifact.inset,result:meshArtifact.result});
    for(const trace of traces.filter(Boolean)){
      const request=traceRequest(trace);
      if(request&&trace.nodes?.length)shared.putTrace(appliedGraph,request,{...trace,path:request.path,status:trace.status==="computing"?"stopped":trace.status});
    }
  }
  function restoreSharedTraces(){
    const available=window.AtlasNumerics?.traces(appliedGraph)??[];
    const selected=Array(6).fill(null),identities=[[],[],[]];
    for(const {request,result} of [...available].reverse()){
      const path=request?.path,edge=Number(path?.index);
      if(path?.kind!=="edge"||!Array.isArray(path.start)||!Number.isInteger(edge)||edge<0||edge>2||!result?.nodes?.length||identities[edge].length>=2)continue;
      const identity=JSON.stringify([Number(path.position),path.start]);if(identities[edge].includes(identity))continue;
      const slot=edge*2+identities[edge].length;identities[edge].push(identity);
      selected[slot]={...result,path,edge,position:path.position,start:path.start,sharedRequest:request};
    }
    const count=selected.filter(Boolean).length;if(!count)return false;
    traces=selected;traceSettings=null;updateEndpoints();scheduleDraw();
    message("cusp-status",`Reused ${count} saved computed edge approach${count===1?"":"es"}, shared by Explore and Limits & Zoom. Original starts, positions and stopping points are retained below. Compute cusp series to replace them with the selected six-series setting.`);
    return true;
  }
  const samePoint=pair=>Array.isArray(pair)&&pair.length===2&&Number(pair[0])===point.r&&Number(pair[1])===point.s;
  const cuspParameters=()=>({position:Number(get("cusp-position").value),target:`1e-${get("cusp-depth").value}`,density:overviewDensity});
  function restoreCuspResults(saved){
    const expected=cuspParameters();
    if(!saved||!Object.keys(expected).every(key=>saved[key]===expected[key])||!Array.isArray(saved.start)||!Array.isArray(saved.traces)||saved.traces.length!==6||!saved.traces.every(trace=>!trace||Array.isArray(trace.nodes)))return false;
    traceSettings={...expected,start:[...saved.start]};traces=saved.traces;
    get("cusp-run").disabled=false;get("cusp-cancel").disabled=true;
    updateEndpoints();scheduleDraw();
    const complete=traces.every(trace=>trace?.status==="ok");
    message("cusp-status",`${complete?"Six saved series":"Saved partial series"} · start (${fmt(saved.start[0])}, ${fmt(saved.start[1])}).${complete?"":" Extend to continue."}`);
    return true;
  }
  function restorePointResult(saved){
    if(!saved?.result||!samePoint(saved.source)||isTransientPoint(saved.result))return false;
    showPoint(saved.result);
    const complete=tau(saved.result.base)&&tau(saved.result.result);
    message("point-status",complete?"Saved period ratios restored.":"Restored a partial point result at this source point. An unresolved branch remains undrawn.",!complete);
    return true;
  }
  let point = nextSource ? api.project(nextSource.r,nextSource.s,inset) : api.project(params.has("r") ? Number(params.get("r")) : active.start[0],
    params.has("s") ? Number(params.get("s")) : active.start[1], inset);
  let pointResult = null, mesh = null, traces = Array(6).fill(null);
  let pointWorker = null, pointBusy = false, pointQueued = null, pointVersion = 0, pointTimer = null;
  let meshWorker = null, meshVersion = 0, traceVersion = 0, traceWorkers = [];
  let pointer = null, redrawFrame = null;
  const fmt = value => Number.isFinite(Number(value)) ? Number(Number(value).toPrecision(9)).toString() : "unresolved";
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
  const tau = branch => branch?.status === "ok" && window.AtlasView.finiteTau(branch.tau) ? branch.tau : null;
  const nodeTau = (node, field) => tau(node?.[field]);
  const makeWorker = () => new Worker("repeated-worker.mjs?v=20261025-tv-maintenance1", {type: "module"});
  const caseRequest = () => ({caseId: active.id, caseConfig: active.caseConfig});

  function message(id, text, error = false) {
    get(id).textContent = text;
    get(id).classList.toggle("error", error);
  }
  function applyLabels(preserveControls=false) {
    const quality=preserveControls?get("quality").value:params.get("quality");
    for (const id of ["a", "b", "c"]) get(id).closest("div").hidden = true;
    get("mesh-run").textContent = "Compute interior image";
    if(qualityRepeated!==active.isRepeated){
      get("quality").innerHTML = active.isRepeated ? '<option value="4">Quick — 15 points</option><option value="8">Standard — 45 points</option><option value="12">Fine — 91 points</option>' : '<option value="20">Precomputed — 231 points</option>';
      qualityRepeated=active.isRepeated;
    }
    get("quality").value = active.isRepeated ? (["4", "8", "12"].includes(quality) ? quality : "8") : "20";
    get("attachment-caption").textContent = active.isRepeated ? "Apply decorations to start the recommended interior image and six cusp approaches. Larger or mixed graphs use lighter sampling to keep the preview responsive. You can choose a finer mesh or deeper cusp check below. Only resolved samples are drawn." : "Apply decorations to load the verified precomputed image and compute its six cusp approaches. You can choose a deeper cusp check below. Dragging uses the selected geometric branch independently.";
    get("base-source-label").textContent = `Base: ${active.baseTitle}`;
    get("attached-source-label").textContent = active.operation === "seed" ? "Same seed · no decoration" : `After decoration: ${active.title}`;
    get("base-heading").textContent = `Base: ${active.baseTitle}`;
    get("result-heading").textContent = active.operation === "seed" ? "Same seed · no decoration" : active.title;
    get("tau-workspace-heading").textContent = "Images in the marked τ chart";
    get("explore-chart-caption").textContent = "This is a marked upper-half-plane chart with Γ(2) reference guides; the points are not moved into a single fundamental domain. Both comparison plots share a horizontal range, which expands to retain branches outside 0 ≤ Re τ ≤ 1. Height is compressed as Im τ / (1 + Im τ). The original τ values appear above. Unresolved samples are omitted; numerical checks are not interval certificates.";
    get("tau-image-caption").textContent = "Shading shows only the computed compact interior. Cusp curves contain numerical samples, stop at the last resolved point, and use no forced endpoint. " + (active.seed==="bouquet"?"The three open source edges approach the raw marked cusps ∞, 0, and 1.":"The raw cusp at ½, when present, is retained in this marking even though it is Γ(2)-equivalent to ∞.");
    get("explore-sampling-note").textContent = "Dark shading samples the compact barycentric interior Kδ with δ = 0.04; its dashed contour is an artificial cutoff. Two numerical series approach each open source edge. Their light ribbon interpolates between computed paths and is not a certified image boundary. Source vertices have direction-dependent finite boundary limits.";
    get("source-domain-note").textContent = "Geometric characters: χ(1) = exp(−2πiR), χ(τ) = exp(2πiS). The translated triangle is an unwrapped lift of ordinary monodromy. Drag either point: the half-period shift is inverted automatically. Every decoration preserves projective monodromy.";
    document.querySelector('label[for="source-r"]').textContent = "Base R";
    document.querySelector('label[for="source-s"]').textContent = "Base S";
    get("source-r").min = "0"; get("source-r").max = String(active.maxR??1);
    get("source-s").min = "0"; get("source-s").max = String(active.top);
    get("source-r").step = "any"; get("source-s").step = "any";
    const cuspNote = get("cusp-sampling-note");
    cuspNote.textContent = `Two approaches per cusp at t and 1 − t show the width of the tail. Apply decorations chooses a recommended mesh and cusp depth for this graph; mixed or larger graphs use lighter overview sampling. Adaptive continuation and arithmetic checks stay enabled. Use the controls for a deeper check, or Limits & Zoom for a focused real-part comparison. Increasing only depth keeps matching checked samples and adds the remaining tail. Research depths are slower and can return partial curves; dragging stays independent.`;
    get("source-region-details").innerHTML = `<summary>Source region and cusp correspondence</summary><p class="caption">The dashed full triangle is excluded. Dragging stays inside its compact barycentric inset, with every barycentric coordinate at least ${inset}. Coordinates use the undeformed seed triangle.</p><div class="edge-cusp-key"><strong>Open source edges → raw marked cusps</strong>${active.edgeLabels.map(label => `<span>${escape(label)}</span>`).join("")}</div><p class="caption">The limit along an open edge is a cusp. Approaching a source vertex can instead give a finite limit depending on direction.</p>`;
    if(!preserveControls){
      get("explore-view").value = params.get("view") === "overlay" ? "overlay" : "separate";
      const p = Number(params.get("cuspT"));
      get("cusp-position").value = p >= .02 && p <= .48 ? String(p) : "0.2";
      get("cusp-depth").innerHTML = '<option value="3">10⁻³ · quick preview</option><option value="5">10⁻⁵ · longer tail</option><option value="8">10⁻⁸</option><option value="10">10⁻¹⁰</option><option value="20">10⁻²⁰ · deeper check</option><optgroup label="Research · slower, may be partial"><option value="40">10⁻⁴⁰</option><option value="80">10⁻⁸⁰</option></optgroup>';
      get("cusp-depth").value = ["3", "5", "8", "10", "20", "40", "80"].includes(params.get("cuspDepth")) ? params.get("cuspDepth") : "8";
    }
  }
  function boundaryPoint(edge, position) {
    if(api.boundaryPoint)return api.boundaryPoint(edge,position);
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
    if (!disposed && redrawFrame === null) redrawFrame = requestAnimationFrame(() => {redrawFrame = null; if(!disposed)redraw();});
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
  function edgeReference(edge,field="result") {
    const info=api.asymptotic(edge,...boundaryPoint(edge,Number(get("cusp-position").value)),field==="base");
    const helper=window.TVCuspPredictions;
    if(!helper||!info)return null;
    const options={weights:field==="base"?[1,0,0,0]:active.weights,chart:info.chart};
    const exact=helper.reference(info,[],options);if(exact)return exact;
    const candidates=[traces[2*edge],traces[2*edge+1]].filter(Boolean).map(trace=>helper.reference(info,(trace.nodes??[]).map(node=>helper.chartValues(node,field,info.chart)),options)).filter(Boolean);
    if(!candidates.length||candidates.some(candidate=>candidate.label!==candidates[0].label))return null;
    return candidates[0];
  }
  const infinityGuides=field=>{
    const ref=edgeReference(0,field);
    return ref?[{re:ref.value,color:colors[field],label:`${ref.status==="proposed"?"Proposed":"Candidate"}: Re τ → ${ref.label}`}]:[];
  };
  function redraw() {
    const overlap = get("explore-view").value === "overlay";
    get("separate-plots").hidden = overlap; get("overlap-plot").hidden = !overlap;
    const layers = field => ({mesh, field, color: colors[field]});
    const paths = field => [...boundaryPaths(field), ...traces.filter(Boolean).map(trace => ({color: colors[field], width: 1.35, radius: 1.5,
      points: trace.nodes.map(node => nodeTau(node, field))}))];
    const bands = field => [0, 1, 2].flatMap(edge => {
      const left = traces[2 * edge], right = traces[2 * edge + 1];
      if (!left || !right || Math.abs(Number(left.position)+Number(right.position)-1)>1e-12 || JSON.stringify(left.path?.start??left.start??null)!==JSON.stringify(right.path?.start??right.start??null)) return [];
      const paired = new Map(right.nodes.map(node => [String(node.epsilon), node]));
      return [{color: colors[field], left: left.nodes.map(node => nodeTau(node, field)), right: left.nodes.map(node => nodeTau(paired.get(String(node.epsilon)), field))}];
    });
    const base = tau(pointResult?.base), result = tau(pointResult?.result);
    const reals=[base,result,...(mesh?.samples??[]).flatMap(p=>[tau(p.base),tau(p.result)]),...traces.filter(Boolean).flatMap(t=>t.nodes.flatMap(n=>[nodeTau(n,"base"),nodeTau(n,"result")]))].filter(Boolean).map(z=>z.re);
    const settings = {height: 520, markedStrip: true,realRange:reals.length?[Math.min(0,...reals),Math.max(1,...reals)]:[0,1]};
    if (overlap) {
      const chart = window.AtlasView.drawDomain(get("overlay-image"), {...settings, references:[...infinityGuides("base"),...infinityGuides("result")], meshLayers: [layers("base"), layers("result")],
        paths: [...paths("base"), ...paths("result")], bands: [...bands("base"), ...bands("result")],
        points: [{tau: base, color: colors.base}, {tau: result, color: colors.result}].filter(p => p.tau),
        arrow: base && result ? {from: base, to: result} : null});
      const distance = chart && base && result ? Math.hypot(...chart.map(base).map((value, index) => value - chart.map(result)[index])) : 0;
      get("overlap-note").textContent = !base || !result ? "The comparison arrow appears after both points resolve." : distance < 12 ? "The two points coincide or are too close to distinguish at this scale." : "The solid arrow compares base τ and decorated τ at the same source character. It is not a continuation path.";
    } else {
      for (const field of ["base", "result"]) window.AtlasView.drawDomain(get(field + "-image"), {...settings, references:infinityGuides(field), meshLayers: [layers(field)],
        paths: paths(field), bands: bands(field), points: [{tau: field === "base" ? base : result, color: colors[field]}].filter(p => p.tau)});
    }
  }
  function showPoint(result) {
    pointResult = result;pointArtifact={source:[point.r,point.s],result};
    // Drag previews are deliberately cheaper; only the refined release result
    // is eligible for reuse by a different page.
    if(![result?.base,result?.result].some(branch=>branch?.precisionDigits===15))window.AtlasNumerics?.putPoint(appliedGraph,pointArtifact);
    for (const field of ["base", "result"]) {
      const z = tau(result?.[field]), output = get(field + "-tau");
      output.setAttribute("aria-busy", "false");
      output.textContent = z ? `τ = ${fmt(z.re)} + ${fmt(z.im)}i` : "Unresolved";
      if (z) output.title = `τ = ${z.re.toPrecision(16)} + ${z.im.toPrecision(16)}i (numerical)`;
      else output.removeAttribute("title");
      get(field + "-diagnostic").textContent = z ? (active.seed==="bouquet" && result[field].precisionDigits===15 ? "Fast preview on the marked graph branch. Release the point for higher-precision refinement." : "Resolved on the selected geometric branch. Numerical result; not an interval certificate.") : "No replacement root is displayed. Move the point to retry.";
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
    if (disposed || pointBusy || !pointQueued) return;
    const job = pointQueued; pointQueued = null;
    if (!pointWorker) {
      pointWorker = makeWorker();
      const worker=pointWorker;
      pointWorker.onmessage = ({data}) => {
        if(disposed||pointWorker!==worker)return;
        if (data.type === "progress") return;
        pointBusy = false;
        if (data.id === pointVersion) {
          if (data.type === "point" && data.result) showPoint(data.result);
          else {clearPoint("Unresolved"); message("point-status", data.message ?? data.error ?? "The selected branch could not be resolved.", true);}
        }
        dispatchPoint();
      };
      pointWorker.onerror = () => {
        if(disposed||pointWorker!==worker)return;
        pointWorker?.terminate(); pointWorker = null; pointBusy = false;
        clearPoint("Unavailable"); message("point-status", "The point solver could not run. Drag to retry.", true); dispatchPoint();
      };
    }
    pointBusy = true; pointWorker.postMessage(job);
  }
  function requestPoint() {
    const shared=window.AtlasNumerics?.getPoint(appliedGraph,[point.r,point.s]);
    if(shared?.result){++pointVersion;pointQueued=null;showPoint(shared.result);message("point-status","Reused the computed marked period ratios at this source point.");return;}
    clearPoint();
    pointQueued = {type: "point", id: ++pointVersion, ...caseRequest(), R: point.r, S: point.s,preview:active.seed==="bouquet" && pointer!==null};
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
    get("mesh-run").disabled = false; get("cancel").disabled = true;
    if (text) message("attachment-status", text);
    const settled=meshCompletion;meshCompletion=null;settled?.resolve();
  }
  function acceptMesh(result, elapsed = 0, cached = false, quality=Number(get("quality").value)) {
    mesh = result;meshArtifact={quality,inset,result};rememberResults();
    const total = result.total??result.points.length, count = field => result.samples.filter(sample => tau(sample?.[field])).length;
    const counts=count("base")===count("result")?`${count("result")}/${total} points`:`base ${count("base")}/${total}, decorated ${count("result")}/${total}`;
    message("attachment-status", `Interior: ${counts} · ${cached ? "saved image" : result.cached ? "precomputed image" : `${elapsed.toFixed(2)} s`}.${!cached&&result.metrics?.reusedSamples?` Kept ${result.metrics.reusedSamples} checked samples.`:""}${result.status!=="ok"?" Partial image. Compute interior image continues the missing samples.":""}`);
    scheduleDraw();
  }
  function restoreMeshResults(exactQuality = false, completeOnly = false) {
    const resolution = Number(get("quality").value), saved=savedResults()?.mesh;
    if(saved?.inset===inset&&Array.isArray(saved.result?.samples)&&(!completeOnly||saved.result.status==="ok")&&(!exactQuality||saved.quality===resolution)){acceptMesh(saved.result,0,true,saved.quality);return true;}
    const sharedMesh=window.AtlasNumerics?.getMesh(appliedGraph,{inset});
    if(sharedMesh?.result&&(!completeOnly||sharedMesh.result.status==="ok")&&(!exactQuality||Number(sharedMesh.resolution)>=resolution)){acceptMesh(sharedMesh.result,0,true,Number(sharedMesh.resolution));return true;}
    return false;
  }
  function computeMesh(event,options={}) {
    event?.preventDefault();
    if(event)++previewEpoch;
    if(!get("mesh-controls").reportValidity())return;
    stopMesh(); const id = meshVersion;
    const resolution = Number(get("quality").value);
    save();
    if(restoreMeshResults(true,Boolean(event)))return Promise.resolve();
    let settle;const promise=new Promise(resolve=>{settle=resolve;});meshCompletion={promise,resolve:settle,quality:resolution};
    get("mesh-run").disabled = true; get("cancel").disabled = false;
    message("attachment-status", `Computing ${active.title} compact image…`);
    const started = performance.now();
    try {
    meshWorker = makeWorker();
    meshWorker.onmessage = ({data}) => {
      if (disposed || data.id !== meshVersion) return;
      if (data.type === "progress") {const p = data.result ?? data; message("attachment-status", `Compact image: ${p.completed}/${p.total} samples…`);}
      else if (data.type === "mesh" && data.result) {const result = data.result; stopMesh(); acceptMesh(result, (performance.now() - started) / 1000,false,resolution);}
      else {stopMesh(); message("attachment-status", data.message ?? data.error ?? "Compact sampling was interrupted. Point dragging is independent.", true);}
    };
    meshWorker.onerror = () => {if (!disposed && id === meshVersion) {stopMesh(); message("attachment-status", "Compact sampling could not run. Point dragging is independent.", true);}};
    const baseMesh=window.AtlasNumerics?.getBaseMesh(appliedGraph,{inset});
    const reuseMesh=meshArtifact?.inset===inset?{...appliedGraph,result:meshArtifact.result}:null;
    meshWorker.postMessage({type: "mesh", id, ...caseRequest(), R: point.r, S: point.s, resolution, inset,...(baseMesh?{baseMesh}:{}),...(reuseMesh?{reuseMesh}:{}),...(options.budgetMs?{budgetMs:options.budgetMs}:{})});
    } catch(error) {stopMesh();message("attachment-status",error.message??"Compact sampling could not start.",true);}
    return promise;
  }
  function updateEndpoints() {
    showAsymptotics();
    get("cusp-endpoints").textContent = traces.map((trace, index) => {
      const edge = Math.floor(index / 2);
      if (!trace) return `${active.edgeLabels[edge]}, series ${index % 2 + 1}: not computed.`;
      const describe = field => {
        const node = trace.nodes.findLast(node => nodeTau(node, field)), z = nodeTau(node, field);
        return z ? `${field === "base" ? "Base" : "Decorated"}: Re τ = ${fmt(z.re)}, Im τ = ${fmt(z.im)} at path ε = ${node.epsilon}` : `${field}: unresolved`;
      };
      return `${active.edgeLabels[edge]}, t = ${trace.position} · ${trace.status}\n${describe("base")}\n${describe("result")}${trace.tailBranchValidation?`\n${trace.tailBranchValidation}.`:""}${trace.message?`\n${trace.message}`:""}`;
    }).join("\n\n");
    window.AtlasUI?.progressTable(get("cusp-progress"),traces,active.edgeLabels,(trace,index)=>{
      const params=new URLSearchParams(location.search);params.set("path","edge"+Math.floor(index/2));
      const start=trace.path?.start??trace.start??[point.r,point.s];
      params.set("t",trace.position);params.set("fromR",start[0]);params.set("fromS",start[1]);
      params.delete("sourceMode");const request=traceRequest(trace);params.set("startEpsilon",request?.startEpsilon??.04);
      if(request){params.set("depth",Math.round(-Math.log10(Number(request.target))));params.set("density",request.samplesPerDecade);params.set("tailRefinement",request.tailRefinement??0);}
      return window.AtlasUI.pageHref("limits.html",params);
    },{scopeKey:appliedGraphKey});
  }
  function showAsymptotics() {
    const element=get("explore-case-asymptotics"),position=Number(get("cusp-position").value);
    const openHelp=new Map([...element.querySelectorAll("details[id]")].map(node=>[node.id,node.open]));
    element.hidden=false;
    const describe=(edge,base)=>{
      const a=api.asymptotic(edge,...boundaryPoint(edge,position),base),variable=a?.chart==="raw"?"τ":"w";
      return a&&Number.isFinite(a.alpha)&&(Number.isFinite(a.C)||Number.isFinite(a.logC))
        ? `Proposed: cusp ${a.cusp}; Re ${variable} → ${a.xLabel??fmt(a.x0)}; Im ${variable} ∼ ${a.logCoefficientLabel??(a.M?`1/(${a.M}π)`:`(${fmt(a.alpha)}/π)`)} log(${Number.isFinite(a.C)?fmt(a.C):`exp(${fmt(a.logC)})`}/δ)`
        : a?.known&&Number.isFinite(a.x0)?`Conjecture: Re ${variable} → ${a.xLabel}; logarithmic slope α = ${fmt(a.alpha)}. The leading constant is not assigned.`:"Mixed decoration: rational candidates are inferred only when the checked tail stabilizes; no full leading formula is asserted.";
    };
    const cards=[0,1,2].map(edge=>{
      const info=api.asymptotic(edge,...boundaryPoint(edge,position)),ref=edgeReference(edge),base=edgeReference(edge,"base"),variable=info?.chart==="raw"?"τ":"w";
      return `<article class="tv-cusp-prediction" data-prediction-status="${ref?.status??"unresolved"}"><span class="tv-prediction-kind">${escape(ref?.kind??"Rational candidate pending")}</span><strong>Edge ${edge+1} · cusp ${escape(info?.cusp??"?")}</strong><span class="tv-rational-number">Re ${variable} → ${ref?window.TVCuspPredictions.fractionHTML(ref):"<span class=\"tv-rational-unresolved\">not resolved yet</span>"}</span><span class="tv-prediction-detail">${escape(info?.chartLabel??"")}${base?` · base: ${escape(base.label)}`:""}</span><details id="cusp-formula-help-${edge}" class="control-help"><summary>Formula &amp; evidence</summary>${window.TVCuspFormulas?.html(active.caseConfig,edge,{compact:true})??""}<p class="caption">${escape(ref?.detail??"Compute or extend the cusp series. A numerical candidate appears only after consecutive resolved samples stabilize near the same simple rational.")}</p></details></article>`;
    }).join("");
    element.innerHTML=`<div class="tv-cusp-predictions" aria-label="Proposed rational cusp limits">${cards}</div><p class="caption">Simulation references, not proved limits.</p><details id="cusp-leading-help" class="control-help"><summary>Leading asymptotics &amp; conventions</summary><p class="caption">The fractions are proposed references for simulation, not proved conclusions. At a finite raw cusp, the rational horizontal value belongs to the displayed cusp-opening coordinate w. Dashed reference lines are independent of the computed curves.</p><p class="caption">At normalized edge position t = ${escape(position)}, away from vertices. The trace parameter ε need not equal edge distance δ.</p><div class="table-scroll"><table><thead><tr><th>Source edge and coordinate</th><th>Base seed</th><th>Selected decoration</th></tr></thead><tbody>${active.edgeLabels.map((label,edge)=>{const a=api.asymptotic(edge,...boundaryPoint(edge,position));return `<tr><th>${escape(label)}${a?`<br><small>${escape(a.chartLabel)}; δ = ${escape(a.deltaLabel)}</small>`:""}</th><td>${escape(describe(edge,true))}</td><td>${escape(describe(edge,false))}</td></tr>`;}).join("")}</tbody></table></div><p class="caption">For the three-loop graphs, the graph proposal is p/(2M), with the marked phase p and resonance M read from the decorated graph. All admissible double-edge decorations use the marked all-count conjecture; the leading constants retain the placement of the decorations. All plotted points remain independently computed roots; the proposed values never replace them.</p><p class="caption"><a href="mathematics.html#count-conjectures">Mathematical background →</a></p></details>`;
    for(const node of element.querySelectorAll("details[id]"))if(openHelp.has(node.id))node.open=openHelp.get(node.id);
  }
  function stopTraces(text) {
    ++traceVersion; traceWorkers.forEach(worker => worker.terminate()); traceWorkers = [];
    traces.forEach(trace => {if (trace?.status === "computing") trace.status = "stopped";});
    get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
    updateEndpoints(); if (text) message("cusp-status", text);
  }
  function computeTraces(event,options={}) {
    event?.preventDefault();
    if(event)++previewEpoch;
    if (!get("cusp-controls").reportValidity()) return;
    const position = Number(get("cusp-position").value);
    if (!(position >= .02 && position <= .48)) return;
    const runStart=[point.r,point.s];
    stopTraces(); const version = traceVersion;
    const target = `1e-${get("cusp-depth").value}`;
    save(); showAsymptotics();
    const retainedResult=trace=>Array.isArray(trace?.nodes)&&(trace.nodes.length>0||/unresolved|failed|unavailable|error/.test(trace.status??""));
    if(!event&&restoreCuspResults(savedResults()?.traces)&&(!options.automatic||traces.every(retainedResult)))return;
    rememberResults();
    const previous=traces.map(trace=>trace?{trace,request:traceRequest(trace)}:null);
    traceSettings={...cuspParameters(),start:runStart.slice()};
    const plans=Array.from({length:6},(_,index)=>{
      const request={path:{kind:"edge",index:Math.floor(index/2),position:index%2?1-position:position,start:runStart.slice()},startEpsilon:.04,target,samplesPerDecade:overviewDensity,tailRefinement:0};
      const keep=options.automatic?previous.find(item=>item?.request?.path?.kind==="edge"&&Number(item.request.path.index)===request.path.index&&Math.abs(Number(item.request.path.position)-Number(request.path.position))<1e-12&&Array.isArray(item.request.path.start)&&item.request.path.start.every((value,j)=>Number(value)===Number(runStart[j]))&&retainedResult(item.trace)):null;
      return {request,keep,extension:keep?null:window.AtlasNumerics?.getExtension(appliedGraph,request)};
    });
    traces=plans.map(({request,keep,extension})=>keep?{...keep.trace,sharedRequest:keep.request}:extension?{nodes:extension.result.nodes,status:"waiting",edge:request.path.index,position:request.path.position,start:runStart.slice(),sharedRequest:request,reusedCount:extension.result.nodes.length}:null);
    updateEndpoints();scheduleDraw();
    get("cusp-run").disabled = true; get("cusp-cancel").disabled = false;
    let nextIndex = 0, completed = 0;
    const poolSize = navigator.hardwareConcurrency > 2 && !(Number.isFinite(navigator.deviceMemory)&&navigator.deviceMemory<=2) ? 2 : 1;
    function finish() {
      if (completed !== 6 || version !== traceVersion) return;
      traceWorkers.forEach(worker => worker.terminate()); traceWorkers = [];
      get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
      const complete = traces.every(trace => trace?.status === "ok");
      rememberResults();
      const reused=traces.reduce((sum,trace)=>sum+(trace?.reusedCount??0),0);
      message("cusp-status",(reused?`Kept ${reused} checked samples. `:"")+(complete?"Six series complete.":`Cusp series partial: ${traces.filter(trace=>trace?.status==="ok").length}/6 complete. Checked samples kept; extend or inspect progress.`));
      updateEndpoints(); scheduleDraw();
    }
    function next(worker) {
      if (version !== traceVersion || nextIndex >= 6) {finish(); return;}
      const index = nextIndex++, edge = Math.floor(index / 2), t = index % 2 ? 1 - position : position;
      const id = version * 8 + index;
      const {request,keep,extension}=plans[index];
      if(keep){traces[index]={...keep.trace,sharedRequest:keep.request};completed++;updateEndpoints();scheduleDraw();next(worker);finish();return;}
      traces[index] = {nodes:extension?.result.nodes??[],edge,position:t,start:runStart.slice(),status:"computing",sharedRequest:request,reusedCount:extension?.result.nodes.length??0};
      const shared=window.AtlasNumerics?.getTrace(appliedGraph,request);
      if(shared?.result?.status==="ok"&&Array.isArray(shared.result.nodes)){
        traces[index]={...shared.result,path:request.path,edge,position:t,start:runStart.slice(),sharedRequest:request};completed++;updateEndpoints();scheduleDraw();next(worker);finish();return;
      }
      worker.onmessage = ({data}) => {
        if (version !== traceVersion || data.id !== id) return;
        if (data.type === "progress") {
          const p = data.result ?? data;
          if(p.solverVersion)Object.assign(traces[index],{solverVersion:p.solverVersion,branchValidation:p.branchValidation,geometryChecks:p.geometryChecks,
            ...(p.bouquetContinuation?{bouquetContinuation:p.bouquetContinuation,tailBranchValidationStatus:p.tailBranchValidationStatus??"continued"}:{})});
          const incoming=Array.isArray(p.nodes)?p.nodes:[...traces[index].nodes,...(p.node?[p.node]:[])];
          traces[index].nodes=extension?window.AtlasNumerics.mergeNodes(extension.result.nodes,incoming):incoming;
          message("cusp-status", `${active.edgeLabels[edge]}, t = ${fmt(t)}: ${p.completed}/${p.total} samples (${completed}/6 series finished).`);
          updateEndpoints(); scheduleDraw();
        } else if (data.type === "trace" && data.result) {
          traces[index] = {...traces[index], ...data.result, status: data.result.status ?? "ok"};
          completed++; updateEndpoints(); scheduleDraw(); next(worker); finish();
        } else {traces[index].status = "unresolved"; traces[index].message=data.message??"The worker did not return a resolved trace."; completed++; next(worker); finish();}
      };
      worker.onerror = () => {
        if (version !== traceVersion || traces[index]?.status!=="computing") return;
        traces[index].status = "unavailable"; completed++; worker.terminate();
        worker.onmessage=null;worker.onerror=null;
        traceWorkers = traceWorkers.filter(w => w !== worker);
        if (nextIndex < 6) addWorker(); finish();
      };
      worker.postMessage({type:"trace",id,...caseRequest(),R:runStart[0],S:runStart[1],...request,...(options.budgetMs?{budgetMs:options.budgetMs}:{}),...(extension?{extendFrom:extension}:{})});
    }
    function addWorker(){
      let worker;
      try{worker=makeWorker();}
      catch(error){
        // A denied replacement must not leave the controls permanently busy.
        // Keep any checked prefixes and exact cached series already available.
        while(nextIndex<6){
          const index=nextIndex++,plan=plans[index],saved=plan.keep?.trace??window.AtlasNumerics?.getTrace(appliedGraph,plan.request)?.result;
          traces[index]=saved?.status==="ok"?{...saved,sharedRequest:plan.request}:{...traces[index],nodes:traces[index]?.nodes??[],status:"unavailable",message:error.message??"The numerical worker could not start."};
          completed++;
        }
        finish();return;
      }
      traceWorkers.push(worker);next(worker);
    }
    for (let i = 0; i < poolSize && nextIndex<6; i++) addWorker();
  }
  function changeNumber(id, direction, large = false) {
    const input = get(id), current = input.valueAsNumber;
    if (!Number.isFinite(current)) return;
    const delta = (id === "cusp-position" ? .01 : .002) * (large ? 5 : 1);
    const value = Math.max(Number(input.min), Math.min(Number(input.max), Number((current + direction * delta).toPrecision(15))));
    if (id === "cusp-position") {++previewEpoch;input.value = String(value); stopTraces("Paired positions changed. Recompute to apply them; the existing numerical curves remain visible."); showAsymptotics();}
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

  function runRecommendedPreview() {
    const key=appliedGraphKey;
    if(autoPreview?.key===key&&autoPreview.token===previewEpoch){if(get("decoration-status"))message("decoration-status","Decorations are applied. The recommended preview is already running; checked results stay visible.");return;}
    const profile=window.AtlasPreview?.recommend({family:active.seed,caseConfig:active.caseConfig,isRepeated:active.isRepeated,operation:active.operation})??{label:"Standard preview",quality:active.isRepeated?8:20,depth:3,density:1,reason:"Adaptive arithmetic checks remain enabled."};
    const token=++previewEpoch;
    const wanted=Math.max(Number(profile.quality),Number(meshArtifact?.quality)||0,Number(meshCompletion?.quality)||0);
    const qualities=[...get("quality").options].map(option=>Number(option.value));
    get("quality").value=String(qualities.includes(wanted)?wanted:qualities.find(value=>value>=wanted)??qualities.at(-1));
    if(!traceWorkers.length)overviewDensity=profile.density;
    let depth=Math.max(profile.depth,traceWorkers.length?Number(get("cusp-depth").value):0);
    if([...get("cusp-depth").options].some(option=>Number(option.value)===depth))get("cusp-depth").value=String(depth);
    const note=get("preview-recommendation");
    if(note)note.textContent=`${profile.label}: up to ${(Number(get("quality").value)+1)*(Number(get("quality").value)+2)/2} interior samples, cusp ε down to 10⁻${String(get("cusp-depth").value).replace(/\d/g,d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d])}, ${overviewDensity} sample${overviewDensity===1?"":"s"} per decade. ${profile.reason} Finer saved images and compatible deeper traces are retained. These settings change sampling; arithmetic precision stays adaptive.`;
    if(get("decoration-status"))message("decoration-status",`${active.title} applied. The interior image runs first, then all six cusp approaches.`);
    save();showAsymptotics();
    const requested=JSON.stringify(cuspParameters()),job={key,token};autoPreview=job;
    get("cusp-cancel").disabled=false;
    if(!traceWorkers.length)message("cusp-status","Six cusp approaches are queued after the interior image. Stop series cancels this queued part.");
    const meshReady=meshWorker&&meshCompletion?meshCompletion.promise:computeMesh(undefined,{budgetMs:profile.meshBudgetMs??profile.budgetMs});
    Promise.resolve(meshReady).then(()=>{
      if(disposed||token!==previewEpoch||key!==appliedGraphKey||requested!==JSON.stringify(cuspParameters()))return;
      if(!traceWorkers.length)computeTraces(undefined,{automatic:true,budgetMs:profile.budgetMs});
    }).catch(error=>{if(!disposed&&token===previewEpoch)message("cusp-status",error.message??"The cusp preview could not start. Computed samples are retained.",true);}).finally(()=>{if(autoPreview===job){autoPreview=null;if(!traceWorkers.length)get("cusp-cancel").disabled=true;}});
  }
  function applyCaseChange() {
    // Keep the page and controls in place. No result from the previous labelled
    // graph may replace the freshly applied case, including an idle point job.
    if(graphKey()===appliedGraphKey){runRecommendedPreview();return;}
    ++previewEpoch;autoPreview=null;
    rememberResults();
    clearTimeout(pointTimer);pointTimer=null;
    pointWorker?.terminate();pointWorker=null;pointBusy=false;pointQueued=null;++pointVersion;
    stopMesh();stopTraces();pointer=null;
    pointResult=null;mesh=null;traces=Array(6).fill(null);
    pointArtifact=null;meshArtifact=null;traceSettings=null;appliedGraphKey=graphKey();appliedGraph=graphDescriptor();
    overviewDensity=traceDensity();
    applyLabels(true);showSource();showAsymptotics();updateEndpoints();
    const saved=savedResults();
    if(!restorePointResult(saved?.point))requestPoint();
    if(!restoreMeshResults())message("attachment-status","Ready to compute the interior.");
    if(!restoreCuspResults(saved?.traces)&&!restoreSharedTraces())message("cusp-status","Decoration changed. Compute cusp series to draw its six approaches, or open Limits & Zoom for one edge.");
    scheduleDraw();
    runRecommendedPreview();
  }

  applyLabels(); showSource(); showAsymptotics();
  // These persistent DOM nodes belonged to the previous family. Clear its
  // series table and endpoints before a compatible saved bundle is restored.
  stopTraces();
  listen(get("mesh-controls"), "submit", computeMesh);
  listen(get("cancel"), "click", () => stopMesh("Compact image computation stopped. Point dragging stays active."));
  listen(get("point-form"), "submit", event => {event.preventDefault(); if (get("point-form").reportValidity()) {move(get("source-r").valueAsNumber, get("source-s").valueAsNumber, true); dispatchPoint();}});
  listen(get("explore-view"), "change", () => {save(); scheduleDraw();});
  listen(get("cusp-controls"), "submit", computeTraces);
  listen(get("cusp-cancel"), "click", () => {++previewEpoch;stopTraces("Stopped. Computed samples remain visible.");});
  listen(get("cusp-position"), "input", () => {++previewEpoch;stopTraces("Paired positions changed. Recompute to apply them."); showAsymptotics();});
  listen(get("quality"), "change",()=>{++previewEpoch;});
  listen(get("cusp-depth"), "change",()=>{++previewEpoch;});
  for (const id of ["source-r", "source-s", "cusp-position"]) {
    for (const [suffix, direction] of [["up", 1], ["down", -1]]) listen(get(`${id}-${suffix}`), "click", event => changeNumber(id, direction, event.shiftKey));
    listen(get(id), "keydown", event => {if (["ArrowUp", "ArrowDown"].includes(event.key)) {event.preventDefault(); changeNumber(id, event.key === "ArrowUp" ? 1 : -1, event.shiftKey);}});
  }
  for (const [id, attached] of [["monodromy-source", false], ["monodromy-attached", true]]) {
    const svg = get(id);
    listen(svg, "pointerdown", event => {if (event.button != null && event.button !== 0) return; pointer = {svg, id: event.pointerId}; svg.setPointerCapture(event.pointerId); pointerPoint(event, svg, attached);});
    listen(svg, "pointermove", event => {if (pointer?.svg === svg && pointer.id === event.pointerId) pointerPoint(event, svg, attached);});
    const finish = event => {if (pointer?.svg === svg && pointer.id === event.pointerId) {pointer = null;if(active.seed==="bouquet")requestPoint();dispatchPoint();save();}};
    listen(svg, "pointerup", finish); listen(svg, "pointercancel", finish);
    listen(svg, "lostpointercapture", () => {if (pointer?.svg === svg) pointer = null;});
    listen(svg, "keydown", event => {
      const step = event.shiftKey ? .01 : .002, d = {ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, step], ArrowDown: [0, -step]}[event.key];
      if (d) {event.preventDefault(); move(point.r + d[0], point.s + d[1], true);}
    });
  }
  listen(window, "resize", scheduleDraw);
  listen(window, "atlas:casechange",applyCaseChange);
  function dispose() {
    if(disposed)return;rememberResults();disposed=true;events.abort();
    clearTimeout(pointTimer);pointTimer=null;
    ++pointVersion;++meshVersion;++traceVersion;++previewEpoch;autoPreview=null;
    meshCompletion?.resolve();meshCompletion=null;
    pointWorker?.terminate();pointWorker=null;meshWorker?.terminate();meshWorker=null;
    traceWorkers.forEach(worker=>worker.terminate());traceWorkers=[];
    pointQueued=null;pointBusy=false;
    if(pointer?.svg.hasPointerCapture?.(pointer.id))pointer.svg.releasePointerCapture(pointer.id);
    pointer=null;
    if(redrawFrame!==null)cancelAnimationFrame(redrawFrame);redrawFrame=null;
  }
  const saved=savedResults();
  if(!restorePointResult(saved?.point))requestPoint();
  if(!restoreMeshResults())message("attachment-status","Ready to compute the interior.");
  const tailsRestored=restoreCuspResults(saved?.traces)||restoreSharedTraces();
  if(!tailsRestored)message("cusp-status","Six approach series · ready to compute.");
  return {dispose,remember:rememberResults,source:()=>({...point})};
  }
})();

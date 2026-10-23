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
// Each source family owns its controller only while it is displayed. Keeping
// the common DOM nodes lets focus, details, graph size and display settings stay.
if (document.getElementById("monodromy-source")) {
  const ids = ["monodromy-source", "monodromy-attached", "base-source-label", "attached-source-label", "base-heading", "result-heading", "tau-workspace-heading", "attachment-caption", "explore-chart-caption", "tau-image-caption", "explore-sampling-note", "source-domain-note", "source-region-details", "quality", "cusp-depth"];
  const templates = new Map(ids.map(id => [id, document.getElementById(id)?.innerHTML]));
  const sourceLabels = ["source-r", "source-s"].map(id => ({id, label: document.querySelector(`label[for="${id}"]`).textContent, min: document.getElementById(id).min, max: document.getElementById(id).max, step: document.getElementById(id).step}));
  const cuspDescription = document.getElementById("cusp-sampling-note").textContent;
  const settingIds = {view:"explore-view",quality:"quality",cuspT:"cusp-position",cuspDepth:"cusp-depth"};
  const preferred = Object.fromEntries(Object.entries(settingIds).map(([key,id]) => [key, query.get(key) ?? document.getElementById(id).value]));
  if(window.SeedCases?.isNew && !query.has("quality"))preferred.quality="8";
  let mounted = null, applied = null, hiddenSource = null, cacheReady = !window.AtlasResults?.ready?.then;
  const mounts = new Map();
  function savePageDisplay() {
    remember();
    window.AtlasResults?.set("explore-display-v1","current",{settings:{...preferred},details:[...document.querySelectorAll("details")].map(node=>node.open)});
  }
  function restorePageDisplay() {
    const saved=window.AtlasResults?.get("explore-display-v1","current");
    if(saved?.settings)Object.assign(preferred,saved.settings);
    return saved;
  }
  function readControls() {return Object.fromEntries(Object.entries(settingIds).map(([key,id]) => [key,document.getElementById(id).value]));}
  function remember() {
    const actual = readControls();
    for (const key of Object.keys(actual)) if (!applied || actual[key] !== applied[key]) preferred[key] = actual[key];
  }
  function rememberDetails() {return [...document.querySelectorAll("details")].map(node => ({node,open:node.open,id:node.id,parent:node.parentElement?.id}));}
  function restoreDetails(states) {for (const state of states) {const node=state.node.isConnected?state.node:state.id?document.getElementById(state.id):state.parent?document.getElementById(state.parent)?.querySelector("details"):null;if(node)node.open=state.open;}}
  function restoreClassical() {
    for (const [id,html] of templates) if(html!=null)document.getElementById(id).innerHTML=html;
    for (const field of sourceLabels) {const input=document.getElementById(field.id);input.min=field.min;input.max=field.max;input.step=field.step;document.querySelector(`label[for="${field.id}"]`).textContent=field.label;}
    for (const id of countNames) document.getElementById(id).closest("div").hidden=false;
    document.getElementById("compute").textContent="Apply & compute";
    document.getElementById("mesh-run").textContent="Compute interior image";
    document.getElementById("cusp-sampling-note").textContent=cuspDescription;
    document.getElementById("explore-case-asymptotics").hidden=true;
  }
  function mount(detail={},initial=false) {
    if(!cacheReady)return;
    if(validCounts(detail.counts??[]))selectedCounts=[...detail.counts];
    const kind=window.SeedCases?.isNew?"tv":"lame", factory=mounts.get(kind);
    if(!factory)return;
    const states=rememberDetails();
    if(kind==="lame")restoreClassical();
    const recommendation=document.getElementById("preview-recommendation");
    if(recommendation)recommendation.textContent="Apply & compute chooses a recommended preview for this graph. Saved results are restored automatically; the controls below refine either image.";
    mounted=factory({source:detail.source,settings:{...preferred},initial});
    applied=readControls();restoreDetails(states);
  }
  window.AtlasExplorer={
    register(kind,factory) {mounts.set(kind,factory);if(!mounted&&kind===(window.SeedCases?.isNew?"tv":"lame"))mount({},true);},
    get source() {return mounted?.source?.()??null;}
  };
  window.addEventListener("atlas:beforefamilychange",event=>{
    remember();event.detail.source=mounted?.source?.()??event.detail.source;
    mounted?.dispose();mounted=null;
  });
  window.addEventListener("atlas:familychange",event=>{
    // Support direct family events as well as the normal before/change pair.
    if(mounted){remember();mounted.dispose();mounted=null;}
    mount(event.detail??{});
  });
  window.addEventListener("pagehide",()=>{savePageDisplay();hiddenSource=mounted?.source?.()??hiddenSource;mounted?.dispose();mounted=null;});
  window.addEventListener("atlas:cachedpage",()=>{
    // Back/Forward may revive this document without running its scripts again.
    // The old workers stay terminated; a fresh controller reads saved data.
    if(!mounted)mount({source:hiddenSource,counts:[...selectedCounts]},false);
  });
  window.addEventListener("atlas:beforepagechange",()=>{savePageDisplay();mounted?.remember?.();});
  const initialDisplay=cacheReady?restorePageDisplay():null;
  window.AtlasExplorer.register("lame",connectExplorer);
  const applySavedDetails=saved=>{if(Array.isArray(saved?.details))[...document.querySelectorAll("details")].forEach((node,i)=>{if(typeof saved.details[i]==="boolean")node.open=saved.details[i];});};
  if(cacheReady)applySavedDetails(initialDisplay);
  else Promise.resolve(window.AtlasResults.ready).catch(()=>null).then(()=>{const saved=restorePageDisplay();cacheReady=true;if(!mounted)mount({},true);applySavedDetails(saved);});
}

function connectExplorer({source:nextSource,settings={},initial=true}={}) {
  const events = new AbortController();
  const listen = (target,type,callback) => target.addEventListener(type,callback,{signal:events.signal});
  let disposed=false;
  const query=new URLSearchParams(location.search);
  for(const [key,value] of Object.entries(settings))if(value!=null)query.set(key,value);
  const get = id => document.getElementById(id);
  const source = get("monodromy-source");
  const attachedSource = get("monodromy-attached");
  get("explore-view").value = query.get("view") === "overlay" ? "overlay" : "separate";
  const initialPair = Number(query.get("cuspT"));
  get("cusp-position").value = initialPair >= .02 && initialPair <= .48 ? String(initialPair) : "0.2";
  get("cusp-depth").value = ["3","5","8","10"].includes(query.get("cuspDepth")) ? query.get("cuspDepth") : "8";
  // Retain numerical data only; the shared bounded cache never keeps workers.
  let appliedCounts = [...selectedCounts], pointArtifact = null, meshArtifact = null, traceSettings = null;
  let overviewDensity=2, overviewRefinement=null, previewEpoch=0, autoPreview=null, meshCompletion=null;
  const graph = () => ({family:"lame",counts:[...appliedCounts]});
  const resultKey = counts => counts.join(",");
  const savedResults = () => window.AtlasResults?.get("explore-lame",resultKey(appliedCounts));
  function rememberResults() {
    if(!window.AtlasResults)return;
    const key=resultKey(appliedCounts), previous=window.AtlasResults.get("explore-lame",key)??{};
    const snapshot={...previous};
    if(pointArtifact)snapshot.point=pointArtifact;
    if(meshArtifact)snapshot.mesh=meshArtifact;
    if(traceSettings&&cuspTraces.some(Boolean))snapshot.traces={...traceSettings,traces:cuspTraces.map(trace=>trace?{...trace,status:trace.status==="computing"?"stopped":trace.status}:null)};
    if(snapshot.point||snapshot.mesh||snapshot.traces)window.AtlasResults.set("explore-lame",key,snapshot);
    publishSharedResults();
  }
  function traceRequest(trace) {
    if(trace.sharedRequest)return trace.sharedRequest;
    if(!traceSettings)return null;
    return {path:trace.path??{kind:"edge",index:trace.edge,position:trace.position},target:traceSettings.target,startEpsilon:".04",samplesPerDecade:traceSettings.density,tailRefinement:traceSettings.tailRefinement};
  }
  function publishSharedResults() {
    const shared=window.AtlasNumerics;if(!shared)return;
    if(pointArtifact)shared.putPoint(graph(),pointArtifact);
    if(meshArtifact)shared.putMesh(graph(),{resolution:meshArtifact.quality,inset:meshArtifact.inset,result:meshArtifact.result});
    for(const trace of cuspTraces.filter(Boolean)){
      const request=traceRequest(trace);
      if(request&&trace.nodes?.length)shared.putTrace(graph(),request,{...trace,path:request.path,status:trace.status==="computing"?"stopped":trace.status});
    }
  }
  function restoreSharedTraces() {
    const available=window.AtlasNumerics?.traces(graph())??[];
    const selected=Array(6).fill(null),positions=[[],[],[]];
    for(const {request,result} of [...available].reverse()){
      const path=request?.path,edge=Number(path?.index);
      if(path?.kind!=="edge"||!Number.isInteger(edge)||edge<0||edge>2||!result?.nodes?.length||positions[edge].length>=2)continue;
      const identity=JSON.stringify([Number(path.position),path.start??null]);
      if(positions[edge].includes(identity))continue;
      const slot=edge*2+positions[edge].length;positions[edge].push(identity);
      selected[slot]={...result,path,edge,position:path.position,sharedRequest:request};
    }
    const count=selected.filter(Boolean).length;if(!count)return false;
    cuspTraces=selected;traceSettings=null;
    updateCuspEndpoints();redraw();
    message("cusp-status",`Reused ${count} saved computed edge approach${count===1?"":"es"}, shared by Explore and Limits & Zoom. Original paths, positions and stopping points are retained below. Compute cusp series to replace them with the selected six-series setting.`);
    return true;
  }
  const samePoint = pair => Array.isArray(pair)&&pair.length===2&&Number(pair[0])===point.r&&Number(pair[1])===point.s;
  const cuspParameters = () => {
    const sum=appliedCounts.reduce((a,b)=>a+b,0);
    return {position:String(Number(get("cusp-position").value)),target:`1e-${get("cusp-depth").value}`,density:overviewDensity,tailRefinement:overviewRefinement??(sum>=8?2:sum>=4?1:0)};
  };
  function restoreCuspResults(saved) {
    const expected=cuspParameters();
    if(!saved||!Object.keys(expected).every(key=>saved[key]===expected[key])||!Array.isArray(saved.traces)||saved.traces.length!==6||!saved.traces.every(trace=>!trace||Array.isArray(trace.nodes)))return false;
    traceSettings=expected;cuspTraces=saved.traces;
    get("cusp-run").disabled=false;get("cusp-cancel").disabled=true;
    updateCuspEndpoints();redraw();
    const complete=cuspTraces.every(trace=>trace?.status==="ok");
    message("cusp-status",complete?"Restored six saved computed cusp approaches; no recomputation needed.":"Restored partial saved cusp approaches. Only checked samples are drawn; Compute or extend series reuses compatible checked paths.");
    return true;
  }
  let cuspWorkers = [], cuspVersion = 0, cuspTraces = Array(6).fill(null);
  let meshWorker = null, meshVersion = 0, currentMesh = null;
  let pointWorker = null, pointVersion = 0, pointBusy = false, queuedPoint = null;
  let activePoint = null, pointTimer = null, pointResult = null, pointerId = null, pointerSource = null;
  let point = nextSource ? projectCompact(nextSource.r,nextSource.s) : projectCompact(query.has("r") ? Number(query.get("r")) : 0.31,
    query.has("s") ? Number(query.get("s")) : 0.29);
  countNames.forEach((name, i) => { get(name).value = selectedCounts[i]; });
  if(["4","8","12","18"].includes(query.get("quality")))get("quality").value=query.get("quality");
  listen(get("quality"), "change",()=>{++previewEpoch;});
  listen(get("cusp-depth"), "change",()=>{++previewEpoch;});
  window.SeedCases?.updateLameGraph?.(selectedCounts);
  countNames.forEach(name=>listen(get(name), "input",()=>{
    const draft=countNames.map(key=>get(key).valueAsNumber);
    if(validCounts(draft))window.SeedCases?.updateLameGraph?.(draft,{draft:draft.some((x,i)=>x!==selectedCounts[i])});
  }));

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
      if (!left || !right || Math.abs(Number(left.position)+Number(right.position)-1)>1e-12 || JSON.stringify(left.path?.start??null)!==JSON.stringify(right.path?.start??null)) return [];
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
    window.AtlasUI?.progressTable(get("cusp-progress"),cuspTraces,
      window.AtlasView.cuspForEdge.map(cusp=>`Cusp ${cusp}`),(trace,index)=>{
        const params=new URLSearchParams(location.search);
        params.set("path","edge"+Math.floor(index/2));params.set("t",trace.position);
        const request=traceRequest(trace),start=request?.path?.start??trace.path?.start;
        if(start){params.set("fromR",start[0]);params.set("fromS",start[1]);params.delete("sourceMode");}
        else{params.delete("fromR");params.delete("fromS");params.delete("r");params.delete("s");params.set("start",".04");params.set("sourceMode","ray");}
        params.set("fromOverview","1");params.set("startEpsilon",request?.startEpsilon??".04");
        if(request){params.set("depth",Math.round(-Math.log10(Number(request.target))));params.set("density",request.samplesPerDecade);params.set("tailRefinement",request.tailRefinement??0);}
        return window.AtlasUI.pageHref("limits.html",params);
      });
  }
  function stopCuspApproaches(text = "Stopped. Only computed samples are retained.") {
    ++cuspVersion; cuspWorkers.forEach(worker => worker.terminate()); cuspWorkers = [];
    cuspTraces.forEach(trace => {if (trace?.status === "computing") trace.status = "stopped";});
    updateCuspEndpoints();
    get("cusp-run").disabled = false; get("cusp-cancel").disabled = true;
    message("cusp-status", text);
  }
  function startCuspApproaches(event, options={}) {
    event?.preventDefault();
    if(event)++previewEpoch;
    if (!get("cusp-controls").reportValidity()) return;
    const position = String(Number(get("cusp-position").value));
    if (!(Number(position) >= .02 && Number(position) <= .48)) return;
    const positions = [position, String(1 - Number(position))];
    const target = `1e-${get("cusp-depth").value}`;
    stopCuspApproaches(); const id = cuspVersion;
    const counts = [...selectedCounts], sum = counts.reduce((a, b) => a + b, 0);
    const tailRefinement = cuspParameters().tailRefinement;
    const retainedResult=trace=>Array.isArray(trace?.nodes)&&(trace.nodes.length>0||/unresolved|failed|unavailable|error/.test(trace.status??""));
    if(!event&&restoreCuspResults(savedResults()?.traces)&&(!options.automatic||cuspTraces.every(retainedResult)))return;
    rememberResults();
    const previous=cuspTraces.map(trace=>trace?{trace,request:traceRequest(trace)}:null);
    traceSettings=cuspParameters();
    const plans=Array.from({length:6},(_,index)=>{
      const request={path:{kind:"edge",index:Math.floor(index/2),position:positions[index%2]},startEpsilon:".04",target,samplesPerDecade:overviewDensity,tailRefinement};
      const keep=options.automatic?previous.find(item=>item?.request?.path?.kind==="edge"&&Number(item.request.path.index)===request.path.index&&Math.abs(Number(item.request.path.position)-Number(request.path.position))<1e-12&&!item.request.path.start&&retainedResult(item.trace)):null;
      return {request,keep,extension:keep?null:window.AtlasNumerics?.getExtension(graph(),request)};
    });
    cuspTraces=plans.map(({request,keep,extension})=>keep?{...keep.trace,sharedRequest:keep.request}:extension?{nodes:extension.result.nodes,status:"waiting",edge:request.path.index,position:request.path.position,sharedRequest:request,reusedCount:extension.result.nodes.length}:null);
    updateCuspEndpoints();redraw();
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
      rememberResults();
      const reused=cuspTraces.reduce((sum,trace)=>sum+(trace?.reusedCount??0),0);
      message("cusp-status", complete ? `${reused?`Kept ${reused} previous checked samples. `:""}Six computed approaches complete for (${counts.join(", ")}), t = ${positions.join(" and ")}. ${tailRefinement ? `${tailRefinement} extra tail samples per series for this attachment. ` : ""}No extrapolated endpoints.` : "Paired cusp computations include partial traces. Checked data and original stopping points are retained; use Compute or extend series for an explicit retry or deeper check. Ribbons stop at their last shared resolved sample.");
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
      const {request,keep,extension}=plans[index];
      if(keep){cuspTraces[index]={...keep.trace,sharedRequest:keep.request};finishJob(slot);return;}
      cuspTraces[index] = {nodes:extension?.result.nodes??[],status:"computing",edge,position:pathPosition,sharedRequest:request,reusedCount:extension?.result.nodes.length??0};
      const shared=window.AtlasNumerics?.getTrace(graph(),request);
      if(shared?.result?.status==="ok"&&Array.isArray(shared.result.nodes)){
        cuspTraces[index]={...shared.result,path:request.path,edge,position:pathPosition,sharedRequest:request};finishJob(slot);return;
      }
      message("cusp-status", `Computing cusp ${window.AtlasView.cuspForEdge[edge]}, series ${series + 1}/2 (${index + 1}/6), t = ${pathPosition}…`);
      try {
        slot.job = {id:jobId,index,edge,position:pathPosition,extension};
        if (!slot.worker) {
          const worker = new Worker("limits-worker.js?v=20261022-auto2");
          slot.worker = worker; cuspWorkers.push(worker);
          worker.onmessage = ({data}) => {
            const job = slot.job;
            if (id !== cuspVersion || slot.worker !== worker || !cuspWorkers.includes(worker) || !job || data.id !== job.id) return;
            const index = job.index, edge = job.edge, series = index % 2;
            if (data.type === "progress") {
              const incoming=Array.isArray(data.nodes)?data.nodes:[...cuspTraces[index].nodes,data.node];
              cuspTraces[index].nodes=job.extension?window.AtlasNumerics.mergeNodes(job.extension.result.nodes,incoming):incoming;
              message("cusp-status", `Cusp ${window.AtlasView.cuspForEdge[edge]}, series ${series + 1}/2 (${index + 1}/6): ${data.completed}/${data.total}; ε ≈ ${Number(data.node.epsilon).toExponential(2)}.`);
              updateCuspEndpoints(); redraw();
            }
            if (data.type === "done") {cuspTraces[index] = {...data.result,edge,position:job.position,sharedRequest:plans[index].request};finishJob(slot);}
            if (data.type === "error") {cuspTraces[index].status = "unresolved"; cuspTraces[index].message=data.message; finishJob(slot);}
          };
          worker.onerror = () => {
            if (id !== cuspVersion || slot.worker !== worker || !slot.job) return;
            cuspTraces[slot.job.index].status = "unavailable";
            closeSlot(slot); finishJob(slot);
          };
        }
        slot.worker.postMessage({type:"trace",id:jobId,counts,...request,...(options.budgetMs?{budgetMs:options.budgetMs}:{}),...(extension?{extendFrom:extension}:{})});
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
    url.searchParams.set("quality",get("quality").value);
    url.searchParams.set("cuspT",get("cusp-position").value);
    url.searchParams.set("cuspDepth",get("cusp-depth").value);
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
    pointArtifact={source:[Number(result.r??point.r),Number(result.s??point.s)],result};
    window.AtlasNumerics?.putPoint(graph(),pointArtifact);
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
    if (disposed || pointBusy || !queuedPoint) return;
    const job = queuedPoint; queuedPoint = null;
    try {
      if (!pointWorker) {
        pointWorker = new Worker("numerics-worker.js?v=20261021-overall-final1");
        pointWorker.onmessage = ({data}) => {
          if(disposed)return;
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
          if(disposed)return;
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
    const shared=window.AtlasNumerics?.getPoint(graph(),[point.r,point.s]);
    if(shared?.result){queuedPoint=null;showPointResult(shared.result);message("point-status","Reused the computed period ratios at this source point.");return;}
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
    get("mesh-run").disabled = false; get("cancel").disabled = true;
    const settled=meshCompletion;meshCompletion=null;settled?.resolve();
  }
  function acceptMesh(mesh, elapsed, cached = false, quality=Number(get("quality").value)) {
    currentMesh = mesh;meshArtifact={quality,inset:.04,result:mesh};rememberResults();
    const count = name => mesh.samples.filter(sample => sample?.[name]?.status === "ok").length;
    message("attachment-status", `Compact mesh: base ${count("base")}/${mesh.points.length}; result ${count("result")}/${mesh.points.length} samples. ` +
      (cached ? "Reused a saved computed image; the sample counts show its actual resolution." : `${elapsed.toFixed(2)} seconds.`) + (mesh.status === "ok" ? "" : " Partial image: unresolved faces are omitted. Compute interior image retries with the normal time budget."));
    redraw();
  }
  function restoreMeshResults(exactQuality = false, completeOnly = false) {
    const resolution = Number(get("quality").value), saved = savedResults()?.mesh;
    if(saved?.inset===.04&&Array.isArray(saved.result?.samples)&&(!completeOnly||saved.result.status==="ok")&&(!exactQuality||saved.quality===resolution)){
      acceptMesh(saved.result,0,true,saved.quality);return true;
    }
    const sharedMesh=window.AtlasNumerics?.getMesh(graph(),{inset:.04});
    if(sharedMesh?.result&&(!completeOnly||sharedMesh.result.status==="ok")&&(!exactQuality||Number(sharedMesh.resolution)>=resolution)){
      acceptMesh(sharedMesh.result,0,true,Number(sharedMesh.resolution));return true;
    }
    return false;
  }
  // Explicit Apply starts one preview sequence. Restoring a family does not.
  function runRecommendedPreview() {
    const key=appliedCounts.join(",");
    if(autoPreview?.key===key&&autoPreview.token===previewEpoch){if(get("decoration-status"))message("decoration-status","Decorations are applied. The recommended preview is already running; checked results stay visible.");return;}
    const profile=window.AtlasPreview?.recommend({family:"lame",counts:[...appliedCounts]})??{label:"Standard preview",quality:8,depth:5,density:1,reason:"Adaptive arithmetic checks remain enabled."};
    const token=++previewEpoch;
    const wanted=Math.max(Number(profile.quality),Number(meshArtifact?.quality)||0,Number(meshCompletion?.quality)||0);
    const qualities=[...get("quality").options].map(option=>Number(option.value));
    get("quality").value=String(qualities.includes(wanted)?wanted:qualities.find(value=>value>=wanted)??qualities.at(-1));
    if(!cuspWorkers.length){overviewDensity=profile.density;overviewRefinement=profile.tailRefinement??0;}
    let depth=Math.max(profile.depth,cuspWorkers.length?Number(get("cusp-depth").value):0);
    if([...get("cusp-depth").options].some(option=>Number(option.value)===depth))get("cusp-depth").value=String(depth);
    const note=get("preview-recommendation");
    if(note)note.textContent=`${profile.label}: up to ${(Number(get("quality").value)+1)*(Number(get("quality").value)+2)/2} interior samples, cusp ε down to 10⁻${String(get("cusp-depth").value).replace(/\d/g,d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d])}, ${overviewDensity} sample${overviewDensity===1?"":"s"} per decade. ${profile.reason} Finer saved images and compatible deeper traces are retained. These settings change sampling; arithmetic precision stays adaptive.`;
    if(get("decoration-status"))message("decoration-status",`Applied (${appliedCounts.join(", ")}) · n = ${order(appliedCounts)}. The interior image runs first, then all six cusp approaches.`);
    saveSelection();
    const requested=JSON.stringify(cuspParameters());
    const job={key,token};autoPreview=job;
    get("cusp-cancel").disabled=false;
    if(!cuspWorkers.length)message("cusp-status","Six cusp approaches are queued after the interior image. Stop series cancels this queued part.");
    const meshReady=meshWorker&&meshCompletion?meshCompletion.promise:computeMesh(undefined,{budgetMs:profile.meshBudgetMs??profile.budgetMs});
    Promise.resolve(meshReady).then(()=>{
      if(disposed||token!==previewEpoch||key!==appliedCounts.join(",")||requested!==JSON.stringify(cuspParameters()))return;
      if(!cuspWorkers.length)startCuspApproaches(undefined,{automatic:true,budgetMs:profile.budgetMs});
    }).catch(error=>{if(!disposed&&token===previewEpoch)message("cusp-status",error.message??"The cusp preview could not start. Computed samples are retained.",true);}).finally(()=>{if(autoPreview===job){autoPreview=null;if(!cuspWorkers.length)get("cusp-cancel").disabled=true;}});
  }
  function applyDecorations(event) {
    event?.preventDefault();
    if (!get("attachment-form").reportValidity()) return;
    const counts = countNames.map(name => get(name).valueAsNumber);
    if (!validCounts(counts) || order(counts) > 21) {
      message(get("decoration-status")?"decoration-status":"attachment-status", "Enter nonnegative whole numbers with a + b + c ≤ 20.", true); return;
    }
    if(event&&counts.every((value,index)=>value===appliedCounts[index])){
      window.SeedCases?.updateLameGraph?.(counts);
      runRecommendedPreview();
      return;
    }
    ++previewEpoch;autoPreview=null;
    rememberResults();
    stopMesh(); ++meshVersion;
    stopCuspApproaches("Choose Compute or extend series to draw the six approaches, or open Limits & Zoom for one focused edge.");
    selectedCounts = counts;appliedCounts=[...counts];
    window.SeedCases?.updateLameGraph?.(counts);
    showSource();
    saveSelection();
    get("result-heading").textContent = `Attachment (${counts.join(", ")}): n = ${order(counts)}`;
    if(get("decoration-status"))message("decoration-status",`Applied (${counts.join(", ")}) · n = ${order(counts)}. Choose the interior mesh or cusp series below.`);
    currentMesh=null;pointResult=null;pointArtifact=null;meshArtifact=null;traceSettings=null;cuspTraces=Array(6).fill(null);
    const saved=savedResults();
    // Old point workers are invalidated even when a restored value needs no job.
    clearTimeout(pointTimer);pointTimer=null;++pointVersion;queuedPoint=null;activePoint=null;pointBusy=false;pointWorker?.terminate();pointWorker=null;
    if(saved?.point?.result&&samePoint(saved.point.source)){
      showPointResult(saved.point.result);
      const complete=saved.point.result.base?.status==="ok"&&saved.point.result.result?.status==="ok";
      message("point-status",complete?"Restored both computed period ratios at this source point.":"Restored a partial point result at this source point. An unresolved branch remains undrawn.",!complete);
    }else requestPoint();
    const tailsRestored=restoreCuspResults(saved?.traces)||restoreSharedTraces();
    if(!tailsRestored){
      updateCuspEndpoints();
      get("cusp-run").textContent="Compute or extend series";
    }
    if(!restoreMeshResults())message("attachment-status", "Decorations applied. Choose Compute interior image to sample the compact image. Cusp series run separately.");
    redraw();
    if(event)runRecommendedPreview();
  }
  function computeMesh(event,options={}) {
    event?.preventDefault();
    if(event)++previewEpoch;
    if (!get("mesh-controls").reportValidity()) return;
    stopMesh(); const id = ++meshVersion;
    const counts=[...appliedCounts], resolution=Number(get("quality").value);
    saveSelection();
    if(restoreMeshResults(true,Boolean(event)))return Promise.resolve();
    let settle;const promise=new Promise(resolve=>{settle=resolve;});meshCompletion={promise,resolve:settle,quality:resolution};
    get("mesh-run").disabled = true; get("cancel").disabled = false;
    message("attachment-status", "Sampling the compact interior image…");
    const started = performance.now();
    try {
      meshWorker = new Worker("numerics-worker.js?v=20261021-overall-final1");
      meshWorker.onmessage = ({data}) => {
        if (disposed || data.id !== meshVersion) return;
        if (data.type === "progress") message("attachment-status", `Compact mesh samples: ${data.completed}/${data.total}…`);
        if (data.type === "done") {
          stopMesh();acceptMesh(data.mesh, (performance.now() - started) / 1000,false,resolution);
        }
        if (data.type === "error") { stopMesh(); message("attachment-status", data.message, true); }
      };
      meshWorker.onerror = () => {
        if (disposed || id !== meshVersion) return;
        stopMesh(); message("attachment-status", "The interior mesh could not run. Point solving is independent.", true);
      };
      meshWorker.postMessage({type: "compute", id, counts, resolution, inset: 0.04,...(options.budgetMs?{budgetMs:options.budgetMs}:{})});
    } catch (error) { stopMesh(); message("attachment-status", error.message, true); }
    return promise;
  }

  listen(get("attachment-form"), "submit", applyDecorations);
  listen(get("mesh-controls"), "submit", computeMesh);
  listen(get("cusp-controls"), "submit", startCuspApproaches);
  listen(get("cusp-cancel"), "click", () => {++previewEpoch;stopCuspApproaches();});
  listen(get("explore-view"), "change", () => {saveSelection(); redraw();});
  listen(get("cancel"), "click", () => {
    ++meshVersion; stopMesh(); message("attachment-status", "Mesh stopped. The draggable point remains active.");
  });
  listen(get("point-form"), "submit", event => {
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
      ++previewEpoch;
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
      listen(get(id + "-" + suffix), "click", event => changeNumber(id, direction, event.shiftKey));
    }
    listen(get(id), "keydown", event => {
      if (!["ArrowUp", "ArrowDown"].includes(event.key)) return;
      event.preventDefault(); changeNumber(id, event.key === "ArrowUp" ? 1 : -1, event.shiftKey);
    });
  }
  listen(get("cusp-position"), "input", () => {
    ++previewEpoch;
    stopCuspApproaches("Paired positions changed. Recompute cusp series to apply them; existing computed samples are retained.");
  });
  for (const svg of [source, attachedSource]) {
    listen(svg, "pointerdown", event => {
      if (event.button != null && event.button !== 0) return;
      pointerId = event.pointerId; pointerSource = svg;
      svg.setPointerCapture(pointerId); pointerPoint(event, svg);
    });
    listen(svg, "pointermove", event => {
      if (pointerSource === svg && pointerId === event.pointerId) pointerPoint(event, svg);
    });
    function finishPointer(event) {
      if (pointerSource !== svg || pointerId !== event.pointerId) return;
      pointerId = null; pointerSource = null; dispatchPoint(); saveSelection();
    }
    listen(svg, "pointerup", finishPointer);
    listen(svg, "pointercancel", finishPointer);
    listen(svg, "lostpointercapture", () => {
      if (pointerSource === svg) { pointerId = null; pointerSource = null; }
    });
    listen(svg, "keydown", event => {
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
  listen(window, "resize", redraw);
  function dispose() {
    if(disposed)return;rememberResults();disposed=true;events.abort();
    clearTimeout(pointTimer);pointTimer=null;
    ++pointVersion;++meshVersion;++cuspVersion;++previewEpoch;autoPreview=null;
    meshCompletion?.resolve();meshCompletion=null;
    pointWorker?.terminate();pointWorker=null;meshWorker?.terminate();meshWorker=null;
    cuspWorkers.forEach(worker=>worker.terminate());cuspWorkers=[];
    queuedPoint=null;activePoint=null;pointBusy=false;
    if(pointerSource&&pointerId!=null&&pointerSource.hasPointerCapture?.(pointerId))pointerSource.releasePointerCapture(pointerId);
    pointerId=null;pointerSource=null;
  }
  showSource(); applyDecorations();
  return {dispose,remember:rememberResults,source:()=>({...point})};
}


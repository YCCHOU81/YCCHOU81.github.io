/* Lamé classic cusp-worker entry. The protocol, checkpoints and depth-extension
 * validation below are unchanged; shared numerical definitions live in workers/. */
"use strict";
importScripts("workers/decimal.js?v=20261025-maintenance1",
  "workers/atlas-math.js?v=20261025-maintenance1",
  "workers/precision-engine.js?v=20261025-maintenance1",
  "workers/trace-initializer.js?v=20261025-maintenance1",
  "workers/approach-engine.js?v=20261025-maintenance1",
  "workers/attachment-anchors.js?v=20261025-maintenance1");
const approachMath=createApproachMath(attachmentMath);
approachMath.setApproachAnchorData(attachmentAnchors);

const TRACE_TIME_BUDGET_MS = 45000;
const ViewDecimal = Decimal.clone({precision: 120, rounding: Decimal.ROUND_HALF_EVEN});

function decorateNode(node) {
  for (const branch of [node.base, node.result]) {
    if (branch?.status !== "ok" || !branch.tauExact) continue;
    const x = new ViewDecimal(branch.tauExact[0]);
    const y = new ViewDecimal(branch.tauExact[1]);
    const denominator0 = x.times(x).plus(y.times(y));
    const u = new ViewDecimal(1).minus(x);
    const denominator1 = u.times(u).plus(y.times(y));
    const halfU = new ViewDecimal(1).minus(x.times(2));
    const denominatorHalf = halfU.times(halfU).plus(y.times(y).times(4));
    const error = new ViewDecimal(branch.numericalErrorEstimate);
    branch.views = {
      raw: [...branch.tauExact],
      zero: [new ViewDecimal(1).minus(x.div(denominator0)).toString(), y.div(denominator0).toString()],
      one: [u.div(denominator1).toString(), y.div(denominator1).toString()],
      half: [x.times(halfU).minus(y.times(y).times(2)).div(denominatorHalf).toString(),y.div(denominatorHalf).toString()]
    };
    branch.viewErrors = {
      raw: error.toString(),
      zero: error.div(denominator0).toString(),
      one: error.div(denominator1).toString(),
      half: error.div(denominatorHalf).toString()
    };
  }
  return node;
}

function validateRequest(data) {
  const counts = data.counts;
  if (!Array.isArray(counts) || counts.length !== 3 ||
      !counts.every(x => Number.isSafeInteger(x) && x >= 0 && x <= 20) ||
      counts.reduce((sum, x) => sum + x, 0) > 20) {
    throw new RangeError("Enter three nonnegative whole numbers with a + b + c <= 20.");
  }
  const path = data.path;
  if (!path || !["edge", "vertex"].includes(path.kind) ||
      !Number.isInteger(path.index) || path.index < 0 || path.index > 2) {
    throw new RangeError("Choose an edge or vertex numbered 0, 1, or 2.");
  }
  const position = String(path.position ?? ".5");
  if (position.length > 60 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(position)) {
    throw new RangeError("Enter a decimal path position strictly between 0 and 1.");
  }
  const p = new Decimal(position);
  if (!p.isFinite() || !p.gt(0) || !p.lt(1)) {
    throw new RangeError("The path position must lie strictly between 0 and 1.");
  }
  const normalizedPath = approachMath.normalizePath({...path,position});
  const target = String(data.target ?? "1e-5");
  if (target.length > 40 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(target)) {
    throw new RangeError("Choose a target epsilon between 0.1 and 1e-10.");
  }
  const targetDecimal = new Decimal(target);
  if (!targetDecimal.isFinite() || targetDecimal.lt("1e-10") || targetDecimal.gt(".1")) {
    throw new RangeError("Choose a target epsilon between 0.1 and 1e-10.");
  }
  const startEpsilon = normalizedPath.start ? String(data.startEpsilon ?? "1") : (data.startEpsilon == null ? null : String(data.startEpsilon));
  if (startEpsilon !== null) {
    if (startEpsilon.length > 40 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(startEpsilon)) {
      throw new RangeError(normalizedPath.start ? "A custom segment starts at epsilon 1." : "Choose a starting epsilon above the target and at most 0.3.");
    }
    const startDecimal = new Decimal(startEpsilon);
    if (!startDecimal.isFinite() || !startDecimal.gt(targetDecimal) || (normalizedPath.start ? !startDecimal.eq(1) : startDecimal.gt(".3"))) {
      throw new RangeError(normalizedPath.start ? "A custom segment starts at epsilon 1." : "Choose a starting epsilon above the target and at most 0.3.");
    }
  }
  const samplesPerDecade = data.samplesPerDecade ?? 2;
  if (!Number.isInteger(samplesPerDecade) || samplesPerDecade < 1 || samplesPerDecade > 4) {
    throw new RangeError("Choose between 1 and 4 samples per decade.");
  }
  const tailRefinement = data.tailRefinement ?? 0;
  if (!Number.isInteger(tailRefinement) || tailRefinement < 0 || tailRefinement > 2) {
    throw new RangeError("Choose between 0 and 2 extra tail samples.");
  }
  return {
    counts: [...counts],
    path: normalizedPath,
    target,
    samplesPerDecade,
    tailRefinement,
    startEpsilon
  };
}

function pathSignature(path){
  return JSON.stringify([path.kind,path.index,new ViewDecimal(path.position).toString(),path.start??null]);
}
function decimalText(value,label){
  if(typeof value!=='string'||value.length>240||!/^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(value))throw new RangeError('The saved '+label+' is not valid decimal data. Compute a new approach.');
  const number=new ViewDecimal(value);
  if(!number.isFinite())throw new RangeError('The saved '+label+' is not finite. Compute a new approach.');
  return number;
}
function validateExtension(extension,request){
  if(!extension||typeof extension!=='object')throw new RangeError('The saved approach is missing. Compute a new approach.');
  const previous=validateRequest(extension.request??{}),result=extension.result;
  const sameStart=(a,b)=>a===null&&b===null||a!==null&&b!==null&&new ViewDecimal(a).eq(b);
  if(JSON.stringify(previous.counts)!==JSON.stringify(request.counts)||pathSignature(previous.path)!==pathSignature(request.path)||
    !sameStart(previous.startEpsilon,request.startEpsilon)||previous.samplesPerDecade!==request.samplesPerDecade||previous.tailRefinement!==request.tailRefinement||
    !new ViewDecimal(request.target).lt(previous.target))throw new RangeError('Only the depth may increase when continuing a saved approach. Compute a new approach for other changes.');
  if(!result||!['ok','partial','aborted','stopped'].includes(result.status)||!Array.isArray(result.nodes)||!result.nodes.length||result.nodes.length>99||
    (result.counts!==undefined&&JSON.stringify(result.counts)!==JSON.stringify(request.counts))||
    (result.path!==undefined&&pathSignature(approachMath.normalizePath(result.path))!==pathSignature(request.path)))
    throw new RangeError('The saved approach does not match this case and path. Compute a new approach.');
  const stopKind=result.stop?.kind??result.stopReason;
  if(stopKind&&/numerical|verification|precision|correction|branch/i.test(stopKind))throw new RangeError('This saved approach stopped at a numerical check and cannot be extended past it. Compute a new approach.');
  if(result.status!=='ok'&&(!['time_budget','cancelled','canceled','aborted','page_change','family_change'].includes(stopKind)&&!(result.status==='stopped'&&!stopKind)||/unresolved|(?:verification|numerical|precision|correction|branch).{0,40}(?:fail|stop|limit)/i.test(result.message??'')))
    throw new RangeError('This saved approach stopped at a numerical check and cannot be extended past it. Compute a new approach.');
  const wanted=result.requestedEpsilons,expectedN=1+request.counts.reduce((a,b)=>a+b,0),nodes=result.nodes;
  const start=new ViewDecimal(previous.startEpsilon??'.1'),oldTarget=new ViewDecimal(previous.target),largestRatio=ViewDecimal.pow(10,new ViewDecimal(1).div(previous.samplesPerDecade)).times('1.00000000000001');
  let preceding=null;
  for(let i=0;i<nodes.length;i++){
    const node=nodes[i];
    if(!node||node.status!=='ok')throw new RangeError('Only a contiguous prefix of independently checked samples can be extended.');
    const epsilon=decimalText(node.epsilon,'epsilon');
    if(!epsilon.gte(oldTarget)||epsilon.gt(start)||(i===0&&!epsilon.eq(start))||(preceding&&!epsilon.lt(preceding)))throw new RangeError('The saved samples are not a continuous ordered approach prefix.');
    if(preceding&&preceding.div(epsilon).gt(largestRatio))throw new RangeError('The saved approach skips a requested interval. Compute a new approach.');
    if(wanted&&(!Array.isArray(wanted)||i>=wanted.length||!epsilon.eq(decimalText(String(wanted[i]),'requested epsilon'))))throw new RangeError('The saved samples do not match their requested approach grid.');
    const point=approachMath.approachCoordinates(request.path,node.epsilon,70),r=decimalText(node.rExact,'source coordinate'),s=decimalText(node.sExact,'source coordinate');
    if(r.minus(point.rExact).abs().gt('1e-65')||s.minus(point.sExact).abs().gt('1e-65')||!Number.isFinite(node.r)||!Number.isFinite(node.s)||Math.abs(node.r-point.r)>1e-14||Math.abs(node.s-point.s)>1e-14||!Array.isArray(node.barycentric)||node.barycentric.length!==3||
      node.barycentric.some((value,j)=>decimalText(value,'barycentric coordinate').minus(point.barycentric[j]).abs().gt('1e-65')))
      throw new RangeError('The saved samples use a different source path. Compute a new approach.');
    for(const [branch,n] of [[node.base,1],[node.result,expectedN]]){
      if(branch?.status!=='ok'||branch.n!==n||branch.validation!=='two_precisions_not_interval_certified'||!Array.isArray(branch.tauExact)||branch.tauExact.length!==2||
        !(branch.precisionAgreementDigits>=17)||!Number.isSafeInteger(branch.precisionDigits)||branch.precisionDigits<60||!Number.isSafeInteger(branch.validationDigits)||branch.validationDigits<branch.precisionDigits+25)
        throw new RangeError('The saved samples lack independent branch checks. Compute a new approach.');
      const re=decimalText(branch.tauExact[0],'tau coordinate'),im=decimalText(branch.tauExact[1],'tau coordinate'),error=decimalText(branch.numericalErrorEstimate,'error estimate');
      if(!re.gt(0)||!re.lt(1)||!im.gt(0)||re.minus('.5').pow(2).plus(im.pow(2)).lte('.25')||error.lt(0)||error.gt('1e-17'))
        throw new RangeError('The saved branch lies outside the marked domain or lacks the required accuracy.');
      if(!Number.isFinite(branch.tau?.re)||!Number.isFinite(branch.tau?.im)||Math.abs(branch.tau.re-re.toNumber())>1e-14*Math.max(1,re.abs().toNumber())||Math.abs(branch.tau.im-im.toNumber())>1e-14*Math.max(1,im.abs().toNumber()))
        throw new RangeError('The saved plotting coordinates do not match their checked decimal values.');
    }
    preceding=epsilon;
  }
  if(result.status==='ok'&&!preceding.eq(oldTarget))throw new RangeError('The saved complete approach did not reach its stated target.');
  if(!new ViewDecimal(request.target).lt(preceding))throw new RangeError('Choose a deeper target to continue this approach.');
  return{request:previous,result};
}
function requestedGrid(request,extension=null){
  const prefix=extension?.result.nodes.map(node=>node.epsilon)??[];
  const out=prefix.length?[...prefix]:request.startEpsilon===null?approachMath.defaultApproachEpsilons(request.target,request.samplesPerDecade):[];
  if(prefix.length||request.startEpsilon!==null){
    const start=new ViewDecimal(prefix.at(-1)??request.startEpsilon),target=new ViewDecimal(request.target);
    if(!prefix.length)out.push(start.toString());
    for(let k=1;k<100;k++){
      const epsilon=start.times(ViewDecimal.pow(10,new ViewDecimal(-k).div(request.samplesPerDecade)));
      if(epsilon.lte(target))break;
      out.push(epsilon.toPrecision(30));
    }
    out.push(target.toString());
  }
  if(request.tailRefinement&&out.length>1){
    const target=new ViewDecimal(out.at(-1)),previous=new ViewDecimal(out.at(-2)),ratio=target.div(previous);
    const extra=Array.from({length:request.tailRefinement},(_,index)=>previous.times(ratio.pow(new ViewDecimal(index+1).div(request.tailRefinement+1))).toPrecision(60));
    out.splice(out.length-1,0,...extra);
  }
  if(out.length>100)throw new RangeError('This extended approach exceeds the 100-sample budget. Compute a new approach.');
  return out;
}

let traceRunning = false;
// Only complete, independently checked traces enter this small worker-local
// cache. Each exact request includes its marking, path, range, and density.
const completedTraceCache = new Map();
const traceCheckpoints=new Map();let checkpointSerial=0;

self.onmessage = async ({data}) => {
  if (!data || data.type !== "trace") return;
  const id = data.id;
  if (traceRunning) {
    self.postMessage({type: "error", id, message: "A trace is already running in this worker."});
    return;
  }
  let request;
  try {
    request = validateRequest(data);
  } catch (error) {
    self.postMessage({type: "error", id, message: error.message});
    return;
  }
  const cacheKey = JSON.stringify([request.counts,request.path.kind,request.path.index,
    new ViewDecimal(request.path.position).toString(),new ViewDecimal(request.target).toString(),
    request.samplesPerDecade,request.startEpsilon===null?null:new ViewDecimal(request.startEpsilon).toString(),
    request.path.start??null,request.tailRefinement]);
  if(completedTraceCache.has(cacheKey)){
    const result=completedTraceCache.get(cacheKey);
    completedTraceCache.delete(cacheKey);completedTraceCache.set(cacheKey,result);
    result.nodes.forEach((node,index)=>self.postMessage({type:'progress',id,node,nodes:result.nodes.slice(0,index+1),completed:index+1,total:result.nodes.length}));
    self.postMessage({type:'done',id,result:{...result,fromCache:true,reusedCount:result.nodes.filter(node=>node.status==='ok').length,newCount:0}});return;
  }
  const saved=data.resumeToken?traceCheckpoints.get(String(data.resumeToken)):null;
  if(data.resumeToken&&(!saved||saved.cacheKey!==cacheKey)){
    self.postMessage({type:'error',id,message:'The saved approach does not match these settings. Compute a new approach.'});return;
  }
  let extension=null,epsilons;
  try{
    if(data.extendFrom){
      if(saved)throw new RangeError('Use either the current live checkpoint or a deeper saved approach, not both.');
      extension=validateExtension(data.extendFrom,request);
    }
    // A same-target live resume must keep the exact grid of that attempt,
    // including earlier refinement nodes retained by a depth extension.
    epsilons=saved?.result?.requestedEpsilons?[...saved.result.requestedEpsilons]:requestedGrid(request,extension);
  }catch(error){self.postMessage({type:'error',id,message:error.message});return;}
  traceRunning = true;
  const started = performance.now();
  const nodes = [];
  const timeBudgetMs=Math.max(250,Math.min(TRACE_TIME_BUDGET_MS,Number(request.budgetMs)||TRACE_TIME_BUDGET_MS));
  const signal = {get aborted() {return performance.now() - started > timeBudgetMs;}};
  try {
    let result = await approachMath.solveApproachPath(request.counts, request.path, epsilons, {
      anchorData: attachmentAnchors,
      resumeState:saved?.result??extension?.result??null,
      revalidateRetainedSeed:Boolean(extension),
      precisionDigits: 70,
      validationDigits: 100,
      signal,
      onProgress(node, completed, total) {
        decorateNode(node);
        nodes.push(node);
        self.postMessage({type: "progress", id, node, nodes:[...nodes], completed, total});
      }
    });
    if (result.status === "aborted") {
      result = {
        ...result,
        status: "partial",
        n: 1 + request.counts.reduce((sum, x) => sum + x, 0),
        completed: result.nodes.filter(node => node.status === "ok").length,
        total: epsilons.length,
        requestedEpsilons: epsilons,
        precisionDigits: 70,
        validationDigits: 100,
        stopReason: "time_budget",
        message: `The ${Math.round(timeBudgetMs/1000)}-second continuation budget was reached. Validated samples are retained; the target may not have been reached.`
      };
    }
    result.nodes.forEach(decorateNode);
    result.tailRefinement = request.tailRefinement;
    result.reusedCount=(saved?.result??extension?.result)?.nodes.filter(node=>node.status==='ok').length??0;
    result.newCount=Math.max(0,result.nodes.filter(node=>node.status==='ok').length-result.reusedCount);
    result.extended=Boolean(extension||saved?.result.extended);
    result.resumed=Boolean(saved);result.attempts=(saved?.attempts??0)+1;
    if(result.status==='partial'){
      const timeBudget=result.stopReason==='time_budget';
      const validated=result.nodes.filter(node=>node.status==='ok');
      result.stopReason=result.stopReason??'verification_stopped';
      result.message=result.message??result.nodes.find(node=>node.status!=='ok')?.result?.message??'Continuation stopped when this branch no longer passed the numerical checks. Validated points are retained.';
      result.stop={kind:timeBudget?'time_budget':'numerical_verification',failedEpsilon:result.stoppedAt??result.requestedEpsilons?.[validated.length]??null,completed:validated.length,total:result.total,retryable:timeBudget,resumeAvailable:timeBudget&&validated.length>0};
      result.resumeAvailable=result.stop.resumeAvailable;
      if(result.resumeAvailable){
        const token=data.resumeToken??('classic-'+Date.now().toString(36)+'-'+(++checkpointSerial));
        traceCheckpoints.set(token,{cacheKey,result,attempts:result.attempts});result.resumeToken=token;
        if(traceCheckpoints.size>6)traceCheckpoints.delete(traceCheckpoints.keys().next().value);
      }
    }else if(data.resumeToken){traceCheckpoints.delete(String(data.resumeToken));}

    if(result.status==='ok'){
      completedTraceCache.set(cacheKey,result);
      if(completedTraceCache.size>12)completedTraceCache.delete(completedTraceCache.keys().next().value);
    }
    self.postMessage({type: "done", id, result});
  } catch (error) {
    if (nodes.length) {
      self.postMessage({type: "done", id, result: {
        status: "partial", counts: request.counts,
        n: 1 + request.counts.reduce((sum, x) => sum + x, 0),
        path: request.path, nodes, completed: nodes.filter(node => node.status === "ok").length,
        total: epsilons.length, requestedEpsilons: epsilons,
        tailRefinement: request.tailRefinement,
        precisionDigits: 70, validationDigits: 100,
        elapsedMs: performance.now() - started,
        reusedCount:(saved?.result??extension?.result)?.nodes.filter(node=>node.status==='ok').length??0,
        newCount:Math.max(0,nodes.filter(node=>node.status==='ok').length-((saved?.result??extension?.result)?.nodes.filter(node=>node.status==='ok').length??0)),extended:Boolean(extension||saved?.result.extended),
        stopReason: "numerical_error", message: error.message, stop:{kind:"numerical_error",completed:nodes.filter(node=>node.status==="ok").length,total:epsilons.length,retryable:false,resumeAvailable:false},resumeAvailable:false,resumed:Boolean(saved),attempts:(saved?.attempts??0)+1
      }});
    } else {
      self.postMessage({type: "error", id, message: error.message});
    }
  } finally {
    traceRunning = false;
  }
};

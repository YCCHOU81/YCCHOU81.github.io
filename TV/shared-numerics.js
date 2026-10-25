/* Share numerical artifacts by marked graph and actual path. Display rounding
   and view settings do not enter the keys; stored decimal strings are untouched. */
(() => {
  "use strict";
  const cache=()=>window.AtlasResults;
  const namespaces={mesh:"numeric-mesh-v1",point:"numeric-point-v1",trace:"numeric-traces-v1",baseMesh:"numeric-base-mesh-v1"};
  const baseMeshVersion="tv-base-mesh-checked-1",viewNames=["raw","zero","one","half"];
  const baseFamily=g=>typeof g==="string"?g:g?.family;
  function checkedBase(branch,family) {
    const precision=Number(branch?.precisionDigits??branch?.precision);
    if(!["double","bouquet"].includes(family)||!branchGood(branch)||!Number.isFinite(precision)||precision<50)return false;
    if(family==="bouquet"&&!(branch.endpointChecked===true&&branch.lowRequestedDigits===18&&branch.highRequestedDigits===28))return false;
    const threshold=family==="double"?1e-10:1e-12;
    return viewNames.every(name=>branch.viewErrors?.[name]!==null&&branch.viewErrors?.[name]!==undefined&&Number.isFinite(Number(branch.viewErrors[name]))&&Number(branch.viewErrors[name])>=0&&Number(branch.viewErrors[name])<=threshold);
  }
  function validBaseMesh(artifact,family,inset) {
    return artifact?.version===baseMeshVersion&&artifact.family===family&&artifact.inset===inset&&
      Number.isInteger(artifact.resolution)&&artifact.resolution>=3&&Number.isFinite(inset)&&inset>=.025&&inset<1/3&&
      Array.isArray(artifact.points)&&Array.isArray(artifact.samples)&&artifact.points.length===artifact.samples.length&&artifact.samples.length>0&&
      artifact.points.every(point=>Number.isFinite(point?.r)&&Number.isFinite(point?.s))&&artifact.samples.every(branch=>checkedBase(branch,family));
  }
  function baseQuality(artifact) {
    return {count:artifact.samples.length,precision:Math.min(...artifact.samples.map(branch=>Number(branch.precisionDigits??branch.precision))),
      error:Math.max(...artifact.samples.flatMap(branch=>viewNames.map(name=>Number(branch.viewErrors[name]))))};
  }
  function putBaseMesh(g,artifact) {
    const family=baseFamily(g),result=artifact?.result;
    if(!cache()||!["double","bouquet"].includes(family)||!Array.isArray(result?.points)||!Array.isArray(result.samples)||result.points.length!==result.samples.length)return false;
    // The worker may cap the requested quality. The returned grid owns the
    // identity, so a requested twelve cannot mislabel an actual six grid.
    const resolution=result.resolution,inset=result.inset,points=[],samples=[];
    for(let index=0;index<result.points.length;index++)if(checkedBase(result.samples[index]?.base,family)){
      points.push(result.points[index]);samples.push(result.samples[index].base);
    }
    const saved={version:baseMeshVersion,family,inset,resolution,points,samples};
    if(!validBaseMesh(saved,family,inset))return false;
    const key=cache().key({version:baseMeshVersion,family,inset,resolution}),old=cache().get(namespaces.baseMesh,key);
    if(validBaseMesh(old,family,inset)){
      const a=baseQuality(old),b=baseQuality(saved);
      if(a.count>b.count||a.count===b.count&&a.precision>=b.precision&&a.error<=b.error)return true;
    }
    return cache().set(namespaces.baseMesh,key,saved);
  }
  function decimal(value) {
    if(value===null||value===undefined)return null;
    const text=String(value).trim(),m=/^([+-]?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?$/i.exec(text);
    if(!m||!(m[2]+(m[3]??"")).length)return text;
    let digits=(m[2]+(m[3]??"")).replace(/^0+/,""),power=Number(m[4]??0)-(m[3]??"").length;
    if(!digits)return "0";
    while(digits.endsWith("0")){digits=digits.slice(0,-1);power++;}
    return (m[1]==="-"?"-":"")+digits+"e"+power;
  }
  function graph(value) {
    return value?.family==="lame" ? {family:"lame",counts:value.counts}
      : {family:value?.family,caseId:value?.caseId,caseConfig:value?.caseConfig??null};
  }
  const graphKey=g=>cache()?.key(graph(g))??JSON.stringify(graph(g));
  function traceKey(request) {
    const path=request.path??{};
    return cache().key({path:{kind:path.kind,index:Number(path.index),position:decimal(path.position),start:Array.isArray(path.start)?path.start.map(decimal):null},
      target:decimal(request.target),startEpsilon:decimal(request.startEpsilon??(path.start?1:null)),
      samplesPerDecade:Number(request.samplesPerDecade??2),tailRefinement:Number(request.tailRefinement??0)});
  }
  const goodCount=result=>(result?.nodes??[]).filter(node=>node.status==="ok").length;
  const branchGood=branch=>branch?.status==="ok"&&Array.isArray(branch.tauExact)&&branch.tauExact.length===2&&branch.tauExact.every(x=>Number.isFinite(Number(x)))&&Number(branch.tauExact[1])>0;
  function extensionPrefix(result) {
    if(!Array.isArray(result?.nodes)||!result.nodes.length||result.tailBranchValidationStatus==="pending")return null;
    const stop=result.stop?.kind??"",failure=String(stop)+" "+String(result.message??"");
    if(result.branchRejected===true||result.stop?.branchRejected===true)return null;
    if(/numerical|verification|precision|correction|branch|error|aborted/i.test(stop))return null;
    if(/fail|reject|mismatch|verification|did not agree/i.test(failure))return null;
    if(!["ok","partial","stopped"].includes(result.status))return null;
    const userStopped=result.status==="stopped"&&!result.stop;
    if(result.status!=="ok"&&!["time_budget","cancelled","canceled","page_change","family_change"].includes(stop)&&!result.cancelled&&!userStopped)return null;
    const prefix=[];let previous=Infinity,stopped=false;
    for(const node of result.nodes){
      const good=node?.status==="ok"&&branchGood(node.base)&&branchGood(node.result);
      if(!good){stopped=true;continue;}
      if(stopped)return null;
      const epsilon=Number(node.epsilon);
      if(!(epsilon>0&&epsilon<previous))return null;
      previous=epsilon;prefix.push(node);
    }
    return prefix.length?prefix:null;
  }
  const meshScore=artifact=>(artifact?.result?.samples??[]).filter(node=>node?.result?.status==="ok").length;
  window.AtlasNumerics=Object.freeze({
    graphKey,traceKey,decimalKey:decimal,
    mergeNodes(prefix,incoming){
      const merged=new Map((incoming??[]).filter(Boolean).map(node=>[decimal(node.epsilon),node]));
      for(const node of prefix??[])merged.set(decimal(node.epsilon),node);
      return [...merged.values()].sort((a,b)=>Number(b.epsilon)-Number(a.epsilon));
    },
    putMesh(g,artifact) {
      if(!cache()||!Array.isArray(artifact?.result?.samples)||!Array.isArray(artifact.result.triangles)||!Array.isArray(artifact.result.points))return false;
      putBaseMesh(g,artifact);
      const key=graphKey(g),old=cache().get(namespaces.mesh,key);
      if(old&&old.inset===artifact.inset&&(meshScore(old)>meshScore(artifact)||meshScore(old)===meshScore(artifact)&&Number(old.resolution)>Number(artifact.resolution)))return true;
      return cache().set(namespaces.mesh,key,artifact);
    },
    getMesh(g,{inset=.04}={}){const saved=cache()?.get(namespaces.mesh,graphKey(g));return saved?.inset===inset?saved:null;},
    getBaseMesh(g,{inset=.04}={}) {
      const family=baseFamily(g);if(!["double","bouquet"].includes(family))return null;
      let best=null;
      for(const {value} of cache()?.list(namespaces.baseMesh)??[])if(validBaseMesh(value,family,inset)&&
        (!best||value.resolution>best.resolution||value.resolution===best.resolution&&value.samples.length>best.samples.length))best=value;
      return best;
    },
    putPoint(g,artifact){return Array.isArray(artifact?.source)&&artifact.source.length===2&&artifact.result?cache()?.set(namespaces.point,graphKey(g),artifact)??false:false;},
    getPoint(g,source){const saved=cache()?.get(namespaces.point,graphKey(g));return Array.isArray(saved?.source)&&cache().key(saved.source.map(decimal))===cache().key(source.map(decimal))?saved:null;},
    putTrace(g,request,result) {
      if(!cache()||!request?.path||!Array.isArray(result?.nodes)||!result.nodes.length)return false;
      const key=graphKey(g),identity=traceKey(request),stored=cache().get(namespaces.trace,key),saved={traces:Array.isArray(stored?.traces)?stored.traces.filter(trace=>trace?.request?.path&&Array.isArray(trace.result?.nodes)):[]};
      const index=saved.traces.findIndex(trace=>traceKey(trace.request)===identity),old=saved.traces[index];
      if(old&&(old.result.status==="ok"&&result.status!=="ok"||goodCount(old.result)>goodCount(result)))return true;
      if(index>=0)saved.traces.splice(index,1);
      saved.traces.push({request,result});saved.traces=saved.traces.slice(-12);
      return cache().set(namespaces.trace,key,saved);
    },
    getTrace(g,request){return this.traces(g).find(trace=>traceKey(trace.request)===traceKey(request))??null;},
    getExtension(g,request){
      const identity=traceKey({...request,target:null}),target=Number(request.target);
      if(!(target>0))return null;
      let best=null,last=Infinity;
      for(const trace of this.traces(g)){
        const priorTarget=Number(trace.request.target);
        if(traceKey({...trace.request,target:null})!==identity||!(priorTarget>=target)||priorTarget===target&&trace.result.status==="ok")continue;
        const prefix=extensionPrefix(trace.result),end=prefix?Number(prefix.at(-1).epsilon):NaN;
        if(!(end>target&&end<last))continue;
        const original=g.family==="lame"?{counts:g.counts}:{seed:g.family,caseId:g.caseId,caseConfig:g.caseConfig??null};
        best={request:{...trace.request,...original},result:{...trace.result,nodes:prefix}};last=end;
      }
      return best;
    },
    traces(g){const saved=cache()?.get(namespaces.trace,graphKey(g));return Array.isArray(saved?.traces)?saved.traces.filter(trace=>trace?.request?.path&&Array.isArray(trace.result?.nodes)):[];}
  });
})();

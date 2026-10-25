import {cachedBaseBranch,verifiedBaseMesh,baseMeshCovers} from './base-mesh-cache.mjs?v=20261024-tv-audit2';
import {reusableMeshSamples} from './compact-mesh-reuse.mjs?v=20261024-tv-audit2';
import * as C from './bouquet-continuation.mjs?v=20261024-tv-audit2';
import {anchor as geometricAnchor} from './bouquet-anchor.mjs?v=20261024-tv-audit2';
import {createBouquetTrace} from './bouquet-trace.mjs?v=20261024-tv-audit2';
import * as F from './generic-fast.mjs';
import {createEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';
import {fixedSchur} from './fixed-schur.mjs';
import verifiedAnchors from './bouquet-verified-anchors.mjs';
import {traceBudget} from './trace-checkpoints.mjs?v=20261024-audit1';
import {checkedCompactMeshRoot,compactMeshCorrect,compactMeshPlotRoot} from './bouquet-mesh-corrector.mjs?v=20261024-tv-audit2';
import {pooledCompactMesh} from './compact-mesh-pool.mjs?v=20261024-tv-audit2';

const engines=new Map(),anchors=new Map(),compactAnchors=new Map(),clock=()=>performance.now(),pause=()=>new Promise(r=>setTimeout(r,0));
const engine=p=>{p=Math.max(35,Math.min(950,Math.ceil(p/5)*5));if(!engines.has(p))engines.set(p,createEvaluator(p));return engines.get(p)};
const CH=[[1,0,0,1],[1,-1,1,0],[0,1,-1,1]],PERM=[[0,1,2,3],[0,2,3,1],[0,3,1,2]];
const inside=(R,S)=>Number(R)>0&&Number(R)<.5&&Number(S)>0&&Number(S)<.5&&Number(R)+Number(S)>.5;
const distance=(a,b)=>Math.hypot(Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]));
const baseOf=c=>({id:'bouquet-base1000',seed:'bouquet',core:'bouquet_B0_T0',weights:[1,0,0,0],shift:[0,0],graph:c.baseGraph??c.graph});
const isBase=c=>c.weights.join(',')==='1,0,0,0';
const checkDeadline=d=>{if(clock()>d)throw Error('The 45-second numerical budget was reached; the last resolved points are retained')};
function validate(c){if(c.seed!=='bouquet'||c.weights?.length!==4||c.weights.some(n=>!Number.isInteger(n)||n<0)||!c.graph)throw Error('The labelled three-loop configuration is incomplete');if(c.weights.reduce((a,b)=>a+b,0)>144||Math.max(...c.weights)>72)throw Error('This graph exceeds the local numerical budget (144 hemispheres / weight 72). Reduce the counts.');return c}
function correct(c,source,guess,{precision=70,digits=18,edge=0,deadline=Infinity,radius=.15}={}){
 checkDeadline(deadline);let got;
 for(const extra of[0,30,60]){
  got=C.correct({config:C.cleanConfig(c),source,guess,edge,precision:Math.max(70,precision)+extra,digits,radius:Math.min(.2,radius),deadline,jetScale:16});
  if(got.converged)return{...got,tau:got.tau,precision:got.precision};
  if(!/precision|rounded|singular|derivative|iteration|decreasing/i.test(got.message))break;
 }
 const error=Error(got?.message??'The local correction needs more precision');error.corrector=got;throw error;
}
function fastCorrect(c,source,guess,{deadline=Infinity,metrics=null}={}){
 checkDeadline(deadline);guess=guess.map(Number);let reference;
 // fixedSchur builds its reference at guess and immediately evaluates that
 // same matrix. Reuse those exact entries; its elimination copies the rows.
 const get=w=>{if(reference&&w[0]===guess[0]&&w[1]===guess[1])return reference;if(metrics)metrics.fastMatrixEvaluations++;const M=F.tvMatrix(c.weights,w,Number(source[1])+c.shift[1],Number(source[0])+c.shift[0],{reduce:false}).M;if(w[0]===guess[0]&&w[1]===guess[1])reference=M;return M;},residual=fixedSchur(F,get,guess);
 let a=guess,b=F.add(a,[1e-5,2e-5]),fa=residual.evaluate(a),fb=residual.evaluate(b),best=null;
 for(let i=0;i<24;i++){
  const step=F.div(F.mul(fb,F.sub(b,a)),F.sub(fb,fa)),err=F.abs(step);if(!Number.isFinite(err))throw Error('A preview correction needs more precision');
  if(!best||err<best.err)best={tau:b,err};
  if(err<2e-8)return{tau:b.map(String),precision:15,correction:String(err),iterations:i+2};
  a=b;fa=fb;b=F.sub(b,step);if(b[1]<=.001||F.abs(F.sub(b,guess))>.15*(1+F.abs(guess)))throw Error('The preview correction left the selected branch');fb=residual.evaluate(b);
 }
 if(best&&best.err<2e-7)return{tau:best.tau.map(String),precision:15,correction:String(best.err),iterations:25};
 throw Error('The preview correction requires higher precision');
}
function branch(tau,diagnostics={}){
 const {tau:ignoredTau,comparisonTauExact,...extra}=diagnostics;
 const E=engine(Math.max(100,(diagnostics.precisionDigits??diagnostics.precision??70)+20)),t=E.c(tau),views={raw:t,zero:E.sub(1,E.div(1,t)),one:E.div(1,E.sub(1,t)),half:E.div(t,E.sub(1,E.mul(2,t)))};
 const lower=comparisonTauExact?E.c(comparisonTauExact):null,other=lower?{raw:lower,zero:E.sub(1,E.div(1,lower)),one:E.div(1,E.sub(1,lower)),half:E.div(lower,E.sub(1,E.mul(2,lower)))}:null;
 const viewErrors=other?Object.fromEntries(Object.keys(views).map(k=>[k,E.abs(E.sub(views[k],other[k])).toString()])):null;
 return{status:'ok',tau:{re:Number(tau[0]),im:Number(tau[1])},tauExact:tau.map(String),views:Object.fromEntries(Object.entries(views).map(([k,v])=>[k,E.strs(v)])),precisionDigits:diagnostics.precision??15,numericalErrorEstimate:viewErrors?.raw??diagnostics.correction??null,...extra,...(viewErrors?{viewErrors,crossPrecisionTauError:viewErrors.raw}:{}),engine:"stable Fourier matching and bounded graph continuation",formulaInputsUsed:false};
}
function cacheAnchor(c,record,{compact=false}={}){const key=C.configKey(c),cache=compact?compactAnchors:anchors;cache.delete(key);cache.set(key,record);while(cache.size>24)cache.delete(cache.keys().next().value);return record;}
async function anchor(c,deadline,{compact=false}={}){
 const key=C.configKey(c),cache=compact?compactAnchors:anchors;if(cache.has(key))return cache.get(key);checkDeadline(deadline);
 const counts=c.graph.stable_edges?.map(e=>e.count??e.q??0)??[],known=counts.every(q=>q===0)?verifiedAnchors.cores[c.core??'bouquet_B0_T0']:verifiedAnchors.counts[(c.core??c.graph.core_id)+':'+counts.join(',')];
 if(known&&known.weights.join(',')===c.weights.join(',')){
  const source=known.source.map(String),got=compact?checkedCompactMeshRoot(C.cleanConfig(c),source,known.tau,{deadline,radius:.08,agreementTolerance:'1e-18'}):await C.repeatedRoot(C.cleanConfig(c),source,known.tau,{precision:80,digits:18,deadline,jetScale:16,radius:.08});
  return cacheAnchor(c,{source,tau:compact?got.tau:got.tauExact,precision:compact?got.precision:got.precisionDigits,accepted:true,checked:true,configKey:key,formulaInputsUsed:false,provenance:known.provenance??'Previously checked labelled graph, independently refined'},{compact});
 }
 // The compact callback changes arithmetic only. Both independent physical
 // meshes, their health tests, and the anchor's precision/repeat policy stay
 // inside geometricAnchor. Its namespace preserves the standard cache path.
 const compactCorrect=options=>{const got=compactMeshCorrect(options.config,options.source,options.guess,options);return got.converged?got:C.correct(options);};
 const got=await geometricAnchor(C.cleanConfig(c),['.29','.31'],{deadline,precision:90,digits:18,...(compact?{correct:compactCorrect,cacheScope:'mesh-fast-v1'}:{})});
 if(!got.accepted){const error=Error(got.message);error.failureKind=got.failureKind;throw error;}
 return cacheAnchor(c,{...got,provenance:'Healthy locally regularized marked Hodge meshes with bounded root agreement'},{compact});
}
function transport(c,from,to,guess,{edge=0,precision=0,digits=25,deadline=Infinity,depth=0,maxDistance=.012,metrics=null}={}){
 checkDeadline(deadline);if(depth>20)throw Error('The branch needs a shorter continuation step');const E=engine(Math.max(65,precision)),dist=distance(from,to);
 const split=()=>{const m=from.map((x,i)=>new E.D(x).plus(to[i]).div(2).toString()),a=transport(c,from,m,guess,{edge,precision,digits,deadline,depth:depth+1,maxDistance,metrics});return transport(c,m,to,a.tau,{edge,precision,digits,deadline,depth:depth+1,maxDistance,metrics})};
 if(dist>maxDistance)return split();
 try{if(metrics)metrics.localSteps++;let out;if(precision)out=correct(c,to,guess,{edge,precision,digits,deadline});else try{out=fastCorrect(c,to,guess,{deadline,metrics})}catch{if(metrics)metrics.fixedFallbacks++;out=correct(c,to,guess,{edge,precision:45,digits:20,deadline})}const a=E.mobius(E.c(guess),CH[edge]),b=E.mobius(E.c(out.tau),CH[edge]);if(E.abs(E.sub(a,b)).gt(E.abs(a).plus(1).times('.13')))throw Error('The branch needs a smaller step');if(metrics)metrics.acceptedLocalSteps++;return out}
 catch(error){if(metrics)metrics.rejectedLocalSteps++;if(depth>=20||clock()>deadline||dist<1e-16)throw error;return split()}
}

export async function bouquetPoint(request,metrics=null){
 const c=validate(request.caseConfig),source=[String(request.R),String(request.S)],begun=clock(),deadline=request.deadline??begun+45000;if(!inside(...source))throw Error('Choose a point strictly inside 0 < R,S < 1/2, R+S > 1/2');const bc=baseOf(c),compact=request._meshChecked===true,cachedPage=compact?cachedBaseBranch(request.baseMesh,'bouquet',source,Math.max(.025,Number(request.inset)||.04)):null,cached=cachedPage??(compact?cachedBaseBranch(request._baseMeshFallback,'bouquet',source,Math.max(.025,Number(request.inset)||.04)):null),ba=cached?null:await anchor(bc,deadline,{compact}),ca=isBase(c)?ba:await anchor(c,deadline,{compact});if(cached&&metrics){metrics.baseCacheHits++;if(!cachedPage)metrics.staticBaseHits++;}
 const go=(cfg,at)=>{if(metrics)metrics.transportCalls++;const got=transport(cfg,at.source,source,at.tau,{deadline,metrics});const high=compact?(isBase(cfg)?checkedCompactMeshRoot(cfg,source,got.tau,{deadline}):compactMeshPlotRoot(cfg,source,got.tau,{deadline})):request.preview?got:correct(cfg,source,got.tau,{precision:80,digits:18,deadline});if(metrics){if(high.endpointChecked){metrics.endpointPrecisionChecks++;metrics.endpointMatrixEvaluations+=high.metrics.matrixEvaluations;}if(high.plotChecked)metrics.plotEndpointChecks++;if(high.plotFallback)metrics.plotFallbacks++;metrics.plotMatrixEvaluations+=high.plotAttemptMatrixEvaluations??0;}cacheAnchor(cfg,{...at,source,tau:high.tau},{compact});return branch(high.tau,{...high,provenance:at.provenance})};
 const base=cached??go(bc,ba),result=isBase(c)?base:go(c,ca);return{base,result,source:source.map(Number),elapsedMs:clock()-begun,precisionDescription:'Marked Hodge branch, fixed Schur matching refinement; numerical estimates, not interval certified.'};
}
export async function bouquetMesh(request){
 const c=validate(request.caseConfig),begun=clock(),deadline=request._meshPoolChild?Math.min(begun+traceBudget(request),begun+Math.max(0,Number(request._meshWallDeadline)-Date.now())):begun+traceBudget(request),total=c.weights.reduce((a,b)=>a+b,0),resolution=Math.max(3,Math.min(Number(request.resolution)||8,total>48?3:total>24?4:12)),inset=Math.max(.025,Number(request.inset)||.04),points=[],triangles=[],map=new Map();
 for(let i=0;i<=resolution;i++)for(let j=0;j<=resolution-i;j++){const b=[inset+(1-3*inset)*i/resolution,inset+(1-3*inset)*j/resolution,inset+(1-3*inset)*(resolution-i-j)/resolution];map.set(i+':'+j,points.length);points.push({r:.5*(b[0]+b[2]),s:.5*(b[1]+b[2])})}
 for(let i=0;i<resolution;i++)for(let j=0;j<resolution-i;j++){const a=map.get(i+':'+j),b=map.get((i+1)+':'+j),d=map.get(i+':'+(j+1));triangles.push([a,b,d]);if(j<resolution-i-1)triangles.push([b,map.get((i+1)+':'+(j+1)),d])}
 const progress=(completed,total)=>{const result={completed,total};if(typeof request.onProgress==='function')request.onProgress(result);if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result});};
 if(!request._meshPoolChild&&!baseMeshCovers(request.baseMesh,'bouquet',points,inset))request={...request,_baseMeshFallback:await verifiedBaseMesh('bouquet',resolution,inset,points)};
 if(!request._meshPoolChild&&request._meshWorkerURL){const pooled=await pooledCompactMesh({request,points,triangles,resolution,inset,deadline,workerURL:request._meshWorkerURL,onProgress:progress});if(pooled)return{...pooled,elapsedMs:clock()-begun};}
 const selected=request._meshPoolChild?request._meshIndices:points.map((_,i)=>i);if(!Array.isArray(selected)||selected.some(i=>!Number.isInteger(i)||i<0||i>=points.length)||new Set(selected).size!==selected.length)throw Error('Invalid compact three-loop source partition');
 const samples=[],pending=new Set(selected),metrics={transportCalls:0,localSteps:0,acceptedLocalSteps:0,rejectedLocalSteps:0,fixedFallbacks:0,fastMatrixEvaluations:0,endpointPrecisionChecks:0,endpointMatrixEvaluations:0,plotEndpointChecks:0,plotFallbacks:0,plotMatrixEvaluations:0,baseCacheHits:0,staticBaseHits:0};let previous=[request.R??.29,request.S??.31],error;
 const reusable=reusableMeshSamples(request,points,resolution,inset);
 for(const i of selected)if(reusable.has(i)){samples[i]=reusable.get(i);pending.delete(i);}
 metrics.reusedSamples=selected.length-pending.size;
 if(metrics.reusedSamples)progress(metrics.reusedSamples,selected.length);
 while(pending.size){let pick=-1,best=Infinity;for(const i of pending){const d=distance(previous,[points[i].r,points[i].s]);if(d<best){pick=i;best=d}}const p=points[pick];try{checkDeadline(deadline);samples[pick]=await bouquetPoint({...request,R:p.r,S:p.s,deadline,preview:true,_meshChecked:true},metrics);pending.delete(pick);previous=[p.r,p.s];if(request._meshPoolChild&&typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:{completed:selected.length-pending.size,total:selected.length,poolChunk:{sourceIndices:[pick],samples:[samples[pick]]}}});}catch(e){error=e.message;break}if((selected.length-pending.size)%4===0)progress(selected.length-pending.size,selected.length);await pause();}
 const retained=points.map((_,i)=>i).filter(i=>samples[i]),remap=new Map(retained.map((i,j)=>[i,j]));return{status:pending.size?retained.length?'partial':'unresolved':'ok',samples:retained.map(i=>samples[i]),points:retained.map(i=>points[i]),triangles:triangles.filter(t=>t.every(v=>remap.has(v))).map(t=>t.map(v=>remap.get(v))),resolution,inset,cached:false,completed:retained.length,total:selected.length,...(request._meshPoolChild?{sourceIndices:retained}:{}),metrics,elapsedMs:clock()-begun,message:error??'Computed compact image; open-edge cusp tails are resolved separately.'};
}

async function anchorAt(c,source,options){
 const at=await anchor(c,options.deadline),E=engine(150);
 if(source.every((x,i)=>new E.D(x).eq(at.source[i])))return{...at,source:source.map(String),configKey:C.configKey(c)};
 const moved=transport(c,at.source,source,at.tau,{deadline:options.deadline}),got=await C.repeatedRoot(C.cleanConfig(c),source,moved.tau,{precision:80,digits:18,deadline:options.deadline,jetScale:16});
 return{...at,source:source.map(String),tau:got.tauExact,precision:got.precisionDigits,configKey:C.configKey(c),accepted:true,checked:true,formulaInputsUsed:false};
}
export const bouquetTrace=createBouquetTrace({validate,baseOf,isBase,branch,anchorAt});
export const VERSION=C.VERSION;

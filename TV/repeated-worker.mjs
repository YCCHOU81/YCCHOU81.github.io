import {pooledCompactMesh} from './numerics/compact-mesh-pool.mjs?v=20261024-tv-audit2';
import {reusableMeshSamples} from './numerics/compact-mesh-reuse.mjs?v=20261024-tv-audit2';
import {cachedBaseBranch,verifiedBaseMesh,baseMeshCovers} from './numerics/base-mesh-cache.mjs?v=20261024-tv-audit2';
import {traceResult as tetraTrace,compactBaseConfig,compactAnchor as tetraAnchor,compactCorrect as tetraCorrect,compactRefine as tetraRefine} from './numerics/tetra-solver.mjs?v=20261024-tv-audit2';
import {boundedFastCorrect as tetraInternalCorrect} from './numerics/tetra-compact-fast.mjs?v=20261024-tv-audit2';
import {traceResult as doubleTrace,compactRefine as doubleRefine,compactCorrect as doubleCompactCorrect} from './numerics/double-solver.mjs?v=20261024-tv-audit2';
import {geometryAnchor as doubleAnchor} from './numerics/double-geometry.mjs?v=20261024-tv-audit2';
import {correct as doubleCorrect} from './numerics/double-corrector.mjs?v=20261024-tv-audit2';
import Decimal from './decimal.mjs';
import {markedHodge} from './numerics/marked-hodge.mjs';
import * as F from './numerics/generic-fast.mjs';
import {createEvaluator} from './numerics/decimal-generic.mjs?v=20261024-tv-audit2';
import {resolveConnectorInfinity} from './numerics/fourier-adaptive.mjs?v=20261024-tetra-dense4';
import {bouquetPoint,bouquetMesh,bouquetTrace} from './numerics/bouquet-solver.mjs?v=20261024-tv-audit2';
import {traceCheckpoint,traceBudget,traceStop,traceProgress,traceExtension,extendTraceGrid} from './numerics/trace-checkpoints.mjs?v=20261024-audit1';

// All caches include the full labelled graph/count key. Ordered TV weights alone
// do not distinguish geometrically different decoration branches.
const records=new Map(),engines=new Map(),clock=()=>performance.now(),pause=()=>new Promise(r=>setTimeout(r,0));
const evaluator=p=>{p=Math.min(970,Math.max(50,Math.ceil(p/10)*10));if(!engines.has(p))engines.set(p,createEvaluator(p));return engines.get(p)};
const emit=(r,node,completed,total,nodes)=>{if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:r.id,result:{node,completed,total,...(nodes?{nodes}:{})}})};
const topOf=c=>c.seed==='double'?.5:1;
const sourceGood=(c,R,S)=>Number(R)>0&&Number(R)<1&&Number(S)>0&&Number(S)<topOf(c)&&Number(R)+Number(S)/topOf(c)>1;
const distance=(a,b)=>Math.hypot(Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]));
const baseConfig=c=>c.seed==='tetra'?compactBaseConfig():({id:c.seed+'-seed',seed:c.seed,weights:c.seed==='double'?[1,1,0,0]:[1,1,1,1],shift:[0,0],connector:{q:0,family:c.seed,isogenyLift:[0,0]}});
const sameWeights=(a,b)=>a.every((v,i)=>v===b[i]);
const isBase=c=>sameWeights(c.weights,baseConfig(c).weights)&&c.shift.every(x=>x%1===0);
const opts=c=>({isogeny:!!c.connector,isogenyLift:c.connector?.isogenyLift??[0,0]});
const chartFor=(c,index)=>index===0?[1,0,0,1]:index===1?[1,-1,1,0]:c.seed==='double'?[1,0,-2,1]:[0,1,-1,1];
function validate(c){if(!c||!['double','tetra'].includes(c.seed)||!Array.isArray(c.weights)||c.weights.length!==4||c.weights.some(x=>!Number.isInteger(x)||x<0))throw Error('Invalid repeated decoration configuration');if(c.weights.reduce((a,b)=>a+b,0)>144||Math.max(...c.weights)>72)throw Error('This labelled graph exceeds the current local numerical budget (144 hemispheres / weight 72). Reduce the counts.');return c}
function characters(E,c,source){return[new E.D(String(source[1])).plus(c.shift[1]),new E.D(String(source[0])).plus(c.shift[0])]}
function decimalCorrect(c,source,guess,{precision=90,chart=[1,0,0,1],digits=30,deadline=Infinity,compact=false}={}){
 if(clock()>deadline)throw Error('The computation time budget was reached');
 if(c.seed==='tetra'){const got=tetraCorrect({config:c,source,guess,precision,chart,digits,radius:Math.min(.12,.075*(1+Math.hypot(...guess.map(Number)))),maxStep:.025,deadline,fast:false});if(!got.converged)throw Error(got.message);return got;}
 if(c.seed==='double'){const got=(compact?doubleCompactCorrect:doubleCorrect)({config:c,source,guess,precision,chart,digits,radius:Math.min(.12,.075*(1+Math.hypot(...guess.map(Number)))),maxStep:.025,deadline,fast:false});if(!got.converged)throw Error(got.message);return got;}
 const E=evaluator(precision),[r,s]=characters(E,c,source),initial=E.mobius(E.c(guess),chart),radius=E.abs(initial).plus(1).times('.22');
 const f=w=>{if(clock()>deadline)throw Error('The computation time budget was reached');return E.determinant(c.weights,E.inverse(w,chart),r,s,opts(c)).value};
 const got=E.solveFunction(f,initial,{digits:Math.min(digits,E.precision-30),maxsteps:35,guard:w=>w[1].gt('.002')&&E.abs(E.sub(w,initial)).lt(radius)}),tau=E.inverse(got.root,chart);
 return{tau:E.strs(tau),correction:E.abs(got.lastStep).toString(),iterations:got.iterations,precision:E.precision};
}
function fastCorrect(c,source,guess,{deadline=Infinity,metrics=null}={}){
 if(c.seed==='tetra'){const local={config:c,source,guess,radius:Math.min(.1,.06*(1+Math.hypot(...guess.map(Number)))),maxStep:.02,deadline};let got=tetraInternalCorrect(local);if(!got.converged&&metrics)metrics.fixedFallbacks++;if(!got.converged)got=tetraCorrect({...local,precision:80,digits:17,fast:true,reuseFastDerivative:true,verifyResidual:false});if(!got.converged)throw Error(got.message);return got;}
 const got=F.solve(c.weights,Number(source[1])+c.shift[1],Number(source[0])+c.shift[0],guess.map(Number),{...opts(c),tol:3e-12,maxsteps:32});
 if(!got.tau.every(Number.isFinite)||got.residual>2e-7||F.abs(F.sub(got.tau,guess.map(Number)))>.15*(1+F.abs(guess.map(Number))))throw Error('Fast correction needs a smaller step or more precision');
 return{tau:got.tau.map(String),correction:String(got.residual),iterations:got.iterations,precision:15};
}
function branch(tau,diagnostics={}){
 const {tau:ignoredTau,...extra}=diagnostics;
 const P=Math.max(70,diagnostics.precision??70),E=evaluator(P),t=E.c(tau),views={raw:t,zero:E.sub(1,E.div(1,t)),one:E.div(1,E.sub(1,t)),half:E.div(t,E.sub(1,E.mul(2,t)))};
 return{status:'ok',tau:{re:Number(tau[0]),im:Number(tau[1])},tauExact:tau.map(String),views:Object.fromEntries(Object.entries(views).map(([k,v])=>[k,E.strs(v)])),precisionDigits:diagnostics.precision??15,numericalErrorEstimate:diagnostics.correction??null,...extra};
}
const dbPromise=typeof indexedDB==='undefined'?Promise.resolve(null):new Promise(resolve=>{try{const r=indexedDB.open('tv-labelled-anchors-v1',1);r.onupgradeneeded=()=>r.result.createObjectStore('anchors');r.onsuccess=()=>resolve(r.result);r.onerror=()=>resolve(null)}catch{resolve(null)}});
async function readAnchor(id){const db=await dbPromise;if(!db)return null;return new Promise(resolve=>{try{const r=db.transaction('anchors').objectStore('anchors').get(id);r.onsuccess=()=>resolve(r.result??null);r.onerror=()=>resolve(null)}catch{resolve(null)}})}
async function saveAnchor(id,v){const db=await dbPromise;if(db)try{db.transaction('anchors','readwrite').objectStore('anchors').put(v,id)}catch{}}
let anchorsPromise;
function verifiedAnchors(){if(!anchorsPromise)anchorsPromise=fetch(new URL('./numerics/verified-anchors.json',import.meta.url)).then(r=>r.ok?r.json():{}).catch(()=>({}));return anchorsPromise}
async function anchor(c,deadline=Infinity){
 if(records.has(c.id))return records.get(c.id);
 if(c.seed==='tetra'){const record=await tetraAnchor(c,deadline);records.set(c.id,record);return record;}
 if(c.seed==='double'&&c.graph){const source=['.58','.31'],got=await doubleAnchor({config:c,source,deadline,correct:doubleCorrect});if(!got.checked||!got.accepted){const budget=got.stopKind==='time_budget';throw Error(budget?'The geometric anchor exhausted the computation budget':got.message);}const record={source,tau:got.tau,precision:got.precision,provenance:'independent healthy marked Hodge meshes with bounded root agreement'};records.set(c.id,record);return record;}
 const source=['.58',c.seed==='double'?'.31':'.62'];let tau,provenance,coarseTau=null;
 if(isBase(c)){tau=c.seed==='double'?['0.2952476880722449','0.3981971675549512']:['0.5904953761444898','0.7963943351099024'];provenance='calibrated marked seed'}
 else{
  const old=await readAnchor(c.id);if(old?.anchorVersion===2&&JSON.stringify(old.weights)===JSON.stringify(c.weights)){tau=old.tau;provenance='cached labelled marked graph with two-mesh branch agreement'}
  if(!tau&&c.connector){const verified=(await verifiedAnchors())[c.seed+':'+c.connector.q];if(verified){tau=verified;provenance='independently checked labelled connector anchor'}}
  if(!tau){if(!c.graph)throw Error('The labelled graph is required to select this branch');const h=markedHodge(c.graph,c.seed,...source.map(Number),{sub:8,rings:5,grades:10});coarseTau=h.tau.map(String);if(clock()>deadline)throw Error('The geometric anchor exhausted the computation budget');const fine=markedHodge(c.graph,c.seed,...source.map(Number),{sub:12,rings:7,grades:14});tau=fine.tau.map(String);provenance='two marked Hodge meshes with independent analytic branch agreement';if(clock()>deadline)throw Error('The finer geometric anchor exhausted the computation budget')}
 }
 const precision=Math.max(90,70+2*Math.max(...c.weights)),high=decimalCorrect(c,source,tau,{precision,digits:35,deadline});
 if(coarseTau){const coarse=decimalCorrect(c,source,coarseTau,{precision,digits:35,deadline}),relative=distance(coarse.tau,high.tau)/(1+Math.hypot(...high.tau.map(Number)));if(relative>1e-8)throw Error('Two geometric meshes selected different analytic branches. Reduce the counts; no image has been substituted.');}
 const record={source,tau:high.tau,precision:high.precision,provenance};records.set(c.id,record);await saveAnchor(c.id,{tau:high.tau,weights:c.weights,anchorVersion:2});return record;
}
function transport(c,from,to,guess,{precision=0,chart=[1,0,0,1],digits=30,deadline=Infinity,maxDistance=.012,depth=0,metrics=null,compact=false}={}){
 if(metrics)metrics.transportCalls++;if(clock()>deadline)throw Error('The computation time budget was reached');if(depth>22)throw Error('The selected labelled branch could not be continued safely');
 const E=precision?evaluator(precision):null,dist=distance(from,to),mid=()=>from.map((v,j)=>E?new E.D(String(v)).plus(String(to[j])).div(2).toString():String((Number(v)+Number(to[j]))/2));
 const recurse=()=>{const m=mid(),one=transport(c,from,m,guess,{precision,chart,digits,deadline,maxDistance,depth:depth+1,metrics,compact});return transport(c,m,to,one.tau,{precision,chart,digits,deadline,maxDistance,depth:depth+1,metrics,compact})};
 if(dist>maxDistance*(1+1e-12)){
  if(c.seed!=='tetra'&&!(compact&&c.seed==='double'))return recurse();
  // Equal short segments avoid rounding every requested step count up to a
  // power of two. Failed local checks still bisect their own segment.
  const steps=Math.ceil(dist/maxDistance);let at=from,got={tau:guess};
  for(let j=1;j<=steps;j++){
   const next=j===steps?to:from.map((v,k)=>E?new E.D(String(v)).plus(new E.D(String(to[k])).minus(v).times(j).div(steps)).toString():String(Number(v)+(Number(to[k])-Number(v))*j/steps));
   got=transport(c,at,next,got.tau,{precision,chart,digits,deadline,maxDistance,depth,metrics,compact});at=next;
  }
  return got;
 }
 try{
  if(metrics)metrics.localSteps++;
  const high=precision?decimalCorrect(c,to,guess,{precision,chart,digits,deadline,compact}):fastCorrect(c,to,guess,{deadline,metrics});
  const A=evaluator(Math.max(70,precision)),a=A.mobius(A.c(guess),chart),b=A.mobius(A.c(high.tau),chart),displacement=A.abs(A.sub(a,b));
  if(displacement.gt(A.abs(a).plus(1).times('.075')))throw Error('Continuation displacement is too large');
  if(c.seed==='tetra'&&dist>1e-10){if(metrics)metrics.reverseChecks++;const local={config:c,source:from,guess:high.tau,radius:Math.min(.14,.025+2*displacement.toNumber()),maxStep:.025,deadline};let reverse=precision?null:tetraInternalCorrect(local);if(!reverse?.converged&&metrics)metrics.fixedFallbacks++;if(!reverse?.converged)reverse=tetraCorrect({...local,precision:Math.max(80,precision),chart,digits:17,fast:true,reuseFastDerivative:true,verifyResidual:false});const tolerance=precision?'1e-11':A.abs(a).plus(1).times('2e-8');if(!reverse.converged||A.abs(A.sub(A.mobius(A.c(reverse.tau),chart),a)).gt(tolerance))throw Error('The reverse continuation did not recover the preceding marked root');}
  if(metrics)metrics.acceptedLocalSteps++;return high;
 }catch(error){if(metrics)metrics.rejectedLocalSteps++;if(clock()>deadline||depth>=22||(!E&&dist<1e-15))throw error;return recurse()}
}
export async function pointResult(request,metrics=null){
 if(request.caseConfig?.seed==='bouquet')return bouquetPoint(request);
 const c=validate(request.caseConfig),source=[String(request.R),String(request.S)],begun=clock();if(!sourceGood(c,...source))throw Error('Choose a point strictly inside the source triangle');const deadline=request.deadline??begun+traceBudget(request),bc=baseConfig(c),cachedPage=request.type==='mesh'?cachedBaseBranch(request.baseMesh,c.seed,source,Math.max(.025,Number(request.inset)||.04)):null,cached=cachedPage??(request.type==='mesh'?cachedBaseBranch(request._baseMeshFallback,c.seed,source,Math.max(.025,Number(request.inset)||.04)):null),b0=cached?null:await anchor(bc,deadline),a0=isBase(c)?b0:await anchor(c,deadline);if(cached&&metrics){metrics.baseCacheHits++;if(!cachedPage)metrics.staticBaseHits++;}
 const go=async(cfg,at)=>{const max=Math.max(...cfg.weights);let got;if(source.every((value,index)=>new Decimal(value).eq(at.source[index]))){got={tau:at.tau,precision:at.precision};}else try{got=transport(cfg,at.source,source,at.tau,{precision:cfg.seed==='tetra'?0:max>6?Math.max(90,70+2*max):0,deadline,maxDistance:cfg.seed==='tetra'?.025:.012,metrics,compact:request.type==='mesh'})}catch(error){got=transport(cfg,at.source,source,at.tau,{precision:Math.max(90,70+2*max),deadline,maxDistance:.006,metrics,compact:request.type==='mesh'})}const refined=cfg.seed==='tetra'?await tetraRefine(cfg,source,got.tau,deadline):cfg.seed==='double'?await doubleRefine(cfg,source,got.tau,deadline,{full:request.type!=='mesh'}):null;if(refined&&metrics)metrics.endpointPrecisionChecks++;records.set(cfg.id,{...at,source,tau:refined?.tauExact??got.tau,precision:refined?.precisionDigits??got.precision});return refined?{...refined,provenance:at.provenance}:branch(got.tau,{...got,provenance:at.provenance})};
 const base=cached??await go(bc,b0),result=isBase(c)?base:await go(c,a0);return{base,result,source:source.map(Number),elapsedMs:clock()-begun,precisionDescription:'Labelled marked anchor and bounded analytic continuation; displayed points receive decimal coordinate checks, and released points include higher-precision repeats. Numerical estimates, not interval certificates.'};
}
export async function meshResult(request){
 if(request.caseConfig?.seed==='bouquet')return bouquetMesh({...request,_meshWorkerURL:String(new URL(import.meta.url))});
 const c=validate(request.caseConfig),begun=clock(),deadline=request._meshPoolChild?Math.min(begun+traceBudget(request),begun+Math.max(0,Number(request._meshWallDeadline)-Date.now())):begun+traceBudget(request),aliases={'double-seed':'double-seed1100','double-pair':'double-pair2200','double-split':'double-split3100','double-sprout':'double-sprout2100','tetra-seed':'tetra-seed1111','tetra-pair':'tetra-pair2211'};
 if(aliases[c.id])try{const all=await fetch(new URL('./seed-cache.json',import.meta.url)).then(r=>r.json()),old=all.cases[aliases[c.id]];if(old)return{...old,status:'ok',samples:old.samples.map(p=>({base:branch(p.base),result:branch(p.result)})),cached:true,resolution:20,inset:.04,elapsedMs:clock()-begun}}catch{}
 const total=c.weights.reduce((a,b)=>a+b,0),resolution=Math.max(3,Math.min(Number(request.resolution)||8,total>48?4:total>24?6:12)),inset=Math.max(.025,Number(request.inset)||.04),top=topOf(c),points=[],triangles=[],map=new Map();
 for(let i=0;i<=resolution;i++)for(let j=0;j<=resolution-i;j++){const b=[inset+(1-3*inset)*i/resolution,inset+(1-3*inset)*j/resolution,inset+(1-3*inset)*(resolution-i-j)/resolution],r=b[0]+b[2],s=top*(b[1]+b[2]);map.set(i+':'+j,points.length);points.push({r,s})}
 for(let i=0;i<resolution;i++)for(let j=0;j<resolution-i;j++){const a=map.get(i+':'+j),b=map.get((i+1)+':'+j),d=map.get(i+':'+(j+1));triangles.push([a,b,d]);if(j<resolution-i-1)triangles.push([b,map.get((i+1)+':'+(j+1)),d])}
 if(c.seed==='double'&&!request._meshPoolChild&&!baseMeshCovers(request.baseMesh,c.seed,points,inset))request={...request,_baseMeshFallback:await verifiedBaseMesh(c.seed,resolution,inset,points)};
 if(['double','tetra'].includes(c.seed)&&!request._meshPoolChild){
  const pooled=await pooledCompactMesh({request,points,triangles,resolution,inset,deadline,workerURL:new URL(import.meta.url),onProgress:(completed,total)=>emit(request,null,completed,total)});
  if(pooled)return{...pooled,elapsedMs:clock()-begun};
 }
 const metrics={transportCalls:0,localSteps:0,acceptedLocalSteps:0,rejectedLocalSteps:0,reverseChecks:0,fixedFallbacks:0,endpointPrecisionChecks:0,...(c.seed==='double'?{baseCacheHits:0,staticBaseHits:0}:{})};
 const selectedIndices=request._meshPoolChild?request._meshIndices:points.map((_,i)=>i);if(!Array.isArray(selectedIndices)||selectedIndices.some(i=>!Number.isInteger(i)||i<0||i>=points.length)||new Set(selectedIndices).size!==selectedIndices.length)throw Error('Invalid compact source partition');const targetCount=selectedIndices.length;
 const reusable=reusableMeshSamples(request,points,resolution,inset),samples=new Array(points.length),pending=new Set(selectedIndices);let completed=0,error,previous=[request.R??.58,request.S??(c.seed==='double'?.31:.62)];
 for(const i of selectedIndices)if(reusable.has(i)){samples[i]=reusable.get(i);pending.delete(i);completed++;}
 metrics.reusedSamples=completed;
 if(completed)emit(request,null,completed,targetCount);
 // Visit neighboring source vertices, rather than jumping across the compact
 // triangle at the end of every row. The returned triangulation is unchanged.
 while(pending.size){let index=-1,best=Infinity;for(const i of pending){const p=points[i],d=distance(previous,[p.r,p.s]);if(d<best){best=d;index=i}}const p=points[index];if(clock()>deadline){error='The compact image time budget was reached';break}try{samples[index]=await pointResult({...request,R:p.r,S:p.s,deadline},metrics);completed++;pending.delete(index);previous=[p.r,p.s];if(request._meshPoolChild&&typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:{completed,total:targetCount,poolChunk:{sourceIndices:[index],samples:[samples[index]]}}})}catch(e){error=e.message;break}if(completed%4===0){emit(request,null,completed,points.length);await pause()}}
 const retained=points.map((_,i)=>i).filter(i=>samples[i]),remap=new Map(retained.map((i,j)=>[i,j]));
 return{status:completed===targetCount?'ok':completed?'partial':'unresolved',samples:retained.map(i=>samples[i]),points:retained.map(i=>points[i]),triangles:triangles.filter(t=>t.every(v=>remap.has(v))).map(t=>t.map(v=>remap.get(v))),resolution,inset,cached:false,completed,total:targetCount,...(request._meshPoolChild?{sourceIndices:retained}:{}),metrics,elapsedMs:clock()-begun,message:error??'Computed image of the compact source triangle; cusp tails are computed separately.'};
}
function makePath(c,request,P){
 const E=evaluator(P),path=request.path??{},start=(path.start??[request.R??'.58',request.S??(c.seed==='double'?'.31':'.62')]).map(String),kind=path.kind??'edge',index=Number(path.index??0),position=String(path.position??'.5'),top=new E.D(topOf(c)),p=new E.D(position);
 if(!sourceGood(c,...start))throw Error('The starting point must be strictly inside the source triangle');if(![0,1,2].includes(index))throw Error('Unknown source boundary');if(kind!=='vertex'&&(!p.gt(0)||!p.lt(1)))throw Error('Choose an open source edge position between zero and one');
 const end=kind==='vertex'?[[new E.D(1),new E.D(0)],[new E.D(0),top],[new E.D(1),top]][index]:[[new E.D(1),top.times(p)],[p,top],[p,top.times(new E.D(1).minus(p))]][index];
 return{start,kind,index,position,end:end.map(String),source:eps=>end.map((v,j)=>v.plus(new E.D(eps).times(new E.D(start[j]).minus(v))).toString())};
}
function checkedBranch(low,high,P,c,index,chart){const E=evaluator(P),error=E.abs(E.sub(E.c(high.tau),E.c(low.tau)));if(error.gt(low.precision===15?'1e-7':'1e-18'))throw Error('Independent precision corrections did not agree');const out=branch(high.tau,{precision:high.precision,correction:high.correction,crossPrecisionTauError:error.toString(),numericalErrorEstimate:error.toString(),precisionAgreementDigits:error.isZero()?Math.min(low.precision,high.precision):Math.max(0,-error.e)});const l=branch(low.tau,{precision:Math.max(P,low.precision)});out.viewErrors={};for(const k of Object.keys(out.views))out.viewErrors[k]=E.abs(E.sub(E.c(out.views[k]),E.c(l.views[k]))).toString();
 if(c?.connector&&index!==undefined){const n=c.connector.q+1,x0=index===0?new E.D(c.seed==='double'?'.25':'.5'):index===1?new E.D(c.seed==='double'?n-1:n).div(c.seed==='double'?n:n+1):c.seed==='double'?new E.D(1).div(4*n).minus('.5'):new E.D(1).div(n+1),a=E.mobius(E.c(low.tau),chart)[0].minus(x0),b=E.mobius(E.c(high.tau),chart)[0].minus(x0);out.horizontalCorrection=b.toString();if(!b.isZero()){const relative=a.minus(b).abs().div(b.abs());if(relative.gt('1e-10'))throw Error('The small horizontal cusp correction did not pass the relative precision check');out.checkedRelativeHorizontalError=relative.toString()}else out.horizontalCorrectionStatus='Zero at both working precisions; no relative accuracy claim';}
 return out}
function fourierBranch(c,source,depth,deadline){
 if(clock()>deadline)throw Error('The computation time budget was reached');
 const E=evaluator(Math.ceil(3*depth+120)),delta=new E.D(1).minus(source[0]).toString(),got=resolveConnectorInfinity({family:c.seed,q:c.connector?.q??0,delta,S:source[1],digits:12,verifyModes:true,maxPrecision:Math.min(900,Math.ceil(4*depth+250+2*Math.max(...c.weights)))});
 if(got.status!=='checked')throw Error('The Fourier cusp did not pass its precision and mode checks');
 const passed=got.checks.filter(x=>x.passed),error=passed.reduce((a,x)=>new E.D(x.absoluteWError??0).gt(a)?new E.D(x.absoluteWError??0):a,new E.D(0)).toString();const out=branch(got.w,{precision:got.precision,correction:error,engine:'deflated Fourier',precisionAgreementDigits:12,horizontalCorrection:got.correction,checkedRelativeHorizontalError:passed.find(x=>x.kind==='precision')?.relativeHorizontalError});out.viewErrors={raw:error};return out;
}
export async function traceResult(request){
 if(request.caseConfig?.seed==='bouquet')return bouquetTrace(request);
 if(request.caseConfig?.seed==='tetra'&&(request.path?.kind??'edge')==='edge'&&!(request.caseConfig.connector&&Number(request.path?.index??0)===0)){return tetraTrace({...request,onProgress:progress=>{if(typeof request.onProgress==='function')request.onProgress(progress);if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:progress});}});}
 if(request.caseConfig?.seed==='double'&&(request.path?.kind??'edge')==='edge'&&!(request.caseConfig.connector&&Number(request.path?.index??0)===0)){return doubleTrace({...request,onProgress:progress=>{if(typeof request.onProgress==='function')request.onProgress(progress);if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:progress});}});}
 const c=validate(request.caseConfig),begun=clock(),deadline=begun+traceBudget(request),target=String(request.target??'1e-10'),depth=-Math.log10(Number(target)),density=Math.max(1,Math.min(5,Number(request.samplesPerDecade)||2));if(!(depth>=1&&depth<=80))throw Error('Choose final epsilon from 0.1 through 1e-80');
 const P=Math.max(100,Math.ceil(3*depth+80+2*Math.max(...c.weights))),A=evaluator(P+40),path=makePath(c,request,P+40),chart=chartFor(c,path.index),first=Math.min(1,Math.max(Number(target),Number(request.startEpsilon)||1)),startDepth=-Math.log10(first),steps=Math.max(1,Math.ceil((depth-startDepth)*density)),bc=baseConfig(c),extension=traceExtension(request,'tv');
 let wanted=extension?extendTraceGrid(extension,target,density,d=>new A.D(10).pow(d).toString()):Array.from({length:steps+1},(_,k)=>k===steps?target:new A.D(10).pow(-startDepth-(depth-startDepth)*k/steps).toString());
 // The serialized prefix replaces anchor initialization, but its last seed is
 // independently corrected again before adding a point on the same branch.
 const {saved,resumed}=traceCheckpoint(request,'tv',()=>extension?({nodes:extension.nodes.slice(0,extension.reusedCount),wanted,oldSource:[extension.last.rExact??extension.last.r,extension.last.sExact??extension.last.s].map(String),previousEpsilon:String(extension.last.epsilon),base:extension.last.base.tauExact.slice(),selected:extension.last.result.tauExact.slice(),intermediateChecked:true,extension:{reusedCount:extension.reusedCount},extensionNeedsVerification:true}):({nodes:[],wanted,oldSource:path.start,previousEpsilon:'1',base:null,selected:null,intermediateChecked:!!c.connector||isBase(c)}));
 wanted=saved.state.wanted??wanted;
 const state=saved.state,nodes=state.nodes;let {oldSource,previousEpsilon,base,selected,intermediateChecked}=state,k=nodes.length;
 const finish=(status,message,stop=null)=>({status,nodes,completed:nodes.length,total:wanted.length,elapsedMs:clock()-begun,path:{kind:path.kind,index:path.index,position:path.position,start:path.start,endpoint:path.end},precisionDescription:'Marked graph branch; adaptive exact-decimal continuation and independent higher-precision checks. Calibrated infinity tails use cancellation-free Fourier jets with doubled modes. Numerical estimates, not interval certified.',message,...traceProgress(saved,resumed,stop)});
 try{
 if(state.extensionNeedsVerification){
  const d=-Math.log10(Number(previousEpsilon)),precision=Math.max(110,Math.ceil(3*d+100+2*Math.max(...c.weights))),E=evaluator(precision+40);
  if(path.source(previousEpsilon).some((x,i)=>new E.D(x).minus(oldSource[i]).abs().gt('1e-20')))throw Error('The saved endpoint failed source-path branch verification; its old points are retained.');
  for(const [cfg,seed] of [[bc,base],[c,selected]]){
   const checked=path.kind==='edge'&&path.index===0&&c.connector&&Number(previousEpsilon)<=.1?fourierBranch(cfg,oldSource,d,deadline).tauExact:decimalCorrect(cfg,oldSource,seed,{precision,chart,digits:Math.max(32,Math.ceil(2*d+25)),deadline}).tau;
   if(E.abs(E.sub(E.mobius(E.c(checked),chart),E.mobius(E.c(seed),chart))).gt('1e-10'))throw Error('The saved trace endpoint failed independent branch revalidation; its old points are retained.');
  }
  state.extensionNeedsVerification=false;
 }
 if(!base){const b0=await anchor(bc,deadline);base=transport(bc,b0.source,oldSource,b0.tau,{precision:100,deadline,maxDistance:.006}).tau;state.base=base;}
 if(!selected){const a0=isBase(c)?null:await anchor(c,deadline);selected=isBase(c)?base:transport(c,a0.source,oldSource,a0.tau,{precision:Math.max(100,75+2*Math.max(...c.weights)),deadline,maxDistance:.006}).tau;state.selected=selected;}
 for(;k<wanted.length;k++){
  const epsilon=wanted[k],d=-Math.log10(Number(epsilon)),source=path.source(epsilon),precision=Math.max(85,Math.ceil(3*d+65+2*Math.max(...c.weights))),digits=Math.max(28,Math.ceil(2*d+18));
  try{
   // Display sampling and continuation sampling are separate. The first two
   // decades need finer source steps for the mixed graph branches.
   const prevDepth=-Math.log10(previousEpsilon),directFourier=path.kind==='edge'&&path.index===0&&c.connector&&Number(epsilon)<=.1,internal=directFourier?1:Math.max(1,Math.ceil((d-prevDepth)*(prevDepth<2?16:5)));
   let bLow,aLow;
   for(let j=1;j<=internal;j++){
    const e=j===internal?epsilon:new A.D(10).pow(-prevDepth-(d-prevDepth)*j/internal).toString(),src=path.source(e),pd=-Math.log10(Number(e)),pr=Math.max(85,Math.ceil(3*pd+65+2*Math.max(...c.weights))),dg=Math.max(28,Math.ceil(2*pd+18));
    if(path.kind==='edge'&&path.index===0&&Number(e)<=.1&&c.connector){
     bLow={tau:base,precision:pr,correction:'0'};
     if(c.connector)aLow={tau:selected,precision:pr,correction:'0'};
     else aLow=transport(c,oldSource,src,selected,{precision:pr,chart,digits:dg,deadline,maxDistance:.006});
    }else{const fast=Number(e)>.001&&Math.max(...c.weights)<=6;bLow=transport(bc,oldSource,src,base,{precision:fast?0:pr,chart,digits:dg,deadline,maxDistance:.006});aLow=isBase(c)?bLow:transport(c,oldSource,src,selected,{precision:fast?0:pr,chart,digits:dg,deadline,maxDistance:.006})}
    base=bLow.tau;selected=aLow.tau;oldSource=src;previousEpsilon=e;Object.assign(state,{base,selected,oldSource,previousEpsilon});
    if(!intermediateChecked&&Number(e)<=.1&&c.graph){
     const h=markedHodge(c.graph,c.seed,...src.map(Number),{sub:8,rings:5,grades:10}),check=decimalCorrect(c,src,h.tau.map(String),{precision:pr,chart,digits:dg,deadline}),difference=Math.hypot(Number(check.tau[0])-Number(selected[0]),Number(check.tau[1])-Number(selected[1]));
     if(difference>.01*(1+Math.hypot(...check.tau.map(Number))))throw Error('The intermediate geometric period rejected this continuation branch; decrease the path length or counts');
     selected=check.tau;aLow=check;intermediateChecked=true;Object.assign(state,{selected,intermediateChecked});
    }
   }
   let b,a;
   if(path.kind==='edge'&&path.index===0&&Number(epsilon)<=.1&&c.connector){b=fourierBranch(bc,source,d,deadline);a=isBase(c)?b:fourierBranch(c,source,d,deadline);base=b.tauExact;selected=a.tauExact}
   else{if(bLow.precision===15)bLow=decimalCorrect(bc,source,base,{precision,chart,digits,deadline});if(!isBase(c)&&aLow.precision===15)aLow=decimalCorrect(c,source,selected,{precision,chart,digits,deadline});const highB=decimalCorrect(bc,source,base,{precision:precision+40,chart,digits,deadline}),highA=isBase(c)?highB:decimalCorrect(c,source,selected,{precision:precision+40,chart,digits,deadline});b=checkedBranch(bLow,highB,precision+40,path.kind==='edge'?bc:null,path.index,chart);a=isBase(c)?b:checkedBranch(aLow,highA,precision+40,path.kind==='edge'?c:null,path.index,chart);base=highB.tau;selected=highA.tau}
   const node={epsilon,r:source[0],s:source[1],rExact:source[0],sExact:source[1],status:'ok',base:b,result:a};nodes.push(node);previousEpsilon=epsilon;Object.assign(state,{base,selected,oldSource,previousEpsilon,intermediateChecked});emit(request,node,k+1,wanted.length,nodes);await pause();
  }catch(error){throw error}
 }
 saved.resumable=false;
 return finish('ok','Actual selected path computed locally. No cusp endpoint was inserted.');
 }catch(error){
  const stop=traceStop(saved,error,wanted[k]??target,nodes.length,wanted.length);
  return finish(nodes.length?'partial':'unresolved',stop.kind==='time_budget'?'Time budget reached. All resolved points and internal continuation steps are saved; Continue fills the remaining path.':'Stopped at the last resolved numerical point: '+error.message,stop);
 }

}
// Keep seed routing separate from numerical continuation. All three operations
// retain their existing validation, caches, precision checks and message shape.
export async function handleRequest(request) {
  const isBouquet = request.caseConfig?.seed === 'bouquet';
  switch (request.type) {
    case 'point':
      return isBouquet ? bouquetPoint(request) : pointResult(request);
    case 'mesh':
      return isBouquet
        ? bouquetMesh({...request, _meshWorkerURL:String(new URL(import.meta.url))})
        : meshResult(request);
    case 'trace':
      return isBouquet ? bouquetTrace(request) : traceResult(request);
    default:
      throw Error('Unknown numerical request');
  }
}
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
  self.onmessage = async event => {
    const request = event.data;
    try {
      self.postMessage({type:request.type, id:request.id, result:await handleRequest(request)});
    } catch (error) {
      self.postMessage({type:'error', id:request.id, message:error.message});
    }
  };
}


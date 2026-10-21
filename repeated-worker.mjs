import Decimal from './decimal.mjs';
import {markedHodge} from './numerics/marked-hodge.mjs';
import * as F from './numerics/generic-fast.mjs';
import {createEvaluator} from './numerics/decimal-generic.mjs';
import {resolveConnectorInfinity} from './numerics/fourier-adaptive.mjs';

// All caches include the full labelled graph/count key. Ordered TV weights alone
// do not distinguish geometrically different decoration branches.
const records=new Map(),engines=new Map(),clock=()=>performance.now(),pause=()=>new Promise(r=>setTimeout(r,0));
const evaluator=p=>{p=Math.min(970,Math.max(50,Math.ceil(p/10)*10));if(!engines.has(p))engines.set(p,createEvaluator(p));return engines.get(p)};
const emit=(r,node,completed,total)=>{if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:r.id,result:{node,completed,total}})};
const topOf=c=>c.seed==='double'?.5:1;
const sourceGood=(c,R,S)=>Number(R)>0&&Number(R)<1&&Number(S)>0&&Number(S)<topOf(c)&&Number(R)+Number(S)/topOf(c)>1;
const distance=(a,b)=>Math.hypot(Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]));
const baseConfig=c=>({id:c.seed+'-seed',seed:c.seed,weights:c.seed==='double'?[1,1,0,0]:[1,1,1,1],shift:[0,0],connector:{q:0,family:c.seed,isogenyLift:[0,0]}});
const sameWeights=(a,b)=>a.every((v,i)=>v===b[i]);
const isBase=c=>sameWeights(c.weights,baseConfig(c).weights)&&c.shift.every(x=>x%1===0);
const opts=c=>({isogeny:!!c.connector,isogenyLift:c.connector?.isogenyLift??[0,0]});
const chartFor=(c,index)=>index===0?[1,0,0,1]:index===1?[1,-1,1,0]:c.seed==='double'?[1,0,-2,1]:[0,1,-1,1];
function validate(c){if(!c||!['double','tetra'].includes(c.seed)||!Array.isArray(c.weights)||c.weights.length!==4||c.weights.some(x=>!Number.isInteger(x)||x<0))throw Error('Invalid repeated decoration configuration');if(c.weights.reduce((a,b)=>a+b,0)>144||Math.max(...c.weights)>72)throw Error('This labelled graph exceeds the current local numerical budget (144 hemispheres / weight 72). Reduce the counts.');return c}
function characters(E,c,source){return[new E.D(String(source[1])).plus(c.shift[1]),new E.D(String(source[0])).plus(c.shift[0])]}
function decimalCorrect(c,source,guess,{precision=90,chart=[1,0,0,1],digits=30,deadline=Infinity}={}){
 if(clock()>deadline)throw Error('The 45-second computation budget was reached');
 const E=evaluator(precision),[r,s]=characters(E,c,source),initial=E.mobius(E.c(guess),chart),radius=E.abs(initial).plus(1).times('.22');
 const f=w=>{if(clock()>deadline)throw Error('The 45-second computation budget was reached');return E.determinant(c.weights,E.inverse(w,chart),r,s,opts(c)).value};
 const got=E.solveFunction(f,initial,{digits:Math.min(digits,E.precision-30),maxsteps:35,guard:w=>w[1].gt('.002')&&E.abs(E.sub(w,initial)).lt(radius)}),tau=E.inverse(got.root,chart);
 return{tau:E.strs(tau),correction:E.abs(got.lastStep).toString(),iterations:got.iterations,precision:E.precision};
}
function fastCorrect(c,source,guess){
 const got=F.solve(c.weights,Number(source[1])+c.shift[1],Number(source[0])+c.shift[0],guess.map(Number),{...opts(c),tol:3e-12,maxsteps:32});
 if(!got.tau.every(Number.isFinite)||got.residual>2e-7||F.abs(F.sub(got.tau,guess.map(Number)))>.15*(1+F.abs(guess.map(Number))))throw Error('Fast correction needs a smaller step or more precision');
 return{tau:got.tau.map(String),correction:String(got.residual),iterations:got.iterations,precision:15};
}
function branch(tau,diagnostics={}){
 const {tau:ignoredTau,...extra}=diagnostics;
 const P=Math.max(70,diagnostics.precision??70),E=evaluator(P),t=E.c(tau),views={raw:t,zero:E.sub(1,E.div(1,t)),one:E.div(1,E.sub(1,t)),half:E.div(t,E.sub(1,E.mul(2,t)))};
 return{status:'ok',tau:{re:Number(tau[0]),im:Number(tau[1])},tauExact:tau.map(String),views:Object.fromEntries(Object.entries(views).map(([k,v])=>[k,E.strs(v)])),precisionDigits:diagnostics.precision??15,numericalErrorEstimate:diagnostics.correction??null,...extra};
}
const dbPromise=typeof indexedDB==='undefined'?Promise.resolve(null):new Promise(resolve=>{try{const r=indexedDB.open('hemisphere-labelled-anchors-v1',1);r.onupgradeneeded=()=>r.result.createObjectStore('anchors');r.onsuccess=()=>resolve(r.result);r.onerror=()=>resolve(null)}catch{resolve(null)}});
async function readAnchor(id){const db=await dbPromise;if(!db)return null;return new Promise(resolve=>{try{const r=db.transaction('anchors').objectStore('anchors').get(id);r.onsuccess=()=>resolve(r.result??null);r.onerror=()=>resolve(null)}catch{resolve(null)}})}
async function saveAnchor(id,v){const db=await dbPromise;if(db)try{db.transaction('anchors','readwrite').objectStore('anchors').put(v,id)}catch{}}
let anchorsPromise;
function verifiedAnchors(){if(!anchorsPromise)anchorsPromise=fetch(new URL('./numerics/verified-anchors.json',import.meta.url)).then(r=>r.ok?r.json():{}).catch(()=>({}));return anchorsPromise}
async function anchor(c,deadline=Infinity){
 if(records.has(c.id))return records.get(c.id);
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
function transport(c,from,to,guess,{precision=0,chart=[1,0,0,1],digits=30,deadline=Infinity,maxDistance=.012,depth=0}={}){
 if(clock()>deadline)throw Error('The 45-second computation budget was reached');if(depth>22)throw Error('The selected labelled branch could not be continued safely');
 const E=precision?evaluator(precision):null,dist=distance(from,to),mid=()=>from.map((v,j)=>E?new E.D(String(v)).plus(String(to[j])).div(2).toString():String((Number(v)+Number(to[j]))/2));
 const recurse=()=>{const m=mid(),one=transport(c,from,m,guess,{precision,chart,digits,deadline,maxDistance,depth:depth+1});return transport(c,m,to,one.tau,{precision,chart,digits,deadline,maxDistance,depth:depth+1})};
 if(dist>maxDistance)return recurse();
 try{
  const high=precision?decimalCorrect(c,to,guess,{precision,chart,digits,deadline}):fastCorrect(c,to,guess);
  const A=evaluator(Math.max(70,precision)),a=A.mobius(A.c(guess),chart),b=A.mobius(A.c(high.tau),chart),displacement=A.abs(A.sub(a,b));
  if(displacement.gt(A.abs(a).plus(1).times('.075')))throw Error('Continuation displacement is too large');
  return high;
 }catch(error){if(clock()>deadline||depth>=22||(!E&&dist<1e-15))throw error;return recurse()}
}
export async function pointResult(request){
 const c=validate(request.caseConfig),source=[String(request.R),String(request.S)];if(!sourceGood(c,...source))throw Error('Choose a point strictly inside the source triangle');const deadline=request.deadline??clock()+45000,bc=baseConfig(c),b0=await anchor(bc,deadline),a0=isBase(c)?b0:await anchor(c,deadline);
 const go=(cfg,at)=>{const max=Math.max(...cfg.weights);let got;try{got=transport(cfg,at.source,source,at.tau,{precision:max>6?Math.max(90,70+2*max):0,deadline,maxDistance:.012})}catch(error){got=transport(cfg,at.source,source,at.tau,{precision:Math.max(90,70+2*max),deadline,maxDistance:.006})}records.set(cfg.id,{...at,source,tau:got.tau,precision:got.precision});return branch(got.tau,{...got,provenance:at.provenance})};
 const base=go(bc,b0),result=isBase(c)?base:go(c,a0);return{base,result,source:source.map(Number),elapsedMs:45000-(deadline-clock()),precisionDescription:'Labelled marked Hodge anchor and analytic branch continuation; numerical, not interval certified.'};
}
export async function meshResult(request){
 const c=validate(request.caseConfig),begun=clock(),deadline=begun+45000,aliases={'double-seed':'double-seed1100','double-pair':'double-pair2200','double-split':'double-split3100','double-sprout':'double-sprout2100','tetra-seed':'tetra-seed1111','tetra-pair':'tetra-pair2211'};
 if(aliases[c.id])try{const all=await fetch(new URL('./seed-cache.json',import.meta.url)).then(r=>r.json()),old=all.cases[aliases[c.id]];if(old)return{...old,status:'ok',samples:old.samples.map(p=>({base:branch(p.base),result:branch(p.result)})),cached:true,resolution:20,inset:.04,elapsedMs:clock()-begun}}catch{}
 const total=c.weights.reduce((a,b)=>a+b,0),resolution=Math.max(3,Math.min(Number(request.resolution)||8,total>48?4:total>24?6:12)),inset=Math.max(.025,Number(request.inset)||.04),top=topOf(c),points=[],triangles=[],map=new Map();
 for(let i=0;i<=resolution;i++)for(let j=0;j<=resolution-i;j++){const b=[inset+(1-3*inset)*i/resolution,inset+(1-3*inset)*j/resolution,inset+(1-3*inset)*(resolution-i-j)/resolution],r=b[0]+b[2],s=top*(b[1]+b[2]);map.set(i+':'+j,points.length);points.push({r,s})}
 for(let i=0;i<resolution;i++)for(let j=0;j<resolution-i;j++){const a=map.get(i+':'+j),b=map.get((i+1)+':'+j),d=map.get(i+':'+(j+1));triangles.push([a,b,d]);if(j<resolution-i-1)triangles.push([b,map.get((i+1)+':'+(j+1)),d])}
 const samples=new Array(points.length),pending=new Set(points.map((_,i)=>i));let completed=0,error,previous=[request.R??.58,request.S??(c.seed==='double'?.31:.62)];
 // Visit neighboring source vertices, rather than jumping across the compact
 // triangle at the end of every row. The returned triangulation is unchanged.
 while(pending.size){let index=-1,best=Infinity;for(const i of pending){const p=points[i],d=distance(previous,[p.r,p.s]);if(d<best){best=d;index=i}}const p=points[index];if(clock()>deadline){error='The 45-second mesh budget was reached';break}try{samples[index]=await pointResult({...request,R:p.r,S:p.s,deadline});completed++;pending.delete(index);previous=[p.r,p.s]}catch(e){error=e.message;break}if(completed%4===0)await pause()}
 const retained=points.map((_,i)=>i).filter(i=>samples[i]),remap=new Map(retained.map((i,j)=>[i,j]));
 return{status:completed===points.length?'ok':completed?'partial':'unresolved',samples:retained.map(i=>samples[i]),points:retained.map(i=>points[i]),triangles:triangles.filter(t=>t.every(v=>remap.has(v))).map(t=>t.map(v=>remap.get(v))),resolution,inset,cached:false,completed,total:points.length,elapsedMs:clock()-begun,message:error??'Computed image of the compact source triangle; cusp tails are computed separately.'};
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
 if(clock()>deadline)throw Error('The 45-second computation budget was reached');
 const E=evaluator(Math.ceil(3*depth+120)),delta=new E.D(1).minus(source[0]).toString(),got=resolveConnectorInfinity({family:c.seed,q:c.connector.q,delta,S:source[1],digits:12,verifyModes:true,maxPrecision:Math.min(900,Math.ceil(4*depth+250+2*Math.max(...c.weights)))});
 if(got.status!=='checked')throw Error('The Fourier cusp did not pass its precision and mode checks');
 const passed=got.checks.filter(x=>x.passed),error=passed.reduce((a,x)=>new E.D(x.absoluteWError??0).gt(a)?new E.D(x.absoluteWError??0):a,new E.D(0)).toString();const out=branch(got.w,{precision:got.precision,correction:error,engine:'deflated Fourier',precisionAgreementDigits:12,horizontalCorrection:got.correction,checkedRelativeHorizontalError:passed.find(x=>x.kind==='precision')?.relativeHorizontalError});out.viewErrors={raw:error};return out;
}
export async function traceResult(request){
 const c=validate(request.caseConfig),begun=clock(),deadline=begun+45000,target=String(request.target??'1e-10'),depth=-Math.log10(Number(target)),density=Math.max(1,Math.min(5,Number(request.samplesPerDecade)||2));if(!(depth>=1&&depth<=80))throw Error('Choose final epsilon from 0.1 through 1e-80');
 const P=Math.max(100,Math.ceil(3*depth+80+2*Math.max(...c.weights))),A=evaluator(P+40),path=makePath(c,request,P+40),chart=chartFor(c,path.index),first=Math.min(1,Math.max(Number(target),Number(request.startEpsilon)||1)),startDepth=-Math.log10(first),steps=Math.max(1,Math.ceil((depth-startDepth)*density)),wanted=Array.from({length:steps+1},(_,k)=>k===steps?target:new A.D(10).pow(-startDepth-(depth-startDepth)*k/steps).toString()),bc=baseConfig(c),nodes=[];let intermediateChecked=!!c.connector||isBase(c);
 let b0=await anchor(bc,deadline),a0=isBase(c)?b0:await anchor(c,deadline),oldSource=path.start,previousEpsilon=1,base=transport(bc,b0.source,oldSource,b0.tau,{precision:100,deadline,maxDistance:.006}).tau,selected=isBase(c)?base:transport(c,a0.source,oldSource,a0.tau,{precision:Math.max(100,75+2*Math.max(...c.weights)),deadline,maxDistance:.006}).tau;
 for(let k=0;k<wanted.length;k++){
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
    base=bLow.tau;selected=aLow.tau;oldSource=src;
    if(!intermediateChecked&&Number(e)<=.1&&c.graph){
     const h=markedHodge(c.graph,c.seed,...src.map(Number),{sub:8,rings:5,grades:10}),check=decimalCorrect(c,src,h.tau.map(String),{precision:pr,chart,digits:dg,deadline}),difference=Math.hypot(Number(check.tau[0])-Number(selected[0]),Number(check.tau[1])-Number(selected[1]));
     if(difference>.01*(1+Math.hypot(...check.tau.map(Number))))throw Error('The intermediate geometric period rejected this continuation branch; decrease the path length or counts');
     selected=check.tau;aLow=check;intermediateChecked=true;
    }
   }
   let b,a;
   if(path.kind==='edge'&&path.index===0&&Number(epsilon)<=.1&&c.connector){b=fourierBranch(bc,source,d,deadline);a=isBase(c)?b:fourierBranch(c,source,d,deadline);base=b.tauExact;selected=a.tauExact}
   else{if(bLow.precision===15)bLow=decimalCorrect(bc,source,base,{precision,chart,digits,deadline});if(!isBase(c)&&aLow.precision===15)aLow=decimalCorrect(c,source,selected,{precision,chart,digits,deadline});const highB=decimalCorrect(bc,source,base,{precision:precision+40,chart,digits,deadline}),highA=isBase(c)?highB:decimalCorrect(c,source,selected,{precision:precision+40,chart,digits,deadline});b=checkedBranch(bLow,highB,precision+40,path.kind==='edge'?bc:null,path.index,chart);a=isBase(c)?b:checkedBranch(aLow,highA,precision+40,path.kind==='edge'?c:null,path.index,chart);base=highB.tau;selected=highA.tau}
   const node={epsilon,r:source[0],s:source[1],rExact:source[0],sExact:source[1],status:'ok',base:b,result:a};nodes.push(node);emit(request,node,k+1,wanted.length);previousEpsilon=Number(epsilon);await pause();
  }catch(error){return{status:nodes.length?'partial':'unresolved',nodes,completed:nodes.length,total:wanted.length,elapsedMs:clock()-begun,path:{kind:path.kind,index:path.index,position:path.position,start:path.start,endpoint:path.end},precisionDescription:'Adaptive exact-decimal path; higher-precision corrections, Fourier mode checks on calibrated infinity branches.',message:'Stopped at the last resolved numerical point: '+error.message}}
 }
 return{status:'ok',nodes,completed:nodes.length,total:wanted.length,elapsedMs:clock()-begun,path:{kind:path.kind,index:path.index,position:path.position,start:path.start,endpoint:path.end},precisionDescription:'Marked graph branch; adaptive exact-decimal continuation and higher-precision checks. Calibrated infinity tails use cancellation-free Fourier jets with doubled modes. Numerical estimates, not interval certified.',message:'Actual selected path computed locally. No cusp endpoint was inserted.'};
}
export async function handleRequest(request){return request.type==='point'?pointResult(request):request.type==='mesh'?meshResult(request):request.type==='trace'?traceResult(request):Promise.reject(Error('Unknown numerical request'))}
if(typeof self!=='undefined'&&typeof self.postMessage==='function')self.onmessage=async event=>{const request=event.data;try{self.postMessage({type:request.type,id:request.id,result:await handleRequest(request)})}catch(error){self.postMessage({type:'error',id:request.id,message:error.message})}};

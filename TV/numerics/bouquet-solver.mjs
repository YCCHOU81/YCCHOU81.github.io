import {markedHodge} from './marked-hodge.mjs?v=20261022-bouquet1';
import * as F from './generic-fast.mjs';
import {createEvaluator} from './decimal-generic.mjs';
import {createEvaluator as createFixedEvaluator} from './bigint-generic.mjs';
import {fixedSchur} from './fixed-schur.mjs';
import verifiedAnchors from './bouquet-verified-anchors.mjs';
import {traceCheckpoint,traceBudget,traceStop,traceProgress,traceExtension,extendTraceGrid} from './trace-checkpoints.mjs?v=20261022-depth2';

const engines=new Map(),anchors=new Map(),clock=()=>performance.now(),pause=()=>new Promise(r=>setTimeout(r,0));
const engine=p=>{p=Math.max(35,Math.min(950,Math.ceil(p/5)*5));if(!engines.has(p))engines.set(p,createEvaluator(p));return engines.get(p)};
const fixedEngines=new Map();
const stableEngine=(p,n)=>{const guard=Math.ceil(2*n*Math.log10(Math.max(2,2*n))+45),P=Math.max(80,Math.min(950,Math.ceil((p+guard)/5)*5));if(!fixedEngines.has(P))fixedEngines.set(P,createFixedEvaluator(P));return fixedEngines.get(P)};
const CH=[[1,0,0,1],[1,-1,1,0],[0,1,-1,1]],PERM=[[0,1,2,3],[0,2,3,1],[0,3,1,2]];
const inside=(R,S)=>Number(R)>0&&Number(R)<.5&&Number(S)>0&&Number(S)<.5&&Number(R)+Number(S)>.5;
const distance=(a,b)=>Math.hypot(Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]));
const baseOf=c=>({id:'bouquet-base1000',seed:'bouquet',weights:[1,0,0,0],shift:[0,0],graph:c.baseGraph??c.graph,cusps:[0,1,2].map(edge=>({edge,M:1,p:1,rational:'1/2',x0:.5,K:'2',logK:Math.log(2),signL:1}))});
const isBase=c=>c.weights.join(',')==='1,0,0,0';
const rat=(E,x)=>String(x).includes('/')?new E.D(String(x).split('/')[0]).div(String(x).split('/')[1]):new E.D(String(x));
const checkDeadline=d=>{if(clock()>d)throw Error('The 45-second numerical budget was reached; the last resolved points are retained')};
function validate(c){if(c.seed!=='bouquet'||c.weights?.length!==4||c.weights.some(n=>!Number.isInteger(n)||n<0)||c.cusps?.length!==3)throw Error('The labelled three-loop configuration is incomplete');if(c.weights.reduce((a,b)=>a+b,0)>144||Math.max(...c.weights)>72)throw Error('This graph exceeds the local numerical budget (144 hemispheres / weight 72). Reduce the counts.');return c}
function correct(c,source,guess,{precision=55,digits=24,edge=0,deadline=Infinity,radius=.24}={}){
 checkDeadline(deadline);const E=stableEngine(precision,Math.max(...c.weights)),chart=CH[edge],initial=E.mobius(E.c(guess),chart),[a,b,cc,d]=chart;
 const r=new E.D(source[1]).plus(c.shift[1]),s=new E.D(source[0]).plus(c.shift[0]),rr=r.times(a).minus(s.times(b)),ss=s.times(d).minus(r.times(cc)),weights=PERM[edge].map(i=>c.weights[i]);
 const residual=fixedSchur(E,w=>{checkDeadline(deadline);return E.tvMatrix(weights,w,rr,ss,{reduce:false}).M},initial);
 const got=E.solveFunction(residual.evaluate,initial,{digits:Math.min(Math.ceil(digits),E.precision-18),maxsteps:28,guard:w=>w[1].gt('.00001')&&E.abs(E.sub(w,initial)).lt(E.abs(initial).plus(1).times(radius))});
 return{tau:E.strs(E.inverse(got.root,chart)),w:E.strs(got.root),precision:E.precision,correction:E.abs(got.lastStep).toString(),iterations:got.iterations};
}
function fastCorrect(c,source,guess,{deadline=Infinity}={}){
 checkDeadline(deadline);guess=guess.map(Number);const residual=fixedSchur(F,w=>F.tvMatrix(c.weights,w,Number(source[1])+c.shift[1],Number(source[0])+c.shift[0],{reduce:false}).M,guess);
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
 const {tau:ignoredTau,...extra}=diagnostics;
 const E=engine(Math.max(65,diagnostics.precision??65)),t=E.c(tau),views={raw:t,zero:E.sub(1,E.div(1,t)),one:E.div(1,E.sub(1,t)),half:E.div(t,E.sub(1,E.mul(2,t)))};
 return{status:'ok',tau:{re:Number(tau[0]),im:Number(tau[1])},tauExact:tau.map(String),views:Object.fromEntries(Object.entries(views).map(([k,v])=>[k,E.strs(v)])),precisionDigits:diagnostics.precision??15,numericalErrorEstimate:diagnostics.correction??null,...extra};
}
async function anchor(c,deadline){
 if(anchors.has(c.id))return anchors.get(c.id);checkDeadline(deadline);if(!c.graph)throw Error('The marked graph is required to identify this branch');
 const counts=c.graph.stable_edges?.map(e=>e.count??e.q??0)??[],old=counts.every(q=>q===0)?verifiedAnchors.cores[c.core??'bouquet_B0_T0']:verifiedAnchors.counts[(c.core??c.graph.core_id)+':'+counts.join(',')];
 if(old&&old.weights.join(',')===c.weights.join(',')){const got=correct(c,old.source.map(String),old.tau,{precision:55,digits:25,deadline}),record={source:old.source.map(String),tau:got.tau,precision:got.precision,provenance:old.provenance};anchors.set(c.id,record);return record;}
 const source=['.29','.31'],coarse=markedHodge(c.graph,'bouquet',.29,.31,{sub:8,rings:5,grades:10}),fine=markedHodge(c.graph,'bouquet',.29,.31,{sub:12,rings:7,grades:14});
 const a=correct(c,source,coarse.tau,{precision:35,digits:20,deadline,radius:.5}),b=correct(c,source,fine.tau,{precision:45,digits:24,deadline,radius:.5});
 if(distance(a.tau,b.tau)>1e-11)throw Error('Two independent marked Hodge meshes selected different roots; no branch has been substituted');
 const record={source,tau:b.tau,precision:b.precision,provenance:'Two marked graph meshes refined by the fixed Schur matching equation'};anchors.set(c.id,record);return record;
}
function transport(c,from,to,guess,{edge=0,precision=0,digits=25,deadline=Infinity,depth=0,maxDistance=.012}={}){
 checkDeadline(deadline);if(depth>20)throw Error('The branch needs a shorter continuation step');const E=engine(Math.max(65,precision)),dist=distance(from,to);
 const split=()=>{const m=from.map((x,i)=>new E.D(x).plus(to[i]).div(2).toString()),a=transport(c,from,m,guess,{edge,precision,digits,deadline,depth:depth+1,maxDistance});return transport(c,m,to,a.tau,{edge,precision,digits,deadline,depth:depth+1,maxDistance})};
 if(dist>maxDistance)return split();
 try{let out;if(precision)out=correct(c,to,guess,{edge,precision,digits,deadline});else try{out=fastCorrect(c,to,guess,{deadline})}catch{out=correct(c,to,guess,{edge,precision:45,digits:20,deadline})}const a=E.mobius(E.c(guess),CH[edge]),b=E.mobius(E.c(out.tau),CH[edge]);if(E.abs(E.sub(a,b)).gt(E.abs(a).plus(1).times('.13')))throw Error('The branch needs a smaller step');return out}
 catch(error){if(depth>=20||clock()>deadline||dist<1e-16)throw error;return split()}
}
function predicted(c,source,edge,E){const data=c.cusps[edge],delta=[new E.D('.5').minus(source[0]),new E.D('.5').minus(source[1]),new E.D(source[0]).plus(source[1]).minus('.5')][edge],t=new E.D(source[edge===0?1:0]),K=data.KLabel?rat(E,data.KLabel):data.K?rat(E,data.K):new E.D(data.logK).exp();return[rat(E,data.rational??`${data.p}/${2*data.M}`),K.times(E.TPI.times(t).sin()).div(delta).ln().div(E.PI.times(data.M))]}
export function bouquetCuspCheck(c,source,edge,depth,deadline=Infinity,previousTau=null){
 const M=c.cusps[edge].M,needed=Math.ceil(2*depth/M+18),P=Math.max(70,Math.ceil(depth+35+1.1*Math.max(...c.weights)),Math.ceil(3*depth/M+45)),E=engine(P),initial=predicted(c,source,edge,E);
 if(initial[1].lte('1.15'))return null;
 // A fixed Schur pivot can be very steep even when two roots are close in
 // the chart. Reuse only a virtually identical initializer, otherwise the
 // proved phase/constant gives the safer local starting point.
 const reuse=previousTau&&E.abs(E.sub(E.mobius(E.c(previousTau),CH[edge]),initial)).lt('.0001'),initialTau=reuse?previousTau:E.strs(E.inverse(initial,CH[edge]));
 const low=correct(c,source,initialTau,{edge,precision:P,digits:needed,deadline,radius:.15}),high=correct(c,source,low.tau,{edge,precision:P+35,digits:needed+3,deadline,radius:.08});
 const H=engine(P+35),error=H.abs(H.sub(H.c(high.w),H.c(low.w))),x0=rat(H,c.cusps[edge].rational??`${c.cusps[edge].p}/${2*c.cusps[edge].M}`),dx=new H.D(high.w[0]).minus(x0),dxlow=new H.D(low.w[0]).minus(x0);
 if(error.gt('1e-12'))throw Error('The cusp corrections did not agree at independent precisions');
 if(dx.abs().gte(new H.D(1).div(4*c.cusps[edge].M)))throw Error('The cusp root is outside the marked sector; continue farther before using the asymptotic initializer');
 const out=branch(high.tau,{precision:high.precision,correction:error.toString(),crossPrecisionTauError:H.abs(H.sub(H.c(high.tau),H.c(low.tau))).toString(),horizontalCorrection:dx.toString(),engine:'fixed Schur elliptic matching; proved marked cusp initializer'});
 const other=branch(low.tau,{precision:P});out.viewErrors=Object.fromEntries(Object.keys(out.views).map(k=>[k,H.abs(H.sub(H.c(out.views[k]),H.c(other.views[k]))).toString()]));
 if(!dx.isZero()){const rel=dx.minus(dxlow).abs().div(dx.abs());if(!rel.lt('1e-8'))throw Error('The tiny horizontal cusp displacement needs more precision');out.checkedRelativeHorizontalError=rel.toString();}
 else out.horizontalCorrectionStatus='Zero at both working precisions; no nonzero relative displacement is claimed';
 return out;
}
export async function bouquetPoint(request){
 const c=validate(request.caseConfig),source=[String(request.R),String(request.S)],begun=clock(),deadline=request.deadline??begun+45000;if(!inside(...source))throw Error('Choose a point strictly inside 0 < R,S < 1/2, R+S > 1/2');const bc=baseOf(c),ba=await anchor(bc,deadline),ca=isBase(c)?ba:await anchor(c,deadline);
 const go=(cfg,at)=>{const got=transport(cfg,at.source,source,at.tau,{deadline});const high=request.preview?got:correct(cfg,source,got.tau,{precision:55,digits:25,deadline});anchors.set(cfg.id,{...at,source,tau:high.tau});return branch(high.tau,{...high,provenance:at.provenance})};
 const base=go(bc,ba),result=isBase(c)?base:go(c,ca);return{base,result,source:source.map(Number),elapsedMs:clock()-begun,precisionDescription:'Marked Hodge branch, fixed Schur matching refinement; numerical estimates, not interval certified.'};
}
export async function bouquetMesh(request){
 const c=validate(request.caseConfig),begun=clock(),deadline=begun+traceBudget(request),total=c.weights.reduce((a,b)=>a+b,0),resolution=Math.max(3,Math.min(Number(request.resolution)||8,total>48?3:total>24?4:12)),inset=Math.max(.025,Number(request.inset)||.04),points=[],triangles=[],map=new Map();
 for(let i=0;i<=resolution;i++)for(let j=0;j<=resolution-i;j++){const b=[inset+(1-3*inset)*i/resolution,inset+(1-3*inset)*j/resolution,inset+(1-3*inset)*(resolution-i-j)/resolution];map.set(i+':'+j,points.length);points.push({r:.5*(b[0]+b[2]),s:.5*(b[1]+b[2])})}
 for(let i=0;i<resolution;i++)for(let j=0;j<resolution-i;j++){const a=map.get(i+':'+j),b=map.get((i+1)+':'+j),d=map.get(i+':'+(j+1));triangles.push([a,b,d]);if(j<resolution-i-1)triangles.push([b,map.get((i+1)+':'+(j+1)),d])}
 const samples=[],pending=new Set(points.map((_,i)=>i));let previous=[request.R??.29,request.S??.31],error;while(pending.size){let pick=-1,best=Infinity;for(const i of pending){const d=distance(previous,[points[i].r,points[i].s]);if(d<best){pick=i;best=d}}const p=points[pick];try{checkDeadline(deadline);samples[pick]=await bouquetPoint({...request,R:p.r,S:p.s,deadline,preview:true});pending.delete(pick);previous=[p.r,p.s]}catch(e){error=e.message;break}await pause()}
 const retained=points.map((_,i)=>i).filter(i=>samples[i]),remap=new Map(retained.map((i,j)=>[i,j]));return{status:pending.size?retained.length?'partial':'unresolved':'ok',samples:retained.map(i=>samples[i]),points:retained.map(i=>points[i]),triangles:triangles.filter(t=>t.every(v=>remap.has(v))).map(t=>t.map(v=>remap.get(v))),resolution,inset,cached:false,completed:retained.length,total:points.length,elapsedMs:clock()-begun,message:error??'Computed compact image; open-edge cusp tails are resolved separately.'};
}
export async function bouquetTrace(request){
 const c=validate(request.caseConfig),begun=clock(),deadline=begun+traceBudget(request),target=String(request.target??'1e-10'),depth=-Math.log10(Number(target)),density=Math.max(1,Math.min(5,Number(request.samplesPerDecade)||2));if(depth<1||depth>80)throw Error('Choose final epsilon from 0.1 through 1e-80');const E=engine(3*depth+90),path=request.path??{},kind=path.kind??'edge',edge=Number(path.index??0),position=new E.D(path.position??'.5'),start=(path.start??['.29','.31']).map(String);if(!inside(...start))throw Error('Choose an interior starting point');
 if(kind!=='edge')throw Error('Choose an open source edge for the proved three-loop cusp law');if(!position.gt(0)||!position.lt(1))throw Error('Choose a position strictly between zero and one');
 const end=[[new E.D('.5'),position.div(2)],[position.div(2),new E.D('.5')],[position.div(2),new E.D('.5').minus(position.div(2))]][edge],source=e=>end.map((v,j)=>v.plus(new E.D(e).times(new E.D(start[j]).minus(v))).toString());
 const first=Math.min(1,Math.max(Number(target),Number(request.startEpsilon)||1)),d0=-Math.log10(first),steps=Math.max(1,Math.ceil((depth-d0)*density)),bc=baseOf(c),extension=traceExtension(request,'bouquet');
 let wanted=extension?extendTraceGrid(extension,target,density,d=>new E.D(10).pow(d).toString()):Array.from({length:steps+1},(_,k)=>k===steps?target:new E.D(10).pow(-d0-(depth-d0)*k/steps).toString());
 const unresolved=e=>{const src=source(e);return{epsilon:e,r:src[0],s:src[1],rExact:src[0],sExact:src[1],status:'unresolved',base:{status:'unresolved'},result:{status:'unresolved'}}};
 const {saved,resumed}=traceCheckpoint(request,'bouquet',()=>extension?({nodes:[...extension.nodes,...wanted.slice(extension.nodes.length).map(unresolved)],wanted,oldSource:[extension.last.rExact??extension.last.r,extension.last.sExact??extension.last.s].map(String),previousDepth:-Math.log10(Number(extension.last.epsilon)),baseAt:extension.last.base.tauExact.slice(),chosenAt:extension.last.result.tauExact.slice(),overlapChecked:extension.result.tailBranchValidationStatus==='matched',tailFirst:false,tailPassComplete:true,nextIndex:extension.reusedCount,recoveries:JSON.parse(JSON.stringify(extension.result.recoveries??[])),extension:{reusedCount:extension.reusedCount},extensionNeedsVerification:true}):({nodes:wanted.map(unresolved),wanted,oldSource:start,previousDepth:0,baseAt:null,chosenAt:null,overlapChecked:false,tailFirst:false,tailPassComplete:false,nextIndex:0,recoveries:[]}));
 wanted=saved.state.wanted??wanted;
 const state=saved.state,nodes=state.nodes;let {oldSource,previousDepth,baseAt,chosenAt,overlapChecked,tailFirst}=state,completed=nodes.filter(n=>n.status==='ok').length,failedEpsilon=target;
 const checkpoint=()=>Object.assign(state,{oldSource,previousDepth,baseAt,chosenAt,overlapChecked,tailFirst});
 const emit=(k)=>{if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:{node:nodes[k],nodes,completed,total:wanted.length}})};
 const keep=(k,b,a,connection)=>{nodes[k]={...nodes[k],status:'ok',base:b,result:a,branchConnection:connection};completed++;emit(k)};
 const hasCuspRoots=()=>nodes.some(n=>n.status==='ok'&&n.result?.engine?.includes('cusp'));
 const finish=(message,stop=null)=>({status:completed===wanted.length&&(!hasCuspRoots()||overlapChecked)?'ok':completed?'partial':'unresolved',nodes,completed,total:wanted.length,elapsedMs:clock()-begun,path:{kind,index:edge,position:position.toString(),start,endpoint:end.map(String)},message,recoveries:state.recoveries,...traceProgress(saved,resumed,stop),tailBranchValidationStatus:overlapChecked?'matched':hasCuspRoots()?'pending':'continued',tailBranchValidation:overlapChecked?'Independent cusp roots matched the continued marked graph branch':hasCuspRoots()?'Locally refined marked-phase tail; geometric continuation overlap is still pending':'Displayed points were continued from the marked graph anchor; no independent cusp shortcut was used',precisionDescription:'Fixed-point theta jets and a fixed Schur residual; independent precision and relative horizontal checks in the cusp. Numerical estimates, not interval certified.'});
 try{
  if(state.extensionNeedsVerification){
   const precision=Math.max(90,Math.ceil(previousDepth+65+Math.max(...c.weights))),V=engine(precision+35);
   if(source(extension?.last.epsilon??nodes[state.nextIndex-1]?.epsilon).some((x,i)=>new V.D(x).minus(oldSource[i]).abs().gt('1e-20')))throw Error('The saved endpoint failed source-path branch verification; its old points are retained.');
   for(const [cfg,seed] of [[bc,baseAt],[c,chosenAt]]){
    const check=correct(cfg,oldSource,seed,{edge,precision,digits:Math.ceil(2*previousDepth+25),deadline});
    if(V.abs(V.sub(V.mobius(V.c(check.tau),CH[edge]),V.mobius(V.c(seed),CH[edge]))).gt('1e-10'))throw Error('The saved three-loop endpoint failed independent branch revalidation; its old points are retained.');
   }
   state.extensionNeedsVerification=false;
  }
  // Larger graphs resolve the requested deep point first. Any unfinished
  // middle slots remain explicit gaps, rather than invented connecting lines.
  const finalSource=source(target),probeFinal=predicted(c,finalSource,edge,E);
  if(!state.tailPassComplete&&c.weights.reduce((a,b)=>a+b,0)>16&&depth>5&&probeFinal[1].gt('1.2')){
   let precedingBase=null,precedingSelected=null;
   for(let k=wanted.length-1;k>=Math.max(0,wanted.length-3);k--){failedEpsilon=wanted[k];if(nodes[k].status==='ok'){precedingBase=nodes[k].base.tauExact;precedingSelected=nodes[k].result.tauExact;continue;}const src=source(wanted[k]),dd=-Math.log10(Number(wanted[k]));if(predicted(c,src,edge,E)[1].lte('1.2'))break;
    const b=bouquetCuspCheck(bc,src,edge,dd,deadline,precedingBase),a=isBase(c)?b:bouquetCuspCheck(c,src,edge,dd,deadline,precedingSelected);
    if(!b||!a)break;keep(k,b,a,'Locally refined marked-phase tail; continuation overlap pending');precedingBase=b.tauExact;precedingSelected=a.tauExact;tailFirst=true;checkpoint();await pause();
   }
  }
  state.tailPassComplete=true;
  if(!baseAt){const ba=await anchor(bc,deadline);baseAt=transport(bc,ba.source,start,ba.tau,{deadline}).tau;checkpoint();}
  if(!chosenAt){const ca=isBase(c)?null:await anchor(c,deadline);chosenAt=isBase(c)?baseAt:transport(c,ca.source,start,ca.tau,{deadline}).tau;checkpoint();}
  for(let k=state.nextIndex;k<wanted.length;k++){
   const epsilon=wanted[k];failedEpsilon=epsilon;const d=-Math.log10(Number(epsilon)),src=source(epsilon);let b,a;
   const probe=predicted(c,src,edge,E);let useCusp=Number(epsilon)<.025&&probe[1].gt('1.2'),postponedCusp=false;
   if(useCusp)try{
    const cached=tailFirst&&nodes[k].status==='ok'?nodes[k]:null;b=cached?.base??bouquetCuspCheck(bc,src,edge,d,deadline);a=cached?.result??(isBase(c)?b:bouquetCuspCheck(c,src,edge,d,deadline));if(!b||!a)throw Error('The terminal cusp initializer is outside its stable chart');
    if(!overlapChecked){const pr=Math.max(65,Math.ceil(d+55+Math.max(...c.weights))),bo=transport(bc,oldSource,src,baseAt,{edge,precision:pr,digits:25,deadline}),ao=isBase(c)?bo:transport(c,oldSource,src,chosenAt,{edge,precision:pr,digits:25,deadline});const difference=E.abs(E.sub(E.mobius(E.c(ao.tau),CH[edge]),E.c(a.views[['raw','zero','one'][edge]]))),baseDifference=E.abs(E.sub(E.mobius(E.c(bo.tau),CH[edge]),E.c(b.views[['raw','zero','one'][edge]])));if(difference.gt('1e-9')||baseDifference.gt('1e-9'))throw Error('The cusp root did not match the independently continued marked graph branch');overlapChecked=true;checkpoint();nodes.forEach((node,i)=>{if(i>=(state.extension?.reusedCount??0)&&node.status==='ok')node.branchConnection='Independent cusp roots matched the continued marked graph branch'});}
    baseAt=b.tauExact;chosenAt=a.tauExact;oldSource=src;previousDepth=d;checkpoint();
   }catch(error){
    checkDeadline(deadline);
    // An asymptotic initializer can become useful later than its height
    // suggests. Its failure must not discard a safely continued graph root.
    // Continue that root instead; the cusp shortcut is retried farther out.
    useCusp=false;postponedCusp=true;
    if(state.recoveries.length<12)state.recoveries.push({epsilon,action:'Continued the marked graph branch; cusp shortcut postponed',reason:error.message});
   }
   if(!useCusp){const fromDepth=previousDepth,n=Math.max(1,Math.ceil((d-fromDepth)*(fromDepth<2?8:3)));let bo,ao;for(let j=1;j<=n;j++){const dd=fromDepth+(d-fromDepth)*j/n,ss=j===n?src:source(new E.D(10).pow(-dd).toString()),pr=Math.max(55,Math.ceil(dd+40+Math.max(...c.weights))),prec=dd>3?pr:0;bo=transport(bc,oldSource,ss,baseAt,{edge,precision:prec,digits:Math.ceil(2*dd+18),deadline});ao=isBase(c)?bo:transport(c,oldSource,ss,chosenAt,{edge,precision:prec,digits:Math.ceil(2*dd+18),deadline});baseAt=bo.tau;chosenAt=ao.tau;oldSource=ss;previousDepth=dd;checkpoint();}const bh=correct(bc,src,baseAt,{edge,precision:Math.max(65,3*d+55),digits:2*d+22,deadline}),ah=isBase(c)?bh:correct(c,src,chosenAt,{edge,precision:Math.max(65,d+55+Math.max(...c.weights)),digits:Math.ceil(2*d/c.cusps[edge].M+22),deadline});if(postponedCusp){const verify=correct(c,src,ah.tau,{edge,precision:Math.max(90,d+80+Math.max(...c.weights)),digits:Math.ceil(2*d/c.cusps[edge].M+25),deadline}),V=engine(Math.max(100,verify.precision)),error=V.abs(V.sub(V.c(verify.w),V.c(ah.w)));if(error.gt('1e-12'))throw Error('Independent corrections rejected the continued branch after a cusp initializer failure');const tauError=V.abs(V.sub(V.c(verify.tau),V.c(ah.tau))).toString();a=branch(verify.tau,{...verify,crossPrecisionTauError:tauError,numericalErrorEstimate:tauError,cuspChartError:error.toString(),recovery:'Asymptotic shortcut postponed; graph continuation independently refined'});const lower=branch(ah.tau,ah);a.viewErrors=Object.fromEntries(Object.keys(a.views).map(key=>[key,V.abs(V.sub(V.c(a.views[key]),V.c(lower.views[key]))).toString()]));chosenAt=verify.tau;}else{a=isBase(c)?null:branch(ah.tau,ah);chosenAt=ah.tau;}b=branch(bh.tau,bh);if(isBase(c)){if(postponedCusp)b=a;else a=b;baseAt=a.tauExact;chosenAt=a.tauExact;}else baseAt=bh.tau;}
   if(nodes[k].status==='ok')completed--;keep(k,b,a,overlapChecked?'Independent cusp roots matched the continued marked graph branch':postponedCusp?'Continued from the marked graph anchor; cusp shortcut postponed':'Continued from the marked graph anchor');previousDepth=d;state.nextIndex=k+1;checkpoint();await pause();
  }
 }catch(error){checkpoint();const stop=traceStop(saved,error,failedEpsilon,completed,wanted.length);return finish(stop.kind==='time_budget'?'Time budget reached. Checked tail points and the continued graph branch are saved; Continue fills unresolved slots.':'Stopped at resolved numerical points, with explicit gaps: '+error.message,stop)}
 checkpoint();saved.resumable=false;
 return finish('Actual numerical points on the marked branch; the exact cusp phase and leading constant initialize independent precision corrections. No cusp endpoint was inserted.');
}



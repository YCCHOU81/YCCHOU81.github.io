// Tetrahedral continuation with independent marked-graph branch checks.
// All period values are corrected roots, selected from independent graph data.
import {correct as decimalCorrect,getEvaluator} from './double-corrector.mjs?v=20261024-tv-audit2';
import {correct as genericCorrect,getGuardedEvaluator,fourierDiagnosticsSince} from './tetra-fixed-corrector.mjs?v=20261024-tv-audit2';
import * as Reduction from './tetra-lame-reduction.mjs?v=20261024-tv-audit2';
import {checkCanonicalTetra} from './tetra-canonical.mjs?v=20261024-tetra1';
const correct=options=>Reduction.eligible(options.config)?Reduction.correct(options):genericCorrect(options);
const compactEndpointCorrect=options=>genericCorrect({...options,compactFourier:true});
import {geometryAnchor,clearGeometryCache} from './tetra-geometry.mjs?v=20261024-tv-audit2';
import baseGraph from './tetra-base-graph.mjs?v=20261024-tetra1';
import {traceExtension,extendTraceGrid} from './tetra-checkpoints.mjs?v=20261024-audit1';

const VERSION='tetra-stable-1';
const anchors=new Map(),jobs=new Map();let serial=0;
const now=()=>performance.now();
const charts=[[1,0,0,1],[1,-1,1,0],[0,1,-1,1]];
const canonical=x=>Array.isArray(x)?x.map(canonical):x&&typeof x==='object'?Object.fromEntries(Object.keys(x).sort().map(k=>[k,canonical(x[k])])):x;
const key=x=>JSON.stringify(canonical(x));
const clone=x=>JSON.parse(JSON.stringify(x));
const caseKey=c=>key({seed:c.seed,weights:c.weights,shift:c.shift,counts:c.counts,graph:c.graph});
const countKeys=['qAB','qAC','qAD','qBC','qBD','qCD'];
const baseConfig=()=>({id:'tetra-seed',seed:'tetra',operation:'repeated',weights:[1,1,1,1],shift:[0,0],top:1,counts:Object.fromEntries(countKeys.map(k=>[k,0])),core:'tetra',graph:baseGraph});
const isBase=c=>c.weights.join(',')==='1,1,1,1'&&c.shift.every(x=>x%1===0);
function validate(c){if(c?.seed!=='tetra'||c.weights?.length!==4||c.weights.some(n=>!Number.isInteger(n)||n<1)||!c.graph)throw Error('Choose a valid labelled tetrahedral graph.');if(c.weights.reduce((a,b)=>a+b,0)>144||Math.max(...c.weights)>72)throw Error('The graph exceeds the existing numerical weight budget.');if(countKeys.some(k=>!Number.isInteger(c.counts?.[k])||c.counts[k]<0))throw Error('The six labelled nonnegative integer decoration counts are required.');const q=c.counts,n=[1+q.qAB+q.qAC+q.qBC,1+q.qAB+q.qAD+q.qBD,1+q.qAC+q.qAD+q.qCD,1+q.qBC+q.qBD+q.qCD],shift=[(q.qAC+q.qBD+q.qAD+q.qBC)%2/2,(q.qAB+q.qCD+q.qAD+q.qBC)%2/2];if(n.some((v,i)=>v!==c.weights[i])||shift.some((v,i)=>v!==c.shift?.[i]))throw Error('The labelled graph counts, weights and character shifts do not match.');const graph=checkCanonicalTetra(c);if(!graph.valid)throw Error('The full labelled ribbon graph does not match its counts: '+graph.reason);return c;}
function ensure(deadline,signal){if(signal?.aborted)throw Error('Calculation cancelled');if(now()>deadline)throw Error('Computation time budget reached');}
function localPrecision(c,depth,extra=0){return Math.min(850,Math.max(80,Math.ceil(2*depth+60+3*Math.max(...c.weights)+extra)));}
function transportPrecision(c,depth){return Math.min(850,Math.max(70,Math.ceil(2*depth+40+3*Math.max(...c.weights))));}
function chartDistance(a,b,chart,precision=120){const E=getEvaluator(precision);return E.abs(E.sub(E.mobius(E.c(a),chart),E.mobius(E.c(b),chart))).toNumber();}
function outputBranch(tau,diagnostics={}){
 const E=getEvaluator(Math.max(80,diagnostics.precision??80)),z=E.c(tau),views={raw:z,zero:E.sub(1,E.div(1,z)),one:E.div(1,E.sub(1,z)),half:E.div(z,E.sub(1,E.mul(2,z)))};
 return {status:'ok',tau:{re:Number(tau[0]),im:Number(tau[1])},tauExact:tau.map(String),views:Object.fromEntries(Object.entries(views).map(([k,v])=>[k,E.strs(v)])),precisionDigits:E.precision,numericalErrorEstimate:diagnostics.correction??null,...diagnostics};
}
async function independentAnchor(config,source,chart,deadline,metrics,signal){
 const id=caseKey(config)+'|'+key(source);
 if(anchors.has(id)){metrics.anchorCacheHits++;return{...anchors.get(id),chart,cacheHit:true};}
 ensure(deadline,signal);
 let selected;
 if(Reduction.eligible(config)){
  const got=await Reduction.geometryAnchor({config,source,chart,deadline,signal});
  metrics.reductionAnchors=(metrics.reductionAnchors??0)+1;
  if(!got.checked||!got.accepted||!got.tau){const e=Error(got.message??'The marked Lamé factor anchor is unresolved.');e.geometry=got;e.failureKind=got.failureKind??'geometry_inconclusive';throw e;}
  selected=got;
 }
 if(!selected&&isBase(config)){
  const guess=['0.5904953761444898','0.7963943351099024'];
  const atDefault=source.every((v,j)=>Math.abs(Number(v)-[.58,.62][j])<1e-14);
  if(atDefault){const a=correct({config,source,guess,chart:[1,0,0,1],precision:80,digits:30,radius:.08,maxStep:.02,deadline});if(!a.converged)throw Error('The calibrated base root could not be refined.');selected={status:'verified',tau:a.tau,precision:a.precision,diagnostics:{kind:'calibrated marked seed'}};}
 }
 if(!selected){
  const got=await geometryAnchor({config,source,chart,deadline,correct:decimalCorrect,signal});
  metrics.geometryCalls++;metrics.geometryMs+=got.elapsedMs??got.milliseconds??0;
  if(!got.tau||!['ok','verified','accepted'].includes(got.status)){
   const budgetFailure=now()>deadline||got.stopKind==='time_budget'||got.attempts?.some(a=>a.status==='time-budget');
   const e=Error(budgetFailure?'Computation time budget reached during independent geometric branch selection':got.message??'The independent geometric branch selector is inconclusive.');e.geometry=got;e.failureKind=budgetFailure?'time_budget':'geometry_inconclusive';throw e;
  }
  selected=got;
 }
 anchors.set(id,selected);while(anchors.size>32)anchors.delete(anchors.keys().next().value);return selected;
}
function addMetrics(metrics,a){metrics.correctorCalls++;metrics.evaluations+=a.evaluations??0;metrics.correctorMs+=a.elapsedMs??0;}
// Compact endpoint repeats may use the same guarded Fourier arithmetic as the
// mixed matrices. The independent marked anchors and all cusp continuations
// retain their original Decimal reduction. A failed fixed-point attempt is
// retried with Decimal in the identical bounded neighbourhood.
function compactReducedCorrect(options){
 const begun=now(),precision=Math.max(Number(options.precision)||80,(Number(options.digits)||25)+30);
 let fixed;
 try{
  const arithmetic=getGuardedEvaluator(precision,{compactFourier:true}),before=arithmetic.getFourierDiagnostics();
  fixed=Reduction.correct({...options,arithmetic});
  fixed.diagnostics={...fixed.diagnostics,arithmetic:'guarded fixed-point Fourier jets; compact marked Lamé factor',fourier:fourierDiagnosticsSince(arithmetic,before)};
  if(fixed.converged)return fixed;
 }catch(error){fixed={message:String(error?.message??error),evaluations:0};}
 if(now()>(options.deadline??Infinity)||fixed.failureKind==='time_budget')return fixed;
 const fallback=Reduction.correct(options);
 fallback.elapsedMs=now()-begun;
 fallback.evaluations+=(fixed.evaluations??0);
 fallback.diagnostics={...fallback.diagnostics,arithmetic:'Decimal fallback after compact marked Lamé factor check',
  fixedPointFallback:{message:fixed.message,precision:fixed.precision,evaluations:fixed.evaluations??0}};
 return fallback;
}
async function verifiedRoot(config,source,guess,chart,depth,deadline,metrics,signal,{fast=false,radius=.03,compactFixedReduction=false}={}){
 const endpointCorrect=compactFixedReduction?(Reduction.eligible(config)?compactReducedCorrect:compactEndpointCorrect):correct;
 let last;
 for(const extra of [0,30,60]){
  ensure(deadline,signal);const precision=localPrecision(config,depth,extra),digits=Math.max(24,Math.ceil(depth+20));
  const low=endpointCorrect({config,source,guess,chart,precision,digits,radius,maxStep:Math.min(.04,radius/2),deadline,fast});addMetrics(metrics,low);last=low;
  if(!low.converged)continue;
  const V=getEvaluator(precision+45),perturbed=V.strs(V.inverse(V.add(V.mobius(V.c(low.tau),chart),['1e-7','2e-7']),chart));
  const high=endpointCorrect({config,source,guess:perturbed,chart,precision:precision+35,digits:digits+10,radius:Math.min(.02,radius),maxStep:.01,deadline,fast:false});addMetrics(metrics,high);last=high;
  if(!high.converged)continue;
  const E=getEvaluator(precision+60),l=E.mobius(E.c(low.tau),chart),h=E.mobius(E.c(high.tau),chart),error=E.abs(E.sub(l,h));
  if(error.gt('1e-18'))continue;
  const lb=outputBranch(low.tau,{precision:low.precision}),hb=outputBranch(high.tau,{precision:high.precision,correction:error.toString(),crossPrecisionChartError:error.toString(),engine:Reduction.eligible(config)?'marked Lamé factor; repeated-precision correction':'adaptive marked predictor-corrector',reduction:Reduction.eligible(config)?Reduction.eligibility(config):null,evaluations:low.evaluations+high.evaluations});
  if(compactFixedReduction){hb.arithmetic=high.diagnostics?.arithmetic;hb.fourier=high.diagnostics?.fourier;}
  hb.viewErrors={};for(const k of Object.keys(hb.views))hb.viewErrors[k]=E.abs(E.sub(E.c(lb.views[k]),E.c(hb.views[k]))).toString();
  return hb;
 }
 const e=Error('Independent precision corrections did not resolve the same local root: '+(last?.message??'no corrected root'));e.corrector=last;
 // Exhausting arithmetic precision does not reject the marked graph. Keep
 // the checked prefix, and let the UI distinguish it from sheet disagreement.
 if(last?.failureKind==='precision'||/working precision|zero local|singular|division by zero/i.test(last?.message??''))e.failureKind='numerical_precision_limit';
 throw e;
}
function pathFor(request){
 const c=request.caseConfig,path=request.path??{},index=Number(path.index??0),position=String(path.position??'.5');
 if(path.kind&&path.kind!=='edge'||![0,1,2].includes(index)||!(Number(position)>0&&Number(position)<1))throw Error('Choose an open source edge.');
 const E=getEvaluator(300),p=new E.D(position),start=(path.start??[request.R??'.58',request.S??'.62']).map(String);
 if(!(Number(start[0])>0&&Number(start[0])<1&&Number(start[1])>0&&Number(start[1])<1&&Number(start[0])+Number(start[1])>1))throw Error('The starting point must lie inside the tetrahedral source triangle.');
 const end=[[new E.D(1),p],[p,new E.D(1)],[p,new E.D(1).minus(p)]][index];
 const source=epsilon=>end.map((v,j)=>v.plus(new E.D(String(epsilon)).times(new E.D(start[j]).minus(v))).toString());
 return{index,position,start,end:end.map(String),source,chart:charts[index],strictTracking:!!request.strictTracking,maxDepthStep:Math.max(.01,Math.min(1,Number(request.maxDepthStep)||1)),sourceDepth:d=>source(new E.D(10).pow(-d).toString())};
}
function predict(history,depth,chart,precision){
 const E=getEvaluator(precision),last=history.at(-1),a=E.mobius(E.c(last.tau),chart);
 if(history.length<2)return{tau:last.tau,predicted:a,displacement:0};
 const prior=history.at(-2),b=E.mobius(E.c(prior.tau),chart),ratio=(depth-last.depth)/(last.depth-prior.depth),move=E.mul(E.sub(a,b),String(ratio)),next=E.add(a,move);
 return{tau:E.strs(E.inverse(next,chart)),predicted:next,displacement:E.abs(move).toNumber()};
}
function chooseGate(config,path){
 if(Math.max(...config.weights)<10)return'.2';
 const clearance=source=>Math.min(Number(source[0]),1-Number(source[0]),Number(source[1]),1-Number(source[1]),Number(source[0])+Number(source[1])-1),goal=Math.min(.12,.6*clearance(path.start));
 // Choose a genuinely interior source, rather than applying the same ray
 // parameter to edges with different normal source distances.
 for(const epsilon of ['.2','.3','.4','.5','.65','.8'])if(clearance(path.source(epsilon))>=goal-1e-14)return epsilon;
 return'.9';
}
async function advance(config,history,target,path,deadline,metrics,signal){
 let step=history.at(-1)?.nextStep??(history.length<2?.035:.1);
 while(history.at(-1).depth<target-1e-12){
  ensure(deadline,signal);const last=history.at(-1),weight=Math.max(...config.weights),initialCap=.06/Math.sqrt(Math.max(1,weight/4)),maxStep=Math.min(path.maxDepthStep??1,last.depth<1?initialCap:last.depth<2?.15:last.depth>3&&last.predictorRatio<.035?1:.5),d=Math.min(target,last.depth+Math.min(step,maxStep));
  const precision=transportPrecision(config,d),E=getEvaluator(precision),guess=predict(history,d,path.chart,precision),source=path.sourceDepth(d),radius=Math.min(.15,Math.max(.025,2*guess.displacement+.012));
  let got=correct({config,source,guess:guess.tau,chart:path.chart,precision,digits:Math.max(17,Math.ceil(d+12)),radius,maxStep:Math.min(.04,radius/2),deadline,fast:d<=8,reuseFastDerivative:true,verifyResidual:false});addMetrics(metrics,got);
  let ratio=Infinity;
  if(got.converged){const z=E.mobius(E.c(got.tau),path.chart),correction=E.abs(E.sub(z,guess.predicted)).toNumber();ratio=correction/Math.max(.003,guess.displacement);}
  if(!got.converged||ratio>.25){
   metrics.rejectedSteps++;step=(d-last.depth)/2;if(step<1e-5)throw Error('The local marked branch could not be continued inside its prediction neighbourhood.');continue;
  }
  // Embedded continuation: two half steps must recover the same local sheet
  // as one full step. Precision agreement alone would miss a branch jump.
  const early=last.depth<1.1||(weight>=7&&last.depth<3),stepIndex=last.stepIndex??0;
  if(last.depth<.12||ratio>.06||(early&&(path.strictTracking||stepIndex%4===0))){
   const mid=(last.depth+d)/2,midGuess=predict(history,mid,path.chart,precision),midSource=path.sourceDepth(mid);
   const one=correct({config,source:midSource,guess:midGuess.tau,chart:path.chart,precision,digits:Math.max(17,Math.ceil(mid+12)),radius:Math.min(radius,.1),maxStep:.025,deadline,fast:mid<=8,reuseFastDerivative:true,verifyResidual:false});addMetrics(metrics,one);
   let two=null;
   if(one.converged){const history2=[...history,{depth:mid,source:midSource,tau:one.tau}],next=predict(history2,d,path.chart,precision);two=correct({config,source,guess:next.tau,chart:path.chart,precision,digits:Math.max(17,Math.ceil(d+12)),radius:Math.min(radius,.1),maxStep:.025,deadline,fast:d<=8,reuseFastDerivative:true,verifyResidual:false});addMetrics(metrics,two);}
   if(!one.converged||!two?.converged||chartDistance(got.tau,two.tau,path.chart,precision)>1e-11){metrics.rejectedSteps++;metrics.embeddedDisagreements=(metrics.embeddedDisagreements??0)+1;step=(d-last.depth)/2;if(step<1e-5)throw Error('The full and half-step continuations could not resolve the same local sheet.');continue;}
   metrics.embeddedChecks=(metrics.embeddedChecks??0)+1;got=two;
  }
  // A reverse correction must recover the preceding local sheet. This check
  // uses the equation and the preceding root, never a proposed cusp center.
  if(last.depth<.1||ratio>.15||(last.depth<1.1&&metrics.acceptedSteps%10===0)){
   const reverseRadius=Math.min(.18,Math.max(.025,2*chartDistance(got.tau,last.tau,path.chart,precision)+.01));
   const reverse=correct({config,source:last.source,guess:got.tau,chart:path.chart,precision,digits:17,radius:reverseRadius,maxStep:.04,deadline,fast:last.depth<=8,reuseFastDerivative:true,verifyResidual:false});addMetrics(metrics,reverse);
   if(!reverse.converged||chartDistance(reverse.tau,last.tau,path.chart,precision)>1e-12){metrics.rejectedSteps++;step=(d-last.depth)/2;if(step<1e-5)throw Error('The reverse branch check could not recover the preceding marked root.');continue;}
   metrics.reverseChecks++;
  }
  const proposed=Math.min(maxStep,(d-last.depth)*(ratio<.025?1.7:ratio<.1?1.3:1));
  const nextStep=d===target&&ratio<.1?Math.max(step,proposed):proposed;
  history.push({depth:d,epsilon:new E.D(10).pow(-d).toString(),source,tau:got.tau,precision:got.precision,predictorRatio:ratio,nextStep,stepIndex:stepIndex+1});
  if(history.length>4)history.shift();metrics.acceptedSteps++;
  step=nextStep;
  await new Promise(resolve=>setTimeout(resolve,0));
 }
 return history.at(-1);
}
export async function traceResult(request){
 const config=validate(request.caseConfig),path=pathFor(request),target=String(request.target??'1e-8'),depth=-Math.log10(Number(target)),density=Math.max(1,Math.min(5,Number(request.samplesPerDecade)||1)),budget=Math.max(50,Math.min(300000,Number(request.budgetMs)||45000)),begun=now(),deadline=begun+budget,gateEpsilon=chooseGate(config,path),gateDepth=-Math.log10(Number(gateEpsilon));
 const first=String(request.startEpsilon??1),startDepth=-Math.log10(Number(first));
 if(!(depth>=1&&depth<=80))throw Error('Choose a final epsilon from 1e-1 through 1e-80.');
 if(!(Number(first)>0&&Number(first)<=1&&startDepth<=depth))throw Error('The starting epsilon must lie between the target and one.');
 const signature=key({version:VERSION,config,path:{index:path.index,position:path.position,start:path.start},target,density,startEpsilon:request.startEpsilon??1,strictTracking:path.strictTracking,maxDepthStep:path.maxDepthStep});
 let job=request.resumeToken?jobs.get(request.resumeToken):null;
 if(request.resumeToken&&(!job||job.signature!==signature||!job.resumable))throw Error('This checkpoint does not match the selected labelled graph and source path.');
 const resumed=!!job;
 if(!job){
  const extension=request.extendFrom?traceExtension(request,'tv'):null;
  if(extension&&(extension.result.solverVersion!==VERSION||!extension.result.geometryChecks?.some(g=>g.kind==='start'&&g.checked&&g.accepted)||(!extension.result.branchValidation?.intermediateVerified&&Number(extension.last.epsilon)<=Number(gateEpsilon))))throw Error('This saved trace needs the updated marked-graph checks. Compute the approach again; display settings are retained.');
  const E=getEvaluator(300),steps=Math.max(1,Math.ceil((depth-startDepth)*density)),wanted=extension?extendTraceGrid(extension,target,density,d=>new E.D(10).pow(d).toString()):Array.from({length:steps+1},(_,i)=>i===0?first:i===steps?target:new E.D(10).pow(-startDepth-(depth-startDepth)*i/steps).toString());
  job={token:`tetra-${++serial}`,signature,nodes:[],history:[],baseHistory:[],geometryChecks:[],auditDone:false,index:0,resumable:false,metrics:{correctorCalls:0,evaluations:0,correctorMs:0,geometryCalls:0,geometryMs:0,anchorCacheHits:0,acceptedSteps:0,rejectedSteps:0,reverseChecks:0}};
  job.wanted=wanted;
  if(extension){
   job.nodes=extension.nodes;job.index=extension.reusedCount;job.auditDone=!!extension.result.branchValidation?.intermediateVerified;job.geometryChecks=clone(extension.result.geometryChecks??[]);job.reusedCount=extension.reusedCount;
   job.history=job.nodes.slice(-3).map(n=>({depth:-Math.log10(Number(n.epsilon)),epsilon:n.epsilon,source:[n.rExact,n.sExact],tau:n.result.tauExact,precision:n.result.precisionDigits}));
   job.baseHistory=job.nodes.slice(-3).map(n=>({depth:-Math.log10(Number(n.epsilon)),epsilon:n.epsilon,source:[n.rExact,n.sExact],tau:n.base.tauExact,precision:n.base.precisionDigits}));
   job.importPending=true;
  }
  jobs.set(job.token,job);while(jobs.size>8)jobs.delete(jobs.keys().next().value);
 }
 const finish=(status,message,stop=null)=>({status,nodes:clone(job.nodes),completed:job.nodes.length,total:job.wanted.length,elapsedMs:now()-begun,message,stop,resumeToken:job.resumable?job.token:null,resumed,reusedCount:job.reusedCount??0,resumeAvailable:job.resumable,extended:!!job.reusedCount,newCount:Math.max(0,job.nodes.length-(job.reusedCount??0)),precisionDescription:'Adaptive marked-graph continuation; healthy independent Hodge meshes, bounded corrections, embedded half steps and perturbed higher-precision repeats. Numerical estimates, not interval certificates.',geometryChecks:clone(job.geometryChecks),metrics:clone(job.metrics),solverVersion:VERSION,branchValidation:{version:VERSION,intermediateVerified:job.auditDone,gateEpsilon,reduction:Reduction.eligible(config)?Reduction.eligibility(config):null},path:{kind:'edge',index:path.index,position:path.position,start:path.start,endpoint:path.end}});
 try{
  ensure(deadline,request.signal);
  const base=baseConfig();
  if(job.importPending){
   const last=job.nodes.at(-1),source=path.source(last.epsilon),E=getEvaluator(300);
   if(source.some((x,j)=>!new E.D(x).eq([last.rExact,last.sExact][j])))throw Error('The imported checkpoint has an inconsistent exact source point.');
   for(const [c,part,history] of [[config,'result',job.history],[base,'base',job.baseHistory]]){
    const refined=await verifiedRoot(c,source,last[part].tauExact,path.chart,-Math.log10(Number(last.epsilon)),deadline,job.metrics,request.signal,{radius:.01});
    if(chartDistance(refined.tauExact,last[part].tauExact,path.chart,180)>1e-18)throw Error('The imported endpoint does not match its independent precision correction.');
    history.at(-1).tau=refined.tauExact;
   }
   job.importPending=false;
  }
  if(!job.history.length){
   const anchor=await independentAnchor(config,path.start,path.chart,deadline,job.metrics,request.signal),b=await independentAnchor(base,path.start,path.chart,deadline,job.metrics,request.signal);
   job.history=[{depth:0,epsilon:'1',source:path.start,tau:anchor.tau,precision:anchor.precision??90}];
   job.baseHistory=[{depth:0,epsilon:'1',source:path.start,tau:b.tau,precision:b.precision??90}];
   job.geometryChecks.push({kind:'start',source:path.start,...anchor});
  }
  const steps=job.wanted.length-1;
  for(;job.index<=steps;job.index++){
   ensure(deadline,request.signal);const epsilon=job.wanted[job.index],d=-Math.log10(Number(epsilon));
   if(d>=gateDepth&&!job.auditDone){
    await advance(config,job.history,gateDepth,path,deadline,job.metrics,request.signal);
    const independent=await independentAnchor(config,path.source(gateEpsilon),path.chart,deadline,job.metrics,request.signal),difference=chartDistance(independent.tau,job.history.at(-1).tau,path.chart,150);
    job.geometryChecks.push({kind:'intermediate',source:path.source(gateEpsilon),epsilon:gateEpsilon,continuedTau:job.history.at(-1).tau,difference,...independent});
    if(difference>1e-12)throw Error('Independent healthy geometric meshes disagree with the continued marked branch.');
    job.auditDone=true;
   }
   const selected=await advance(config,job.history,d,path,deadline,job.metrics,request.signal),baseSelected=isBase(config)?selected:await advance(base,job.baseHistory,d,path,deadline,job.metrics,request.signal),source=path.source(epsilon);
   // The undecorated comparison and selected curve are identical in this
   // case. Share the checked continuation rather than solving it twice.
   if(isBase(config))job.baseHistory=clone(job.history);
   const result=await verifiedRoot(config,source,selected.tau,path.chart,d,deadline,job.metrics,request.signal,{radius:.025}),baseResult=isBase(config)?result:await verifiedRoot(base,source,baseSelected.tau,path.chart,d,deadline,job.metrics,request.signal,{radius:.025});
   Object.assign(job.history.at(-1),{tau:result.tauExact,source,epsilon});Object.assign(job.baseHistory.at(-1),{tau:baseResult.tauExact,source,epsilon});
   const node={status:'ok',epsilon,r:Number(source[0]),s:Number(source[1]),rExact:source[0],sExact:source[1],base:baseResult,result};job.nodes.push(node);
   if(typeof request.onProgress==='function')request.onProgress({node,completed:job.nodes.length,total:steps+1,solverVersion:VERSION,branchValidation:{version:VERSION,intermediateVerified:job.auditDone},geometryChecks:job.geometryChecks});
  }
  job.resumable=false;return finish('ok','Computed marked roots with healthy independent graph checks; no cusp endpoint was inserted.');
 }catch(error){
  const cancelled=request.signal?.aborted,budgetFailure=error.failureKind==='time_budget'||/budget|deadline/i.test(error.message),inconclusive=error.failureKind==='geometry_inconclusive',precisionLimit=error.failureKind==='numerical_precision_limit';job.resumable=budgetFailure||cancelled||inconclusive;
  const stop={kind:cancelled?'cancelled':budgetFailure?'time_budget':inconclusive?'geometry_inconclusive':precisionLimit?'numerical_precision_limit':'verification',message:error.message,resumeAvailable:job.resumable,retryable:job.resumable,branchRejected:!job.resumable&&!precisionLimit,...(error.geometry?{geometry:error.geometry}:{}),...(error.corrector?{corrector:error.corrector}:{})};
  return finish(job.nodes.length?'partial':'unresolved',error.message,stop);
 }
}
export async function pointResult(request){
 const config=validate(request.caseConfig),source=[String(request.R),String(request.S)],begun=now(),deadline=request.deadline??begun+(request.budgetMs??45000),metrics={anchorCacheHits:0,geometryCalls:0,geometryMs:0,correctorCalls:0,evaluations:0,correctorMs:0};
 const a=await independentAnchor(config,source,[1,0,0,1],deadline,metrics,request.signal),b=await independentAnchor(baseConfig(),source,[1,0,0,1],deadline,metrics,request.signal);
 return{source:source.map(Number),result:await verifiedRoot(config,source,a.tau,[1,0,0,1],0,deadline,metrics,request.signal),base:await verifiedRoot(baseConfig(),source,b.tau,[1,0,0,1],0,deadline,metrics,request.signal),metrics,elapsedMs:now()-begun,solverVersion:VERSION};
}
export function clearSolverCaches(){anchors.clear();jobs.clear();clearGeometryCache();Reduction.clearCache();}

// Compact meshes reuse a single independent anchor, followed by short bounded
// continuations in the worker. A Hodge mesh is not rebuilt at every vertex.
export const compactBaseConfig=baseConfig;
export const compactCorrect=correct;
export async function compactRefine(config,source,guess,deadline=Infinity){
 const metrics={correctorCalls:0,evaluations:0,correctorMs:0};
 return verifiedRoot(config,source,guess,[1,0,0,1],0,deadline,metrics,undefined,{radius:.02,compactFixedReduction:true});
}
export async function compactAnchor(config,deadline=Infinity){
 const source=['.58','.62'],metrics={anchorCacheHits:0,geometryCalls:0,geometryMs:0};
 const got=await independentAnchor(validate(config),source,[1,0,0,1],deadline,metrics);
 return{source,tau:got.tau,precision:got.precision??90,provenance:Reduction.eligible(config)?'independently marked classical Lamé factor anchor':'independent healthy marked Hodge meshes with bounded root agreement',solverVersion:VERSION};
}

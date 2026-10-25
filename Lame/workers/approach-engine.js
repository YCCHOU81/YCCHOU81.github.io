/* Marked log-epsilon continuation and independent two-precision validation.
 * One implementation serves Limits and Overlap with their existing initializer
 * and precision policy. All continuation, acceptance and stopping checks remain. */
function createApproachMath(attachmentMath,{reuseAcceptedEvaluation=false,exposePrecisionEngine=false}={}) {
const {solveAttachment,validateCounts,setAnchorData,MAX_N}=attachmentMath;
const createPrecisionEngine=createLamePrecisionEngine({reuseAcceptedEvaluation});
const clock=()=>globalThis.performance?.now?.()??Date.now();
let storedAnchors=null;
function setApproachAnchorData(data){setAnchorData(data);storedAnchors=data;return{status:'ok',maxN:Math.min(MAX_N,data.maxN)};}

function normalizePath(path={}){
 const kind=path.kind??path.mode??'edge',index=Number(kind==='vertex'?(path.vertex??path.index??0):(path.edge??path.index??0));
 const position=String(kind==='vertex'?(path.direction??path.position??'.5'):(path.position??'.5'));
 if(!['edge','vertex'].includes(kind)||![0,1,2].includes(index))throw new RangeError('Choose an edge or vertex numbered 0, 1, or 2.');
 const t=new Decimal(position);if(!t.gt(0)||!t.lt(1))throw new RangeError('The path position must lie strictly between zero and one.');
 let start;
 if(path.start!==undefined){
  if(!Array.isArray(path.start)||path.start.length!==2)throw new RangeError('Choose an interior starting point with two coordinates.');
  start=path.start.map(value=>{
   const text=String(value);
   if(text.length>60||!/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(text))throw new RangeError('Use finite decimal coordinates for the starting point.');
   const coordinate=new Decimal(text);if(!coordinate.isFinite())throw new RangeError('Use finite decimal coordinates for the starting point.');
   return coordinate.toString();
  });
  const r=new Decimal(start[0]),s=new Decimal(start[1]);
  if(!r.gt(0)||!s.gt(0)||!r.lt('.5')||!s.lt('.5')||!r.plus(s).gt('.5'))throw new RangeError('The starting point must lie strictly inside the base monodromy triangle.');
 }
 return{kind,index,position,...(kind==='edge'?{edge:index}:{vertex:index,direction:position}),...(start?{start}:{})};
}
function coordinates(engine,path,epsilon){
 const {d}=engine,e=d(epsilon),t=d(path.position),b=[];
 if(!e.gt(0)||(path.start?e.gt(1):!e.lt(1)))throw new RangeError(path.start?'Segment epsilon must be positive and at most one.':'Approach epsilon must lie strictly between zero and one.');
 const other=[0,1,2].filter(i=>i!==path.index);
 if(path.start){
  const r=d(path.start[0]),s=d(path.start[1]),start=[d(1).minus(s.times(2)),r.plus(s).times(2).minus(1),d(1).minus(r.times(2))];
  const boundary=[d(0),d(0),d(0)];
  if(path.kind==='edge'){boundary[other[0]]=t;boundary[other[1]]=d(1).minus(t);}else boundary[path.index]=d(1);
  for(let i=0;i<3;i++)b[i]=boundary[i].times(d(1).minus(e)).plus(start[i].times(e));
 }else{
  b[path.index]=path.kind==='edge'?e:d(1).minus(e);
  const remaining=path.kind==='edge'?d(1).minus(e):e;
  b[other[0]]=remaining.times(t);b[other[1]]=remaining.times(d(1).minus(t));
 }
 return{epsilon:e,r:d(1).minus(b[2]).div(2),s:d(1).minus(b[0]).div(2),barycentric:b};
}
function approachCoordinates(path,epsilon,precisionDigits=70){
 const p=createPrecisionEngine(precisionDigits),point=coordinates(p,normalizePath(path),epsilon);
 return{epsilon:point.epsilon.toString(),r:point.r.toNumber(),s:point.s.toNumber(),
  rExact:point.r.toString(),sExact:point.s.toString(),barycentric:point.barycentric.map(v=>v.toString())};
}
function defaultApproachEpsilons(target='1e-6',samplesPerDecade=3){
 const last=new Decimal(target);if(!last.gt(0)||last.gt('.1')||last.lt('1e-10'))throw new RangeError('Choose a target epsilon from 0.1 through 1e-10.');
 const end=-Math.log10(last.toNumber()),out=[];
 for(let exponent=1;exponent<end-1e-12;exponent+=1/samplesPerDecade)
  out.push(Decimal.pow(10,-exponent).toPrecision(18));
 out.push(last.toString());return out;
}

/** Follow one geometrically labeled sheet to actual points approaching a face.
 * Epsilon is barycentric; the ordinary UI's source margin is deliberately not
 * used here. Returned coordinates remain in the fixed original F0 marking. */
async function solveApproachPath(counts,path,epsilons,options={}){
 const started=clock(),data=options.anchorData??storedAnchors;
 if(!data)throw new TypeError('Load the geometric attachment anchors before tracing an approach.');
 setApproachAnchorData(data);const checked=validateCounts(counts,data);
 if(checked.status!=='ok')throw new RangeError(checked.message);
 path=normalizePath(path);epsilons=epsilons??defaultApproachEpsilons(options.target??'1e-6');
 if(!Array.isArray(epsilons)||!epsilons.length||epsilons.length>100)throw new RangeError('Provide between one and 100 approach samples.');
 const precisionDigits=Math.max(60,Math.min(110,options.precisionDigits??70));
 const validationDigits=Math.max(precisionDigits+25,Math.min(160,options.validationDigits??100));
 const p=createPrecisionEngine(precisionDigits),v=createPrecisionEngine(validationDigits);
 const requested=epsilons.map(e=>p.d(String(e)));
 for(let i=0;i<requested.length;i++)if(!requested[i].gte('1e-10')||requested[i].gt(path.start?'1':'.3')||(i&&requested[i].gte(requested[i-1])))
  throw new RangeError(path.start?'Segment epsilons must strictly decrease from at most 1 to at least 1e-10.':'Epsilons must strictly decrease from at most 0.3 to at least 1e-10.');
 const key=counts.join(','),entry=data.anchors[key],baseEntry=data.anchors['0,0,0'];
 if(!entry||!baseEntry)throw new RangeError('This geometric sheet has no validated anchor.');
 const shift=[p.d(counts[0]+counts[2]).div(2),p.d(counts[1]+counts[2]).div(2)];
 let baseTau=p.C(...(baseEntry.highPrecisionTau??baseEntry.tau)),selectedTau=p.C(...(entry.highPrecisionTau??entry.tau));
 let currentR=p.d(String(data.r)),currentS=p.d(String(data.s)),currentEpsilon=null,steps=0,totalIterations=0;
 let stateLogEpsilon=null,previousLogEpsilon=null,previousBaseTau=null,previousSelectedTau=null;
 const retained=options.resumeState?.nodes?.filter(node=>node.status==='ok')??[];
 const nodes=[...retained],cancelled=()=>Boolean(options.signal?.aborted),yieldUI=()=>new Promise(resolve=>setTimeout(resolve,0));
 async function correctAt(point,initial=false){
  let baseSeed=baseTau,selectedSeed=selectedTau;
  const nextLogEpsilon=point.epsilon?.ln()??null;
  if(!initial&&nextLogEpsilon&&stateLogEpsilon&&previousLogEpsilon&&!stateLogEpsilon.eq(previousLogEpsilon)){
   // A local secant predictor in log(epsilon) reduces Newton work near a
   // cusp. The accepted root still obeys the original continuation guard.
   const ratio=nextLogEpsilon.minus(stateLogEpsilon).div(stateLogEpsilon.minus(previousLogEpsilon));
   if(ratio.gt(0)&&ratio.lte(3)){
    const predict=(tau,previous)=>{
     let move=p.K(p.S(tau,previous),ratio),distance=p.abs(move);
     if(distance.gt('.12'))move=p.K(move,p.d('.12').div(distance));
     return p.A(tau,move);
    };
    baseSeed=predict(baseTau,previousBaseTau);selectedSeed=predict(selectedTau,previousSelectedTau);
   }
  }
  const base=p.correct(baseSeed,point.r,point.s,1,{tolerance:'1e-27',trustRadius:initial?.2:.3});
  if(base.status!=='ok')return{status:base.status,branch:'base'};
  const selected=checked.n===1?base:p.correct(selectedSeed,point.r.plus(shift[0]),point.s.plus(shift[1]),checked.n,{tolerance:'1e-27',trustRadius:initial?.2:.3});
  if(selected.status!=='ok')return{status:selected.status,branch:'result'};
  if(p.abs(p.S(base.tau,baseTau)).gt('.16')||p.abs(p.S(selected.tau,selectedTau)).gt('.16'))return{status:'step_too_large'};
  if(nextLogEpsilon){previousLogEpsilon=stateLogEpsilon;stateLogEpsilon=nextLogEpsilon;previousBaseTau=baseTau;previousSelectedTau=selectedTau;}
  baseTau=base.tau;selectedTau=selected.tau;currentR=point.r;currentS=point.s;
  totalIterations+=base.iterations+(checked.n===1?0:selected.iterations);steps++;await yieldUI();return{status:'ok'};
 }
 if(retained.length){
  const last=retained.at(-1);baseTau=p.asComplex(last.base.tauExact);selectedTau=p.asComplex(last.result.tauExact);
  currentR=p.d(last.rExact);currentS=p.d(last.sExact);currentEpsilon=p.d(last.epsilon);stateLogEpsilon=currentEpsilon.ln();
  if(options.revalidateRetainedSeed){
   // A transported cache entry supplies data, never a live numerical state.
   // Recheck its final root at both precisions before using it to continue;
   // retained output strings themselves remain exactly as originally checked.
   const point=coordinates(p,path,currentEpsilon),oldBase=baseTau,oldSelected=selectedTau;
   const base=p.correct(oldBase,point.r,point.s,1,{tolerance:'1e-27',trustRadius:.025});
   const selected=checked.n===1?base:p.correct(oldSelected,point.r.plus(shift[0]),point.s.plus(shift[1]),checked.n,{tolerance:'1e-27',trustRadius:.025});
   if(base.status!=='ok'||selected.status!=='ok'||p.abs(p.S(base.tau,oldBase)).gt('1e-17')||p.abs(p.S(selected.tau,oldSelected)).gt('1e-17'))
    throw new Error('The saved endpoint no longer passes the branch correction checks. Compute a new approach.');
   baseTau=base.tau;selectedTau=selected.tau;
   const baseCheck=await checkedBranch(baseTau,1,0,0),selectedCheck=checked.n===1?baseCheck:await checkedBranch(selectedTau,checked.n,shift[0].toString(),shift[1].toString());
   if(baseCheck.status!=='ok'||selectedCheck.status!=='ok')throw new Error('The saved endpoint no longer passes the independent precision check. Compute a new approach.');
  }
  if(retained.length>1){const prior=retained.at(-2);previousLogEpsilon=p.d(prior.epsilon).ln();previousBaseTau=p.asComplex(prior.base.tauExact);previousSelectedTau=p.asComplex(prior.result.tauExact);}
  steps=last.continuationSteps??0;
  for(let i=0;i<retained.length;i++)if(options.onProgress)await options.onProgress(retained[i],i+1,requested.length);
 }else{
 const first=coordinates(p,path,requested[0]);
 // The existing double engine supplies only an initial root on this exact
 // labeled sheet; multiprecision correction and comparison accept every node.
 const warm=solveAttachment(first.r.toNumber(),first.s.toNumber(),counts,null,data);
 if(warm.status==='ok'){
  baseTau=p.asComplex(warm.base.tau);selectedTau=p.asComplex(warm.result.tau);
  const result=await correctAt(first,true);if(result.status!=='ok')throw new Error('The high-precision initial correction failed: '+result.status);
 }else{
  const originR=currentR,originS=currentS,distance=Math.hypot(first.r.minus(originR).toNumber(),first.s.minus(originS).toNumber());
  let fraction=0,increment=Math.min(1,.012/Math.max(distance,.001));
  for(let attempts=0;fraction<1&&attempts<300;attempts++){
   if(cancelled())return{status:'aborted',counts,path,nodes,elapsedMs:clock()-started};
   const next=Math.min(1,fraction+increment),point={r:originR.plus(first.r.minus(originR).times(next)),s:originS.plus(first.s.minus(originS).times(next))};
   const result=await correctAt(point,true);
   if(result.status==='ok'){fraction=next;increment=Math.min(1-fraction,increment*1.3);}else{increment/=2;if(increment<1e-6)throw new Error('The geometric anchor could not reach the requested approach path.');}
  }
  if(fraction<1)throw new Error('The initial continuation exceeded its safe step budget.');
 }
 }
 async function checkedBranch(tau,n,dr,ds){
  const point=coordinates(v,path,currentEpsilon.toString());
  const initial=v.C(tau[0].toString(),tau[1].toString());
  const corrected=v.correct(initial,point.r.plus(dr),point.s.plus(ds),n,{tolerance:'1e-38',trustRadius:.025,maxIterations:16});
  if(corrected.status!=='ok')return{status:'precision_limited',message:'The independent higher-precision correction did not converge.'};
  const error=v.abs(v.S(corrected.tau,initial));
  if(error.gt('1e-17'))return{status:'precision_limited',crossPrecisionTauError:error.toString(),message:'The two working precisions did not agree to the required accuracy.'};
  // An unchanged higher-precision iterate is not an exact root. Include its
  // remaining Newton correction so zero inter-precision difference cannot
  // create an inflated accuracy label. This remains a numerical estimate.
  const numericalErrorEstimate=v.D.max(error,corrected.correction,v.d(10).pow(-precisionDigits+15));
  const agreementDigits=Math.max(0,Math.min(precisionDigits-15,Math.floor(-numericalErrorEstimate.log(10).toNumber())));
  return{status:'ok',n,tau:{re:corrected.tau[0].toNumber(),im:corrected.tau[1].toNumber()},
   tauExact:corrected.tau.map(z=>z.toPrecision(Math.min(validationDigits,60))),precisionDigits,validationDigits,
   precisionAgreementDigits:agreementDigits,crossPrecisionTauError:error.toString(),numericalErrorEstimate:numericalErrorEstimate.toString(),
   correctionEstimate:corrected.correction.toString(),iterations:corrected.iterations,validation:'two_precisions_not_interval_certified'};
 }
 for(let index=retained.length;index<requested.length;index++){
  if(cancelled())return{status:'aborted',counts:[...counts],path,nodes,elapsedMs:clock()-started};
  const target=requested[index];let failure=null;
  if(currentEpsilon!==null){
   let log=currentEpsilon.ln(),targetLog=target.ln(),increment=p.d('.8');
   for(let attempts=0;log.gt(targetLog)&&attempts<500;attempts++){
    if(cancelled())return{status:'aborted',counts:[...counts],path,nodes,elapsedMs:clock()-started};
    const next=p.D.max(targetLog,log.minus(increment)),point=coordinates(p,path,next.exp());
    const corrected=await correctAt(point);
    if(corrected.status==='ok'){log=next;increment=p.D.min('1.2',increment.times('1.25'));}
    else{increment=increment.div(2);if(increment.lt('1e-5')){failure=corrected;break;}}
   }
   if(log.gt(targetLog)&&!failure)failure={status:'continuation_limit'};
  }
  currentEpsilon=target;const point=coordinates(p,path,target);
  let base,result;
  if(failure){base={status:'unresolved'};result={status:failure.status,message:'Continuation stopped before this node; no other root was substituted.'};}
  else{
   base=await checkedBranch(baseTau,1,0,0);
   result=checked.n===1?{...base}:await checkedBranch(selectedTau,checked.n,shift[0].toString(),shift[1].toString());
  }
  const node={epsilon:target.toString(),r:point.r.toNumber(),s:point.s.toNumber(),rExact:point.r.toString(),sExact:point.s.toString(),
   barycentric:point.barycentric.map(z=>z.toString()),base,result,status:base.status==='ok'&&result.status==='ok'?'ok':'unresolved',
   elapsedMs:clock()-started,continuationSteps:steps};
  nodes.push(node);if(options.onProgress)await options.onProgress(node,index+1,requested.length);await yieldUI();
  if(node.status!=='ok')break;
 }
 return{status:nodes.length===requested.length&&nodes.every(n=>n.status==='ok')?'ok':'partial',counts:[...counts],n:checked.n,path,nodes,
  requestedEpsilons:requested.map(z=>z.toString()),completed:nodes.filter(n=>n.status==='ok').length,total:requested.length,
  stoppedAt:nodes.at(-1)?.status==='ok'?null:nodes.at(-1)?.epsilon,precisionDigits,validationDigits,elapsedMs:clock()-started,
  continuationSteps:steps,iterations:totalIterations,method:'geometric_label_continuation_and_two_decimal_precisions',
  epsilonDefinition:path.start?'Segment fraction: weights=(1-epsilon)*boundaryWeights+epsilon*startWeights; epsilon=1 is the selected start and epsilon=0 is the ideal boundary endpoint.':'Normalized barycentric face or corner distance; weights=(1-2s,2(r+s)-1,1-2r).'};
}

const api={setApproachAnchorData,defaultApproachEpsilons,solveApproachPath,normalizePath,approachCoordinates};
if(exposePrecisionEngine)api.createPrecisionEngine=createPrecisionEngine;
return api;
}

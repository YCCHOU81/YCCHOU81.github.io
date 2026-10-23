/* Preview-only cusp trajectories. All plotted nodes are corrected numerical roots.
 * The geometric anchor selects the branch. Small steps in log(epsilon) carry
 * it; an independently corrected higher-precision root checks every output.
 */
"use strict";
const previewTailTrace = (()=>{
 const finite=z=>z&&Number.isFinite(z.re)&&Number.isFinite(z.im)&&z.im>0;
 const yieldUI=()=>new Promise(resolve=>setTimeout(resolve,0));
 const cusp=['infinity','one','zero'];
 const epsilonKey=x=>Number(x).toPrecision(12);
 function grid(depth=5,density=1,inset=.06){
  depth=Math.max(3,Math.min(20,Math.round(Number(depth)||5)));
  density=Math.max(1,Math.min(3,Math.round(Number(density)||1)));
  inset=Number(inset);if(!(inset>=.009&&inset<1/3))throw RangeError('Use a barycentric inset between .009 and 1/3.');
  // The first tail sample is a genuine root on the compact mesh boundary.
  // Short transition steps replace the old uncomputed .06-to-.04 gap, while
  // retaining .04 and every historical logarithmic node for cache reuse.
  const joinEnd=Math.min(.04,inset),values=[];
  if(inset>joinEnd)for(let i=0;i<4;i++)values.push(Number(inset+(joinEnd-inset)*i/4).toPrecision(15));
  values.push(String(joinEnd));
  for(let exponent=2;exponent<depth-1e-10;exponent+=1/density){const epsilon=10**(-exponent);if(epsilon<joinEnd-1e-14)values.push(String(epsilon));}
  values.push('1e-'+depth);return values;
 }
 async function trace(counts,edge,t,epsilons,options={}){
  const started=performance.now(),n=1+counts.reduce((a,b)=>a+b,0),path={kind:'edge',edge,index:edge,position:String(t)};
  const requestedDepth=Math.max(...epsilons.map(epsilon=>-Math.log10(Number(epsilon)))),precisionDigits=requestedDepth>10?70:45,validationDigits=precisionDigits>45?100:70;
  const p=previewApproach.createPrecisionEngine(precisionDigits),v=previewApproach.createPrecisionEngine(validationDigits),shift=[(counts[0]+counts[2])/2,(counts[1]+counts[2])/2];
  const checked=previewTailAttachment.validateCounts(counts);if(checked.status!=='ok')throw RangeError(checked.message);
  const existing=new Map((options.existingTail?.nodes??[]).filter(node=>node.status==='ok'&&finite(node.tau)&&node.tauExact?.length===2).map(node=>[epsilonKey(node.epsilon),node]));
  const nodes=[],budget=options.timeBudgetMs??45000;
  let state=null,previousState=null,resolved=0,reused=0,stopReason=null,steps=0;
  const withinBudget=()=>performance.now()-started<=budget;
  const coordinates=epsilon=>previewApproach.approachCoordinates(path,epsilon,70);
  function asState(node){return{log:Math.log(Number(node.epsilon)),tau:p.asComplex(node.tauExact),epsilon:node.epsilon};}
  async function initial(epsilon){
   const source=coordinates(epsilon),ordinary=previewTailAttachment.solveAttachment(source.r,source.s,counts,null);
   if(ordinary.result?.status==='ok'&&finite(ordinary.result.tau))return{tau:p.asComplex(ordinary.result.tau),log:Math.log(Number(epsilon)),epsilon};
   // Only this same geometric sheet may supply a fallback initial root.
   const slow=await previewApproach.solveApproachPath(counts,path,[epsilon],{anchorData:attachmentAnchors,precisionDigits:60,validationDigits:85,
    signal:{get aborted(){return!withinBudget();}}});
   const node=slow.nodes?.[0];
   if(node?.status!=='ok'||!node.result?.tauExact)throw Error('The geometric branch could not reach the first tail sample.');
   return{tau:p.asComplex(node.result.tauExact),log:Math.log(Number(epsilon)),epsilon};
  }
  function correct(seed,epsilon,trust=.12){
   const source=coordinates(epsilon),r=p.d(source.rExact).plus(String(shift[0])),s=p.d(source.sExact).plus(String(shift[1]));
   return p.correct(seed,r,s,n,{tolerance:precisionDigits>45?'1e-30':'1e-20',trustRadius:trust,maxIterations:20});
  }
  function verify(tau,epsilon){
   const source=coordinates(epsilon),refined=correct(tau,epsilon,.025);
   if(refined.status!=='ok')return{status:refined.status,source};
   const candidate=v.correct(refined.tau.map(x=>x.toString()),v.d(source.rExact).plus(String(shift[0])),v.d(source.sExact).plus(String(shift[1])),n,
    {tolerance:precisionDigits>45?'1e-40':'1e-30',trustRadius:.025,maxIterations:18});
   if(candidate.status!=='ok')return{status:candidate.status,source};
   const discrepancy=v.abs(v.S(candidate.tau,v.asComplex(refined.tau.map(x=>x.toString()))));
   const displacement=v.abs(v.S(candidate.tau,v.asComplex(tau.map(x=>x.toString()))));
   if(discrepancy.gt('1e-12')||displacement.gt('.005'))return{status:'precision_limited',source};
   const error=v.D.max(discrepancy,candidate.correction,v.d('1e-30'));
   return{status:'ok',source,tau:{re:candidate.tau[0].toNumber(),im:candidate.tau[1].toNumber()},tauExact:candidate.tau.map(x=>x.toPrecision(55)),
    errorEstimate:error.toString(),comparisonTauError:discrepancy.toNumber(),precisionDigits,validationDigits,
    method:'marked_log_epsilon_continuation_and_two_fourier_precisions'};
  }
  for(const epsilon of epsilons){
   const saved=existing.get(epsilonKey(epsilon));
   if(saved){previousState=state;state=asState(saved);nodes.push({...saved,reused:true});resolved++;reused++;
    if(options.onNode)await options.onNode(nodes.at(-1),nodes.length,epsilons.length);continue;}
   if(!withinBudget()){stopReason='time_budget';break;}
   let failure=null;
   try{
    if(!state)state=await initial(epsilon);
    else{
     const targetLog=Math.log(Number(epsilon));let increment=.45,attempts=0;
     while(state.log>targetLog+1e-13&&attempts++<400){
      if(!withinBudget()){failure='time_budget';break;}
      const nextLog=Math.max(targetLog,state.log-increment),nextEpsilon=String(Math.exp(nextLog));let seed=state.tau;
      if(previousState&&Math.abs(state.log-previousState.log)>1e-12){
       const ratio=(nextLog-state.log)/(state.log-previousState.log),move=p.K(p.S(state.tau,previousState.tau),String(ratio));
       if(ratio>0&&ratio<=3&&p.abs(move).lt('.08'))seed=p.A(seed,move);
      }
      const result=correct(seed,nextEpsilon);
      if(result.status==='ok'&&p.abs(p.S(result.tau,state.tau)).lt('.08')){
       previousState=state;state={tau:result.tau,log:nextLog,epsilon:nextEpsilon};steps++;increment=Math.min(.6,increment*1.2);
      }else{increment/=2;if(increment<1e-5){failure=result.status==='ok'?'continuation_guard':result.status;break;}}
      if(steps%4===0)await yieldUI();
     }
     if(state.log>targetLog+1e-13&&!failure)failure='continuation_limit';
    }
    let value=failure?{status:failure,source:coordinates(epsilon)}:verify(state.tau,epsilon);
    const source=value.source,node={epsilon:String(epsilon),r:source.r,s:source.s,source:{r:source.r,s:source.s},rExact:source.rExact,sExact:source.sExact,
     barycentric:source.barycentric,tau:value.status==='ok'?value.tau:null,status:value.status==='ok'?'ok':'unresolved',solverStatus:value.status,
     ...(value.status==='ok'?{tauExact:value.tauExact,errorEstimate:value.errorEstimate,comparisonTauError:value.comparisonTauError,
       precisionDigits,validationDigits,method:value.method}:{}),elapsedMs:performance.now()-started,continuationSteps:steps};
    nodes.push(node);if(node.status==='ok'){resolved++;state.tau=p.asComplex(node.tauExact);}else stopReason=failure??'numerical_verification';
    if(options.onNode)await options.onNode(node,nodes.length,epsilons.length);
    if(node.status!=='ok')break;
   }catch(error){
    stopReason=error.message;const source=coordinates(epsilon),node={epsilon:String(epsilon),r:source.r,s:source.s,source:{r:source.r,s:source.s},tau:null,status:'unresolved',solverStatus:'initialization_failed',message:error.message};
    nodes.push(node);if(options.onNode)await options.onNode(node,nodes.length,epsilons.length);break;
   }
   await yieldUI();
  }
  return{edge,t,position:t,path,cusp:cusp[edge],nodes,status:resolved===epsilons.length?'ok':'partial',resolved,total:epsilons.length,
   transitionInset:Number(options.inset??epsilons[0]),transitionVersion:'join-v1',
   requestedEpsilons:[...epsilons],reused,stopReason,elapsedMs:performance.now()-started,continuationSteps:steps,
   method:'marked_log_epsilon_continuation_and_two_fourier_precisions',numericalNotIntervalCertified:true};
 }
 return{grid,trace};
})();

// Independent marked-graph selection. No cusp formulas enter this module.
import{hemisphereMesh,markedHodge,meshHealth}from'./bouquet-hodge.mjs?v=20261024-bouquet-stable1';
import{correct as defaultCorrect,configKey}from'./bouquet-continuation.mjs?v=20261024-tv-audit2';
import{createEvaluator}from'./decimal-generic.mjs?v=20261024-tv-audit2';
const E=createEvaluator(200),now=()=>performance.now(),cache=new Map(),partialCache=new Map();
export const VERSION='bouquet-local-cone-hodge-2';
const error=(a,b)=>E.abs(E.sub(E.c(a),E.c(b)));
const clone=x=>JSON.parse(JSON.stringify(x));
const normalizePlan=o=>({sub:o.sub??6,rings:o.rings??4,grades:o.grades??8,regularize:!!o.regularize});
const savePartial=(key,state)=>{partialCache.delete(key);partialCache.set(key,clone(state));while(partialCache.size>16)partialCache.delete(partialCache.keys().next().value);};
const record=(state,attempt)=>{state.attempts.push(attempt);const budgets=state.attempts.filter(a=>a.status==='time-budget');if(budgets.length>16){const remove=new Set(budgets.slice(0,budgets.length-16));state.attempts=state.attempts.filter(a=>!remove.has(a));}};
export async function anchor(config,source=['.29','.31'],{deadline=now()+45000,correct=defaultCorrect,plans,precision=100,digits=25,jetScale=16,fresh=false,cacheScope}={}){
 source=source.map(String);const begun=now(),fullKey=configKey(config),policyPlans=(plans??[{sub:8,rings:5,grades:8,regularize:true},{sub:12,rings:7,grades:12,regularize:true},{sub:16,rings:9,grades:16,regularize:true},{sub:20,rings:11,grades:20,regularize:true}]).map(normalizePlan),uniquePlans=[...new Map(policyPlans.map(p=>[JSON.stringify(p),p])).values()],key=JSON.stringify([VERSION,fullKey,source,precision,digits,jetScale,uniquePlans,...(cacheScope===undefined?[]:[{cacheScope}])]);
 if(fresh){cache.delete(key);partialCache.delete(key);}
 if(cache.has(key)){const result=clone(cache.get(key));cache.delete(key);cache.set(key,clone(result));return{...result,cacheHit:true,milliseconds:now()-begun};}
 const resumed=partialCache.has(key),state=resumed?clone(partialCache.get(key)):{nextPlan:0,attempts:[],candidates:[],budgetStops:0,elapsedMilliseconds:0},reusedPlans=state.candidates.length;
 if(resumed){partialCache.delete(key);partialCache.set(key,clone(state));}
 let budgetFailure=false;
 for(let pi=state.nextPlan;pi<uniquePlans.length;pi++){
  const option=uniquePlans[pi],planId=JSON.stringify(option);state.nextPlan=pi;
  if(now()>deadline){budgetFailure=true;break;}const attempt={option,planIndex:pi};
  try{
   let t=now();const mesh=hemisphereMesh(config.graph,'bouquet',...source.map(Number),option);attempt.meshMilliseconds=now()-t;attempt.health=meshHealth(mesh);if(!attempt.health.healthy){attempt.status='invalid-mesh';record(state,attempt);state.nextPlan=pi+1;continue;}
   t=now();const h=markedHodge(config.graph,'bouquet',...source.map(Number),{...option,mesh});attempt.hodge=h;attempt.hodgeMilliseconds=now()-t;
   if(now()>deadline)throw Error('Geometry budget reached after marked Hodge computation');
   const radius=Math.min(.16,.055*(1+Math.hypot(...h.tau))),got=correct({config,source,guess:h.tau,precision,digits,radius,maxStep:.025,jetScale,deadline});attempt.refinement=got;
   if(!got.converged){attempt.status=/budget|deadline/i.test(got.message)?'time-budget':'refinement-inconclusive';record(state,attempt);if(attempt.status==='time-budget'){budgetFailure=true;break;}state.nextPlan=pi+1;continue;}
   attempt.hodgeDisplacement=error(got.tau,h.tau).toString();
   // A second independent mesh can confirm an already repeated HIGH root.
   // We return that high root, so no duplicate high solve is needed here.
   // An unmatched low root is never retained as a partial candidate.
   for(const previous of state.candidates.filter(c=>c.planId!==planId)){
    const agreement=error(got.tau,previous.tau),rawDifference=Math.hypot(...h.tau.map((v,i)=>v-previous.hodge[i]));
    if(agreement.gt('1e-19')||rawDifference>.12*(1+Math.hypot(...h.tau)))continue;
    attempt.status='independent-consensus';attempt.agreementWithVerifiedHigh=agreement.toString();record(state,attempt);state.nextPlan=pi+1;
    const elapsed=now()-begun,result={status:'accepted',accepted:true,checked:true,certified:false,source,tau:previous.tau,precision:previous.precision,configKey:fullKey,formulaInputsUsed:false,geometryVersion:VERSION,
      cacheHit:false,resumedGeometry:resumed,reusedGeometryPlans:reusedPlans,budgetStops:state.budgetStops,consensus:{meshOptions:[previous.option,option],hodgeDifference:rawDifference,rootDifference:agreement.toString(),perturbedRepeatDifference:previous.localRepeatDifference,localRepeatDifferences:[previous.localRepeatDifference],returnedHighRootPlan:previous.option,independentConfirmationPrecision:got.precision},attempts:clone(state.attempts),milliseconds:elapsed,totalMilliseconds:state.elapsedMilliseconds+elapsed};
    partialCache.delete(key);cache.set(key,clone(result));while(cache.size>12)cache.delete(cache.keys().next().value);return clone(result);
   }
   // Only a locally repeated, healthy single-mesh candidate is resumable.
   // It is still not an accepted graph anchor until another mesh agrees.
   const guess=E.strs(E.add(E.c(got.tau),['1e-7','2e-7'])),check=correct({config,source,guess,precision:precision+35,digits:digits+10,radius:.025,maxStep:.01,jetScale,deadline});attempt.perturbedRepeat=check;
   if(!check.converged){attempt.status=/budget|deadline/i.test(check.message)?'time-budget':'repeat-inconclusive';record(state,attempt);if(attempt.status==='time-budget'){budgetFailure=true;break;}state.nextPlan=pi+1;continue;}
   const repeatError=error(check.tau,got.tau);attempt.perturbedRepeatDifference=repeatError.toString();
   if(repeatError.gt('1e-19')){attempt.status='repeat-disagreement';record(state,attempt);state.nextPlan=pi+1;continue;}
   attempt.status='verified-candidate';record(state,attempt);
   const candidate={tau:check.tau,precision:check.precision,option,planId,hodge:h.tau,health:attempt.health,localRepeatDifference:repeatError.toString()};
   // A repeated physical mesh can never supply the second independent vote.
   if(!state.candidates.some(c=>c.planId===planId))state.candidates.push(candidate);
   state.nextPlan=pi+1;
   for(const previous of state.candidates.filter(c=>c.planId!==planId)){
    const agreement=error(candidate.tau,previous.tau),rawDifference=Math.hypot(...candidate.hodge.map((v,i)=>v-previous.hodge[i]));
    if(agreement.gt('1e-19')||rawDifference>.12*(1+Math.hypot(...candidate.hodge)))continue;
    const elapsed=now()-begun,result={status:'accepted',accepted:true,checked:true,certified:false,source,tau:candidate.tau,precision:candidate.precision,configKey:fullKey,formulaInputsUsed:false,geometryVersion:VERSION,
      cacheHit:false,resumedGeometry:resumed,reusedGeometryPlans:reusedPlans,budgetStops:state.budgetStops,consensus:{meshOptions:[previous.option,candidate.option],hodgeDifference:rawDifference,rootDifference:agreement.toString(),perturbedRepeatDifference:repeatError.toString(),localRepeatDifferences:[previous.localRepeatDifference,candidate.localRepeatDifference]},attempts:clone(state.attempts),milliseconds:elapsed,totalMilliseconds:state.elapsedMilliseconds+elapsed};
    partialCache.delete(key);cache.set(key,clone(result));while(cache.size>12)cache.delete(cache.keys().next().value);return clone(result);
   }
  }catch(e){attempt.status=/budget|deadline/i.test(String(e))?'time-budget':'geometry-inconclusive';attempt.message=String(e);record(state,attempt);if(attempt.status==='time-budget'){budgetFailure=true;break;}state.nextPlan=pi+1;}
  savePartial(key,state);
 }
 const elapsed=now()-begun;state.elapsedMilliseconds+=elapsed;if(budgetFailure){state.budgetStops++;savePartial(key,state);}else partialCache.delete(key);
 return{status:'inconclusive',accepted:false,checked:false,certified:false,source,configKey:fullKey,formulaInputsUsed:false,geometryVersion:VERSION,attempts:clone(state.attempts),milliseconds:elapsed,totalMilliseconds:state.elapsedMilliseconds,resumedGeometry:resumed,reusedGeometryPlans:reusedPlans,completedGeometryPlans:state.candidates.length,nextGeometryPlan:state.nextPlan,budgetStops:state.budgetStops,failureKind:budgetFailure?'time_budget':'geometry_inconclusive',stopKind:budgetFailure?'time_budget':'verification',retryable:true,branchRejected:false,message:budgetFailure?'Independent geometric selection reached its budget; completed verified meshes are retained':'Independent healthy graph refinements did not yet agree'};
}
export const clearCache=()=>{cache.clear();partialCache.clear();};
export const cacheStats=()=>({version:VERSION,accepted:cache.size,partial:partialCache.size,acceptedLimit:12,partialLimit:16});

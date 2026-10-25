// Validated, adaptive marked-graph branch selection. Used by the TV numerical worker.
import {hemisphereMesh,markedHodge} from './double-marked-hodge.mjs?v=20261023-double1';
import {hemisphereMesh as coneMesh,markedHodge as coneHodge,meshHealth as coneHealth} from './bouquet-hodge.mjs?v=20261024-bouquet-stable1';
import {createEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';

export const GEOMETRY_VERSION='tetra-health-bounded-v3-medium-cone';
const cache=new Map(),evaluators=new Map(),clock=()=>performance.now();
const engine=p=>{if(!evaluators.has(p))evaluators.set(p,createEvaluator(p));return evaluators.get(p);};
const abs=z=>Math.hypot(...z.map(Number));
const difference=(a,b)=>Math.hypot(Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]));
const identity=[1,0,0,1];
const sourceArray=source=>Array.isArray(source)?source.map(String):[String(source.R),String(source.S)];
const isDeadlineError=error=>/budget|deadline|time limit/i.test(String(error));

export function meshHealth(mesh){
 let zeroAreas=0,nonFinite=0,minArea=Infinity,minRelativeArea=Infinity,maxCotangent=0;
 for(const f of mesh.faces){
  const [a,b,c]=f.map(i=>mesh.xy[i]),u=[b[0]-a[0],b[1]-a[1]],v=[c[0]-a[0],c[1]-a[1]];
  const area=Math.abs(u[0]*v[1]-u[1]*v[0]);
  if(!Number.isFinite(area)){nonFinite++;continue;}if(area===0){zeroAreas++;continue;}
  minArea=Math.min(minArea,area);const scale=Math.max(u[0]**2+u[1]**2,v[0]**2+v[1]**2,(u[0]-v[0])**2+(u[1]-v[1])**2);
  minRelativeArea=Math.min(minRelativeArea,area/scale);
  for(const [p,q,r] of [[a,b,c],[b,c,a],[c,a,b]]){
   const x=[p[0]-r[0],p[1]-r[1]],y=[q[0]-r[0],q[1]-r[1]],cot=(x[0]*y[0]+x[1]*y[1])/area;
   if(!Number.isFinite(cot))nonFinite++;else maxCotangent=Math.max(maxCotangent,Math.abs(cot));
  }
 }
 return{healthy:zeroAreas===0&&nonFinite===0,vertices:mesh.nv,triangles:mesh.faces.length,
  zeroAreas,nonFinite,minArea:Number.isFinite(minArea)?minArea:null,minRelativeArea:Number.isFinite(minRelativeArea)?minRelativeArea:null,maxCotangent};
}

export function checkedHemisphereMesh(graph,seed,R,S,option,{deadline=Infinity}={}){
 const requested={...option},attempts=[],grades=[option.grades,...[40,32,24,16,12].filter(g=>g<option.grades)];
 for(const grade of [...new Set(grades)]){
  if(clock()>deadline)return{status:'inconclusive',failureKind:'time_budget',stopKind:'time_budget',timeBudgetReached:true,retryable:true,message:'Geometry time budget reached before the next mesh.',requested,attempts};
  const actual={...option,grades:grade},begun=clock();
  try{const regularized=seed==='tetra'&&actual.regularize,mesh=(regularized?coneMesh:hemisphereMesh)(graph,seed,R,S,actual),health=regularized?coneHealth(mesh):meshHealth(mesh);
   attempts.push({option:actual,health,milliseconds:clock()-begun});
   if(health.healthy)return{status:'healthy',requested,option:actual,mesh,health,attempts};
  }catch(error){attempts.push({option:actual,error:String(error),milliseconds:clock()-begun});}
 }
 return{status:'inconclusive',message:'No nondegenerate floating-point mesh was found; no geometric root is inferred.',requested,attempts};
}

function plans(config,source){
 const R=Number(source[0]),S=Number(source[1]),n=Math.max(...config.weights),top=config.seed==='double'?.5:1;
 const distance=Math.min(1-R,1-S/top,R+S/top-1);
 // Tetrahedral cones need deep radial grading before angular refinement.
 // Low grades may give healthy meshes whose Hodge guesses lie in a different
 // analytic basin; two high-grade mesh sizes give a genuinely local gate.
 if(config.seed==='tetra'){
  // A high-order cone needs its own local coordinate. This keeps the same
  // physical gluing and integral marking, while avoiding deep globally
  // rounded triangles. Two independently refined meshes remain required.
  if(n>12)return distance>.15?
   [{sub:6,rings:4,grades:8,regularize:true},{sub:8,rings:5,grades:12,regularize:true},{sub:12,rings:7,grades:16,regularize:true},{sub:16,rings:9,grades:20,regularize:true}]:
   [{sub:8,rings:5,grades:12,regularize:true},{sub:12,rings:7,grades:16,regularize:true},{sub:16,rings:9,grades:20,regularize:true},{sub:20,rings:11,grades:24,regularize:true}];
  // These moderate cone orders are calibrated against the former deeply
  // graded selector and a third independent physical mesh. Local cone
  // coordinates need fewer radial layers; all health, trust-neighbourhood,
  // marking and two-mesh consensus checks below remain unchanged.
  if(n>=4&&n<=8)return[
   {sub:4,rings:3,grades:6,regularize:true},
   {sub:6,rings:4,grades:8,regularize:true},
   {sub:8,rings:5,grades:12,regularize:true},
   {sub:12,rings:7,grades:16,regularize:true}];
  const grades=Math.min(48,6*n+12);
  return distance>.15?[{sub:8,rings:5,grades},{sub:10,rings:7,grades},{sub:12,rings:8,grades},{sub:16,rings:10,grades}]:
   [{sub:10,rings:7,grades},{sub:12,rings:8,grades},{sub:16,rings:10,grades},{sub:20,rings:12,grades}];
 }
 if(distance>.15)return[{sub:8,rings:5,grades:Math.max(10,Math.min(18,2*n+2))},
  {sub:12,rings:7,grades:Math.max(14,Math.min(24,3*n+2))},
  {sub:16,rings:10,grades:Math.min(48,6*n+12)},{sub:20,rings:12,grades:Math.min(48,6*n+12)}];
 return[{sub:10,rings:6,grades:Math.max(16,Math.min(32,3*n))},
  {sub:12,rings:8,grades:Math.max(24,Math.min(40,4*n+4))},
  {sub:16,rings:10,grades:Math.min(48,6*n+12)},
  {sub:20,rings:12,grades:Math.min(48,6*n+12)},
  {sub:24,rings:14,grades:Math.min(48,6*n+12)}];
}

function cacheKey(config,source,chart){
 return JSON.stringify([GEOMETRY_VERSION,config.id,source,config.weights,config.shift,
  config.graph.vertices,config.graph.pairs,config.graph.length_coefficients,config.graph.face_representatives]);
}
function store(key,value){cache.delete(key);cache.set(key,value);while(cache.size>24)cache.delete(cache.keys().next().value);}
export function clearGeometryCache(){cache.clear();}
export function geometryCacheStats(){return{version:GEOMETRY_VERSION,entries:cache.size,limit:24};}

function comparison(result,guess,chart){
 if(!guess)return{...result,status:'accepted',accepted:true,comparison:null};
 const E=engine(90),a=E.mobius(E.c(result.tau),chart),b=E.mobius(E.c(guess),chart),error=E.abs(E.sub(a,b)),limit=E.abs(a).plus(1).times('1e-8'),agrees=error.lte(limit);
 return{...result,status:agrees?'accepted':'disagreement',accepted:agrees,comparison:{agrees,chartError:error.toString(),tolerance:limit.toString(),
  explanation:agrees?'The independently selected geometric root agrees with the continued root.':'Two independent geometric refinements agree with each other and disagree with the continued root.'}};
}

/**
 * Independent marked-graph root selection. `guess` is only compared AFTER a
 * consensus is obtained; it is never passed to a geometric initializer.
 * `correct` may be sync or async and has the corrector.mjs structured contract.
 */
export async function geometryAnchor({config,source,chart=identity,guess,deadline=clock()+30000,correct,meshPlans}={}){
 const begun=clock();source=sourceArray(source);chart=chart.map(Number);
 if(!config?.graph||!['double','tetra'].includes(config.seed))return{status:'inconclusive',accepted:false,checked:false,branchRejected:false,retryable:false,verificationOutcome:'inconclusive',message:'A full labelled double-edge or tetrahedral graph is required.',attempts:[]};
 const [R,S]=source.map(Number),top=config.seed==='double'?.5:1,distance=Math.min(R,1-R,S,top-S,R+S/top-1);
 if(!(distance>2e-5))return{status:'inconclusive',accepted:false,checked:false,branchRejected:false,retryable:true,verificationOutcome:'inconclusive',message:'The source point is too close to the boundary for the finite-difference Hodge marking. Use an interior geometric gate.',attempts:[]};
 const key=cacheKey(config,source,chart),cached=cache.get(key);
 if(cached){cache.delete(key);cache.set(key,cached);return comparison({...cached,chart,cacheHit:true,elapsedMs:clock()-begun},guess,chart);}
 if(!correct)correct=(await import('./double-corrector.mjs?v=20261024-tv-audit2')).correct;
 const highCone=config.seed==='tetra'&&Math.max(...config.weights)>12;
 const attempts=[],candidates=[],seen=new Set(),precision=highCone?Math.max(70,50+2*Math.max(...config.weights)):Math.max(70,Math.ceil((40+4*Math.max(...config.weights))/10)*10);
 for(const requested of meshPlans??plans(config,source)){
  if(clock()>deadline)break;
  const checked=checkedHemisphereMesh(config.graph,config.seed,R,S,requested,{deadline});
  const row={requested,healthAttempts:checked.attempts};attempts.push(row);
  if(checked.status!=='healthy'){row.status=checked.failureKind==='time_budget'?'time-budget':'invalid-mesh';row.message=checked.message;if(checked.failureKind==='time_budget')break;continue;}
  if(clock()>deadline){row.status='time-budget';row.message='Geometry time budget reached after mesh health checks.';break;}
  row.option=checked.option;row.health=checked.health;
  const meshId=JSON.stringify(checked.option);if(seen.has(meshId)){row.status='duplicate-healthy-mesh';continue;}seen.add(meshId);
  try{
   const coneCentered=config.seed==='tetra'&&!!checked.option.regularize;
   const h=(coneCentered?coneHodge:markedHodge)(config.graph,config.seed,R,S,{...checked.option,mesh:checked.mesh});
   row.hodge={tau:h.tau,vertices:h.vertices,iterations:h.iterations,markError:h.markError,milliseconds:h.milliseconds};
   if(clock()>deadline){row.status='time-budget';break;}
   // Raw-coordinate trust neighbourhood prevents an opened chart from making
   // a physically distant root appear admissible. The requested chart is used
   // only later to compare independently selected roots with the trace.
   const radius=Math.min(.065,.035*(1+abs(h.tau))),maxStep=.004*(1+abs(h.tau));
   const got=await correct({config,source,guess:h.tau.map(String),chart:identity,precision,digits:coneCentered&&highCone?22:30,radius,maxStep,deadline});
   row.refinement={tau:got.tau,converged:got.converged,precision:got.precision,iterations:got.iterations,evaluations:got.evaluations,
    correction:got.correction,estimatedError:got.estimatedError,retryable:got.retryable,failureKind:got.failureKind,message:got.message,radius,maxStep};
   if(!got.converged||!got.tau){row.status=got.failureKind==='time_budget'?'time-budget':'refinement-inconclusive';if(got.failureKind==='time_budget')break;continue;}
   const rawShift=difference(got.tau,h.tau);
   if(rawShift>radius*(1+1e-8)){row.status='outside-hodge-neighbourhood';continue;}
   row.status='candidate';row.hodgeDisplacement=rawShift;
   const candidate={tau:got.tau.map(String),hodge:h.tau,option:checked.option,precision:got.precision,attempt:attempts.length-1};
   candidates.push(candidate);
   for(const other of candidates.slice(0,-1)){
    const exactError=engine(90).abs(engine(90).sub(engine(90).c(candidate.tau),engine(90).c(other.tau))),rootError=exactError.toNumber(),threshold=1e-14*(1+abs(candidate.tau));
    const rawMeshDifference=difference(candidate.hodge,other.hodge),rawScale=1+abs(candidate.hodge);
    // Consensus requires two different healthy meshes, bounded correction on
    // each, and raw periods in the same local geometric neighbourhood.
    if(rootError>threshold||rawMeshDifference>.08*rawScale)continue;
    const result={status:'accepted',accepted:true,checked:true,certified:false,tau:candidate.tau,
    precision:Math.min(candidate.precision??precision,other.precision??precision),source,chart,
    coneCentered:!!checked.option.regularize,formulaInputsUsed:false,
     geometryVersion:GEOMETRY_VERSION,cacheHit:false,consensus:{meshOptions:[other.option,candidate.option],
      rootDifference:rootError,rootDifferenceExact:exactError.toString(),threshold,rawMeshDifference,independentInitializers:true},
     attempts,elapsedMs:clock()-begun,message:'Two independent healthy marked Hodge meshes refine within their neighbourhoods to the same analytic root.'};
    store(key,result);return comparison(result,guess,chart);
   }
  }catch(error){row.status=isDeadlineError(error)?'time-budget':'geometric-inconclusive';row.message=String(error);if(isDeadlineError(error))break;}
 }
 const timeBudgetReached=clock()>deadline||attempts.some(a=>a.status==='time-budget'||a.refinement?.failureKind==='time_budget');
 return{status:'inconclusive',accepted:false,checked:false,certified:false,source,chart,geometryVersion:GEOMETRY_VERSION,
  failureKind:timeBudgetReached?'time_budget':'geometric_inconclusive',stopKind:timeBudgetReached?'time_budget':'verification',timeBudgetReached,retryable:true,branchRejected:false,verificationOutcome:'inconclusive',
  candidates:candidates.map(({tau,option,precision})=>({tau,option,precision})),attempts,elapsedMs:clock()-begun,
  message:timeBudgetReached?'Geometry computation time budget reached before independent local consensus. Resume with the saved numerical prefix.':'The geometric selector did not obtain independent local consensus. This is not evidence that the continued branch is wrong.'};
}

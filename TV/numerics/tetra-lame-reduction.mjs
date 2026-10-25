// Marked multiplication-by-two reduction for opposite edge counts.
// Branch anchors are actual classical geometric anchors, never cusp guesses.
import copiedAnchors from './tetra-lame-anchors.mjs?v=20261024-tetra1';
import {checkCanonicalTetra} from './tetra-canonical.mjs?v=20261024-tetra1';
export {canonicalGraph as canonicalGraphForCounts,checkCanonicalTetra as validateCanonicalGraph} from './tetra-canonical.mjs?v=20261024-tetra1';
import {correct as boundedCorrect, getEvaluator} from './double-corrector.mjs?v=20261024-tv-audit2';

export const REDUCTION_VERSION='tetra-lame-opposite-pairs-v2';
const keys=['qAB','qAC','qAD','qBC','qBD','qCD'];
const identity=[1,0,0,1],clock=()=>performance.now(),cache=new Map(),LIMIT=24;
const asSource=x=>Array.isArray(x)?x.map(String):[String(x.R),String(x.S)];
const failure=(message,extra={})=>({status:'inconclusive',accepted:false,checked:false,converged:false,retryable:false,message,...extra});
function decimalParts(value){
 const match=String(value).trim().match(/^([+-]?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?$/i);
 if(!match||!(match[2]||match[3]))throw Error('A finite decimal source coordinate is required.');
 const coefficient=BigInt((match[1]==='-'?'-':'')+(match[2]||'0')+(match[3]||''));
 return{coefficient,exponent:Number(match[4]||0)-(match[3]||'').length};
}
/** Exact finite-decimal division; never rounds a tiny boundary distance. */
export function halfDecimal(value){let{coefficient,exponent}=decimalParts(value);if(coefficient%2n===0n)coefficient/=2n;else{coefficient*=5n;exponent--;}return coefficient.toString()+'e'+exponent;}
function canonicalDecimal(value){let{coefficient,exponent}=decimalParts(value);if(coefficient===0n)return'0';while(coefficient%10n===0n){coefficient/=10n;exponent++;}return coefficient.toString()+'e'+exponent;}
function addInteger(value,k){let{coefficient,exponent}=decimalParts(value);if(exponent>=0)return(coefficient*10n**BigInt(exponent)+BigInt(k)).toString();return(coefficient+BigInt(k)*10n**BigInt(-exponent)).toString()+'e'+exponent;}
/** Eligibility uses the labelled graph counts, never only the TV weights. */
export function eligibility(config){
 if(config?.seed!=='tetra'||!config.counts)return{eligible:false,reason:'A full tetrahedral configuration with six labelled counts is required.'};
 const values=keys.map(k=>config.counts[k]);
 if(!values.every(v=>Number.isSafeInteger(v)&&v>=0))return{eligible:false,reason:'Every labelled tetrahedral count must be a nonnegative integer.'};
 const [a,b,c]=values;
 if(values[5]!==a||values[4]!==b||values[3]!==c)return{eligible:false,reason:'Every pair of opposite labelled edges must have equal counts; equal weights alone are insufficient.'};
 const n=1+a+b+c,anchorKey=[a,b,c].join(','),q=a===b&&b===c?a:null;
 if(n>21||!copiedAnchors.anchors[anchorKey])return{eligible:false,q,reason:'The validated classical geometric anchor table covers indices n <= 21.'};
 if(!Array.isArray(config.weights)||config.weights.length!==4||config.weights.some(v=>v!==n))return{eligible:false,reason:'The labelled counts and ordered TV weights disagree.'};
 if(!Array.isArray(config.shift)||config.shift.length!==2||config.shift.some(v=>v!==0))return{eligible:false,reason:'The opposite-pair family requires the canonical ordinary TV shifts to equal [0,0] exactly.'};
 const canonical=checkCanonicalTetra(config);if(!canonical.valid)return{eligible:false,reason:canonical.reason,field:canonical.field};
 return{eligible:true,q,n,anchorKey,keys:[...keys],counts:[a,b,c],weights:[n,n,n,n],lift:[(a+c)%2,(b+c)%2],version:REDUCTION_VERSION,
  identity:'Sum over the four half-period translates of wp(z) equals 4 wp(2z).',
  marking:'The modulus tau and homology basis are unchanged. Both Lamé Bloch multipliers square under u=2z.',
  graphIdentification:'Opposite labelled tetrahedral counts (a,b,c); corresponding marked classical decoration (a,b,c).'};
}
export const eligible=config=>eligibility(config).eligible;
export function character(config,R,S){
 const info=eligibility(config);if(!info.eligible)return info;
 const [a,b,c]=info.counts,r=halfDecimal(addInteger(S,a+c)),s=halfDecimal(addInteger(R,b+c));
 return{...info,r,s,analyticTV:[String(S),String(R)],analyticLame:[r,s],geometricClassical:[halfDecimal(S),halfDecimal(R)],accessoryScale:4};
}
function reduced(config){const info=eligibility(config);return{...config,id:'lame-reduction-'+info.anchorKey,weights:[info.n,0,0,0],shift:[info.lift[1]/2,info.lift[0]/2],connector:undefined};}
function reducedSource(source){return asSource(source).map(halfDecimal);}

/** Same bounded-corrector contract as corrector.mjs. `guess` is raw tau.
 * All corrections, trust radii, and movement checks use the supplied opened
 * chart. Precision agreement validates arithmetic, not branch identity.
 */
export function factorCorrect(options={}){
 const info=eligibility(options.config);
 if(!info.eligible)return failure(info.reason,{method:REDUCTION_VERSION,eligible:false,failureKind:'ineligible',tau:options.guess});
 const got=boundedCorrect({...options,config:reduced(options.config),source:reducedSource(options.source)});
 return{...got,method:REDUCTION_VERSION,reduction:info,
  diagnostics:{...got.diagnostics,reducedDimension:Math.ceil(info.n/2),fullDimension:4*Math.ceil(info.n/2),
   characterDictionary:character(options.config,...asSource(options.source)),conjecturalAsymptoticDataUsed:false}};
}
export const correct=factorCorrect;

/** Independent, unfactored full-TV solve near an already selected root. */
export function verifyTV(options={}){
 const info=eligibility(options.config);if(!info.eligible)return failure(info.reason,{failureKind:'ineligible'});
 const E=getEvaluator(Math.max(options.precision??80,options.digits??30)),selected=E.c(options.guess),chart=options.chart??identity;
 // An independent solve must actually move: a residual check at the same
 // rounded input can otherwise return that input at both precisions.
 const w=E.mobius(selected,chart),probe=options.probe??['1e-7','2e-7'],initial=E.strs(E.inverse(E.add(w,E.c(probe)),chart));
 const config={...options.config,connector:undefined},got=boundedCorrect({...options,guess:initial,config,fast:false,
  radius:options.radius??'.005',maxStep:options.maxStep??'.001'});
 if(got.converged){const F=getEvaluator(Math.max(got.precision,80));got.selectedRootDifference=F.abs(F.sub(F.c(got.tau),F.c(options.guess))).toString();}
 return{...got,method:'independent unfactored full-TV correction from a perturbed start',reduction:info,probe};
}

function cacheKey(config,source){const i=eligibility(config);return JSON.stringify([REDUCTION_VERSION,i.anchorKey,asSource(source).map(canonicalDecimal)]);}
function save(key,value){cache.delete(key);cache.set(key,value);while(cache.size>LIMIT)cache.delete(cache.keys().next().value);}
export const clearCache=()=>cache.clear();
export const cacheStats=()=>({entries:cache.size,limit:LIMIT,version:REDUCTION_VERSION});

function compare(result,guess,chart){
 if(!guess)return result;
 const E=getEvaluator(100),a=E.mobius(E.c(result.tau),chart),b=E.mobius(E.c(guess),chart),d=E.abs(E.sub(a,b)),t=E.abs(a).plus(1).times('1e-8'),agrees=d.lte(t);
 return{...result,status:agrees?'accepted':'disagreement',accepted:agrees,comparison:{agrees,chartError:d.toString(),tolerance:t.toString(),explanation:agrees?'The independently marked classical Lamé branch agrees.':'The independently transported marked classical Lamé branch disagrees.'}};
}

/** Actual marked classical anchor, followed by bounded source continuation.
 * guess is used only in the final comparison, never to select this branch.
 * The fixed-table seeds and transported points are cached separately.
 */
export async function geometryAnchor({config,source,chart=identity,guess,deadline=clock()+60000,verifyFullTV=true}={}){
 const begun=clock(),info=eligibility(config);source=asSource(source??['.58','.62']);
 if(!info.eligible)return failure(info.reason,{failureKind:'ineligible',elapsedMs:clock()-begun});
 const [R,S]=source.map(Number);if(!(R>0&&R<1&&S>0&&S<1&&R+S>1))return failure('An interior tetrahedral source point is required.',{failureKind:'invalid_source'});
 const reusable=value=>value&&(!verifyFullTV||value.consensus?.fullTV),key=cacheKey(config,source),hit=reusable(cache.get(key))?cache.get(key):null;
 if(hit)return compare({...hit,chart,cacheHit:true,elapsedMs:clock()-begun},guess,chart);
 const initialSource=[copiedAnchors.sourceR,copiedAnchors.sourceS],seedKey=cacheKey(config,initialSource),precision=Math.max(80,Math.ceil((55+4*info.n)/10)*10),digits=32;
 let current=reusable(cache.get(seedKey))?cache.get(seedKey):null,attempts=[];
 if(!current){
  const seed=copiedAnchors.anchors[info.anchorKey].highPrecisionTau;
  const low=factorCorrect({config,source:initialSource,guess:seed,chart:identity,precision,digits,radius:'.015',maxStep:'.002',deadline,fast:false});
  attempts.push({phase:'classical-anchor-refine',...low});
  if(!low.converged)return failure('The copied geometric anchor did not refine in its bounded neighbourhood.',{attempts,failureKind:low.failureKind,retryable:low.retryable,elapsedMs:clock()-begun});
  const high=factorCorrect({config,source:initialSource,guess:low.tau,chart:identity,precision:precision+30,digits:45,radius:'.003',maxStep:'.0005',deadline,fast:false});
  attempts.push({phase:'classical-anchor-independent-precision',...high});
  const E=getEvaluator(precision+30),error=high.converged?E.abs(E.sub(E.c(low.tau),E.c(high.tau))):new E.D(Infinity);
  if(!high.converged||error.gt('1e-25'))return failure('The classical anchor did not pass the repeated precision check.',{attempts,failureKind:high.failureKind??'precision',retryable:high.retryable,elapsedMs:clock()-begun});
  let full=null;
  if(verifyFullTV){
   full=verifyTV({config,source:initialSource,guess:high.tau,chart:identity,precision:precision+20,digits:35,deadline});
   attempts.push({phase:'unfactored-TV-anchor-check',...full});
   if(!full.converged||new E.D(full.selectedRootDifference??Infinity).gt('1e-25'))return failure('The marked reduced anchor did not pass its full-TV check.',{attempts,failureKind:full.failureKind??'verification',retryable:full.retryable,elapsedMs:clock()-begun});
  }
  current={status:'accepted',accepted:true,checked:true,certified:false,tau:high.tau,source:initialSource,chart:identity,
   precision:high.precision,method:REDUCTION_VERSION,geometryVersion:REDUCTION_VERSION,reduction:info,
   consensus:{independentInitializers:false,exactPullback:true,classicalGeometricAnchor:true,rootDifferenceExact:error.toString(),fullTV:full?{precision:full.precision,difference:full.selectedRootDifference}:null},
   sourceAnchor:{classicalCounts:[...info.counts],classicalSource:copiedAnchors.classicalSource,sourceHash:copiedAnchors.source_sha256},
   attempts:[...attempts],message:'The selected classical geometric anchor passes the exact pullback and independent precision checks.'};
  save(seedKey,current);
 }
 // Prefer the nearest already checked point with the same labelled q. A
 // caller's continued root is deliberately excluded from this selection.
 const E=getEvaluator(Math.max(precision,100));
 let start=current;
 for(const value of cache.values())if(value.reduction.anchorKey===info.anchorKey&&value.checked){
  const d=value.source.reduce((sum,x,j)=>sum+(Number(x)-Number(source[j]))**2,0),old=start.source.reduce((sum,x,j)=>sum+(Number(x)-Number(source[j]))**2,0);
  if(d<old)start=value;
 }
 let tau=start.tau,position=start.source,steps=0,rejections=0,fraction=0;
 const distance=Math.hypot(...source.map((x,j)=>Number(x)-Number(position[j]))),origin=position;
 let increment=distance?Math.min(1,.018/(1+.04*(info.n-1))/distance):1;
 while(fraction<1&&distance>0){
  if(clock()>deadline)return failure('Marked classical continuation time budget reached.',{failureKind:'time_budget',retryable:true,attempts,tau,source:position,elapsedMs:clock()-begun});
  const next=Math.min(1,fraction+increment),target=source.map((x,j)=>new E.D(origin[j]).plus(new E.D(x).minus(origin[j]).times(String(next))).toString());
  const got=factorCorrect({config,source:target,guess:tau,chart:identity,precision,digits,radius:'.035',maxStep:'.004',deadline,fast:true});
  const movement=got.converged?E.abs(E.sub(E.c(got.tau),E.c(tau))).toNumber():Infinity;
  if(!got.converged||movement>.015){rejections++;increment/=2;if(increment<1e-7)return failure('No sufficiently small marked Lamé continuation step was found.',{failureKind:got.failureKind??'local_convergence',retryable:true,attempts,tau,source:position,steps,rejections,elapsedMs:clock()-begun});continue;}
  tau=got.tau;position=target;fraction=next;steps++;increment=Math.min(1-fraction,increment*(movement<.006?1.25:1));
 }
 if(distance>0){
  const high=factorCorrect({config,source,guess:tau,chart:identity,precision:precision+30,digits:45,radius:'.003',maxStep:'.0005',deadline,fast:false});
  attempts.push({phase:'transported-final-precision',...high});
  const error=high.converged?E.abs(E.sub(E.c(tau),E.c(high.tau))):new E.D(Infinity);
  if(!high.converged||error.gt('1e-25'))return failure('The independently transported marked Lamé point did not pass repeated precision.',{attempts,failureKind:high.failureKind??'precision',retryable:high.retryable,elapsedMs:clock()-begun});
  tau=high.tau;
 }
 const result={...current,tau,source,chart,cacheHit:false,attempts,continuation:{steps,rejections,independentStart:start.source},elapsedMs:clock()-begun,
  message:distance?'Independently transported from an actual marked classical Lamé anchor using the exact pullback.':current.message};
 save(key,result);return compare(result,guess,chart);
}
export const anchor=geometryAnchor;

// Research public entry point: infinity cusp for the two repeated-connector
// families. Inputs/outputs are decimal strings. The known marked branch and
// zero isogeny lift are retained. This does not select other graph branches.
import Decimal from '../decimal.mjs';
import{createFourierEvaluator}from'./fourier-kronecker.mjs?v=20261024-tetra-dense4';
const evaluatorCache=new Map();

export function resolveConnectorInfinity({family='double',q=8,delta,S,epsilon,
 startR='.58',startS=family==='double'?'.31':'.62',endS=family==='double'?'.25':'.5',
 digits=20,verifyModes=true,maxPrecision=900}={}){
 if(!['double','tetra'].includes(family)||!Number.isInteger(q)||q<0)throw Error('A connector family and nonnegative integer count are required');
 if(epsilon===undefined&&(delta===undefined||S===undefined))throw Error('Provide epsilon for the marked straight path, or exact decimal delta and S');
 const A=Decimal.clone({precision:maxPrecision}),target=new A(10).pow(-digits),n=q+1;
 const Probe=Decimal.clone({precision:40});
 const depth=epsilon!==undefined?-new Probe(String(epsilon)).logarithm(10).toNumber():-new Probe(String(delta)).logarithm(10).toNumber();
 let P=Math.ceil(Math.max(n<=5?70:90,2*depth+(n<=5?55:100)+2*n)/20)*20,previous=null,high=null,checks=[];
 const begun=Date.now();
 function one(precision,modesMultiplier=1){
  const key=precision+':'+modesMultiplier;if(!evaluatorCache.has(key))evaluatorCache.set(key,createFourierEvaluator(precision,{bigint:true,modesMultiplier}));
  const E=evaluatorCache.get(key),t=Date.now();
  let dd,ss;
  if(epsilon!==undefined){const ee=new E.D(String(epsilon));dd=ee.times(new E.D(1).minus(String(startR)));ss=new E.D(String(endS)).plus(ee.times(new E.D(String(startS)).minus(String(endS))))}
  else{dd=new E.D(String(delta));ss=new E.D(String(S))}
  if(dd.lte(0)||dd.gte('.1'))throw Error('The Fourier cusp entry point requires 0 < delta < .1');
  // Request enough digits to resolve the small horizontal displacement. A
  // tolerance tied only to P-70 makes the first low-cost stage stop too early
  // and needlessly repeats the whole calculation at a third precision.
  const toleranceDigits=Math.min(precision-30,Math.ceil(2*depth+digits+20));
  const result=E.solve(family,n,dd,ss,E.prediction(family,n,dd,ss),{tol:new E.D(10).pow(-toleranceDigits)});
  return{precision,modesMultiplier,w:E.strs(result.w),correction:result.correction.toString(),delta:dd.toString(),S:ss.toString(),evaluations:result.evaluations,milliseconds:Date.now()-t};
 }
 function compare(a,b){
  const dx=new A(a.w[0]).minus(b.w[0]).abs(),dy=new A(a.w[1]).minus(b.w[1]).abs(),den=new A(b.correction).abs();
  const rel=den.isZero()?new A(Infinity):new A(a.correction).minus(b.correction).abs().div(den);
  return{absoluteWError:A.max(dx,dy).toSignificantDigits(12).toString(),relativeHorizontalError:rel.toSignificantDigits(12).toString(),passed:A.max(dx,dy).lt(target)&&rel.lt(target)};
 }
 while(P+80<=maxPrecision){
  try{previous=high??one(P);high=one(P+80)}catch(error){checks.push({kind:'precision retry',precision:P,error:String(error)});P+=80;high=null;continue}
  const check={kind:'precision',...compare(previous,high),lowerPrecision:previous.precision,higherPrecision:high.precision};checks.push(check);
  if(check.passed){
   if(verifyModes){const doubled=one(high.precision,2),modeCheck={kind:'Fourier modes doubled',...compare(high,doubled),precision:high.precision};checks.push(modeCheck);if(!modeCheck.passed){P+=80;high=doubled;continue}high=doubled}
   return{status:'checked',family,q,n,digits,source:epsilon!==undefined?'marked straight path':'fixed S',epsilon:epsilon===undefined?undefined:String(epsilon),...high,checks,totalMilliseconds:Date.now()-begun};
  }
  P+=80;
 }
 return{status:'precision limit',family,q,n,digits,last:high,checks,totalMilliseconds:Date.now()-begun};
}

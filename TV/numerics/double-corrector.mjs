// Bounded local TV correction; Decimal is the default arithmetic.
import * as F from './generic-fast.mjs';
import {createEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';

const engines=new Map();
export function getEvaluator(precision=80){const p=Math.max(50,Math.min(970,Math.ceil(Number(precision)/10)*10));if(!engines.has(p))engines.set(p,createEvaluator(p));return engines.get(p);}
const now=()=>performance.now();
const fmobius=(z,[a,b,c,d])=>F.div(F.add(F.mul(a,z),b),F.add(F.mul(c,z),d));
const finverse=(z,[a,b,c,d])=>F.div(F.sub(F.mul(d,z),b),F.sub(a,F.mul(c,z)));
const opts=config=>({isogeny:!!config.connector,isogenyLift:config.connector?.isogenyLift??[0,0]});

function detFast(M,scales){
 let product=[1,0];const N=M.length;
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=F.mul(1/scales[i],M[i][j]);
 for(let j=0;j<N;j++){
  let pivot=j;for(let i=j+1;i<N;i++)if(F.abs(M[i][j])>F.abs(M[pivot][j]))pivot=i;
  if(!F.abs(M[pivot][j]))return[0,0];if(pivot!==j){[M[pivot],M[j]]=[M[j],M[pivot]];product=F.neg(product);}
  product=F.mul(product,M[j][j]);for(let i=j+1;i<N;i++){const ratio=F.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=F.sub(M[i][k],F.mul(ratio,M[j][k]));}
 }
 return product;
}
function detDecimal(E,M,scales){
 let product=E.c(1);const N=M.length;
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=E.div(M[i][j],scales[i]);
 for(let j=0;j<N;j++){
  let pivot=j;for(let i=j+1;i<N;i++)if(E.norm2(M[i][j]).gt(E.norm2(M[pivot][j])))pivot=i;
  if(E.norm2(M[pivot][j]).isZero())return E.c(0);if(pivot!==j){[M[pivot],M[j]]=[M[j],M[pivot]];product=E.neg(product);}
  product=E.mul(product,M[j][j]);for(let i=j+1;i<N;i++){const ratio=E.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=E.sub(M[i][k],E.mul(ratio,M[j][k]));}
 }
 return product;
}
function frozenFast(config,source,initial,chart,counter,deadline){
 const r=Number(source[1])+config.shift[1],s=Number(source[0])+config.shift[0];
 if(config.connector)return w=>{if(now()>deadline)throw Error('Local corrector deadline reached');counter.fast++;return F.determinant(config.weights,finverse(w,chart),r,s,opts(config)).value;};
 const red=F.reduceModular(config.weights,finverse(initial,chart),r,s),matrix=red.matrix;
 const first=F.tvMatrix(red.weights,red.tau,red.r,red.s,{reduce:false}).M;
 const scales=first.map(row=>Math.max(...row.map(F.abs)));if(scales.some(x=>!x||!Number.isFinite(x)))throw Error('Fast frozen row scaling is singular');
 let firstValue=detFast(first,scales);counter.fast++;
 const f=w=>{if(now()>deadline)throw Error('Local corrector deadline reached');if(firstValue&&F.abs(F.sub(w,initial))===0){const value=firstValue;firstValue=null;return value;}counter.fast++;return detFast(F.tvMatrix(red.weights,fmobius(finverse(w,chart),matrix),red.r,red.s,{reduce:false}).M,scales);};
 Object.assign(f,{frame:matrix,scales,weights:red.weights,characters:[red.r,red.s]});return f;
}
function frozenDecimal(E,config,source,initial,chart,counter,deadline){
 const r=new E.D(String(source[1])).plus(config.shift[1]),s=new E.D(String(source[0])).plus(config.shift[0]);
 if(config.connector)return w=>{if(now()>deadline)throw Error('Local corrector deadline reached');counter.decimal++;return E.determinant(config.weights,E.inverse(w,chart),r,s,opts(config)).value;};
 const red=E.reduceModular(config.weights,E.inverse(initial,chart),r,s),matrix=red.matrix;
 const first=E.tvMatrix(red.weights,red.tau,red.r,red.s,{reduce:false}).M;
 const scales=E.rowScales?E.rowScales(first):first.map(row=>E.D.max(...row.map(E.abs)));if(scales.some(x=>x.isZero()))throw Error('Frozen row scaling is singular');
 let firstValue=detDecimal(E,first,scales);counter.decimal++;
 const f=w=>{if(now()>deadline)throw Error('Local corrector deadline reached');if(firstValue&&E.abs(E.sub(w,initial)).isZero()){const value=firstValue;firstValue=null;return value;}counter.decimal++;return detDecimal(E,E.tvMatrix(red.weights,E.mobius(E.inverse(w,chart),matrix),red.r,red.s,{reduce:false}).M,scales);};
 Object.assign(f,{frame:matrix,scales,weights:red.weights,characters:[red.r,red.s]});return f;
}
function fastPolish(config,source,initial,chart,{radius,maxStep,deadline,counter,maxIterations=16,linearize=false}){
 const f=frozenFast(config,source,initial,chart,counter,deadline),scale=1+F.abs(initial),h=[1e-6*scale,2e-6*scale];let z=initial,value=f(z),derivative=F.div(F.sub(f(F.add(z,h)),value),h),iterations=0,best=initial;
 const history=[],finish=(w,iterations,converged)=>{let linearization=null;if(linearize&&f.frame){const h=1e-5*(1+F.abs(w)),derivative=F.div(F.sub(f(F.add(w,h)),f(F.sub(w,h))),2*h);if(Number.isFinite(F.abs(derivative))&&F.abs(derivative)>0)linearization={frame:f.frame,scales:f.scales,weights:f.weights,characters:f.characters,derivative};}return{w,iterations,converged,history,linearization};};
 for(;iterations<maxIterations;iterations++){
  let step=F.div(value,derivative),length=F.abs(step);if(!Number.isFinite(length))break;
  if(length<2e-10)return finish(z,iterations,true);if(length>maxStep)step=F.mul(maxStep/length,step);
  let accepted=false,factor=1,next,nextValue;
  for(let k=0;k<9;k++){
   next=F.sub(z,F.mul(factor,step));if(next[1]>.002&&F.abs(F.sub(next,initial))<=radius){nextValue=f(next);if(Number.isFinite(F.abs(nextValue))&&F.abs(nextValue)<F.abs(value)){accepted=true;break;}}factor/=2;
  }
  history.push({step:length,factor,accepted});if(!accepted)break;
  derivative=F.div(F.sub(nextValue,value),F.sub(next,z));z=next;value=nextValue;best=z;
 }
 return finish(best,iterations,false);
}

/** Local correction only: a graph-selected anchor or checked predecessor is required.
 * radius/maxStep are absolute distances in the supplied opened chart.
 * estimatedError is the last Newton/secant step, never an interval bound.
 * An optional arithmetic evaluator changes number representation only; the
 * same neighbourhood, convergence target and residual checks are applied.
 */
export function correct({config,source,guess,chart=[1,0,0,1],precision=80,digits=25,radius,maxStep,deadline=Infinity,fast=true,maxIterations=28,verifyResidual=true,reuseFastDerivative=false,arithmetic=null}={}){
 const begun=now(),E=arithmetic??getEvaluator(Math.max(precision,digits+30)),initial=E.mobius(E.c(guess),chart),initialSize=E.abs(initial).plus(1),R=new E.D(radius??initialSize.times('.075')),H=new E.D(maxStep??initialSize.times('.02')),tolerance=new E.D(10).pow(-digits),counter={fast:0,decimal:0,verification:0};
 const diagnostics={engine:'frozen modular basis and row scales; bounded secant with optional double prepolish',fast:null,history:[]};let current=initial,lastStep=null,residual=null,iterations=0,rebuilds=0;
 const finish=(converged,message)=>({tau:E.strs(E.inverse(current,chart)),w:E.strs(current),converged,iterations,evaluations:counter.fast+counter.decimal+counter.verification,evaluationCounts:{...counter},correction:lastStep?E.abs(lastStep).toString():null,estimatedError:lastStep?E.abs(lastStep).toString():null,precision:E.precision,elapsedMs:now()-begun,retryable:!converged,failureKind:converged?null:/deadline/i.test(message)?'time_budget':/precision|singular|derivative/i.test(message)?'precision':'local_convergence',message,residual:residual?.toString()??null,diagnostics});
 try{
  if(!(initial[1].gt('.002')&&R.gt(0)&&H.gt(0)))throw Error('Invalid local correction neighbourhood');
  if(now()>deadline)throw Error('Local corrector deadline reached');
  if(fast){try{const got=fastPolish(config,source,initial.map(Number),chart,{radius:R.toNumber(),maxStep:H.toNumber(),deadline,counter,linearize:reuseFastDerivative});diagnostics.fast=got;const candidate=E.c(got.w.map(String));if(candidate[1].gt('.002')&&E.abs(E.sub(candidate,initial)).lte(R))current=candidate;}catch(error){diagnostics.fast={message:String(error)};}}
  const start=current,f=frozenDecimal(E,config,source,start,chart,counter,deadline);let value=f(current),derivative,usedFastDerivative=false;
  const lin=diagnostics.fast?.linearization;
  if(reuseFastDerivative&&lin&&f.frame&&JSON.stringify(lin.frame)===JSON.stringify(f.frame)&&JSON.stringify(lin.weights)===JSON.stringify(f.weights)&&lin.characters.every((x,i)=>Math.abs(x-Number(f.characters[i]))<1e-8)){
   let factor=new E.D(1);for(let j=0;j<lin.scales.length;j++)factor=factor.times(String(lin.scales[j])).div(f.scales[j]);derivative=E.mul(E.c(lin.derivative.map(String)),factor);usedFastDerivative=true;diagnostics.fastDerivativeReused=true;
  }else{const perturb=E.D.min(new E.D('1e-6').times(E.abs(current).plus(1)),R.div(20)),h=[perturb,perturb.times(2)];derivative=E.div(E.sub(f(E.add(current,h)),value),h);}
  for(;iterations<maxIterations;iterations++){
   if(now()>deadline)throw Error('Local corrector deadline reached');
   if(E.norm2(value).isZero())throw Error('Exact zero at working precision; increase precision');
   if(E.norm2(derivative).isZero())throw Error('Zero local secant derivative');
   let step=E.div(value,derivative),length=E.abs(step);lastStep=step;residual=E.abs(value);
   if(length.lt(tolerance)&&(!usedFastDerivative||iterations>0)){
    if(verifyResidual){counter.verification++;const r=new E.D(String(source[1])).plus(config.shift[1]),s=new E.D(String(source[0])).plus(config.shift[0]),check=E.determinant(config.weights,E.inverse(current,chart),r,s,opts(config));diagnostics.productionResidual=E.abs(check.value).toString();diagnostics.minRelativePivot=check.minRelativePivot?.toString();}
    return finish(true,'Finite selected root corrected inside the local chart neighbourhood.');
   }
   if(length.gt(H))step=E.mul(step,H.div(length));let accepted=false,factor=new E.D(1),next,nextValue;
   for(let k=0;k<12;k++){
    next=E.sub(current,E.mul(step,factor));if(next[1].gt('.002')&&E.abs(E.sub(next,initial)).lte(R)){nextValue=f(next);if(E.abs(nextValue).lt(E.abs(value))){accepted=true;break;}}factor=factor.div(2);
   }
   diagnostics.history.push({iteration:iterations,step:length.toString(),factor:factor.toString(),accepted});
   if(!accepted){
    if(rebuilds++>=2)throw Error('No decreasing bounded local step');
    const a=new E.D(10).pow(-Math.min(digits+5,Math.floor(E.precision/2))),h2=[a,a.times(2)];derivative=E.div(E.sub(f(E.add(current,h2)),value),h2);continue;
   }
   derivative=E.div(E.sub(nextValue,value),E.sub(next,current));current=next;value=nextValue;
  }
  return finish(false,'Local corrector iteration budget reached.');
 }catch(error){diagnostics.failure=String(error);return finish(false,String(error));}
}

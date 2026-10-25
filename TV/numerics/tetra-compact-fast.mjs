// Binary64 correction for INTERNAL compact-source transport only.
// The caller must retain bounded source steps, reverse checks, Decimal
// fallback, and perturbed higher-precision verification of displayed points.
import * as F from './generic-fast.mjs';
import * as Reduction from './tetra-lame-reduction.mjs?v=20261024-tv-audit2';
import {checkCanonicalTetra} from './tetra-canonical.mjs?v=20261024-tetra1';

const now=()=>performance.now(),identity=[1,0,0,1];
const finite=z=>Array.isArray(z)&&z.length===2&&z.every(Number.isFinite);
const mobius=(z,[a,b,c,d])=>F.div(F.add(F.mul(a,z),b),F.add(F.mul(c,z),d));


function model(config,source){
 const check=checkCanonicalTetra(config);if(!check.valid)throw Error('The full labelled tetrahedral graph is invalid: '+check.reason);
 if(config.weights.reduce((a,b)=>a+b,0)>144||Math.max(...config.weights)>72)throw Error('The compact graph exceeds the existing local numerical weight budget');
 const reduction=Reduction.eligibility(config);
 if(reduction.eligible)return{weights:[reduction.n,0,0,0],
  r:Number(Reduction.halfDecimal(source[1]))+reduction.lift[0]/2,
  s:Number(Reduction.halfDecimal(source[0]))+reduction.lift[1]/2,
  basis:identity,reduction,method:'frozen binary64 marked Lamé factor'};
 const r=Number(source[1])+config.shift[1],s=Number(source[0])+config.shift[0];
 // Preserve the existing exact connector factor, then freeze its modular
 // frame and row scales. No dynamically normalized F.solve call is used.
 if(config.connector&&config.weights[0]===config.weights[1]&&config.weights[2]===config.weights[3]){
  const [u,v]=config.connector.isogenyLift??[0,0];
  if(config.weights[0]===config.weights[2])return{weights:[config.weights[0],0,0,0],r:(r+u)/2,s:(s+v)/2,basis:identity,reduction:null,method:'frozen binary64 connector Lamé factor'};
  if(config.weights[2]===0)return{weights:[config.weights[0],0,0,0],r,s:(s+v)/2,basis:[2,0,0,1],reduction:null,method:'frozen binary64 connector factor'};
  return{weights:[config.weights[0],config.weights[2],0,0],r:(s+v)/2,s:-r,basis:[0,-1,2,0],reduction:null,method:'frozen binary64 connector factor'};
 }
 return{weights:config.weights,r,s,basis:identity,reduction:null,method:'frozen binary64 full-TV determinant'};
}

function frozenFast(data,initial,counter,deadline){
 const red=F.reduceModular(data.weights,mobius(initial,data.basis),data.r,data.s),frame=red.matrix;
 const get=tau=>{
  if(now()>deadline)throw Error('Compact internal correction deadline reached');
  counter.evaluations++;
  return F.tvMatrix(red.weights,mobius(mobius(tau,data.basis),frame),red.r,red.s,{reduce:false}).M;
 };
 const raw=get(initial),N=raw.length,rows=Array.from({length:N},(_,i)=>i),cols=[...rows];
 const rs=raw.map(row=>Math.max(...row.map(F.abs))),M=raw.map((row,i)=>row.map(z=>F.mul(1/rs[i],z))),cs=Array.from({length:N},(_,j)=>Math.max(...M.map(row=>F.abs(row[j]))));
 if([...rs,...cs].some(x=>!Number.isFinite(x)||!x))throw Error('Binary64 frozen equilibration is singular');
 const inverse=rs.map(r=>cs.map(c=>(1/r)/c));
 if(inverse.some(row=>row.some(w=>!Number.isFinite(w)||w===0)))throw Error('Binary64 frozen normalizer is unsafe');
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=F.mul(inverse[i][j],raw[i][j]);
 // Freeze a complete-pivot leading minor. The final Schur entry has the
 // same zeros as the local determinant, without multiplying tiny pivots.
 // Row/column scales and pivot order never change during this correction.
 const healthy=z=>finite(z)&&F.abs(z)>1e-13;
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=F.abs(M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const v=F.abs(M[i][k]);if(v>best){p=i;q=k;best=v;}}
  if(!(Number.isFinite(best)&&best>1e-13))throw Error('Binary64 reference minor needs higher precision');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]];}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]];}
  if(!healthy(M[j][j]))throw Error('Binary64 reference pivot is unsafe');
  for(let i=j+1;i<N;i++){const r=F.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=F.sub(M[i][k],F.mul(r,M[j][k]));}
 }
 let first=M[N-1][N-1];
 const f=tau=>{
  if(first&&F.abs(F.sub(tau,initial))===0){const v=first;first=null;return v;}
  const raw=get(tau),A=rows.map(i=>cols.map(j=>F.mul(inverse[i][j],raw[i][j])));
  for(let j=0;j<N-1;j++){
   if(!healthy(A[j][j]))throw Error('Binary64 fixed minor became unsafe; use higher precision');
   for(let i=j+1;i<N;i++){const r=F.div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=F.sub(A[i][k],F.mul(r,A[j][k]));}
  }
  const v=A[N-1][N-1];if(!finite(v))throw Error('Binary64 Schur residual is nonfinite');return v;
 };
 Object.assign(f,{frame,scales:rs,columnScales:cs,weights:red.weights,characters:[red.r,red.s],residual:'fixed complete-pivot Schur'});return f;
}

/** Frozen-frame fastPolish from double-corrector.mjs, without its Decimal
 * polish. This is a local transport helper, never a displayed-point solver.
 */
export function boundedFastCorrect({config,source,guess,radius,maxStep=.02,deadline=Infinity,maxIterations=16}={}){
 const begun=now(),initial=Array.isArray(guess)?guess.map(Number):[],counter={evaluations:0},history=[];
 let current=initial,iterations=0,lastStep=null,residual=null,data=null;
 const finish=(converged,message)=>({tau:current.map(String),w:current.map(String),converged,status:converged?'ok':'fallback',
  iterations,evaluations:counter.evaluations,precision:15,elapsedMs:now()-begun,
  correction:lastStep===null?null:String(lastStep),estimatedError:lastStep===null?null:String(lastStep),
  residual:residual===null?null:String(residual),retryable:!converged,
  failureKind:converged?null:/deadline/i.test(message)?'time_budget':'binary64_fallback',message,
  diagnostics:{method:(data?.method??'frozen binary64 compact correction')+'; fixed complete-pivot Schur residual',internalOnly:true,requiresDisplayedPrecisionCheck:true,
   reduction:data?.reduction??null,history}});
 try{
  radius=Number(radius??Math.min(.1,.06*(1+F.abs(initial))));maxStep=Number(maxStep);
  maxIterations=Math.min(32,Math.max(1,Math.floor(Number(maxIterations)||16)));
  if(now()>deadline)throw Error('Compact internal correction deadline reached');
  if(!finite(initial)||initial[1]<=.002||!(radius>0)||!(maxStep>0)||!Number.isFinite(radius)||!Number.isFinite(maxStep))throw Error('Invalid bounded compact correction neighbourhood');
  if(!Array.isArray(source)||source.length!==2||source.some(x=>!Number.isFinite(Number(x)))||!(Number(source[0])>0&&Number(source[0])<1&&Number(source[1])>0&&Number(source[1])<1&&Number(source[0])+Number(source[1])>1))throw Error('An interior tetrahedral source point is required');
  data=model(config,source);
  const f=frozenFast(data,initial,counter,deadline),scale=1+F.abs(initial),h=[1e-6*scale,2e-6*scale];
  let value=f(current),derivative=F.div(F.sub(f(F.add(current,h)),value),h);
  for(;iterations<maxIterations;iterations++){
   if(now()>deadline)throw Error('Compact internal correction deadline reached');
   if(!finite(value)||!finite(derivative)||!F.abs(derivative))throw Error('The binary64 local derivative is unresolved');
   if(!F.abs(value))throw Error('The binary64 determinant rounded to zero; use Decimal correction');
   let step=F.div(value,derivative),length=F.abs(step);lastStep=length;residual=F.abs(value);
   if(!Number.isFinite(length))throw Error('The binary64 local step is not finite');
   if(length<2e-10)return finish(true,'Internal bounded binary64 correction; displayed points still require Decimal precision checks.');
   if(length>maxStep)step=F.mul(maxStep/length,step);
   let accepted=false,factor=1,next,nextValue;
   for(let k=0;k<9;k++){
    next=F.sub(current,F.mul(factor,step));
    if(finite(next)&&next[1]>.002&&F.abs(F.sub(next,initial))<=radius){
     nextValue=f(next);
     if(finite(nextValue)&&F.abs(nextValue)<F.abs(value)){accepted=true;break;}
    }
    factor/=2;
   }
   history.push({step:length,factor,accepted});
   if(!accepted)return finish(false,'No decreasing bounded binary64 step; use Decimal correction');
   derivative=F.div(F.sub(nextValue,value),F.sub(next,current));current=next;value=nextValue;
  }
  return finish(false,'The bounded binary64 iteration budget was reached; use Decimal correction');
 }catch(error){return finish(false,String(error?.message??error));}
}

// Compact-image correction only. The labelled graph still supplies its anchor;
// cusp continuation and graph/mesh acceptance use their existing solvers.
import {createEvaluator} from './bouquet-fixed.mjs?v=20261024-tv-audit2';
import {createEvaluator as decimalEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';
import * as F from './generic-fast.mjs';

const engines=new Map(),now=()=>performance.now();
const ensure=deadline=>{if(now()>deadline)throw Error('The compact mesh numerical budget was reached');};
function engine(precision,compact,decimal=false){
 // Precision is the total requested working precision. The fixed arithmetic
 // already retains its intrinsic guard; correction targets add their guard
 // once in compactMeshCorrect, rather than again at evaluator construction.
 const p=Math.max(50,Math.min(950,Math.ceil(precision/10)*10)),key=p+':'+compact+':'+decimal;
 if(!engines.has(key))engines.set(key,decimal?decimalEvaluator(p):createEvaluator(p,{jetScale:compact?1:16,axisShortcuts:true,cacheMatrixConstants:true,compactFourier:compact,cacheFourierConstants:true,cacheSourceConstants:true}));
 return engines.get(key);
}
const squared=(E,z)=>{const a=E.c(z);return typeof a[0].n==='bigint'?a[0].n*a[0].n+a[1].n*a[1].n:E.norm2(a);};
const greater=(a,b)=>typeof a==='bigint'?a>b:a.gt(b),zero=a=>typeof a==='bigint'?a===0n:a.isZero();
function maximum(E,row){let index=0,best=squared(E,row[0]);for(let j=1;j<row.length;j++){const n=squared(E,row[j]);if(greater(n,best)){best=n;index=j;}}return E.abs(row[index]);}
function frozen(E,config,source,initial,deadline,metrics){
 const r=new E.D(String(source[1])).plus(config.shift[1]),s=new E.D(String(source[0])).plus(config.shift[0]),red=E.reduceModular(config.weights,initial,r,s);
 const get=w=>{ensure(deadline);metrics.matrixEvaluations++;return E.tvMatrix(red.weights,E.mobius(w,red.matrix),red.r,red.s,{reduce:false}).M;};
 const raw=get(initial),N=raw.length,rows=Array.from({length:N},(_,i)=>i),cols=[...rows],rs=raw.map(row=>maximum(E,row));
 if(rs.some(x=>x.isZero()))throw Error('Compact reference row rounded to zero');
 const M=raw.map((row,i)=>row.map(z=>E.div(z,rs[i]))),cs=Array.from({length:N},(_,j)=>maximum(E,M.map(row=>row[j])));
 if(cs.some(x=>x.isZero()))throw Error('Compact reference column rounded to zero');
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=E.div(M[i][j],cs[j]);
 const pivots=[];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=squared(E,M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const n=squared(E,M[i][k]);if(greater(n,best)){p=i;q=k;best=n;}}
  if(zero(best))throw Error('Compact reference minor rounded singular');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]];}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]];}
  pivots.push(E.abs(M[j][j]));
  for(let i=j+1;i<N;i++){const ratio=E.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=E.sub(M[i][k],E.mul(ratio,M[j][k]));}
 }
 let first=M[N-1][N-1];
 return{dimension:N,frame:red.matrix,referenceMinPivot:pivots.length?E.D.min(...pivots).toString():null,evaluate:w=>{
  if(first&&zero(squared(E,E.sub(w,initial)))){const value=first;first=null;return value;}
  const raw=get(w),A=rows.map(i=>cols.map(j=>E.div(E.div(raw[i][j],rs[i]),cs[j])));
  for(let j=0;j<N-1;j++){
   const pivot=E.abs(A[j][j]);if(pivot.isZero()||pivot.lt(pivots[j].times('1e-10')))throw Error('Compact frozen Schur minor lost its local pivot health');
   for(let i=j+1;i<N;i++){const ratio=E.div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=E.sub(A[i][k],E.mul(ratio,A[j][k]));}
  }
  return A[N-1][N-1];
 }};
}

export function compactMeshCorrect(config,source,guess,{precision=80,digits=18,deadline=Infinity,radius=.15,maxStep=.025,compact=true,derivativeProbe=null,decimal=false}={}){
 const begun=now(),useCompact=!decimal&&compact&&Math.max(...config.weights)<=9,E=engine(Math.max(precision,digits+30),useCompact,decimal),initial=E.c(guess),R=new E.D(radius),H=new E.D(maxStep),tolerance=new E.D(10).pow(-digits),metrics={matrixEvaluations:0,iterations:0,lineTrials:0,restarts:0};
 let current=initial,value,derivative,lastStep=null,model;
 const finish=(converged,message)=>({converged,tau:E.strs(current),precision:E.precision,requestDigits:digits,correction:lastStep?.toString()??null,pivotHealthy:converged,formulaInputsUsed:false,dimension:model?.dimension,frame:model?.frame,referenceMinPivot:model?.referenceMinPivot,arithmetic:decimal?'Decimal generic Fourier fallback':useCompact?'guarded compact fixed Fourier jets':'scaled stable fixed Fourier jets',metrics,milliseconds:now()-begun,message});
 try{
  ensure(deadline);if(!initial[1].gt('.00001')||!R.gt(0)||!H.gt(0))throw Error('Invalid compact local correction neighbourhood');
  const rebuild=()=>{model=frozen(E,config,source,current,deadline,metrics);value=model.evaluate(current);const h=E.D.min(new E.D('1e-6').times(E.abs(current).plus(1)),R.div(30)),candidate=derivativeProbe?E.c(derivativeProbe):null,delta=candidate?E.sub(candidate,current):null,probe=delta&&E.abs(delta).gt('1e-20')&&E.abs(delta).lte(R)?delta:[h,h.times(2)];derivative=E.div(E.sub(model.evaluate(E.add(current,probe)),value),probe);};
  rebuild();
  for(let i=0;i<28;i++){
   ensure(deadline);metrics.iterations++;
   if(zero(squared(E,value)))throw Error('Compact residual rounded to zero; increase working precision');
   if(zero(squared(E,derivative)))throw Error('Compact local derivative rounded to zero; increase working precision');
   let step=E.div(value,derivative),length=E.abs(step);lastStep=length;
   if(length.lt(tolerance))return finish(true,'Compact root corrected inside a healthy bounded Schur neighbourhood');
   if(length.gt(H))step=E.mul(step,H.div(length));
   let accepted=false,factor=new E.D(1),next,nextValue;
   for(let k=0;k<10;k++){
    ensure(deadline);metrics.lineTrials++;next=E.sub(current,E.mul(step,factor));
    if(next[1].gt('.00001')&&E.abs(E.sub(next,initial)).lte(R)){
     try{nextValue=model.evaluate(next);if(E.abs(nextValue).lt(E.abs(value))){accepted=true;break;}}catch{}
    }
    factor=factor.div(2);
   }
   if(!accepted){if(metrics.restarts++>=2)throw Error('No decreasing bounded compact Schur step');rebuild();continue;}
   derivative=E.div(E.sub(nextValue,value),E.sub(next,current));current=next;value=nextValue;
  }
  throw Error('Compact local iteration budget reached');
 }catch(error){return finish(false,String(error.message??error));}
}

export function checkedCompactMeshRoot(config,source,guess,{deadline=Infinity,radius=.15,agreementTolerance='1e-12'}={}){
 if(!(Number(agreementTolerance)>0&&Number(agreementTolerance)<=1e-12))throw Error('Compact agreement tolerance must be positive and at most1e-12');
 const begun=now();let last;
 for(const extra of [0,30,60,100]){
  const decimal=extra===100,compact=extra<60;
  ensure(deadline);const low=compactMeshCorrect(config,source,guess,{precision:50+extra,digits:18,deadline,radius,compact,decimal});last=low;if(!low.converged)continue;
  const E=engine(low.precision+30,false,decimal),probe=E.strs(E.add(E.c(low.tau),['1e-7','2e-7']));
  // Both high-precision values are freshly evaluated in the perturbed model.
  // The lower root is a useful derivative probe, never an accepted high value.
  const high=compactMeshCorrect(config,source,probe,{precision:low.precision+30,digits:28,deadline,radius:Math.min(.025,radius),derivativeProbe:low.tau,compact,decimal});last=high;if(!high.converged)continue;
  const views=t=>({raw:t,zero:E.sub(1,E.div(1,t)),one:E.div(1,E.sub(1,t)),half:E.div(t,E.sub(1,E.mul(2,t)))}),a=views(E.c(low.tau)),b=views(E.c(high.tau));
  const errors=Object.fromEntries(Object.keys(a).map(k=>[k,E.abs(E.sub(a[k],b[k])).toString()]));
  if(Object.values(errors).some(e=>new E.D(e).gt(agreementTolerance)))continue;
  return{...high,comparisonTauExact:low.tau,lowPrecision:low.precision,lowRequestedDigits:18,highRequestedDigits:28,precisionAgreementTolerance:String(agreementTolerance),crossPrecisionTauError:errors.raw,numericalErrorEstimate:E.D.max(new E.D(high.correction),new E.D(errors.raw)).toString(),viewErrors:errors,highPrecisionProvenance:'independently perturbed fresh compact correction',endpointChecked:true,milliseconds:now()-begun,metrics:{matrixEvaluations:low.metrics.matrixEvaluations+high.metrics.matrixEvaluations,iterations:low.metrics.iterations+high.metrics.iterations},formulaInputsUsed:false};
 }
 const error=Error('Compact endpoint precision repeats did not agree: '+(last?.message??'unresolved'));error.corrector=last;throw error;
}

// Plot-level checks strengthen the original binary64 compact preview. They
// never replace graph anchoring, standard point checks, or cusp refinement.
const finiteComplex=z=>Array.isArray(z)&&z.length===2&&z.every(Number.isFinite);
const fastMobius=(z,[a,b,c,d])=>F.div(F.add(F.mul(a,z),b),F.add(F.mul(c,z),d));
function plotModel(config,source,initial,deadline,metrics){
 metrics.models=(metrics.models??0)+1;
 const red=F.reduceModular(config.weights,initial,Number(source[1])+config.shift[1],Number(source[0])+config.shift[0]);
 if(!finiteComplex(red.tau)||red.tau[1]<.8||red.tau[1]>8)throw Error('Compact plot modular height requires guarded precision');
 const get=w=>{ensure(deadline);metrics.matrixEvaluations++;const M=F.tvMatrix(red.weights,fastMobius(w,red.matrix),red.r,red.s,{reduce:false}).M;if(!M.length||M.some(row=>row.length!==M.length||row.some(z=>!finiteComplex(z))))throw Error('Nonfinite compact plot matrix');return M;};
 const raw=get(initial),N=raw.length,rows=Array.from({length:N},(_,i)=>i),cols=[...rows],rs=raw.map(row=>Math.max(...row.map(F.abs)));
 if(rs.some(x=>!(x>0&&Number.isFinite(x))))throw Error('Unresolved compact plot row scale');
 const M=raw.map((row,i)=>row.map(z=>F.div(z,rs[i]))),cs=Array.from({length:N},(_,j)=>Math.max(...M.map(row=>F.abs(row[j]))));
 if(cs.some(x=>!(x>0&&Number.isFinite(x))))throw Error('Unresolved compact plot column scale');
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=F.div(M[i][j],cs[j]);
 const pivots=[];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=F.abs(M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const value=F.abs(M[i][k]);if(value>best){p=i;q=k;best=value;}}
  if(!(best>1e-12&&Number.isFinite(best)))throw Error('Unresolved compact plot reference minor');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]];}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]];}
  pivots.push(F.abs(M[j][j]));
  for(let i=j+1;i<N;i++){const ratio=F.div(M[i][j],M[j][j]);if(!finiteComplex(ratio))throw Error('Unsafe compact plot pivot division');for(let k=j+1;k<N;k++)M[i][k]=F.sub(M[i][k],F.mul(ratio,M[j][k]));}
 }
 let first=M[N-1][N-1];
 return w=>{
  let value;
  if(first&&w[0]===initial[0]&&w[1]===initial[1]){value=first;first=null;}
  else{
   const raw=get(w),A=rows.map(i=>cols.map(j=>F.div(F.div(raw[i][j],rs[i]),cs[j])));
   for(let j=0;j<N-1;j++){
    const pivot=F.abs(A[j][j]);if(!(Number.isFinite(pivot)&&pivot>1e-12&&pivot>=pivots[j]*1e-10))throw Error('Compact plot frozen minor lost pivot health');
    for(let i=j+1;i<N;i++){const ratio=F.div(A[i][j],A[j][j]);if(!finiteComplex(ratio))throw Error('Unsafe compact plot pivot division');for(let k=j+1;k<N;k++)A[i][k]=F.sub(A[i][k],F.mul(ratio,A[j][k]));}
   }
   value=A[N-1][N-1];
  }
  if(!finiteComplex(value)||F.abs(value)===0)throw Error('Compact plot residual unresolved at binary64 precision');
  return value;
 };
}
function plotCorrect(config,source,guess,{deadline,radius,metrics}){
 const initial=guess.map(Number);if(!finiteComplex(initial)||initial[1]<=.001||!(radius>0&&Number.isFinite(radius)))throw Error('Invalid compact plot neighbourhood');
 const model=plotModel(config,source,initial,deadline,metrics),h=Math.min(1e-6*(1+F.abs(initial)),radius/30),probe=[h,2*h];
 let current=initial,value=model(current),derivative=F.div(F.sub(model(F.add(current,probe)),value),probe);
 for(let i=0;i<28;i++){
  ensure(deadline);if(!finiteComplex(derivative)||F.abs(derivative)===0)throw Error('Unresolved compact plot derivative');
  let step=F.div(value,derivative),length=F.abs(step);if(!finiteComplex(step)||!Number.isFinite(length))throw Error('Nonfinite compact plot correction');
  if(length<=3e-12)return{tau:current,correction:length,iterations:i+1};
  if(length>.025)step=F.mul(.025/length,step);
  let accepted=false,factor=1,next,nextValue;
  for(let k=0;k<10;k++){
   ensure(deadline);next=F.sub(current,F.mul(factor,step));
   if(finiteComplex(next)&&next[1]>.001&&F.abs(F.sub(next,initial))<=radius){nextValue=model(next);if(F.abs(nextValue)<F.abs(value)){accepted=true;break;}}
   factor/=2;
  }
  if(!accepted)throw Error('No decreasing bounded compact plot step');
  derivative=F.div(F.sub(nextValue,value),F.sub(next,current));current=next;value=nextValue;
 }
 throw Error('Compact plot correction iteration budget reached');
}
const plotViews=t=>({raw:t,zero:F.sub(1,F.div(1,t)),one:F.div(1,F.sub(1,t)),half:F.div(t,F.sub(1,F.mul(2,t)))});
function propagatedPlotCorrection(t,error,view){
 if(view==='raw')return error;
 const denominator=view==='zero'?F.abs(t):view==='one'?F.abs(F.sub(1,t)):F.abs(F.sub(1,F.mul(2,t))),margin=denominator-(view==='half'?2:1)*error;
 return margin>0?error/(denominator*margin):Infinity;
}
export function compactMeshPlotRoot(config,source,guess,{deadline=Infinity,radius=.15}={}){
 const begun=now(),metrics={matrixEvaluations:0};
 try{
  if(Math.max(...config.weights)>9)throw Error('Compact plot derivative degree requires guarded precision');
  const initial=guess.map(Number),neighbourhood=Math.min(radius,.13*(1+F.abs(initial))),perturb=[1e-6,2e-6];
  const a=plotCorrect(config,source,initial,{deadline,radius:neighbourhood,metrics}),b=plotCorrect(config,source,F.add(initial,perturb),{deadline,radius:neighbourhood-F.abs(perturb),metrics});
  const av=plotViews(a.tau),bv=plotViews(b.tau),errors={};
  for(const name of Object.keys(av)){
   if(!finiteComplex(av[name])||!finiteComplex(bv[name]))throw Error('Unsafe compact plot view');
   const error=F.abs(F.sub(av[name],bv[name]))+propagatedPlotCorrection(a.tau,a.correction,name)+propagatedPlotCorrection(b.tau,b.correction,name)+64*Number.EPSILON*(1+F.abs(av[name])+F.abs(bv[name]));
   if(!(Number.isFinite(error)&&error<=1e-9))throw Error('Independent compact plot views did not agree');errors[name]=String(error);
  }
  return{converged:true,tau:b.tau.map(String),precision:15,correction:String(b.correction),numericalErrorEstimate:errors.raw,viewErrors:errors,plotComparisonTauExact:a.tau.map(String),crossModelTauError:String(F.abs(F.sub(a.tau,b.tau))),plotChecked:true,twoFreshModels:metrics.models===2,pivotHealthy:true,plotCorrectionTarget:'3e-12',plotErrorTolerance:'1e-9',plotAttemptMatrixEvaluations:metrics.matrixEvaluations,metrics,milliseconds:now()-begun,formulaInputsUsed:false,arithmetic:'two independently frozen binary64 Schur models',message:'Bounded compact plot roots agree in every displayed chart; numerical estimates, not interval bounds'};
 }catch(error){
  const high=checkedCompactMeshRoot(config,source,guess,{deadline,radius});
  return{...high,plotFallback:true,plotFallbackReason:String(error.message??error),plotAttemptMatrixEvaluations:metrics.matrixEvaluations,milliseconds:now()-begun};
 }
}

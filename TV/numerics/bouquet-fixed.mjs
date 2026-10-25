// Stable browser evaluator. theta/exp Taylor variable is scaled by jetScale.
// Matrix relation: M_scaled[(i,l),(j,h)]=jetScale^(l+h+1)*M_original.
// Browser-safe fixed-point Bloch jets for the same elliptic matching matrix.
// Precision includes a derivative guard; independent precision comparisons
// remain required. These are numerical estimates, not interval certificates.
import {createFixed} from './bigint-fixed.mjs?v=20261024-tetra-dense4';

export function createEvaluator(precision=80,{jetScale=1,thetaMode='stable',axisShortcuts=false,cacheMatrixConstants=false,compactFourier=false,cacheFourierConstants=false,cacheSourceConstants=false}={}) {
 if(!(Number(jetScale)>0&&Number.isFinite(Number(jetScale))))throw Error('A positive finite derivative scale is required');
 if(!Number.isInteger(precision)||precision<25||precision>970)throw Error('Decimal precision must be an integer from 25 to 970');
 const D=createFixed(precision);
 const Z=new D(0),O=new D(1),PI=D.acos(-1),TPI=PI.times(2),JS=new D(jetScale);
 const d=x=>x?.constructor===D?x:new D(x);
 const c=x=>Array.isArray(x)?[d(x[0]),d(x[1])]:[d(x),Z];
 const add=(a,b)=>{a=c(a);b=c(b);return[a[0].plus(b[0]),a[1].plus(b[1])]};
 const neg=a=>{a=c(a);return[a[0].neg(),a[1].neg()]};
 const sub=(a,b)=>add(a,neg(b));
 const minusOne=O.neg(),scale=(x,k)=>k.eq(O)?x:k.eq(minusOne)?x.neg():x.times(k);
 // The optional branches remove operations involving exact zeros and units.
 // They preserve the fixed-point truncation of the general formulas.
 const mul=(a,b)=>{a=c(a);b=c(b);
  if(axisShortcuts){
   if(b[1].isZero())return[scale(a[0],b[0]),scale(a[1],b[0])];
   if(b[0].isZero())return[scale(a[1],b[1]).neg(),scale(a[0],b[1])];
   if(a[1].isZero())return[scale(b[0],a[0]),scale(b[1],a[0])];
   if(a[0].isZero())return[scale(b[1],a[1]).neg(),scale(b[0],a[1])];
  }
  return[a[0].times(b[0]).minus(a[1].times(b[1])),a[0].times(b[1]).plus(a[1].times(b[0]))]};
 const div=(a,b)=>{a=c(a);b=c(b);
  if(axisShortcuts){
   if(b[1].isZero()){
    if(b[0].eq(O))return[a[0],a[1]];
    if(b[0].eq(minusOne))return[a[0].neg(),a[1].neg()];
    return[a[0].div(b[0]),a[1].div(b[0])];
   }
   if(b[0].isZero()){
    if(b[1].eq(O))return[a[1],a[0].neg()];
    if(b[1].eq(minusOne))return[a[1].neg(),a[0]];
    return[a[1].div(b[1]),a[0].neg().div(b[1])];
   }
  }
  if(b[0].abs().gte(b[1].abs())){const r=b[1].div(b[0]),v=b[0].plus(b[1].times(r));return[a[0].plus(a[1].times(r)).div(v),a[1].minus(a[0].times(r)).div(v)]}const r=b[0].div(b[1]),v=b[1].plus(b[0].times(r));return[a[0].times(r).plus(a[1]).div(v),a[1].times(r).minus(a[0]).div(v)]};
 const abs=a=>{a=c(a);if(axisShortcuts){if(a[0].isZero())return a[1].abs();if(a[1].isZero())return a[0].abs();}const hi=D.max(a[0].abs(),a[1].abs()),lo=D.min(a[0].abs(),a[1].abs());return hi.isZero()?Z:hi.times(O.plus(lo.div(hi).pow(2)).sqrt())};
 // Only comparisons and zero tests use norm2; the norm has the same order
 // and avoids squaring a resolved tiny residual below the fixed-point scale.
 const norm2=abs;
 const exp=a=>{a=c(a);const mag=a[0].exp();return[mag.times(a[1].cos()),mag.times(a[1].sin())]};
 const strs=a=>a.map(x=>x.toString());

 function reduceModular(weights,tau,r,s) {
  weights=[...weights];tau=c(tau);r=new D(r);s=new D(s);let matrix=[1,0,0,1];
  for(let steps=0;steps<100;steps++) {
   const m=tau[0].round().toNumber();
   if(m) {tau[0]=tau[0].minus(m);r=r.plus(s.times(m));
    if(Math.abs(m)%2)weights=[weights[0],weights[1],weights[3],weights[2]];
    matrix=[matrix[0]-m*matrix[2],matrix[1]-m*matrix[3],matrix[2],matrix[3]];
   }
   if(norm2(tau).gte(O.minus(new D(10).pow(-precision+8))))break;
   tau=div(-1,tau);[r,s]=[s,r.neg()];weights=[weights[0],weights[2],weights[1],weights[3]];
   matrix=[-matrix[2],-matrix[3],matrix[0],matrix[1]];
  }
  r=r.minus(r.round());s=s.minus(s.round());
  return{weights,tau,r,s,matrix};
 }

 // Taylor coefficients in v=2*pi*z; omit the common q^(1/4).
 function thetaJetsLegacy(z,tau,N) {
  z=c(z);const y=tau[1].toNumber(),zi=Math.abs(z[1].toNumber());
  const count=Math.max(4,Math.ceil(Math.sqrt((precision*Math.LN10+N*Math.log(N+2)+N*Math.log(Math.max(1,Number(jetScale)))+70)/(Math.PI*y))+zi/y)+3);
  const out=Array.from({length:N+1},()=>c(0));
  // exp terms use recurrence rather than repeated trigonometric functions.
  const modeScale=JS.pow(N);
  let qp=c(modeScale),qratio=exp(mul([Z,PI.times(2)],tau));
  let ep=exp(mul([Z,PI],z)),em=div(1,ep),epRatio=exp(mul([Z,TPI],z)),emRatio=div(1,epRatio);
  const qbase=qratio;
  for(let n=0;n<count;n++) {
   const k=new D(n).plus('.5'),factor=mul([Z,n%2?1:-1],qp);let fp=c(1),fm=c(1);
   for(let d=0;d<=N;d++) {
    out[d]=add(out[d],div(mul(factor,sub(mul(fp,ep),mul(fm,em))),modeScale));
    fp=mul(fp,[Z,k.times(JS).div(d+1)]);fm=mul(fm,[Z,k.neg().times(JS).div(d+1)]);
   }
   qp=mul(qp,qratio);qratio=mul(qratio,qbase);ep=mul(ep,epRatio);em=mul(em,emRatio);
  }
  return out;
 }
 // The derivative coefficients are independent of tau and z. Keep the
 // exact original recurrence, bounded by the compact degree and mode budget.
 const fourierConstants=new Map(),fourierStats={recurrenceCalls:0,directCalls:0,guardFallbackCalls:0,coefficientCacheHits:0,coefficientCacheMisses:0};
 const getFourierDiagnostics=()=>({...fourierStats,coefficientCacheModes:fourierConstants.size});
 function coefficients(n,N){
  if(!cacheFourierConstants||N>18||n>=128)return null;
  let item=fourierConstants.get(n);
  if(item)fourierStats.coefficientCacheHits++;else fourierStats.coefficientCacheMisses++;
  if(!item){item={k:new D(n).plus('.5'),plus:[c(1)],minus:[c(1)]};fourierConstants.set(n,item);}
  for(let order=item.plus.length-1;order<N;order++){
   item.plus.push(mul(item.plus[order],[Z,item.k.times(JS).div(order+1)]));
   item.minus.push(mul(item.minus[order],[Z,item.k.neg().times(JS).div(order+1)]));
  }
  return item;
 }
 // Combine the Fourier damping and the two traveling-wave amplitudes
 // BEFORE fixed-point rounding. Inverting exp(2*pi*i*z) is unsafe when
 // Im(z) is large; splitting a tiny q factor from a huge wave is unsafe too.
 // The extra common derivative scale stays inside exp until after the
 // Taylor coefficient multiplier, preserving terms amplified by derivatives.
 function thetaJetsStable(z,tau,N) {
  fourierStats.directCalls++;z=c(z);const y=tau[1].toNumber(),zi=Math.abs(z[1].toNumber());
  const count=Math.max(4,Math.ceil(Math.sqrt((precision*Math.LN10+N*Math.log(N+2)+N*Math.log(Math.max(1,Number(jetScale)))+70)/(Math.PI*y))+zi/y)+3);
  const out=Array.from({length:N+1},()=>c(0)),modeScale=JS.pow(N),logScale=JS.ln().times(N);
  for(let n=0;n<count;n++) {
   const k=new D(n).plus('.5'),harmonic=2*n+1,frequency=n*(n+1),sign=n%2?1:-1;
   const exponentP=add(mul([Z,PI],add(mul(frequency,tau),mul(harmonic,z))),logScale);
   const exponentM=add(mul([Z,PI],sub(mul(frequency,tau),mul(harmonic,z))),logScale);
   const ep=exp(exponentP),em=exp(exponentM),factor=c([0,sign]);
   const cached=coefficients(n,N);let fp=c(1),fm=c(1);
   for(let order=0;order<=N;order++) {
    if(cached){fp=cached.plus[order];fm=cached.minus[order];}
    out[order]=add(out[order],div(mul(factor,sub(mul(fp,ep),mul(fm,em))),modeScale));
    if(!cached){fp=mul(fp,[Z,k.times(JS).div(order+1)]);fm=mul(fm,[Z,k.neg().times(JS).div(order+1)]);}
   }
  }
  return out;
 }
 // For bounded compact images, multiply the combined Fourier
 // amplitudes, never a tiny damping factor and a large separate wave.
 function thetaJetsCompact(z,tau,N) {
  z=c(z);tau=c(tau);
  if(N>18||!JS.eq(1)||tau[1].lt('.8')||tau[1].gt(8)||z[1].abs().gt(tau[1])){fourierStats.guardFallbackCalls++;return thetaJetsStable(z,tau,N);}
  fourierStats.recurrenceCalls++;
  const y=tau[1].toNumber(),zi=Math.abs(z[1].toNumber());
  const count=Math.max(4,Math.ceil(Math.sqrt((precision*Math.LN10+N*Math.log(N+2)+N*Math.log(Math.max(1,Number(jetScale)))+70)/(Math.PI*y))+zi/y)+3);
  const out=Array.from({length:N+1},()=>c(0)),modeScale=JS.pow(N),logScale=JS.ln().times(N);
  let ep=exp(add(mul([Z,PI],z),logScale)),em=exp(add(mul([Z,PI],neg(z)),logScale));
  let rp=exp(mul([Z,TPI],add(tau,z))),rm=exp(mul([Z,TPI],sub(tau,z)));
  const qr=exp(mul([Z,TPI],tau));
  for(let n=0;n<count;n++) {
   const k=new D(n).plus('.5'),factor=c([0,n%2?1:-1]);
   const cached=coefficients(n,N);let fp=c(1),fm=c(1);
   for(let order=0;order<=N;order++) {
    if(cached){fp=cached.plus[order];fm=cached.minus[order];}
    out[order]=add(out[order],div(mul(factor,sub(mul(fp,ep),mul(fm,em))),modeScale));
    if(!cached){fp=mul(fp,[Z,k.times(JS).div(order+1)]);fm=mul(fm,[Z,k.neg().times(JS).div(order+1)]);}
   }
   ep=mul(ep,rp);em=mul(em,rm);rp=mul(rp,qr);rm=mul(rm,qr);
  }
  return out;
 }
 const thetaJets=(z,tau,N)=>(thetaMode==='legacy'?thetaJetsLegacy:compactFourier?thetaJetsCompact:thetaJetsStable)(z,tau,N);
 function polyDiv(a,b,N) {
  const out=[];
  for(let k=0;k<=N;k++) {
   let v=a[k]??c(0);
   for(let j=1;j<=Math.min(k,b.length-1);j++)v=sub(v,mul(b[j],out[k-j]));
   out.push(div(v,b[0]));
  }
  return out;
 }
 function polyMul(a,b,N) {
  const out=[];
  for(let k=0;k<=N;k++) {
   let v=c(0);for(let j=Math.max(0,k-b.length+1);j<=Math.min(k,a.length-1);j++)v=add(v,mul(a[j],b[k-j]));
   out.push(v);
  }
  return out;
 }
 const choose=(n,k)=>{let v=O;for(let i=1;i<=k;i++)v=v.times(n-k+i).div(i);return v};
 // Optional compact-mesh caches contain only already computed exact entries.
 // Source phases do not depend on tau; binomial coefficients depend on neither
 // tau nor the source. Keep both caches bounded for long-lived workers.
 const sourceConstants=new Map(),matrixCoefficients=new Map(),matrixCacheStats={sourceHits:0,sourceMisses:0,phaseHits:0,phaseMisses:0,coefficientHits:0,coefficientMisses:0};
 const getMatrixCacheDiagnostics=()=>({...matrixCacheStats,sourcePairs:sourceConstants.size,coefficients:matrixCoefficients.size});
 function tvMatrix(weights,tau,r,s,{reduce=true}={}) {
  tau=c(tau);r=new D(r);s=new D(s);
  if(!tau[1].gt(0))throw Error('Upper half-plane required');
  if(reduce)({weights,tau,r,s}=reduceModular(weights,tau,r,s));
  else{r=r.minus(r.round());s=s.minus(s.round())}
  const sites=[];for(let i=0;i<4;i++)for(let h=weights[i]-1;h>=0;h-=2)sites.push([i,h]);
  const degree=2*Math.max(0,...sites.map(x=>x[1])),N=degree+1;
  const alpha=add(r,mul(s,tau)),zero=thetaJets(0,tau,N+1),num=thetaJets(alpha,tau,N),den=zero.slice(1);
  const expjet=[c(1)];for(let d=1;d<=N;d++)expjet.push(mul(expjet[d-1],[Z,s.times(JS).div(d)]));
  const factor=div(zero[1],num[0]),diag=polyMul(expjet,polyDiv(num,den,N),N).map(x=>mul(x,factor));
  const bits=[[0,0],[1,0],[0,1],[1,1]],points=[c(0),c('.5'),mul('.5',tau),mul('.5',add(1,tau))],jets=new Map();
  // A matrix has only four possible period translations. Their phases and
  // its repeated binomial coefficients are independent of the matrix entry.
  // Reuse exact computed values; no Fourier term or digit is omitted.
  let sourceEntry=null;
  if(cacheSourceConstants){
   const sourceKey=r.toString()+':'+s.toString();sourceEntry=sourceConstants.get(sourceKey);
   if(sourceEntry)matrixCacheStats.sourceHits++;
   else{matrixCacheStats.sourceMisses++;sourceEntry={phases:new Map()};sourceConstants.set(sourceKey,sourceEntry);while(sourceConstants.size>64)sourceConstants.delete(sourceConstants.keys().next().value);}
  }
  const phases=sourceEntry?.phases??(cacheMatrixConstants?new Map():null),coefficients=cacheSourceConstants?matrixCoefficients:cacheMatrixConstants?new Map():null;
  const coefficient=(n,k)=>{
   if(!coefficients)return choose(n,k);
   const key=n+':'+k;
   if(!coefficients.has(key)){if(cacheSourceConstants)matrixCacheStats.coefficientMisses++;coefficients.set(key,choose(n,k));if(cacheSourceConstants&&coefficients.size>256)coefficients.delete(coefficients.keys().next().value);}
   else if(cacheSourceConstants)matrixCacheStats.coefficientHits++;
   return coefficients.get(key);
  };
  const M=sites.map(([i,l])=>sites.map(([j,h])=>{
   const d=l+h;if(i===j)return mul(coefficient(d,l),diag[d+1]);
   const du=bits[i][0]-bits[j][0],dv=bits[i][1]-bits[j][1],hu=(du+2)%2,hv=(dv+2)%2;
   const key=hu+2*hv,periodU=(du-hu)/2,periodV=(dv-hv)/2;
   if(!jets.has(key)) {
    const z=points[key],scalar=mul(exp(mul([Z,TPI.times(s)],z)),factor);
    jets.set(key,polyMul(expjet,polyDiv(thetaJets(add(z,alpha),tau,degree),thetaJets(z,tau,degree),degree),degree).map(x=>mul(x,scalar)));
   }
   const phaseKey=periodU+2*periodV;
   let periodFactor=phases?.get(phaseKey);
   if(!periodFactor){if(cacheSourceConstants)matrixCacheStats.phaseMisses++;const phase=TPI.times(s.times(periodU).minus(r.times(periodV)));periodFactor=exp([Z,phase]);phases?.set(phaseKey,periodFactor);}
   else if(cacheSourceConstants)matrixCacheStats.phaseHits++;
   return mul(coefficient(d,l),mul(periodFactor,jets.get(key)[d]));
  }));
  return{M,sites,weights,tau,r,s};
 }
 function determinant(weights,tau,r,s,options={}) {
  if(options.isogeny&&weights[0]===weights[1]&&weights[2]===weights[3]) {
   const [u,v]=options.isogenyLift??[0,0];r=new D(r);s=new D(s);
   if(weights[0]===weights[2])return determinant([weights[0],0,0,0],tau,r.plus(u).div(2),s.plus(v).div(2),{...options,isogeny:false});
   if(weights[2]===0)return determinant([weights[0],0,0,0],mul(2,tau),r,s.plus(v).div(2),{...options,isogeny:false});
   return determinant([weights[0],weights[2],0,0],div(-1,mul(2,tau)),s.plus(v).div(2),r.neg(),{...options,isogeny:false});
  }
  const{M}=tvMatrix(weights,tau,r,s,options),N=M.length;let product=c(1),minRelativePivot=O;
  for(let i=0;i<N;i++) {
   let scale=N===1?O:D.max(...M[i].map(abs));
   if(scale.isZero())return{value:c(0),dimension:N,minRelativePivot:Z};
   for(let j=0;j<N;j++)M[i][j]=div(M[i][j],scale);
  }
  for(let j=0;j<N;j++) {
   let p=j;for(let i=j+1;i<N;i++)if(norm2(M[i][j]).gt(norm2(M[p][j])))p=i;
   const pivot=abs(M[p][j]);minRelativePivot=D.min(minRelativePivot,pivot);
   if(pivot.isZero())return{value:c(0),dimension:N,minRelativePivot:Z};
   if(p!==j){[M[p],M[j]]=[M[j],M[p]];product=neg(product)}
   product=mul(product,M[j][j]);
   for(let i=j+1;i<N;i++) {
    const ratio=div(M[i][j],M[j][j]);
    for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]));
   }
  }
  return{value:product,dimension:N,minRelativePivot};
 }
 function solveFunction(f,guess,{digits=precision-25,maxsteps=35,perturb='1e-6',guard=()=>true}={}) {
  let a=c(guess),b=add(a,[perturb,new D(perturb).times(2)]),fa=f(a),fb=f(b),evaluations=2;
  const tolerance=new D(10).pow(-digits);let lastStep=null;
  for(let i=0;i<maxsteps;i++) {
   if(norm2(fb).isZero())throw Error('Exact zero determinant at finite precision: increase precision');
   const denominator=sub(fb,fa);if(norm2(denominator).isZero())throw Error('Zero secant denominator');
   const step=div(mul(fb,sub(b,a)),denominator);lastStep=step;
   if(abs(step).lt(tolerance))return{root:b,iterations:i+2,evaluations,residual:abs(fb),lastStep};
   a=b;fa=fb;b=sub(b,step);
   if(!b.every(x=>x.isFinite())||!guard(b))throw Error('Unsafe secant iterate');
   fb=f(b);evaluations++;
  }
  throw Error('No secant convergence: '+strs(lastStep));
 }
 const mobius=(z,[a,b,c0,d])=>div(add(mul(a,z),b),add(mul(c0,z),d));
 const inverse=(z,[a,b,c0,d])=>div(sub(mul(d,z),b),sub(a,mul(c0,z)));
 return{D,precision,getFourierDiagnostics,getMatrixCacheDiagnostics,jetScale:Number(jetScale),thetaMode,debugThetaJets:thetaJets,PI,TPI,c,add,sub,neg,mul,div,exp,abs,norm2,strs,reduceModular,tvMatrix,determinant,solveFunction,mobius,inverse};
}


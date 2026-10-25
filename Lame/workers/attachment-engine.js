/* Interior branch continuation and mesh sampling. The two factories retain
 * their original private state; callers select the existing precision mode. */
function createClassicalPrecision(options={}) {
const createPrecisionEngine=createLamePrecisionEngine(options);
const schurPredictor=(()=>{// Standalone Bloch-jet determinant. No website imports or mutations.
const C=x=>typeof x==='number'?[x,0]:x;
const add=(...xs)=>xs.reduce((a,b)=>{b=C(b);return[a[0]+b[0],a[1]+b[1]]},[0,0]);
const neg=a=>[-a[0],-a[1]];
const sub=(a,b)=>add(a,neg(C(b)));
const mul=(a,b)=>{a=C(a);b=C(b);return[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]]};
const div=(a,b)=>{a=C(a);b=C(b);const d=b[0]**2+b[1]**2;return[(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d]};
const abs=a=>Math.hypot(...a);
const exp=a=>{const r=Math.exp(a[0]);return[r*Math.cos(a[1]),r*Math.sin(a[1])]};

function reduceModular(weights,tau,r,s){
 weights=[...weights];tau=[...tau];let matrix=[1,0,0,1];
 for(let steps=0;steps<60;steps++){
  const m=Math.round(tau[0]);
  if(m){
   tau[0]-=m;r+=m*s;
   if(Math.abs(m)%2)weights=[weights[0],weights[1],weights[3],weights[2]];
   matrix=[matrix[0]-m*matrix[2],matrix[1]-m*matrix[3],matrix[2],matrix[3]];
  }
  if(abs(tau)>=1-1e-13)break;
  tau=div(-1,tau);[r,s]=[s,-r];weights=[weights[0],weights[2],weights[1],weights[3]];
  matrix=[-matrix[2],-matrix[3],matrix[0],matrix[1]];
 }
 r-=Math.round(r);s-=Math.round(s);
 return {weights,tau,r,s,matrix};
}

// Taylor coefficients in v=2*pi*z; theta common q^(1/4) omitted.
function thetaJets(z,tau,N){
 z=C(z);const count=Math.max(4,Math.ceil(Math.sqrt(42/(Math.PI*tau[1]))+Math.abs(z[1])/tau[1])+3);
 const out=Array.from({length:N+1},()=>[0,0]);
 for(let n=0;n<count;n++){
  const k=n+.5,q=exp(mul([0,Math.PI*n*(n+1)],tau));
  const ep=exp(mul([0,2*Math.PI*k],z)),em=div(1,ep),factor=mul([0,-((-1)**n)],q);
  let fp=[1,0],fm=[1,0];
  for(let d=0;d<=N;d++){
   out[d]=add(out[d],mul(factor,sub(mul(fp,ep),mul(fm,em))));
   fp=mul(fp,[0,k/(d+1)]);fm=mul(fm,[0,-k/(d+1)]);
  }
 }
 return out;
}
function polyDiv(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=a[k]??[0,0];
  for(let j=1;j<=Math.min(k,b.length-1);j++)v=sub(v,mul(b[j],out[k-j]));
  out.push(div(v,b[0]));
 }
 return out;
}
function polyMul(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=[0,0];for(let j=Math.max(0,k-b.length+1);j<=Math.min(k,a.length-1);j++)v=add(v,mul(a[j],b[k-j]));
  out.push(v);
 }
 return out;
}
function choose(n,k){let v=1;for(let i=1;i<=k;i++)v=v*(n-k+i)/i;return v}

function tvMatrix(weights,tau,r,s,{reduce=true}={}){
 if(!(tau[1]>0))throw Error('Upper half-plane required');
 if(reduce)({weights,tau,r,s}=reduceModular(weights,tau,r,s));
 else {r-=Math.round(r);s-=Math.round(s)}
 const sites=[];
 for(let i=0;i<4;i++)for(let h=weights[i]-1;h>=0;h-=2)sites.push([i,h]);
 const degree=2*Math.max(0,...sites.map(x=>x[1])),N=degree+1;
 const alpha=add(r,mul(s,tau)),zero=thetaJets(0,tau,N+1),num=thetaJets(alpha,tau,N),den=zero.slice(1);
 const expjet=[[1,0]];
 for(let d=1;d<=N;d++)expjet.push(mul(expjet[d-1],[0,s/d]));
 const factor=div(zero[1],num[0]),diag=polyMul(expjet,polyDiv(num,den,N),N).map(x=>mul(x,factor));
 const bits=[[0,0],[1,0],[0,1],[1,1]],points=[[0,0],[.5,0],mul(.5,tau),mul(.5,add(1,tau))],jets=new Map();
 const M=sites.map(([i,l])=>sites.map(([j,h])=>{
  const d=l+h;
  if(i===j)return mul(choose(d,l),diag[d+1]);
  const du=bits[i][0]-bits[j][0],dv=bits[i][1]-bits[j][1],hu=(du+2)%2,hv=(dv+2)%2;
  const key=hu+2*hv,pointIndex=key,periodU=(du-hu)/2,periodV=(dv-hv)/2;
  if(!jets.has(key)){
   const z=points[pointIndex],scalar=mul(exp(mul([0,2*Math.PI*s],z)),factor);
   const jet=polyMul(expjet,polyDiv(thetaJets(add(z,alpha),tau,degree),thetaJets(z,tau,degree),degree),degree).map(x=>mul(x,scalar));
   jets.set(key,jet);
  }
  const periodFactor=exp([0,2*Math.PI*(s*periodU-r*periodV)]);
  return mul(choose(d,l),mul(periodFactor,jets.get(key)[d]));
 }));
 return {M,sites,weights,tau,r,s};
}
function determinant(weights,tau,r,s,options={}){
 // This shortcut selects a square-root character lift. It is NOT the full
 // zero locus of the larger-period equation. Enable only after the graph
 // branch and its lift have been calibrated; full generic is the default.
 if(options.isogeny && weights[0]===weights[1] && weights[2]===weights[3]){
  const [u,v]=options.isogenyLift??[0,0];
  if(weights[0]===weights[2])return determinant([weights[0],0,0,0],tau,(r+u)/2,(s+v)/2,{...options,isogeny:false});
  if(weights[2]===0)return determinant([weights[0],0,0,0],mul(2,tau),r,(s+v)/2,{...options,isogeny:false});
  return determinant([weights[0],weights[2],0,0],div(-1,mul(2,tau)),(s+v)/2,-r,{...options,isogeny:false});
 }
 const {M}=tvMatrix(weights,tau,r,s,options),N=M.length;
 let product=[1,0],minRelativePivot=1;
 // Row normalization keeps intermediate products in the double range.
 for(let i=0;i<N;i++){
  // A 1x1 row is the residual itself: dividing by its absolute value would
  // erase its zero and make secant iteration impossible.
  const scale=N===1?1:Math.max(...M[i].map(abs));
  if(!scale||!Number.isFinite(scale))return {value:[NaN,NaN],minRelativePivot:0,dimension:N};
  for(let j=0;j<N;j++)M[i][j]=mul(1/scale,M[i][j]);
 }
 for(let j=0;j<N;j++){
  let p=j;for(let i=j+1;i<N;i++)if(abs(M[i][j])>abs(M[p][j]))p=i;
  const pivot=abs(M[p][j]);minRelativePivot=Math.min(minRelativePivot,pivot);
  if(pivot<1e-300)return {value:[0,0],minRelativePivot:0,dimension:N};
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];product=neg(product)}
  product=mul(product,M[j][j]);
  for(let i=j+1;i<N;i++){
   const ratio=div(M[i][j],M[j][j]);
   for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]));
  }
 }
 return {value:product,minRelativePivot,dimension:N};
}

function solve(weights,r,s,guess,{tol=1e-11,maxsteps=32,isogeny=false,isogenyLift=[0,0]}={}){
 let a=guess,b=add(guess,[1e-6,2e-6]),fa=determinant(weights,a,r,s,{isogeny,isogenyLift}).value,fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 for(let i=0;i<maxsteps;i++){
  const denominator=sub(fb,fa);if(abs(denominator)===0)break;
  const step=div(mul(fb,sub(b,a)),denominator);
  if(abs(step)<tol)return {tau:b,iterations:i+2,residual:abs(fb)};
  a=b;fa=fb;b=sub(b,step);
  if(!Number.isFinite(abs(b))||b[1]<=.005)throw Error('Unsafe iterate');
  fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 }
 throw Error('No convergence');
}

// Fixed row/column equilibration and complete-pivot Schur residual.
// Unlike renormalizing a determinant at every iterate, this preserves a
// holomorphic residual and resolves the formerly unstable count-32 anchor.
function fixedSchur(E,get,initial){
 const abs=E.abs,div=E.div,mul=E.mul,sub=E.sub;
 const ref=get(initial),N=ref.length;if(!N)throw Error('Empty elliptic matching matrix');
 const greater=(a,b)=>a?.gt?a.gt(b):a>b;
 const maximum=a=>a.reduce((m,x)=>greater(abs(x),m)?abs(x):m,abs(a[0]));
 const scales=ref.map(maximum),M=ref.map((row,i)=>row.map(x=>div(x,scales[i])));
 const cs=Array.from({length:N},(_,j)=>maximum(M.map(row=>row[j])));
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=div(M[i][j],cs[j]);
 const rows=Array.from({length:N},(_,i)=>i),cols=[...rows];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=abs(M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const v=abs(M[i][k]);if(greater(v,best)){p=i;q=k;best=v}}
  if(!greater(best,0))throw Error('The Schur complement has a singular fixed pivot');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]]}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]]}
  for(let i=j+1;i<N;i++){const ratio=div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]))}
 }
 return{dimension:N,evaluate:z=>{
  const raw=get(z),A=rows.map(i=>cols.map(j=>div(div(raw[i][j],scales[i]),cs[j])));
  for(let j=0;j<N-1;j++)for(let i=j+1;i<N;i++){const ratio=div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=sub(A[i][k],mul(ratio,A[j][k]))}
  return A[N-1][N-1];
 }};
}

function predict(initial,n,r,s){
 try{
  const seed=[...initial],reduced=reduceModular([n,0,0,0],seed,r,s),[a,b,c,d]=reduced.matrix;
  const get=z=>tvMatrix([n,0,0,0],div(add(mul(a,z),b),add(mul(c,z),d)),reduced.r,reduced.s,{reduce:false}).M;
  const schur=fixedSchur({abs,div,mul,sub},get,seed);
  let tau=seed;
  for(let iteration=0;iteration<40;iteration++){
   const value=schur.evaluate(tau),h=1e-5;
   const derivative=div(sub(schur.evaluate(add(tau,[h,0])),schur.evaluate(sub(tau,[h,0]))),[2*h,0]);
   let step=div(value,derivative);if(!Number.isFinite(abs(step)))return null;
   if(abs(step)<3e-10)return tau;
   if(abs(step)>.06)step=mul(step,.06/abs(step));
   let accepted=false;
   for(let back=0;back<16;back++){
    const trial=sub(tau,mul(step,2**-back));
    if(!(trial[0]>0&&trial[0]<1&&trial[1]>0&&abs(sub(trial,.5))>.5)||abs(sub(trial,seed))>.24)continue;
    if(abs(schur.evaluate(trial))<abs(value)){tau=trial;accepted=true;break;}
   }
   if(!accepted)return null;
  }
  return null;
 }catch{return null;}
}

return{predict};})();
function transport(r,s,counts,previous,data,id){
 const started=performance.now(),n=1+counts.reduce((sum,x)=>sum+x,0),key=counts.join(',');
 const entry=data.anchors[key];
 if(!entry)return{id,n,status:'missing_anchor',message:'This labelled sheet has no geometric anchor.'};
 const p=createPrecisionEngine(45),v=createPrecisionEngine(70);
 const warm=previous?.status==='ok'&&previous?.result?.status==='ok'&&previous.counts?.every((x,i)=>x===counts[i]);
 const originR=warm?previous.r:data.r,originS=warm?previous.s:data.s;
 let tau=p.asComplex(warm?(previous.result.tauExact??previous.result.tau):(entry.highPrecisionTau??entry.tau));
 const dr=p.d(counts[0]+counts[2]).div(2),ds=p.d(counts[1]+counts[2]).div(2);
 const distance=Math.hypot(r-originR,s-originS);let fraction=0,increment=distance?Math.min(1,.02/distance):1,steps=0,iterations=0,last,priorTau=null,priorFraction=null;
 for(let attempts=0;fraction<1&&attempts<200;attempts++){
  if(performance.now()-started>4000)return{id,n,status:'time_budget',continuationSteps:steps,message:'The point budget was reached. The same labelled sheet can be retried.'};
  const next=Math.min(1,fraction+increment),rr=p.d(String(originR)).plus(p.d(String(r)).minus(String(originR)).times(String(next))).plus(dr),ss=p.d(String(originS)).plus(p.d(String(s)).minus(String(originS)).times(String(next))).plus(ds);
  let seed=tau;
  if(priorTau&&fraction>priorFraction){
   const ratio=(next-fraction)/(fraction-priorFraction);let movement=p.K(p.S(tau,priorTau),p.d(String(ratio)));
   if(ratio>0&&ratio<=3&&p.abs(movement).lt('.12'))seed=p.A(tau,movement);
  }
  if(distance*(next-fraction)>.0005){
   const proposed=schurPredictor.predict(seed.map(x=>x.toNumber()),n,rr.toNumber(),ss.toNumber());
   if(proposed&&p.abs(p.S(p.asComplex(proposed),tau)).lt('.16'))seed=p.asComplex(proposed);
  }
  last=p.correct(seed,rr,ss,n,{tolerance:'1e-18',trustRadius:.23,maxIterations:20});
  iterations+=last.iterations??0;
  if(last.status==='ok'){
   const movement=p.abs(p.S(last.tau,tau)).toNumber();
   if(movement>.16&&increment>1e-5){increment/=2;continue;}
   priorTau=tau;priorFraction=fraction;tau=last.tau;fraction=next;steps++;increment=Math.min(1-fraction,increment*(last.iterations>7?.8:1.4));
  }else{increment/=2;if(increment<1e-6)break;}
 }
 if(fraction!==1||last?.status!=='ok')return{id,n,status:last?.status??'unresolved',continuationSteps:steps,iterations,message:'This labelled branch did not pass the correction checks; no other root was substituted.'};
 const vr=v.d(String(r)).plus(dr.toString()),vs=v.d(String(s)).plus(ds.toString());
 const check=v.correct(tau.map(x=>x.toString()),vr,vs,n,{tolerance:'1e-28',trustRadius:.025,maxIterations:16});
 if(check.status!=='ok')return{id,n,status:'precision_limited',message:'The independent higher-precision point correction did not converge.'};
 const difference=v.abs(v.S(check.tau,v.asComplex(tau.map(x=>x.toString()))));
 if(difference.gt('2e-9'))return{id,n,status:'precision_limited',comparisonTauError:difference.toNumber(),message:'The two working precisions did not agree.'};
 const estimate=v.D.max(difference,check.correction,v.d('1e-30'));
 return{id,n,status:'ok',counts:[...counts],sheetKey:key,target:{r:r+dr.toNumber(),s:s+ds.toNumber()},tau:{re:check.tau[0].toNumber(),im:check.tau[1].toNumber()},tauExact:check.tau.map(x=>x.toPrecision(55)),correctionEstimate:estimate.toNumber(),comparisonTauError:difference.toNumber(),numericalErrorEstimate:estimate.toString(),iterations:iterations+check.iterations,continuationSteps:steps,method:'marked_continuation/scaled_multiprecision_45_70',validation:'two_precisions_not_interval_certified'};
}

return{transport};
}

function createAttachmentEngine(classicalPrecision) {
const {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh} = atlasMath;
/* Validated finite-range attachment engine. Geometry chooses the sheet;
 * numerical continuation transports it. No global root reseeding is used. */

// Smith division avoids squaring the tiny recurrence denominators. Such
// squares can underflow long before the quotient itself loses precision.
function div(a,b){
 if(Math.abs(b.re)>=Math.abs(b.im)){
  const q=b.im/b.re,d=b.re+b.im*q;
  return c((a.re+a.im*q)/d,(a.im-a.re*q)/d);
 }
 const q=b.re/b.im,d=b.im+b.re*q;
 return c((a.re*q+a.im)/d,(a.im*q-a.re)/d);
}

const MAX_N = 21;
const MAX_ATTACHMENTS = MAX_N - 1;
const DEFAULT_SOURCE={r:.31,s:.29};
const BASE_TAU=c(.5904953761444916,.7963943351098978);
let seedData=null;
const anchorCache=new Map();
const finite=z=>z&&Number.isFinite(z.re)&&Number.isFinite(z.im);
const zcopy=z=>Array.isArray(z)?c(z[0],z[1]):c(z.re,z.im);
const keyOf=counts=>counts.join(',');
const sameCounts=(a,b)=>Array.isArray(a)&&Array.isArray(b)&&a.length===3&&a.every((v,i)=>v===b[i]);
const pow=(z,n)=>{if(n<0)return div(c(1),pow(z,-n));let out=c(1),base=z;while(n){if(n&1)out=mul(out,base);base=mul(base,base);n=Math.floor(n/2);}return out;};

function validateCounts(counts,anchorData=seedData){
 if(!Array.isArray(counts)||counts.length!==3||!counts.every(v=>Number.isInteger(v)&&v>=0))
  return {status:'invalid_counts',message:'Enter three nonnegative whole numbers.'};
 const n=1+counts.reduce((a,b)=>a+b,0);
 const available=Math.min(MAX_N,anchorData?.maxN??anchorData?.max_n_completed??MAX_N);
 if(n>available)return {status:'unsupported_index',n,maxN:available,
  message:`This release validates at most ${available-1} attachments in total.`};
 return {status:'ok',n,maxN:available};
}

function setAnchorData(data){
 if(!data?.anchors||!Number.isFinite(data.r)||!Number.isFinite(data.s))
  throw new TypeError('Attachment anchors must include r, s and an anchors dictionary.');
 seedData=data;anchorCache.clear();
 return {status:'ok',maxN:Math.min(MAX_N,data.maxN??data.max_n_completed??MAX_N),
  anchorCount:Object.keys(data.anchors).length};
}

function ellipticBasis(tau,r,s){
 const e=evaluate(tau,r,s,1);
 let e1=evaluate(tau,.5,0,1).wp,e2=evaluate(tau,0,.5,1).wp,e3=evaluate(tau,.5,.5,1).wp;
 // The premodular polynomial is symmetric in the three half-period values.
 // Choose the largest Legendre denominator; a nearly coalescing pair makes
 // the rational recurrence grow although the underlying polynomial is finite.
 const d12=abs(sub(e2,e1)),d13=abs(sub(e3,e1)),d23=abs(sub(e3,e2));
 if(d13>=d12&&d13>=d23)[e2,e3]=[e3,e2];
 else if(d23>d12)[e1,e2,e3]=[e2,e3,e1];
 const D=sub(e2,e1),t=div(sub(e3,e1),D);
 const natural=Math.max(1,abs(e.z),Math.sqrt(abs(e.wp)),Math.cbrt(abs(e.wpPrime)),Math.sqrt(abs(D)));
 // The equation has weight n(n+1)/2. Normalize its weight-one scale BEFORE
 // the recurrence: large intermediate H and squared division denominators
 // otherwise overflow even while the final modulus is well conditioned.
 // This multiplies the equation by a nonzero real factor, preserving roots.
 const normalizedSize=1,inv=normalizedSize/natural,inv2=inv*inv;
 return {...e,z:scale(e.z,inv),wp:scale(e.wp,inv2),wpPrime:scale(e.wpPrime,inv2*inv),
  e1:scale(e1,inv2),e2:scale(e2,inv2),e3:scale(e3,inv2),D:scale(D,inv2),t,natural:normalizedSize,
  weightScale:natural/normalizedSize};
}
function leading(j,D){const k=j%2===0?j*(j+2)/4:(j+1)**2/4;return scale(pow(D,-k),2**(-2*k));}

function recurrence(b,n){
 const Q={'-2':c(1),'-1':c(1),0:b.z};
 let G=div(b.wpPrime,b.D),R=add(mul(div(sub(b.wp,b.e1),b.D),b.z),scale(G,.5));
 let suspect=false,minRelativePrior=Math.min(1,abs(b.z)/b.natural);
 for(let j=1;j<n;j++){
  const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1];
  const phi=add(mul(R,q2),scale(G,(j-1)/2));
  const product=mul(mul(q3,q2),q1),u=sub(phi,product),v=sub(phi,mul(b.t,product));
  const H=mul(mul(phi,u),v),Hp=add(add(mul(phi,u),mul(phi,v)),mul(u,v));
  const nextG=div(H,mul(G,pow(q3,3)));
  const nextQ=sub(div(mul(q3,nextG),G),scale(div(Hp,mul(G,pow(q3,2))),j/2));
  Q[j]=nextQ;
  if(j<n-1){
   const first=mul(phi,nextQ),second=scale(mul(q3,nextG),(j+1)/2),numerator=add(first,second);
   const cancellation=abs(numerator)/Math.max(1e-300,abs(first)+abs(second));
   if(cancellation<1e-10&&minRelativePrior<1e-7)suspect=true;
   R=div(numerator,mul(q3,q1));
   const degree=(j+1)*(j+2)/2;
   const typical=abs(leading(j,b.D))*b.natural**degree;
   minRelativePrior=Math.min(minRelativePrior,abs(nextQ)/Math.max(1e-300,typical));
  }
  G=nextG;
 }
 if(n>=3&&minRelativePrior<1e-11)suspect=true;
 const value=div(Q[n-1],leading(n-1,b.D));
 return {value,suspect: suspect||!finite(value),minRelativePrior};
}

function explicitTerms(b,n){
 const Z=b.z,p=b.wp,pp=b.wpPrime,p2=mul(p,p),pp2=mul(pp,pp);
 const g2=scale(add(add(mul(b.e1,b.e2),mul(b.e1,b.e3)),mul(b.e2,b.e3)),-4);
 const g3=scale(mul(mul(b.e1,b.e2),b.e3),4);
 if(n===3)return[pow(Z,6),scale(mul(p,pow(Z,4)),-15),scale(mul(pp,pow(Z,3)),-20),
  mul(sub(scale(g2,27/4),scale(p2,45)),pow(Z,2)),scale(mul(mul(p,pp),Z),-12),scale(pp2,-5/4)];
 if(n===4)return[pow(Z,10),scale(mul(p,pow(Z,8)),-45),scale(mul(pp,pow(Z,7)),-120),
  mul(sub(scale(g2,399/4),scale(p2,630)),pow(Z,6)),scale(mul(mul(p,pp),pow(Z,5)),-504),
  scale(mul(sub(sub(scale(pow(p,3),280),scale(mul(g2,p),49)),scale(g3,115)),pow(Z,4)),-15/4),
  scale(mul(mul(sub(scale(g2,11),scale(p2,24)),pp),pow(Z,3)),15),
  scale(mul(add(add(sub(scale(pow(p,4),140),scale(mul(g2,p2),245)),scale(mul(g3,p),190)),scale(mul(g2,g2),21)),pow(Z,2)),-9/4),
  scale(mul(mul(add(sub(scale(pow(p,3),40),scale(mul(g2,p),163)),scale(g3,125)),pp),Z),-1),
  scale(mul(sub(scale(g2,25),scale(p2,3)),pp2),3/4)];
 return null;
}
function safelyAwayFromLattice(tau,r,s){
 const rr=r-Math.round(r),ss=s-Math.round(s);
 let shortest=Math.min(1,abs(tau),abs(sub(tau,c(1)))),distance=Infinity;
 for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)distance=Math.min(distance,abs(add(c(rr-i),scale(tau,ss-j))));
 return distance>.2*shortest;
}

/** Polynomial projection cancels removable recurrence singularities. It is a
 * fallback, not a universal replacement for the pointwise recurrence. */
function project(b,n,radiusRatio=.025){
 const degree=n*(n+1)/2,M=2**Math.ceil(Math.log2(degree+1)),radius=radiusRatio*b.natural;
 let sr=0,si=0,cr=0,ci=0,maxNode=0;
 for(let k=0;k<M;k++){
  const angle=2*Math.PI*(k+.1732050807568877)/M;
  const x=add(b.z,c(radius*Math.cos(angle),radius*Math.sin(angle)));
  const value=recurrence({...b,z:x},n).value;
  if(!finite(value))throw new RangeError('A polynomial-projection node lost precision.');
  maxNode=Math.max(maxNode,abs(value));
  let y=value.re-cr,t=sr+y;cr=(t-sr)-y;sr=t;
  y=value.im-ci;t=si+y;ci=(t-si)-y;si=t;
 }
 return {value:c(sr/M,si/M),maxNode,nodes:M,radius};
}

function evaluateAttachmentEquation(tau,r,s,n,options={}){
 if(!Number.isInteger(n)||n<1||n>MAX_N)throw new RangeError(`Supported indices are one through ${MAX_N}.`);
 tau=zcopy(tau);
 if(n<=2){const e=evaluate(tau,r,s,n);return {...e,method:n===1?'Hecke':'pole_subtracted_explicit',normalizer:e.residual/Math.max(e.relativeResidual,1e-300)};}
 const b=ellipticBasis(tau,r,s),raw=recurrence(b,n),degree=n*(n+1)/2;
 let value=raw.value,method='CKL_recurrence',normalizer=Math.max(1,b.natural**degree),projection=null;
 const explicit=explicitTerms(b,n),safe=safelyAwayFromLattice(tau,r,s);
 if(options.mode==='projection'){
  projection=project(b,n,options.radiusRatio??.025);value=projection.value;method='polynomial_projection';
 }else if(raw.suspect){
  if(explicit&&safe){value=explicit.reduce(add,c());normalizer=Math.max(1,explicit.reduce((a,v)=>a+abs(v),0));method='explicit_removed_division';}
  else{projection=project(b,n);value=projection.value;method='polynomial_projection';}
 }
 if(!finite(value))throw new RangeError('The equation could not be evaluated reliably.');
 return {value,residual:abs(value),relativeResidual:abs(value)/normalizer,normalizer,
  method,projectionNodes:projection?.nodes??0,recurrenceCondition:raw.minRelativePrior};
}

function differential(tau,r,s,n){
 const h=1e-5*Math.max(1,abs(tau));
 const plus=evaluateAttachmentEquation(c(tau.re+h,tau.im),r,s,n).value;
 const minus=evaluateAttachmentEquation(c(tau.re-h,tau.im),r,s,n).value;
 return scale(sub(plus,minus),1/(2*h));
}
function finalCheck(tau,r,s,n,e,derivative){
 if(n<=2)return{status:'ok',comparisonTauError:0};
 const b=ellipticBasis(tau,r,s),terms=explicitTerms(b,n);
 let difference;
 if(terms&&safelyAwayFromLattice(tau,r,s))difference=abs(sub(e.value,terms.reduce(add,c())));
 else{
  const a=project(b,n,.0125),d=project(b,n,.025);
  difference=Math.max(abs(sub(e.value,a.value)),abs(sub(a.value,d.value)));
 }
 const comparisonTauError=difference/Math.max(1e-300,abs(derivative));
 return {status:comparisonTauError<2e-8?'ok':'precision_limited',comparisonTauError};
}

function correct(tau,r,s,n,trustRadius=.2){
 tau=zcopy(tau);const seed=zcopy(tau);let evaluations=0;
 try{
  if(!inF0(tau))return{status:'outside_chart',tau,iterations:0};
  let e=evaluateAttachmentEquation(tau,r,s,n);evaluations++;
  for(let iteration=0;iteration<28;iteration++){
   const derivative=differential(tau,r,s,n);evaluations+=2;
   if(!finite(derivative)||abs(derivative)<1e-250)return{status:'small_derivative',tau,iterations:iteration};
   const step=div(e.value,derivative),correctionEstimate=abs(step);
   if(correctionEstimate<2e-10){
    const checked=finalCheck(tau,r,s,n,e,derivative);
    return {...e,...checked,tau,iterations:iteration,evaluations,correctionEstimate};
   }
   let accepted=false;
   for(let back=0;back<22;back++){
    const trial=sub(tau,scale(step,2**-back));
    if(!inF0(trial)||abs(sub(trial,seed))>trustRadius)continue;
    const next=evaluateAttachmentEquation(trial,r,s,n);evaluations++;
    if(next.residual<e.residual*(1-.00001*2**-back)){
     tau=trial;e=next;accepted=true;break;
    }
   }
   if(!accepted)return{status:'correction_failed',tau,iterations:iteration+1,evaluations};
  }
  return{status:'iteration_limit',tau,iterations:28,evaluations};
 }catch(error){return{status:'precision_limited',tau,iterations:0,evaluations,message:error.message};}
}

function initialAnchor(counts,data){
 const key=keyOf(counts);
 if(anchorCache.has(key))return anchorCache.get(key);
 if(counts.every(x=>x===0))return{r:DEFAULT_SOURCE.r,s:DEFAULT_SOURCE.s,tau:BASE_TAU,status:'validated_anchor'};
 const entry=data?.anchors?.[key];
 if(!entry)return null;
 const tau=entry.tau??entry.analyticTau??entry.fine?.tau;
 if(!tau)return null;
 const answer={r:data.r,s:data.s,tau:zcopy(tau),status:entry.status??'geometric_anchor'};
 anchorCache.set(key,answer);return answer;
}

function transport(r,s,counts,previous,data,id){
 const transportStarted=globalThis.performance?.now?.()??Date.now();
 const n=1+counts.reduce((a,b)=>a+b,0),dr=(counts[0]+counts[2])/2,ds=(counts[1]+counts[2])/2;
 if(n>12)return classicalPrecision.transport(r,s,counts,previous,data,id);
 const prior=id==='base'?previous?.base:(sameCounts(previous?.counts,counts)?previous?.result:null);
 const anchor=prior?.status==='ok'&&sourceStatus(previous.r,previous.s)==='ok'?
  {r:previous.r,s:previous.s,tau:prior.tau}:initialAnchor(counts,data);
 const target={r:r+dr,s:s+ds};
 if(!anchor)return{id,n,target,status:'missing_anchor',message:'The geometric sheet anchor is unavailable.'};
 let tau=zcopy(anchor.tau),fraction=0,steps=0,iterations=0,evaluations=0,last=null;
 const distance=Math.hypot(r-anchor.r,s-anchor.s);
 let increment=distance?Math.min(1,.022/distance):1;
 for(let attempts=0;fraction<1&&attempts<220;attempts++){
  if((globalThis.performance?.now?.()??Date.now())-transportStarted>350){
   last={status:'precision_limited'};break;
  }
  const next=Math.min(1,fraction+increment),rr=anchor.r+next*(r-anchor.r),ss=anchor.s+next*(s-anchor.s);
  last=correct(tau,rr+dr,ss+ds,n,Math.min(.23,Math.max(.035,.8*abs(tau))));
  iterations+=last.iterations??0;evaluations+=last.evaluations??0;
  if(last.status==='ok'){
   const movement=abs(sub(last.tau,tau));
   if(movement>.16&&increment>3e-5){increment*=.5;continue;}
   tau=last.tau;fraction=next;steps++;
   increment=Math.min(1-fraction,increment*(last.iterations>7?.8:1.35));
  }else{increment*=.5;if(increment<1e-5)break;}
 }
 if(fraction!==1||last?.status!=='ok')return{id,n,target,status:last?.status??'unresolved',
  iterations,evaluations,continuationSteps:steps,message:'This sheet could not be transported reliably; no other root was substituted.'};
 return {...last,id,n,target,tau,iterations,evaluations,continuationSteps:steps,
  counts:[...counts],sheetKey:keyOf(counts),method:`marked_continuation/${last.method}`};
}

/** Input counts refer to CWW's clockwise edge labels L1,L2,L3. */
function solveAttachment(r,s,counts,previousState=null,anchorData=null){
 const started=globalThis.performance?.now?.()??Date.now();
 const elapsed=()=> (globalThis.performance?.now?.()??Date.now())-started;
 if(anchorData&&anchorData!==seedData)setAnchorData(anchorData);
 const checked=validateCounts(counts,seedData),source=sourceStatus(r,s);
 if(checked.status!=='ok'||source!=='ok')return{r,s,counts:Array.isArray(counts)?[...counts]:counts,
  n:checked.n,status:checked.status!=='ok'?checked.status:source,branches:[],elapsedMs:elapsed(),
  message:checked.message??'Choose a point slightly inside the source triangle.'};
 let base=transport(r,s,[0,0,0],previousState,seedData,'base');
 if(base.status!=='ok'&&previousState?.base?.status==='ok')
  base={...transport(r,s,[0,0,0],null,seedData,'base'),restartedFromAnchor:true};
 let result=counts.every(x=>x===0)?{...base,id:'attached',counts:[0,0,0],sheetKey:'0,0,0'}:
  transport(r,s,counts,previousState,seedData,'attached');
 // A failed long warm path may restart only from this same sheet's fixed
 // geometric anchor. It must never borrow another sheet's analytic root.
 if(result.status!=='ok'&&sameCounts(previousState?.counts,counts)&&previousState?.result?.status==='ok')
  result={...transport(r,s,counts,null,seedData,'attached'),restartedFromAnchor:true};
 return{r,s,counts:[...counts],n:checked.n,status:base.status==='ok'&&result.status==='ok'?'ok':'partial',
  base,result,branches:[base,result],maxN:checked.maxN,elapsedMs:elapsed(),
  validation:{kind:'numerical_not_interval_certified',sourceMargin:.004,
   equationDegree:checked.n*(checked.n+1)/2,geometricAnchorKey:keyOf(counts)}};
}

async function sampleAttachment(options={}){
 if(options.anchorData)setAnchorData(options.anchorData);
 const counts=options.counts??[0,0,0],checked=validateCounts(counts),meshStarted=performance.now();
 if(checked.status!=='ok')return{status:checked.status,message:checked.message,counts,points:[],triangles:[],samples:[]};
 const mesh=buildSourceMesh(options.resolution??16,options.inset??.025),samples=new Array(mesh.points.length);
 // Alternate rows so continuation never jumps across the source triangle.
 const order=[];let cursor=0;
 for(let row=0;row<=mesh.resolution;row++){
  const indices=Array.from({length:mesh.resolution-row+1},()=>cursor++);
  order.push(...(row%2?indices.reverse():indices));
 }
 let previous=options.previousState??null,completed=0;
 for(const index of order){
  if(performance.now()-meshStarted>45000)break;
  if(options.signal?.aborted)return{...mesh,counts:[...counts],n:checked.n,samples,status:'aborted',completed};
  const p=mesh.points[index],sample=solveAttachment(p.r,p.s,counts,previous);
  samples[index]=sample;completed++;
  if(sample.status==='ok')previous=sample;
  if(options.onSample)await options.onSample(sample,index,completed,mesh.points.length);
  // Every point yields to queued interactive jobs in a worker.
  await new Promise(resolve=>setTimeout(resolve,0));
 }
 return{...mesh,counts:[...counts],n:checked.n,samples,completed,
  status:completed===mesh.points.length&&samples.every(s=>s?.status==='ok')?'ok':'partial',
  stopReason:completed<mesh.points.length?'time_budget':null,
  message:completed<mesh.points.length?'The 45-second mesh budget was reached. Checked points are retained; choose the quick mesh for a complete preview.':null,
  metadata:{sampledInteriorOnly:true,sourceCoordinateSideGap:mesh.inset/2,
   maximumRelativeResidual:Math.max(...samples.flatMap(s=>s.branches.filter(b=>b.status==='ok'&&Number.isFinite(b.relativeResidual)).map(b=>b.relativeResidual)))}};
}

return {MAX_N,validateCounts,setAnchorData,solveAttachment,sampleAttachment};
}

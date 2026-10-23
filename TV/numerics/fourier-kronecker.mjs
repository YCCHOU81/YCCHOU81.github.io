// Research-only, cancellation-free deep infinity cusp jets for connector
// families.  No website imports.  Same exact Bloch determinant as generic
// theta quotient, after the known isogeny; odd trigonometric jets removed
// analytically before evaluating the even derivative constraints.
import Decimal from '../decimal.mjs';
import {createFixed}from './bigint-fixed.mjs';

export function createFourierEvaluator(precision=80,{bigint=false,modesMultiplier=1}={}){
 const D=bigint?createFixed(precision):Decimal.clone({precision,rounding:Decimal.ROUND_HALF_EVEN,maxE:9e15,minE:-9e15});
 const Z=new D(0),O=new D(1),PI=D.acos(-1),TP=PI.times(2),HALF=new D('.5');
 const c=x=>Array.isArray(x)?x.map(v=>new D(v)):[new D(x),Z];
 const add=(a,b)=>{a=c(a);b=c(b);return[a[0].plus(b[0]),a[1].plus(b[1])]};
 const neg=a=>{a=c(a);return a.map(v=>v.neg())};
 const sub=(a,b)=>add(a,neg(b));
 const mul=(a,b)=>{a=c(a);b=c(b);return[a[0].times(b[0]).minus(a[1].times(b[1])),a[0].times(b[1]).plus(a[1].times(b[0]))]};
 const div=(a,b)=>{a=c(a);b=c(b);if(b[0].abs().gte(b[1].abs())){const ratio=b[1].div(b[0]),den=b[0].plus(b[1].times(ratio));return[a[0].plus(a[1].times(ratio)).div(den),a[1].minus(a[0].times(ratio)).div(den)]}const ratio=b[0].div(b[1]),den=b[1].plus(b[0].times(ratio));return[a[0].times(ratio).plus(a[1]).div(den),a[1].times(ratio).minus(a[0]).div(den)]};
 const scale=(a,b)=>mul(a,b),abs=a=>{a=c(a);const hi=D.max(a[0].abs(),a[1].abs()),lo=D.min(a[0].abs(),a[1].abs());return hi.isZero()?Z:hi.times(O.plus(lo.div(hi).pow(2)).sqrt())};
 const exp=a=>{a=c(a);const v=a[0].exp();return[v.times(a[1].cos()),v.times(a[1].sin())]};
 const phase=a=>exp([Z,TP.times(a)]);
 const ipow=k=>[[O,Z],[Z,O],[O.neg(),Z],[Z,O.neg()]][((k%4)+4)%4];
 function factorial(n){let v=O;for(let j=2;j<=n;j++)v=v.times(j);return v}
 function binomial(n,k){let v=O;for(let j=1;j<=k;j++)v=v.times(n-k+j).div(j);return v}
 // v/(2 sin(v/2)) = sum b_(2j) v^(2j). This reciprocal has no
 // small-nome or small-character subtraction.
 const cscCache=new Map();
 function cscCoeffs(N){
  if(cscCache.has(N))return cscCache.get(N);
  const sinc=[],b=[O];
  for(let j=0;j<=N;j++)sinc.push(new D(j%2?-1:1).div(new D(2).pow(2*j).times(factorial(2*j+1))));
  for(let j=1;j<=N;j++){let v=Z;for(let k=1;k<=j;k++)v=v.minus(sinc[k].times(b[j-k]));b.push(v)}
  cscCache.set(N,b);return b;
 }
 function matrix(family,n,w,delta,S){
  w=c(w);delta=new D(delta);S=new D(S);const tau=scale(w,2),s=HALF.minus(delta.div(2)),epsilon=delta.div(-2),r=S.plus(new D(n-1).div(2));
  const qCache=new Map();
  const qpow=p=>{p=new D(p);const key=p.toString();if(!qCache.has(key))qCache.set(key,exp([TP.times(tau[1]).times(p).neg(),TP.times(tau[0]).times(p)]));return qCache.get(key)};
  const degree=2*(n-1),fact=[O],b=cscCoeffs(n);
  for(let k=1;k<=degree+1;k++)fact.push(fact[k-1].times(k));
  const powerCache=new Map();
  const powers=t=>{const key=t.toString();if(!powerCache.has(key)){const out=[O];for(let d=1;d<=degree;d++)out.push(out[d-1].times(t));powerCache.set(key,out)}return powerCache.get(key)};
  const spow=powers(s);
  const phaseR=phase(r),E=mul(phaseR,qpow(s)),phaseInv=phase(r.neg()),epsPower=[O];
  for(let j=1;j<=degree+1;j++)epsPower.push(epsPower[j-1].times(epsilon));
  const digitPerNome=TP.times(tau[1]).toNumber()/Math.LN10;
  // Every omitted term is bounded geometrically; two extra modes cover
  // degree-dependent polynomial factors. Conservative for tested n<=65.
  const modes=Math.max(4,Math.ceil(2*(precision+20+degree)/digitPerNome)+4)*modesMultiplier;
  const es=div(E,sub(1,E)),modesData=[],crossData=[],baseCross=[];
  for(let m=1;m<=modes;m++){
   const am=mul(phaseR,qpow(s.plus(m))),bm=mul(phaseInv,qpow(new D(m).minus(s)));
   modesData.push([m,div(am,sub(1,am)),div(bm,sub(1,bm)),powers(s.plus(m)),powers(s.minus(m))]);
   if(family==='tetra'){
    const a=div(mul(phaseR,qpow(s.plus(m).times('1.5'))),sub(1,am));
    const b=div(mul(phaseInv,qpow(new D(m).minus(s).div(2))),sub(1,bm));
    crossData.push([m,a,b,powers(s.plus(m)),powers(s.minus(m))]);
   }
  }
  if(family==='tetra')for(let j=0;j<=modes;j++)baseCross.push([j,qpow(s.plus(j).div(2)),powers(s.plus(j))]);
  const crossConstant=family==='tetra'?div(mul(phaseR,qpow(s.times('1.5'))),sub(1,E)):null;
  function evenJet(d){
   let base=c(0);
   for(let j=0;j<=d/2;j++){const k=d+1-2*j;const coefficient=b[j].times(epsPower[k]).div(fact[k]);base=add(base,mul(ipow(k),coefficient))}
   let tail=scale(es,spow[d]);
   for(const [m,a,b,pp,pm]of modesData)tail=add(tail,sub(scale(a,pp[d]),scale(b,pm[d])));
   // Derivatives, rather than Taylor coefficients, avoid factorial-sized
   // underflow. This is a harmless row/column rescaling of the jet matrix.
   return add(scale(base,fact[d]),mul(ipow(d-1),tail));
  }
  function crossJet(d){
   let sum=c(0);
   for(const[j,z,pp]of baseCross)sum=add(sum,scale(z,pp[d]));
   sum=add(sum,scale(crossConstant,spow[d]));
   for(const[m,a,b,pp,pm]of crossData){
    sum=add(sum,sub(scale(a,pp[d]),scale(b,pm[d])));
   }
   return mul(ipow(d-1),sum);
  }
  const orders=[];for(let h=n-1;h>=0;h-=2)orders.push(h);
  const even=new Map(),cross=new Map();
  const ev=d=>{if(!even.has(d))even.set(d,evenJet(d));return even.get(d)};
  const cr=d=>{if(!cross.has(d))cross.set(d,crossJet(d));return cross.get(d)};
  const M=orders.map(l=>orders.map(h=>ev(l+h)));
  if(family==='tetra'){
   for(let i=0;i<orders.length;i++)M[i].push(mul(phaseR,cr(orders[i])));
   M.push([...orders.map(h=>cr(h)),ev(0)]);
  }
  return {M,modes,dimension:M.length};
 }
 function determinant(family,n,w,delta,S){
  const {M,modes,dimension}=matrix(family,n,w,delta,S);let value=c(1),minPivot=O;
  // A scalar row is the vanishing residual. Normalizing it would erase its
  // magnitude and prevent the n=1 and n=2 double-seed secant from converging.
  if(dimension>1)for(const row of M){const norm=D.max(...row.map(abs));for(let j=0;j<row.length;j++)row[j]=div(row[j],norm)}
  for(let j=0;j<dimension;j++){
   let p=j;for(let k=j+1;k<dimension;k++)if(abs(M[k][j]).gt(abs(M[p][j])))p=k;
   if(p!==j){[M[p],M[j]]=[M[j],M[p]];value=neg(value)}
   const pivot=M[j][j];minPivot=D.min(minPivot,abs(pivot));value=mul(value,pivot);
   if(abs(pivot).isZero())return{value:c(0),minPivot,dimension,modes};
   // The product of harmless nonzero pivots may be 10^(-O(n^2)).
   // Remove that real positive scale before fixed-point underflow.
   // Zero locations and residual phases are unchanged. As with row
   // normalization, this residual is smooth but is not holomorphic.
   if(j<dimension-1)value=div(value,abs(value));
   for(let i=j+1;i<dimension;i++){const a=div(M[i][j],pivot);for(let k=j+1;k<dimension;k++)M[i][k]=sub(M[i][k],mul(a,M[j][k]))}
  }
  return{value,minPivot,dimension,modes};
 }
 function solve(family,n,delta,S,guess,{tol=new D(10).pow(-precision+15),maxsteps=30}={}){
  let a=c(guess),b=add(a,['.000001','.000002']),fa=determinant(family,n,a,delta,S).value,fb=determinant(family,n,b,delta,S).value,evaluations=2;
  for(let k=0;k<maxsteps;k++){
   if(abs(fb).isZero())throw Error('Zero residual at finite precision; repeat at higher precision');
   const step=div(mul(fb,sub(b,a)),sub(fb,fa));
   if(abs(step).lt(tol))return{w:b,correction:b[0].minus(family==='double'?'.25':'.5'),iterations:k+2,evaluations,residual:abs(fb)};
   a=b;fa=fb;b=sub(b,step);if(b[1].lte(0))throw Error('Outside upper half-plane');fb=determinant(family,n,b,delta,S).value;evaluations++;
  }
  throw Error('No convergence');
 }
 function prediction(family,n,delta,S){const P=Decimal.clone({precision:Math.min(precision,1000)}),pi=P.acos(-1),C=family==='double'?new P(2*n*(n+1)).times(pi.times(2).times(String(S)).sin()):new P(8*n*(n+1)).sqrt().times(pi.times(String(S)).sin());return[new D(family==='double'?'.25':'.5'),new D(C.div(String(delta)).ln().div(pi).times(family==='double'?'.5':1).toString())]}
 return{D,c,add,sub,mul,div,abs,exp,matrix,determinant,solve,prediction,strs:a=>a.map(v=>v.toString())};
}

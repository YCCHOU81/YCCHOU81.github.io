/* Shared stable Fourier/Decimal evaluator. Each owner keeps its own precision cache.
 * reuseAcceptedEvaluation preserves the existing Overlap optimization; it is
 * deliberately off for Explore and Limits, whose arithmetic is unchanged. */
function createLamePrecisionEngine({reuseAcceptedEvaluation=false}={}) {
const precisionEngineCache=new Map();

function createPrecisionEngine(precision=75){
 if(precisionEngineCache.has(precision))return precisionEngineCache.get(precision);
 const D=Decimal.clone({precision,rounding:Decimal.ROUND_HALF_EVEN}),d=x=>new D(x);
 const C=(x=0,y=0)=>[d(x),d(y)],A=(a,b)=>[a[0].plus(b[0]),a[1].plus(b[1])],S=(a,b)=>[a[0].minus(b[0]),a[1].minus(b[1])];
 const M=(a,b)=>[a[0].times(b[0]).minus(a[1].times(b[1])),a[0].times(b[1]).plus(a[1].times(b[0]))];
 const K=(a,k)=>[a[0].times(k),a[1].times(k)],norm=a=>a[0].times(a[0]).plus(a[1].times(a[1])),abs=a=>norm(a).sqrt();
 const V=(a,b)=>K([a[0].times(b[0]).plus(a[1].times(b[1])),a[1].times(b[0]).minus(a[0].times(b[1]))],d(1).div(norm(b)));
 // Preserve the two existing arithmetic paths, including Overlap's unused-square optimization.
 const P=reuseAcceptedEvaluation ? (z,n)=>{if(n<0)return V(C(1),P(z,-n));let o=C(1),a=z;while(n){if(n%2)o=M(o,a);n=Math.floor(n/2);if(n)a=M(a,a);}return o;} : (z,n)=>{if(n<0)return V(C(1),P(z,-n));let o=C(1),a=z;while(n){if(n%2)o=M(o,a);a=M(a,a);n=Math.floor(n/2);}return o;};
 const pi=D.acos(-1),twoPi=pi.times(2),one=C(1),threshold=d(10).pow(-precision+12),thresholdSquared=threshold.times(threshold);
 const exp=z=>{const v=z[0].exp();return[v.times(z[1].cos()),v.times(z[1].sin())];};
 const imaginaryExp=(real,imag)=>exp([d(imag).times(twoPi).neg(),d(real).times(twoPi)]);
 const asComplex=z=>Array.isArray(z)?C(z[0],z[1]):C(z.re,z.im);
 const inChart=t=>t[0].gt(0)&&t[0].lt(1)&&t[1].gt(0)&&norm(S(t,C('.5'))).gt('.25');
 function reduce(tau,r,s){
  let t=asComplex(tau),rr=d(r),ss=d(s),factor=one,ma=d(1),mb=d(0),mc=d(0),md=d(1);
  for(let k=0;k<100;k++){
   const integer=t[0].toNearest(1);t=[t[0].minus(integer),t[1]];rr=rr.plus(integer.times(ss));
   ma=ma.minus(integer.times(mc));mb=mb.minus(integer.times(md));
   rr=rr.minus(rr.toNearest(1));ss=ss.minus(ss.toNearest(1));
   if(norm(t).gte(d(1).minus(threshold)))return{t,r:rr,s:ss,factor,matrixC:mc};
   factor=M(t,factor);t=V(C(-1),t);[rr,ss]=[ss,rr.neg()];
   [ma,mb,mc,md]=[mc.neg(),md.neg(),ma,mb];
  }
  throw Error('Could not reduce the elliptic evaluation chart.');
 }
 function basis(tau,r,s,n=8,withDerivative=false){
  const chart=reduce(tau,r,s),{t}=chart; r=chart.r;s=chart.s;
  const ar=r.plus(s.times(t[0])),ai=s.times(t[1]),sign=ai.gte(0)?1:-1;
  const w=imaginaryExp(ar.times(sign),ai.abs()),den=S(w,one),quot=V(A(w,one),den);
  const cot=K([quot[1].neg(),quot[0]],sign),csc=K(V(w,M(den,den)),-4);
  let sin=C(),psum=C(),ppsum=C(),eisen2=C();
  let dsin=C(),dpsum=C(),dppsum=C();
  const timesI=(z,k)=>K([z[1].neg(),z[0]],k);
  const q=imaginaryExp(t[0],t[1]);
  const up=sign===1?M(q,w):V(q,w),down=sign===1?V(q,w):M(q,w);
  let qk=q,uk=up,vk=down;
  for(let k=1;k<400;k++){
   const denominator=S(one,qk),coeff=V(qk,denominator),plus=V(uk,denominator),minus=V(vk,denominator);
   const dif=S(plus,minus),sine=K([dif[1],dif[0].neg()],'.5'),cos=K(A(plus,minus),'.5');
   sin=A(sin,sine);psum=A(psum,K(S(coeff,cos),k));ppsum=A(ppsum,K(sine,k*k));
   eisen2=A(eisen2,K(coeff,k));
   if(withDerivative){
    // Differentiate the convergent Fourier series in the reduced chart.
    // r and s are fixed there; q' = 2 pi i q and exp(2 pi i z)' = 2 pi i s exp(2 pi i z).
    const dcoeff=timesI(V(coeff,denominator),twoPi.times(k));
    const dplus=timesI(M(plus,A(C(d(1).plus(s)),coeff)),twoPi.times(k));
    const dminus=timesI(M(minus,A(C(d(1).minus(s)),coeff)),twoPi.times(k));
    const dsine=K([S(dplus,dminus)[1],S(dplus,dminus)[0].neg()],'.5');
    const dcos=K(A(dplus,dminus),'.5');
    dsin=A(dsin,dsine);dpsum=A(dpsum,K(S(dcoeff,dcos),k));dppsum=A(dppsum,K(dsine,k*k));
   }
   // The equivalent squared tail test avoids three arbitrary-precision
   // square roots per Fourier term without weakening the stopping threshold.
   if(k>5&&D.max(norm(plus),norm(minus),norm(coeff)).times(k**10).lt(thresholdSquared))break;
   if(k===399)throw Error('High-precision Fourier series did not converge.');
   qk=M(qk,q);uk=M(uk,up);vk=M(vk,down);
  }
  const z=A(A(K(cot,pi),K(sin,pi.times(4))),C(0,twoPi.times(s)));
  const wp=A(K(S(csc,C(d(1).div(3))),pi.pow(2)),K(psum,pi.pow(2).times(8)));
  const pp=A(K(M(cot,csc),pi.pow(3).times(-2)),K(ppsum,pi.pow(3).times(16)));
  const eta1=K(S(one,K(eisen2,24)),pi.pow(2).div(3));
  let derivative;
  if(withDerivative){
   const dcot=K(csc,pi.times(s).neg()),dcsc=K(M(cot,csc),twoPi.times(s).neg());
   derivative={z:A(K(dcot,pi),K(dsin,pi.times(4))),
    wp:A(K(dcsc,pi.pow(2)),K(dpsum,pi.pow(2).times(8))),
    pp:A(K(A(M(dcot,csc),M(cot,dcsc)),pi.pow(3).times(-2)),K(dppsum,pi.pow(3).times(16)))};
  }
  if(n<=2)return{z,wp,pp,eta1,factor:chart.factor,matrixC:chart.matrixC,derivative};
  let t3=one,t4=one,t2sum=one;
  let dt3=C(),dt4=C(),dt2sum=C();
  const sqrt=z=>{const radius=abs(z);return[radius.plus(z[0]).div(2).sqrt(),radius.minus(z[0]).div(2).sqrt().times(z[1].isNeg()?-1:1)];};
  const halfq=sqrt(q),eighthq=sqrt(sqrt(halfq));
  for(let k=1;k<40;k++){
   const term=P(halfq,k*k),term2=P(halfq,k*(k+1));t3=A(t3,K(term,2));t4=A(t4,K(term,k%2?-2:2));
   t2sum=A(t2sum,term2);
   if(withDerivative){
    const dt=timesI(term,pi.times(k*k));dt3=A(dt3,K(dt,2));dt4=A(dt4,K(dt,k%2?-2:2));
    dt2sum=A(dt2sum,timesI(term2,pi.times(k*(k+1))));
   }
   if(k>3&&norm(term).lt(thresholdSquared))break;
  }
  const t2=K(M(eighthq,t2sum),2);
  const e1=K(A(P(t3,4),P(t4,4)),pi.pow(2).div(3));
  const e2=K(A(P(t2,4),P(t3,4)),pi.pow(2).div(-3));
  const e3=K(S(P(t2,4),P(t4,4)),pi.pow(2).div(3));
  if(withDerivative){
   const dt2=K(A(M(timesI(eighthq,pi.div(4)),t2sum),M(eighthq,dt2sum)),2);
   const p2=K(M(P(t2,3),dt2),4),p3=K(M(P(t3,3),dt3),4),p4=K(M(P(t4,3),dt4),4);
   derivative.e1=K(A(p3,p4),pi.pow(2).div(3));
   derivative.e2=K(A(p2,p3),pi.pow(2).div(-3));
   derivative.e3=K(S(p2,p4),pi.pow(2).div(3));derivative.D=S(derivative.e2,derivative.e1);
  }
  return{z,wp,pp,e1,e2,e3,D:S(e2,e1),factor:chart.factor,matrixC:chart.matrixC,derivative};
 }
 function fromBasis(b,n){
  if(n===1)return b.z;
  if(n===2)return S(S(P(b.z,3),K(M(b.wp,b.z),3)),b.pp);
  const T=V(S(b.e3,b.e1),b.D),Q={'-2':one,'-1':one,0:b.z};
  let G=V(b.pp,b.D),R=A(M(V(S(b.wp,b.e1),b.D),b.z),K(G,'.5'));
  for(let j=1;j<n;j++){
   const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1],phi=A(M(R,q2),K(G,d(j-1).div(2))),Y=M(M(q3,q2),q1),u=S(phi,Y),v=S(phi,M(T,Y));
   const H=M(M(phi,u),v),Hp=A(A(M(phi,u),M(phi,v)),M(u,v));
   const nextG=V(H,M(G,P(q3,3))),nextQ=S(V(M(q3,nextG),G),K(V(Hp,M(G,P(q3,2))),d(j).div(2)));
   Q[j]=nextQ;
   if(j<n-1)R=V(A(M(phi,nextQ),K(M(q3,nextG),d(j+1).div(2))),M(q3,q1));
   G=nextG;
  }
  const j=n-1,k=j%2?(j+1)**2/4:j*(j+2)/4,lead=K(P(b.D,-k),d(2).pow(-2*k));
  return V(Q[n-1],lead);
 }
 function value(tau,r,s,n){const b=basis(tau,r,s,n);return n>12?fromBasis(b,n):V(fromBasis(b,n),P(b.factor,n*(n+1)/2));}
 function valueAndDerivative(tau,r,s,n){
  if(n===1)return baseValueAndDerivative(tau,r,s);
  const b=basis(tau,r,s,n,true),zero=C();
  // First-order complex arithmetic propagates the exact series derivative
  // through the same premodule recurrence used for the equation itself.
  const constant=z=>[z,zero],jet=k=>[b[k],b.derivative[k]],
   ja=(a,b)=>[A(a[0],b[0]),A(a[1],b[1])],js=(a,b)=>[S(a[0],b[0]),S(a[1],b[1])],
   jm=(a,b)=>[M(a[0],b[0]),A(M(a[1],b[0]),M(a[0],b[1]))],
   jk=(a,k)=>[K(a[0],k),K(a[1],k)],
   jv=(a,b)=>{const q=V(a[0],b[0]);return[q,V(S(a[1],M(q,b[1])),b[0])];},
   jp=(a,n)=>[P(a[0],n),K(M(P(a[0],n-1),a[1]),n)];
  const z=jet('z'),wp=jet('wp'),pp=jet('pp');let f;
  if(n===2)f=js(js(jp(z,3),jk(jm(wp,z),3)),pp);
  else{
   const e1=jet('e1'),e3=jet('e3'),delta=jet('D'),T=jv(js(e3,e1),delta),Q={'-2':constant(one),'-1':constant(one),0:z};
   let G=jv(pp,delta),R=ja(jm(jv(js(wp,e1),delta),z),jk(G,'.5'));
   for(let j=1;j<n;j++){
    const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1],phi=ja(jm(R,q2),jk(G,d(j-1).div(2))),Y=jm(jm(q3,q2),q1),u=js(phi,Y),v=js(phi,jm(T,Y));
    const H=jm(jm(phi,u),v),Hp=ja(ja(jm(phi,u),jm(phi,v)),jm(u,v));
    const nextG=jv(H,jm(G,jp(q3,3))),nextQ=js(jv(jm(q3,nextG),G),jk(jv(Hp,jm(G,jp(q3,2))),d(j).div(2)));
    Q[j]=nextQ;
    if(j<n-1)R=jv(ja(jm(phi,nextQ),jk(jm(q3,nextG),d(j+1).div(2))),jm(q3,q1));
    G=nextG;
   }
   const j=n-1,k=j%2?(j+1)**2/4:j*(j+2)/4,lead=jk(jp(delta,-k),d(2).pow(-2*k));
   f=jv(Q[n-1],lead);
  }
  const weight=n*(n+1)/2;
  if(n>12)return{value:f[0],derivative:V(f[1],P(b.factor,2))};
  return{value:V(f[0],P(b.factor,weight)),
   derivative:S(V(f[1],P(b.factor,weight+2)),K(V(f[0],P(b.factor,weight+1)),d(weight).times(b.matrixC)))};
 }
 function baseValueAndDerivative(tau,r,s){
  const b=basis(tau,r,s,1),sum=A(b.pp,K(M(A(b.wp,b.eta1),b.z),2));
  const reducedDerivative=K([sum[1].neg(),sum[0]],d(1).div(pi.times(4)));
  const derivative=S(V(reducedDerivative,P(b.factor,3)),K(V(b.z,P(b.factor,2)),b.matrixC));
  return{value:V(b.z,b.factor),derivative};
 }
 function correct(seed,r,s,n,{maxIterations=18,trustRadius=.3,tolerance='1e-26'}={}){
  let tau=asComplex(seed),atTau=null;const start=tau,tol=d(tolerance);
  for(let iteration=0;iteration<maxIterations;iteration++){
   let f,df;
   {const exact=atTau??valueAndDerivative(tau,r,s,n);f=exact.value;df=exact.derivative;atTau=null;}
   const step=V(f,df),error=abs(step);
   if(!error.isFinite())return{status:'precision_limited',iterations:iteration};
   if(error.lt(tol))return{status:'ok',tau,iterations:iteration,correction:error};
   let accepted=false;
   for(let back=0;back<18;back++){
    const trial=S(tau,K(step,d(2).pow(-back)));
    if(!inChart(trial)||abs(S(trial,start)).gt(trustRadius))continue;
    // Overlap already retained this derivative; the other pages retain their original value-only trial.
    const trialEvaluation=reuseAcceptedEvaluation&&back===0?valueAndDerivative(trial,r,s,n):null;
    if(abs(trialEvaluation?.value??value(trial,r,s,n)).lt(abs(f))){tau=trial;atTau=trialEvaluation;accepted=true;break;}
   }
   if(!accepted)return{status:'correction_failed',iterations:iteration,tau};
  }
  return{status:'iteration_limit',tau};
 }
 const engine={D,d,C,A,S,M,V,K,P,abs,inChart,basis,value,valueAndDerivative,baseValueAndDerivative,correct,precision,asComplex};
 precisionEngineCache.set(precision,engine);return engine;
}
return createPrecisionEngine;
}

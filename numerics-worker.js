/* Hemisphere Atlas: self-contained classic numerical worker.
 * Bundled without changing the numerical equations, continuation, or branch
 * labels in the checked atlas engine. Geometric anchors select the sheet.
 * Supported order: n = 1 + a + b + c <= 12. Interior samples only.
 * Interactive point requests stay in the compact source triangle with all
 * barycentric weights >= 0.04; the monodromy boundary is never sampled.
 * Results are numerical approximations, not interval-certified values.
 */
"use strict";

const atlasMath = (() => {
/* Elliptic monodromy atlas. No external dependencies.
 * Numerical roots, not interval certificates. Parameter continuation preserves
 * the chosen F0 chart; modular transformations are used ONLY inside evaluation.
 */

const PI = Math.PI;
const TAU = 2 * PI;
const SOURCE_MARGIN = 0.004;
const ROOT_TOL = 1e-13;
const MAX_Y = 60;
const c = (re = 0, im = 0) => ({ re, im });
const add = (a, b) => c(a.re + b.re, a.im + b.im);
const sub = (a, b) => c(a.re - b.re, a.im - b.im);
const mul = (a, b) => c(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const scale = (a, k) => c(a.re * k, a.im * k);
const abs = a => Math.hypot(a.re, a.im);
const div = (a, b) => {
  const d = b.re * b.re + b.im * b.im;
  return c((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
};
const finite = z => Number.isFinite(z.re) && Number.isFinite(z.im);
const polar = (m, a) => c(m * Math.cos(a), m * Math.sin(a));
const asComplex = z => Array.isArray(z) ? c(z[0], z[1]) : c(z.re, z.im);
const centered = x => x - Math.round(x);

/** Standard SL2 reduction for evaluation, with its covariant character action.
 * If t' = gamma t, then r'=a r-b s, s'=-c r+d s, and Z'=factor Z.
 * This function NEVER changes a solver's root or displayed marking.
 */
function evaluationChart(tau, r, s) {
  let t = asComplex(tau), rr = r, ss = s, factor = c(1, 0);
  if (!finite(t) || !(t.im > 0)) throw new RangeError('The modulus must lie in the upper half-plane.');
  for (let k = 0; k < 48; k++) {
    const integer = Math.round(t.re);
    if (integer) {
      t.re -= integer;
      rr += integer * ss;
    }
    rr = centered(rr);
    ss = centered(ss);
    if (abs(t) >= 1 - 2e-15) return { tau: t, r: rr, s: ss, factor, transformations: k };
    factor = mul(t, factor);
    t = div(c(-1, 0), t);
    [rr, ss] = [ss, -rr];
  }
  throw new RangeError('Elliptic evaluation could not find a stable chart.');
}

/** Fourier evaluation after reduction. Each exponential has nonpositive real
 * part, including near cusps, so no large theta terms are formed and cancelled.
 */
function reducedElliptic(t, r, s) {
  const ar = r + s * t.re, ai = s * t.im;
  const sign = ai >= 0 ? 1 : -1;
  const w = polar(Math.exp(-TAU * Math.abs(ai)), sign * TAU * ar);
  const denominator = sub(w, c(1, 0));
  if (abs(denominator) < 2e-15) throw new RangeError('The elliptic argument is too close to a lattice pole.');
  const quotient = div(add(w, c(1, 0)), denominator);
  const cot = c(-sign * quotient.im, sign * quotient.re);
  const csc2 = scale(div(w, mul(denominator, denominator)), -4);
  let sinSum = c(), pSum = c(), ppSum = c(), eisen2 = c(), eisen4 = c(), eisen6 = c();
  let last = Infinity, terms = 0;
  for (let k = 1; k <= 32; k++) {
    const qk = polar(Math.exp(-TAU * k * t.im), TAU * k * t.re);
    const den = sub(c(1, 0), qk);
    const coeff = div(qk, den);
    const plus = div(polar(Math.exp(-TAU * k * t.im * (1 + s)),
      TAU * k * (r + (1 + s) * t.re)), den);
    const minus = div(polar(Math.exp(-TAU * k * t.im * (1 - s)),
      TAU * k * (-r + (1 - s) * t.re)), den);
    const difference = sub(plus, minus);
    const sine = c(difference.im / 2, -difference.re / 2);
    const cosine = scale(add(plus, minus), .5);
    sinSum = add(sinSum, sine);
    pSum = add(pSum, scale(sub(coeff, cosine), k));
    ppSum = add(ppSum, scale(sine, k * k));
    eisen2 = add(eisen2, scale(coeff, k));
    eisen4 = add(eisen4, scale(coeff, k ** 3));
    eisen6 = add(eisen6, scale(coeff, k ** 5));
    last = Math.max(abs(plus), abs(minus), abs(coeff)) * k ** 5;
    terms = k;
    if (k >= 5 && last < 3e-18) break;
  }
  if (!(last < 1e-14)) throw new RangeError('The elliptic Fourier series did not converge.');
  const z = add(add(scale(cot, PI), scale(sinSum, 4 * PI)), c(0, TAU * s));
  const wp = add(scale(sub(csc2, c(1 / 3, 0)), PI * PI), scale(pSum, 8 * PI * PI));
  const wpPrime = add(scale(mul(cot, csc2), -2 * PI ** 3), scale(ppSum, 16 * PI ** 3));
  const result = { z, wp, wpPrime, terms };
  const argument = c(ar, ai);
  if (abs(argument) < .22) {
    // Remove the universal Laurent pole BEFORE forming Z^3-3 wp Z-wp'.
    // Direct polynomial evaluation loses digits when both character entries
    // approach integers, despite a misleadingly small relative residual.
    const eta1 = scale(sub(c(1, 0), scale(eisen2, 24)), PI * PI / 3);
    const g2 = scale(add(c(1, 0), scale(eisen4, 240)), 4 * PI ** 4 / 3);
    const g3 = scale(sub(c(1, 0), scale(eisen6, 504)), 8 * PI ** 6 / 27);
    const coefficients = [c(), scale(g2, 1 / 20), scale(g3, 1 / 28)];
    for (let m = 3; m <= 14; m++) {
      let sum = c();
      for (let i = 1; i <= m - 2; i++) sum = add(sum, mul(coefficients[i], coefficients[m - 1 - i]));
      coefficients.push(scale(sum, 3 / ((m - 2) * (2 * m + 3))));
    }
    let H = sub(c(0, TAU * s), mul(eta1, argument)), P = c(), Pp = c();
    const a2 = mul(argument, argument);
    let evenPower = a2, oddPower = argument;
    for (let m = 1; m <= 14; m++) {
      P = add(P, mul(coefficients[m], evenPower));
      Pp = add(Pp, scale(mul(coefficients[m], oddPower), 2 * m));
      H = sub(H, scale(mul(coefficients[m], mul(evenPower, argument)), 1 / (2 * m + 1)));
      evenPower = mul(evenPower, a2);
      oddPower = mul(oddPower, a2);
    }
    const H2 = mul(H, H);
    const stableTerms = [scale(div(H2, argument), 3), mul(H2, H),
      scale(div(P, argument), -3), scale(mul(P, H), -3), scale(Pp, -1)];
    result.stableZ2 = stableTerms.reduce(add, c());
    result.stableNormalization = Math.max(1, stableTerms.reduce((sum, term) => sum + abs(term), 0));
  }
  return result;
}

/** Return Z or Z^(2), together with elliptic values and a cancellation-aware
 * relative residual. tau can be {re,im} or [re,im].
 */
function evaluate(tau, r, s, n = 1) {
  if (n !== 1 && n !== 2) throw new RangeError('This engine implements indices one and two.');
  const chart = evaluationChart(tau, r, s);
  const e = reducedElliptic(chart.tau, chart.r, chart.s);
  const z3 = mul(mul(e.z, e.z), e.z), pz = scale(mul(e.wp, e.z), 3);
  const fReduced = n === 1 ? e.z : (e.stableZ2 ?? sub(sub(z3, pz), e.wpPrime));
  const normalization = n === 1 ? Math.max(1, abs(e.z), TAU) :
    (e.stableNormalization ?? Math.max(1, abs(z3) + abs(pz) + abs(e.wpPrime)));
  const f2 = mul(chart.factor, chart.factor), f3 = mul(f2, chart.factor);
  const value = div(fReduced, n === 1 ? chart.factor : f3);
  const result = { value, z: div(e.z, chart.factor), wp: div(e.wp, f2),
    wpPrime: div(e.wpPrime, f3), residual: abs(value),
    relativeResidual: abs(fReduced) / normalization,
    terms: e.terms, transformations: chart.transformations };
  if (!finite(value)) throw new RangeError('Elliptic evaluation lost numerical precision.');
  return result;
}

const thetaZ = (tau, r, s) => evaluate(tau, r, s, 1).value;
const Z2 = (tau, r, s) => evaluate(tau, r, s, 2).value;
const cauchyZ = thetaZ;

function inF0(tau, strict = true) {
  const t = asComplex(tau), epsilon = strict ? 2e-9 : -2e-12;
  return finite(t) && t.re > epsilon && t.re < 1 - epsilon &&
    t.im > 2e-6 && t.im < MAX_Y &&
    (t.re - .5) ** 2 + t.im ** 2 > .25 + epsilon;
}

function sourceStatus(r, s, margin = SOURCE_MARGIN) {
  if (![r, s].every(Number.isFinite)) return 'invalid_parameters';
  if (!(r > 0 && s > 0 && r < .5 && s < .5 && r + s > .5)) return 'outside_source_triangle';
  if (r < margin || s < margin || .5 - r < margin || .5 - s < margin || r + s - .5 < margin)
    return 'near_boundary';
  return 'ok';
}


function buildSourceMesh(resolution = 14, inset = .025) {
  const n = Math.max(2, Math.min(80, Math.round(resolution)));
  if (!(inset >= .009 && inset < 1 / 3)) throw new RangeError('Use a barycentric inset between .009 and 1/3.');
  const points = [], triangles = [], lookup = new Map();
  for (let i = 0; i <= n; i++) for (let j = 0; j <= n - i; j++) {
    const a = inset + (1 - 3 * inset) * i / n;
    const b = inset + (1 - 3 * inset) * j / n;
    const d = 1 - a - b;
    lookup.set(`${i},${j}`, points.length);
    points.push({ r: .5 * (a + b), s: .5 * (b + d), barycentric: [a, b, d] });
  }
  const at = (i, j) => lookup.get(`${i},${j}`);
  for (let i = 0; i < n; i++) for (let j = 0; j < n - i; j++) {
    triangles.push([at(i, j), at(i + 1, j), at(i, j + 1)]);
    if (i + j < n - 1) triangles.push([at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
  }
  return { points, triangles, resolution: n, inset };
}


return {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh};
})();

const atlasAttachment = (() => {
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

const MAX_N = 12;
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
 const counts=options.counts??[0,0,0],checked=validateCounts(counts);
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
  if(options.signal?.aborted)return{...mesh,counts:[...counts],n:checked.n,samples,status:'aborted',completed};
  const p=mesh.points[index],sample=solveAttachment(p.r,p.s,counts,previous);
  samples[index]=sample;completed++;
  if(sample.status==='ok')previous=sample;
  if(options.onSample)await options.onSample(sample,index,completed,mesh.points.length);
  // Every point yields to queued interactive jobs in a worker.
  await new Promise(resolve=>setTimeout(resolve,0));
 }
 return{...mesh,counts:[...counts],n:checked.n,samples,completed,
  status:samples.every(s=>s?.status==='ok')?'ok':'partial',
  metadata:{sampledInteriorOnly:true,sourceCoordinateSideGap:mesh.inset/2,
   maximumRelativeResidual:Math.max(...samples.flatMap(s=>s.branches.filter(b=>b.status==='ok').map(b=>b.relativeResidual)))}};
}

return {MAX_N,validateCounts,setAnchorData,solveAttachment,sampleAttachment};
})();

const attachmentAnchors = {"r":0.31,"s":0.29,"maxN":12,"anchors":{"0,0,0":{"tau":[0.5904953761444915,0.7963943351098978]},"0,0,1":{"tau":[0.7515095138734595,0.6524433021803057]},"0,1,0":{"tau":[0.35448685090813264,0.6879997451009816]},"1,0,0":{"tau":[0.5823785379137855,1.1341826530853316]},"0,0,2":{"tau":[0.8103860604950234,0.5797764171757713]},"0,1,1":{"tau":[0.5400832084805506,0.646138229333322]},"0,2,0":{"tau":[0.2661530927667447,0.6170110529155463]},"1,0,1":{"tau":[0.7554851567438122,0.9044222083961427]},"1,1,0":{"tau":[0.33641698617555044,0.9361692214053354]},"2,0,0":{"tau":[0.5802309037096132,1.3519173892371166]},"0,0,3":{"tau":[0.841728094023958,0.5344069082972436]},"0,1,2":{"tau":[0.6299910352020642,0.6065667355063957]},"0,2,1":{"tau":[0.4336528965362644,0.6170485958998836]},"0,3,0":{"tau":[0.2194990271876052,0.5699249827074794]},"1,0,2":{"tau":[0.8190763243004178,0.788631148844644]},"1,1,1":{"tau":[0.5345911839088004,0.8413672802455601]},"1,2,0":{"tau":[0.24508363041652606,0.820541713936849]},"2,0,1":{"tau":[0.7566877827649808,1.0780300107149037]},"2,1,0":{"tau":[0.330728071930688,1.1089889089128258]},"3,0,0":{"tau":[0.5793596671311855,1.5133680043418907]},"0,0,4":{"tau":[0.8615254401911367,0.5026248506350188]},"0,1,3":{"tau":[0.6842316944717046,0.575021459683328]},"0,2,2":{"tau":[0.5267688263467021,0.6000381495242691]},"0,3,1":{"tau":[0.36881306450740986,0.5894886050882322]},"0,4,0":{"tau":[0.19030803642972227,0.5361316608784794]},"1,0,3":{"tau":[0.8526909974897778,0.7164178524340136]},"1,1,2":{"tau":[0.6314478387366457,0.7723654650396221]},"1,2,1":{"tau":[0.42284498422152705,0.7812601291251633]},"1,3,0":{"tau":[0.1973984260185703,0.7458993031656952]},"2,0,2":{"tau":[0.8220118840367037,0.9366554849150815]},"2,1,1":{"tau":[0.5326043104373224,0.9845389328712051]},"2,2,0":{"tau":[0.23787172685968883,0.9673218781297337]},"3,0,1":{"tau":[0.7572248859548885,1.2123248746321915]},"3,1,0":{"tau":[0.32814204362632093,1.2429897592333141]},"4,0,0":{"tau":[0.5789215292300838,1.64185870386297]},"0,0,5":{"tau":[0.87532736448948,0.47873158903313423]},"0,1,4":{"tau":[0.7209441947002367,0.5497539960286885]},"0,2,3":{"tau":[0.587659408841123,0.5813159221898846]},"0,3,2":{"tau":[0.45853298654130986,0.5864119195544156]},"0,4,1":{"tau":[0.32481896292121665,0.5658832421781945]},"0,5,0":{"tau":[0.17012729074846442,0.5104168701229547]},"1,0,4":{"tau":[0.8737146117761868,0.6659979038828044]},"1,1,3":{"tau":[0.6897768161485427,0.7213488596924934]},"1,2,2":{"tau":[0.5228133583527507,0.7416723275535352]},"1,3,1":{"tau":[0.3550011865750481,0.7333908261815947]},"1,4,0":{"tau":[0.16791151078300529,0.6931050520647739]},"2,0,3":{"tau":[0.856600950061261,0.8473150080377062]},"2,1,2":{"tau":[0.632157377106861,0.8969694583022713]},"2,2,1":{"tau":[0.41854260474439214,0.905436867652474]},"2,3,0":{"tau":[0.18945832036306007,0.875279890371322]},"3,0,2":{"tau":[0.8234166474721952,1.0536888952271188]},"3,1,1":{"tau":[0.5316145465016989,1.0991956202308]},"3,2,0":{"tau":[0.23439333176450308,1.0838565533760522]},"4,0,1":{"tau":[0.7575144695144952,1.3223411987468239]},"4,1,0":{"tau":[0.3267340104740559,1.3528628092637536]},"5,0,0":{"tau":[0.5786704237611745,1.7486344886564826]},"0,0,6":{"tau":[0.8855883749676,0.459893337978774]},"0,1,5":{"tau":[0.7476547390629024,0.5290741282185445]},"0,2,4":{"tau":[0.6310714701953587,0.5638367265401152]},"0,3,3":{"tau":[0.5204211943977096,0.5770842731764217]},"0,4,2":{"tau":[0.40965039211323434,0.5716922698326926]},"0,5,1":{"tau":[0.29281381542373114,0.5458775986357257]},"0,6,0":{"tau":[0.1552308493390186,0.49000353533730145]},"1,0,5":{"tau":[0.8882196634076015,0.6282589367010324]},"1,1,4":{"tau":[0.7290769007922019,0.6820042681833915]},"1,2,3":{"tau":[0.5881532493915438,0.7075651684523052]},"1,3,2":{"tau":[0.4510576512137605,0.7117896923823168]},"1,4,1":{"tau":[0.3092118506224589,0.695216006203829]},"1,5,0":{"tau":[0.1477519141514894,0.6533558004801285]},"2,0,4":{"tau":[0.8782133522035309,0.7844578175952746]},"2,1,3":{"tau":[0.6922991189465529,0.8329147799736228]},"2,2,2":{"tau":[0.5212230759002172,0.8511681141462283]},"2,3,1":{"tau":[0.34917986947047064,0.8442598735593233]},"2,4,0":{"tau":[0.15960719886127595,0.8098782385880396]},"3,0,3":{"tau":[0.8585464147633195,0.9521664394926173]},"3,1,2":{"tau":[0.6325417195306986,0.9985719839741891]},"3,2,1":{"tau":[0.4162831369003377,1.0068608466821334]},"3,3,0":{"tau":[0.18548711432422435,0.979494178808416]},"4,0,2":{"tau":[0.8242124660453741,1.151180299071085]},"4,1,1":{"tau":[0.5310399425071837,1.1953234039601657]},"4,2,0":{"tau":[0.2324122394543271,1.1810929728116042]},"5,0,1":{"tau":[0.7576893917453332,1.415707662218239]},"5,1,0":{"tau":[0.32587834922456627,1.4461479825119008]},"6,0,0":{"tau":[0.578513211933155,1.8400026346529397]},"0,0,7":{"tau":[0.8935695170631885,0.4445238550546686]},"0,1,6":{"tau":[0.7680816619446168,0.5117859055568753]},"0,2,5":{"tau":[0.6638087929376337,0.5481312045758051]},"0,3,4":{"tau":[0.5662502583558813,0.5660674216966612]},"0,4,3":{"tau":[0.47038431580611356,0.5691191152884885]},"0,5,2":{"tau":[0.37270851191302595,0.5575773347273169]},"0,6,1":{"tau":[0.2683647278635283,0.52879748719793]},"0,7,0":{"tau":[0.14371472207625965,0.4732798229549187]},"1,0,6":{"tau":[0.8988946785558809,0.5986474112765273]},"1,1,5":{"tau":[0.757515289092806,0.6505645672211018]},"1,2,4":{"tau":[0.6345895968123779,0.6786954967391343]},"1,3,3":{"tau":[0.5173362701255009,0.6896700070199766]},"1,4,2":{"tau":[0.39982624715111464,0.6851433733448178]},"1,5,1":{"tau":[0.276090469234877,0.6641428959307979]},"1,6,0":{"tau":[0.13302346902452356,0.6220753749174971]},"2,0,5":{"tau":[0.8930918423482282,0.7371942141027298]},"2,1,4":{"tau":[0.7328497447335313,0.7837854959610272]},"2,2,3":{"tau":[0.5885851008708306,0.8064731758072106]},"2,3,2":{"tau":[0.447737001043542,0.8104257477661121]},"2,4,1":{"tau":[0.3023745099048237,0.7961194599484742]},"2,5,0":{"tau":[0.13926952084402813,0.7604944279863987]},"3,0,4":{"tau":[0.8805111459535212,0.88014891120271]},"3,1,3":{"tau":[0.6936929117459214,0.9249420837481607]},"3,2,2":{"tau":[0.5203604410850786,0.9420472610477201]},"3,3,1":{"tau":[0.3460120862317194,0.9359838200637809]},"3,4,0":{"tau":[0.1553506593979241,0.9048369846987677]},"4,0,3":{"tau":[0.8596823845481332,1.040458892557427]},"4,1,2":{"tau":[0.6327740811462541,1.0849285948795253]},"4,2,1":{"tau":[0.41492179738690005,1.0931249845327629]},"4,3,0":{"tau":[0.18315877131009872,1.0674498166149522]},"5,0,2":{"tau":[0.8247117287312125,1.2350094503071851]},"5,1,1":{"tau":[0.5306736884551937,1.2782904040218237]},"5,2,0":{"tau":[0.23116478113594885,1.264772313747665]},"6,0,1":{"tau":[0.7578034945743788,1.4968886398046868]},"6,1,0":{"tau":[0.3253179805035468,1.5272780873967173]},"7,0,0":{"tau":[0.5784082889971616,1.9198626778550154]},"0,0,8":{"tau":[0.8999886735985638,0.43165749057710767]},"0,1,7":{"tau":[0.7842843731758058,0.4970654360110489]},"0,2,6":{"tau":[0.6895018910795871,0.5341374020413163]},"0,3,5":{"tau":[0.6017919435613742,0.5549293341084078]},"0,4,4":{"tau":[0.5166510006237209,0.5631702795678296]},"0,5,3":{"tau":[0.43147276876014107,0.5599209902678184]},"0,6,2":{"tau":[0.34369136170007636,0.5445324732829234]},"0,7,1":{"tau":[0.24900061019025993,0.5140511474271879]},"0,8,0":{"tau":[0.1345010488513675,0.4592428537479833]},"1,0,7":{"tau":[0.9071185345899363,0.5746069889852019]},"1,1,6":{"tau":[0.7791416217273306,0.6247271194210997]},"1,2,5":{"tau":[0.6694686003680538,0.6541175224491951]},"1,3,4":{"tau":[0.566291131542106,0.6690524526333079]},"1,4,3":{"tau":[0.464760493774037,0.6715543759066558]},"1,5,2":{"tau":[0.36126522625152746,0.6618033415718381]},"1,6,1":{"tau":[0.25093182248719836,0.6383140187936306]},"1,7,0":{"tau":[0.1217443668562202,0.5966418365607582]},"2,0,6":{"tau":[0.9040110611158669,0.7000080185294845]},"2,1,5":{"tau":[0.762178227429426,0.7446696257600844]},"2,2,4":{"tau":[0.6365164099219819,0.7694924688840932]},"2,3,3":{"tau":[0.5160314298136841,0.7793499753523696]},"2,4,2":{"tau":[0.3952342014435136,0.7754820407072519]},"2,5,1":{"tau":[0.2685551858462416,0.7572438938373648]},"2,6,0":{"tau":[0.12446556844497446,0.7215668307205353]},"3,0,5":{"tau":[0.8956278040791181,0.8256929702307054]},"3,1,4":{"tau":[0.734983197185718,0.8684124065360601]},"3,2,3":{"tau":[0.5888624815053458,0.8894315983125628]},"3,3,2":{"tau":[0.4458530141184615,0.8932607935389507]},"3,4,1":{"tau":[0.2985581051723347,0.8803522880710194]},"3,5,0":{"tau":[0.13484413406375953,0.8481938901081237]},"4,0,4":{"tau":[0.8818814555448774,0.9613304545714055]},"4,1,3":{"tau":[0.694557975951264,1.0038902962802296]},"4,2,2":{"tau":[0.5198248176004105,1.0202785499707214]},"4,3,1":{"tau":[0.3440525661399386,1.0147704677679885]},"4,4,0":{"tau":[0.15280393431861738,0.9856218132307234]},"5,0,3":{"tau":[0.8604126867429136,1.1170602611367253]},"5,1,2":{"tau":[0.632926040787829,1.1602697414617]},"5,2,1":{"tau":[0.41402886114936516,1.1684116115104175]},"5,3,0":{"tau":[0.18165721310513627,1.1438482089549935]},"6,0,2":{"tau":[0.8250473662854825,1.308670545276684]},"6,1,1":{"tau":[0.5304247761233722,1.351368916988763]},"6,2,0":{"tau":[0.23032393881576027,1.338337097595142]},"7,0,1":{"tau":[0.7578821950635308,1.568741781822955]},"7,1,0":{"tau":[0.3249304151620226,1.599097138180921]},"8,0,0":{"tau":[0.5783347884054936,1.9907964524404396]},"0,0,9":{"tau":[0.9052863407274091,0.42066837982167543]},"0,1,8":{"tau":[0.7975003961237308,0.4843360211540446]},"0,2,7":{"tau":[0.7102810744526495,0.5216545187363899]},"0,3,6":{"tau":[0.6302881075104035,0.5442438663912446]},"0,4,5":{"tau":[0.5533301666697178,0.5558690621223308]},"0,5,4":{"tau":[0.47722162225118536,0.5579177435088031]},"0,6,3":{"tau":[0.4002268745047575,0.5505193291289712]},"0,7,2":{"tau":[0.32022070355933313,0.5326224962894985]},"0,8,1":{"tau":[0.2332316901885675,0.5011748764210534]},"0,9,0":{"tau":[0.12693194518398745,0.44723328196147555]},"1,0,8":{"tau":[0.9136737045620708,0.5545792734711206]},"1,1,7":{"tau":[0.7962008484973406,0.6030158039280128]},"1,2,6":{"tau":[0.6967284862987426,0.6329664899109504]},"1,3,5":{"tau":[0.6041306363884229,0.6503436645774036]},"1,4,4":{"tau":[0.5141192369559593,0.6572773685404725]},"1,5,3":{"tau":[0.42400586746770635,0.6544141023426473]},"1,6,2":{"tau":[0.331101689698656,0.6413589740543262]},"1,7,1":{"tau":[0.23111429379408077,0.6164446751774848]},"1,8,0":{"tau":[0.1127983490862128,0.5754377482669324]},"2,0,7":{"tau":[0.912397230479871,0.6697708427064374]},"2,1,6":{"tau":[0.7844552780175824,0.712614200135365]},"2,2,5":{"tau":[0.6725187345711167,0.7384467602617243]},"2,3,4":{"tau":[0.5665518369109653,0.7517816841966645]},"2,4,3":{"tau":[0.462099501317581,0.7540950926399174]},"2,5,2":{"tau":[0.35572927723620884,0.7455407913806367]},"2,6,1":{"tau":[0.2429065728888235,0.7251063477057245]},"2,7,0":{"tau":[0.11317042601141826,0.6898883177182566]},"3,0,6":{"tau":[0.9067124091344732,0.7826787101347856]},"3,1,5":{"tau":[0.7648628516050452,0.8233650732281796]},"3,2,4":{"tau":[0.6376872271367137,0.8462169074874232]},"3,3,3":{"tau":[0.5152900929628823,0.8554069855901278]},"3,4,2":{"tau":[0.39255422302925636,0.8519925774317841]},"3,5,1":{"tau":[0.2642680884070581,0.8354776252924652]},"3,6,0":{"tau":[0.11994084773381948,0.8033902657430243]},"4,0,5":{"tau":[0.8971639809366251,0.9011761716805341]},"4,1,4":{"tau":[0.7363335146569225,0.9414993163681471]},"4,2,3":{"tau":[0.589046567364497,0.9614434107222732]},"4,3,2":{"tau":[0.4446503039921923,0.9652054073582298]},"4,4,1":{"tau":[0.2961510986030259,0.9532245778049533]},"4,5,0":{"tau":[0.1321565057912849,0.9232359840162704]},"5,0,4":{"tau":[0.8827780402701756,1.0322177292958412]},"5,1,3":{"tau":[0.6951374123938552,1.073296201030138]},"5,2,2":{"tau":[0.5194643875961141,1.0892017355913224]},"5,3,1":{"tau":[0.3427400369916056,1.0840791608266542]},"5,4,0":{"tau":[0.1511330994806695,1.0562653963191864]},"6,0,3":{"tau":[0.8609137234808509,1.1848760155011235]},"6,1,2":{"tau":[0.6330312822893928,1.2272136268889113]},"6,2,1":{"tau":[0.4134077837747312,1.2353204330404501]},"6,3,0":{"tau":[0.18062458321973887,1.2115308098835826]},"7,0,2":{"tau":[0.8252846473115993,1.3744349481415081]},"7,1,1":{"tau":[0.5302474343356233,1.4167199794724403]},"7,2,0":{"tau":[0.22972834561982114,1.404035548910081]},"8,0,1":{"tau":[0.7579388342521736,1.6332140920183909]},"8,1,0":{"tau":[0.3246509413617829,1.6635454197877955]},"9,0,0":{"tau":[0.5782813053022888,2.054602837847386]},"0,0,10":{"tau":[0.9097487527704582,0.4111306827603229]},"0,1,9":{"tau":[0.8085205030604694,0.4731841000501837]},"0,2,8":{"tau":[0.727484913591236,0.5104695998726397]},"0,3,7":{"tau":[0.6537225575156654,0.5342070175203936]},"0,4,6":{"tau":[0.5832557827619747,0.5481280598591305]},"0,5,5":{"tau":[0.5141321200052162,0.5537657245792036]},"0,6,4":{"tau":[0.4449961846412119,0.551605558803055]},"0,7,3":{"tau":[0.3745100535182204,0.5413557437284982]},"0,8,2":{"tau":[0.3007934566085034,0.5217792175341178]},"0,9,1":{"tau":[0.22010452946094078,0.4898137606886881]},"0,10,0":{"tau":[0.12058200288812793,0.4367980258823749]},"1,0,9":{"tau":[0.9190384211101494,0.5375538951023546]},"1,1,8":{"tau":[0.8100412482120526,0.5844412215281697]},"1,2,7":{"tau":[0.7186835859387825,0.6145586210627302]},"1,3,6":{"tau":[0.6343628118808657,0.6334855626949449]},"1,4,5":{"tau":[0.5531333364818894,0.6433577571540114]},"1,5,4":{"tau":[0.4727572835445658,0.645028601098364]},"1,6,3":{"tau":[0.3913913054166512,0.6385772047483168]},"1,7,2":{"tau":[0.3068027654273528,0.6233524874653934]},"1,8,1":{"tau":[0.21505973900580513,0.5976350819356401]},"1,9,0":{"tau":[0.1055076289568248,0.557406650482888]},"2,0,8":{"tau":[0.9190606076711629,0.644560676162795]},"2,1,7":{"tau":[0.8020005852137946,0.6857419783728032]},"2,2,6":{"tau":[0.7006387670582919,0.7119871111940945]},"2,3,5":{"tau":[0.6056123636276993,0.7274684775197592]},"2,4,4":{"tau":[0.5130170313391375,0.7337170748264283]},"2,5,3":{"tau":[0.4202900175632987,0.7312134305812743]},"2,6,2":{"tau":[0.3248540459501834,0.7196891904709856]},"2,7,1":{"tau":[0.22273975015311986,0.6980096115053703]},"2,8,0":{"tau":[0.10424399065456832,0.6634698710594106]},"3,0,7":{"tau":[0.9152155571535485,0.7476036605229398]},"3,1,6":{"tau":[0.7875582513635886,0.7864254368915393]},"3,2,5":{"tau":[0.674390219373184,0.8101019886409881]},"3,3,4":{"tau":[0.5667527982850816,0.8224381075576594]},"3,4,3":{"tau":[0.4605137349559376,0.8246607204324721]},"3,5,2":{"tau":[0.3524306517251559,0.8169135208837309]},"3,6,1":{"tau":[0.23827210887070233,0.7984044507900027]},"3,7,0":{"tau":[0.10858981240248783,0.7668390401184514]},"4,0,6":{"tau":[0.9083686362080391,0.853471128910256]},"4,1,5":{"tau":[0.7665877408192001,0.8916706712153726]},"4,2,4":{"tau":[0.6384606435127549,0.9132281048057832]},"4,3,3":{"tau":[0.5148113446942287,0.9219781522126247]},"4,4,2":{"tau":[0.3908091918996164,0.918884894011992]},"4,5,1":{"tau":[0.26152358969543366,0.9035272110690888]},"4,6,0":{"tau":[0.11716117880608924,0.8737082193732995]},"5,0,5":{"tau":[0.8981825574879845,0.9674058441093298]},"5,1,4":{"tau":[0.7372529540960431,1.006118644389196]},"5,2,3":{"tau":[0.5891746214659662,1.0253215576534365]},"5,3,2":{"tau":[0.4438246155489609,1.0290425985006992]},"5,4,1":{"tau":[0.29451315034278686,1.0177125624401901]},"5,5,0":{"tau":[0.13037027880383903,0.989190627497009]},"6,0,4":{"tau":[0.883402440379279,1.0953251469013674]},"6,1,3":{"tau":[0.6955470786137932,1.1353624094097583]},"6,2,2":{"tau":[0.5192082355101879,1.1509252351605825]},"6,3,1":{"tau":[0.34181112754341164,1.1460818676471523]},"6,4,0":{"tau":[0.14996694760821253,1.1192104404698613]},"7,0,3":{"tau":[0.8612740853110933,1.2458070785396733]},"7,1,2":{"tau":[0.6331074021030786,1.2875138975974518]},"7,2,1":{"tau":[0.4129566717215887,1.2955967017886751]},"7,3,0":{"tau":[0.17988051687742626,1.2723691650393238]},"8,0,2":{"tau":[0.8254589686702685,1.4338726835189464]},"8,1,1":{"tau":[0.5301163985559172,1.4758532559781143]},"8,2,0":{"tau":[0.22929014552917062,1.4634259603517479]},"9,0,1":{"tau":[0.7579809867273274,1.6916950098860004]},"9,1,0":{"tau":[0.32444264395471506,1.7220087308402454]},"10,0,0":{"tau":[0.578241177463598,2.112585796401073]},"0,0,11":{"tau":[0.9135705198512508,0.40274319084995835]},"0,1,10":{"tau":[0.8178748914036359,0.4633057764560679]},"0,2,9":{"tau":[0.7419995615345809,0.5003921913931301]},"0,3,8":{"tau":[0.6733854903201331,0.5248610171571945]},"0,4,7":{"tau":[0.6082149762714262,0.5403819371013867]},"0,5,6":{"tau":[0.5446816420698254,0.5485497483853404]},"0,6,5":{"tau":[0.4816333320402924,0.5500281308875302]},"0,7,4":{"tau":[0.4180892754139823,0.544884048450245]},"0,8,3":{"tau":[0.3529239050607367,0.5326145495488299]},"0,9,2":{"tau":[0.2844107691180769,0.5118959676022632]},"0,10,1":{"tau":[0.2089795013447162,0.4796957634301651]},"0,11,0":{"tau":[0.11516327836048905,0.42761431079597856]},"1,0,10":{"tau":[0.923522083580442,0.5228436920813326]},"1,1,9":{"tau":[0.8215231642780433,0.5683143809096617]},"1,2,8":{"tau":[0.7367880463090198,0.59836942505751]},"1,3,7":{"tau":[0.6591386056231198,0.6182916746489376]},"1,4,6":{"tau":[0.5848640159391224,0.6301805299640165]},"1,5,5":{"tau":[0.5119820863111095,0.6349826256776173]},"1,6,4":{"tau":[0.43905263991071825,0.6330082494050312]},"1,7,3":{"tau":[0.36463837528057746,0.6240648804722684]},"1,8,2":{"tau":[0.28676838049023473,0.6073809424113679]},"1,9,1":{"tau":[0.2017602409055539,0.5812406655903369]},"1,10,0":{"tau":[0.09943620012791543,0.5418265137711081]},"2,0,9":{"tau":[0.9244965945055537,0.6231238419640506]},"2,1,8":{"tau":[0.816210312738608,0.6628008359309334]},"2,2,7":{"tau":[0.7232639602197353,0.6891250629829302]},"2,3,6":{"tau":[0.6368087234473495,0.7059664572287357]},"2,4,5":{"tau":[0.5532820844409287,0.7148409867471551]},"2,5,4":{"tau":[0.4705583079566005,0.7163735393510486]},"2,6,3":{"tau":[0.38684796399943494,0.7106298176971176]},"2,7,2":{"tau":[0.300010590519598,0.6971539737455481]},"2,8,1":{"tau":[0.2064339173398105,0.6747839892977446]},"2,9,0":{"tau":[0.09699457093004263,0.6410066821773549]},"3,0,8":{"tau":[0.9219625862618696,0.718300170447209]},"3,1,7":{"tau":[0.8054263691316562,0.7554462829596336]},"3,2,6":{"tau":[0.7030651717204842,0.7794188051572261]},"3,3,5":{"tau":[0.6065802114192296,0.793688038869172]},"3,4,4":{"tau":[0.5123726481422266,0.7995062578759802]},"3,5,3":{"tau":[0.4180178046003677,0.7972728654869734]},"3,6,2":{"tau":[0.3210707964935976,0.7867852240665781]},"3,7,1":{"tau":[0.21784528541401077,0.7671610731254896]},"3,8,0":{"tau":[0.09963559115192062,0.7363017420446022]},"4,0,7":{"tau":[0.916960169040881,0.8144531131952817]},"4,1,6":{"tau":[0.7895758107902655,0.8507391923521915]},"4,2,5":{"tau":[0.6756412810453973,0.8729810884061157]},"4,3,4":{"tau":[0.5668984611129122,0.884639177811809]},"4,4,3":{"tau":[0.45946160702678335,0.8868103076086719]},"4,5,2":{"tau":[0.35024985048842067,0.8796329618867217]},"4,6,1":{"tau":[0.23526991624878854,0.8624294725194341]},"4,7,0":{"tau":[0.10575014763779851,0.8331886558415479]},"5,0,6":{"tau":[0.909478336130215,0.9158150216727714]},"5,1,5":{"tau":[0.7677772064188196,0.952327367648329]},"5,2,4":{"tau":[0.6390028604706922,0.9729788345009572]},"5,3,3":{"tau":[0.5144782751688971,0.9814199027712702]},"5,4,2":{"tau":[0.38959268539149444,0.9785620709554678]},"5,5,1":{"tau":[0.25963277795822404,0.9640250154428804]},"5,6,0":{"tau":[0.11529505721141707,0.935752735570142]},"6,0,5":{"tau":[0.898900147925036,1.0266196244194499]},"6,1,4":{"tau":[0.7379121304054201,1.0641875888152226]},"6,2,3":{"tau":[0.5892673990558824,1.0828543465832914]},"6,3,2":{"tau":[0.44322851898129456,1.0865484268975043]},"6,4,1":{"tau":[0.29333838206637347,1.0756945523849286]},"6,5,0":{"tau":[0.12910937368416608,1.0482188756958317]},"7,0,4":{"tau":[0.8838573865685705,1.1523009263231543]},"7,1,3":{"tau":[0.6958486702907646,1.1915745597379706]},"7,2,2":{"tau":[0.5190187547853987,1.206884256897688]},"7,3,1":{"tau":[0.3411264129888459,1.2022500080874907]},"7,4,0":{"tau":[0.14911578112600785,1.1760716818208807]},"8,0,3":{"tau":[0.8615427773268887,1.3011769323276752]},"8,1,2":{"tau":[0.6331643609690645,1.3424113761332421]},"8,2,1":{"tau":[0.4126178490351031,1.350476962940132]},"8,3,0":{"tau":[0.17932493075464095,1.3276715574532627]},"9,0,2":{"tau":[0.8255909930696255,1.4881190103136215]},"9,1,1":{"tau":[0.5300167236508954,1.5298685448852112]},"9,2,0":{"tau":[0.22895789502527408,1.5176370820909544]},"10,0,1":{"tau":[0.7580132224444402,1.745212327495684]},"10,1,0":{"tau":[0.32428317208656704,1.7755127431036775]},"11,0,0":{"tau":[0.578210300862536,2.1657209292992188]}}};

atlasAttachment.setAnchorData(attachmentAnchors);

const POINT_INSET = 0.04;
const POINT_ROUNDOFF = 1e-12;
const pointStateCache = new Map();

function compactSourcePoint(r, s, inset = POINT_INSET) {
  if (![r, s].every(Number.isFinite)) {
    throw new RangeError("The source coordinates must be finite numbers.");
  }
  if (inset !== POINT_INSET) {
    throw new RangeError("Interactive points use the fixed compact inset 0.04.");
  }

  const weights = [1 - 2 * s, 2 * (r + s) - 1, 1 - 2 * r];

  if (weights.some(weight => weight < POINT_INSET - POINT_ROUNDOFF)) {
    throw new RangeError(
      "Choose a point in the compact source triangle: " +
      "all three barycentric weights must be at least 0.04."
    );
  }

  // Pointer projection can undershoot a compact edge by floating roundoff.
  // Restore just that tiny excess; do not move points outside the accepted K.
  if (weights.some(weight => weight < POINT_INSET)) {
    const clamped = weights.map(weight => Math.max(POINT_INSET, weight));
    const largest = clamped.indexOf(Math.max(...clamped));
    clamped[largest] -= clamped.reduce((sum, weight) => sum + weight, 0) - 1;
    return {r: (1 - clamped[2]) / 2, s: (1 - clamped[0]) / 2};
  }

  return {r, s};
}

self.onmessage = async ({data}) => {
  if (data.type === "solve") {
    const id = data.id;

    try {
      const checked = atlasAttachment.validateCounts(data.counts);
      if (checked.status !== "ok") throw new RangeError(checked.message);

      const {r, s} = compactSourcePoint(data.r, data.s, data.inset);
      const key = data.counts.join(",");
      const previous = pointStateCache.get(key) ?? null;
      const result = atlasAttachment.solveAttachment(r, s, data.counts, previous);

      // Only a complete numerical success may seed the next point. The
      // original fixed geometric anchors remain available for solver retries.
      if (result.status === "ok") {
        pointStateCache.delete(key);
        pointStateCache.set(key, result);

        if (pointStateCache.size > 8) {
          pointStateCache.delete(pointStateCache.keys().next().value);
        }
      }

      self.postMessage({type: "point", id, result});
    } catch (error) {
      self.postMessage({
        type: "pointError", id,
        message: error instanceof Error ? error.message : String(error)
      });
    }

    return;
  }

  if (data.type !== "compute") return;
  const id = data.id;
  try {
    const checked = atlasAttachment.validateCounts(data.counts);
    if (checked.status !== "ok") throw new RangeError(checked.message);
    const mesh = await atlasAttachment.sampleAttachment({
      counts: data.counts,
      resolution: data.resolution ?? 12,
      inset: data.inset ?? 0.04,
      onSample: (_sample, _index, completed, total) => {
        if (completed % 4 === 0 || completed === total) {
          self.postMessage({type: "progress", id, completed, total});
        }
      }
    });
    self.postMessage({type: "done", id, mesh});
  } catch (error) {
    self.postMessage({
      type: "error", id,
      message: error instanceof Error ? error.message : String(error)
    });
  }
};

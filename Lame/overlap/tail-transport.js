/* Preview-only boundary transport: exact same marked branch engine, with its UI-only source margin removed. */
const previewTailAttachment=(function(atlasMath){
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
})();


return atlasAttachment;
})({...atlasMath,sourceStatus(r,s){return Number.isFinite(r)&&Number.isFinite(s)&&r>0&&s>0&&r<.5&&s<.5&&r+s>.5?'ok':'outside_source_triangle'}});
previewTailAttachment.setAnchorData(attachmentAnchors);

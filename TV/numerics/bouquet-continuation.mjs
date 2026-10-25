// Browser continuation core. Cusp predictions are never solver inputs.
import {createEvaluator as fixedEvaluator} from './bigint-generic.mjs?v=20261024-tetra-dense4';
import {createEvaluator as decimalEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';
import {createEvaluator as conditionedEvaluator} from './bouquet-fixed.mjs?v=20261024-tv-audit2';
import {compactPrepolish} from './bouquet-compact.mjs?v=20261024-bouquet-stable1';

export const VERSION='tv-bouquet-stable-1';
const engines=new Map(),clock=()=>performance.now(),pause=()=>new Promise(r=>setTimeout(r,0));
export function engine(p=100,jetScale=0){p=Math.max(50,Math.min(950,10*Math.ceil(p/10)));const key=p+':'+jetScale;if(!engines.has(key))engines.set(key,jetScale?conditionedEvaluator(p,{jetScale}):fixedEvaluator(p));return engines.get(key);}
const sourceEngine=decimalEvaluator(350),CH=[[1,0,0,1],[1,-1,1,0],[0,1,-1,1]],PERM=[[0,1,2,3],[0,2,3,1],[0,3,1,2]];
const clone=x=>JSON.parse(JSON.stringify(x));
export const cleanConfig=c=>({id:c.id,seed:'bouquet',core:c.core,weights:[...c.weights],shift:[...c.shift],graph:clone(c.graph)});
const sorted=x=>Array.isArray(x)?x.map(sorted):x&&typeof x==='object'?Object.fromEntries(Object.keys(x).sort().map(k=>[k,sorted(x[k])])):x;
export const configKey=c=>JSON.stringify(sorted(cleanConfig(c)));
const ensure=d=>{if(clock()>d)throw Error('The numerical time budget was reached');};
const inside=s=>Number(s[0])>0&&Number(s[0])<.5&&Number(s[1])>0&&Number(s[1])<.5&&Number(s[0])+Number(s[1])>.5;
const max=(E,a)=>a.reduce((m,x)=>E.abs(x).gt(m)?E.abs(x):m,new E.D(0));

/** Freeze the modular frame, row/column scales and complete-pivot minor.
 * The last pivot is the local Schur residual. Reference evaluation is reused.
 * A minor whose pivot shrinks drastically is rejected rather than mistaken
 * for a new zero. Rebuilding happens between holomorphic local solves only.
 */
function residual(E,config,source,initial,edge,deadline,metrics){
 const [a,b,c,d]=CH[edge],r=new E.D(String(source[1])).plus(config.shift[1]),s=new E.D(String(source[0])).plus(config.shift[0]);
 const rr=r.times(a).minus(s.times(b)),ss=s.times(d).minus(r.times(c)),weights=PERM[edge].map(i=>config.weights[i]);
 const red=E.reduceModular(weights,initial,rr,ss),get=w=>{ensure(deadline);metrics.matrixEvaluations++;return E.tvMatrix(red.weights,E.mobius(w,red.matrix),red.r,red.s,{reduce:false}).M;};
 const raw=get(initial),N=raw.length;if(!N)throw Error('Empty matching matrix');
 const rows=Array.from({length:N},(_,i)=>i),cols=[...rows],rs=raw.map(row=>max(E,row));
 if(rs.some(x=>x.isZero()))throw Error('Zero reference row; increase working precision');
 const first=raw.map((row,i)=>row.map(z=>E.div(z,rs[i]))),cs=Array.from({length:N},(_,j)=>max(E,first.map(row=>row[j])));
 if(cs.some(x=>x.isZero()))throw Error('Zero reference column; increase working precision');
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)first[i][j]=E.div(first[i][j],cs[j]);
 const pivots=[];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=E.abs(first[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const v=E.abs(first[i][k]);if(v.gt(best)){p=i;q=k;best=v;}}
  if(best.isZero())throw Error('Reference minor rounded singular; increase working precision');
  if(p!==j){[first[p],first[j]]=[first[j],first[p]];[rows[p],rows[j]]=[rows[j],rows[p]];}
  if(q!==j){for(let i=0;i<N;i++)[first[i][q],first[i][j]]=[first[i][j],first[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]];}
  pivots.push(E.abs(first[j][j]));
  for(let i=j+1;i<N;i++){const t=E.div(first[i][j],first[j][j]);for(let k=j+1;k<N;k++)first[i][k]=E.sub(first[i][k],E.mul(t,first[j][k]));}
 }
 let firstValue=first[N-1][N-1];
 const evaluate=w=>{
  if(firstValue&&E.abs(E.sub(w,initial)).isZero()){const v=firstValue;firstValue=null;return v;}
  const raw=get(w),M=rows.map(i=>cols.map(j=>E.div(E.div(raw[i][j],rs[i]),cs[j])));
  for(let j=0;j<N-1;j++){
   const v=E.abs(M[j][j]);if(v.isZero()||v.lt(pivots[j].times('1e-10')))throw Error('Frozen Schur minor lost its local pivot health');
   for(let i=j+1;i<N;i++){const t=E.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=E.sub(M[i][k],E.mul(t,M[j][k]));}
  }
  return M[N-1][N-1];
 };
 return{evaluate,dimension:N,frame:red.matrix,referenceMinPivot:pivots.length?E.D.min(...pivots).toString():null};
}

/** A local bounded correction. Precision is TOTAL, with no blanket n-dependent
 * guard. Arithmetic adequacy is determined by perturbed precision repeats.
 * A selected graph root or checked predecessor must supply the initial value.
 */
export function correct({config,source,guess,edge=0,precision=100,digits=25,radius=.08,maxStep=.025,deadline=Infinity,maxIterations=28,maxRestarts=2,jetScale=16,fast=true}={}){
 const begun=clock(),E=engine(Math.max(precision,digits+30),jetScale),initial=E.mobius(E.c(guess),CH[edge]),R=new E.D(radius),H=new E.D(maxStep),tolerance=new E.D(10).pow(-digits),metrics={matrixEvaluations:0,iterations:0,restarts:0,lineTrials:0},history=[];
 let current=initial,lastStep=null,lastValue=null,model=null;
 const finish=(converged,message)=>({converged,status:converged?'ok':'unresolved',tau:E.strs(E.inverse(current,CH[edge])),w:E.strs(current),precision:E.precision,jetScale,source:source.map(String),edge,configKey:configKey(config),requestDigits:digits,pivotHealthy:converged,correction:lastStep?.toString()??null,residual:lastValue?.toString()??null,dimension:model?.dimension,frame:model?.frame,referenceMinPivot:model?.referenceMinPivot,milliseconds:clock()-begun,metrics,history,message,formulaInputsUsed:false});
 try{
  ensure(deadline);if(!initial[1].gt('.00001')||!R.gt(0)||!H.gt(0))throw Error('Invalid local chart neighbourhood');
  let value,derivative;
  const R0=new E.D(String(source[0])),S0=new E.D(String(source[1])),clearance=E.D.min(new E.D('.5').minus(R0),new E.D('.5').minus(S0),R0.plus(S0).minus('.5'));
  if(fast&&clearance.gt('1e-8')){const pre=compactPrepolish({config,source,guess:E.strs(initial),edge,radius:R.toNumber(),maxStep:H.toNumber(),deadline});metrics.prepolish=pre;if(pre.converged){const candidate=E.c(pre.w);if(candidate[1].gt('.00001')&&E.abs(E.sub(candidate,initial)).lte(R))current=candidate;}}
  const rebuild=()=>{
   ensure(deadline);model=residual(E,config,source,current,edge,deadline,metrics);value=model.evaluate(current);
   const h=E.D.min(new E.D('1e-6').times(E.abs(current).plus(1)),R.div(30)),probe=[h,h.times(2)];
   derivative=E.div(E.sub(model.evaluate(E.add(current,probe)),value),probe);
  };
  rebuild();
  for(let iteration=0;iteration<maxIterations;iteration++){
   ensure(deadline);metrics.iterations++;
   if(E.abs(value).isZero())throw Error('Residual rounded to zero; increase working precision');
   if(E.abs(derivative).isZero())throw Error('Local derivative rounded to zero; increase working precision');
   let step=E.div(value,derivative),length=E.abs(step);lastStep=length;lastValue=E.abs(value);
   if(length.lt(tolerance))return finish(true,'Root corrected inside a bounded healthy Schur neighbourhood');
   if(length.gt(H))step=E.mul(step,H.div(length));
   let factor=new E.D(1),accepted=false,next,nextValue,reason=null;
   for(let trial=0;trial<10;trial++){
    ensure(deadline);metrics.lineTrials++;next=E.sub(current,E.mul(step,factor));
    if(next[1].gt('.00001')&&E.abs(E.sub(next,initial)).lte(R)){
     try{nextValue=model.evaluate(next);if(E.abs(nextValue).lt(E.abs(value))){accepted=true;break;}}catch(error){reason=String(error);}
    }
    factor=factor.div(2);
   }
   history.push({iteration,step:length.toString(),factor:factor.toString(),accepted,reason});
   if(!accepted){
    if(metrics.restarts>=maxRestarts)throw Error('No decreasing bounded Schur step');
    metrics.restarts++;rebuild();continue;
   }
   derivative=E.div(E.sub(nextValue,value),E.sub(next,current));current=next;value=nextValue;
  }
  throw Error('Bounded correction iteration budget reached');
 }catch(error){const message=String(error.message??error);return finish(false,/Fixed division by zero/i.test(message)?'Arithmetic denominator rounded to zero; increase working precision ('+message+')':message);}
}

export function chartError(a,b,edge=0,p=180){const E=engine(p);return E.abs(E.sub(E.mobius(E.c(a),CH[edge]),E.mobius(E.c(b),CH[edge])));}
export async function repeatedRoot(config,source,guess,{edge=0,depth=0,precision,digits=25,deadline=Infinity,radius=.04,jetScale=16,lowCandidate}={}){
 let last;const base=precision??Math.max(80,Math.ceil(depth+55));
 for(const extra of [0,30,60,100]){
  ensure(deadline);const p=Math.min(900,base+extra),canReuse=extra===0&&lowCandidate?.converged===true&&lowCandidate.pivotHealthy===true&&lowCandidate.formulaInputsUsed===false&&lowCandidate.configKey===configKey(config)&&lowCandidate.edge===edge&&lowCandidate.jetScale===jetScale&&lowCandidate.precision>=p&&lowCandidate.requestDigits>=digits&&Array.isArray(lowCandidate.source)&&lowCandidate.source.length===2&&lowCandidate.source.every((x,i)=>new sourceEngine.D(x).eq(source[i]))&&lowCandidate.correction!==null&&new sourceEngine.D(lowCandidate.correction).lt(new sourceEngine.D(10).pow(-digits)),low=canReuse?lowCandidate:correct({config,source,guess,edge,precision:p,digits,radius,deadline,jetScale,fast:false});last=low;if(!low.converged)continue;
  const highPrecision=Math.max(p+30,low.precision+30),E=engine(highPrecision+20),probe=E.strs(E.inverse(E.add(E.c(low.w),['1e-7','2e-7']),CH[edge])),high=correct({config,source,guess:probe,edge,precision:highPrecision,digits:digits+8,radius:Math.min(.025,radius),deadline,jetScale,fast:false});last=high;if(!high.converged)continue;
  const error=chartError(low.tau,high.tau,edge,p+80);if(error.gt('1e-18'))continue;
  return{status:'ok',tauExact:high.tau,comparisonTauExact:low.tau,w:high.w,precisionDigits:high.precision,crossPrecisionChartError:error.toString(),correction:high.correction,milliseconds:(canReuse?0:low.milliseconds)+high.milliseconds,evaluations:(canReuse?0:low.metrics.matrixEvaluations)+high.metrics.matrixEvaluations,lowPrecisionReused:canReuse,lowPrecision:low.precision,lowRequestedDigits:low.requestDigits,highPrecisionProvenance:'independently perturbed fresh correction',formulaInputsUsed:false};
 }
 const error=Error('Perturbed precision repeats did not resolve one root: '+(last?.message??'no solve'));error.corrector=last;throw error;
}

function makePath({edge=0,position='.37',start=['.29','.31']}){
 const E=sourceEngine,t=new E.D(String(position));start=start.map(String);if(![0,1,2].includes(edge)||!t.gt(0)||!t.lt(1)||!inside(start))throw Error('Choose an open source edge and interior start');
 const end=[[new E.D('.5'),t.div(2)],[t.div(2),new E.D('.5')],[t.div(2),new E.D('.5').minus(t.div(2))]][edge];
 const source=epsilon=>end.map((x,j)=>x.plus(new E.D(String(epsilon)).times(new E.D(start[j]).minus(x))).toString());
 return{edge,position:String(position),start,end:end.map(String),source,atDepth:d=>source(new E.D(10).pow(-d).toString())};
}
function predict(history,d,edge,p){const E=engine(p),last=history.at(-1),w=E.mobius(E.c(last.tau),CH[edge]);if(history.length<2)return{tau:last.tau,w,displacement:0};const before=history.at(-2),a=E.mobius(E.c(before.tau),CH[edge]),move=E.mul(E.sub(w,a),String((d-last.depth)/(last.depth-before.depth)));return{tau:E.strs(E.inverse(E.add(w,move),CH[edge])),w:E.add(w,move),displacement:E.abs(move).toNumber()};}
function correctionOrThrow(o){const got=correct(o);if(!got.converged){const e=Error(got.message);e.corrector=got;throw e;}return got;}
async function advance(config,history,target,path,deadline,metrics,{strictTracking=false,jetScale=16,targetPrecision=0,targetDigits=18}={}){
 if(history.length===1&&target>history[0].depth&&!strictTracking){
  // Obtain an actual local tangent from a nearby source behind the compact
  // anchor. This avoids repeatedly rejecting a zero-order first predictor.
  // The root is provisional and is never exposed as a displayed sample.
  const at=history[0];let bootstrapped=false;
  for(const d of[-.01,-.005]){
   const source=path.atDepth(d);if(!inside(source))continue;ensure(deadline);
   const got=correct({config,source,guess:at.tau,edge:path.edge,precision:80,digits:18,radius:.04,maxStep:.015,deadline,jetScale});metrics.corrections++;metrics.evaluations+=got.metrics.matrixEvaluations;
   if(got.converged){history.unshift({depth:d,source,tau:got.tau,precision:got.precision});at.nextStep=.035;metrics.bootstrap={direction:'backward',depth:d,source,precision:got.precision};bootstrapped=true;break;}
  }
  if(!bootstrapped){
   const d=Math.min(.001,target/4),source=path.atDepth(d),got=correctionOrThrow({config,source,guess:at.tau,edge:path.edge,precision:80,digits:18,radius:.04,maxStep:.015,deadline,jetScale});metrics.corrections++;metrics.evaluations+=got.metrics.matrixEvaluations;history.push({depth:d,source,tau:got.tau,precision:got.precision,nextStep:Math.min(.02,target-d)});metrics.bootstrap={direction:'forward',depth:d,source,precision:got.precision};
  }
 }
 let step=history.at(-1).nextStep??.035;
 while(history.at(-1).depth<target-1e-12){
  ensure(deadline);const last=history.at(-1);let cap=last.depth<1?.07:last.depth<2?.16:last.depth<3?.35:1;
  if(last.depth>=3&&last.predictorRatio<.04&&history.length>1){const before=history.at(-2),speed=chartError(last.tau,before.tau,path.edge,120).toNumber()/(last.depth-before.depth);cap=Math.max(1,Math.min(4,.08/Math.max(.005,speed)));}
  const d=Math.min(target,last.depth+Math.min(step,cap)),finishing=Math.abs(d-target)<1e-12,internalPrecision=d<12?50:d<20?60:Math.max(70,Math.ceil(d+30)),p=Math.max(internalPrecision,finishing?targetPrecision:0),digits=finishing?Math.max(18,targetDigits):18,source=path.atDepth(d),guess=predict(history,d,path.edge,p+30),radius=Math.min(.15,Math.max(.025,2*guess.displacement+.012));
  let got,ratio=Infinity;
  const run=(at,from,extra={})=>{let last;for(const more of[0,30,60]){const z=correct({config,source:at,guess:from,edge:path.edge,precision:p+more,digits,radius,deadline,jetScale,...extra});metrics.corrections++;metrics.evaluations+=z.metrics.matrixEvaluations;if(z.metrics.prepolish){metrics.prepolishMilliseconds=(metrics.prepolishMilliseconds??0)+z.metrics.prepolish.milliseconds;metrics.prepolishEvaluations=(metrics.prepolishEvaluations??0)+z.metrics.prepolish.evaluations;metrics.prepolishCalls=(metrics.prepolishCalls??0)+1;}last=z;if(z.converged)return z;if(!/precision|rounded|singular|derivative|no decreasing|iteration budget/i.test(z.message))break;}const e=Error(last.message);e.corrector=last;throw e;};
  try{
   got=run(source,guess.tau);const E=engine(p+30),c=E.abs(E.sub(E.c(got.w),guess.w)).toNumber();ratio=c/Math.max(.003,guess.displacement);if(ratio>.25)throw Error('Prediction neighbourhood disagrees');
   if(metrics.steps===0||(strictTracking&&last.depth<.12)||ratio>.06||(last.depth<1.1&&(strictTracking||metrics.steps%4===0))||(d-last.depth>1.1&&metrics.steps%4===0)){
    const mid=(last.depth+d)/2,g1=predict(history,mid,path.edge,p+30),one=run(path.atDepth(mid),g1.tau,{radius:Math.min(.1,radius)}),g2=predict([...history,{depth:mid,tau:one.tau}],d,path.edge,p+30),two=run(source,g2.tau,{radius:Math.min(.1,radius)});metrics.halfChecks++;
    if(chartError(got.tau,two.tau,path.edge,p+50).gt('1e-11'))throw Error('Full and two-half continuations disagree');got=two;
   }
   if(metrics.steps===0||(strictTracking&&last.depth<.12)||ratio>.15||metrics.steps%12===0){
    const R=Math.min(.2,.025+2*chartError(got.tau,last.tau,path.edge,p+30).toNumber()),back=run(last.source,got.tau,{radius:R});metrics.reverseChecks++;if(chartError(back.tau,last.tau,path.edge,p+50).gt('1e-11'))throw Error('Reverse continuation disagrees');
   }
  }catch(error){metrics.rejected++;step=(d-last.depth)/2;if(step<1e-5||clock()>deadline)throw error;continue;}
  const proposed=Math.min(cap,(d-last.depth)*(ratio<.025?1.65:ratio<.1?1.3:1)),nextStep=d===target&&ratio<.1?Math.max(step,proposed):proposed;history.push({depth:d,source,tau:got.tau,precision:got.precision,nextStep,predictorRatio:ratio,corrected:got});if(history.length>4)history.shift();metrics.steps++;step=nextStep;await pause();
 }
 return history.at(-1);
}
const jobs=new Map();let serial=0;
export async function trace(request={}){
 const original=request.config??request.caseConfig,config=cleanConfig(original),edge=Number(request.edge??request.path?.index??0),position=String(request.position??request.path?.position??'.37'),start=request.start??request.path?.start??['.29','.31'],path=makePath({edge,position,start}),depth=Number(request.depth??-Math.log10(Number(request.target??'1e-40'))),begun=clock(),deadline=begun+(request.budgetMs??45000),density=Math.max(1,Math.min(5,Number(request.density??request.samplesPerDecade)||1)),gateEpsilon=String(request.gateEpsilon??'.5'),gateDepth=-Math.log10(Number(gateEpsilon));
 if(!(depth>=0&&depth<=80))throw Error('Choose a depth from zero through80');
 if(!Number.isFinite(Number(gateEpsilon))||!(Number(gateEpsilon)>0&&Number(gateEpsilon)<1))throw Error('The independent graph gate epsilon must lie strictly between zero and one');
 const trackingPolicy={strictTracking:!!request.strictTracking,jetScale:request.jetScale??16,digits:request.digits??25,gateRequired:!!request.gateProvider,gateEpsilon},wanted=request.depths??Array.from({length:Math.ceil(depth*density)+1},(_,i)=>Math.min(depth,i/density));
 if(!Array.isArray(wanted)||!wanted.length||wanted.some((d,i)=>!Number.isFinite(d)||d<0||d>depth||(i>0&&d<=wanted[i-1]))||wanted[0]!==0||Math.abs(wanted.at(-1)-depth)>1e-12)throw Error('Requested depths must increase strictly from zero through the target depth');
 const signature=JSON.stringify(sorted({version:VERSION,config,path:{edge,position,start:path.start},wanted,trackingPolicy}));
 let job=request.resumeToken?jobs.get(request.resumeToken):null;if(request.resumeToken&&(!job||job.signature!==signature||!job.resumable))throw Error('The checkpoint does not match this full labelled graph and path');
 const resumed=!!job;
 if(!job){job={token:'bouquet-numeric-'+(++serial),signature,nodes:[],history:[],checks:[],index:0,resumable:false,metrics:{evaluations:0,corrections:0,steps:0,rejected:0,halfChecks:0,reverseChecks:0},imported:false,gateChecked:!request.gateProvider};
  if(request.extendFrom){const saved=request.extendFrom.result??request.extendFrom,oldRequest=request.extendFrom.request??{},oldKey=saved.caseKey??configKey(oldRequest.config??oldRequest.caseConfig);if(saved.solverVersion!==VERSION||JSON.stringify(sorted(saved.trackingPolicy))!==JSON.stringify(sorted(trackingPolicy)))throw Error('The saved prefix uses another continuation/precision policy');if(!saved.checks?.some(c=>c.kind==='start'&&c.accepted===true&&c.checked===true&&c.formulaInputsUsed===false&&c.configKey===configKey(config)&&Array.isArray(c.source)&&c.source.length===2&&c.source.every((x,i)=>new sourceEngine.D(x).eq(path.start[i]))))throw Error('The saved prefix lacks an accepted formula-free graph anchor');if(oldKey!==configKey(config)||saved.path?.edge!==edge||String(saved.path?.position)!==position||JSON.stringify(saved.path?.start)!==JSON.stringify(path.start))throw Error('The saved trace belongs to another graph/source path');if(saved.status!=='ok'&&!saved.stop?.retryable)throw Error('Only a checked or retryable prefix can be extended');const nodes=clone(saved.nodes??[]);if(!nodes.length||nodes.length>wanted.length||nodes.some(n=>n.status!=='ok'||!Number.isFinite(n.depth)))throw Error('The saved prefix has unresolved, invalid or excessive samples');for(let j=0;j<nodes.length;j++)if(Math.abs(nodes[j].depth-wanted[j])>1e-12)throw Error('The old samples are not a prefix of the requested grid');if(request.gateProvider&&saved.gateChecked&&!saved.checks?.some(c=>c.kind==='intermediate'&&c.accepted===true&&c.checked===true&&c.formulaInputsUsed===false&&c.configKey===configKey(config)&&Array.isArray(c.source)&&c.source.length===2&&c.source.every((x,i)=>new sourceEngine.D(x).eq(path.source(gateEpsilon)[i]))))throw Error('The saved prefix lacks its claimed checked intermediate graph gate');job.nodes=nodes;job.index=nodes.length;job.reusedCount=nodes.length;job.checks=clone(saved.checks??[]);job.history=nodes.slice(-3).map(n=>({depth:n.depth,source:[n.rExact,n.sExact],tau:n.tauExact,precision:n.precisionDigits}));job.imported=true;job.gateChecked=!!saved.gateChecked||!request.gateProvider;}
  jobs.set(job.token,job);while(jobs.size>8)jobs.delete(jobs.keys().next().value);
 }
 const finish=(status,message,stop=null)=>({status,message,nodes:clone(job.nodes),completed:job.nodes.length,total:wanted.length,elapsedMs:clock()-begun,solverVersion:VERSION,caseKey:configKey(config),path:{edge,position,start:path.start,endpoint:path.end},trackingPolicy,checks:clone(job.checks),gateChecked:job.gateChecked,metrics:clone(job.metrics),formulaInputsUsed:false,resumeToken:job.resumable?job.token:null,resumeAvailable:job.resumable,resumed,reusedCount:job.reusedCount??0,newCount:job.nodes.length-(job.reusedCount??0),stop});
 try{
  ensure(deadline);
  if(job.imported){const last=job.nodes.at(-1),source=path.atDepth(last.depth);if(source.some((s,i)=>!new sourceEngine.D(s).eq([last.rExact,last.sExact][i])))throw Error('Saved source coordinates disagree with the exact path');const z=await repeatedRoot(config,source,last.tauExact,{edge,depth:last.depth,deadline,digits:request.digits??25,jetScale:request.jetScale??16});if(chartError(z.tauExact,last.tauExact,edge,180).gt('1e-18'))throw Error('Saved endpoint failed perturbed revalidation');job.history.at(-1).tau=z.tauExact;job.imported=false;}
  if(!job.history.length){const anchor=request.anchor??(request.anchorProvider?await request.anchorProvider(config,path.start,{deadline,correct}):null);if(anchor?.accepted!==true||anchor.checked!==true||anchor.formulaInputsUsed!==false)throw Error('A formula-free accepted and checked graph anchor is required');if(anchor.configKey!==configKey(config))throw Error('The supplied graph anchor belongs to another configuration');if(!Array.isArray(anchor.source)||anchor.source.length!==2||anchor.source.some((s,i)=>!new sourceEngine.D(s).eq(path.start[i])))throw Error('The supplied anchor has another source point');const tau=anchor.tau??anchor.tauExact;job.checks.push({kind:'start',...anchor});job.history=[{depth:0,source:path.start,tau,precision:anchor.precision}];}
  for(;job.index<wanted.length;job.index++){
   if(!job.gateChecked&&wanted[job.index]>=gateDepth){const at=await advance(config,job.history,gateDepth,path,deadline,job.metrics,request),source=path.source(gateEpsilon),gate=await request.gateProvider(config,source,{deadline,correct});job.checks.push({kind:'intermediate',epsilon:gateEpsilon,...gate});if(!gate.accepted){const budget=gate.failureKind==='time_budget'||gate.stopKind==='time_budget'||clock()>deadline||/budget|deadline/i.test(gate.message??''),e=Error(budget?'Independent graph gate reached its budget':'Independent graph gate is inconclusive');e.failureKind=budget?'time_budget':'geometry_inconclusive';throw e;}if(gate.checked!==true||gate.formulaInputsUsed!==false||gate.configKey!==configKey(config)||!Array.isArray(gate.source)||gate.source.length!==2||gate.source.some((s,i)=>!new sourceEngine.D(s).eq(source[i])))throw Error('The intermediate graph gate lacks matching checked formula-free metadata');if(chartError(gate.tau??gate.tauExact,at.tau,edge,180).gt('1e-10'))throw Error('Independent graph gate disagrees with the continued marked branch');job.gateChecked=true;}
   const d=wanted[job.index],outputPrecision=request.digits<=18?Math.max(70,Math.ceil(d+40)):Math.max(80,Math.ceil(d+55)),at=await advance(config,job.history,d,path,deadline,job.metrics,{...request,targetPrecision:outputPrecision,targetDigits:request.digits??25}),source=path.atDepth(d),root=await repeatedRoot(config,source,at.tau,{edge,depth:d,deadline,digits:request.digits??25,jetScale:request.jetScale??16,precision:outputPrecision,lowCandidate:at.corrected}),epsilon=new sourceEngine.D(10).pow(-d).toString();
   const node={status:'ok',depth:d,epsilon,rExact:source[0],sExact:source[1],...root};job.nodes.push(node);job.history.at(-1).tau=root.tauExact;const snapshot=finish('partial','Checked prefix while calculation is running',{kind:'in_progress',retryable:true,branchRejected:false});snapshot.resumeToken=null;snapshot.resumeAvailable=false;request.onProgress?.({node,completed:job.nodes.length,total:wanted.length,continuation:snapshot});await pause();
  }
  job.resumable=false;return finish('ok','All roots continued from the accepted graph anchor; no asymptotic initializer used');
 }catch(error){const budget=error.failureKind==='time_budget'||/budget|deadline/i.test(String(error)),inconclusive=error.failureKind==='geometry_inconclusive',numericalPrecision=/precision|rounded|Fixed division by zero/i.test(String(error)),retryable=budget||inconclusive;job.resumable=retryable;return finish(job.nodes.length?'partial':'unresolved',String(error.message??error),{kind:budget?'time_budget':inconclusive?'geometry_inconclusive':numericalPrecision?'numerical_precision_limit':'verification',retryable,resumeAvailable:retryable,branchRejected:!retryable&&!numericalPrecision,...(error.corrector?{corrector:error.corrector}:{})});}
}
export function clearCaches(){jobs.clear();engines.clear();}

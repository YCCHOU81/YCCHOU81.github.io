// Adapt the bounded numerical core to the site's paired base/decoration API.
// The bounded browser checkpoint includes both branches and exact path/grid.
import * as C from './bouquet-continuation.mjs?v=20261024-tv-audit2';
import {anchor as geometricAnchor} from './bouquet-anchor.mjs?v=20261024-tv-audit2';
import {createEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';
import {traceCheckpoint,traceBudget,traceExtension,extendTraceGrid} from './trace-checkpoints.mjs?v=20261024-audit1';
const E=createEvaluator(350),clock=()=>performance.now(),clone=x=>x==null?null:JSON.parse(JSON.stringify(x));
const inside=start=>Number(start[0])>0&&Number(start[0])<.5&&Number(start[1])>0&&Number(start[1])<.5&&Number(start[0])+Number(start[1])>.5;

export function createBouquetTrace({validate,baseOf,isBase,branch,anchorAt}){
 return async function bouquetTrace(request){
  const c=validate(request.caseConfig),begun=clock(),deadline=begun+traceBudget(request),target=String(request.target??'1e-10'),depth=-Math.log10(Number(target));
  if(!(depth>=1&&depth<=80))throw Error('Choose final epsilon from 0.1 through 1e-80');
  const path=request.path??{},edge=Number(path.index??0),position=String(path.position??'.5'),start=(path.start??['.29','.31']).map(String),t=new E.D(position);
  if((path.kind??'edge')!=='edge'||![0,1,2].includes(edge)||!t.gt(0)||!t.lt(1)||!inside(start))throw Error('Choose an open edge, an interior start, and a position strictly between zero and one');
  const density=Math.max(1,Math.min(5,Number(request.samplesPerDecade)||2)),first=Math.min(1,Math.max(Number(target),Number(request.startEpsilon)||1)),d0=-Math.log10(first),steps=Math.max(1,Math.ceil((depth-d0)*density)),bc=baseOf(c);
  const extension=traceExtension(request,'bouquet');
  const wanted=extension?extendTraceGrid(extension,target,density,d=>new E.D(10).pow(d).toString()):depth-d0<1e-12?[target]:Array.from({length:steps+1},(_,i)=>i===steps?target:new E.D(10).pow(-d0-(depth-d0)*i/steps).toString());
  if(extension&&extension.result.bouquetContinuation?.version!==C.VERSION)throw Error('This saved three-loop trace uses the earlier evaluator. Compute it once with the stable solver before extending it.');
  const {saved,resumed}=traceCheckpoint(request,'bouquet',()=>{
   const continuation=extension?.result.bouquetContinuation,limit=extension?-Math.log10(Number(extension.last.epsilon)):Infinity;
   const prefix=record=>{const copy=clone(record);if(copy&&copy.nodes.some(n=>n.depth>limit+1e-12)){copy.nodes=copy.nodes.filter(n=>n.depth<=limit+1e-12);copy.completed=copy.total=copy.nodes.length;copy.status='ok';copy.stop=null;copy.resumeToken=null;}return copy;};
   return{wanted,base:prefix(continuation?.base),selected:prefix(continuation?.selected),baseLive:false,selectedLive:false,
    extension:extension?{reusedCount:extension.reusedCount}:null,oldNodes:extension?.nodes.slice(0,extension.reusedCount)??[],nodes:[]};
  });
  const state=saved.state,grid=state.wanted,displayDepths=grid.map(v=>-Math.log10(Number(v))),oldDepths=state.selected?.nodes?.map(n=>n.depth)??state.base?.nodes?.map(n=>n.depth)??[];
  const coreDepths=oldDepths.length?[...oldDepths,...displayDepths.filter(d=>d>oldDepths.at(-1)+1e-12)]:[...new Set([0,...displayDepths])];
  const publicNodes=()=>grid.map((epsilon,i)=>{
   if(i<state.oldNodes.length)return state.oldNodes[i];
   const d=displayDepths[i],find=trace=>trace?.nodes?.find(n=>Math.abs(n.depth-d)<1e-12),b=find(state.base),a=isBase(c)?b:find(state.selected),n=a??b;
   return{epsilon,r:n?.rExact,s:n?.sExact,rExact:n?.rExact,sExact:n?.sExact,status:b&&a?'ok':'unresolved',
    base:b?branch(b.tauExact,b):{status:'unresolved'},result:a?branch(a.tauExact,a):{status:'unresolved'},
    branchConnection:b&&a?'Continued from the checked marked graph with bounded steps and perturbed precision checks':null};
  });
  const emit=()=>{state.nodes=publicNodes();const good=state.nodes.filter(n=>n.status==='ok');if(typeof self!=='undefined'&&self.postMessage)self.postMessage({type:'progress',id:request.id,result:{node:good.at(-1),nodes:state.nodes,completed:good.length,total:grid.length,solverVersion:C.VERSION,bouquetContinuation:{version:C.VERSION,base:clone(state.base),selected:clone(isBase(c)?state.base:state.selected)}}});};
  const run=async(cfg,role)=>{
   const old=state[role];if(old?.status==='ok'&&Math.abs(old.nodes.at(-1).depth-depth)<1e-12)return;
   const live=state[role+'Live'],options={config:C.cleanConfig(cfg),edge,position,start,depth,depths:coreDepths,digits:18,jetScale:16,budgetMs:Math.max(1,deadline-clock()),
    anchorProvider:(config,source,settings)=>anchorAt(cfg,source,settings),
    gateProvider:(config,source,settings)=>geometricAnchor(config,source,{...settings,precision:80,digits:18}),gateEpsilon:'.5',
    ...(old?.resumeAvailable&&live?{resumeToken:old.resumeToken}:old?.nodes?.length?{extendFrom:{request:{config:C.cleanConfig(cfg)},result:{...old,resumeToken:null}}}:{}),
    onProgress:p=>{state[role]=p.continuation??{...(state[role]??{}),nodes:[...(state[role]?.nodes??[]),p.node]};emit();}};
   const result=await C.trace(options);state[role]=result;state[role+'Live']=true;emit();
  };
  let stop=null,errorMessage;
  try{
   if(clock()>deadline)throw Error('The numerical time budget was reached');
   await run(bc,'base');
   if(state.base.status!=='ok')stop=state.base.stop;
   else if(isBase(c))state.selected=state.base;
   else{
    if(clock()>deadline)throw Error('The numerical time budget was reached');
    await run(c,'selected');if(state.selected.status!=='ok')stop=state.selected.stop;
   }
  }catch(error){errorMessage=String(error.message??error);const budget=/budget|deadline/i.test(errorMessage),precision=/precision|rounded|division by zero/i.test(errorMessage);stop={kind:budget?'time_budget':precision?'numerical_precision_limit':'numerical',retryable:budget,resumeAvailable:budget,branchRejected:false,message:errorMessage};}
  state.nodes=publicNodes();const completed=state.nodes.filter(n=>n.status==='ok').length,all=completed===grid.length&&state.base?.status==='ok'&&(isBase(c)||state.selected?.status==='ok');
  saved.resumable=!all&&!!stop?.retryable;
  const end=[[new E.D('.5'),t.div(2)],[t.div(2),new E.D('.5')],[t.div(2),new E.D('.5').minus(t.div(2))]][edge];
  return{status:all?'ok':completed?'partial':'unresolved',nodes:clone(state.nodes),completed,total:grid.length,elapsedMs:clock()-begun,
   path:{kind:'edge',index:edge,position,start,endpoint:end.map(String)},solverVersion:C.VERSION,formulaInputsUsed:false,
   stop:all?null:stop,resumeAvailable:saved.resumable,resumeToken:saved.resumable?saved.token:null,resumed,attempts:saved.attempts,
   ...(state.extension?{extended:true,reusedCount:state.extension.reusedCount,newCount:Math.max(0,completed-state.extension.reusedCount)}:{}),
   bouquetContinuation:{version:C.VERSION,base:clone(state.base),selected:clone(isBase(c)?state.base:state.selected)},
   tailBranchValidationStatus:'continued',tailBranchValidation:'Every displayed point continues the marked graph branch; no conjectural cusp initializer is used',
   precisionDescription:'Practical 18-digit corrections with a perturbed higher-precision repeat; stable Fourier evaluation and adaptive working precision. Precision limits keep the accepted prefix.',
   message:all?'Computed graph-selected numerical trajectory. The conjectural cusp formula is a comparison, not a substitute for roots.':saved.resumable?'Time budget reached. Checked points and internal continuation are saved; Continue calculation resumes the same trace.':errorMessage??state.selected?.message??state.base?.message??'Stopped at the checked numerical prefix.'};
 };
}

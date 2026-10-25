// Independent, bounded compact-mesh jobs. Every child keeps the numerical
// anchor, continuation, reverse, and endpoint checks of the serial worker.
import {reusableMeshSamples} from './compact-mesh-reuse.mjs?v=20261024-tv-audit2';
const clock=()=>performance.now();

export function meshPoolSize(total,{supported=typeof Worker==='function',cores=globalThis.navigator?.hardwareConcurrency,memory=globalThis.navigator?.deviceMemory}={}){
 if(total<45||!supported||!Number.isFinite(cores)||cores<2||Number.isFinite(memory)&&memory<=2)return 0;
 return Math.min(cores<6?2:3,Number.isFinite(memory)&&memory<=4?2:3);
}

// Equal-cardinality strips in the original source coordinates give connected
// regions while avoiding one child crossing the whole source triangle.
export function meshPartitions(points,count){
 if(!Number.isInteger(count)||count<2||count>3||points.length<count)throw Error('Invalid compact worker count');
 const order=points.map((p,i)=>({i,r:p.r,s:p.s})).sort((a,b)=>a.r-b.r||a.s-b.s||a.i-b.i);
 return Array.from({length:count},(_,k)=>order.slice(Math.floor(k*order.length/count),Math.floor((k+1)*order.length/count)).map(p=>p.i));
}

export async function pooledCompactMesh({request,points,triangles,resolution,inset,deadline,workerURL,onProgress,signal=request.signal,
 size=meshPoolSize(points.length),workerFactory=url=>new Worker(url,{type:'module'})}={}){
 if(size<2||request._meshPoolChild)return null;
 const begun=clock(),remaining=deadline-begun;
 if(!(remaining>0))return null;
 const samples=reusableMeshSamples(request,points,resolution,inset),initialCount=samples.size,missing=points.map((_,i)=>i).filter(i=>!samples.has(i));
 if(missing.length===1)return null;
 size=Math.min(size,missing.length||size);
 const partitions=missing.length?meshPartitions(missing.map(i=>points[i]),size).map(part=>part.map(i=>missing[i])):[],wallDeadline=Date.now()+remaining,workers=[],done=new Set(),terminated=new Set(),errors=[],metrics={reusedSamples:samples.size};
 let timer,finished=false,startupFailed=false,numericallyStopped=false;
 return new Promise(resolve=>{
  const stopWorker=index=>{if(terminated.has(index))return;const worker=workers[index];if(!worker)return;terminated.add(index);worker.onmessage=null;worker.onerror=null;try{worker.terminate();}catch{}};
  const cleanup=()=>{
   clearTimeout(timer);signal?.removeEventListener?.('abort',aborted);
   workers.forEach((_,index)=>stopWorker(index));
  };
  const finish=(message,{allowSerial=false}={})=>{
   if(finished)return;finished=true;cleanup();
   // Some browsers expose Worker but refuse nested module workers only after
   // construction. Retry serially inside the original deadline in that case.
   // A numerical rejection, cancellation or timeout must never take this path.
   if(allowSerial&&!numericallyStopped&&samples.size===initialCount&&!signal?.aborted&&clock()<deadline){resolve(null);return;}
   const retained=[...samples.keys()].sort((a,b)=>a-b),remap=new Map(retained.map((v,i)=>[v,i])),completed=retained.length;
   resolve({status:completed===points.length?'ok':completed?'partial':'unresolved',samples:retained.map(i=>samples.get(i)),points:retained.map(i=>points[i]),
    triangles:triangles.filter(t=>t.every(i=>remap.has(i))).map(t=>t.map(i=>remap.get(i))),resolution,inset,cached:false,completed,total:points.length,
    metrics,elapsedMs:clock()-begun,pooled:true,poolWorkers:size,poolCompletedWorkers:done.size,
    message:completed===points.length?'Computed image of the compact source triangle; cusp tails are computed separately.':message??errors.find(Boolean)??'The compact image time budget was reached'});
  };
  const aborted=()=>finish('Compact image calculation cancelled');
  const merge=(result,allowed)=>{
   const ids=result.sourceIndices,values=result.samples;
   if(!Array.isArray(ids)||!Array.isArray(values)||ids.length!==values.length)throw Error('A compact worker returned malformed point indices');
   for(let j=0;j<ids.length;j++){
    const id=ids[j],sample=values[j],p=points[id];
    if(!allowed.has(id)||!p||sample?.base?.status!=='ok'||sample?.result?.status!=='ok'||sample?.source?.[0]!==p.r||sample?.source?.[1]!==p.s)
     throw Error('A compact worker returned a point outside its labelled source region');
    samples.set(id,sample);
   }
  };
  const completedWorker=(index,message,{startupFailure=false}={})=>{if(done.has(index))return;done.add(index);stopWorker(index);if(message)errors[index]=message;startupFailed ||= startupFailure;numericallyStopped ||= !!message&&!startupFailure;if(done.size===size)finish(undefined,{allowSerial:startupFailed});};
  try{
   if(signal?.aborted){finish('Compact image calculation cancelled');return;}
   signal?.addEventListener?.('abort',aborted,{once:true});
   if(!missing.length){finish();return;}
   if(samples.size)onProgress?.(samples.size,points.length);
   for(let index=0;index<size;index++){
    const indices=partitions[index],allowed=new Set(indices),worker=workerFactory(workerURL);workers.push(worker);
    const id=String(request.id??'mesh')+'-part-'+index;
    worker.onmessage=event=>{
     if(finished||done.has(index))return;
     const message=event.data;if(message?.id!==id)return;
     try{
      if(message.type==='progress'){
       if(message.result?.poolChunk)merge(message.result.poolChunk,allowed);
       onProgress?.(samples.size,points.length);return;
      }
      if(message.type==='error'){completedWorker(index,message.message??'A compact worker stopped');return;}
      if(message.type!=='mesh')return;
      merge(message.result,allowed);
      for(const [name,value]of Object.entries(message.result.metrics??{}))if(Number.isFinite(value))metrics[name]=(metrics[name]??0)+value;
      onProgress?.(samples.size,points.length);completedWorker(index,message.result.status==='ok'?null:message.result.message);
     }catch(error){completedWorker(index,error.message);}
    };
    worker.onerror=event=>completedWorker(index,event?.message??'A compact worker could not start',{startupFailure:true});
    worker.postMessage({...request,id,_meshPoolChild:true,_meshIndices:indices,_meshWallDeadline:wallDeadline,budgetMs:Math.max(1,wallDeadline-Date.now()),signal:undefined});
   }
   if(!finished)timer=setTimeout(()=>finish('The compact image time budget was reached'),Math.max(0,deadline-clock()));
  }catch(error){
   finish(error.message,{allowSerial:true});
  }
 });
}

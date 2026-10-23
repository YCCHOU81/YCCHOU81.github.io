/* Separate Lamé comparison preview. solver.js is an unchanged atlas worker.
 * One labelled branch per job; cancellation is Worker.terminate().
 */
"use strict";
importScripts("solver.js","tail-solver.js","tail-transport.js","tail-trace.js","symmetry.js");

const previewFiniteTau = branch => branch?.status === "ok" &&
  Number.isFinite(branch.tau?.re) && Number.isFinite(branch.tau?.im) && branch.tau.im > 0;
const previewTau = branch => previewFiniteTau(branch) ? {re:branch.tau.re, im:branch.tau.im} : null;
const previewYield = () => new Promise(resolve => setTimeout(resolve, 0));

self.onmessage = async ({data}) => {
  if (data?.type !== "run") return;
  const id = data.id, started = performance.now();
  try {
    const counts = Array.isArray(data.counts) ? [...data.counts] : null;
    const checked = atlasAttachment.validateCounts(counts);
    if (checked.status !== "ok") throw new RangeError(checked.message);
    const resolution = Math.max(2, Math.min(24, Math.round(Number(data.resolution) || 4)));
    const inset = Number.isFinite(data.inset) ? data.inset : .06;
    const mesh = atlasMath.buildSourceMesh(resolution, inset);
    const timeBudgetMs = Number.isFinite(data.timeBudgetMs) ? Math.max(5000,Math.min(180000,data.timeBudgetMs)) : 45000;
    const existing = data.existingRecord;
    const reuseInterior = existing?.counts?.every((v,i)=>v===counts[i]) && existing.resolution===mesh.resolution &&
      Math.abs(existing.inset-mesh.inset)<1e-12 && existing.tau?.length===mesh.points.length &&
      existing.tau.every(z=>Number.isFinite(z?.re)&&Number.isFinite(z?.im)&&z.im>0) &&
      Number.isFinite(existing.reference?.re)&&Number.isFinite(existing.reference?.im)&&existing.reference.im>0;
    const tau = reuseInterior?existing.tau.map(z=>({...z})):new Array(mesh.points.length).fill(null),
      pointDiagnostics = reuseInterior?(existing.pointDiagnostics??[]).map(p=>({...p,reused:true})):new Array(mesh.points.length);
    const diagnostics = {numericalNotIntervalCertified:true, methods:{}, statuses:{},
      maximumComparisonTauError:0, maximumCorrectionEstimate:0, unresolved:[],
      branchPolicy:"same labelled geometric anchor and continuation; no substitute roots"};
    self.postMessage({type:"start",id,counts,points:mesh.points,triangles:mesh.triangles,
      resolution:mesh.resolution,inset:mesh.inset,total:mesh.points.length});
    let referenceSample = reuseInterior?null:atlasAttachment.solveAttachment(1/3, 1/3, counts, null);
    const reference = reuseInterior?{...existing.reference}:previewTau(referenceSample.result);
    diagnostics.referenceStatus = reuseInterior?"ok":referenceSample.result?.status ?? referenceSample.status;
    diagnostics.referenceElapsedMs = referenceSample?.elapsedMs??0;
    diagnostics.referenceComparisonTauError=referenceSample?.result?.comparisonTauError??null;
    diagnostics.referenceCorrectionEstimate=referenceSample?.result?.correctionEstimate??null;
    diagnostics.reusedInterior=Boolean(reuseInterior);
    let previous = referenceSample?.status === "ok" ? referenceSample : null;
    self.postMessage({type:"reference", id, counts, reference, status:diagnostics.referenceStatus,
      elapsedMs:performance.now()-started});

    // Begin at the nearest mesh point to the reference and then walk to the
    // nearest unvisited point. All mesh vertices remain identical in every case.
    const remaining = new Set(reuseInterior?[]:mesh.points.map((_,i)=>i));
    let cursor = {r:1/3,s:1/3}, completed = reuseInterior?tau.length:0, resolved = reuseInterior?tau.length:0;
    if(reuseInterior)for(let index=0;index<tau.length;index++)self.postMessage({type:"sample",id,counts,index,source:mesh.points[index],tau:tau[index],status:"ok",completed:tau.length,total:tau.length,reused:true,elapsedMs:performance.now()-started});
    while (remaining.size) {
      if(performance.now()-started>timeBudgetMs)break;
      let index = -1, shortest = Infinity;
      for (const candidate of remaining) {
        const p = mesh.points[candidate], d = (p.r-cursor.r)**2+(p.s-cursor.s)**2;
        if (d < shortest) {shortest=d;index=candidate;}
      }
      remaining.delete(index);
      const p = mesh.points[index];
      const sample = atlasAttachment.solveAttachment(p.r,p.s,counts,previous);
      const branch = sample.result, value = previewTau(branch);
      tau[index]=value; completed++; if(value)resolved++;
      const status=branch?.status??sample.status;
      pointDiagnostics[index]={status,elapsedMs:sample.elapsedMs,
        comparisonTauError:branch?.comparisonTauError??null,
        correctionEstimate:branch?.correctionEstimate??null,method:branch?.method??null};
      diagnostics.statuses[status]=(diagnostics.statuses[status]??0)+1;
      if(branch?.method)diagnostics.methods[branch.method]=(diagnostics.methods[branch.method]??0)+1;
      if(Number.isFinite(branch?.comparisonTauError))diagnostics.maximumComparisonTauError=Math.max(diagnostics.maximumComparisonTauError,branch.comparisonTauError);
      if(Number.isFinite(branch?.correctionEstimate))diagnostics.maximumCorrectionEstimate=Math.max(diagnostics.maximumCorrectionEstimate,branch.correctionEstimate);
      if(!value)diagnostics.unresolved.push({index,r:p.r,s:p.s,status});
      // A partially failed BASE solve never contaminates the labelled result.
      // The existing solver only accepts complete successful states as warm starts.
      if(sample.status==="ok")previous=sample;
      cursor=p;
      self.postMessage({type:"sample",id,counts,index,source:{r:p.r,s:p.s},tau:value,
        status,completed,total:mesh.points.length,elapsedMs:performance.now()-started});
      self.postMessage({type:"progress",id,completed,total:mesh.points.length,resolved});
      await previewYield();
    }
    const record={counts,n:checked.n,points:mesh.points.map(({r,s,barycentric})=>({r,s,barycentric})),
      triangles:mesh.triangles,tau,reference,referenceSource:{r:1/3,s:1/3},
      resolution:mesh.resolution,inset:mesh.inset,
      status:reference&&resolved===tau.length?"ok":"partial",resolved,total:tau.length,
      completed,stopReason:remaining.size?"time_budget":null,
      elapsedMs:performance.now()-started,diagnostics:{...diagnostics,uncomputed:remaining.size,timeBudgetMs},pointDiagnostics};
    const depth=Math.max(3,Math.min(20,Math.round(Number(data.depth??data.tailDepth)||5))),density=Math.max(1,Math.min(3,Math.round(Number(data.density??data.tailDensity)||1)));
    const tailStarted=performance.now(),tailBudgetMs=Math.max(5000,Math.min(180000,Number(data.tailBudgetMs)||90000)),epsilons=previewTailTrace.grid(depth,density,mesh.inset);
    record.tails=[];record.tailDepth=depth;record.tailDensity=density;record.transitionVersion='join-v1';
    const stabilizer=LameSymmetry.words.filter(word=>LameSymmetry.applyCounts(counts,word).every((v,i)=>v===counts[i]));
    for(let edge=0;edge<3;edge++)for(const t of [.2,.8]){
      const tailIndex=record.tails.length,existingTail=existing?.tails?.find(tail=>tail.edge===edge&&Math.abs(Number(tail.t??tail.position)-t)<1e-12);
      self.postMessage({type:"tail-start",id,counts,tailIndex,edge,t,total:epsilons.length,depth,density});
      const derivedFrom=record.tails.flatMap(tail=>stabilizer.map(word=>({tail,word,mapped:LameSymmetry.applyEdge(tail.edge,tail.t,word)})))
        .find(item=>item.tail.status==="ok"&&item.mapped.edge===edge&&Math.abs(item.mapped.t-t)<1e-10);
      let tail;
      if(derivedFrom){
        const v=previewApproach.createPrecisionEngine(70),word=derivedFrom.word;
        tail=LameSymmetry.applyTail(derivedFrom.tail,word);tail.edge=edge;tail.t=t;tail.position=t;tail.path={kind:"edge",edge,index:edge,position:String(t)};
        tail.nodes=derivedFrom.tail.nodes.map(node=>{
          let value=v.asComplex(node.tauExact);
          for(const operation of word)value=operation==="N"?[v.d(1).minus(value[0]),value[1]]:v.V(v.C(1),v.S(v.C(1),value));
          const source=previewApproach.approachCoordinates(tail.path,node.epsilon,70);
          return{...node,...source,source:{r:source.r,s:source.s},tau:{re:value[0].toNumber(),im:value[1].toNumber()},tauExact:value.map(z=>z.toPrecision(55)),
            symmetryDerived:true,symmetryWord:word,sourceTauExact:node.tauExact,inheritedRepresentativeError:true};
        });
        tail.symmetryDerived=true;tail.elapsedMs=0;
        for(let index=0;index<tail.nodes.length;index++)self.postMessage({type:"tail-sample",id,counts,tailIndex,edge,t,node:tail.nodes[index],completed:index+1,total:epsilons.length,symmetryDerived:true,elapsedMs:performance.now()-started});
      }else{
        const remainingBudget=Math.max(0,tailBudgetMs-(performance.now()-tailStarted));
        tail=await previewTailTrace.trace(counts,edge,t,epsilons,{existingTail,inset:mesh.inset,timeBudgetMs:remainingBudget,
          onNode(node,tailCompleted,tailTotal){self.postMessage({type:"tail-sample",id,counts,tailIndex,edge,t,node,completed:tailCompleted,total:tailTotal,elapsedMs:performance.now()-started});}});
      }
      record.tails.push(tail);
      self.postMessage({type:"tail-done",id,counts,tailIndex,edge,t,tail,elapsedMs:performance.now()-started});
      await previewYield();
    }
    record.tailStatus=record.tails.every(tail=>tail.status==="ok")?"ok":"partial";
    record.tailResolved=record.tails.reduce((sum,tail)=>sum+tail.resolved,0);
    record.tailTotal=record.tails.reduce((sum,tail)=>sum+tail.total,0);
    record.tailElapsedMs=performance.now()-tailStarted;
    record.elapsedMs=performance.now()-started;
    record.diagnostics.tailPolicy="genuine marked log-epsilon continuation; every plotted node corrected at 45/70 Fourier digits (70/100 beyond depth10); gaps retained; no manufactured endpoints";
    self.postMessage({type:"done",id,record});
  } catch(error) {
    self.postMessage({type:"error",id,message:error instanceof Error?error.message:String(error),elapsedMs:performance.now()-started});
  }
};

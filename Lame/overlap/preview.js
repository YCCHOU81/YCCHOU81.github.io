"use strict";
(() => {
  const $=id=>document.getElementById(id),M=ComparisonMath,V=ViewMath,C=ContourMath;
  const INSET=.06,VERSION="classical-normalized-20261023-v1",TAIL_VERSION="tails-v2";
  const query=new URLSearchParams(location.search),records=new Map(),cards=new Map(),visibleCards=new Set(),hiddenCases=new Set();
  const palette=["#355e86","#9b5639","#557d62","#745f91","#8d7530","#397f86","#986376","#65708c"];
  let n=4,resolution=4,depth=5,density=1,mode="affine",layout="overlay",selected="1,1,1";
  let triples=[],jobs=[],slots=[],epoch=0,running=false,started=0,elapsed=0,totalJobs=0,doneJobs=0,reused=0,meshReused=0;
  let drawn=false,orderTimer=null,prepared=new Map(),pointSelection=null,publishTimer=0,legendElements=new Map();
  const fmt=value=>Number.isFinite(value)?Number(value.toPrecision(8)).toString():"—";
  const key=M.key,label=c=>`(${c.join(", ")})`,canonical=c=>LameSymmetry.canonical(c);
  const canonicalKey=c=>key(canonical(c)),expected=()=> (resolution+1)*(resolution+2)/2;
  const recommended=order=>order<=8?4:order<=12?3:2;
  const baseId=c=>`${VERSION}|${key(c)}|${resolution}|${INSET}`;
  const fullId=c=>`${baseId(c)}|${TAIL_VERSION}|${depth}|${density}`;
  const finite=z=>z&&Number.isFinite(z.re)&&Number.isFinite(z.im);
  const tauText=z=>finite(z)?`${Number(z.re.toPrecision(15))} ${z.im<0?"−":"+"} ${Number(Math.abs(z.im).toPrecision(15))}i`:"unresolved";
  const superscript=value=>String(value).split("").map(c=>({"-":"⁻","0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹"}[c]??c)).join("");

  // Bounded overlap-atlas storage, separate from the single-case atlas cache.
  const cache=(()=>{
    const memory=new Map(),MAX_BYTES=8*1024*1024,MAX_ITEMS=600;let database=null,bytes=0;
    const remember=(id,record,when=Date.now())=>{
      const size=2*(id.length+JSON.stringify(record).length),old=memory.get(id);if(old)bytes-=old.size;
      memory.delete(id);memory.set(id,{record,size,when});bytes+=size;const removed=[];
      while(bytes>MAX_BYTES||memory.size>MAX_ITEMS){const first=memory.keys().next().value;bytes-=memory.get(first).size;memory.delete(first);removed.push(first);}return removed;
    };
    const ready=new Promise(resolve=>{
      if(typeof indexedDB==="undefined"){resolve();return;}let settled=false;
      const finish=()=>{if(!settled){settled=true;resolve();}},timer=setTimeout(finish,1500);
      try{const request=indexedDB.open("lame-atlas-overlap-v1",1);
        request.onupgradeneeded=()=>request.result.createObjectStore("meshes",{keyPath:"key"});
        request.onerror=request.onblocked=()=>{clearTimeout(timer);finish();};
        request.onsuccess=()=>{database=request.result;const get=database.transaction("meshes","readonly").objectStore("meshes").getAll();
          get.onsuccess=()=>{for(const row of get.result.sort((a,b)=>a.when-b.when))if(row?.record?.counts&&Array.isArray(row.record.tau)&&!memory.has(row.key))remember(row.key,row.record,row.when);
            const stale=get.result.filter(row=>!memory.has(row.key));if(stale.length){try{const store=database.transaction("meshes","readwrite").objectStore("meshes");for(const row of stale)store.delete(row.key);}catch(_){}}
            clearTimeout(timer);finish();};get.onerror=()=>{clearTimeout(timer);finish();};};
      }catch(_){clearTimeout(timer);finish();}
    });
    return{ready,get(id){const row=memory.get(id);if(!row)return null;memory.delete(id);memory.set(id,row);return row.record;},
      put(id,record){const removed=remember(id,record);if(!database)return;try{const store=database.transaction("meshes","readwrite").objectStore("meshes");store.put({key:id,record,when:Date.now()});for(const old of removed)store.delete(old);}catch(_){}},get persistent(){return!!database;}};
  })();

  function settings(){const order=Number($("order-n").value);if(!Number.isInteger(order)||order<1||order>21)return null;
    return{n:order,resolution:$("mesh-quality").value==="auto"?recommended(order):Number($("mesh-quality").value),depth:$("tail-depth").value==="auto"?(order<=8?5:3):Number($("tail-depth").value),density:Number($("tail-density").value)};
  }
  function saveURL(){const url=new URL(location.href);for(const [name,value]of Object.entries({n,quality:$("mesh-quality").value,depth:$("tail-depth").value,density,coordinate:mode,layout,case:selected}))url.searchParams.set(name,String(value));url.searchParams.delete("mode");history.replaceState(null,"",url);}
  function hasInterior(record){return!!record?.reference&&record.tau?.length===expected()&&record.tau.every(finite);}
  function tailDepth(record){return Number(record?.tailDepth??record?.depth??0);}
  function requestedEpsilons(){const values=[INSET,.055,.05,.045,.04];for(let exponent=2;exponent<depth-1e-10;exponent+=1/density)values.push(10**(-exponent));values.push(10**(-depth));return values.map(value=>Number(value).toPrecision(11));}
  function tailReady(tail){if(!tail)return false;const available=new Set((tail.nodes??[]).filter(node=>finite(node.tau)&&node.status==="ok").map(node=>Number(node.epsilon).toPrecision(11)));return requestedEpsilons().every(epsilon=>available.has(epsilon));}
  function hasTails(record){return record?.tails?.filter(Boolean).length===6&&record.tails.filter(Boolean).every(tailReady);}
  function cacheRecord(record){if(!record?.counts||!hasInterior(record))return;cache.put(baseId(record.counts),record);if(hasTails(record))cache.put(fullId(record.counts),record);}
  function stop(announce=true){epoch++;clearTimeout(publishTimer);for(const slot of slots){slot.worker.terminate();if(slot.job){const record=records.get(key(slot.job.counts));if(record){record.status=hasInterior(record)?"ok":"stopped";record.tailStatus="stopped";cacheRecord(record);publish(record);}}}
    slots=[];jobs=[];if(running)elapsed=performance.now()-started;running=false;$("stop-compute").disabled=true;$("compute-all").disabled=false;
    if(announce)$("run-status").textContent="Stopped. Computed interiors and tail samples are kept; Compute or extend resumes the case.";
    updateProgress();updateCards();requestDraw();
  }
  function buildCards(){intersection.disconnect();cards.clear();visibleCards.clear();legendElements.clear();$("case-grid").replaceChildren();$("case-legend").replaceChildren();
    const fragment=document.createDocumentFragment(),legend=document.createDocumentFragment();
    triples.forEach((counts,index)=>{const id=key(counts),color=palette[index%palette.length],card=document.createElement("article");card.className="case-card";card.dataset.key=id;card.id=`case-${id.replaceAll(",","-")}`;
      const button=document.createElement("button");button.type="button";button.className="case-select";button.setAttribute("aria-label",`Select attachment ${label(counts)}`);button.innerHTML=`<span class="case-title">${label(counts)}</span><span class="case-state">Waiting</span>`;
      const canvas=document.createElement("canvas");canvas.className="case-canvas";canvas.width=1;canvas.height=1;canvas.setAttribute("aria-label",`Computed interior and tails of attachment ${label(counts)}`);
      const meta=document.createElement("div");meta.className="case-meta";const progress=document.createElement("progress");progress.className="case-progress";progress.max=expected();progress.value=0;
      card.append(button,canvas,meta,progress);fragment.append(card);cards.set(id,{card,button,canvas,meta,progress,state:button.lastElementChild,color});button.addEventListener("click",()=>choose(id));canvas.addEventListener("click",()=>choose(id));intersection.observe(card);
      const chip=document.createElement("button");chip.type="button";chip.className="case-chip";chip.style.setProperty("--case-color",color);chip.textContent=label(counts);chip.title="Select this attachment; Shift-click to hide or show it.";chip.setAttribute("aria-label",`Select ${label(counts)}; shift-click hides it`);chip.addEventListener("click",event=>{if(event.shiftKey){hiddenCases.has(id)?hiddenCases.delete(id):hiddenCases.add(id);prepared.clear();updateCards();requestDraw();}else choose(id);});legend.append(chip);legendElements.set(id,chip);
    });$("case-grid").append(fragment);$("case-legend").append(legend);$("empty-state").hidden=true;updateCards();
  }
  const intersection=new IntersectionObserver(entries=>{for(const entry of entries){const id=entry.target.dataset.key,card=cards.get(id);if(!card)continue;if(entry.isIntersecting){visibleCards.add(id);drawCard(id);}else{visibleCards.delete(id);card.canvas.width=1;card.canvas.height=1;}}},{rootMargin:"160px 0px"});

  async function start(event){event?.preventDefault();clearTimeout(orderTimer);if(!$("compare-form").reportValidity())return;const next=settings();if(!next)return;
    stop(false);const token=epoch;({n,resolution,depth,density}=next);records.clear();prepared.clear();hiddenCases.clear();pointSelection=null;triples=M.enumerate(n);
    if(!triples.some(c=>key(c)===selected))selected=key(triples.reduce((a,b)=>{const score=c=>c.reduce((sum,v)=>sum+(v-(n-1)/3)**2,0);return score(a)<score(b)?a:b;}));
    overview.manual=false;comparison.manual=false;buildCards();saveURL();$("run-status").textContent="Restoring computed interiors and tails…";await cache.ready;if(token!==epoch)return;
    const independent=new Map();for(const counts of triples)independent.set(canonicalKey(counts),canonical(counts));totalJobs=independent.size;doneJobs=0;reused=0;meshReused=0;jobs=[];
    const selectedCanonical=canonicalKey(selected.split(",").map(Number));
    const ordered=[...independent.values()].sort((a,b)=>key(a)===selectedCanonical?-1:key(b)===selectedCanonical?1:Math.max(...b)-Math.max(...a));
    for(const counts of ordered){const old=cache.get(fullId(counts))??cache.get(baseId(counts));if(old){publish(old);if(hasInterior(old)&&hasTails(old)){doneJobs++;reused++;continue;}if(hasInterior(old))meshReused++;}jobs.push({counts,existingRecord:old??null});}
    started=performance.now();elapsed=0;running=jobs.length>0;$("compute-all").disabled=running;$("stop-compute").disabled=!running;updateProgress();updateCards();showSelection();requestDraw();
    if(!running){finish();return;}
    const count=Math.min(jobs.length,Math.max(1,Math.min(3,(navigator.hardwareConcurrency||4)-1)));
    for(let i=0;i<count;i++){const worker=new Worker("normalized-worker.js?v=3"),slot={worker,job:null,sequence:0};slots.push(slot);worker.onmessage=event=>receive(slot,event.data,token);worker.onerror=()=>{if(token===epoch&&slot.job)failJob(slot,"Numerical worker stopped.");};dispatch(slot,token);}
  }
  function dispatch(slot,token){if(token!==epoch)return;const job=jobs.shift();if(!job){slot.job=null;if(slots.every(s=>!s.job))finish();return;}slot.job=job;slot.sequence++;
    const id=key(job.counts),old=job.existingRecord,record=old?{...old,counts:job.counts,status:hasInterior(old)?"ok":"computing",tailStatus:"computing"}:{counts:job.counts,n,resolution,inset:INSET,points:[],triangles:[],tau:Array(expected()).fill(null),reference:null,total:expected(),resolved:0,completed:0,status:"computing",tails:[],tailStatus:"computing"};
    records.set(id,record);updateCard(id);slot.worker.postMessage({type:"run",id:`${token}:${slot.sequence}:${id}`,counts:job.counts,resolution,inset:INSET,depth,density,existingRecord:old,timeBudgetMs:90000,tailBudgetMs:depth>10?180000:90000});
  }
  function receive(slot,data,token){if(token!==epoch||!slot.job||data.id!==`${token}:${slot.sequence}:${key(slot.job.counts)}`)return;const id=key(slot.job.counts),record=records.get(id);
    if(data.type==="start"){Object.assign(record,{points:data.points,triangles:data.triangles,total:data.total});if(!hasInterior(record))record.tau=Array(data.total).fill(null);}
    else if(data.type==="reference"){record.reference=data.reference;prepared.clear();}
    else if(data.type==="sample"){record.tau[data.index]=data.tau;record.completed=data.completed;record.resolved=record.tau.filter(finite).length;if(record.resolved===record.total)record.status="ok";prepared.clear();updateCard(id);requestDraw();}
    else if(data.type==="tail-start"){record.tails??=[];const saved=record.tails[data.tailIndex];record.tails[data.tailIndex]={...saved,edge:data.edge,t:data.t??data.position,position:data.t??data.position,nodes:saved?.nodes?.slice()??[],total:data.total,status:"computing"};record.tailStatus="computing";}
    else if(data.type==="tail-sample"){const tail=record.tails[data.tailIndex];if(tail){const epsilon=Number(data.node.epsilon).toPrecision(12),index=tail.nodes.findIndex(node=>Number(node.epsilon).toPrecision(12)===epsilon);if(index<0)tail.nodes.push(data.node);else tail.nodes[index]=data.node;tail.nodes.sort((a,b)=>Number(b.epsilon)-Number(a.epsilon));tail.resolved=tail.nodes.filter(node=>finite(node.tau)).length;tail.completed=data.completed;}prepared.clear();throttledPublish(record);requestDraw();}
    else if(data.type==="tail-done"){record.tails[data.tailIndex]=data.tail;record.tailDepth=depth;record.tailDensity=density;cacheRecord(record);publish(record);updateCards();updateProgress();showSelection();requestDraw();}
    else if(data.type==="done"){const result=data.record;cacheRecord(result);publish(result);doneJobs++;slot.job=null;updateCards();updateProgress();showSelection();requestDraw();dispatch(slot,token);}
    else if(data.type==="error")failJob(slot,data.message);
  }
  function throttledPublish(record){if(publishTimer)return;publishTimer=setTimeout(()=>{publishTimer=0;if(records.get(key(record.counts))===record){publish(record);updateCards();updateProgress();}},120);}
  function failJob(slot,message){const record=records.get(key(slot.job.counts));record.status=hasInterior(record)?"ok":"error";record.tailStatus="partial";record.message=message;cacheRecord(record);publish(record);doneJobs++;slot.job=null;updateCards();updateProgress();requestDraw();dispatch(slot,epoch);}
  function publish(record){for(const counts of triples){if(canonicalKey(counts)!==canonicalKey(record.counts))continue;try{const derived=key(counts)===key(record.counts)?record:LameSymmetry.derive(record,counts);records.set(key(counts),derived);}catch(_){/* A not-yet-built source mesh is not guessed. */}}prepared.clear();$("download-data").disabled=!records.size;}
  function finish(){if(running)elapsed=performance.now()-started;running=false;for(const slot of slots)slot.worker.terminate();slots=[];jobs=[];$("compute-all").disabled=false;$("stop-compute").disabled=true;updateCards();updateProgress();requestDraw();showSelection();
    const complete=[...records.values()].filter(hasTails).length;$("run-status").textContent=complete===triples.length?`All ${complete} interiors and six-series tails are ready.`:`${complete}/${triples.length} complete tail sets. Other resolved samples remain visible; Compute or extend retries the gaps.`;
  }
  function clippedRecord(record){const tails=(record.tails??[]).filter(Boolean).map(tail=>({...tail,nodes:(tail.nodes??[]).filter(node=>{const e=Number(node.epsilon);return e>=10**(-depth)*(1-1e-9);})}));return{...record,tails};}
  function viewRecord(id,coordinate=mode){const token=`${coordinate}|${id}`;if(prepared.has(token))return prepared.get(token);const record=records.get(id);if(!record)return null;const result=V.prepare(clippedRecord(record),{mode:coordinate});prepared.set(token,result);return result;}
  function updateCard(id){const card=cards.get(id),record=records.get(id);if(!card)return;const count=record?.resolved??0,total=record?.total??expected(),tails=record?.tails?.filter(Boolean)??[];
    card.card.classList.toggle("is-selected",id===selected);card.button.setAttribute("aria-pressed",String(id===selected));card.card.classList.toggle("is-partial",!!record&&(!hasInterior(record)||!hasTails(record))&&!running);
    card.state.textContent=hasTails(record)?"Ready":record?.tailStatus==="computing"?"Computing tails":hasInterior(record)?"Interior ready":record?.status==="error"?"Unresolved":"Waiting";
    card.meta.textContent=`${count}/${total} interior · ${tails.reduce((sum,t)=>sum+(t.nodes??[]).filter(node=>finite(node.tau)).length,0)} tail points`;
    card.progress.max=6;card.progress.value=tails.filter(tailReady).length;card.progress.hidden=hasTails(record);const chip=legendElements.get(id);if(chip){chip.classList.toggle("is-selected",id===selected);chip.classList.toggle("is-hidden",hiddenCases.has(id));chip.setAttribute("aria-pressed",String(id===selected));}
  }
  function updateCards(){for(const id of cards.keys())updateCard(id);}
  function updateProgress(){const interiors=[...records.values()].filter(hasInterior).length,tails=[...records.values()].filter(hasTails).length;
    $("summary-count").textContent=`${interiors}/${triples.length} interiors · ${tails}/${triples.length} tail sets`;$("gallery-count").textContent=`${triples.length} attachments`;
    $("run-progress").max=Math.max(totalJobs,1);$("run-progress").value=doneJobs;$("runtime").textContent=started?`${((running?performance.now()-started:elapsed)/1000).toFixed(1)} s`:"";
    if(running)$("run-status").textContent=`Computing interiors and six cusp approaches · ${doneJobs}/${totalJobs} independent cases complete.`;
    $("cache-note").textContent=`${totalJobs} independent computations cover ${triples.length} ordered triples. ${reused?`${reused} complete cases restored. `:""}${meshReused?`${meshReused} saved interiors reused. `:""}Overlap Atlas cache: at most 8 MB.`;
    const nodes=[...records.values()].flatMap(record=>clippedRecord(record).tails.flatMap(tail=>tail.nodes??[])),resolved=nodes.filter(node=>finite(node.tau)).length;
    $("tail-status").textContent=`${resolved} resolved tail points · target ε = 10${superscript(-depth)} · two paths at each of the three cusps. ${nodes.length>resolved?`${nodes.length-resolved} unresolved points are left as gaps.`:"Dashed curves stop at computed endpoints; sample markers are optional."}`;
  }
  function choose(id){if(!cards.has(id))return;const old=selected;selected=id;pointSelection=null;hiddenCases.delete(id);const priority=canonicalKey(id.split(",").map(Number));jobs.sort((a,b)=>canonicalKey(a.counts)===priority?-1:canonicalKey(b.counts)===priority?1:0);updateCard(old);updateCard(id);showSelection();saveURL();requestDraw();}
  function showSelection(){const record=records.get(selected),counts=selected.split(",").map(Number),value=pointSelection?.id===selected?pointSelection.tau:record?.reference;
    $("selection-info").textContent=`Attachment ${label(counts)} · n = ${n}\n${pointSelection?.id===selected?"Selected sample":"Reference"}: τ = ${tauText(value)}${pointSelection?.epsilon?`\nε = ${pointSelection.epsilon}`:""}\n${(record?.tails??[]).filter(Boolean).filter(tailReady).length}/6 completed cusp series`;
  }

  // A single Euclidean scale is used on both axes. Zooming never changes Im τ alone.
  function createPlot(canvas){const plot={canvas,cx:.5,cy:1,scale:100,manual:false,map:null,unmap:null,geometry:[],samples:[],width:0,height:0};let drag=null,moved=false;
    canvas.addEventListener("wheel",event=>{event.preventDefault();if(!plot.unmap)return;const box=canvas.getBoundingClientRect(),x=event.clientX-box.left,y=event.clientY-box.top,point=plot.unmap(x,y);plot.scale=Math.max(1e-7,Math.min(1e13,plot.scale*Math.exp(-event.deltaY*.0015)));plot.cx=point.re-(x-plot.width/2)/plot.scale;plot.cy=point.im+(y-plot.height/2)/plot.scale;plot.manual=true;requestDraw();},{passive:false});
    canvas.addEventListener("pointerdown",event=>{drag={x:event.clientX,y:event.clientY,cx:plot.cx,cy:plot.cy};moved=false;canvas.setPointerCapture(event.pointerId);});
    canvas.addEventListener("pointermove",event=>{const box=canvas.getBoundingClientRect(),x=event.clientX-box.left,y=event.clientY-box.top;if(drag){const dx=event.clientX-drag.x,dy=event.clientY-drag.y;if(Math.hypot(dx,dy)>3)moved=true;plot.cx=drag.cx-dx/plot.scale;plot.cy=drag.cy+dy/plot.scale;plot.manual=true;requestDraw();}else if(plot===comparison)inspectOverlap(plot,x,y);});
    canvas.addEventListener("pointerup",event=>{drag=null;if(moved)return;const box=canvas.getBoundingClientRect(),x=event.clientX-box.left,y=event.clientY-box.top;let nearest=null,best=144;for(const item of plot.samples){const d=(item.x-x)**2+(item.y-y)**2;if(d<best){nearest=item;best=d;}}if(nearest){choose(nearest.id);pointSelection=nearest;showSelection();requestDraw();return;}const hits=hitCases(plot,x,y);if(hits.length)choose(hits.includes(selected)?hits[(hits.indexOf(selected)+1)%hits.length]:hits[0]);});
    canvas.addEventListener("pointercancel",()=>{drag=null;});canvas.addEventListener("dblclick",()=>{plot.manual=false;requestDraw();});return plot;
  }
  const overview=createPlot($("overview")),comparison=createPlot($("comparison-plot"));
  function canvasContext(canvas){const width=canvas.clientWidth||400,height=canvas.clientHeight||400,dpr=Math.min(devicePixelRatio||1,1.7);if(canvas.width!==Math.round(width*dpr)||canvas.height!==Math.round(height*dpr)){canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);}const ctx=canvas.getContext("2d");ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);return{ctx,width,height};}
  function boundsOf(views,coordinate){let minRe=Infinity,maxRe=-Infinity,minIm=Infinity,maxIm=-Infinity;const include=z=>{if(finite(z)){minRe=Math.min(minRe,z.re);maxRe=Math.max(maxRe,z.re);minIm=Math.min(minIm,z.im);maxIm=Math.max(maxIm,z.im);}};
    for(const view of views){if(!view)continue;if($("show-interior").checked)for(const z of view.points??[])include(z);if($("show-tails").checked)for(const tail of view.tails??[])for(const z of tail.values??[])include(z);}
    if(!Number.isFinite(minRe))return{minRe:0,maxRe:1,minIm:0,maxIm:2};if(coordinate==="raw"||coordinate==="cuspinfty"){minRe=Math.min(0,minRe);maxRe=Math.max(1,maxRe);minIm=0;maxIm=Math.max(1.5,maxIm);}if(coordinate==="disk"||coordinate==="common-disk"){minRe=Math.min(-1,minRe);maxRe=Math.max(1,maxRe);minIm=Math.min(-1,minIm);maxIm=Math.max(1,maxIm);}return{minRe,maxRe,minIm,maxIm};
  }
  function setupPlot(plot,views,coordinate){const result=canvasContext(plot.canvas),{width,height}=result;plot.width=width;plot.height=height;
    if(!plot.manual){const box=boundsOf(views,coordinate),dx=Math.max(1e-8,box.maxRe-box.minRe),dy=Math.max(1e-8,box.maxIm-box.minIm);plot.cx=(box.minRe+box.maxRe)/2;plot.cy=(box.minIm+box.maxIm)/2;plot.scale=Math.min((width-85)/dx,(height-70)/dy)*.9;}
    plot.map=z=>[width/2+(z.re-plot.cx)*plot.scale,height/2-(z.im-plot.cy)*plot.scale];plot.unmap=(x,y)=>({re:plot.cx+(x-width/2)/plot.scale,im:plot.cy-(y-height/2)/plot.scale});plot.geometry=[];plot.samples=[];return result;
  }
  function tickStep(scale){const rough=65/scale,power=10**Math.floor(Math.log10(rough)),ratio=rough/power;return(ratio<2?1:ratio<5?2:5)*power;}
  function axes(ctx,plot,coordinate){const {width,height,map,unmap}=plot,step=tickStep(plot.scale),lo=unmap(40,height-28),hi=unmap(width-14,20);ctx.fillStyle="#fcfdfb";ctx.fillRect(0,0,width,height);ctx.font="10px system-ui";ctx.lineWidth=.6;ctx.strokeStyle="#e6eae9";ctx.fillStyle="#71808a";
    let count=0;for(let x=Math.ceil(lo.re/step)*step;x<=hi.re&&count++<70;x+=step){const p=map({re:x,im:0})[0];ctx.beginPath();ctx.moveTo(p,15);ctx.lineTo(p,height-28);ctx.stroke();ctx.textAlign="center";ctx.fillText(axisNumber(x,step),p,height-12);}
    count=0;for(let y=Math.ceil(lo.im/step)*step;y<=hi.im&&count++<70;y+=step){const p=map({re:0,im:y})[1];ctx.beginPath();ctx.moveTo(40,p);ctx.lineTo(width-10,p);ctx.stroke();ctx.textAlign="left";ctx.fillText(axisNumber(y,step),4,p+3);}
    const desc=V.describe(coordinate);ctx.textAlign="right";ctx.fillText(desc.xLabel,width-12,height-2);ctx.textAlign="left";ctx.fillText(desc.yLabel,5,12);
    ctx.save();ctx.beginPath();ctx.rect(40,18,width-50,height-47);ctx.clip();
    if(coordinate==="raw"||coordinate==="cuspinfty"||coordinate==="cusp0"||coordinate==="cusp1"){const heightTop=unmap(0,18).im,offset=coordinate==="cusp0"?-1:0;ctx.strokeStyle="#9baea3";ctx.lineWidth=1;ctx.setLineDash([]);ctx.beginPath();ctx.moveTo(...map({re:offset,im:heightTop}));ctx.lineTo(...map({re:offset,im:0}));for(let j=0;j<=180;j++){const x=j/180;ctx.lineTo(...map({re:x+offset,im:Math.sqrt(x*(1-x))}));}ctx.lineTo(...map({re:1+offset,im:heightTop}));ctx.stroke();}
    if(coordinate==="disk"||coordinate==="common-disk"||coordinate.startsWith("q")){ctx.strokeStyle="#b8c4c0";ctx.beginPath();const center=map({re:0,im:0});ctx.arc(...center,plot.scale,0,2*Math.PI);ctx.stroke();}
    if(coordinate==="common-disk"){const reference={re:.5,im:Math.sqrt(3)/2},boundary=[];for(let j=1;j<=150;j++){const y=Math.tan((j/151)*Math.PI/2);boundary.push(V.transform({re:0,im:y},reference,"common-disk"));}boundary.push({re:1,im:0});for(let j=150;j>=1;j--){const y=Math.tan((j/151)*Math.PI/2);boundary.push(V.transform({re:1,im:y},reference,"common-disk"));}for(let j=179;j>=1;j--){const x=j/180;boundary.push(V.transform({re:x,im:Math.sqrt(x*(1-x))},reference,"common-disk"));}polyline(ctx,plot,boundary,"#a9b9af",1,true);ctx.fillStyle="#63776b";ctx.font="11px system-ui";for(const [text,z]of [["∞",{re:1,im:0}],["0",{re:-.5,im:Math.sqrt(3)/2}],["1",{re:-.5,im:-Math.sqrt(3)/2}]]){const p=map(z);ctx.fillText(text,p[0]+5,p[1]-5);}}
    ctx.restore();
  }
  function axisNumber(value,step){if(Math.abs(value)<step*1e-7)value=0;if(Math.abs(value)>=1e5||(value!==0&&Math.abs(value)<1e-4))return value.toExponential(1);return Number(value.toFixed(Math.max(0,-Math.floor(Math.log10(step))+1))).toString();}
  function polyline(ctx,plot,values,color,width,dashed=false){ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dashed?[4,3]:[]);ctx.beginPath();let active=false;for(const value of values){if(!finite(value)){active=false;continue;}const p=plot.map(value);if(active)ctx.lineTo(...p);else ctx.moveTo(...p);active=true;}ctx.stroke();ctx.setLineDash([]);}
  function curvePixels(values,plot,closed=false){return C.smoothSegments(values.map(plot.map),{closed});}
  function segmentPixels(start,segment,steps=8){const points=[];for(let j=1;j<=steps;j++){const t=j/steps,u=1-t;points.push([u*u*u*start[0]+3*u*u*t*segment.cp1[0]+3*u*t*t*segment.cp2[0]+t*t*t*segment.end[0],u*u*u*start[1]+3*u*u*t*segment.cp1[1]+3*u*t*t*segment.cp2[1]+t*t*t*segment.end[1]]);}return points;}
  function flattenCurve(curve){if(!curve.start)return[];const pixels=[curve.start];let start=curve.start;for(const segment of curve.segments){pixels.push(...segmentPixels(start,segment));start=segment.end;}return pixels;}
  function appendArea(path,plot,pixels,id){if(pixels.length<3||pixels.some(p=>!Number.isFinite(p[0])||!Number.isFinite(p[1])))return;let area=0;for(let i=0,j=pixels.length-1;i<pixels.length;j=i++)area+=pixels[j][0]*pixels[i][1]-pixels[i][0]*pixels[j][1];const ordered=area<0?pixels.slice().reverse():pixels;
    path.moveTo(...ordered[0]);for(let i=1;i<ordered.length;i++)path.lineTo(...ordered[i]);path.closePath();const xs=ordered.map(p=>p[0]),ys=ordered.map(p=>p[1]);plot.geometry.push({id,points:ordered,minX:Math.min(...xs),maxX:Math.max(...xs),minY:Math.min(...ys),maxY:Math.max(...ys)});
  }
  function compactArcPixels(curve,arc){if(!curve?.start||arc.startIndex<0||arc.endIndex<0||!arc.complete)return null;const count=curve.segments.length,starts=[curve.start,...curve.segments.slice(0,-1).map(segment=>segment.end)];if(arc.startIndex>=count||arc.endIndex>=count)return null;
    const pixels=[starts[arc.startIndex]];let index=arc.startIndex,guard=0;while(index!==arc.endIndex&&guard++<=count){if(arc.direction>0){pixels.push(...segmentPixels(starts[index],curve.segments[index]));index=(index+1)%count;}else{const previous=(index-1+count)%count,segment=curve.segments[previous];pixels.push(...segmentPixels(starts[index],{cp1:segment.cp2,cp2:segment.cp1,end:starts[previous]}));index=previous;}}return index===arc.endIndex?pixels:null;
  }
  function strokeCurve(ctx,plot,values,color,width,closed=false){if(values.length<2)return;const curve=curvePixels(values,plot,closed);if(!curve.start)return;ctx.beginPath();ctx.moveTo(...curve.start);for(const segment of curve.segments)ctx.bezierCurveTo(...segment.cp1,...segment.cp2,...segment.end);if(curve.closed)ctx.closePath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.lineJoin="round";ctx.lineCap="butt";ctx.setLineDash([7,4]);ctx.stroke();ctx.setLineDash([]);}
  function prepareDrawing(plot,id,view,{card=false}={}){const record=records.get(id);if(!record||!view)return null;const selectedCase=id===selected,color=selectedCase?"#bb4138":cards.get(id)?.color??palette[0],width=selectedCase?2:1.05;
    const showInterior=$("show-interior").checked,showTails=$("show-tails").checked,showSamples=$("show-samples").checked,contour=C.contours(record,view,{showInterior,showTails}),areaPath=new Path2D();
    const compactCurve=contour.compactBoundary?curvePixels(contour.compactBoundary,plot,true):null;
    if(showInterior){if(compactCurve&&view.points.length&&view.points.every(finite))appendArea(areaPath,plot,flattenCurve(compactCurve),id);else for(const triangle of record.triangles??[]){const vertices=triangle.map(i=>view.points[i]);if(vertices.every(finite))appendArea(areaPath,plot,vertices.map(plot.map),id);}}
    const tails=view.tails??[];
    if(showTails&&$("show-ribbons").checked){for(let edge=0;edge<3;edge++){const indexed=tails.map((tail,index)=>({tail,index})).filter(item=>item.tail.edge===edge).sort((a,b)=>Number(a.tail.t??a.tail.position)-Number(b.tail.t??b.tail.position));if(indexed.length<2)continue;const [left,right]=indexed,lookup=new Map(right.tail.nodes.map(node=>[Number(node.epsilon).toPrecision(12),node.value])),runs=[];let run=[];
      for(const node of left.tail.nodes){const other=lookup.get(Number(node.epsilon).toPrecision(12));if(!finite(node.value)||!finite(other)){if(run.length)runs.push(run);run=[];}else run.push({left:node.value,right:other,epsilon:Number(node.epsilon)});}if(run.length)runs.push(run);
      for(const pairRun of runs){if(pairRun.length<2)continue;const leftPixels=flattenCurve(curvePixels(pairRun.map(pair=>pair.left),plot)),rightPixels=flattenCurve(curvePixels(pairRun.map(pair=>pair.right),plot));let cap=[];if(Math.abs(pairRun[0].epsilon-INSET)<1e-8){const arc=contour.baseArcs.find(arc=>arc.edge===edge&&arc.tailIndices[0]===left.index&&arc.tailIndices[1]===right.index);if(arc&&compactCurve)cap=compactArcPixels(compactCurve,arc)?.slice().reverse()??[];}appendArea(areaPath,plot,[...leftPixels,...rightPixels.reverse(),...cap],id);}
    }}
    // Keep one union path per attachment: its own mesh and ribbon joins count once.
    return {areaPath,color,details(ctx){
    ctx.globalAlpha=selectedCase?1:.85;for(let i=0;i<contour.paths.length;i++)strokeCurve(ctx,plot,contour.paths[i],color,width,contour.pathClosed[i]);ctx.globalAlpha=1;
    if(showTails)for(const tail of tails){const values=tail.values??tail.nodes.map(node=>node.value);for(let j=0;j<values.length;j++){const value=values[j];if(!finite(value))continue;const p=plot.map(value);if(p[0]<35||p[0]>plot.width||p[1]<15||p[1]>plot.height-25)continue;if(showSamples){ctx.beginPath();ctx.arc(...p,selectedCase?1.8:1.15,0,2*Math.PI);ctx.fillStyle=color;ctx.fill();}if(!card)plot.samples.push({id,x:p[0],y:p[1],tau:tail.nodes[j].tau,epsilon:tail.nodes[j].epsilon});}
      if(showSamples){const last=values.findLast(finite);if(last){ctx.beginPath();ctx.arc(...plot.map(last),selectedCase?3.5:2.2,0,2*Math.PI);ctx.fillStyle="#fff";ctx.fill();ctx.strokeStyle=color;ctx.lineWidth=1;ctx.stroke();}}}
    if(selectedCase&&finite(view.reference)){const p=plot.map(view.reference);if(showSamples){ctx.beginPath();ctx.arc(...p,3,0,2*Math.PI);ctx.fillStyle=color;ctx.fill();}if(!card)plot.samples.push({id,x:p[0],y:p[1],tau:record.reference});}
    if(pointSelection?.id===id){const value=V.transform(pointSelection.tau,record.reference,view.mode,{scale:view.scale});if(finite(value)){ctx.beginPath();ctx.arc(...plot.map(value),4.5,0,2*Math.PI);ctx.strokeStyle=color;ctx.lineWidth=2;ctx.stroke();}}
    }};
  }
  function fillDrawings(ctx,drawings){if(!drawings.length)return;const union=new Path2D();for(const drawing of drawings)union.addPath(drawing.areaPath);
    // The base tint is drawn once for the union, independently of overlap count.
    ctx.fillStyle=drawings.length>1||layout==="overlay"?"#3e6483":drawings[0].color;ctx.globalAlpha=Number($("image-opacity").value);ctx.fill(union,"nonzero");ctx.globalAlpha=1;
  }
  function renderPlot(plot,coordinate,ids){const views=ids.map(id=>viewRecord(id,coordinate)),{ctx,width,height}=setupPlot(plot,views,coordinate);axes(ctx,plot,coordinate);ctx.save();ctx.beginPath();ctx.rect(40,18,width-50,height-47);ctx.clip();
    const drawings=ids.filter(id=>id!==selected).concat(ids.includes(selected)?[selected]:[]).map(id=>prepareDrawing(plot,id,viewRecord(id,coordinate))).filter(Boolean);
    fillDrawings(ctx,drawings);OverlapShading.paint(ctx,drawings.map(drawing=>drawing.areaPath),{width,height,strength:Number($("overlap-darkness").value)});
    for(const drawing of drawings)drawing.details(ctx);ctx.restore();}
  function drawCard(id){if(layout!=="gallery"||!visibleCards.has(id))return;const card=cards.get(id);if(!card)return;const {ctx,width,height}=canvasContext(card.canvas),view=viewRecord(id),views=triples.map(c=>viewRecord(key(c))).filter(Boolean),bounds=boundsOf(views,mode);
    const scale=Math.min((width-30)/Math.max(1e-8,bounds.maxRe-bounds.minRe),(height-25)/Math.max(1e-8,bounds.maxIm-bounds.minIm))*.88,cx=(bounds.minRe+bounds.maxRe)/2,cy=(bounds.minIm+bounds.maxIm)/2;
    const plot={canvas:card.canvas,width,height,scale,map:z=>[width/2+(z.re-cx)*scale,height/2-(z.im-cy)*scale],geometry:[],samples:[]};ctx.fillStyle="#fff";ctx.fillRect(0,0,width,height);const drawing=prepareDrawing(plot,id,view,{card:true});if(drawing){fillDrawings(ctx,[drawing]);drawing.details(ctx);}}
  function pointInPolygon(x,y,points){let inside=false;for(let i=0,j=points.length-1;i<points.length;j=i++){const a=points[i],b=points[j];if((a[1]>y)!==(b[1]>y)&&x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])inside=!inside;}return inside;}
  function hitCases(plot,x,y){const hits=new Set();for(const item of plot.geometry)if(!hits.has(item.id)&&x>=item.minX&&x<=item.maxX&&y>=item.minY&&y<=item.maxY&&pointInPolygon(x,y,item.points))hits.add(item.id);return[...hits];}
  function inspectOverlap(plot,x,y){if(!plot.unmap)return;const point=plot.unmap(x,y),hits=hitCases(plot,x,y),desc=V.describe(mode);$("overlap-info").textContent=`${desc.xLabel} = ${fmt(point.re)}, ${desc.yLabel} = ${fmt(point.im)} · ${hits.length} sampled region${hits.length===1?"":"s"} here${hits.length?`: ${hits.slice(0,6).map(id=>`(${id})`).join(" ")}${hits.length>6?" …":""}`:""}. Numerical mesh/ribbon overlap, not a certified domain count.`;}
  function requestDraw(){if(drawn)return;drawn=true;requestAnimationFrame(()=>{drawn=false;const ids=triples.map(key).filter(id=>records.has(id)&&!hiddenCases.has(id));renderPlot(overview,"raw",ids);renderPlot(comparison,mode,layout!=="overlay"?ids.filter(id=>id===selected):ids);if(layout==="gallery")for(const id of visibleCards)drawCard(id);});}
  function updateDisplay(){mode=$("coordinate-view").value;layout=$("display-layout").value;prepared.clear();comparison.manual=false;$("comparison").dataset.layout=layout;$("gallery-panel").hidden=layout!=="gallery";$("comparison-title").textContent=layout==="overlay"?"All images overlaid":layout==="selected"?"Selected image":"Selected image · gallery below";const desc=V.describe(mode);$("coordinate-help").textContent=desc.description;$("comparison-caption").textContent=`${desc.label}. Smooth dashed interpolation through computed points; no extrapolated cusp ends. Wheel or +/− to zoom; drag to pan. Darker shading indicates more overlapping sampled images.`;$("overview-caption").textContent="Original Re τ and Im τ, with equal units on both axes. No height compression. Zoom or pan to inspect narrow cusp tails.";saveURL();requestDraw();showSelection();}
  for(const [prefix,plot]of [["overview",overview],["comparison",comparison]]){const fit=$(prefix==="overview"?"fit-overview":"fit-comparison");fit?.addEventListener("click",()=>{plot.manual=false;requestDraw();});$(`${prefix}-zoom-in`)?.addEventListener("click",()=>{plot.scale*=1.6;plot.manual=true;requestDraw();});$(`${prefix}-zoom-out`)?.addEventListener("click",()=>{plot.scale/=1.6;plot.manual=true;requestDraw();});}
  $("reset-overview")?.addEventListener("click",()=>{overview.cx=.5;overview.cy=1;overview.scale=Math.min((overview.width-80)/1.4,(overview.height-60)/2);overview.manual=true;requestDraw();});
  $("compare-form").addEventListener("submit",start);$("stop-compute").addEventListener("click",()=>stop());
  $("order-n").addEventListener("input",()=>{clearTimeout(orderTimer);if(settings())orderTimer=setTimeout(()=>start(),650);});$("order-n").addEventListener("change",()=>{clearTimeout(orderTimer);const next=settings();if(next&&next.n!==n)start();});
  $("mesh-quality").addEventListener("change",()=>start());
  for(const id of ["tail-depth","tail-density"])$(id).addEventListener("change",()=>{const next=settings();if(next&&next.depth>10){stop(false);depth=next.depth;density=next.density;prepared.clear();updateCards();updateProgress();saveURL();requestDraw();$("run-status").textContent="Deeper tails use more computation. Click Compute images & tails to extend the saved paths; Stop keeps resolved points.";}else start();});
  for(const id of ["coordinate-view","display-layout"])$(id).addEventListener("change",updateDisplay);
  for(const id of ["show-interior","show-tails","show-ribbons","show-samples"])$(id).addEventListener("input",()=>requestDraw());
  for(const id of ["image-opacity","overlap-darkness"])$(id).addEventListener("input",()=>{$(`${id}-value`).textContent=`${Math.round(Number($(id).value)*100)}%`;requestDraw();});
  $("show-selected").addEventListener("click",()=>{layout="selected";$("display-layout").value=layout;updateDisplay();});
  $("download-data").addEventListener("click",()=>{const data={description:"Classical Lamé Overlap Atlas; original tau values with computed paired cusp approaches, no forced endpoints",n,resolution,inset:INSET,depth,density,coordinate:mode,records:[...records.values()].map(clippedRecord)};const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:"application/json"}));const a=document.createElement("a");a.href=url;a.download=`lame-n${n}-images-and-tails.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
  new ResizeObserver(()=>requestDraw()).observe($("comparison"));setInterval(()=>{if(running)updateProgress();},600);window.addEventListener("pagehide",()=>stop(false));
  if(Number.isInteger(Number(query.get("n")))&&Number(query.get("n"))>=1&&Number(query.get("n"))<=21)$("order-n").value=query.get("n");
  if(["auto","2","3","4","6","10"].includes(query.get("quality")))$("mesh-quality").value=query.get("quality");
  if(["auto","3","5","8","10","15","20"].includes(query.get("depth")))$("tail-depth").value=query.get("depth");if(["1","2","3"].includes(query.get("density")))$("tail-density").value=query.get("density");
  if(["raw","affine","equal-affine","disk","common-disk","cusp0","cusp1","cuspinfty","q0","q1","qinfty"].includes(query.get("coordinate")))$("coordinate-view").value=query.get("coordinate");
  if(["overlay","selected","gallery"].includes(query.get("layout")))$("display-layout").value=query.get("layout");if(/^\d+,\d+,\d+$/.test(query.get("case")??""))selected=query.get("case");
  updateDisplay();start();
})();

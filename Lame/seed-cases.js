"use strict";

/* Lamé attachment controls and linked dual graphs.
 * This atlas contains one family. Numerical work belongs to the workers;
 * this module only renders the case, moves controls, and links graph focus.
 * Source-chart helpers remain available to the page controllers.
 */
window.SeedCases = (() => {
  const query = new URLSearchParams(location.search);
  const countNames = ["a", "b", "c"];
  const title = "Classical Lamé · attachments (a,b,c)";
  const shift = [0, 0];
  const top = .5, maxR = .5;
  const active = {
    id: "lame-attachment", seed: "lame", operation: "attachment",
    weights: [1, 0, 0, 0], shift, top, maxR, caseConfig: null,
    isRepeated: false, title, baseTitle: title, baseId: "lame-seed",
    start: [.31, .29], sourceVertices: [[.5, 0], [0, .5], [.5, .5]],
    // Compatibility metadata; classical plots use AtlasView's marked chart.
    edgeLabels: ["R → 1 · cusp ∞", "S → ½ · cusp 0", "R + 2S → 1 · cusp ½"]
  };
  const escape = value => String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll('"', "&quot;");
  const isLimitsPage = () => Boolean(document.getElementById("limits-controls"));
  const compactVertices = inset => active.sourceVertices.map(v => v.map((x,j) => (1-3*inset)*x+inset*active.sourceVertices.reduce((sum,w)=>sum+w[j],0)));
  function project(R,S,inset=.04) {
    const vertices=compactVertices(inset);
    if (![R,S].every(Number.isFinite)) return {r:active.start[0],s:active.start[1]};
    if (R <= maxR*(1-inset) && S <= top*(1-inset) && R/maxR+S/top >= 1+inset) return {r:R,s:S};
    let best=null, distance=Infinity;
    vertices.forEach((a,i)=>{
      const b=vertices[(i+1)%3],dx=b[0]-a[0],dy=b[1]-a[1];
      const t=Math.max(0,Math.min(1,((R-a[0])*dx+(S-a[1])*dy)/(dx*dx+dy*dy)));
      const r=a[0]+t*dx,s=a[1]+t*dy,d=(r-R)**2+(s-S)**2;
      if(d<distance){best={r,s};distance=d;}
    });
    return best;
  }
  function map(R,S,attached=false) {return {r:R+(attached?shift[0]:0),s:S+(attached?shift[1]:0)};}
  const scale = 300;
  function svgPoint(R,S,attached=false) {const p=map(R,S,attached);return [48+scale*p.r,350-scale*p.s];}
  function svgInverse(x,y,attached=false) {return {r:(x-48)/scale-(attached?shift[0]:0),s:(350-y)/scale-(attached?shift[1]:0)};}
  function renderMonodromy(svg,{point,attached=false,title="Monodromy triangle"}={}) {
    const color=attached?"#b55e3c":"#285f65", xy=p=>svgPoint(...p,attached).join(",");
    const xMax=maxR+shift[0],yMax=Math.max(1,top+shift[1]);
    const ticks=max=>Array.from({length:Math.floor(max*2)+1},(_,i)=>i/2);
    const number=n=>({0:"0",.5:"½",1:"1",1.5:"1½"}[n]??String(n));
    svg.setAttribute("viewBox","0 0 400 400");
    svg.innerHTML=`<title>${escape(title)}</title><desc>Dashed full source triangle; shaded compact region with each barycentric coordinate at least 0.04. Coordinates are lifted ordinary characters.</desc>
      <path d="M48 ${350-scale*yMax-12} V350 H${48+scale*xMax+12}" fill="none" stroke="#a5afa2" stroke-width="1.4"/>
      ${ticks(xMax).map(n=>`<text x="${48+scale*n}" y="376" text-anchor="middle">${number(n)}</text>`).join("")}
      ${ticks(yMax).filter(n=>n>0).map(n=>`<text x="37" y="${356-scale*n}" text-anchor="end">${number(n)}</text>`).join("")}
      <text x="${48+scale*xMax+18}" y="355">R</text><text x="38" y="${350-scale*yMax-18}">S</text>
      <polygon class="source-full" points="${active.sourceVertices.map(xy).join(" ")}" fill="${color}06" stroke="${color}" stroke-dasharray="6 5" stroke-width="1.8"/>
      <polygon class="source-compact ${attached?"attached-source-compact":""}" points="${compactVertices(.04).map(xy).join(" ")}" fill="${color}20" stroke="${color}" stroke-width="1.4"/>
      ${point?`<circle class="seed-source-halo" cx="${svgPoint(...point,attached)[0]}" cy="${svgPoint(...point,attached)[1]}" r="14" fill="${color}20"/><circle class="seed-source-dot" cx="${svgPoint(...point,attached)[0]}" cy="${svgPoint(...point,attached)[1]}" r="7" fill="${color}" stroke="white" stroke-width="2"/>`:""}`;
  }

  const graphPalette = ["#2c6b70", "#b15b31", "#6972a8"];
  let selectedGraphControl = null;
  const linkTokens = value => String(value ?? "").split(/\s+/).filter(Boolean);
  const linkAttributes=(ids,name,count,{keyboard=true,related=[]}={})=>ids?.length?` class="graph-control-link" data-graph-controls="${escape(ids.join(" "))}"${related.length?` data-graph-related="${escape(related.join(" "))}"`:""} data-edge="${escape(name)}" role="button" tabindex="${keyboard?0:-1}" aria-controls="${escape(ids.join(" "))}" aria-pressed="false" aria-label="${escape(name)}: ${count} edge${count===1?"":"s"}. Select its decoration control"`:"";
  function linkedPath(d,{name,m,color,controls=[],related=[],width=2.8,title=""}){
    const path=`<path class="graph-chain" d="${d}" fill="none" stroke="${color}" stroke-width="${width}"><title>${escape(title||`${name} ×${m}`)}</title></path>`;
    if(!controls.length)return path;
    return `<g${linkAttributes(controls,name,m,{related})}><path class="graph-chain-halo" d="${d}" fill="none" stroke="#e9bc63" stroke-width="10" opacity="0" pointer-events="none"/>${path}<path class="graph-chain-hit" d="${d}" fill="none" stroke="transparent" stroke-width="14" pointer-events="stroke"/></g>`;
  }
  function controlLabel(control) {
    const label = document.querySelector(`label[for="${control.id}"]`);
    return (label?.textContent ?? control.name ?? control.id).replace(/\s+/g, " ").trim();
  }
  function graphControlVisible(control) {
    return Boolean(control && !control.disabled && !control.closest("[hidden]"));
  }
  function refreshGraphControlLinks() {
    const graph = document.getElementById("case-graph");
    if (!graph) return;
    const selected = document.getElementById(selectedGraphControl);
    if (!graphControlVisible(selected)) selectedGraphControl = null;
    for (const element of document.querySelectorAll(".graph-control-selected")) {
      element.classList.remove("graph-control-selected");
    }
    for (const element of graph.querySelectorAll("[data-graph-controls]")) {
      const controls = [...linkTokens(element.dataset.graphControls), ...linkTokens(element.dataset.graphRelated)];
      const match = Boolean(selectedGraphControl && controls.includes(selectedGraphControl));
      element.classList.toggle("graph-link-selected", match);
      element.setAttribute("aria-pressed", String(match));
    }
    if (selectedGraphControl) {
      const field = selected.closest(".repeat-local") ?? selected.closest(".repeat-count-field") ?? selected.parentElement;
      field?.classList.add("graph-control-selected");
    }
    const status = document.getElementById("graph-control-feedback");
    if (status) {
      const message = selectedGraphControl
        ? `Selected ${controlLabel(selected)}. Apply to update τ.`
        : "Click an edge or count to highlight its control.";
      if (status.textContent !== message) status.textContent = message;
    }
  }
  function selectGraphControl(control,{focus=false}={}){
    if(!graphControlVisible(control))return;
    selectedGraphControl=control.id;refreshGraphControlLinks();
    if(focus){
      control.focus?.({preventScroll:true});
      const rect=control.getBoundingClientRect?.();
      if(rect&&(rect.top<110||rect.bottom>window.innerHeight-24))control.scrollIntoView?.({block:"center",behavior:window.matchMedia?.("(prefers-reduced-motion: reduce)")?.matches?"auto":"smooth"});
    }
  }
  function activateGraphLink(link){
    const control=linkTokens(link.dataset.graphControls).map(id=>document.getElementById(id)).find(graphControlVisible);
    if(control)selectGraphControl(control,{focus:true});
  }
  function linkedFormControl(target) {
    const control = target?.closest?.("input, select");
    if (!control) return null;
    const prefix = isLimitsPage() ? "limit-" : "";
    return countNames.some(name => control.id === prefix + name) ? control : null;
  }
  function clearGraphControl(){
    if(!selectedGraphControl)return;
    selectedGraphControl=null;refreshGraphControlLinks();
  }
  document.addEventListener("click",event=>{
    const link=event.target?.closest?.("#case-graph [data-graph-controls]");
    if(link){event.preventDefault();activateGraphLink(link);return;}
    const label=event.target?.closest?.("label[for]");
    const control=linkedFormControl(event.target)??linkedFormControl(label?document.getElementById(label.getAttribute("for")):null);
    if(control){selectGraphControl(control);return;}
    clearGraphControl();
  });
  document.addEventListener("keydown",event=>{
    if(event.key==="Escape"){clearGraphControl();return;}
    const link=event.target?.closest?.("#case-graph [data-graph-controls]");
    if(link&&["Enter"," "].includes(event.key)){event.preventDefault();activateGraphLink(link);}
  });
  for(const type of ["focusin","pointerdown","input","change"])document.addEventListener(type,event=>{const control=linkedFormControl(event.target);if(control)selectGraphControl(control);});
  const sampleCurve=(a,b,controls=[])=>Array.from({length:49},(_,i)=>{
    const t=i/48,u=1-t;
    if(controls.length===1)return [0,1].map(k=>u*u*a[k]+2*u*t*controls[0][k]+t*t*b[k]);
    if(controls.length===2)return [0,1].map(k=>u*u*u*a[k]+3*u*u*t*controls[0][k]+3*u*t*t*controls[1][k]+t*t*t*b[k]);
    return [0,1].map(k=>u*a[k]+t*b[k]);
  });
  function edgeCountLabels(entries,{vertices=[],occupied=[],box=[0,-5,255,130]}={}){
    const overlap=(a,b,pad=0)=>Math.max(0,Math.min(a[0]+a[2]/2+pad,b[0]+b[2]/2)-Math.max(a[0]-a[2]/2-pad,b[0]-b[2]/2))*Math.max(0,Math.min(a[1]+a[3]/2+pad,b[1]+b[3]/2)-Math.max(a[1]-a[3]/2-pad,b[1]-b[3]/2));
    const pointDistance=(p,a,b)=>{const dx=b[0]-a[0],dy=b[1]-a[1],t=Math.max(0,Math.min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/(dx*dx+dy*dy||1)));return Math.hypot(p[0]-a[0]-t*dx,p[1]-a[1]-t*dy);};
    const curveDistance=(p,points)=>Math.min(...points.slice(1).map((b,i)=>pointDistance(p,points[i],b)));
    const reserved=[...occupied,...vertices.map(p=>[...p,15,15])],out=[];
    // Short central chains get their space first; the large outer loops have
    // many safe choices. The count is the actual number of suppressed edges.
    const length=e=>e.points.slice(1).reduce((s,p,i)=>s+Math.hypot(p[0]-e.points[i][0],p[1]-e.points[i][1]),0);
    const ordered=entries.map((e,index)=>({...e,index})).sort((a,b)=>length(a)-length(b));
    for(const e of ordered){
      const text=`×${e.m}`,width=10+text.length*5.8,height=15;let best=null;
      for(const t of [.5,.35,.65,.2,.8]){
        const i=Math.round(t*(e.points.length-1)),target=e.points[i],before=e.points[Math.max(0,i-1)],after=e.points[Math.min(e.points.length-1,i+1)],dx=after[0]-before[0],dy=after[1]-before[1],norm=Math.hypot(dx,dy)||1;
        for(const offset of [12,-12,20,-20,30,-30,42,-42]){
          const p=[target[0]-dy/norm*offset,target[1]+dx/norm*offset],rect=[...p,width,height];
          if(p[0]-width/2<box[0]+3||p[0]+width/2>box[0]+box[2]-3||p[1]-height/2<box[1]+3||p[1]+height/2>box[1]+box[3]-3)continue;
          let score=Math.abs(offset)*.28+Math.abs(t-.5)*20;
          for(const r of reserved)score+=overlap(rect,r,3)*100;
          for(const other of entries){const d=curveDistance(p,other.points);score+=Math.max(0,11-d)**2*(other===entries[e.index]?1:4);}
          if(!best||score<best.score)best={p,target,rect,score};
        }
      }
      if(!best)continue;
      reserved.push(best.rect);const [x,y]=best.p,[tx,ty]=best.target,color=e.color??"#496351";
      const linked=Boolean(e.controls?.length),attrs=linked?linkAttributes(e.controls,e.name,e.m,{keyboard:false,related:e.related}):` data-edge="${escape(e.name)}" pointer-events="none"`;
      out.push(`<g${attrs.replace('class="graph-control-link"','class="graph-edge-count graph-control-link"')}${linked?"":' class="graph-edge-count"'} data-count="${e.m}"><title>${escape(e.name)}: ${e.m} edge${e.m===1?"":"s"} in this chain${linked?"; select its decoration control":""}</title><path d="M${tx.toFixed(2)} ${ty.toFixed(2)} L${x.toFixed(2)} ${y.toFixed(2)}" fill="none" stroke="${color}" stroke-width=".85" opacity=".8"/><rect x="${(x-width/2).toFixed(2)}" y="${(y-height/2).toFixed(2)}" width="${width}" height="${height}" rx="4" fill="#fcfcf7" stroke="${color}" stroke-width=".8"/><text x="${x.toFixed(2)}" y="${(y+3.4).toFixed(2)}" text-anchor="middle" fill="${color}" font-family="system-ui,sans-serif" font-size="10" font-weight="650">${text}</text></g>`);
    }
    return out.join("");
  }
  function graphKey(items){
    return `<ul class="graph-chain-key" aria-label="Labelled chains and edge counts">${items.map((item,i)=>{const content=`<span class="graph-key-stroke" style="--chain-color:${graphPalette[i%graphPalette.length]}" aria-hidden="true"></span><span><strong>${escape(item.name)} ×${item.m}</strong>${item.faces?`<small>F${item.faces.join(" / F")}</small>`:""}</span>`;return item.controls?.length?`<li><button type="button" class="graph-key-link graph-control-link" data-graph-controls="${escape(item.controls.join(" "))}"${item.related?.length?` data-graph-related="${escape(item.related.join(" "))}"`:""} data-edge="${escape(item.name)}" aria-controls="${escape(item.controls.join(" "))}" aria-pressed="false" aria-label="${escape(item.name)}: select its decoration control">${content}</button></li>`:`<li>${content}</li>`;}).join("")}</ul>`;
  }
  // Classical edge attachments grow bridges in the exterior face only.
  // The terminal loops remain one edge: face degrees are (2n+1,1,1,1).
  // This differs from inserting pairs on a loop, which changes two weights.
  function lameGraph(counts=[0,0,0],{interactive=true}={}){
    const c=Array.isArray(counts)?counts:[counts.a,counts.b,counts.c];
    const values=c.map(x=>Math.max(0,Math.trunc(Number(x)||0))),center=[144,104];
    const angles=[-Math.PI*2/3,0,Math.PI*2/3];
    let paths="",nodes=`<circle cx="${center[0]}" cy="${center[1]}" r="4.6" fill="#285f65" stroke="white"/>`,faces='<text x="267" y="24">F₀</text>';const countEntries=[],vertices=[center],occupied=[[267,20,19,14]];
    for(let i=0;i<3;i++){
      const name=["a","b","c"][i],ids=interactive?[(isLimitsPage()?"limit-":"")+name]:[];
      const u=[Math.cos(angles[i]),Math.sin(angles[i])],v=[-u[1],u[0]],length=values[i]?41:0;
      const end=center.map((x,k)=>x+length*u[k]),p=end.map((x,k)=>x+39*u[k]);
      if(values[i]){
        paths+=linkedPath(`M${center} L${end}`,{name:name+" arm",m:values[i],color:graphPalette[i],controls:ids,width:3,title:`${name} arm: ${values[i]} edges, both sides in F₀`});
        nodes+=`<circle cx="${end[0]}" cy="${end[1]}" r="4.2" fill="#285f65" stroke="white"/>`;
        vertices.push(end);countEntries.push({name:name+" arm",m:values[i],points:sampleCurve(center,end),color:graphPalette[i],controls:ids});
      }
      const c1=end.map((x,k)=>x+80*u[k]-38*v[k]),c2=end.map((x,k)=>x+80*u[k]+38*v[k]);
      const loopIDs=values[i]?[]:ids;
      paths+=linkedPath(`M${end} C${c1} ${c2} ${end}`,{name:`F${i+1} loop`,m:1,color:"#496351",controls:loopIDs,width:2.4,title:`Single-edge loop bounding F${i+1}. ${values[i]?"The terminal loop remains one edge.":`Select ${name} to extend the arm at this loop; the loop remains one edge.`}`});
      faces+=`<text x="${p[0]}" y="${p[1]+4}">F${i+1}</text>`;
      occupied.push([p[0],p[1],19,14]);countEntries.push({name:`F${i+1} loop`,m:1,points:sampleCurve(end,end,[c1,c2]),color:"#496351",controls:loopIDs});
    }
    const items=values.map((m,i)=>({name:["a","b","c"][i]+(m?" arm":" unchanged"),m,controls:interactive?[(isLimitsPage()?"limit-":"")+["a","b","c"][i]]:[]}));
    return `<svg class="labelled-dual-svg lame-dual-svg" viewBox="25 -10 270 230" role="${interactive?"group":"img"}" aria-label="Classical Lamé attachment (${values.join(", ")}), ${1+values.reduce((a,b)=>a+b,0)} hemispheres, weights (${1+values.reduce((a,b)=>a+b,0)},0,0,0)${interactive?". Select an arm, or an unchanged loop, to focus its attachment count":""}"><g>${paths}${nodes}</g><g class="case-face-labels" text-anchor="middle">${faces}</g>${edgeCountLabels(countEntries,{vertices,occupied,box:[25,-10,270,230]})}</svg>${graphKey(items)}`;
  }
  let lameGraphState={counts:[query.get("a")??0,query.get("b")??0,query.get("c")??0],draft:false};
  function updateLameGraph(counts,{draft=false}={}){
    const values=(Array.isArray(counts)?counts:[counts.a,counts.b,counts.c]).map(x=>Math.max(0,Math.trunc(Number(x)||0)));
    lameGraphState={counts:values,draft};
    const target=document.getElementById("lame-result-graph"),caption=document.getElementById("lame-graph-caption");
    if(target)target.innerHTML=lameGraph(values);
    if(caption)caption.innerHTML=`<strong>${draft?"Preview":"Selected"}: (${values.join(", ")})</strong> · n = ${1+values.reduce((a,b)=>a+b,0)} · weights (${1+values.reduce((a,b)=>a+b,0)}, 0, 0, 0)${draft?" · apply to update τ":""}`;
    refreshGraphControlLinks();
  }

  function renderSupport() {
    const support = document.getElementById("case-support");
    if (!support) return;
    const instruction = isLimitsPage()
      ? "Apply decorations first, choose the source path, then use Compute approach."
      : "Apply &amp; compute starts the interior image and cusp series. Refine either preview with the controls below.";
    support.innerHTML = `<details id="lame-numerical-range-help" class="case-rules control-help"><summary>Numerical range &amp; speed</summary><p>Nonnegative whole numbers with a + b + c ≤ 20; n = 1 + a + b + c.</p><p class="caption">These are drawing limits, not mathematical restrictions. Start with small counts for a faster comparison. ${instruction}</p></details>`;
  }
  function renderSummary() {
    renderSupport();
    const summary = document.getElementById("case-summary");
    if (!summary) return;
    const wasOpen = summary.querySelector("details")?.open;
    summary.innerHTML = `<p class="case-active-title"><span class="case-badge">CLASSICAL LAMÉ</span> <strong>${escape(active.title)}</strong></p><p class="case-essential">Attachment counts: <strong>a + b + c ≤ 20</strong>.</p><details id="lame-attachment-help" class="case-rules control-help"><summary>How attachments work</summary><p>The labelled graph changes with these counts. Each count extends the corresponding arm; the terminal loops stay single. Compare all triples of one order in the Overlap Atlas.</p><p><a href="mathematics.html" data-page="mathematics.html">Mathematical background →</a></p></details>`;
    if (wasOpen) summary.querySelector("details").open = true;
  }
  function renderCaseGraph() {
    const graph = document.getElementById("case-graph");
    if (!graph) return;
    const expanded = graph.classList.contains("case-graph-expanded");
    graph.innerHTML = `<div class="case-graph-toolbar"><strong>Labelled dual graphs</strong><button type="button" id="case-graph-enlarge" aria-expanded="false" aria-controls="case-graph-comparison">Enlarge graphs</button></div><div id="case-graph-comparison" class="case-graph-comparison"><figure class="case-graph-base">${lameGraph([0,0,0], {interactive:false})}<figcaption>Base n = 1 · three-loop graph</figcaption></figure><span class="case-graph-arrow" aria-hidden="true">→</span><figure class="case-graph-result"><div id="lame-result-graph">${lameGraph(lameGraphState.counts)}</div><figcaption id="lame-graph-caption"><strong>Selected:</strong> (${lameGraphState.counts.join(", ")})</figcaption></figure></div><details id="lame-dual-graph-help" class="case-rules control-help"><summary>Read the dual graph</summary><p class="case-graph-note">Dots are vertices; F labels identify faces. <strong>Each nearby ×m badge gives the total number of edges in that chain</strong>; ×1 is a single edge. Intermediate bivalent vertices are suppressed. a, b, c count bridge edges towards F₁, F₂, F₃; the terminal loops stay single. The drawing records topology, not edge lengths.</p><p>Click an edge, count badge, or key to highlight its attachment control. Focus a control to highlight its edge. Counts a, b, c extend the arms towards F₁, F₂, F₃; each terminal loop stays one edge. The base graph is a fixed reference.</p></details><p id="graph-control-feedback" class="graph-control-feedback" role="status" aria-live="polite">Click an edge or count to highlight its control.</p>`;
    const toggle = document.getElementById("case-graph-enlarge");
    if (expanded) {
      toggle.setAttribute("aria-expanded", "true");
      toggle.textContent = "Compact view";
    }
    toggle.addEventListener("click", () => {
      const isExpanded = graph.classList.toggle("case-graph-expanded");
      toggle.setAttribute("aria-expanded", String(isExpanded));
      toggle.textContent = isExpanded ? "Compact view" : "Enlarge graphs";
    });
    refreshGraphControlLinks();
  }
  function placeAttachmentControls() {
    const host = document.getElementById("lame-attachment-settings");
    if (!host) return;
    for (const id of ["attachment-form", "decoration-status"]) {
      const element = document.getElementById(id);
      if (element) { element.hidden = false; host.appendChild(element); }
    }
    const legacy = document.getElementById("limits-legacy-attachment");
    if (legacy) legacy.hidden = false;
  }
  function mount() {
    const selector = document.getElementById("case-selector");
    if (!selector) return;
    if (!document.getElementById("atlas-seed")) {
      selector.innerHTML = `<div><label for="atlas-seed">Seed family</label><select id="atlas-seed"><option value="lame" selected>${escape(title)}</option></select></div><div><label for="atlas-decoration">Decoration</label><select id="atlas-decoration"><option value="attachment" selected>Hemisphere attachments (a,b,c)</option></select></div>`;
    }
    document.getElementById("atlas-seed").value = "lame";
    document.getElementById("atlas-decoration").value = "attachment";
    renderSummary();
    renderCaseGraph();
    placeAttachmentControls();
    document.documentElement.dataset.seed = "lame";
  }
  function boundaryPoint(edge, position) {
    return edge === 0 ? [.5 * position, .5]
      : edge === 1 ? [.5 * position, .5 * (1 - position)]
      : [.5, .5 * (1 - position)];
  }

  mount();
  updateLameGraph(lameGraphState.counts);
  return {
    isNew: false, active, project, map, svgPoint, svgInverse, renderMonodromy,
    asymptotic: () => null, boundaryPoint, updateLameGraph,
    // Compatibility with the common page-controller interface. This separate
    // atlas cannot commit a TV decoration or switch to an unavailable family.
    commitCase: () => false,
    previewGraph: () => {},
    switchFamily: family => family === "lame",
    graphRenderers: {lameGraph},
    drawDomain: (canvas, options) => window.AtlasView.drawDomain(canvas, {...options, markedStrip:true})
  };
})();

"use strict";
(() => {
  const get=id=>document.getElementById(id);
  const escape=value=>String(value).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
  const seedNames={bouquet:"Three-loop dual seed",double:"Double-edge dual seed",tetra:"Tetrahedral seed"};
  const doubleDefaults={qt:0,qu:0,ql:0,qm:0,modeA:"arm",countA:0,modeB:"arm",countB:0};
  const states={bouquet:{counts:{bouquetCore:"bouquet_B0_T0",bq0:0,bq1:0,bq2:0,bq3:0,bq4:0,bq5:0},chain:0},double:{counts:{...doubleDefaults},chain:0},tetra:{counts:Object.fromEntries(["qAB","qAC","qAD","qBC","qBD","qCD"].map(k=>[k,0])),chain:0}};
  const params=new URLSearchParams(location.search);
  let seed=["bouquet","double","tetra"].includes(params.get("seed"))?params.get("seed"):"double";
  let config,scene,selected=null,selectedChain=null,mode="edge",source=null,activeCycle="a",lengthDraft=false,arcRequest=null,arcResult=null;
  const developmentWidget=GraphDevelopment.create({onPole:i=>{setMode("pole");selectObject(`p${i}`);},onArc:request=>selectArc(request),onViewChange:()=>syncEdgeArcLinks()});
  const arcView=GraphArcView.create({onClear:()=>clearSelection(),onEdge:id=>selectObject(id),onArc:request=>selectArc(request)});
  const shiftWidget=GraphShift.create({onSource:p=>{lengthDraft=false;source=projectSource(...p);drawCharacters();updateExploreLink();}});
  const metricWidget=GraphMetric.create({onEdge:id=>{setMode("edge");selectObject(id);},onHemisphere:id=>{setMode("hemisphere");selectObject(id);},onCycle:cycle=>{activeCycle=cycle;arcRequest=null;arcResult=null;updateHolonomy();}});
  const maxHemispheres=16;
  const imported=GraphPageState.read(params,seed,states[seed].counts);
  states[seed].counts=imported.counts;
  if(imported.notice){
    const note=get("graph-import-note");note.hidden=false;
    note.append(document.createTextNode(imported.notice+" "));
    const link=document.createElement("a");link.href=new URL("../explore.html"+location.search,location.href).href;
    link.textContent="Open the requested full case in Explore →";note.append(link);
  }
  if(imported.source)source=projectSource(...imported.source);
  function makeConfig(){return seed==="bouquet"?BouquetModel.make("repeated",states[seed].counts):DecorationModel.make(seed,"repeated",states[seed].counts);}
  function countKey(index){
    if(seed==="bouquet")return `bq${index}`;
    const e=config.graph.stable_edges[index];
    if(seed==="tetra")return {"0,1":"qAB","0,2":"qAC","1,2":"qAD","0,3":"qBC","1,3":"qBD","2,3":"qCD"}[e.faces.slice().sort((a,b)=>a-b).join(",")];
    const darts=new Set(e.chain);
    return darts.has("t")||darts.has("T")?"qt":darts.has("u")||darts.has("U")?"qu":darts.has("l+")||darts.has("l-")?"ql":"qm";
  }
  function uiChainName(index){
    if(seed==="bouquet")return `ℓ${index+1}`;
    const key=countKey(index);
    return seed==="tetra"?key.slice(1):{qt:"t",qu:"u",ql:"ℓ",qm:"m"}[key];
  }
  function syncControls(){
    const container=get("decoration-controls");
    if(container.dataset.seed!==seed){
      container.innerHTML=config.graph.stable_edges.map((e,i)=>`<div class="decoration-control" data-control-chain="${i}"><label for="chain-count-${i}">${escape(uiChainName(i))}<small>p${e.faces[0]} ↔ p${e.faces[1]}</small></label><input id="chain-count-${i}" data-count-chain="${i}" type="number" min="0" max="3" step="1" aria-label="Pairs on ${escape(uiChainName(i))}"></div>`).join("");
      container.dataset.seed=seed;
    }
    config.graph.stable_edges.forEach((e,i)=>{const input=get(`chain-count-${i}`),q=states[seed].counts[countKey(i)];if(document.activeElement!==input)input.value=String(q);input.max=String(Math.max(0,Math.min(3,q+Math.floor((maxHemispheres-config.hemispheres)/2))));});
  }
  function updateMetric(){
    if(!config||!scene||!source)return;
    metricWidget.update({config,source,scene,selected,mode,chainIndex:selectedChain});
    const data=GraphMetric.metric(seed,...source),length=d=>config.graph.length_coefficients[d].reduce((s,n,i)=>s+n*data.normalized[i],0);
    get("length-z-label").textContent=`${data.names[2]} / π`;
    if(!lengthDraft)for(const [i,id]of ["length-x","length-y","length-z"].entries())if(document.activeElement!==get(id))get(id).value=String(Number(data.normalized[i].toFixed(8)));
    scene.metadata.edges.forEach(e=>document.querySelectorAll(`.graph-view [data-edge="${e.id}"]`).forEach(el=>el.setAttribute("aria-description",`Spherical length ${length(e.darts[0]).toFixed(4)}π`)));
    config.graph.stable_edges.forEach((e,i)=>{const out=document.querySelector(`[data-control-chain="${i}"] small`);if(out)out.textContent=`p${e.faces[0]} ↔ p${e.faces[1]} · first ${length(e.start).toFixed(3)}π`;});
    const out=get("edge-length-readout"),edge=scene.metadata.edges.find(e=>e.id===selected);
    if(out&&edge)out.textContent=`Actual spherical length: ${length(edge.darts[0]).toFixed(5)}π. Complementary segments along its chain have length ${(2-length(edge.darts[0])).toFixed(5)}π.`;
    updateHolonomy();
  }
  let holonomyKey="",holonomy;
  function updateHolonomy(){
    const key=`${config.id}|${source.join(",")}`;
    if(key!==holonomyKey){holonomy=GraphHolonomy.compute(config,source);holonomyKey=key;}
    const f=n=>Number(n.toFixed(6)).toString(),h=holonomy;
    get("holonomy-check").textContent=h.error<1e-9?"Numerical check: graph products agree with the character formulas to phase error below 10⁻⁹.":"The graph products and character formulas disagree; this case needs checking.";
    const edge=scene.metadata.edges.find(e=>e.id===selected)||scene.metadata.edges.find(e=>e.chainIndex===selectedChain)||scene.metadata.edges[0],t=h.transitions.find(t=>edge.darts.includes(t.d));
    const math=(body,block=true)=>`<math xmlns="http://www.w3.org/1998/Math/MathML"${block?' display="block"':''}>${body}</math>`;
    const index=(symbol,i)=>`<msub><mi>${symbol}</mi><mn>${i}</mn></msub>`;
    const word=key=>h.selectedwords[key].terms.map(term=>{const pair=`<mrow><mo stretchy="false">(</mo>${index("γ",term.i)}${index("γ",3)}<mo stretchy="false">)</mo></mrow>`;return term.power===1?pair:`<msup>${pair}<mn>${term.power}</mn></msup>`;}).join("")||'<mi>I</mi>';
    const datum=(symbol,subscript,value)=>`<msub><mi>${symbol}</mi>${subscript}</msub><mo>≈</mo><mn>${f(value)}</mn><mi>π</mi>`;
    get("gluing-detail").innerHTML=`<p><strong>${escape(edge.label)}: H<sub>${t.from+1}</sub> ↔ H<sub>${t.to+1}</sub>.</strong> For side length ℓ<sub>d</sub> and starting boundary angles α<sub>d</sub>, α<sub>d̄</sub>, reversing the edge gives</p>${math('<msub><mi>J</mi><mi>d</mi></msub><mo stretchy="false">(</mo><mi>ζ</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><msup><mi>e</mi><mrow><mi>i</mi><mo stretchy="false">(</mo><msub><mi>α</mi><mi>d</mi></msub><mo>+</mo><msub><mi>α</mi><mover><mi>d</mi><mo>¯</mo></mover></msub><mo>+</mo><msub><mi>ℓ</mi><mi>d</mi></msub><mo stretchy="false">)</mo></mrow></msup><mi>ζ</mi></mfrac>')}`;
    get("period-word-detail").innerHTML=`<p>The selected ordinary monodromy words for the marked periods are</p>${math(`<msub><mi>M</mi><mn>1</mn></msub><mo>=</mo>${word("one")}<mo>,</mo>`)}${math(`<msub><mi>M</mi><mi>τ</mi></msub><mo>=</mo>${word("tau")}<mo>.</mo>`)}`;
    const c=z=>math(`<mn>${f(z.re).replace("-","−")}</mn><mo>${z.im<0?"−":"+"}</mo><mn>${f(Math.abs(z.im))}</mn><mi>i</mi>`,false);
    const scientific=x=>{if(!x)return math('<mn>0</mn>',false);const [coefficient,exponent]=x.toExponential(1).split("e");return math(`<mn>${coefficient}</mn><mo>×</mo><msup><mn>10</mn><mn>${String(Number(exponent)).replace("-","−")}</mn></msup>`,false);};
    get("holonomy-table").innerHTML=`<caption>Ordinary character values from the graph products</caption><thead><tr><th>Period</th><th>χ, numerical value</th><th>Phase discrepancy (rad)</th></tr></thead><tbody><tr><td>1</td><td>${c(h.multipliers.one)}</td><td>${scientific(h.error1)}</td></tr><tr><td>τ</td><td>${c(h.multipliers.tau)}</td><td>${scientific(h.errorTau)}</td></tr></tbody>`;
    arcResult=arcRequest?GraphArcEdges.resolve({config,holonomy:h,scene,request:arcRequest}):null;
    developmentWidget.update({config,source,holonomy:h,cycle:activeCycle,arcRequest});
    shiftWidget.update({config,source,chainIndex:selectedChain});
    syncShiftControls();
    arcView.update({result:arcResult,scene,config,request:arcRequest});
    syncEdgeArcLinks();
  }
  function syncEdgeArcLinks(){
    if(!scene||!holonomy)return;
    const edge=mode==="edge"&&selected?scene.metadata.edges.find(item=>item.id===selected):null;
    const matches=[];
    if(edge)for(const link of developmentWidget.visibleArcRequests()){
      const result=GraphArcEdges.resolve({config,holonomy,scene,request:link.request});
      const term=result.terms.find(item=>item.edgeId===edge.id);
      if(term)matches.push({...link,result,coefficient:term.coefficient});
    }
    developmentWidget.highlightEdgeMatches(matches,edge?.label);
    const oldNote=get("selection-detail").querySelector(".edge-to-arcs-note");
    if(!edge){oldNote?.remove();return;}
    arcView.update({scene,config,edgeLink:{edge,matches}});
    const note=oldNote??document.createElement("span");note.className="edge-to-arcs-note";
    note.innerHTML=matches.length?`${matches.length} matching arc sum${matches.length===1?"":"s"} for ${escape(edge.label)}. <a href="#arc-selection-status">View highlighted arcs ↓</a>`:`No displayed arc sum uses ${escape(edge.label)} in this continuation. <a href="#arc-selection-status">Inspect another period or step ↓</a>`;
    if(!oldNote)get("selection-detail").append(note);
  }
  function syncShiftControls(){
    const choice=get("shift-chain-choice"),signature=`${seed}|${config.graph.stable_edges.map(e=>e.id).join(",")}`;
    if(choice.dataset.signature!==signature){choice.innerHTML=config.graph.stable_edges.map((e,i)=>`<option value="${i}">${escape(uiChainName(i))} · F${e.faces[0]} / F${e.faces[1]}</option>`).join("");choice.dataset.signature=signature;}
    const i=selectedChain??states[seed].chain,q=states[seed].counts[countKey(i)];choice.value=String(i);
    get("shift-chain-count").textContent=`${q} pair${q===1?"":"s"}`;
    get("shift-add-pair").disabled=q>=3||config.hemispheres+2>maxHemispheres;
    get("shift-remove-pair").disabled=q===0;
    get("shift-case-label").textContent=`${seedNames[seed]} · weights (${config.weights.join(", ")}) · ${config.hemispheres} hemispheres`;
  }
  get("shift-chain-choice").addEventListener("change",()=>{states[seed].chain=Number(get("shift-chain-choice").value);setMode("edge");selectedChain=states[seed].chain;refreshSelection();});
  for(const [id,step]of [["shift-add-pair",1],["shift-remove-pair",-1]])get(id).addEventListener("click",()=>{const top=get("shift-title").getBoundingClientRect().top,i=Number(get("shift-chain-choice").value),key=countKey(i);states[seed].chain=i;states[seed].counts[key]+=step;selectedChain=i;render();window.scrollBy(0,get("shift-title").getBoundingClientRect().top-top);});
  function render(){
    try{
      config=makeConfig();scene=GraphGeometry.build(config);
      const m=scene.metadata,N=config.hemispheres;
      if(m.edges.length!==N+2||m.poles.length!==4||m.hemispheres.length!==N)throw new Error("Graph cell counts do not match the hemisphere model.");
      get("original-graph").innerHTML=scene.primalSvg;get("dual-graph").innerHTML=scene.dualSvg;
      document.querySelectorAll(".graph-view .graph-edge").forEach(el=>{const chain=m.chains.find(c=>c.id===el.dataset.chain),hit=el.querySelector(".edge-hit");if(hit)hit.style.strokeWidth=chain?.multiplicity<=3?"15":"7";});
      get("original-count").textContent=`${N} hemisphere${N===1?"":"s"}`;
      get("dual-count").textContent=`${N} dot${N===1?"":"s"}`;
      get("case-description").textContent=`${seedNames[seed]} · weights (${config.weights.join(", ")}) · ${N} downstairs hemisphere${N===1?"":"s"}.`;
      states[seed].chain=Math.min(states[seed].chain,config.graph.stable_edges.length-1);syncControls();
      get("count-summary").innerHTML=`<table><thead><tr><th>Sphere decomposition</th><th>Vertices</th><th>Edges</th><th>Faces</th></tr></thead><tbody><tr><td>Original graph</td><td>4</td><td>${N+2}</td><td>${N}</td></tr><tr><td>Dual graph</td><td>${N}</td><td>${N+2}</td><td>4</td></tr></tbody></table><p>Both satisfy vertices − edges + faces = 2. Pair insertion adds 2 edges and 2 hemisphere regions.</p>`;
      selected=null;setMode(mode,true);drawCharacters();updateExploreLink();
      get("setup-status").textContent="";
    }catch(error){get("setup-status").textContent=`This graph could not be displayed: ${error.message}`;console.error(error);}
  }
  const modeAttribute={edge:"data-edge",hemisphere:"data-hemisphere",pole:"data-pole"};
  function setMode(next,preserveChain=false,preserveArc=false){
    if(!preserveArc){arcRequest=null;arcResult=null;}
    mode=next;selected=null;if(!preserveChain)selectedChain=null;get("graph-row").dataset.mode=mode;
    document.querySelectorAll('input[name="object-mode"]').forEach(el=>el.checked=el.value===mode);
    document.querySelectorAll("[data-mode-choice]").forEach(el=>el.dataset.active=String(el.dataset.modeChoice===mode));
    get("interaction-hint").textContent={edge:"Click an edge to link its decoration, its matching edge and its red arc sums. Click elsewhere to clear.",hemisphere:"Select a region H or its dual vertex H. Its spherical boundary is shown in Step 2.",pole:"Select pᵢ or face Fᵢ to read its cone angle."}[mode];
    const active=modeAttribute[mode];
    document.querySelectorAll("#graph-row .graph-view [data-edge],#graph-row .graph-view [data-hemisphere],#graph-row .graph-view [data-pole]").forEach(el=>{
      const enabled=el.hasAttribute(active);
      el.dataset.interactive=String(enabled);el.setAttribute("tabindex",enabled?"0":"-1");
      if(enabled){el.setAttribute("role","button");el.setAttribute("aria-pressed","false");const id=el.getAttribute(active);el.setAttribute("aria-label",objectLabel(mode,id));}
      else{el.removeAttribute("role");el.removeAttribute("aria-pressed");}
    });
    refreshSelection();
  }
  function objectLabel(kind,id){
    const list=kind==="edge"?scene.metadata.edges:kind==="pole"?scene.metadata.poles:scene.metadata.hemispheres;
    const item=list.find(v=>v.id===id);
    return kind==="edge"?`Select ${item?.label??id} and its corresponding edge`:kind==="pole"?`Select marked point ${item?.label??id} and its dual face`:`Select hemisphere ${item?.label??id} and its dual vertex`;
  }
  function names(items,ids,prefix){return ids.map(id=>items.find(v=>v.id===id||v.index===id)?.label??`${prefix}${id}`).join(" and ");}
  function selectObject(id){
    arcRequest=null;arcResult=null;
    selected=id;
    if(mode==="edge"){
      const e=scene.metadata.edges.find(v=>v.id===id);
      if(e){states[seed].chain=e.chainIndex;selectedChain=e.chainIndex;syncControls();}
    }
    refreshSelection();
  }
  function refreshSelection(){
    const attr=modeAttribute[mode];
    document.querySelectorAll(".graph-view .is-selected").forEach(el=>el.classList.remove("is-selected"));
    document.querySelectorAll(".is-chain-selected").forEach(el=>el.classList.remove("is-chain-selected"));
    document.querySelectorAll(".decoration-control").forEach(el=>el.classList.toggle("is-selected",Number(el.dataset.controlChain)===selectedChain));
    if(selectedChain!==null)document.querySelectorAll(`.graph-view [data-chain="chain${selectedChain}"]`).forEach(el=>el.classList.add("is-chain-selected"));
    document.querySelectorAll(".graph-view [aria-pressed]").forEach(el=>el.setAttribute("aria-pressed","false"));
    if(selected)document.querySelectorAll(`.graph-view [${attr}="${selected}"]`).forEach(el=>{el.classList.add("is-selected");el.setAttribute("aria-pressed","true");});
    updateMetric();
    if(arcResult){get("selection-detail").innerHTML=`<strong>${escape(arcResult.title)} · ${arcResult.terms.length} contributing edges.</strong><span>Solid red edges add their length; dashed red edges subtract. Matching edge labels appear in both drawings.</span>`;return;}
    if(!selected){get("selection-detail").innerHTML=selectedChain!==null?`<strong>Chain ${escape(uiChainName(selectedChain))} · ${config.graph.stable_edges[selectedChain].count} pair insertions.</strong><span>Both drawings show its ${config.graph.stable_edges[selectedChain].multiplicity} edges. Each pair creates two new hemisphere regions.</span>`:`<strong>${{edge:"Select an edge or a decoration count.",hemisphere:"Select a hemisphere region or a dual vertex.",pole:"Select a marked vertex or a labelled dual face."}[mode]}</strong><span>Click elsewhere to clear the highlight. Step 2 assigns actual spherical lengths to these edges.</span>`;return;}
    const m=scene.metadata;
    if(mode==="edge"){
      const e=m.edges.find(v=>v.id===selected);if(!e)return;
      const ps=names(m.poles,e.poles,"p"),hs=names(m.hemispheres,e.hemispheres,"H");
      const same=e.hemispheres[0]===e.hemispheres[1];
      get("selection-detail").innerHTML=`<strong>${escape(e.label)} · chain ${escape(uiChainName(e.chainIndex))}: original edge ↔ dual edge.</strong><span>Endpoints ${escape(ps)}; adjacent hemispheres ${escape(hs)}.${same?" A bridge in the original graph becomes a dual loop.":""} The corresponding count is highlighted on the left.</span><span id="edge-length-readout"></span>`;
    }else if(mode==="hemisphere"){
      const h=m.hemispheres.find(v=>v.id===selected);if(!h)return;
      get("selection-detail").innerHTML=`<strong>${escape(h.label)}: one hemisphere region ↔ one dual vertex.</strong><span>Its boundary has ${h.degree} edge incidences, counted with repetition.${h.exterior?" This is the unbounded region of the planar drawing; on the sphere it is an ordinary hemisphere face.":""} The dual vertex has the same degree.</span>`;
    }else{
      const p=m.poles.find(v=>v.id===selected);if(!p)return;
      get("selection-detail").innerHTML=`<strong>${escape(p.label)} ↔ dual face F${p.index}.</strong><span>Degree ${p.degree} = 2 × ${config.weights[p.index]} + 1. The upstairs cone angle is ${2*p.degree}π.${config.weights[p.index]===0?" This is a smooth marked point on the torus.":""} Both sides of a bridge count in this boundary.</span>`;
    }
    updateMetric();
  }
  function selectArc(request){arcRequest=JSON.stringify(arcRequest)===JSON.stringify(request)?null:request;arcResult=null;selected=null;selectedChain=null;setMode("edge",false,true);}
  function clearSelection(){selected=null;selectedChain=null;arcRequest=null;arcResult=null;refreshSelection();}
  get("graph-row").addEventListener("click",event=>{const node=event.target.closest(`[${modeAttribute[mode]}]`);if(node){event.preventDefault();selectObject(node.getAttribute(modeAttribute[mode]));}else clearSelection();});
  get("graph-row").addEventListener("keydown",event=>{if(event.key==="Escape"){clearSelection();return;}if(!["Enter"," "].includes(event.key))return;const node=event.target.closest(`[${modeAttribute[mode]}]`);if(node){event.preventDefault();selectObject(node.getAttribute(modeAttribute[mode]));}});
  document.addEventListener("click",event=>{const controls="#graph-row,.comparison-toolbar,.dictionary-grid,.reading-jumps,#selection-detail,.decoration-controls,.setup input,.setup select,.setup button,.setup label,.metric-lab svg,.metric-lab input,.metric-lab label,.metric-lab select,.metric-lab button,.length-entry,.development-section svg,.development-section input,.development-section button,#halfturn-guide,#halfturn-steps,#arc-selection-status,#arc-linked-graphs,.shift-section,.calculation-details";if(event.composedPath().some(node=>node instanceof Element&&node.matches(controls)))return;clearSelection();});
  document.querySelectorAll('input[name="object-mode"]').forEach(el=>el.addEventListener("change",()=>setMode(el.value)));
  document.querySelectorAll("[data-mode-choice]").forEach(el=>el.addEventListener("click",()=>setMode(el.dataset.modeChoice)));
  get("clear-selection").addEventListener("click",clearSelection);
  document.addEventListener("keydown",event=>{if(event.key==="Escape"&&(arcRequest||selected||selectedChain!==null))clearSelection();});
  get("larger-drawings").addEventListener("change",()=>get("graph-row").classList.toggle("is-enlarged",get("larger-drawings").checked));
  const focusChain=event=>{const input=event.target.closest("[data-count-chain]");if(!input)return;states[seed].chain=Number(input.dataset.countChain);setMode("edge");selectedChain=states[seed].chain;refreshSelection();};
  get("decoration-controls").addEventListener("focusin",focusChain);
  get("decoration-controls").addEventListener("click",focusChain);
  get("decoration-controls").addEventListener("input",event=>{const input=event.target.closest("[data-count-chain]");if(!input)return;const index=Number(input.dataset.countChain),value=Number(input.value);if(input.value===""||!Number.isInteger(value)||value<0||value>Number(input.max)){get("setup-status").textContent=`Choose a whole number from 0 to ${input.max}.`;return;}states[seed].counts[countKey(index)]=value;selectedChain=index;render();});
  get("graph-seed").value=seed;
  get("graph-seed").addEventListener("change",()=>{seed=get("graph-seed").value;source=null;selectedChain=null;lengthDraft=false;get("length-status").textContent="";render();});
  get("reset-graph").addEventListener("click",()=>{for(const key of Object.keys(states[seed].counts))if(!["bouquetCore","modeA","modeB"].includes(key))states[seed].counts[key]=0;render();});

  function triangleData(){const x=seed==="bouquet"?.5:1,y=seed==="tetra"?1:.5;return {x,y,vertices:[[x,0],[0,y],[x,y]]};}
  function projectSource(R,S){
    const {x,y,vertices}=triangleData(),inset=.06,vs=vertices.map(v=>[(1-3*inset)*v[0]+inset*2*x,(1-3*inset)*v[1]+inset*2*y]);
    if(R<=x*(1-inset)&&S<=y*(1-inset)&&R/x+S/y>=1+inset)return [R,S];
    let best=null,d=Infinity;
    vs.forEach((a,i)=>{const b=vs[(i+1)%3],dx=b[0]-a[0],dy=b[1]-a[1],t=Math.max(0,Math.min(1,((R-a[0])*dx+(S-a[1])*dy)/(dx*dx+dy*dy))),p=[a[0]+t*dx,a[1]+t*dy],dist=(p[0]-R)**2+(p[1]-S)**2;if(dist<d){best=p;d=dist;}});
    return best;
  }
  let characterScale=1;
  function drawCharacters(){
    const t=triangleData();source??={bouquet:[.31,.34],double:[.62,.29],tetra:[.58,.62]}[seed].slice();
    get("seed-character").setAttribute("viewBox",seed==="double"?"0 85 320 190":"0 0 320 280");
    const limitX=t.x,limitY=t.y;characterScale=Math.min(235/(limitX+.08),190/(limitY+.06));
    const number=n=>Number.isInteger(n)?String(n):n===.5?"½":n===1.5?"1½":String(n);
    for(const id of ["seed-character"]){
      const color="#315d83",xy=p=>[42+characterScale*p[0],228-characterScale*p[1]],inside=t.vertices.map(p=>[(1-.18)*p[0]+.12*t.x,(1-.18)*p[1]+.12*t.y]),dot=xy(source);
      get(id).innerHTML=`<path d="M42 ${228-characterScale*limitY-9} V228 H${42+characterScale*limitX+10}" fill="none" stroke="#9cae9d" stroke-width="1.2"/>${Array.from({length:Math.floor(limitX*2)+1},(_,i)=>i/2).map(n=>`<text x="${42+characterScale*n}" y="249" text-anchor="middle" font-size="12">${number(n)}</text>`).join("")}${Array.from({length:Math.floor(limitY*2)},(_,i)=>(i+1)/2).map(n=>`<text x="33" y="${232-characterScale*n}" text-anchor="end" font-size="12">${number(n)}</text>`).join("")}<text x="${53+characterScale*limitX}" y="232" font-size="13">R</text><text x="33" y="${214-characterScale*limitY}" font-size="13">S</text><polygon points="${t.vertices.map(p=>xy(p).join(",")).join(" ")}" fill="${color}08" stroke="${color}" stroke-width="1.4" stroke-dasharray="5 4"/><polygon points="${inside.map(p=>xy(p).join(",")).join(" ")}" fill="${color}15" stroke="${color}" stroke-width="1.2"/><circle cx="${dot[0]}" cy="${dot[1]}" r="11" fill="${color}25"/><circle class="character-dot" cx="${dot[0]}" cy="${dot[1]}" r="5.5" fill="${color}"/>`;
    }
    get("source-inequality").textContent=seed==="bouquet"?"Source C: 0 < R, S < ½ and R + S > ½.":seed==="double"?"Source D: 0 < R < 1, 0 < S < ½ and R + 2S > 1.":"Source T: 0 < R, S < 1 and R + S > 1.";
    for(const [i,id]of ["source-r","source-s"].entries())if(document.activeElement!==get(id))get(id).value=String(Number(source[i].toFixed(6)));
    updateMetric();
  }
  for(const id of ["seed-character"]){
    const svg=get(id);let dragging=false;
    const move=event=>{const p=svg.createSVGPoint();p.x=event.clientX;p.y=event.clientY;const q=p.matrixTransform(svg.getScreenCTM().inverse());lengthDraft=false;source=projectSource((q.x-42)/characterScale,(228-q.y)/characterScale);drawCharacters();updateExploreLink();};
    svg.addEventListener("pointerdown",event=>{if(event.button!==0)return;event.preventDefault();dragging=true;svg.setPointerCapture(event.pointerId);move(event);});
    svg.addEventListener("pointermove",event=>{if(dragging)move(event);});
    svg.addEventListener("pointerup",()=>dragging=false);svg.addEventListener("pointercancel",()=>dragging=false);
    svg.addEventListener("keydown",event=>{const directions={ArrowLeft:[-1,0],ArrowRight:[1,0],ArrowUp:[0,1],ArrowDown:[0,-1]},d=directions[event.key];if(!d)return;event.preventDefault();const t=triangleData(),step=event.shiftKey?.04:.01;source=projectSource(source[0]+d[0]*t.x*step,source[1]+d[1]*t.y*step);drawCharacters();updateExploreLink();});
  }
  function editCoordinates(event){
    if(get("source-r").value===""||get("source-s").value==="")return;
    const R=Number(get("source-r").value),S=Number(get("source-s").value);if(!Number.isFinite(R)||!Number.isFinite(S))return;
    const projected=projectSource(R,S);
    if(event.type==="input"&&(Math.abs(R-projected[0])>1e-10||Math.abs(S-projected[1])>1e-10))return;
    lengthDraft=false;source=projected;drawCharacters();updateExploreLink();
    if(event.type==="change")for(const [i,id]of ["source-r","source-s"].entries())get(id).value=String(Number(source[i].toFixed(6)));
  }
  for(const id of ["source-r","source-s"])for(const type of ["input","change"])get(id).addEventListener(type,editCoordinates);
  for(const id of ["length-x","length-y","length-z"])get(id).addEventListener("input",()=>{lengthDraft=true;get("length-status").textContent="Apply the lengths to update all diagrams.";});
  get("apply-lengths").addEventListener("click",()=>{
    const values=["length-x","length-y","length-z"].map(id=>Number(get(id).value));
    if(values.some(n=>!Number.isFinite(n)||n<=0)){get("length-status").textContent="All three lengths must be strictly positive.";return;}
    const [x,y,z]=values,perimeter=seed==="bouquet"?2*(x+y+z):seed==="double"?x+y+2*z:x+y+z;
    if(Math.abs(perimeter-2)>1e-7){get("length-status").textContent=`The hemisphere perimeter would be ${perimeter.toFixed(6)}π. It must equal 2π; the input has not been rescaled.`;return;}
    const R=seed==="double"?(y+2*z)/2:(y+z)/2,S=seed==="double"?(x+y)/4:(x+z)/2;
    const triangle=triangleData();if(!(R>0&&S>0&&R<triangle.x&&S<triangle.y&&R/triangle.x+S/triangle.y>1)){get("length-status").textContent="These lengths do not define a positive metric in this seed family.";return;}
    lengthDraft=false;source=[R,S];drawCharacters();updateExploreLink();get("length-status").textContent="Valid spherical metric. Its developed points, both characters and triangle shift are updated below.";
  });
  function updateExploreLink(){const url=new URL("../explore.html",location.href);url.searchParams.set("seed",seed);url.searchParams.set("decoration","repeated");for(const [key,value]of Object.entries(states[seed].counts))url.searchParams.set(key,String(value));url.searchParams.set("r",String(source[0]));url.searchParams.set("s",String(source[1]));url.searchParams.set("autocompute","0");get("open-explore").href=url.href;}
  if(!window.GraphGeometry){get("setup-status").textContent="The graph drawing module is missing.";return;}
  render();
})();

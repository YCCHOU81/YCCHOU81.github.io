"use strict";

/* TV seed/decorations controls and labelled metric graphs.
 * All three families share one case-description path and the operation registry
 * in AtlasFamily. R,S are geometric characters; the analytic evaluator receives
 * (r,s)=(S,R). No modular parking or inversion is applied here.
 * Family events and public graph helpers retain their existing protocols.
 */
window.SeedCases = (() => {
  const query = new URLSearchParams(location.search);
  const allowedSeeds = window.AtlasFamily.allowedSeeds;
  let seed = allowedSeeds.includes(query.get("seed")) ? query.get("seed") : window.AtlasFamily.defaultSeed;
  let operations = window.AtlasFamily.operations[seed];
  let operation = operations.includes(query.get("decoration")) ? query.get("decoration") : operations[0];
  const names = {bouquet:"Three-loop seed · all admissible decorations",double:"Double-edge seed · (1,1,0,0)",tetra:"Tetrahedral seed · (1,1,1,1)"};
  const opNames = {seed:"Seed only",pair:"Insert a pair on F₀ / F₁",split:"Split the F₂ loop-dart block",sprout:"Sprout the F₂ monogon",repeated:"Repeated decorations · choose counts"};
  let caseConfig=seed==="bouquet"?BouquetModel.make(operation,BouquetModel.parse(operation,query)):DecorationModel.make(seed,operation,DecorationModel.parse(seed,operation,query));
  let shift = caseConfig?.shift??[0,0];
  let top = seed === "tetra" ? 1 : .5;
  let maxR = seed==="bouquet" ? .5 : 1;
  const active = caseMetadata(caseConfig);

  // One description feeds initial render, decoration changes, and family restore.
  // Model-specific construction stays in DecorationModel and BouquetModel.
  function caseMetadata(config) {
    const familyName = seed === "double" ? "Double-edge" : seed === "bouquet" ? "Three-loop" : "Tetrahedron";
    const decorationName = seed === "bouquet" && operation === "repeated"
      ? bouquetLabel(BouquetModel.cores[config.core]) : opNames[operation];
    const expected = config && !config.connector && operation === "repeated" ? "expected " : "";
    const edgeLabels = seed === "bouquet"
      ? ["R → ½ · cusp ∞", "S → ½ · cusp 0", "R + S → ½ · cusp 1"]
      : [
          `R → 1 · ${expected}cusp ∞`,
          `S → ${top === 1 ? "1" : "½"} · ${expected}cusp 0`,
          `${top === 1 ? "R + S" : "R + 2S"} → 1 · ${expected}cusp ${top === 1 ? "1" : "½"}`
        ];
    return {
      id:config?.id ?? `${seed}-${operation}`, seed, operation,
      weights:config?.weights ?? [1,0,0,0], shift, top, maxR, caseConfig:config,
      isRepeated:operation === "repeated" || seed === "bouquet",
      title:`${familyName} · ${decorationName} · (${config.weights.join(",")})`,
      baseTitle:seed === "bouquet" ? "Three-loop seed · (1,0,0,0)" : names[seed],
      baseId:`${seed}-seed`,
      start:seed === "bouquet" ? [.29,.31] : config?.start ?? (seed === "tetra" ? [.58,.62] : [.58,.31]),
      sourceVertices:config?.sourceVertices ?? [[maxR,0],[0,top],[maxR,top]],
      edgeLabels
    };
  }

  const escape = value => String(value).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
  const isLimitsPage=()=>Boolean(document.getElementById("limits-controls"));
  const applyButtonLabel=()=>isLimitsPage()?"Apply decorations":"Apply &amp; compute";
  const applyInstruction=()=>isLimitsPage()?"Apply decorations updates the case. Then choose the source path below and use Compute approach.":"Apply &amp; compute updates the case and starts both recommended previews.";
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
  let scale=Math.min(300/(maxR+shift[0]),300/(seed==="bouquet"?top+shift[1]:Math.max(1,top+shift[1])));
  function svgPoint(R,S,attached=false) {const p=map(R,S,attached);return [48+scale*p.r,350-scale*p.s];}
  function svgInverse(x,y,attached=false) {return {r:(x-48)/scale-(attached?shift[0]:0),s:(350-y)/scale-(attached?shift[1]:0)};}
  function renderMonodromy(svg,{point,attached=false,title="Monodromy triangle"}={}) {
    const color=attached?"#b55e3c":"#285f65", xy=p=>svgPoint(...p,attached).join(",");
    const xMax=maxR+shift[0],yMax=seed==="bouquet"?top+shift[1]:Math.max(1,top+shift[1]);
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
  const labels={"1/4":.25,"1/8":.125,"1/6":1/6,"1/2":.5,"-1/4":-.25,"-3/8":-.375,"-1/8":-.125,"-1/6":-1/6,"2/3":2/3,"1/3":1/3,"0":0};
  const tables={
    "double-seed":[["1/4",.5,4],["0",2,2],["-1/4",.5,4]],
    "double-pair":[["1/4",.5,12],["1/2",1,12],["-3/8",.25,24]],
    "double-split":[["1/8",.25,40],["0",2,2*Math.sqrt(6)],["-1/8",.25,40]],
    "double-sprout":[["1/6",1/3,12],["0",2,2*Math.sqrt(3)],["-1/6",1/3,12]],
    "tetra-seed":[["1/2",1,4],["1/2",1,4],["1/2",1,4]],
    "tetra-pair":[["1/2",1,4*Math.sqrt(3)],["2/3",2/3,12],["1/3",2/3,12]]
  };
  function rational(a,b=1){const gcd=(x,y)=>y?gcd(y,x%y):Math.abs(x),d=gcd(a,b)||1;a/=d;b/=d;return {value:a/b,label:b===1?String(a):`${a}/${b}`,numerator:a,denominator:b};}
  function asymptotic(edge,R,S,base=false) {
    if (![0,1,2].includes(edge) || false) return null;
    if(seed==="bouquet")return {...BouquetModel.asymptotic(caseConfig,edge,R,S,base),proved:false,proposal:true,countFormula:!base?window.TVCuspFormulas?.lookup(caseConfig,edge):null,
      detail:"Graph-based asymptotic proposal for this labelled three-loop graph and fixed counts, away from source vertices. This simulation displays the proposed value for comparison; the numerical points are independently refined roots."};
    let row=tables[base?active.baseId:active.id]?.[edge],rat=null;
    const countFormula=!base?window.TVCuspFormulas?.lookup(caseConfig,edge):null;
    if((seed==="double"||seed==="tetra")&&countFormula){rat=countFormula;row=[rat.label,countFormula.alpha,countFormula.constantFactor];}
    if(!row&&!base&&caseConfig?.connector){
      const n=caseConfig.connector.q+1,tetra=seed==="tetra";
      let bin=1;for(let j=1;j<=n;j++)bin=bin*(n+j)/j;
      rat=tetra?[rational(1,2),rational(n,n+1),rational(1,n+1)][edge]:[rational(1,4),rational(n-1,n),rational(1-2*n,4*n)][edge];
      const exponent=tetra?[1,2/(n+1),2/(n+1)][edge]:[.5,2/n,1/(2*n)][edge];
      const coefficient=tetra?[Math.sqrt(8*n*(n+1)),2*bin,2*bin][edge]:[2*n*(n+1),n*bin,2*n*bin][edge];
      row=[rat.label,exponent,coefficient];
    }
    if(!row&&countFormula){rat=countFormula;row=[rat.label,countFormula.alpha,countFormula.constantFactor];}
    const chart=edge===0?"raw":edge===1?"zero":seed==="double"?"half":"one";
    const sinArgument=edge===0?(seed==="double"?2*Math.PI*S:Math.PI*S):Math.PI*R;
    const cusp=edge===0?"∞":edge===1?"0":seed==="double"?"½":"1";
    const chartLabel={raw:"τ",zero:"w = 1 − 1/τ",half:"w = τ/(1 − 2τ)",one:"w = 1/(1 − τ)"}[chart];
    if(!row)return {known:false,cusp,chart,chartLabel,x0:null,xLabel:"computed from the selected branch",alpha:null,C:null,
      deltaLabel:["1 − R",top===1?"1 − S":"½ − S",top===1?"R + S − 1":"R + 2S − 1"][edge],
      detail:`The selected chart targets cusp ${cusp}; the displayed tail is computed on this labelled decorated graph. A closed leading coefficient is not asserted for this mixed decoration.`};
    if(!rat){const parts=row[0].split("/").map(Number);rat=rational(parts[0],parts[1]??1);}
    const alphaParts=countFormula?.alphaLabel?.split("/").map(Number),logCoefficientLabel=alphaParts?(alphaParts.length===2?`${alphaParts[0]}/(${alphaParts[1]}π)`:`${alphaParts[0]}/π`):null;
    const logK=countFormula?.logConstantFactor??(row[2]>0?Math.log(row[2]):null),logC=Number.isFinite(logK)?logK+Math.log(Math.sin(sinArgument)):null;
    return {known:true,proved:false,proposal:true,cusp,chart,chartLabel,countFormula,M:countFormula?.M,p:countFormula?.p,alphaLabel:countFormula?.alphaLabel,logCoefficientLabel,KLabel:countFormula?.constantLabel,logK,logC,x0:rat.value,xLabel:rat.label,xNumerator:rat.numerator,xDenominator:rat.denominator,numerator:rat.numerator,denominator:rat.denominator,alpha:row[1],C:Number.isFinite(logC)?Math.exp(logC):null,
      deltaLabel:["1 − R",top===1?"1 − S":"½ − S",top===1?"R + S − 1":"R + 2S − 1"][edge],
      constantFactor:row[2],sineLabel:edge===0?(seed==="double"?"sin(2πS)":"sin(πS)"):"sin(πR)",
      detail:countFormula?`${countFormula.name}: ${countFormula.parameters} Proposed Re ${chart==="raw"?"τ":"w"} → ${rat.label}, with logarithmic slope α = ${row[1]}. ${Number.isFinite(logK)?"":"The leading constant C has not been assigned for this mixed family. "}${countFormula.evidence}`:"Graph-based asymptotic proposal on the selected geometric branch, away from source vertices; displayed for simulation, without an exact-limit proof claim."};
  }
  function graphSVG(family,op="seed") {
    const tint="#b55e3c",green="#285f65",stroke="#496351";
    const path=(d,changed=false)=>`<path d="${d}" fill="none" stroke="${changed?tint:stroke}" stroke-width="${changed?3:2.3}"/>`;
    const node=(x,y,changed=false)=>`<circle cx="${x}" cy="${y}" r="${changed?4:4.5}" fill="${changed?tint:green}" stroke="#fcfcf7" stroke-width="1"/>`;
    let paths="",nodes="",labels="";
    if(family==="double"){
      const detached=["split","sprout"].includes(op);
      paths+=path("M86 60 Q130 22 178 60",op==="pair")+path("M86 60 Q130 102 178 60");
      paths+=path(detached?"M42 60 C4 8 4 112 42 60":"M86 60 C30 0 30 120 86 60");
      // One monogon lies inside the connector bigon. Putting both outside
      // changes the face degrees and would depict a different ribbon graph.
      paths+=path("M178 60 C135 38 135 82 178 60");
      nodes+=node(86,60)+node(178,60);
      if(detached){paths+=path("M42 60 L86 60",true);nodes+=node(42,60);if(op==="split")nodes+=node(64,60,true);}
      if(op==="pair")nodes+=node(116,43,true)+node(148,43,true);
      labels=`<text x="130" y="14">F₀</text><text x="120" y="76">F₁</text><text x="${detached?21:54}" y="65">F₂</text><text x="157" y="64">F₃</text>`;
    }else if(family==="tetra"){
      paths=path("M42 98 L212 98 L127 10 Z")+path("M42 98 L127 64 L212 98")+path("M127 10 L127 64");
      if(op==="pair")paths+=path("M42 98 L212 98",true);
      nodes=node(42,98)+node(212,98)+node(127,10)+node(127,64);
      if(op==="pair")nodes+=node(98,98,true)+node(155,98,true);
      labels='<text x="127" y="86">F₁</text><text x="230" y="26">F₀</text><text x="94" y="57">F₂</text><text x="160" y="57">F₃</text>';
    }else{
      paths=path("M127 65 C50 -20 20 100 127 65")+path("M127 65 C210 -20 240 100 127 65")+path("M127 65 C40 136 210 136 127 65");nodes=node(127,65);
      labels='<text x="127" y="10">F₀</text><text x="72" y="46">F₁</text><text x="185" y="46">F₂</text><text x="127" y="102">F₃</text>';
    }
    const geometry=family==="double"?[
      {name:"t",a:[86,60],b:[178,60],c:[[130,22]]},{name:"u",a:[86,60],b:[178,60],c:[[130,102]]},
      {name:"ℓ",a:[86,60],b:[86,60],c:[[30,0],[30,120]]},{name:"m",a:[178,60],b:[178,60],c:[[135,38],[135,82]]}
    ]:family==="tetra"?[
      {name:"AB",a:[42,98],b:[212,98]},{name:"AC",a:[42,98],b:[127,10]},{name:"AD",a:[42,98],b:[127,64]},
      {name:"BC",a:[212,98],b:[127,10]},{name:"BD",a:[212,98],b:[127,64]},{name:"CD",a:[127,10],b:[127,64]}
    ]:[
      {name:"ℓ₁",a:[127,65],b:[127,65],c:[[50,-20],[20,100]]},
      {name:"ℓ₂",a:[127,65],b:[127,65],c:[[210,-20],[240,100]]},
      {name:"ℓ₃",a:[127,65],b:[127,65],c:[[40,136],[210,136]]}
    ];
    const occupied=[...labels.matchAll(/<text x="([\d.-]+)" y="([\d.-]+)"/g)].map(m=>[Number(m[1]),Number(m[2])-4,19,14]);
    const vertices=[...nodes.matchAll(/<circle cx="([\d.-]+)" cy="([\d.-]+)"/g)].map(m=>[Number(m[1]),Number(m[2])]);
    const countLabels=edgeCountLabels(geometry.map(e=>({name:e.name,m:1,points:sampleCurve(e.a,e.b,e.c),color:stroke})),{occupied,vertices});
    return `<svg viewBox="0 -5 255 130" role="img" aria-label="${escape(names[family])}, ${escape(opNames[op]??op)} dual graph; all base edges have count one"><g>${paths}${nodes}</g><g text-anchor="middle" class="case-face-labels">${labels}</g>${countLabels}</svg>`;
  }
  // Full chain names stay in the key. Short, colour-matched count badges are
  // placed beside each chain, away from the marked faces and core vertices.
  const graphPalette=["#2c6b70","#b15b31","#6972a8","#927019","#8c527e","#387a52"];
  // The marked result graph and the form share explicit control IDs. These
  // links only move focus/selection: they never edit a value or apply a case.
  let selectedGraphControl=null,graphLinkScope="";
  const linkTokens=value=>String(value??"").split(/\s+/).filter(Boolean);
  const controlIDs=(family,edge,config)=>{
    if(config.operation!=="repeated")return [];
    if(family==="tetra"){
      const key={"0,1":"qAB","0,2":"qAC","1,2":"qAD","0,3":"qBC","1,3":"qBD","2,3":"qCD"}[edge.faces.slice().sort((a,b)=>a-b).join(",")];
      return key?["repeat-"+key]:[];
    }
    const ds=new Set(edge.chain),has=(...names)=>names.some(name=>ds.has(name));
    const key=has("t","T")?"qt":has("u","U")?"qu":has("l+","l-")?"ql":has("m+","m-")?"qm":null;
    return key?["repeat-"+key]:["repeat-count"+edge.id[0],"repeat-mode"+edge.id[0]];
  };
  const linkAttributes=(ids,name,count,{keyboard=true,related=[]}={})=>ids?.length?` class="graph-control-link" data-graph-controls="${escape(ids.join(" "))}"${related.length?` data-graph-related="${escape(related.join(" "))}"`:""} data-edge="${escape(name)}" role="button" tabindex="${keyboard?0:-1}" aria-controls="${escape(ids.join(" "))}" aria-pressed="false" aria-label="${escape(name)}: ${count} edge${count===1?"":"s"}. Select its decoration control"`:"";
  function linkedPath(d,{name,m,color,controls=[],related=[],width=2.8,title=""}){
    const path=`<path class="graph-chain" d="${d}" fill="none" stroke="${color}" stroke-width="${width}"><title>${escape(title||`${name} ×${m}`)}</title></path>`;
    if(!controls.length)return path;
    return `<g${linkAttributes(controls,name,m,{related})}><path class="graph-chain-halo" d="${d}" fill="none" stroke="#e9bc63" stroke-width="10" opacity="0" pointer-events="none"/>${path}<path class="graph-chain-hit" d="${d}" fill="none" stroke="transparent" stroke-width="14" pointer-events="stroke"/></g>`;
  }
  function controlLabel(control){
    const label=document.querySelector(`label[for="${control.id}"]`);
    const text=(label?.textContent??control.name??control.id).replace(/\s+/g," ").trim();
    return control.id.startsWith("repeat-count")?`Count at vertex ${control.id.at(-1)}`:text;
  }
  function graphControlVisible(control){return Boolean(control&&!control.disabled&&!control.closest("[hidden]"));}
  function refreshGraphControlLinks(){
    const graph=document.getElementById("case-graph");if(!graph)return;
    const scope=`${seed}:${operation}:${seed==="bouquet"?caseConfig?.core:""}`;
    if(scope!==graphLinkScope){selectedGraphControl=null;graphLinkScope=scope;}
    const selected=document.getElementById(selectedGraphControl);
    if(!graphControlVisible(selected))selectedGraphControl=null;
    for(const element of document.querySelectorAll(".graph-control-selected"))element.classList.remove("graph-control-selected");
    for(const element of graph.querySelectorAll("[data-graph-controls]")){
      const match=Boolean(selectedGraphControl&&[...linkTokens(element.dataset.graphControls),...linkTokens(element.dataset.graphRelated)].includes(selectedGraphControl));
      element.classList.toggle("graph-link-selected",match);element.setAttribute("aria-pressed",String(match));
    }
    if(selectedGraphControl){
      const field=selected.closest(".repeat-local")??selected.closest(".repeat-count-field")??selected.parentElement;
      field?.classList.add("graph-control-selected");
    }
    const status=document.getElementById("graph-control-feedback");
    if(status){
      const note="";
      const message=selectedGraphControl?`Selected ${controlLabel(selected)}.${note} Apply decorations to update τ.`:(operation==="repeated")?"Click an edge or count to highlight its control.":"Fixed preset · choose Repeated decorations to edit chains.";
      if(status.textContent!==message)status.textContent=message;
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
  function linkedFormControl(target){
    const control=target?.closest?.("input, select");if(!control)return null;
    const id=control.id;
    
    if(operation!=="repeated")return null;
    return (seed==="bouquet"?/^bouquet-bq[0-5]$/:/^repeat-(?:q(?:t|u|l|m|AB|AC|AD|BC|BD|CD)|(?:mode|count)[AB])$/).test(id)?control:null;
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
  function faceText(polygons,preferred,callouts,occupied=[]){
    const safe=(p,poly)=>[[-7,-6],[7,-6],[-7,5],[7,5]].every(([x,y])=>BouquetLayout.inside([p[0]+x,p[1]+y],poly));
    return [0,1,2,3].map(i=>{
      if(i===0){occupied.push([callouts[0][0],callouts[0][1]-4,19,14]);return `<text x="${callouts[0][0]}" y="${callouts[0][1]}" text-anchor="middle">F₀</text>`;}
      const p=preferred[i]??BouquetLayout.interiorLabel(polygons[i]);
      if(safe(p,polygons[i])){occupied.push([p[0],p[1]-1,19,14]);return `<text x="${p[0]}" y="${p[1]+3}" text-anchor="middle">F${i}</text>`;}
      const q=callouts[i];
      occupied.push([q[0],q[1]-8,19,14]);
      return `<path class="graph-face-leader" d="M${q[0]} ${q[1]} L${p[0]} ${p[1]}"/><circle class="graph-face-target" cx="${p[0]}" cy="${p[1]}" r="1.8"/><text x="${q[0]}" y="${q[1]-4}" text-anchor="middle">F${i}</text>`;
    }).join("");
  }
  function repeatedGraph(config){
    const g=config.graph,verts=g.vertices.filter(v=>v.length!==2),positions=new Map(),family=config.seed,curves=new Map();
    const vertexPoint=v=>{
      if(family==="tetra")return {A:[42,108],B:[220,108],C:[131,10],D:[131,70]}[v[0][0]];
      if(v.includes("t")&&v.includes("u"))return[81,70];
      if(v.includes("T")&&v.includes("U"))return[192,70];
      if(v.includes("l+")&&v.includes("l-"))return[31,70];
      if(v.includes("m+")&&v.includes("m-"))return[156,70];
      if(v.includes("t"))return[81,49];if(v.includes("u"))return[81,91];
      if(v.includes("T"))return[192,49];if(v.includes("U"))return[192,91];
      return[131,70];
    };
    for(const v of verts){const p=vertexPoint(v);for(const d of v)positions.set(d,p);}
    let paths="";const items=[],counts=[],occupied=[];
    for(const [ei,e] of g.stable_edges.entries()){
      const a=positions.get(e.start),b=positions.get(e.end),name=e.id;
      let d,physical=name,controls=[];
      if(family==="tetra")d=`M${a} L${b}`;
      else if(name==="t"||name==="T") {physical="t";controls=[[136,6]];d=`M${a} Q${controls[0]} ${b}`;}
      else if(name==="u"||name==="U"){physical="u";controls=[[136,135]];d=`M${a} Q${controls[0]} ${b}`;}
      else if(name==="l+"||name==="l-"){
        physical="ℓ";
        controls=a[0]===b[0]&&a[1]===b[1]?[[a[0]-55,a[1]-60],[a[0]-55,a[1]+60]]:[[12,70]];
        d=`M${a} ${controls.length===2?"C":"Q"}${controls.map(p=>p.join(",")).join(" ")} ${b}`;
      }else if(name==="m+"||name==="m-"){
        physical="m";
        controls=a[0]===b[0]&&a[1]===b[1]?[[a[0]-46,a[1]-35],[a[0]-46,a[1]+35]]:[[148,70]];
        d=`M${a} ${controls.length===2?"C":"Q"}${controls.map(p=>p.join(",")).join(" ")} ${b}`;
      }else {
        const end=name[0],cross=config.counts["mode"+end]==="cross";
        physical=end+(cross?" cross":" arm");
        if(cross){controls=[[end==="A"?111:235,70]];d=`M${a} Q${controls[0]} ${b}`;}
        else d=`M${a} L${b}`;
      }
      const pts=sampleCurve(a,b,controls);curves.set(e.start,pts);curves.set(e.end,pts.slice().reverse());
      const ids=controlIDs(family,e,config),related=family==="double"&&config.operation==="repeated"&&physical==="ℓ"&&Number(config.counts.countA)===0?["repeat-countA","repeat-modeA"]:family==="double"&&config.operation==="repeated"&&physical==="m"&&Number(config.counts.countB)===0?["repeat-countB","repeat-modeB"]:[];
      paths+=linkedPath(d,{name:physical,m:e.multiplicity,color:graphPalette[ei%graphPalette.length],controls:ids,related,title:`${physical} ×${e.multiplicity}; F${e.faces.join(" / F")}`});
      items.push({name:physical,m:e.multiplicity,faces:e.faces,controls:ids,related});
      counts.push({name:physical,m:e.multiplicity,points:pts,color:graphPalette[ei%graphPalette.length],controls:ids,related});
    }
    const polygons=g.faces.map(f=>f.filter(d=>curves.has(d)).flatMap(d=>curves.get(d).slice(0,-1)));
    const nodes=verts.map(v=>{const p=vertexPoint(v);return `<circle cx="${p[0]}" cy="${p[1]}" r="4.2" fill="#285f65" stroke="white"/>`;}).join("");
    const preferred=Object.fromEntries([1,2,3].map(i=>[i,BouquetLayout.interiorLabel(polygons[i])]));
    const callouts=family==="tetra"?{0:[250,13],1:[132,137],2:[35,22],3:[230,22]}:{0:[136,3],1:[133,151],2:[-17,40],3:[255,70]};
    const faces=faceText(polygons,preferred,callouts,occupied);
    if(family==="tetra")occupied.push([28,115,13,14],[229,115,13,14],[128,-4,13,14],[141,74,13,14]);
    const countLabels=edgeCountLabels(counts,{vertices:verts.map(vertexPoint),occupied,box:family==="tetra"?[10,-16,260,160]:[-30,-12,310,178]});
    const vertexLabels=family==="tetra"?'<g class="graph-vertex-label"><text x="28" y="119">A</text><text x="229" y="119">B</text><text x="128" y="0">C</text><text x="141" y="78">D</text></g>':"";
    return `<svg class="labelled-dual-svg" viewBox="${family==="tetra"?"10 -16 260 160":"-30 -12 310 178"}" role="${config.operation==="repeated"?"group":"img"}" aria-label="Labelled ${escape(family)} dual graph; each nearby number counts edges in that chain${config.operation==="repeated"?". Select an edge to focus its decoration control":""}"><g>${paths}${nodes}</g><g text-anchor="middle" class="case-face-labels">${faces}</g>${vertexLabels}${countLabels}</svg>${graphKey(items)}`;
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
  function updateLameGraph() { /* Compatibility no-op: this folder exposes TV only. */ }
  function repeatedControls(){
    const c=caseConfig.counts;
    const field=(key,label)=>`<div class="repeat-count-field"><label for="repeat-${key}">${label.split(" · ").map((part,i)=>`<span class="${i?"repeat-face-pair":"repeat-edge-name"}">${escape(part)}</span>`).join(" ")}</label><input id="repeat-${key}" name="${key}" type="number" min="0" max="32" step="1" value="${c[key]}"></div>`;
    const local=end=>`<div class="repeat-local"><div><label for="repeat-mode${end}">At vertex ${end}</label><select id="repeat-mode${end}" name="mode${end}"><option value="arm"${c["mode"+end]==="arm"?" selected":""}>Arm · 0 unchanged</option><option value="cross"${c["mode"+end]==="cross"?" selected":""}>Cross split · ≥1</option></select></div>${field("count"+end,"Count")}</div>`;
    const controls=seed==="double"?`<div class="repeat-count-grid">${field("qt","Pairs on t · F₀ / F₁")}${field("qu","Pairs on u · F₀ / F₁")}${field("ql","Pairs on ℓ · F₀ / F₂")}${field("qm","Pairs on m · F₁ / F₃")}</div><div class="repeat-local-grid">${local("A")}${local("B")}</div>`:
      `<div class="repeat-count-grid">${[["qAB","AB · F₀ / F₁"],["qAC","AC · F₀ / F₂"],["qAD","AD · F₁ / F₂"],["qBC","BC · F₀ / F₃"],["qBD","BD · F₁ / F₃"],["qCD","CD · F₂ / F₃"]].map(([k,label])=>field(k,"Pairs on "+label)).join("")}</div>`;
    return `<form id="repeat-settings" class="repeat-settings"${operation==="repeated"?"":" hidden"}><p class="repeat-instruction">Edit counts to preview the graph, then apply.</p>${controls}<div class="repeat-actions"><button type="submit" class="button" id="repeat-apply">${applyButtonLabel()}</button><button type="button" id="repeat-reset" class="button secondary">Reset counts</button><span class="caption" id="repeat-feedback" role="status"></span></div><details class="case-rules control-help"><summary>Count rules &amp; limits</summary><p>${applyInstruction()} Each pair adds 2 hemispheres. Start with counts 0–3 for a faster first comparison.</p>${seed==="double"?'<p>An arm count k adds k hemispheres: odd counts start with a sprout, positive even counts with a loop-block split, followed by pair insertions. A cross-split count c adds 2c hemispheres and must be at least 1. Choose one local type at each vertex; arm and cross split cannot both be active there.</p>':""}<p>Whole-number counts 0–32 per control; total N ≤ 144 hemispheres; every face weight nᵢ ≤ 72. These bounds apply together. They are numerical limits, not mathematical restrictions. Larger counts and mixed decorations can be slower.</p></details></form>`;
  }
  function bouquetLabel(core) {
    if(!core)return "Undecorated core";
    const sprouts=core.sprouts??[],splits=core.diagonals??[];
    return `${sprouts.length?`Sprout ${sprouts.map(i=>`F${i}`).join(", ")}`:"No sprouts"} · ${splits.length?`split ${splits.map(d=>d.join("–")).join(", ")}`:"unsplit center"}`;
  }
  const bouquetCores=()=>Object.values(BouquetModel.cores);
  const sproutKey=core=>(core.sprouts??[]).join("");
  function bouquetEdgeLabel(edge,index,core) {
    const loop=edge.id.match(/^l([123])/),bridge=edge.id.match(/^b([123])/);
    if(loop)return `ℓ${loop[1]}`;
    if(bridge)return `b${bridge[1]}`;
    return `h${core.stable_edges.slice(0,index+1).filter(e=>!/^l[123]|^b[123]/.test(e.id)).length}`;
  }
  function bouquetControls() {
    const core=BouquetModel.cores[caseConfig.core];
    const subsets=[[],[1],[2],[3],[1,2],[1,3],[2,3],[1,2,3]];
    return `<form id="bouquet-settings" class="repeat-settings"${operation==="repeated"?"":" hidden"}><p class="repeat-instruction">Preview the graph, then apply.</p><div class="bouquet-core-grid"><div><label for="bouquet-sprouts">Sprouted monogons</label><select id="bouquet-sprouts">${subsets.map(s=>`<option value="${s.join("")}"${s.join("")===sproutKey(core)?" selected":""}>${s.length?s.map(i=>`F${i}`).join(", "):"None"}</option>`).join("")}</select></div><div><label for="bouquet-core">Admissible center split</label><select id="bouquet-core" name="bouquetCore"></select></div></div><div class="bouquet-core-preview" id="bouquet-core-preview"></div><div id="bouquet-count-fields" class="repeat-count-grid"></div><div class="repeat-actions"><button type="submit" class="button" id="bouquet-apply">${applyButtonLabel()}</button><button type="button" id="bouquet-reset" class="button secondary">Reset pair counts</button><span class="caption" id="bouquet-feedback" role="status"></span></div><details class="case-rules control-help"><summary>Decoration rules &amp; limits</summary><p>Choose sprouts, a center split, and pair counts. Only admissible choices are listed. ${applyInstruction()} Start with counts 0–3 for a faster first comparison.</p><p>Whole-number pair counts 0–32 on each chain; total N ≤ 144 hemispheres; every face weight nᵢ ≤ 72. These bounds apply together. One pair adds 2 hemispheres. The 88 admissible cores and their pair insertions are mathematical families; the caps only limit numerical drawing.</p></details></form>`;
  }
  function bouquetSplitSVG(core) {
    const darts=core.external_darts??[],d=darts.length,cx=92,cy=74,r=50;
    const point=(i,extra=0)=>{const a=-Math.PI/2+(i+.5)*2*Math.PI/d;return[cx+(r+extra)*Math.cos(a),cy+(r+extra)*Math.sin(a)];};
    const gap=i=>{const a=-Math.PI/2+(i+1)*2*Math.PI/d;return[cx+r*Math.cos(a),cy+r*Math.sin(a)];};
    const edges=Array.from({length:d},(_,i)=>gap(i).join(",")).join(" ");
    return `<svg class="bouquet-split-svg" viewBox="0 0 185 150" role="img" aria-label="Noncrossing center split; gap indices label the polygon"><polygon points="${edges}" fill="#285f6508" stroke="#8d9b8b" stroke-width="1.4"/>${(core.diagonals??[]).map(([i,j])=>`<path d="M${gap(i)} L${gap(j)}" stroke="#b55e3c" stroke-width="2.4"/>`).join("")}${darts.map((dart,i)=>{const p=point(i,15),q=gap(i),label=dart.replace(/^l([123])([+-])$/,"ℓ$1$2");return `<text x="${p[0]}" y="${p[1]+3}" text-anchor="middle">${escape(label)}</text><circle cx="${q[0]}" cy="${q[1]}" r="2.6" fill="#496351"/><text class="bouquet-gap-label" x="${q[0]+5}" y="${q[1]+3}">${i}</text>`;}).join("")}</svg><details class="case-rules control-help"><summary>Read this split</summary><p class="caption">Chord endpoints use the numbered gaps. Sprouts replace a pair of loop darts by one bridge dart. Only noncrossing splits are selectable.</p></details>`;
  }
  function bouquetGraph(config) {
    const core=BouquetModel.cores[config.core],vs=core.vertices.filter(v=>v.length!==2),index=new Map(vs.flatMap((v,i)=>v.map(d=>[d,i]))),alpha=new Map(core.pairs.flatMap(([a,b])=>[[a,b],[b,a]]));
    const centers=[[90,33],[222,95],[90,157]],middle=[142,95],direction=i=>{const p=centers[i-1],length=Math.hypot(p[0]-middle[0],p[1]-middle[1]);return[(p[0]-middle[0])/length,(p[1]-middle[1])/length];};
    const edges=core.stable_edges.map((e,i)=>({...e,index:i,start:e.chain[0],end:alpha.get(e.chain.at(-1))}));
    const adjacent=vs.map(()=>[]),anchors=vs.map(()=>[]);
    for(const e of edges){const a=index.get(e.start),b=index.get(e.end),loop=e.id.match(/^l([123])/);if(loop){const u=direction(Number(loop[1]));anchors[a].push([middle[0]+45*u[0],middle[1]+45*u[1]]);anchors[b].push([middle[0]+45*u[0],middle[1]+45*u[1]]);}else if(a!==b){adjacent[a].push(b);adjacent[b].push(a);}}
    const positions=vs.map(()=>middle.slice());
    for(let step=0;step<220;step++)for(let i=0;i<vs.length;i++)for(let k=0;k<2;k++)positions[i][k]=(adjacent[i].reduce((s,j)=>s+positions[j][k],0)+anchors[i].reduce((s,p)=>s+2*p[k],0))/(adjacent[i].length+2*anchors[i].length||1);
    let paths="";const facePositions=BouquetLayout.facePositions(core,edges,positions,index),items=[],counts=[],occupied=[];
    for(const e of edges){const a=positions[index.get(e.start)],b=positions[index.get(e.end)],loop=e.id.match(/^l([123])/),q=config.edgeCounts?.[e.id]??0,m=e.epsilon+2*q;let d,controls=[];
      if(loop){const i=Number(loop[1]),u=direction(i),v=[-u[1],u[0]],c1=[a[0]+78*u[0]-37*v[0],a[1]+78*u[1]-37*v[1]],c2=[b[0]+78*u[0]+37*v[0],b[1]+78*u[1]+37*v[1]];controls=[c1,c2];d=`M${a} C${c1} ${c2} ${b}`;}
      else d=`M${a} L${b}`;
      const name=bouquetEdgeLabel(e,e.index,core);
      const ids=config.operation==="repeated"?[`bouquet-bq${e.index}`]:[];
      paths+=linkedPath(d,{name,m,color:graphPalette[e.index%graphPalette.length],controls:ids,title:`${name} ×${m}; F${e.faces.join(" / F")}`});
      items.push({name,m,faces:e.faces,controls:ids});
      counts.push({name,m,points:sampleCurve(a,b,controls),color:graphPalette[e.index%graphPalette.length],controls:ids});
    }
    const polygons=BouquetLayout.facePolygons(core,edges,positions,index),faces=faceText(polygons,facePositions,{0:[260,12],1:[45,8],2:[280,95],3:[45,190]},occupied);
    return `<svg class="bouquet-dual-svg labelled-dual-svg" viewBox="0 -12 300 218" role="${config.operation==="repeated"?"group":"img"}" aria-label="Labelled decorated three-loop dual graph; each nearby number counts edges in that chain. Select an edge to focus its pair-count control"><g>${paths}${positions.map((p,i)=>`<circle cx="${p[0]}" cy="${p[1]}" r="4.2" fill="#285f65" stroke="white"><title>Core vertex ${i+1}</title></circle>`).join("")}</g><g class="case-face-labels">${faces}</g>${edgeCountLabels(counts,{vertices:positions,occupied,box:[0,-12,300,218]})}</svg>${graphKey(items)}`;
  }
  let syncBouquetControls=null;
  function mountBouquetControls() {
    const form=document.getElementById("bouquet-settings");if(!form)return;
    const sprouts=document.getElementById("bouquet-sprouts"),select=document.getElementById("bouquet-core"),fields=document.getElementById("bouquet-count-fields"),feedback=document.getElementById("bouquet-feedback");
    const read=()=>{const values={bouquetCore:select.value};for(let i=0;i<6;i++)values[`bq${i}`]=document.getElementById(`bouquet-bq${i}`)?.value??0;return values;};
    const budget=()=>{const next=BouquetModel.make("repeated",read()),ok=next.numericalBudget.supported;feedback.textContent=ok?`${next.hemispheres} / 144 hemispheres · largest weight ${Math.max(...next.weights)} / 72 · shift (${next.shift.map(v=>v===.5?"½":"0").join(", ")})`:`${next.hemispheres} hemispheres · reduce counts: total ≤144 and each weight ≤72.`;feedback.classList.toggle("repeat-invalid",!ok);document.getElementById("bouquet-apply").disabled=!ok;
      if(form.checkValidity?.()!==false)previewGraph(next);
    };
    const updateCore=(retain=false)=>{selectedGraphControl=null;const core=BouquetModel.cores[select.value];fields.innerHTML=core.stable_edges.map((e,i)=>`<div class="repeat-count-field"><label for="bouquet-bq${i}"><span class="repeat-edge-name">Pairs on ${escape(bouquetEdgeLabel(e,i,core))}</span> <span class="repeat-face-pair">F${e.faces.join(" / F")}</span></label><input id="bouquet-bq${i}" type="number" name="bq${i}" min="0" max="32" step="1" value="${retain?caseConfig.counts[`bq${i}`]??0:0}"><small class="caption">Core: ${e.epsilon} edge${e.epsilon===1?"":"s"}</small></div>`).join("");document.getElementById("bouquet-core-preview").innerHTML=bouquetSplitSVG(core);budget();};
    const updateChoices=(initial=false)=>{const list=bouquetCores().filter(c=>sproutKey(c)===sprouts.value);select.innerHTML=list.map((c,i)=>`<option value="${escape(c.name)}">${(c.diagonals??[]).length?`Split ${c.diagonals.map(d=>d.join("–")).join(", ")}`:"Unsplit center"} · core ${i+1}/${list.length}</option>`).join("");if(initial&&list.some(c=>c.name===caseConfig.core))select.value=caseConfig.core;updateCore(initial);};
    sprouts.addEventListener("change",()=>updateChoices());select.addEventListener("change",()=>updateCore());form.addEventListener("input",budget);
    form.addEventListener("submit",event=>{event.preventDefault();if(!form.reportValidity())return;const config=BouquetModel.make("repeated",read());if(!config.numericalBudget.supported){budget();return;}commitCase(config,"repeated");});
    document.getElementById("bouquet-reset").addEventListener("click",()=>{for(const input of fields.querySelectorAll("input"))input.value="0";budget();});
    syncBouquetControls=()=>{sprouts.value=sproutKey(BouquetModel.cores[caseConfig.core]);updateChoices(true);};
    updateChoices(true);
  }

  function renderSupport(){
    const support=document.getElementById("case-support");
    if(!support)return;
    const limits=`${seed==="bouquet"?"Whole-number pair counts":"Whole-number counts"} 0–32 per control; total N ≤ 144 hemispheres; every face weight nᵢ ≤ 72. All three bounds apply together.`;
    support.innerHTML=`<details class="case-rules control-help"><summary>Numerical range &amp; speed</summary><p>${limits}</p><p class="caption">These are drawing limits, not mathematical restrictions. Start with small counts for a faster comparison. ${isLimitsPage()?"Apply decorations first, choose the source path, then use Compute approach.":"Apply &amp; compute starts the interior image and cusp series. Refine either preview with the controls below."}</p></details>`;
  }
  function renderSummary(){
    renderSupport();
    const summary=document.getElementById("case-summary");
    const badge=(seed==="double"?"DOUBLE-EDGE":seed==="bouquet"?"THREE-LOOP":"TETRAHEDRON");
    const shiftLabel=shift.map(x=>x===.5?"½":"0").join(", ");
    const familyRule=seed==="tetra"?"The tetrahedral seed admits pair insertion on its six labelled edges; splitting and sprouting are unavailable.":seed==="bouquet"?`All eight sprout subsets and their noncrossing center splits give 88 admissible cores. Pair insertions preserve admissibility. Selected core: ${escape(caseConfig.core)}.`:"At each original vertex choose either an arm or a cross split; these local choices are mutually exclusive. Pairs may be inserted on each labelled chain.";
    const wasOpen=summary?.querySelector("details")?.open;
    if(summary)summary.innerHTML=`<p class="case-active-title"><span class="case-badge">${badge}</span> <strong>${escape(active.title)}</strong></p>${`<p class="case-essential"><strong>${caseConfig.hemispheres} / 144 hemispheres</strong> · largest face weight ${Math.max(...caseConfig.weights)} / 72 · character shift (${shiftLabel}) modulo 1.</p><details class="case-rules control-help"><summary>Case rules &amp; limits</summary><p>${familyRule}</p><p>Weights (n₀, n₁, n₂, n₃) match faces F₀, F₁, F₂, F₃; their sum N is the total number of hemispheres. Numerical drawing uses integer counts 0–32 per control, N ≤ 144, and each nᵢ ≤ 72. These are simultaneous computation limits.</p><p>The full labelled graph selects the numerical branch. Equal weights or monodromy signs need not give the same marked τ-image. Smaller counts give faster first comparisons; unresolved samples are never replaced by a theoretical endpoint.</p></details>`}`;
    if(wasOpen && summary?.querySelector("details"))summary.querySelector("details").open=true;
  }
  function renderCaseGraph(){
    const graph=document.getElementById("case-graph");
    if(graph){
      const compare=(operation!=="seed");
      const expanded=graph.classList.contains("case-graph-expanded");
      graph.innerHTML=`<div class="case-graph-toolbar"><strong>Labelled dual ${compare?"graphs":"graph"}</strong><button type="button" id="case-graph-enlarge" aria-expanded="false" aria-controls="case-graph-comparison">Enlarge ${compare?"graphs":"graph"}</button></div><div id="case-graph-comparison" class="case-graph-comparison${compare?"":" case-graph-single"}"><figure class="case-graph-base">${graphSVG(seed,"seed")}<figcaption>${escape(active.baseTitle)}</figcaption></figure>${compare?`<span class="case-graph-arrow" aria-hidden="true">→</span><figure class="case-graph-result">${`<div id="case-result-graph">${seed==="bouquet"?bouquetGraph(caseConfig):repeatedGraph(caseConfig)}</div>`}<figcaption${' id="case-result-caption"'}><strong>Selected:</strong> ${escape((seed==="bouquet"?bouquetLabel(BouquetModel.cores[caseConfig.core]):opNames[operation]))}</figcaption></figure>`:""}</div><p id="graph-control-feedback" class="graph-control-feedback" role="status" aria-live="polite">${(operation==="repeated")?"Click an edge or count to highlight its control.":"Fixed preset · choose Repeated decorations to edit chains."}</p><details class="case-rules control-help"><summary>Read the dual graph</summary><p class="case-graph-note">Dots are vertices; F labels identify faces. <strong>Each nearby ×m badge gives the total number of edges in that chain</strong>; ×1 is a single edge. Intermediate bivalent vertices are suppressed. ${(compare?"Badge colours match the chains and the key below. ":"")}The drawing records topology, not edge lengths.</p><p>Click an edge, count badge, or key to highlight its control. Focus a control to highlight its edge. The base graph is a fixed reference.</p></details>`;
      const toggle=document.getElementById("case-graph-enlarge");
      if(expanded){toggle.setAttribute("aria-expanded","true");toggle.textContent="Compact view";}
      document.getElementById("case-graph-enlarge").addEventListener("click",event=>{
        const expanded=graph.classList.toggle("case-graph-expanded");
        event.currentTarget.setAttribute("aria-expanded",String(expanded));
        event.currentTarget.textContent=expanded?"Compact view":`Enlarge ${compare?"graphs":"graph"}`;
      });
      refreshGraphControlLinks();
    }
  }
  function previewGraph(config,{draft=config.id!==active.id}={}){
    
    const target=document.getElementById("case-result-graph"),caption=document.getElementById("case-result-caption");
    if(!target)return;
    target.innerHTML=seed==="bouquet"?bouquetGraph(config):repeatedGraph(config);
    if(caption)caption.innerHTML=`<strong>${draft?"Preview":"Selected"}:</strong> weights (${config.weights.join(", ")})${seed==="bouquet"?` · ${escape(bouquetLabel(BouquetModel.cores[config.core]))}`:""}${draft?" · apply to update τ":""}`;
    refreshGraphControlLinks();
  }
  function commitCase(config,nextOperation=operation){
    if((!config?.numericalBudget?.supported))return false;
    const previousId=active.id,previousOperation=operation,x=window.scrollX,y=window.scrollY;
    operation=nextOperation;caseConfig=config;shift=config.shift;
    scale=Math.min(300/(maxR+shift[0]),300/(seed==="bouquet"?top+shift[1]:Math.max(1,top+shift[1])));
    Object.assign(active,caseMetadata(config));
    const url=new URL(location.href);url.searchParams.set("seed",seed);url.searchParams.set("decoration",operation);
    for(const [key,value]of Object.entries(config.counts))url.searchParams.set(key,value);
    history.replaceState(null,"",url);
    const form=document.getElementById(seed==="bouquet"?"bouquet-settings":"repeat-settings");if(form)form.hidden=operation!=="repeated";
    document.getElementById("atlas-decoration").value=operation;
    if(previousOperation!==operation&&seed==="bouquet")syncBouquetControls?.();
    renderSummary();
    if(previousOperation!==operation)renderCaseGraph();
    previewGraph(config,{draft:false});
    const feedback=document.getElementById(seed==="bouquet"?"bouquet-feedback":"repeat-feedback");
    if(feedback)feedback.textContent="Decoration applied. The labelled graph and numerical case are updated.";
    window.dispatchEvent(new CustomEvent("atlas:casechange",{detail:{previousId,active}}));
    window.AtlasUI?.refreshNavigation();
    requestAnimationFrame(()=>{if(typeof window.scrollTo==="function")window.scrollTo(x,y);});
    return true;
  }

  const familyStates=new Map();
  let selectorsBound=false;
  function placeAttachmentControls() {
    if (!document.getElementById("lame-attachment-settings")) return;
    // Keep legacy HTML hosts hidden without moving them on each TV switch.
    for (const id of ["attachment-form", "decoration-status", "limits-legacy-attachment"]) {
      const element = document.getElementById(id);
      if (element) element.hidden = true;
    }
  }
  function legacyCounts(){
    const current=new URLSearchParams(location.search);
    const read=prefix=>["a","b","c"].map((key,i)=>{
      const field=document.getElementById(prefix+key);return field?Number(field.value):Number(current.get(key)??[1,0,0][i]);
    });
    const values=(!document.getElementById("limits-controls"))
      ? ["a","b","c"].map((key,i)=>Number(current.get(key)??[1,0,0][i]))
      : read(document.getElementById("limit-a")?"limit-":"");
    return values.every(x=>Number.isSafeInteger(x)&&x>=0)&&values.reduce((a,b)=>a+b,0)<=20?values:[1,0,0];
  }
  function switchFamily(family){
    if(!allowedSeeds.includes(family)||!Object.hasOwn(names,family))return false;
    if(family===seed)return true;
    const previousSeed=seed,x=window.scrollX,y=window.scrollY,retainedLegacyCounts=legacyCounts();
    const before={previousSeed,nextSeed:family,source:null};
    window.dispatchEvent(new CustomEvent("atlas:beforefamilychange",{detail:before}));
    const current=new URLSearchParams(location.search);
    const source=before.source??{r:Number(current.get(document.getElementById("limits-controls")?"fromR":"r")??active.start[0]),s:Number(current.get(document.getElementById("limits-controls")?"fromS":"s")??active.start[1])};
    const normalized={r:source.r/maxR,s:source.s/top};
    familyStates.set(seed,{operation,config:caseConfig,counts:null});
    const saved=familyStates.get(family);
    seed=family;
    operations = window.AtlasFamily.operations[seed];
    operation=saved?.operation??"repeated";
    if(!operations.includes(operation))operation=operations[0];
    top=seed==="tetra"?1:.5;maxR=seed==="bouquet" ? .5 : 1;
    caseConfig=(saved?.config??(seed==="bouquet"?BouquetModel.make(operation,BouquetModel.preset(operation)):DecorationModel.make(seed,operation,DecorationModel.preset(seed,operation))));
    shift=caseConfig?.shift??[0,0];
    scale=Math.min(300/(maxR+shift[0]),300/(seed==="bouquet"?top+shift[1]:Math.max(1,top+shift[1])));
    const counts=saved?.counts??retainedLegacyCounts;
    const validCounts=counts.every(v=>Number.isSafeInteger(v)&&v>=0)&&counts.reduce((a,b)=>a+b,0)<=20?counts:[1,0,0];
    Object.assign(active,caseMetadata(caseConfig));
    const nextSource=project(normalized.r*maxR,normalized.s*top,.04);
    const url=new URL(location.href);url.searchParams.set("seed",seed);url.searchParams.set("decoration",operation);
    for(const key of [...DecorationModel.queryKeys,...BouquetModel.queryKeys])url.searchParams.delete(key);
    if(caseConfig)for(const [key,value]of Object.entries(caseConfig.counts))url.searchParams.set(key,value);
    for(const [i,key]of ["a","b","c"].entries())url.searchParams.set(key,validCounts[i]);
    url.searchParams.set("r",String(nextSource.r));url.searchParams.set("s",String(nextSource.s));
    if(document.getElementById("limits-controls")||url.searchParams.has("fromR")){url.searchParams.set("fromR",String(nextSource.r));url.searchParams.set("fromS",String(nextSource.s));}
    url.searchParams.delete("fromOverview");url.searchParams.delete("start");
    history.replaceState(null,"",url);
    for(const [i,key]of ["a","b","c"].entries())for(const prefix of ["","limit-"]){const input=document.getElementById(prefix+key);if(input)input.value=String(validCounts[i]);}
    
    mount();
    
    window.dispatchEvent(new CustomEvent("atlas:familychange",{detail:{previousSeed,seed,active,source:nextSource,counts:validCounts}}));
    window.AtlasUI?.refreshNavigation();
    const selector=document.getElementById("atlas-seed");if(selector)selector.focus?.({preventScroll:true});
    requestAnimationFrame(()=>{if(typeof window.scrollTo==="function")window.scrollTo(x,y);});
    return true;
  }
  function mount(){
    const selector=document.getElementById("case-selector");if(!selector)return;
    const existing=document.getElementById("atlas-seed");
    if(!existing){
    selector.innerHTML=`<div><label for="atlas-seed">Seed family</label><select id="atlas-seed">${Object.entries(names).filter(([key])=>allowedSeeds.includes(key)).map(([key,label])=>`<option value="${key}"${key===seed?" selected":""}>${escape(label)}</option>`).join("")}</select></div><div><label for="atlas-decoration">Decoration</label><select id="atlas-decoration">${operations.map(key=>`<option value="${key}"${key===operation?" selected":""}>${escape(seed==="bouquet"&&key==="repeated"?"All admissible decorations · choose core and counts":opNames[key])}</option>`).join("")}</select></div>`;
    }
    const familySelect=document.getElementById("atlas-seed"),decorationSelect=document.getElementById("atlas-decoration");
    familySelect.value=seed;
    decorationSelect.innerHTML=operations.map(key=>`<option value="${key}"${key===operation?" selected":""}>${escape(seed==="bouquet"&&key==="repeated"?"All admissible decorations · choose core and counts":opNames[key])}</option>`).join("");
    for(const id of ["repeat-settings","bouquet-settings"]){const previous=document.getElementById(id);if(previous)previous.remove();}
    syncBouquetControls=null;
    selector.insertAdjacentHTML("beforeend",seed==="bouquet"?bouquetControls():repeatedControls());
    renderSummary();renderCaseGraph();
    const change=family=>{
      if(family===seed&&true){
        const next=document.getElementById("atlas-decoration").value;
        const config=seed==="bouquet"?BouquetModel.make(next,next==="repeated"?caseConfig.counts:BouquetModel.parse(next,new URLSearchParams())):DecorationModel.make(seed,next,next==="repeated"?caseConfig.counts:DecorationModel.preset(seed,next));
        if(seed!=="bouquet")for(const [key,value]of Object.entries(config.counts)){const input=document.getElementById("repeat-"+key);if(input)input.value=value;}
        if(seed==="double")for(const end of ["A","B"]){const input=document.getElementById("repeat-count"+end);if(input)input.min=config.counts["mode"+end]==="cross"?"1":"0";}
        commitCase(config,next);return;
      }
      switchFamily(family);
    };
    if(!selectorsBound){
      document.getElementById("atlas-seed").addEventListener("change",event=>change(event.target.value));
      document.getElementById("atlas-decoration").addEventListener("change",()=>change(seed));
      selectorsBound=true;
    }
    if(seed==="bouquet")mountBouquetControls();
    const repeatForm=document.getElementById("repeat-settings");
    if(repeatForm){
      const keys=seed==="double"?DecorationModel.doubleKeys:DecorationModel.tetraKeys;
      const updateMinimum=end=>{const input=document.getElementById("repeat-count"+end),cross=document.getElementById("repeat-mode"+end).value==="cross";input.min=cross?"1":"0";if(cross&&Number(input.value)<1)input.value="1";};
      const readControls=()=>Object.fromEntries(keys.map(key=>[key,document.getElementById("repeat-"+key).value]));
      const updateBudget=()=>{
        const config=DecorationModel.make(seed,"repeated",readControls()),feedback=document.getElementById("repeat-feedback"),supported=config.numericalBudget.supported;
        feedback.textContent=supported?`${config.hemispheres} / 144 hemispheres · largest weight ${Math.max(...config.weights)} / 72`:`${config.hemispheres} hemispheres · outside numerical budget. Reduce counts to total ≤144 and each weight ≤72.`;
        feedback.classList.toggle("repeat-invalid",!supported);document.getElementById("repeat-apply").disabled=!supported;
        if(repeatForm.checkValidity?.()!==false)previewGraph(config);
      };
      if(seed==="double")for(const end of ["A","B"]){updateMinimum(end);document.getElementById("repeat-mode"+end).addEventListener("change",()=>{updateMinimum(end);updateBudget();});}
      repeatForm.addEventListener("input",updateBudget);updateBudget();
      repeatForm.addEventListener("submit",event=>{
        event.preventDefault();if(!repeatForm.reportValidity())return;
        const next=readControls();
        const config=DecorationModel.make(seed,"repeated",next);
        if(!config.numericalBudget.supported){updateBudget();return;}
        commitCase(config,"repeated");
      });
      document.getElementById("repeat-reset").addEventListener("click",()=>{
        const zero=DecorationModel.preset(seed,"seed");for(const key of keys)document.getElementById("repeat-"+key).value=zero[key];
        if(seed==="double")for(const end of ["A","B"])updateMinimum(end);
        updateBudget();document.getElementById("repeat-feedback").textContent="Counts reset. Graph preview updated; Apply to update τ.";
      });
    }
    placeAttachmentControls();
    document.documentElement.dataset.seed=seed;
  }
  mount();
  
  function boundaryPoint(edge,position){if(seed==="bouquet")return edge===0?[.5,.5*position]:edge===1?[.5*position,.5]:[.5*position,.5*(1-position)];return edge===0?[1,top*position]:edge===1?[position,top]:[position,top*(1-position)];}
  return {get isNew(){return true;},active,project,map,svgPoint,svgInverse,renderMonodromy,asymptotic,boundaryPoint,
    updateLameGraph,commitCase,previewGraph,switchFamily,
    graphRenderers:{lameGraph,repeatedGraph,bouquetGraph},
    drawDomain:(canvas,options)=>window.AtlasView.drawDomain(canvas,{...options,markedStrip:true})};
})();

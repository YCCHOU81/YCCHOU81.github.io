"use strict";

/* Labelled metric cells. R,S are geometric characters; the analytic evaluator
   receives (r,s)=(S,R). No modular parking or inversion is applied here. */
window.SeedCases = (() => {
  const query = new URLSearchParams(location.search);
  const seed = ["double", "tetra"].includes(query.get("seed")) ? query.get("seed") : "lame";
  const operations = seed === "double" ? ["seed", "pair", "split", "sprout", "repeated"] : seed === "tetra" ? ["seed", "pair", "repeated"] : ["attachment"];
  const operation = operations.includes(query.get("decoration")) ? query.get("decoration") : operations[0];
  const names = {lame:"Classical Lamé · three-loop seed",double:"Double-edge seed · (1,1,0,0)",tetra:"Tetrahedral seed · (1,1,1,1)"};
  const opNames = {seed:"Seed only",pair:"Insert a pair on F₀ / F₁",split:"Split the F₂ loop-dart block",sprout:"Sprout the F₂ monogon",repeated:"Repeated decorations · choose counts",attachment:"Hemisphere attachments (a,b,c)"};
  const weights = seed === "double" ? {seed:[1,1,0,0],pair:[2,2,0,0],split:[3,1,0,0],sprout:[2,1,0,0]} : {seed:[1,1,1,1],pair:[2,2,1,1]};
  const caseConfig=seed!=="lame"?DecorationModel.make(seed,operation,DecorationModel.parse(seed,operation,query)):null;
  const shift = caseConfig?.shift??[0,0];
  const top = seed === "tetra" ? 1 : .5;
  const active = {id:caseConfig?.id??`${seed}-${operation}`,seed,operation,weights:caseConfig?.weights??[1,0,0,0],shift,top,caseConfig,isRepeated:operation==="repeated",
    title:seed === "lame" ? names.lame : `${seed === "double" ? "Double-edge" : "Tetrahedron"} · ${opNames[operation]} · (${caseConfig.weights.join(",")})`,
    baseTitle:names[seed],baseId:`${seed}-seed`,start:seed === "tetra" ? [.58,.62] : [.58,.31],
    sourceVertices:[[1,0],[0,top],[1,top]],
    edgeLabels:[`R → 1 · ${caseConfig&&!caseConfig.connector&&operation==="repeated"?"expected ":""}cusp ∞`,`S → ${top === 1 ? "1" : "½"} · ${caseConfig&&!caseConfig.connector&&operation==="repeated"?"expected ":""}cusp 0`,`${top === 1 ? "R + S" : "R + 2S"} → 1 · ${caseConfig&&!caseConfig.connector&&operation==="repeated"?"expected ":""}cusp ${top === 1 ? "1" : "½"}`]};

  const escape = value => String(value).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
  const compactVertices = inset => active.sourceVertices.map(v => v.map((x,j) => (1-3*inset)*x+inset*active.sourceVertices.reduce((sum,w)=>sum+w[j],0)));
  function project(R,S,inset=.04) {
    const vertices=compactVertices(inset);
    if (![R,S].every(Number.isFinite)) return {r:active.start[0],s:active.start[1]};
    if (R <= 1-inset && S <= top*(1-inset) && R+S/top >= 1+inset) return {r:R,s:S};
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
  const scale=Math.min(300/(1+shift[0]),300/Math.max(1,top+shift[1]));
  function svgPoint(R,S,attached=false) {const p=map(R,S,attached);return [48+scale*p.r,350-scale*p.s];}
  function svgInverse(x,y,attached=false) {return {r:(x-48)/scale-(attached?shift[0]:0),s:(350-y)/scale-(attached?shift[1]:0)};}
  function renderMonodromy(svg,{point,attached=false,title="Monodromy triangle"}={}) {
    const color=attached?"#b55e3c":"#285f65", xy=p=>svgPoint(...p,attached).join(",");
    const xMax=1+shift[0],yMax=Math.max(1,top+shift[1]);
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
    if (![0,1,2].includes(edge) || seed==="lame") return null;
    let row=tables[base?active.baseId:active.id]?.[edge],rat=null;
    if(!row&&!base&&caseConfig?.connector){
      const n=caseConfig.connector.q+1,tetra=seed==="tetra";
      let bin=1;for(let j=1;j<=n;j++)bin=bin*(n+j)/j;
      rat=tetra?[rational(1,2),rational(n,n+1),rational(1,n+1)][edge]:[rational(1,4),rational(n-1,n),rational(1-2*n,4*n)][edge];
      const exponent=tetra?[1,2/(n+1),2/(n+1)][edge]:[.5,2/n,1/(2*n)][edge];
      const coefficient=tetra?[Math.sqrt(8*n*(n+1)),2*bin,2*bin][edge]:[2*n*(n+1),n*bin,2*n*bin][edge];
      row=[rat.label,exponent,coefficient];
    }
    const chart=edge===0?"raw":edge===1?"zero":seed==="double"?"half":"one";
    const sinArgument=edge===0?(seed==="double"?2*Math.PI*S:Math.PI*S):Math.PI*R;
    const cusp=edge===0?"∞":edge===1?"0":seed==="double"?"½":"1";
    const chartLabel={raw:"τ",zero:"w = 1 − 1/τ",half:"w = τ/(1 − 2τ)",one:"w = 1/(1 − τ)"}[chart];
    if(!row)return {known:false,cusp,chart,chartLabel,x0:null,xLabel:"computed from the selected branch",alpha:null,C:null,
      deltaLabel:["1 − R",top===1?"1 − S":"½ − S",top===1?"R + S − 1":"R + 2S − 1"][edge],
      detail:`The selected chart targets cusp ${cusp}; the displayed tail is computed on this labelled decorated graph. A closed leading coefficient is not asserted for this mixed decoration.`};
    if(!rat){const parts=row[0].split("/").map(Number);rat=rational(parts[0],parts[1]??1);}
    return {known:true,cusp,chart,chartLabel,x0:rat.value,xLabel:rat.label,xNumerator:rat.numerator,xDenominator:rat.denominator,numerator:rat.numerator,denominator:rat.denominator,alpha:row[1],C:row[2]*Math.sin(sinArgument),
      deltaLabel:["1 − R",top===1?"1 − S":"½ − S",top===1?"R + S − 1":"R + 2S − 1"][edge],
      constantFactor:row[2],sineLabel:edge===0?(seed==="double"?"sin(2πS)":"sin(πS)"):"sin(πR)",
      detail:"Leading estimate on the selected geometric branch, away from source vertices; numerical branch matching is not an interval certificate."};
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
    return `<svg viewBox="0 0 255 120" role="img" aria-label="${escape(names[family])}, ${escape(opNames[op]??op)} dual graph"><g>${paths}${nodes}</g><g text-anchor="middle" class="case-face-labels">${labels}</g></svg>`;
  }
  function repeatedGraph(config){
    const g=config.graph,verts=g.vertices.filter(v=>v.length!==2),positions=new Map();
    const family=config.seed;
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
    let paths="",annotations="";
    for(const e of g.stable_edges){
      const a=positions.get(e.start),b=positions.get(e.end),name=e.id;
      let d,label=[(a[0]+b[0])/2,(a[1]+b[1])/2-5],physical=name;
      if(family==="tetra"){
        physical=name;d=`M${a} L${b}`;
        if(name==="AB")label[1]+=17;
        if(name==="AC"||name==="AD")label[0]-=13;
        if(name==="BC"||name==="BD")label[0]+=13;
        if(name==="CD")label[0]+=16;
      }else if(name==="t"||name==="T") {physical="t";d=`M${a} Q136 6 ${b}`;label=[136,21];}
      else if(name==="u"||name==="U"){physical="u";d=`M${a} Q136 135 ${b}`;label=[136,132];}
      else if(name==="l+"||name==="l-"){
        physical="ℓ";
        d=a[0]===b[0]&&a[1]===b[1]?`M${a} C${a[0]-55} ${a[1]-60},${a[0]-55} ${a[1]+60},${b}`:`M${a} Q12 70 ${b}`;
        label=[Math.max(17,a[0]-44),72];
      }else if(name==="m+"||name==="m-"){
        physical="m";
        d=a[0]===b[0]&&a[1]===b[1]?`M${a} C${a[0]-46} ${a[1]-35},${a[0]-46} ${a[1]+35},${b}`:`M${a} Q148 70 ${b}`;
        label=[a[0]-31,73];
      }else {
        const end=name[0],cross=config.counts["mode"+end]==="cross";
        physical=end+(cross?" cross":" arm");
        if(cross){d=`M${a} Q${end==="A"?111:235} 70 ${b}`;label=[end==="A"?109:231,73];}
        else {d=`M${a} L${b}`;label=[(a[0]+b[0])/2,(a[1]+b[1])/2-8];}
      }
      const changed=e.multiplicity!==1||name.includes("_");
      paths+=`<path d="${d}" fill="none" stroke="${changed?"#b55e3c":"#496351"}" stroke-width="${changed?2.8:2.2}"/>`;
      annotations+=`<text class="repeat-graph-label" x="${label[0]}" y="${label[1]}" text-anchor="middle">${escape(physical)}${e.multiplicity>1?` ×${e.multiplicity}`:""}</text>`;
    }
    const nodes=verts.map(v=>{const p=vertexPoint(v);return `<circle cx="${p[0]}" cy="${p[1]}" r="4.2" fill="#285f65" stroke="white"/>`;}).join("");
    const faces=family==="tetra"?'<text x="241" y="24">F₀</text><text x="131" y="96">F₁</text><text x="99" y="64">F₂</text><text x="164" y="64">F₃</text>':'<text x="136" y="6">F₀</text><text x="131" y="91">F₁</text><text x="45" y="54">F₂</text><text x="175" y="60">F₃</text>';
    return `<svg viewBox="0 0 265 145" role="img" aria-label="Labelled repeated-decoration graph; chain labels give the number of edges"><g>${paths}${nodes}</g><g text-anchor="middle" class="case-face-labels">${faces}</g><g>${annotations}</g></svg>`;
  }
  function repeatedControls(){
    const c=caseConfig.counts;
    const field=(key,label)=>`<div><label for="repeat-${key}">${label}</label><input id="repeat-${key}" name="${key}" type="number" min="0" max="32" step="1" value="${c[key]}"></div>`;
    const local=end=>`<div class="repeat-local"><div><label for="repeat-mode${end}">At vertex ${end}</label><select id="repeat-mode${end}" name="mode${end}"><option value="arm"${c["mode"+end]==="arm"?" selected":""}>Arm · 0 unchanged</option><option value="cross"${c["mode"+end]==="cross"?" selected":""}>Cross split · ≥1</option></select></div>${field("count"+end,"Count")}</div>`;
    const controls=seed==="double"?`<div class="repeat-count-grid">${field("qt","Pairs on t · F₀ / F₁")}${field("qu","Pairs on u · F₀ / F₁")}${field("ql","Pairs on ℓ · F₀ / F₂")}${field("qm","Pairs on m · F₁ / F₃")}</div><div class="repeat-local-grid">${local("A")}${local("B")}</div>`:
      `<div class="repeat-count-grid">${[["qAB","AB · F₀ / F₁"],["qAC","AC · F₀ / F₂"],["qAD","AD · F₁ / F₂"],["qBC","BC · F₀ / F₃"],["qBD","BD · F₁ / F₃"],["qCD","CD · F₂ / F₃"]].map(([k,label])=>field(k,"Pairs on "+label)).join("")}</div>`;
    return `<form id="repeat-settings" class="repeat-settings"${operation==="repeated"?"":" hidden"}><p class="caption">Choose labelled counts, then apply. One pair adds two hemispheres. A chain marked ×m contains m edges; intermediate vertices are suppressed in this diagram.</p>${controls}${seed==="double"?'<p class="caption">An arm count k adds k hemispheres: odd counts start with a sprout, positive even counts with a loop-block split, followed by pair insertions. A cross-split count c adds 2c hemispheres. Choose one local type at each vertex.</p>':""}<div class="repeat-actions"><button type="submit" class="button" id="repeat-apply">Apply decoration counts</button><button type="button" id="repeat-reset" class="button secondary">Reset counts</button><span class="caption" id="repeat-feedback" role="status"></span></div><p class="caption">Counts 0–32 per control. Numerical drawing supports at most 144 hemispheres and each face weight at most 72. Mixed configurations can take longer; only resolved points are shown.</p></form>`;
  }
  function mount(){
    const selector=document.getElementById("case-selector");if(!selector)return;
    selector.innerHTML=`<div><label for="atlas-seed">Seed family</label><select id="atlas-seed">${Object.entries(names).map(([key,label])=>`<option value="${key}"${key===seed?" selected":""}>${escape(label)}</option>`).join("")}</select></div><div><label for="atlas-decoration">Decoration</label><select id="atlas-decoration">${operations.map(key=>`<option value="${key}"${key===operation?" selected":""}>${escape(opNames[key])}</option>`).join("")}</select></div>${seed!=="lame"?repeatedControls():""}`;
    const summary=document.getElementById("case-summary");
    if(summary)summary.innerHTML=`<p class="case-active-title"><span class="case-badge">${seed==="lame"?"CLASSICAL LAMÉ":seed==="double"?"DOUBLE-EDGE":"TETRAHEDRON"}</span> <strong>${escape(active.title)}</strong></p><p class="caption">${seed==="lame"?"Choose a, b, c below. The graph identifies the three-loop base; the monodromy and τ plots compare the chosen attachment.":`${operation==="seed"?"Compare the seed with itself, or choose an admissible decoration above.":operation==="repeated"?`${caseConfig.hemispheres} hemispheres · ordinary character shift (${shift.map(x=>x===.5?"½":"0").join(", ")}) modulo 1. Each labelled count is retained when selecting the numerical branch.`:"The orange edges and vertices mark the selected operation."} Weights are ordered (n₀,n₁,n₂,n₃), matching faces F₀,F₁,F₂,F₃. ${seed==="tetra"?"This seed admits pair insertion; splitting and sprouting are unavailable.":operation==="repeated"?"At each original vertex, choose either an arm or a cross split. These choices are mutually exclusive.":""}`}</p>`;
    const graph=document.getElementById("case-graph");if(graph)graph.innerHTML=`<figure>${graphSVG(seed,"seed")}<figcaption>${escape(seed==="lame"?"Three-loop base":names[seed])}</figcaption></figure>${seed!=="lame"&&operation!=="seed"?`<span class="case-graph-arrow" aria-hidden="true">→</span><figure>${operation==="repeated"?repeatedGraph(caseConfig):graphSVG(seed,operation)}<figcaption>${escape(opNames[operation])}</figcaption></figure>`:""}<p class="case-graph-note">Dual graph · dots are vertices; F labels identify faces. ${operation==="repeated"?"×m gives the chain's edge multiplicity. ":""}Drawings encode topology, not edge lengths.</p>`;
    const change=family=>{
      const url=new URL(location.href);url.searchParams.set("seed",family);
      url.searchParams.set("decoration",family===seed?document.getElementById("atlas-decoration").value:family==="lame"?"attachment":"seed");
      if(family===seed&&url.searchParams.get("decoration")==="repeated"&&operation!=="repeated")for(const [k,v] of Object.entries(caseConfig.counts))url.searchParams.set(k,v);
      if(family!==seed)for(const key of ["r","s","fromR","fromS","path","t","start",...DecorationModel.queryKeys])url.searchParams.delete(key);
      location.href=url.href;
    };
    document.getElementById("atlas-seed").addEventListener("change",event=>change(event.target.value));
    document.getElementById("atlas-decoration").addEventListener("change",()=>change(seed));
    const repeatForm=document.getElementById("repeat-settings");
    if(repeatForm){
      const keys=seed==="double"?DecorationModel.doubleKeys:DecorationModel.tetraKeys;
      const updateMinimum=end=>{const input=document.getElementById("repeat-count"+end),cross=document.getElementById("repeat-mode"+end).value==="cross";input.min=cross?"1":"0";if(cross&&Number(input.value)<1)input.value="1";};
      const readControls=()=>Object.fromEntries(keys.map(key=>[key,document.getElementById("repeat-"+key).value]));
      const updateBudget=()=>{
        const config=DecorationModel.make(seed,"repeated",readControls()),feedback=document.getElementById("repeat-feedback"),supported=config.numericalBudget.supported;
        feedback.textContent=supported?`${config.hemispheres} hemispheres · weights (${config.weights.join(", ")})`:`${config.hemispheres} hemispheres · outside numerical budget. Reduce counts to total ≤144 and each weight ≤72.`;
        feedback.classList.toggle("repeat-invalid",!supported);document.getElementById("repeat-apply").disabled=!supported;
      };
      if(seed==="double")for(const end of ["A","B"]){updateMinimum(end);document.getElementById("repeat-mode"+end).addEventListener("change",()=>{updateMinimum(end);updateBudget();});}
      repeatForm.addEventListener("input",updateBudget);updateBudget();
      repeatForm.addEventListener("submit",event=>{
        event.preventDefault();if(!repeatForm.reportValidity())return;
        const next=readControls();
        const config=DecorationModel.make(seed,"repeated",next),url=new URL(location.href);
        if(!config.numericalBudget.supported){updateBudget();return;}
        url.searchParams.set("decoration","repeated");for(const [key,value]of Object.entries(config.counts))url.searchParams.set(key,value);
        document.getElementById("repeat-feedback").textContent="Applying the labelled decoration…";location.href=url.href;
      });
      document.getElementById("repeat-reset").addEventListener("click",()=>{
        const zero=DecorationModel.preset(seed,"seed");for(const key of keys)document.getElementById("repeat-"+key).value=zero[key];
        if(seed==="double")for(const end of ["A","B"])updateMinimum(end);
        updateBudget();document.getElementById("repeat-feedback").textContent="Counts reset. Apply to update the graphs.";
      });
    }
    document.documentElement.dataset.seed=seed;
  }
  mount();
  return {isNew:seed!=="lame",active,project,map,svgPoint,svgInverse,renderMonodromy,asymptotic,
    drawDomain:(canvas,options)=>window.AtlasView.drawDomain(canvas,{...options,markedStrip:true})};
})();

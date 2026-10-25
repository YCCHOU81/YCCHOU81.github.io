"use strict";

/* A labelled topological drawing, not a conformal or spherical metric model.
   The two drawings are built from the same exact ribbon map. Inserting two
   dual vertices along a stable edge inserts two parallel primal edges. */
globalThis.GraphGeometry = (() => {
  const W=720,H=500,ink="#394444",paper="#fbfbf8",regionFill="#f1f2ee";
  const sub=(a,b)=>a.map((x,i)=>x-b[i]),add=(a,b)=>a.map((x,i)=>x+b[i]),mul=(a,s)=>a.map(x=>x*s),lerp=(a,b,t)=>a.map((x,i)=>x*(1-t)+b[i]*t);
  const fmt=n=>Number(n.toFixed(3));
  const xy=p=>p.map(fmt).join(",");
  const esc=x=>String(x).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const point=(c,t)=>{const u=1-t;return c[0].map((_,k)=>u*u*u*c[0][k]+3*u*u*t*c[1][k]+3*u*t*t*c[2][k]+t*t*t*c[3][k]);};
  const path=c=>`M${xy(c[0])} C${xy(c[1])} ${xy(c[2])} ${xy(c[3])}`;
  const reverse=c=>c.slice().reverse();
  function split(c,t){const a=lerp(c[0],c[1],t),b=lerp(c[1],c[2],t),d=lerp(c[2],c[3],t),e=lerp(a,b,t),f=lerp(b,d,t),g=lerp(e,f,t);return[[c[0],a,e,g],[g,f,d,c[3]]];}
  function interval(c,a,b){if(a===0)return split(c,b)[0];if(b===1)return split(c,a)[1];return split(split(c,b)[0],a/b)[1];}
  function straight(a,b,bend=0){const d=sub(b,a),n=Math.hypot(...d),v=[-d[1]/n,d[0]/n];return[a,add(lerp(a,b,1/3),mul(v,bend)),add(lerp(a,b,2/3),mul(v,bend)),b];}
  function samples(c,n=32){return Array.from({length:n+1},(_,i)=>point(c,i/n));}
  function signedArea(poly){let out=0;for(let i=0;i<poly.length;i++){const a=poly[i],b=poly[(i+1)%poly.length];out+=a[0]*b[1]-b[0]*a[1];}return out/2;}
  function inside(p,poly){let yes=false;for(let i=0,j=poly.length-1;i<poly.length;j=i++){const a=poly[i],b=poly[j];if((a[1]>p[1])!==(b[1]>p[1])&&p[0]<(b[0]-a[0])*(p[1]-a[1])/(b[1]-a[1])+a[0])yes=!yes;}return yes;}
  const outer="M18,18H702V482H18Z";
  function templates(seed){
    if(seed==="bouquet"){
      const hp=[[360,270]],pp=[[360,245],[360,55],[635,405],[85,405]],dual={};
      ["l1+","l2+","l3+"].forEach((id,i)=>{const a=hp[0],angle=(-90+i*120)*Math.PI/180,u=[Math.cos(angle),Math.sin(angle)],v=[-u[1],u[0]];dual[id]=[a,add(add(a,mul(u,230)),mul(v,-100)),add(add(a,mul(u,230)),mul(v,100)),a];});
      return{pp,hp,dual,primal:{},primalExterior:0,dualExterior:0,poleLabels:[[330,254],[360,30],[666,421],[54,421]],hLabels:[[615,75]]};
    }
    if(seed==="double"){
      const hp=[[200,250],[520,250]],pp=[[150,240],[480,240],[350,240],[625,370]];
      return{pp,hp,dual:{t:[hp[0],[280,65],[440,65],hp[1]],u:[hp[0],[280,435],[440,435],hp[1]],"l+":[hp[0],[460,100],[460,400],hp[0]],"m+":[hp[1],[710,40],[710,460],hp[1]]},
        primal:{t:straight(pp[1],pp[0],165),u:straight(pp[0],pp[1],165),"l+":straight(pp[0],pp[2]),"m+":straight(pp[1],pp[3])},
        primalExterior:1,dualExterior:1,poleLabels:[[122,242],[505,229],[350,216],[654,388]],hLabels:[[420,240],[620,65]]};
    }
    if(seed==="tetra"){
      const hp=[[360,65],[90,420],[630,420],[360,300]],pp=[[70,420],[650,420],[360,45],[360,275]],dual={};
      const idx={A:0,B:1,C:2,D:3};for(const id of["AB","AC","AD","BC","BD","CD"])dual[id]=straight(hp[idx[id[0]]],hp[idx[id[1]]]);
      return{pp,hp,dual,primal:{},primalExterior:0,dualExterior:0,poleLabels:[[48,443],[674,443],[360,22],[388,277]],hLabels:[[640,65],[363,369],[250,239],[470,239]]};
    }
    throw new Error("The graph lab supports the three TV base seeds.");
  }
  function concatenate(darts,dartGeometry){let d="",poly=[];for(const dart of darts){const c=dartGeometry.get(dart);if(!c)throw new Error("Missing edge geometry");d+=(d?` C${xy(c[1])} ${xy(c[2])} ${xy(c[3])}`:path(c));poly.push(...samples(c).slice(0,-1));}return{path:d+"Z",poly,area:signedArea(poly)};}
  function basePrimal(graph,t,faceOf){const map=new Map();for(const e of graph.stable_edges){const a=t.pp[faceOf.get(e.start)],b=t.pp[faceOf.get(e.end)];let c=t.primal[e.id]||straight(a,b);if(Math.hypot(...sub(c[0],a))>1)c=reverse(c);map.set(e.start,c);map.set(e.end,reverse(c));}return map;}
  function label(content,p,attrs="",size=15){return`<text x="${fmt(p[0])}" y="${fmt(p[1])}" text-anchor="middle" dominant-baseline="middle" class="graph-label" font-size="${size}" ${attrs}>${esc(content)}</text>`;}
  function hitAttrs(kind,id,description){return`data-${kind}="${id}" tabindex="0" role="button" aria-label="${esc(description)}"`;}
  function build(config){
    const graph=config.graph,seed=config.seed,t=templates(seed);
    if((seed==="bouquet"&&config.core!=="bouquet_B0_T0")||(seed==="double"&&config.core!=="double_none_none")||graph.stable_edges.some(e=>e.epsilon!==1||e.count>3))throw new Error("Graph lab scope: base cores with 0–3 pair insertions per stable edge.");
    const alpha=new Map(graph.pairs.flatMap(([a,b])=>[[a,b],[b,a]])),vertexOf=new Map(graph.vertices.flatMap((ds,i)=>ds.map(d=>[d,i]))),faceOf=new Map(graph.faces.flatMap((ds,i)=>ds.map(d=>[d,i])));
    const metadata={seed,hemispheres:graph.vertices.map((v,i)=>({id:`h${i}`,index:i,label:`H${i+1}`,degree:v.length,exterior:i===t.primalExterior})),poles:graph.faces.map((f,i)=>({id:`p${i}`,index:i,label:`p${i}`,degree:f.length,weight:config.weights[i]})),edges:[],chains:[],counts:{hemispheres:graph.vertices.length,edges:graph.pairs.length,poles:graph.faces.length,euler:graph.vertices.length-graph.pairs.length+graph.faces.length}};
    const primal=new Map(),dual=new Map(),hPositions=t.hp.map(p=>p.slice()),hPrimal=t.hLabels.map(p=>p.slice()),hCandidates=new Map(),edgeRows=[];
    // Core face boundaries determine a common orientation. The extra strips
    // must expand toward the start hemisphere in the labelled cyclic order.
    const coreMap=basePrimal(graph,t,faceOf),coreVertices=graph.vertices.slice(0,t.hp.length),coreFaces=coreVertices.map(ds=>concatenate(ds,coreMap));
    const orient=Math.sign(coreFaces.filter((_,i)=>i!==t.primalExterior).reduce((s,f)=>s+f.area,0))||1;
    let cursor=0;
    graph.stable_edges.forEach((chain,ci)=>{
      const m=chain.multiplicity,q=chain.count||0,chainId=`chain${ci}`,color=ink,a=t.pp[faceOf.get(chain.start)],b=t.pp[faceOf.get(chain.end)],base=coreMap.get(chain.start),delta=sub(b,a),length=Math.hypot(...delta),normal=[-delta[1]/length,delta[0]/length],spread=seed==="bouquet"?80:seed==="double"?(["l+","m+"].includes(chain.id)?60:50):({AB:45,AC:35,AD:35,BC:20,BD:20,CD:30}[chain.id]),dBase=t.dual[chain.id];
      if(!dBase)throw new Error(`Missing seed edge ${chain.id}`);
      const row={id:chainId,index:ci,label:chain.id.replace("l1+","ℓ1").replace("l2+","ℓ2").replace("l3+","ℓ3").replace("l+","ℓ").replace("m+","m"),count:q,multiplicity:m,edges:[],color};
      metadata.chains.push(row);
      const startIndex=vertexOf.get(chain.start),endIndex=vertexOf.get(chain.end);let dCurve=dBase;
      if(Math.hypot(...sub(dCurve[0],t.hp[startIndex]))>1)dCurve=reverse(dCurve);
      for(let j=0;j<m;j++){
        const ei=cursor++,pair=graph.pairs[ei],[d,D]=pair,id=`e${ei}`,offset=m===1?0:orient*spread*(1-2*j/(m-1));
        const pCurve=base.map((p,k)=>k===0||k===3?p:add(p,mul(normal,offset))),dc=interval(dCurve,j/m,(j+1)/m);
        primal.set(d,pCurve);primal.set(D,reverse(pCurve));dual.set(d,dc);dual.set(D,reverse(dc));
        const edge={id,index:ei,chain:chainId,chainIndex:ci,label:`E${ei+1}`,poles:[faceOf.get(d),faceOf.get(D)],hemispheres:[vertexOf.get(d),vertexOf.get(D)],color,darts:[d,D]};metadata.edges.push(edge);row.edges.push(id);edgeRows.push({edge,pCurve,dCurve:dc});
        hPositions[vertexOf.get(d)]=dc[0];hPositions[vertexOf.get(D)]=dc[3];
        if(j<m-1){const nextOffset=orient*spread*(1-2*(j+1)/(m-1)),mean=(offset+nextOffset)/2,tt=.2+.6*(j+.5)/(m-1),mid=base.map((p,k)=>k===0||k===3?p:add(p,mul(normal,mean))),hi=vertexOf.get(D);hPrimal[hi]=point(mid,tt);hCandidates.set(hi,Array.from({length:61},(_,k)=>({t:.18+.64*k/60,p:point(mid,.18+.64*k/60)})).sort((a,b)=>Math.abs(a.t-tt)-Math.abs(b.t-tt)).map(x=>x.p));}
      }
    });
    if(cursor!==graph.pairs.length)throw new Error("Edge partition mismatch");
    const primalFaces=graph.vertices.map(ds=>concatenate(ds,primal)),dualFaces=graph.faces.map(ds=>concatenate(ds,dual));
    // Labels may slide along their own bigon, but never move into a different
    // hemisphere. Keep their estimated text boxes apart at the diagram scale.
    const box=(p,text,size)=>({x:p[0]-text.length*size*.34-3,y:p[1]-size*.55-2,w:text.length*size*.68+6,h:size*1.1+4});
    const overlap=(a,b)=>a.x<b.x+b.w&&b.x<a.x+a.w&&a.y<b.y+b.h&&b.y<a.y+a.h;
    const occupied=[...t.poleLabels.map((p,i)=>box(p,`p${i}`,17)),...metadata.hemispheres.filter(h=>h.index<t.hp.length).map(h=>box(hPrimal[h.index],h.label,17))];
    let labelCollisions=0;
    for(const h of metadata.hemispheres)if(hCandidates.has(h.index)){
      const candidates=hCandidates.get(h.index).filter(p=>inside(p,primalFaces[h.index].poly));
      const chosen=candidates.find(p=>!occupied.some(b=>overlap(box(p,h.label,12),b)))||candidates[0]||hPrimal[h.index];
      if(occupied.some(b=>overlap(box(chosen,h.label,12),b)))labelCollisions++;
      hPrimal[h.index]=chosen;occupied.push(box(chosen,h.label,12));
    }
    const polePositions=dualFaces.map((f,i)=>{
      if(i===t.dualExterior)return[66,65];
      const seedPoint=seed==="bouquet"?lerp(t.hp[0],point(t.dual[`l${i}+`],.5),.62):seed==="double"?(i===0?[455,250]:i===2?[295,250]:[625,250]):i===1?[270,273]:i===2?[450,273]:[360,379];
      return seedPoint;
    });
    // Use one stable-chain label in each picture. Its location is chosen among
    // nearby points on the actual chain, keeping the existing object labels and
    // vertices clear. A subdivided chain keeps one name rather than a name on
    // every inserted edge.
    const chainLabels={primal:[],dual:[]},chainLabelCollisions={primal:0,dual:0};
    for(const mode of ["primal","dual"]){
      const original=mode==="primal",nodes=original?t.pp:hPositions;
      const used=original?[...occupied]:[
        ...metadata.hemispheres.map(h=>box([hPositions[h.index][0],hPositions[h.index][1]-16],h.label,h.index<t.hp.length?15:12)),
        ...metadata.poles.map(p=>box(polePositions[p.index],`F${p.index}`,17))
      ];
      for(const chain of metadata.chains){
        const text=chain.label+(chain.multiplicity>1?` ×${chain.multiplicity}`:""),rows=edgeRows.filter(r=>r.edge.chain===chain.id);
        const curves=original?[rows[0].pCurve,rows[rows.length-1].pCurve]:[t.dual[graph.stable_edges[chain.index].id]];
        const candidates=[];
        for(const c of curves)for(const tt of [.5,.38,.62,.26,.74]){
          const p=point(c,tt),tangent=sub(point(c,Math.min(.999,tt+.002)),point(c,Math.max(.001,tt-.002))),len=Math.hypot(...tangent),normal=[-tangent[1]/len,tangent[0]/len];
          for(const off of [-18,18,-27,27,-36,36,-48,48,-60,60]){
            const pos=add(p,mul(normal,off)),bb=box(pos,text,14),collisions=used.reduce((n,b)=>n+Number(overlap(bb,b)),0),nearNode=nodes.reduce((n,node)=>n+Number(Math.hypot(...sub(pos,node))<29),0);
            if(bb.x<28||bb.x+bb.w>W-28||bb.y<28||bb.y+bb.h>H-28)continue;
            candidates.push({pos,box:bb,score:collisions*10000+nearNode*3000+Math.abs(tt-.5)*20+Math.abs(off)});
          }
        }
        candidates.sort((a,b)=>a.score-b.score);
        const chosen=candidates[0];
        if(!chosen)throw new Error("No space for a stable-chain label");
        if(used.some(b=>overlap(chosen.box,b)))chainLabelCollisions[mode]++;
        chainLabels[mode].push({chain:chain.id,edge:chain.edges[0],label:text,position:chosen.pos});used.push(chosen.box);
      }
    }
    function svg(mode){
      const original=mode==="primal",regions=original?primalFaces:dualFaces,ext=original?t.primalExterior:t.dualExterior;
      const title=original?"Hemisphere boundary graph":"Dual hemisphere graph";
      let out=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="correspondence-svg ${mode}-svg" role="group" aria-label="${title}"><rect x="18" y="18" width="684" height="464" rx="12" fill="${paper}" stroke="#d9ddda"/>`;
      const regionIndices=[ext,...regions.map((_,i)=>i).filter(i=>i!==ext)];
      for(const i of regionIndices){const face=regions[i],kind=original?"hemisphere":"pole",id=(original?"h":"p")+i,description=original?`Hemisphere H${i+1}; dual vertex H${i+1}`:`Dual face F${i}; original marked point p${i}`;
        out+=`<path d="${i===ext?outer:""}${face.path}" fill="${regionFill}" fill-opacity=".52" fill-rule="evenodd" stroke="none" class="graph-region ${original?"hemisphere-region":"pole-region"}" ${hitAttrs(kind,id,description)}><title>${description}${i===ext?" (exterior region; the drawing is on a sphere)":""}</title></path>`;
      }
      for(const {edge,pCurve,dCurve}of edgeRows){const c=original?pCurve:dCurve,chain=metadata.chains[edge.chainIndex],desc=`${edge.label}; chain ${chain.label}; original p${edge.poles[0]}–p${edge.poles[1]}, dual H${edge.hemispheres[0]+1}–H${edge.hemispheres[1]+1}`;
        out+=`<g class="graph-edge" data-chain="${edge.chain}" ${hitAttrs("edge",edge.id,desc)}><title>${esc(desc)}</title><path d="${path(c)}" class="edge-hit" fill="none" stroke="transparent" stroke-width="15"/><path d="${path(c)}" class="edge-ink edge-line" fill="none" stroke="${edge.color}" stroke-width="2.8" stroke-linecap="round"/></g>`;
      }
      if(original){
        for(const h of metadata.hemispheres){const p=hPrimal[h.index];out+=`<g class="hemisphere-label" ${hitAttrs("hemisphere",h.id,`Hemisphere ${h.label}; dual vertex ${h.label}`)}>${label(h.label,p,"",h.index<t.hp.length?17:12)}</g>`;}
        for(const pole of metadata.poles){const p=t.pp[pole.index];out+=`<g class="graph-pole" ${hitAttrs("pole",pole.id,`Marked point ${pole.label}; dual face F${pole.index}`)}><circle class="node-circle" cx="${p[0]}" cy="${p[1]}" r="7.5" fill="white" stroke="${ink}" stroke-width="2"/>${label(pole.label,t.poleLabels[pole.index],`fill="${ink}"`,17)}</g>`;}
      }else{
        for(const p of metadata.poles)out+=`<g class="dual-face-label" ${hitAttrs("pole",p.id,`Dual face F${p.index}; original marked point ${p.label}`)}>${label(`F${p.index}`,polePositions[p.index],`fill="${ink}"`,17)}</g>`;
        for(const h of metadata.hemispheres){const p=hPositions[h.index],r=h.index<t.hp.length?8:5;out+=`<g class="graph-hemisphere" ${hitAttrs("hemisphere",h.id,`Dual vertex ${h.label}; original hemisphere ${h.label}`)}><circle class="node-circle" cx="${fmt(p[0])}" cy="${fmt(p[1])}" r="${r}" fill="white" stroke="${ink}" stroke-width="1.6"/>${label(h.label,[p[0],p[1]-16],"",h.index<t.hp.length?15:12)}</g>`;}
      }
      for(const item of chainLabels[mode])out+=`<g class="graph-chain-label" data-chain="${item.chain}" ${hitAttrs("edge",item.edge,`Decoration chain ${item.label}`)}><title>Decoration chain ${esc(item.label)}</title>${label(item.label,item.position,`fill="${ink}" paint-order="stroke" stroke="${paper}" stroke-width="5" stroke-linejoin="round"`,14)}</g>`;
      out+="</svg>";return out;
    }
    metadata.positions={primalPoles:t.pp,primalHemispheres:hPrimal,dualHemispheres:hPositions,dualPoles:polePositions,chainLabels};
    metadata.layout={labelCollisions,chainLabelCollisions};
    metadata.geometry={primalFaces:primalFaces.map(f=>({area:f.area,polygon:f.poly})),dualFaces:dualFaces.map(f=>({area:f.area,polygon:f.poly})),edges:edgeRows.map(({edge,pCurve,dCurve})=>({id:edge.id,primal:pCurve,dual:dCurve}))};
    return{primalSvg:svg("primal"),dualSvg:svg("dual"),metadata};
  }
  return{build,point,samples,inside,signedArea};
})();

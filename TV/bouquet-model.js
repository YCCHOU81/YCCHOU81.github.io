"use strict";

/* Exact canonical three-loop decorations.  The graph, including its cyclic
   dart order and face labels, is retained to select the analytic branch. */
globalThis.BouquetModel = (() => {
  const cores=globalThis.BouquetData;
  if(!cores)throw new Error("The three-loop core registry was not loaded");
  const defaultCore="bouquet_B0_T0";
  const countKeys=Array.from({length:6},(_,i)=>`bq${i}`);
  const queryKeys=["bouquetCore",...countKeys];
  const numeric=x=>typeof x==="string"&&x.includes("/")?x.split("/").map(Number).reduce((a,b)=>a/b):Number(x);
  const integer=(x,fallback=0)=>Number.isSafeInteger(Number(x))&&Number(x)>=0?Math.min(32,Number(x)):fallback;
  const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);
  function rational(a,b){const d=gcd(a,b)||1;a/=d;b/=d;if(b<0){a=-a;b=-b;}return {value:a/b,label:b===1?String(a):`${a}/${b}`,numerator:a,denominator:b};}
  function preset(operation="seed") {
    const out={bouquetCore:defaultCore,...Object.fromEntries(countKeys.map(k=>[k,0]))};
    if(operation==="pair")out.bq0=1;
    if(operation==="sprout")out.bouquetCore="bouquet_B1_T0";
    if(operation==="split")out.bouquetCore="bouquet_B0_T03";
    return out;
  }
  function parse(operation,query){
    if(operation&&typeof operation.get==="function"){query=operation;operation="repeated";}
    const out=preset(operation);
    if(operation!=="repeated")return out;
    if(query?.get("bouquetCore")&&cores[query.get("bouquetCore")])out.bouquetCore=query.get("bouquetCore");
    countKeys.forEach(k=>out[k]=integer(query?.get(k)));
    return out;
  }
  function faceCycles(graph){
    const alpha=new Map(graph.pairs.flatMap(([a,b])=>[[a,b],[b,a]]));
    const sigma=new Map(graph.vertices.flatMap(v=>v.map((d,i)=>[d,v[(i+1)%v.length]])));
    const faces=[],seen=new Set();
    for(const [a,b]of graph.pairs)for(const first of[a,b])if(!seen.has(first)){
      const f=[];let d=first;
      do{if(seen.has(d)||!sigma.has(alpha.get(d)))throw new Error("Invalid three-loop ribbon face");seen.add(d);f.push(d);d=sigma.get(alpha.get(d));}while(d!==first);
      faces.push(f);
    }
    if(faces.length!==4)throw new Error("The three-loop graph must have four labelled faces");
    return graph.face_representatives.map(d=>faces.find(f=>f.includes(d)));
  }
  function expand(core,edgeCounts,id){
    const weights=core.weights.slice(),eta=core.shift_twice.slice();
    const alpha=new Map(core.pairs.flatMap(([a,b])=>[[a,b],[b,a]]));
    const vertices=core.vertices.filter(v=>v.length!==2).map(v=>v.slice());
    const stableDarts=new Set(vertices.flat());
    const graph={vertices,pairs:[],lengths:{},length_coefficients:{},
      face_representatives:core.faces.map(f=>f.find(d=>stableDarts.has(d))),
      P:core.P.map(numeric),name:id,seed:"bouquet",source:"classical",core_id:core.name,
      weights,shift:[],faces:[],stable_edges:[],counts:{...edgeCounts},
      sprouts:core.sprouts.slice(),diagonals:core.diagonals.map(v=>v.slice()),external_darts:core.external_darts.slice()};
    core.stable_edges.forEach((e,ei)=>{
      const q=edgeCounts[e.id]||0,m=e.epsilon+2*q;
      weights.forEach((_,i)=>weights[i]+=q*e.weight_per_pair[i]);
      eta.forEach((_,i)=>eta[i]=(eta[i]+q*e.shift_twice_per_pair[i])%2);
      const start=e.chain[0],end=alpha.get(e.chain[e.chain.length-1]);
      const ds=[start,...Array.from({length:m-1},(_,j)=>`qe${ei}_${j+1}p`)];
      const Ds=[...Array.from({length:m-1},(_,j)=>`qe${ei}_${j}m`),end];
      const first=e.first_length.map(numeric);
      for(let j=0;j<m;j++){
        const d=ds[j],D=Ds[j],v=j%2?graph.P.map((p,i)=>p-first[i]):first.slice();
        graph.pairs.push([d,D]);
        graph.lengths[d]=v.map((x,i)=>x?`(${x})*${["x","y","z"][i]}`:null).filter(Boolean).join("+")||"0";
        graph.length_coefficients[d]=v;graph.length_coefficients[D]=v.slice();
        if(j<m-1)graph.vertices.push([D,ds[j+1]]);
      }
      graph.stable_edges.push({...e,count:q,q,multiplicity:m,start,end});
    });
    graph.shift=eta.map(n=>n/2);graph.shift_twice=eta;graph.faces=faceCycles(graph);
    if(graph.faces.some((f,i)=>!f||(f.length-1)/2!==weights[i]))throw new Error("Three-loop face-degree check failed");
    if(vertices.length!==weights.reduce((a,b)=>a+b,0))throw new Error("Three-loop hemisphere count failed");
    return graph;
  }
  function spectrum(a,b){
    const out=[],d=Math.abs(a-b);
    for(let h=1;h<=d;h++)out.push(h);
    for(let h=d+2;h<=a+b;h+=2)out.push(h);
    return out;
  }
  const bigGcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){const t=a%b;a=b;b=t;}return a;};
  function residue(weights,edge,M){
    const i=edge+1,other=[1,2,3].filter(j=>j!==i);
    const A=spectrum(weights[0],weights[i]),B=spectrum(weights[other[0]],weights[other[1]]);
    if(!A.includes(M)||B.includes(M))throw new Error("Invalid three-loop cusp resonance");
    let numerator=BigInt(2*M),denominator=1n,signL=1,logK=Math.log(2*M);
    for(const h of[...A,...B])if(h!==M){
      const aa=M+h,bb=M-h;logK+=Math.log(aa)-Math.log(Math.abs(bb));
      if(bb<0)signL*=-1;
      numerator*=BigInt(aa);denominator*=BigInt(Math.abs(bb));
      const d=bigGcd(numerator,denominator);numerator/=d;denominator/=d;
    }
    const absLabel=denominator===1n?String(numerator):`${numerator}/${denominator}`;
    return {signL,logK,K:Math.exp(logK),L:(signL<0?"-":"")+absLabel,KLabel:absLabel,
      residueNumerator:(signL<0?"-":"")+String(numerator),residueDenominator:String(denominator)};
  }
  function cuspData(core,edgeCounts,weights,shift){
    return core.cusp_affine.map((row,edge)=>{
      const q=core.stable_edges.map(e=>edgeCounts[e.id]||0),affine=a=>a[0]+q.reduce((s,n,i)=>s+n*a[i+1],0);
      const M=affine(row.M_affine),p=affine(row.p_affine),x=rational(p,2*M),res=residue(weights,edge,M);
      const eta=[shift[1],shift[0]+shift[1],shift[0]][edge],sigma=Math.round(2*eta)%2?-1:1;
      if(((p%4)+4)%4!==(sigma*res.signL===1?1:3))throw new Error("Three-loop marked phase congruence failed");
      return {edge,M,p,x0:x.value,xLabel:x.label,rational:x.label,numerator:x.numerator,denominator:x.denominator,
        xNumerator:x.numerator,xDenominator:x.denominator,alpha:1/M,eta,sigma,...res};
    });
  }
  let baseConfig=null;
  function make(operation="repeated",counts={}){
    if(operation&&typeof operation==="object"){counts=operation;operation="repeated";}
    const given=operation==="repeated"?counts:preset(operation);
    const core=cores[given.bouquetCore]||cores[defaultCore];
    const normal={bouquetCore:core.name,...Object.fromEntries(countKeys.map((k,i)=>[k,i<core.stable_edges.length?integer(given[k]):0]))};
    const edgeCounts=Object.fromEntries(core.stable_edges.map((e,i)=>[e.id,normal[`bq${i}`]]));
    const id=operation==="repeated"?`bouquet-repeat-${core.name}-${core.stable_edges.map((e,i)=>normal[`bq${i}`]).join(",")}`:`bouquet-${operation}`;
    const graph=expand(core,edgeCounts,id);
    const config={id,seed:"bouquet",operation,weights:graph.weights.slice(),shift:graph.shift.slice(),top:.5,
      counts:normal,core:core.name,edgeCounts,graph,hemispheres:graph.vertices.length,
      sourceVertices:[[.5,0],[0,.5],[.5,.5]],start:[1/3,1/3],
      numericalBudget:{totalLimit:144,individualLimit:72,supported:graph.vertices.length<=144&&Math.max(...graph.weights)<=72}};
    config.cusps=cuspData(core,edgeCounts,config.weights,config.shift);
    if(baseConfig)config.baseGraph=baseConfig.graph;
    return config;
  }
  function asymptotic(config,edge,R,S,base=false){
    if(![0,1,2].includes(edge))return null;
    const c=(base?baseConfig:config).cusps[edge],t=edge===0?S:R;
    const sine=Math.sin(2*Math.PI*t),logC=c.logK+Math.log(Math.max(0,sine));
    return {...c,known:true,proved:true,cusp:["∞","0","1"][edge],chart:["raw","zero","one"][edge],
      chartLabel:["τ","w = 1 − 1/τ","w = 1/(1 − τ)"][edge],
      C:Math.exp(logC),logC,constantFactor:c.K,constantLabel:c.KLabel,
      deltaLabel:["½ − R","½ − S","R + S − ½"][edge],sineLabel:edge===0?"sin(2πS)":"sin(2πR)",
      detail:"Proved for the canonical labelled graph and fixed finite counts, away from source vertices. The plotted points are independently refined roots; the asymptotic formula is not substituted for their computation."};
  }
  const coreList=Object.values(cores);
  coreList.forEach(c=>{
    c.label=`Sprouts ${c.sprouts.length?c.sprouts.join(", "):"none"} · ${c.diagonals.length?"splits "+c.diagonals.map(e=>e.join("–")).join(", "):"unsplit centre"}`;
  });
  baseConfig=make("seed");baseConfig.baseGraph=baseConfig.graph;
  return {cores,coreList,queryKeys,countKeys,preset,parse,make,asymptotic,faceCycles,spectrum,residue,
    defaultCore,baseConfig,baseGraph:baseConfig.graph,sourceVertices:[[.5,0],[0,.5],[.5,.5]]};
})();

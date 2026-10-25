"use strict";

/* Conjectural count formulas for annotations only. The numerical workers do
   not import this file; no plotted point or solver initialization is changed. */
window.TVCuspFormulas = (() => {
  const esc=s=>String(s).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
  const n=x=>`<mn>${x}</mn>`,v=x=>`<mi>${x}</mi>`,op=x=>`<mo>${x}</mo>`;
  const row=(...a)=>`<mrow>${a.join("")}</mrow>`,frac=(a,b)=>`<mfrac>${a}${b}</mfrac>`;
  const par=a=>row(op("("),a,op(")")),sub=(a,b)=>`<msub>${v(a)}${n(b)}</msub>`;
  const math=a=>`<math xmlns="http://www.w3.org/1998/Math/MathML">${a}</math>`;
  function rational(a,b) {
    let x=Math.abs(a),y=b;while(y){const t=x%y;x=y;y=t;}a/=x||1;b/=x||1;
    return {numerator:a,denominator:b,value:a/b,label:b===1?String(a):`${a}/${b}`};
  }
  function choose(a,b){let x=1;for(let j=1;j<=b;j++)x=x*(a-b+j)/j;return x;}
  const zero=(c,keys)=>keys.every(k=>Number(c[k]??0)===0);
  const chart=(seed,i)=>["raw","zero",seed==="double"?"half":"one"][i];
  const doubleCache=new Map(),tetraCache=new Map(),bouquetCache=new Map();
  const gcdBig=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){const t=a%b;a=b;b=t;}return a||1n;};
  function exactFraction(a,b=1n) {
    if(b===0n)throw new Error("A cusp product has an unexpected second pole");
    if(b<0n){a=-a;b=-b;}const g=gcdBig(a,b);return {a:a/g,b:b/g};
  }
  const times=(x,y)=>exactFraction(x.a*y.a,x.b*y.b);
  const absolute=x=>({a:x.a<0n?-x.a:x.a,b:x.b});
  const exactLabel=x=>x.b===1n?String(x.a):`${x.a}/${x.b}`;
  const exactHTML=x=>x.b===1n?n(x.a):frac(n(x.a),n(x.b));
  function logBig(x){const s=String(x);return Math.log(Number(s.slice(0,16)))+(s.length-Math.min(16,s.length))*Math.LN10;}
  const logFraction=x=>logBig(x.a)-logBig(x.b);
  function sqrtBig(x) {
    if(x<2n)return x;
    let a=1n<<BigInt(Math.ceil(x.toString(2).length/2));
    for(;;){const b=(a+x/a)/2n;if(b>=a)return a;a=b;}
  }
  function constantData(x,squareRoot=false) {
    x=absolute(x);const logValue=logFraction(x)/(squareRoot?2:1);
    let exact=x,radical=false;
    if(squareRoot){const a=sqrtBig(x.a),b=sqrtBig(x.b);if(a*a===x.a&&b*b===x.b)exact={a,b};else radical=true;}
    const label=radical?`sqrt(${exactLabel(x)})`:exactLabel(exact);
    return {value:Math.exp(logValue),logValue,label,html:radical?`<msqrt>${exactHTML(x)}</msqrt>`:exactHTML(exact),
      numerator:String(exact.a),denominator:String(exact.b),squareRoot:radical,
      squaredLabel:squareRoot?exactLabel(x):null};
  }
  function capPoles(a,b) {
    const h=[],d=Math.abs(a-b);
    for(let j=1;j<=d;j++)h.push(j);
    for(let j=d+2;j<=a+b;j+=2)h.push(j);
    return h;
  }
  function capProduct(a,b,M,residue=false) {
    const poles=capPoles(a,b);
    if(residue&&!poles.includes(M))throw new Error("The proposed cusp index is not a cap pole");
    if(!residue&&poles.includes(M))throw new Error("The second cap is resonant at the proposed cusp index");
    let result=exactFraction(BigInt(residue?2*M:1));
    for(const h of poles)if(h!==M)result=times(result,exactFraction(BigInt(M+h),BigInt(M-h)));
    return result;
  }
  function doubleCounts(config) {
    const c=config.counts??{},out={};
    for(const key of ["kA","kB","cA","cB","qt","qu","ql","qm"]){
      const value=c[key]??0;if(typeof value==="boolean"||!Number.isSafeInteger(Number(value))||Number(value)<0)return null;
      out[key]=Number(value);
    }
    for(const side of ["A","B"]){
      if(Object.hasOwn(c,"mode"+side)||Object.hasOwn(c,"count"+side)){
        const mode=c["mode"+side]??"arm",value=c["count"+side]??0,count=Number(value);
        if(!["arm","cross"].includes(mode)||typeof value==="boolean"||!Number.isSafeInteger(count)||count<0||(mode==="cross"&&!count))return null;
        out["k"+side]=mode==="arm"?count:0;out["c"+side]=mode==="cross"?count:0;
      }
      if(out["k"+side]&&out["c"+side])return null;
    }
    return out;
  }
  function doubleData(config) {
    const counts=doubleCounts(config);if(!counts)return null;
    const key=JSON.stringify(counts);if(doubleCache.has(key))return doubleCache.get(key);
    const {kA,kB,cA,cB,qt:T,qu:U,ql:l,qm:m}=counts,K=kA+kB,C=cA+cB,Q=l+m;
    const weights=[1+kA+cB+T+U+l,1+kB+cA+T+U+m,cA+l,cB+m];
    // The website admits at most 144 hemispheres. This arithmetic guard also
    // keeps malformed, hand-written configurations from blocking the UI.
    if(weights.reduce((a,b)=>a+b,0)>10000)return null;
    const M=[2+K+C+Q+2*U,1+C+T+U,2+K+C+Q+2*T],p=[1+2*(U+C),T-U,-1-2*(T+C)];
    const x=[rational(p[0],2*M[0]),rational(p[1],M[1]),rational(p[2],2*M[2])];
    const alpha=[rational(1,M[0]),rational(2,M[1]),rational(1,M[2])];
    const [a,b,c,d]=weights;
    let constants;
    try {
      constants=[constantData(times(capProduct(a,b,M[0],true),capProduct(c,d,M[0]))),
        // These are two separate simple residues. Their product has a double
        // pole; taking the residue of that product would give the wrong law.
        constantData(times(capProduct(a,c,M[1],true),capProduct(d,b,M[1],true)),true),
        constantData(times(capProduct(a,b,M[2],true),capProduct(c,d,M[2])))];
    }catch{return null;}
    const result={counts,totals:{K,C,Q,T,U},weights,M,p,x,alpha,constants};
    doubleCache.set(key,result);if(doubleCache.size>32)doubleCache.delete(doubleCache.keys().next().value);
    return result;
  }
  function tetraData(config) {
    const source=config.counts??{},counts={};
    for(const key of ["qAB","qAC","qAD","qBC","qBD","qCD"]){
      const value=source[key]??0;
      if(typeof value==="boolean"||!Number.isSafeInteger(Number(value))||Number(value)<0)return null;
      counts[key]=Number(value);
    }
    const key=JSON.stringify(counts);if(tetraCache.has(key))return tetraCache.get(key);
    const {qAB:AB,qAC:AC,qAD:AD,qBC:BC,qBD:BD,qCD:CD}=counts;
    const A=AB+CD,B=AC+BD,C=AD+BC;
    const weights=[1+AB+AC+BC,1+AB+AD+BD,1+AC+AD+CD,1+BC+BD+CD];
    if(weights.reduce((a,b)=>a+b,0)>10000)return null;
    const M=[B+C+2,A+C+2,A+B+2],p=[C+1,A+1,B+1];
    const x=p.map((v,i)=>rational(v,M[i])),alpha=M.map(m=>rational(2,m));
    const [a,b,c,d]=weights;
    let constants;
    try{
      // Each factor is a separate SIMPLE residue. The product of the cap
      // functions has a double pole; its residue is not this leading scale.
      constants=[constantData(times(capProduct(a,b,M[0],true),capProduct(c,d,M[0],true)),true),
        constantData(times(capProduct(a,c,M[1],true),capProduct(d,b,M[1],true)),true),
        constantData(times(capProduct(a,d,M[2],true),capProduct(b,c,M[2],true)),true)];
    }catch{return null;}
    const result={counts,totals:{A,B,C},weights,M,p,x,alpha,constants};
    tetraCache.set(key,result);if(tetraCache.size>32)tetraCache.delete(tetraCache.keys().next().value);
    return result;
  }
  function bouquetData(config,core) {
    if(!core||!Array.isArray(config.weights)||config.weights.length!==4||
      config.weights.some(x=>!Number.isSafeInteger(x)||x<0)||config.weights.reduce((a,b)=>a+b,0)>10000)return null;
    const counts=core.stable_edges.map(e=>Number(config.edgeCounts?.[e.id]??config.graph?.counts?.[e.id]??0));
    if(counts.some(q=>!Number.isSafeInteger(q)||q<0))return null;
    const key=JSON.stringify([config.core,counts,config.weights]);
    if(bouquetCache.has(key))return bouquetCache.get(key);
    const affine=a=>a[0]+counts.reduce((sum,q,i)=>sum+q*a[i+1],0);
    const M=core.cusp_affine.map(a=>affine(a.M_affine)),p=core.cusp_affine.map(a=>affine(a.p_affine));
    if(M.some(m=>!Number.isSafeInteger(m)||m<1)||p.some(v=>!Number.isSafeInteger(v)))return null;
    let constants;
    try{constants=M.map((m,edge)=>{
      const i=edge+1,[j,k]=[1,2,3].filter(v=>v!==i);
      return constantData(times(capProduct(config.weights[0],config.weights[i],m,true),capProduct(config.weights[j],config.weights[k],m)));
    });}catch{return null;}
    const result={counts,M,p,constants};bouquetCache.set(key,result);
    if(bouquetCache.size>32)bouquetCache.delete(bouquetCache.keys().next().value);
    return result;
  }
  function affine(a) {
    let out=a[0]!==0?n(a[0]):"";
    a.slice(1).forEach((coefficient,i)=>{
      if(!coefficient)return;
      const sign=coefficient<0?"−":"+",c=Math.abs(coefficient);
      out+=(out||coefficient<0?op(sign):"")+(c===1?"":n(c))+sub("q",i+1);
    });
    return row(out||n(0));
  }
  function lookup(config,edge) {
    if(!config||![0,1,2].includes(edge))return null;
    const c=config.counts??{},seed=config.seed;
    if(seed==="double") {
      const data=doubleData(config);if(!data)return null;
      const {K,C,Q,T,U}=data.totals,constant=data.constants[edge],a=data.alpha[edge];
      const sum=(...terms)=>row(...terms.flatMap((term,i)=>i?[op("+"),term]:[term]));
      const common=sum(n(2),v("K"),v("C"),v("Q"));
      const formulas=[frac(sum(n(1),row(n(2),par(sum(v("C"),v("U"))))),row(n(2),par(sum(common,row(n(2),v("U")))))),
        frac(row(v("T"),op("−"),v("U")),sum(n(1),v("C"),v("T"),v("U"))),
        frac(row(op("−"),par(sum(n(1),row(n(2),par(sum(v("C"),v("T"))))))),row(n(2),par(sum(common,row(n(2),v("T"))))))];
      return {...data.x[edge],family:"double-all-decorations",name:"Double-edge count conjecture",conjectural:true,
        chart:chart(seed,edge),counts:data.counts,totals:data.totals,weights:data.weights,M:data.M[edge],p:data.p[edge],
        allIndices:data.M.slice(),allPhases:data.p.slice(),alpha:a.value,alphaRational:a,alphaLabel:a.label,
        constantFactor:constant.value,logConstantFactor:constant.logValue,constantLabel:constant.label,
        constantHTML:math(constant.html),constantExact:{numerator:constant.numerator,denominator:constant.denominator,squareRoot:constant.squareRoot},
        constantSquaredLabel:constant.squaredLabel,
        formulaHTML:math(row(sub("x",edge),op("="),formulas[edge])),
        parameters:`K = kA + kB = ${K}; C = cA + cB = ${C}; Q = qℓ + qm = ${Q}; T = qt = ${T}; U = qu = ${U}. Arms: (kA,kB) = (${data.counts.kA},${data.counts.kB}); cross counts: (cA,cB) = (${data.counts.cA},${data.counts.cB}).`,
        evidence:"Supported by 22 checked source paths in 8 labelled decorated graphs. The law for arbitrary counts, including its leading constant, remains conjectural. These are full marked limits, without reduction modulo one; computed roots do not use these formulas.",
        anchor:"count-conjectures"};
    }
    if(seed==="tetra") {
      const data=tetraData(config);if(!data)return null;
      const {A,B,C}=data.totals,constant=data.constants[edge],a=data.alpha[edge];
      const formulas=[frac(row(v("C"),op("+"),n(1)),row(v("B"),op("+"),v("C"),op("+"),n(2))),
        frac(row(v("A"),op("+"),n(1)),row(v("A"),op("+"),v("C"),op("+"),n(2))),
        frac(row(v("B"),op("+"),n(1)),row(v("A"),op("+"),v("B"),op("+"),n(2)))];
      return {...data.x[edge],family:"tetra-all-decorations",name:"Tetrahedral count conjecture",conjectural:true,chart:chart(seed,edge),
        counts:data.counts,totals:data.totals,weights:data.weights,M:data.M[edge],p:data.p[edge],
        allIndices:data.M.slice(),allPhases:data.p.slice(),alpha:a.value,alphaRational:a,alphaLabel:a.label,
        constantFactor:constant.value,logConstantFactor:constant.logValue,constantLabel:constant.label,
        constantHTML:math(constant.html),constantExact:{numerator:constant.numerator,denominator:constant.denominator,squareRoot:constant.squareRoot},
        constantSquaredLabel:constant.squaredLabel,
        formulaHTML:math(row(sub("x",edge),op("="),formulas[edge])),
        parameters:`A = qAB + qCD = ${A}; B = qAC + qBD = ${B}; C = qAD + qBC = ${C}. These are totals on opposite edge pairs, using the displayed graph labels.`,
        evidence:"Supported by 14 deep and 3 shallow marked approaches across 9 decorated graphs, with 72 independent determinant-coefficient experiments. Opposite-edge totals determine the proposed rational limits and logarithmic rates; the positive leading constants retain the individual counts. The law for arbitrary counts remains conjectural; computed roots do not use these formulas.",anchor:"tetra-all-counts"};
    }
    if(seed==="bouquet") {
      const core=window.BouquetModel?.cores[config.core],a=core?.cusp_affine?.[edge],data=bouquetData(config,core);
      if(!a||!data)return null;
      const M=data.M[edge],p=data.p[edge],constant=data.constants[edge],alpha=rational(1,M);
      const edgeLabel=(e,i)=>/^l[123]/.test(e.id)?`ℓ${e.id[1]}`:/^b[123]/.test(e.id)?`b${e.id[1]}`:`h${core.stable_edges.slice(0,i+1).filter(x=>!/^l[123]|^b[123]/.test(x.id)).length}`;
      return {...rational(p,2*M),family:"bouquet-affine",name:"Three-loop cusp conjecture",conjectural:true,chart:chart(seed,edge),
        M,p,weights:config.weights.slice(),allIndices:data.M.slice(),allPhases:data.p.slice(),
        alpha:alpha.value,alphaRational:alpha,alphaLabel:alpha.label,
        constantFactor:constant.value,logConstantFactor:constant.logValue,constantLabel:constant.label,
        constantHTML:math(constant.html),constantExact:{numerator:constant.numerator,denominator:constant.denominator,squareRoot:false},
        deltaLabel:["1/2 − R","1/2 − S","R + S − 1/2"][edge],
        formulaHTML:math(row(sub("x",0),op("="),frac(sub("P",edge+1),row(n(2),sub("M",edge+1))),op(","),sub("M",edge+1),op("="),affine(a.M_affine),op(","),sub("P",edge+1),op("="),affine(a.p_affine))),
        parameters:`Core ${config.core.replace("bouquet_","")}. `+core.stable_edges.map((e,i)=>`q${i+1} = ${data.counts[i]} on ${edgeLabel(e,i)}`).join("; ")+`. Thus M${edge+1} = ${M}, P${edge+1} = ${p}. The full marked value is retained, including negative values and values greater than one.`,
        evidence:"Numerically checked on 12 source paths across 10 mixed graphs, including sprouts, nested splits, trapped components and counts (32,1,0). The center and cap-residue constant are conjectural comparison references; they do not supply the root guesses.",anchor:"three-loop-cusps"};
    }
    return null;
  }
  function html(config,edge,{compact=false}={}) {
    const f=lookup(config,edge);if(!f)return "";
    const currentSearch=typeof location!=="undefined"?location.search:window.location?.search;
    const formulaHref=esc(`mathematics.html${typeof currentSearch==="string"?currentSearch:""}#${f.anchor}`);
    const variable=f.chart==="raw"?"τ":"w",coordinate={raw:"τ",zero:"w = 1 − 1/τ",half:"w = τ/(1 − 2τ)",one:"w = 1/(1 − τ)"}[f.chart];
    const index=["double-all-decorations","tetra-all-decorations"].includes(f.family)?edge:0;
    const value=math(row(sub("x",index),op("="),f.denominator===1?n(f.numerator):frac(n(f.numerator),n(f.denominator))));
    const parameters=compact&&f.totals?`Totals (${f.family==="tetra-all-decorations"?"A,B,C":"K,C,Q,T,U"}) = (${Object.values(f.totals).join(", ")}).`:f.parameters;
    let leading="";
    if(!compact&&["double-all-decorations","tetra-all-decorations","bouquet-affine"].includes(f.family)){
      const a=f.alphaRational,alpha=math(row(sub("α",edge),op("="),a.denominator===1?n(a.numerator):frac(n(a.numerator),n(a.denominator))));
      const symbol=f.family==="double-all-decorations"?"𝒞":"K";
      const showExact=f.constantLabel.length<=(f.family==="tetra-all-decorations"?40:70);
      const constant=showExact?f.constantHTML:math(row(sub(symbol,edge),op("≈"),`<mtext>${esc(Number.isFinite(f.constantFactor)?f.constantFactor.toExponential(6):`exp(${f.logConstantFactor.toPrecision(8)})`)}</mtext>`));
      leading=`<p class="caption">Logarithmic slope: ${alpha}. Positive leading constant ${showExact?math(row(sub(symbol,edge),op("="))):""}${constant}. Its finite-product formula is given in Mathematics.</p>`;
      if(f.family==="bouquet-affine"){
        const normal=x=>`<mi mathvariant="normal">${x}</mi>`;
        const law=math(row(normal("Im"),v(variable),op("="),frac(n(1),row(v("π"),sub("M",edge+1))),normal("log"),par(frac(row(sub("K",edge+1),normal("sin"),par(row(v("π"),v("t")))),sub("δ",edge+1))),op("+"),normal("o"),par(n(1))));
        leading+=`<div class="tv-count-equation">${law}</div><p class="caption">Here δ = ${esc(f.deltaLabel)} is the actual boundary distance, and t runs from 0 to 1 along the selected edge. The endpoint coordinate is t/2, so the factor is sin(πt). K is the exact rational cap residue. Trapped hemispheres affect M; the complete marked turn affects P. Large M can make convergence slow even for very small ε.</p>`;
      }
    }
    return `<div class="tv-count-formula" data-formula-family="${f.family}"><strong>${esc(f.name)}</strong><p class="caption">Edge ${edge+1} · x${index===0?"₀":index===1?"₁":"₂"} = limiting Re ${variable}${f.chart==="raw"?"":`; ${coordinate}`}.</p><div class="tv-count-equation">${f.formulaHTML}</div><p>${esc(parameters)}</p><p class="tv-count-evaluated">With these counts: ${value}${compact?` <a href="${formulaHref}">Formula and checks →</a>`:""}</p>${leading}${compact?"":`<p class="caption">${esc(f.evidence)} <a href="${formulaHref}">Formula and checks →</a></p>`}</div>`;
  }
  return {lookup,html,rational};
})();

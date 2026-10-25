"use strict";

/* Actual hemisphere-chart holonomy.  No elliptic solver, numerical period
   fitting, or decoration sign table is used to calculate the multipliers. */
globalThis.GraphHolonomy = (() => {
  const add = (a,b) => a.map((x,i) => x+b[i]);
  const scale = (a,s) => a.map(x => x*s);
  const value = (a,R,S) => a[0]+a[1]*R+a[2]*S;
  const integer = x => {
    const n=Math.round(x);
    if (Math.abs(x-n)>1e-8 || !Number.isSafeInteger(n))
      throw Error("The marked holonomy coefficients are not safe integers.");
    return n;
  };
  const mod = (a,m) => ((a%m)+m)%m;
  const gcd = (a,b) => { a=Math.abs(a);b=Math.abs(b);while(b){const c=a%b;a=b;b=c;}return a; };
  function extended(a,b) {
    let x=1,y=0,X=0,Y=1;
    while(b){const q=Math.floor(a/b);[a,b]=[b,a-q*b];[x,X]=[X,x-q*X];[y,Y]=[Y,y-q*Y];}
    return {g:a,x,y};
  }
  function congruence(a,b,m) {
    if(m===1) return {r:0,m:1};
    a=mod(a,m);b=mod(b,m);
    const g=gcd(a,m);
    if(b%g) return null;
    const M=m/g;
    return {r:M===1?0:mod(extended(a/g,M).x*(b/g),M),m:M};
  }
  function combine(a,b) {
    if(!a||!b) return null;
    const correction=congruence(a.m,b.r-a.r,b.m);
    if(!correction) return null;
    const M=a.m*correction.m;
    return {r:mod(a.r+a.m*correction.r,M),m:M};
  }
  function periodWord(rows,target) {
    let best=null;
    for(let p=0;p<3;p++) for(let q=p+1;q<3;q++) {
      const r=3-p-q;
      const [ap,bp]=rows[p].slice(1),[aq,bq]=rows[q].slice(1),[ar,br]=rows[r].slice(1);
      const D=ap*bq-aq*bp;
      if(!D) continue;
      const [A,B]=target;
      const x0=A*bq-aq*B,x1=-ar*bq+aq*br;
      const y0=ap*B-A*bp,y1=-ap*br+ar*bp;
      const z=combine(congruence(x1,-x0,Math.abs(D)),congruence(y1,-y0,Math.abs(D)));
      if(!z) continue;
      const xb=integer((x0+x1*z.r)/D),yb=integer((y0+y1*z.r)/D);
      const xs=integer(x1*z.m/D),ys=integer(y1*z.m/D);
      const candidates=new Set([0]);
      // The L1 norm of the solution is convex in its free integer parameter.
      // Its minimum lies next to a zero of one of these three affine terms.
      for(const [v,s] of [[z.r,z.m],[xb,xs],[yb,ys]]) if(s) {
        const k=-v/s;
        for(const j of [Math.floor(k)-1,Math.floor(k),Math.ceil(k),Math.ceil(k)+1]) candidates.add(j);
      }
      for(const k of candidates) {
        const exponents=[0,0,0];
        exponents[p]=xb+k*xs;exponents[q]=yb+k*ys;exponents[r]=z.r+k*z.m;
        const cost=exponents.reduce((s,n)=>s+Math.abs(n),0);
        if(!best||cost<best.cost) best={exponents,cost};
      }
    }
    if(!best) throw Error("The graph characters do not span the marked torus periods.");
    const affine=rows.reduce((sum,row,i)=>add(sum,scale(row,best.exponents[i])),[0,0,0]);
    if(affine[1]!==target[0]||affine[2]!==target[1]) throw Error("Marked word reconstruction failed.");
    const terms=best.exponents.flatMap((power,i)=>power?[{i,power,gamma:[i,3]}]:[]);
    const text=terms.length?terms.map(t=>`(γ${t.i}γ3)${t.power===1?"":`^(${t.power})`}`).join(" · "):"I";
    return {exponents:best.exponents,terms,text,affinePi:affine};
  }
  function canonicalTree(config) {
    const g=config.graph,seed=config.seed;
    // Keep the based puncture paths fixed when the displayed base cores are
    // subdivided by hemisphere pairs. A fresh shortest-path tree can change
    // an individual developed Q_i by a period rotation, obscuring the fact
    // that pair insertion only changes ordinary lift signs.
    const ids=seed==="bouquet"&&config.core==="bouquet_B0_T0"
      ? ["l1+","l2+","l3+"]
      : seed==="double"&&config.core==="double_none_none"
        ? ["t","l+","u","m+"]
        : seed==="tetra"&&config.core==="tetra"
          ? ["AB","AD","AC","BC","BD","CD"] : null;
    if(!ids||g.stable_edges?.length!==ids.length||
      g.stable_edges.some((e,i)=>e.id!==ids[i]||e.chain.length!==1||
        !Number.isInteger(e.multiplicity)||e.multiplicity<1||e.multiplicity%2!==1)) return null;
    const baseTree=new Set(seed==="double"?["t"]:seed==="tetra"?["AB","AC","AD"]:[]);
    const tree=new Set();let offset=0;
    for(const edge of g.stable_edges) {
      const amount=edge.multiplicity-(baseTree.has(edge.id)?0:1);
      for(const pair of g.pairs.slice(offset,offset+amount)) tree.add(pair[0]);
      offset+=edge.multiplicity;
    }
    return offset===g.pairs.length&&tree.size===g.vertices.length-1?tree:null;
  }
  function compute(config,point) {
    const started=globalThis.performance?.now?.()??Date.now();
    const [R,S]=point,g=config.graph,seed=config.seed;
    if(!g||!["bouquet","double","tetra"].includes(seed)) throw Error("Unsupported metric ribbon graph.");
    const basis=seed==="bouquet"?[[1,-2,0],[1,0,-2],[-1,2,2]]:
      seed==="double"?[[2,-2,0],[-2,2,4],[1,0,-2]]:
      [[2,-2,0],[2,0,-2],[-2,2,2]];
    const lengthAffine=d=>basis.reduce((sum,row,i)=>add(sum,scale(row,g.length_coefficients[d][i])),[0,0,0]);
    const alpha={},vertex={},lengths={},perimeters=[];
    g.vertices.forEach((darts,v)=>{
      let sum=[0,0,0];
      for(const d of darts){alpha[d]=sum;vertex[d]=v;lengths[d]=lengthAffine(d);sum=add(sum,lengths[d]);}
      if(sum.some((x,i)=>Math.abs(x-(i===0?2:0))>1e-9)) throw Error("A hemisphere does not have perimeter 2π.");
      perimeters.push({index:v,darts:darts.slice(),length:2*Math.PI});
    });
    const tree=canonicalTree(config),adjacency=g.vertices.map(()=>[]),transitions=[];
    for(const [d,D] of g.pairs){
      const affine=add(add(alpha[d],alpha[D]),lengths[d]);
      const edge={d,mate:D,from:vertex[d],to:vertex[D],affinePi:affine,
        lengthPi:value(lengths[d],R,S),alphaPi:value(alpha[d],R,S),mateAlphaPi:value(alpha[D],R,S),cPi:value(affine,R,S)};
      edge.c=Math.PI*edge.cPi;edge.length=Math.PI*edge.lengthPi;
      edge.index=transitions.length;
      edge.inChartTree=tree?tree.has(d):false;
      transitions.push(edge);
      if(!tree||tree.has(d)) {
        adjacency[edge.from].push({to:edge.to,affine,d,index:edge.index});
        adjacency[edge.to].push({to:edge.from,affine,d:D,index:edge.index});
      }
    }
    // Any tree is valid.  The affine word reconstruction below recovers the
    // fixed marking instead of confusing a newly based loop with a period.
    const charts=[{sign:1,offset:[0,0,0],parent:null,viaDart:null,viaTransition:null}],queue=[0];
    for(const v of queue) for(const edge of adjacency[v]) if(!charts[edge.to]){
      charts[edge.to]={sign:-charts[v].sign,offset:add(charts[v].offset,scale(edge.affine,charts[v].sign)),
        parent:v,viaDart:edge.d,viaTransition:edge.index};
      transitions[edge.index].inChartTree=true;
      queue.push(edge.to);
    }
    if(queue.length!==g.vertices.length) throw Error("The hemisphere graph is disconnected.");
    for(const chart of charts) chart.offsetPi=value(chart.offset,R,S);
    const axes=g.face_representatives.map((d,i)=>{
      const chart=charts[vertex[d]],affine=add(scale(alpha[d],chart.sign),chart.offset),phasePi=value(affine,R,S);
      const chartPath=[];let v=vertex[d];
      while(charts[v].parent!==null){chartPath.unshift(charts[v].viaTransition);v=charts[v].parent;}
      return {i,d,vertex:vertex[d],weight:g.weights[i],affinePi:affine,phasePi,phase:Math.PI*phasePi,
        localAlphaPi:value(alpha[d],R,S),chartSign:chart.sign,chartOffsetPi:chart.offsetPi,chartPath};
    });
    const rows=[0,1,2].map(i=>{
      const row=add(axes[i].affinePi,scale(axes[3].affinePi,-1));
      row[0]+=g.weights[i]+g.weights[3];
      return row.map(integer);
    });
    const one=periodWord(rows,[-2,0]),tau=periodWord(rows,[0,2]);
    for(const word of [one,tau]){word.phasePi=value(word.affinePi,R,S);word.phase=Math.PI*word.phasePi;}
    const phase1=one.phase,phaseTau=tau.phase;
    const expected={phase1:-2*Math.PI*(R+config.shift[0]),phaseTau:2*Math.PI*(S+config.shift[1])};
    const distance=(a,b)=>Math.abs(Math.atan2(Math.sin(a-b),Math.cos(a-b)));
    const error1=distance(phase1,expected.phase1),errorTau=distance(phaseTau,expected.phaseTau);
    const complex=a=>({re:Math.cos(a),im:Math.sin(a)});
    const finished=globalThis.performance?.now?.()??Date.now();
    return {phase1,phaseTau,expected,error:Math.max(error1,errorTau),error1,errorTau,
      selectedwords:{one,tau},transitions,axes,charts,perimeters,
      lengthsPi:Object.fromEntries(Object.entries(lengths).map(([d,a])=>[d,value(a,R,S)])),
      multipliers:{one:complex(phase1),tau:complex(phaseTau)},
      projective:{one:complex(2*phase1),tau:complex(2*phaseTau)},
      metricBasisPi:basis.map(a=>value(a,R,S)),elapsedMs:finished-started,
      chartTreeMode:tree?"Stable base-core tree, continued along the inserted chains":"Spanning tree with affine period reconstruction",
      method:"Hemisphere gluing cocycle and local TV exponent lifts; affine integer period reconstruction"};
  }
  return {compute};
})();

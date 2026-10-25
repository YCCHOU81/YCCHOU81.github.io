/* Display contours through resolved samples; this module performs no solves.
 * A compact source edge is replaced by two open cusp sides only when both
 * sides have an actual, resolved sample on that compact edge (epsilon=inset).
 * Missing samples split paths. Neither a missing sample nor the cusp is joined
 * by a fabricated segment. Pixel-space PCHIP interpolation is shape preserving.
 */
(function (root) {
  "use strict";
  const finite = z => !!z && Number.isFinite(z.re) && Number.isFinite(z.im);
  const equal = (a,b) => finite(a) && finite(b) && a.re===b.re && a.im===b.im;
  function barycentric(point) {
    if(Array.isArray(point?.barycentric) && point.barycentric.length===3) {
      const b=point.barycentric.map(Number);
      if(b.every(Number.isFinite))return b;
    }
    return Number.isFinite(point?.r)&&Number.isFinite(point?.s)
      ?[1-2*point.s,2*(point.r+point.s)-1,1-2*point.r]:null;
  }
  function edgePosition(b,edge) {
    const others=[0,1,2].filter(i=>i!==edge),denominator=b[others[0]]+b[others[1]];
    return denominator>0?b[others[0]]/denominator:NaN;
  }
  function append(target,value) {
    const z=finite(value)?value:null;
    if(z&&equal(target.at(-1),z))return;
    if(z===null&&target.at(-1)===null)return;
    target.push(z);
  }
  function splitLinear(values) {
    const paths=[];let current=[];
    for(const value of values) {
      if(!finite(value)){if(current.length)paths.push(current);current=[];}
      else append(current,value);
    }
    if(current.length)paths.push(current);
    return paths;
  }
  function splitCycle(values) {
    const firstGap=values.findIndex(z=>!finite(z));
    if(firstGap<0) {
      const path=[];for(const value of values)append(path,value);
      if(path.length>1&&!equal(path[0],path.at(-1)))path.push(path[0]);
      return {paths:path.length?[path]:[],pathClosed:path.length?[true]:[]};
    }
    const rotated=[...values.slice(firstGap+1),...values.slice(0,firstGap),null];
    const paths=splitLinear(rotated);
    return {paths,pathClosed:paths.map(()=>false)};
  }
  function contours(record,view,options={}) {
    const inset=Number(record?.inset??.06),tolerance=Number(options.joinTolerance??1e-8);
    const showInterior=options.showInterior!==false,showTails=options.showTails!==false;
    const edges=[[],[],[]],joinCandidates=[[],[],[]],tails=[],joins=[];
    const points=record?.points??[],values=view?.points??[];
    for(let index=0;index<points.length;index++) {
      const b=barycentric(points[index]);if(!b)continue;
      for(let edge=0;edge<3;edge++)if(Math.abs(b[edge]-inset)<=tolerance) {
        edges[edge].push({t:edgePosition(b,edge),value:finite(values[index])?values[index]:null,sourceIndex:index});
      }
    }
    for(const edge of edges)edge.sort((a,b)=>a.t-b.t);
    for(let tailIndex=0;tailIndex<(view?.tails??[]).length;tailIndex++) {
      const tail=view.tails[tailIndex];if(!tail)continue;
      const edge=Number(tail.edge??tail.path?.index),t=Number(tail.t??tail.position??tail.path?.position);
      if(!Number.isInteger(edge)||edge<0||edge>2||!Number.isFinite(t))continue;
      const sourceNodes=tail.nodes??[],ordered=sourceNodes.map((node,index)=>({node,index,
        epsilon:Number(node.epsilon),value:finite(node.value)?node.value:finite(tail.values?.[index])?tail.values[index]:null}))
        .filter(item=>Number.isFinite(item.epsilon)).sort((a,b)=>b.epsilon-a.epsilon);
      let join=ordered.findIndex(item=>Math.abs(item.epsilon-inset)<=tolerance && finite(item.value) &&
        (item.node.status==null||item.node.status==="ok"));
      if(join>=0) {
        const b=barycentric(ordered[join].node.source??ordered[join].node);
        if(b&&(Math.abs(b[edge]-inset)>tolerance||Math.abs(edgePosition(b,edge)-t)>tolerance))join=-1;
      }
      const row={edge,t,tailIndex,joinIndex:join,
        values:(join<0?ordered:ordered.slice(join)).map(item=>item.value)};
      tails.push(row);
      if(join>=0&&edges[edge].length&&t>=edges[edge][0].t-tolerance&&t<=edges[edge].at(-1).t+tolerance) {
        const value=ordered[join].value,entry={edge,t,value,epsilon:ordered[join].epsilon,tailIndex};
        joins.push(entry);joinCandidates[edge].push({...entry,tail:row});
      }
    }
    const withJoins=(edge)=>{
      const nodes=edges[edge].filter(node=>!joinCandidates[edge].some(join=>Math.abs(join.t-node.t)<=tolerance));
      nodes.push(...joinCandidates[edge]);return nodes.sort((a,b)=>a.t-b.t);
    };
    // The source cycle is b_0=inset increasing t, then b_2 increasing,
    // then b_1 decreasing. Adjacent edges share the same computed vertex.
    const order=[{edge:0,reverse:false},{edge:2,reverse:false},{edge:1,reverse:true}];
    const compact=[];
    for(const {edge,reverse}of order) {
      const nodes=withJoins(edge);if(nodes.length<2)append(compact,null);if(reverse)nodes.reverse();
      for(const node of nodes)append(compact,node.value);
    }
    const boundaryComplete=edges.every(edge=>edge.length>=2)&&compact.length>=3&&compact.every(finite);
    const compactBoundary=boundaryComplete?splitCycle(compact).paths[0]:null;
    const boundaryIndex=value=>compactBoundary?compactBoundary.findIndex(z=>z===value||equal(z,value)):-1;
    const baseArcs=[];
    for(let edge=0;edge<3;edge++) {
      const candidates=joinCandidates[edge].slice().sort((a,b)=>a.t-b.t);
      if(candidates.length<2||candidates.at(-1).t-candidates[0].t<=tolerance)continue;
      const lower=candidates[0],upper=candidates.at(-1),nodes=withJoins(edge)
        .filter(node=>node.t>=lower.t-tolerance&&node.t<=upper.t+tolerance);
      baseArcs.push({edge,tailIndices:[lower.tailIndex,upper.tailIndex],points:nodes.map(node=>node.value),
        startIndex:boundaryIndex(lower.value),endIndex:boundaryIndex(upper.value),direction:edge===1?-1:1,
        complete:nodes.every(node=>finite(node.value))});
    }
    if(!showInterior) {
      const paths=showTails?tails.flatMap(tail=>splitLinear(tail.values)):[];
      return {paths,pathClosed:paths.map(()=>false),closed:false,compactBoundary,baseArcs,joins,splicedEdges:[],fallbackEdges:[]};
    }
    const cycle=[],splicedEdges=[],fallbackEdges=[],used=new Set();
    for(const {edge,reverse}of order) {
      const candidates=joinCandidates[edge].slice().sort((a,b)=>a.t-b.t),nodes=edges[edge],part=[];
      if(showTails&&candidates.length>=2&&candidates.at(-1).t-candidates[0].t>tolerance) {
        const lower=candidates[0],upper=candidates.at(-1);
        for(const node of nodes)if(node.t<lower.t-tolerance)append(part,node.value);
        for(const value of lower.tail.values)append(part,value);
        append(part,null);
        for(const value of upper.tail.values.slice().reverse())append(part,value);
        for(const node of nodes)if(node.t>upper.t+tolerance)append(part,node.value);
        used.add(lower.tailIndex);used.add(upper.tailIndex);splicedEdges.push(edge);
      } else {
        if(showTails&&tails.some(tail=>tail.edge===edge))fallbackEdges.push(edge);
        if(nodes.length<2)append(part,null);
        for(const node of withJoins(edge))append(part,node.value);
      }
      if(reverse)part.reverse();for(const value of part)append(cycle,value);
    }
    const result=splitCycle(cycle);
    if(showTails)for(const tail of tails)if(!used.has(tail.tailIndex)) {
      const paths=splitLinear(tail.values);result.paths.push(...paths);result.pathClosed.push(...paths.map(()=>false));
    }
    return {...result,closed:result.pathClosed.length===1&&result.pathClosed[0],compactBoundary,baseArcs,joins,splicedEdges,fallbackEdges};
  }
  // Nonuniform PCHIP derivatives in x and y separately. The common chord
  // parameter gives the same tangent direction on both sides of every sample.
  // Componentwise limiting keeps each cubic in its adjacent-sample rectangle.
  function derivative(left,right,hLeft,hRight) {
    if(left===0||right===0||Math.sign(left)!==Math.sign(right))return 0;
    const w1=2*hRight+hLeft,w2=hRight+2*hLeft;
    return(w1+w2)/(w1/left+w2/right);
  }
  function endpoint(first,second,hFirst,hSecond) {
    let slope=((2*hFirst+hSecond)*first-hFirst*second)/(hFirst+hSecond);
    if(Math.sign(slope)!==Math.sign(first))slope=0;
    else if(Math.sign(first)!==Math.sign(second)&&Math.abs(slope)>3*Math.abs(first))slope=3*first;
    return slope;
  }
  function smoothSegments(pixelPoints,options={}) {
    const points=[];
    for(const point of pixelPoints??[]) {
      const p=Array.isArray(point)?point:[point?.re,point?.im];
      // A caller must split unresolved gaps before interpolation.
      if(!Number.isFinite(p[0])||!Number.isFinite(p[1]))return{start:null,segments:[],closed:false};
      if(!points.length||p[0]!==points.at(-1)[0]||p[1]!==points.at(-1)[1])points.push([p[0],p[1]]);
    }
    let closed=Boolean(options.closed);
    if(points.length>1&&points[0][0]===points.at(-1)[0]&&points[0][1]===points.at(-1)[1]){points.pop();closed=true;}
    if(points.length<2)return{start:points[0]??null,segments:[],closed:false};
    if(points.length===2)closed=false;
    const length=points.length,count=closed?length:length-1,h=[],slopes=[];
    for(let i=0;i<count;i++) {
      const next=points[(i+1)%length],distance=Math.hypot(next[0]-points[i][0],next[1]-points[i][1]);
      h.push(Math.max(1e-12,Math.sqrt(distance)));slopes.push([(next[0]-points[i][0])/h[i],(next[1]-points[i][1])/h[i]]);
    }
    const tangents=points.map(()=>[0,0]);
    for(let i=0;i<length;i++)for(let dimension=0;dimension<2;dimension++) {
      if(closed){const previous=(i-1+count)%count;tangents[i][dimension]=derivative(slopes[previous][dimension],slopes[i][dimension],h[previous],h[i]);}
      else if(length===2)tangents[i][dimension]=slopes[0][dimension];
      else if(i===0)tangents[i][dimension]=endpoint(slopes[0][dimension],slopes[1][dimension],h[0],h[1]);
      else if(i===length-1)tangents[i][dimension]=endpoint(slopes[count-1][dimension],slopes[count-2][dimension],h[count-1],h[count-2]);
      else tangents[i][dimension]=derivative(slopes[i-1][dimension],slopes[i][dimension],h[i-1],h[i]);
    }
    const segments=[];
    for(let i=0;i<count;i++) {
      const next=(i+1)%length,p=points[i],q=points[next],cp1=[],cp2=[];
      for(let d=0;d<2;d++) {
        const low=Math.min(p[d],q[d]),high=Math.max(p[d],q[d]);
        cp1.push(Math.max(low,Math.min(high,p[d]+h[i]*tangents[i][d]/3)));
        cp2.push(Math.max(low,Math.min(high,q[d]-h[i]*tangents[next][d]/3)));
      }
      segments.push({cp1,cp2,end:[...q]});
    }
    return{start:[...points[0]],segments,closed};
  }
  root.ContourMath=Object.freeze({contours,smoothSegments,barycentric,edgePosition});
})(typeof globalThis!=="undefined"?globalThis:this);

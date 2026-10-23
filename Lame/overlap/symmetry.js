/* Exact relabeling of the three Lamé triangle edges, independent of numerics.
 * λ=(1−2s,1−2r,2r+2s−1), attached counts=(a,b,c).
 * C: (λ,a,b,c) -> cyclic order (3,1,2); τ -> 1/(1−τ).
 * N: interchange edges 2 and 3 and reverse orientation; τ -> 1−conj(τ).
 * Words are applied from left to right. No roots are replaced or reselected.
 */
(function (root) {
  "use strict";
  const words = Object.freeze(["", "C", "CC", "N", "NC", "NCC"]);
  function operations(word) {
    if (typeof word !== "string" || !/^[CN]*$/.test(word)) throw new TypeError("Use a word in C and N.");
    return word;
  }
  function applyCounts(counts, word = "") {
    let value = Array.from(counts);
    if (value.length !== 3 || !value.every(n => Number.isInteger(n) && n >= 0)) {
      throw new TypeError("Three nonnegative attachment counts are required.");
    }
    for (const operation of operations(word)) {
      const [a, b, c] = value;
      value = operation === "C" ? [c, a, b] : [a, c, b];
    }
    return value;
  }
  function applySource(point, word = "") {
    if (!point || !Number.isFinite(point.r) || !Number.isFinite(point.s)) return null;
    let {r, s} = point;
    for (const operation of operations(word)) {
      [r, s] = operation === "C" ? [s, 1 - r - s] : [1 - r - s, s];
    }
    return {r, s};
  }
  function applyTau(tau, word = "") {
    if (!tau) return null;
    let re = Array.isArray(tau) ? tau[0] : tau.re;
    let im = Array.isArray(tau) ? tau[1] : tau.im;
    if (!Number.isFinite(re) || !Number.isFinite(im) || im <= 0) return null;
    for (const operation of operations(word)) {
      if (operation === "N") {
        re = 1 - re;
      } else {
        const x = 1 - re, y = -im;
        const magnitude = Math.max(Math.abs(x), Math.abs(y));
        const xx = x / magnitude, yy = y / magnitude;
        const denominator = xx * xx + yy * yy;
        re = (xx / denominator) / magnitude;
        im = (-yy / denominator) / magnitude;
      }
    }
    return Number.isFinite(re) && Number.isFinite(im) && im > 0 ? {re, im} : null;
  }
  function cuspName(cusp) {
    if (cusp === 0 || cusp === "0" || cusp === "zero") return "zero";
    if (cusp === 1 || cusp === "1" || cusp === "one") return "one";
    if (cusp === Infinity || cusp === "∞" || cusp === "infinity" || cusp === "infty") return "infinity";
    throw new RangeError("A cusp is 0, 1, or infinity.");
  }
  function applyCusp(cusp, word = "") {
    let value = cuspName(cusp);
    for (const operation of operations(word)) {
      value = operation === "C" ? ({zero:"one",one:"infinity",infinity:"zero"})[value]
        : ({zero:"one",one:"zero",infinity:"infinity"})[value];
    }
    return value;
  }
  // Existing atlas edge coordinates are b=(1−2s,2(r+s)−1,1−2r).
  // At edge k, b_k=0 and the first other coordinate is the position t.
  function edgeBoundary(edge, position) {
    if (!Number.isInteger(edge) || edge<0 || edge>2 || !Number.isFinite(position) || position<0 || position>1) {
      throw new RangeError("Use an edge from 0 to 2 and a position from 0 to 1.");
    }
    const other=[0,1,2].filter(index=>index!==edge), b=[0,0,0];
    b[other[0]]=position;b[other[1]]=1-position;
    return {r:(1-b[2])/2,s:(1-b[0])/2};
  }
  function applyEdge(edge, position, word = "") {
    const source=applySource(edgeBoundary(edge,position),word);
    const b=[1-2*source.s,2*(source.r+source.s)-1,1-2*source.r];
    // Tracking a vertex would be ambiguous; tail positions lie strictly inside
    // an edge. Keep the edge label itself under exact barycentric permutations.
    let index=edge;
    for(const operation of operations(word))index=operation==="C" ? [2,0,1][index] : [0,2,1][index];
    const other=[0,1,2].filter(k=>k!==index), denominator=b[other[0]]+b[other[1]];
    return {edge:index,t:Math.max(0,Math.min(1,b[other[0]]/denominator)),source};
  }
  function applyTail(tail, word = "") {
    const oldEdge=Number(tail.edge ?? tail.path?.index), oldPosition=Number(tail.t ?? tail.position ?? tail.path?.position);
    const mapped=applyEdge(oldEdge,oldPosition,word);
    const nodes=(tail.nodes || []).map(node=>{
      const source=applySource(node.source || node,word);
      const result={...node,source,...(source || {}),tau:applyTau(node.tau,word),symmetryDerived:true};
      // Original decimal strings must not masquerade as transformed values.
      // Numerical chart changes here are binary64 display transformations.
      if(node.tauExact!=null){result.sourceTauExact=node.tauExact;result.tauExact=null;}
      if(node.sourceExact!=null){result.originalSourceExact=node.sourceExact;result.sourceExact=null;}
      for(const field of ["rExact","sExact","barycentricExact"]){
        if(node[field]!=null){result[`original${field[0].toUpperCase()}${field.slice(1)}`]=node[field];result[field]=null;}
      }
      if(node.barycentric!=null){
        result.originalBarycentric=node.barycentric;
        result.barycentric=source ? [1-2*source.s,2*(source.r+source.s)-1,1-2*source.r] : null;
      }
      if("value" in result)delete result.value;
      return result;
    });
    const cusp=applyCusp(tail.cusp ?? ["infinity","one","zero"][oldEdge],word);
    return {...tail,edge:mapped.edge,t:mapped.t,position:mapped.t,cusp,nodes,
      path:{...tail.path,kind:"edge",index:mapped.edge,position:String(mapped.t)},
      ...(tail.start ? {start:applySource(tail.start,word)} : {}),
      symmetry:{word,sourceEdge:oldEdge,sourcePosition:oldPosition,exactRelabeling:true},
      inheritedRepresentativeError:true};
  }
  function canonical(counts) {
    return applyCounts(counts).sort((a, b) => a - b);
  }
  function orbit(counts) {
    const found = new Set();
    const result = [];
    for (const word of words) {
      const transformed = applyCounts(counts, word);
      const key = transformed.join(",");
      if (found.has(key)) continue;
      found.add(key);
      const inverseWord = Array.from(word).reverse().map(operation => operation === "C" ? "CC" : "N").join("");
      result.push({counts: transformed, key, word, label: word || "Identity",
        orientation: word.includes("N") ? -1 : 1,
        mapSource: point => applySource(point, word),
        inverseSource: point => applySource(point, inverseWord),
        mapTau: tau => applyTau(tau, word)});
    }
    return result;
  }
  // The uniform inset triangular lattice is invariant under every edge
  // permutation. Return targetIndex[sourceIndex], throwing rather than guessing
  // if a caller supplied a mesh that is not invariant.
  function meshPermutation(points, word = "", tolerance = 5e-11) {
    operations(word);
    if (!Number.isFinite(tolerance) || tolerance <= 0) throw new RangeError("Use a positive source-mesh tolerance.");
    const quantize = p => `${Math.round(p.r / tolerance)},${Math.round(p.s / tolerance)}`;
    const lookup = new Map(points.map((p, i) => [quantize(p), i]));
    const seen = new Set();
    return points.map(point => {
      const transformed = applySource(point, word);
      let index = lookup.get(quantize(transformed));
      if (index === undefined || Math.hypot(points[index].r - transformed.r, points[index].s - transformed.s) > tolerance) {
        index = points.findIndex(p => Math.hypot(p.r - transformed.r, p.s - transformed.s) <= tolerance);
      }
      if (index < 0 || index === undefined || seen.has(index)) {
        throw new RangeError("The source mesh is not invariant under this relabeling.");
      }
      seen.add(index);
      return index;
    });
  }
  function derive(record, targetCounts = record.counts) {
    if (Array.isArray(targetCounts) && targetCounts.length === 3 && targetCounts.every((n, i) => n === record.counts[i])) return record;
    const word = typeof targetCounts === "string" ? operations(targetCounts) :
      orbit(record.counts).find(item => item.key === Array.from(targetCounts).join(","))?.word;
    if (word === undefined) throw new RangeError("The requested counts are not a permutation of the computed branch.");
    const permutation = meshPermutation(record.points, word);
    const tau = new Array(permutation.length).fill(null);
    const pointDiagnostics = new Array(permutation.length).fill(null);
    for (let i = 0; i < permutation.length; i += 1) {
      const index = permutation[i];
      tau[index] = applyTau(record.tau[i], word);
      if (record.pointDiagnostics?.[i]) pointDiagnostics[index] = {...record.pointDiagnostics[i], symmetryDerived: true};
    }
    const diagnostics = {...record.diagnostics};
    if (Array.isArray(diagnostics.unresolved)) {
      diagnostics.unresolved = diagnostics.unresolved.map(failure => {
        const result = {...failure};
        if (Number.isInteger(failure.index) && failure.index >= 0 && failure.index < permutation.length) {
          result.index = permutation[failure.index];
          // Use the actual common mesh node, avoiding a tiny affine-rounding
          // discrepancy between exported failure coordinates and point arrays.
          result.r = record.points[result.index].r;
          result.s = record.points[result.index].s;
        } else if (Number.isFinite(failure.r) && Number.isFinite(failure.s)) {
          Object.assign(result, applySource(failure, word));
        }
        return result;
      });
    }
    return {...record, counts: applyCounts(record.counts, word), tau,
      ...(Array.isArray(record.tails) ? {tails:record.tails.map(tail=>applyTail(tail,word))} : {}),
      reference: applyTau(record.reference, word),
      referenceSource: applySource(record.referenceSource || {r: 1/3, s: 1/3}, word),
      pointDiagnostics, elapsedMs: 0, sourceElapsedMs: record.elapsedMs,
      symmetry: {word, from: Array.from(record.counts), sourceCounts: Array.from(record.counts), exactRelabeling: true},
      diagnostics: {...diagnostics,
        inheritedRepresentativeError: true,
        symmetry: "Exact permutation of triangle edges and transported marking. Numerical errors are inherited from the representative chart, not independently recomputed; unresolved samples remain unresolved."}};
  }
  root.LameSymmetry = Object.freeze({words, applyCounts, applySource, applyTau, applyCusp,edgeBoundary,applyEdge,applyTail,
    canonical, orbit, meshPermutation, derive});
  root.ComparisonSymmetry = root.LameSymmetry;
})(typeof globalThis !== "undefined" ? globalThis : this);

/* Display coordinates only. Measured τ samples and their tails stay unchanged.
 * Affine maps use one real scale in both directions. A cusp chart is a genuine
 * Möbius coordinate, never a compression of the imaginary coordinate alone.
 */
(function (root) {
  "use strict";
  const commonReference = Object.freeze({re: .5, im: Math.sqrt(3)/2});
  const aliases = Object.freeze({tau: "raw", centered: "affine", equal: "equal-affine", commonDisk: "common-disk",
    cuspInfinity: "cuspinfty", qInfinity: "qinfty", cuspZero: "cusp0", cuspOne: "cusp1"});
  const modes = Object.freeze(["raw", "affine", "equal-affine", "common-disk", "disk", "cusp0", "cusp1", "cuspinfty", "q0", "q1", "qinfty"]);
  const descriptions = Object.freeze({
    raw: {label: "Original τ", xLabel: "Re τ", yLabel: "Im τ", preservesOverlap: true,
      description: "Original upper-half-plane coordinates, with a linear imaginary axis. The positions and overlaps are unchanged."},
    affine: {label: "Centered affine", xLabel: "Re Z", yLabel: "Im Z", preservesOverlap: false,
      description: "Z = (τ − τ*) / Im τ*. This is a translation and a uniform scaling for each attachment. It compares recentered shapes; relative positions and overlaps between different attachments are changed."},
    "equal-affine": {label: "Equal interior RMS", xLabel: "Re Z", yLabel: "Im Z", preservesOverlap: false,
      description: "Z = (τ − τ*) / R, where R is the RMS distance of the interior samples from τ*. One uniform scaling is also applied to every tail. It compares recentered shapes, not the original overlaps."},
    "common-disk": {label: "Common disk · overlaps preserved", xLabel: "Re ζ", yLabel: "Im ζ", preservesOverlap: true,
      description: "One fixed disk coordinate for every attachment: actual overlaps are preserved, and all three cusp tails fit in one picture. Cusp labels mark limits; dashed curves stop at computed endpoints. Sample points can be shown separately. The formula is in Coordinate details below."},
    disk: {label: "Centered disk", xLabel: "Re ζ", yLabel: "Im ζ", preservesOverlap: false,
      description: "ζ = (τ − τ*) / (τ − conjugate(τ*)). Each attachment has its own center. Infinity goes to the unit circle; this optional view compares shapes rather than original positions or overlaps."},
    cusp0: {label: "Cusp 0 chart", xLabel: "Re w", yLabel: "Im w", preservesOverlap: true,
      description: "w = −1/τ. Cusp 0 moves to infinity. The same one-to-one coordinate is used for every attachment, so overlaps are preserved."},
    cusp1: {label: "Cusp 1 chart", xLabel: "Re w", yLabel: "Im w", preservesOverlap: true,
      description: "w = 1/(1 − τ). Cusp 1 moves to infinity. The same one-to-one coordinate is used for every attachment, so overlaps are preserved."},
    cuspinfty: {label: "Cusp ∞ chart", xLabel: "Re w = Re τ", yLabel: "Im w = Im τ", preservesOverlap: true,
      description: "w = τ. Cusp infinity is already at infinity. Both coordinate axes remain linear and the original overlaps are unchanged."},
    q0: {label: "Cusp 0 exponential", xLabel: "Re q₀", yLabel: "Im q₀", preservesOverlap: "strip",
      description: "q₀ = exp(πi(−1/τ)). The cusp becomes the origin. This coordinate is periodic under w → w + 2, so overlap is faithful only within a width-two strip. Very deep values can underflow at floating-point precision."},
    q1: {label: "Cusp 1 exponential", xLabel: "Re q₁", yLabel: "Im q₁", preservesOverlap: "strip",
      description: "q₁ = exp(πi/(1 − τ)). The cusp becomes the origin. This coordinate is periodic under w → w + 2, so overlap is faithful only within a width-two strip. Very deep values can underflow at floating-point precision."},
    qinfty: {label: "Cusp ∞ exponential", xLabel: "Re q∞", yLabel: "Im q∞", preservesOverlap: "strip",
      description: "q∞ = exp(πiτ). The cusp becomes the origin. This coordinate is periodic under τ → τ + 2, so overlap is faithful only within a width-two strip. Very deep values can underflow at floating-point precision."}
  });
  function modeName(mode = "raw") {
    const name = aliases[mode] || mode;
    if (!modes.includes(name)) throw new RangeError("Unknown comparison coordinate.");
    return name;
  }
  function complex(value) {
    if (!value) return null;
    const re = Array.isArray(value) ? value[0] : value.re;
    const im = Array.isArray(value) ? value[1] : value.im;
    return Number.isFinite(re) && Number.isFinite(im) ? {re, im} : null;
  }
  function upper(value) { const z = complex(value); return z && z.im > 0 ? z : null; }
  function reciprocal(z) {
    const magnitude = Math.max(Math.abs(z.re), Math.abs(z.im));
    if (!(magnitude > 0) || !Number.isFinite(magnitude)) return null;
    const x = z.re / magnitude, y = z.im / magnitude, d = x*x + y*y;
    const result = {re: (x/d)/magnitude, im: (-y/d)/magnitude};
    return complex(result);
  }
  function cuspChart(tau, cusp = "infinity") {
    const z = upper(tau);
    if (!z) return null;
    if (cusp === "infinity" || cusp === "infty" || cusp === Infinity || cusp === "∞") return z;
    if (cusp === "zero" || cusp === "0" || cusp === 0) {
      const inverse = reciprocal(z);
      return inverse ? {re: -inverse.re, im: -inverse.im} : null;
    }
    if (cusp === "one" || cusp === "1" || cusp === 1) return reciprocal({re: 1-z.re, im: -z.im});
    throw new RangeError("A cusp is 0, 1, or infinity.");
  }
  function inverseCusp(w, cusp = "infinity") {
    const z = upper(w);
    if (!z) return null;
    if (cusp === "infinity" || cusp === "infty" || cusp === Infinity || cusp === "∞") return z;
    if (cusp === "zero" || cusp === "0" || cusp === 0) return cuspChart(z, 0);
    if (cusp === "one" || cusp === "1" || cusp === 1) {
      const inverse = reciprocal(z);
      return inverse ? {re: 1-inverse.re, im: -inverse.im} : null;
    }
    throw new RangeError("A cusp is 0, 1, or infinity.");
  }
  function disk(tau, reference) {
    const t = upper(tau), p = upper(reference);
    if (!t || !p) return null;
    if (root.ComparisonMath?.disk) return root.ComparisonMath.disk(t, p);
    const scale = Math.max(Math.abs(t.re), Math.abs(p.re), t.im, p.im);
    const x = t.re/scale - p.re/scale, y = t.im/scale-p.im/scale, v = t.im/scale+p.im/scale;
    const d = x*x+v*v;
    return complex({re: (x*x+y*v)/d, im: x*(y-v)/d});
  }
  function exponential(w) {
    if (!w) return null;
    const radius = Math.exp(-Math.PI*w.im);
    const angle = Math.PI*(w.re % 2);
    return {re: radius*Math.cos(angle), im: radius*Math.sin(angle)};
  }
  function transform(tau, reference, mode = "raw", options = {}) {
    const name = modeName(mode), t = upper(tau);
    if (!t) return null;
    if (name === "raw" || name === "cuspinfty") return t;
    if (name === "cusp0" || name === "cusp1") return cuspChart(t, name === "cusp0" ? 0 : 1);
    if (name[0] === "q") return exponential(cuspChart(t, name === "q0" ? 0 : name === "q1" ? 1 : "infinity"));
    if (name === "common-disk") return disk(t, commonReference);
    if (name === "disk") return disk(t, reference);
    const p = upper(reference);
    if (!p) return null;
    const divisor = name === "affine" ? p.im : options.absoluteScale ?? p.im*(options.scale ?? 1);
    if (!Number.isFinite(divisor) || divisor <= 0) return null;
    return complex({re: (t.re-p.re)/divisor, im: (t.im-p.im)/divisor});
  }
  function rms(points) {
    const resolved = points.filter(Boolean);
    if (!resolved.length) return null;
    const largest = resolved.reduce((m,z) => Math.max(m,Math.abs(z.re),Math.abs(z.im)),0);
    if (largest === 0) return 0;
    return largest*Math.sqrt(resolved.reduce((sum,z) => sum+(z.re/largest)**2+(z.im/largest)**2,0)/resolved.length);
  }
  function bounds(points) {
    let minRe=Infinity, maxRe=-Infinity, minIm=Infinity, maxIm=-Infinity, maxAbs=0, count=0;
    for (const z of points) {
      if (!z || !Number.isFinite(z.re) || !Number.isFinite(z.im)) continue;
      minRe=Math.min(minRe,z.re);maxRe=Math.max(maxRe,z.re);minIm=Math.min(minIm,z.im);maxIm=Math.max(maxIm,z.im);
      maxAbs=Math.max(maxAbs,Math.hypot(z.re,z.im));count+=1;
    }
    return count ? {minRe,maxRe,minIm,maxIm,maxAbs,width:maxRe-minRe,height:maxIm-minIm,count} : null;
  }
  function prepare(record, options = {}) {
    const mode = modeName(options.mode), raw = Array.isArray(record?.tau) ? record.tau : [];
    const originalReference = upper(record?.reference), affine = raw.map(t => transform(t,originalReference,"affine"));
    const affineRms = rms(affine), absoluteRms = originalReference && affineRms !== null ? originalReference.im*affineRms : null;
    const scale = mode === "equal-affine" ? affineRms : 1;
    const coordinateOptions = {scale,absoluteScale:absoluteRms};
    const map = value => transform(value,originalReference,mode,coordinateOptions);
    const points = raw.map(map);
    let tailResolved=0,tailTotal=0,underflow=0;
    const tails = (record?.tails || []).map(tail => {
      const originalNodes = tail.nodes || tail.points || [];
      const nodes = originalNodes.map(node => {
        const originalTau = node?.tau ?? node;
        const value = map(originalTau);
        tailTotal+=1;if(value) tailResolved+=1;
        const qUnderflow = !!value && mode[0] === "q" && value.re === 0 && value.im === 0;
        if(qUnderflow)underflow+=1;
        return {...node,tau:originalTau,value,qUnderflow};
      });
      return {...tail,nodes,values:nodes.map(node=>node.value),points:nodes.map(node=>node.value)};
    });
    const resolved=points.filter(Boolean).length,total=points.length;
    const tailValues=tails.flatMap(tail=>tail.values);
    const interiorComplete=total>0 && resolved===total && record?.complete!==false && !record?.cancelled && record?.status!=="partial";
    const tailsComplete=tails.length>0 && tailTotal>0 && tailResolved===tailTotal && tails.every(t=>t.status!=="partial" && t.status!=="unresolved" && t.status!=="stopped");
    const scaleAvailable=mode!=="equal-affine" || (absoluteRms!==null && absoluteRms>0 && Number.isFinite(absoluteRms));
    return {points,tails,reference:map(originalReference),originalReference,scale,rms:affineRms,absoluteRms,
      bounds:bounds([...points,...tailValues]),interiorBounds:bounds(points),tailBounds:bounds(tailValues),
      mode,label:descriptions[mode].label,description:descriptions[mode].description,
      complete:interiorComplete && (tails.length===0 || tailsComplete) && scaleAvailable,
      interiorComplete,tailsComplete,resolved,total,tailResolved,tailTotal,underflow,
      normalizationLabel:!scaleAvailable ? "Interior scale unavailable" : descriptions[mode].label,
      message:!scaleAvailable ? "The common interior mesh has not resolved a nonzero RMS scale; its tails are retained in the original data." :
        underflow ? "Very deep exponential coordinates underflow to zero; the original τ samples remain available in linear cusp charts." : ""};
  }
  function describe(mode = "raw") {return {...descriptions[modeName(mode)]};}
  root.ViewMath = Object.freeze({modes,commonReference,modeName,complex,upper,reciprocal,cuspChart,inverseCusp,transform,prepare,rms,bounds,describe});
  root.LameViewMath = root.ViewMath;
})(typeof globalThis !== "undefined" ? globalThis : this);

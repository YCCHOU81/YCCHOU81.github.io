"use strict";

// Shared marked charts. Trajectories contain computed samples only.
window.AtlasView = (() => {
  const sourceVertices = [[328, 312], [328, 32], [48, 32]];
  const sourceXY = (r, s) => [48 + 560 * r, 312 - 560 * s];
  const sourceRS = (x, y) => ({r: (x - 48) / 560, s: (312 - y) / 560});
  const sourceAtlasViewBox = "0 0 400 400";
  const sourceAtlasXY = (r, s) => [40 + 320 * r, 360 - 320 * s];
  const sourceAtlasRS = (x, y) => ({r: (x - 40) / 320, s: (360 - y) / 320});
  function monodromyMap(r, s, counts = [0, 0, 0]) {
    const p = (counts[0] + counts[2]) % 2, q = (counts[1] + counts[2]) % 2;
    // Ordinary characters are taken modulo integer shifts and simultaneous
    // inversion. For Delta_3 use the inverse representative (1-R,1-S), so
    // the four cells tile the triangle R>=0, S>=0, R+S<=1.
    const topLeft = p === 1 && q === 0;
    return {r: topLeft ? .5 - r : q ? 1 - r - p / 2 : r,
      s: topLeft ? 1 - s : q ? .5 - s : s,
      index: [[0, 2], [3, 1]][p][q], shiftR: (counts[0] + counts[2]) / 2,
      shiftS: (counts[1] + counts[2]) / 2, reflected: Boolean(p || q), p, q};
  }
  function monodromyInverse(r, s, counts = [0, 0, 0]) {
    const p = (counts[0] + counts[2]) % 2, q = (counts[1] + counts[2]) % 2;
    if (p === 1 && q === 0) return {r: .5 - r, s: 1 - s};
    return {r: q ? 1 - p / 2 - r : r, s: q ? .5 - s : s};
  }
  const monodromyVertices = (counts = [0, 0, 0]) => [[.5, 0], [.5, .5], [0, .5]].map(([r, s]) => {
    const point = monodromyMap(r, s, counts); return [point.r, point.s];
  });
  function cuspRational(counts = [0, 0, 0]) {
    let numerator = 2 * counts[2] + 1, denominator = 2 * (counts[1] + counts[2] + 1);
    const gcd = (a, b) => b ? gcd(b, a % b) : a;
    const divisor = gcd(numerator, denominator); numerator /= divisor; denominator /= divisor;
    return {numerator, denominator, value: numerator / denominator, label: `${numerator}/${denominator}`, status: "predicted"};
  }
  const finiteTau = z => z && Number.isFinite(z.re) && Number.isFinite(z.im) && z.im >= 0;
  const cuspForEdge = ["∞", "1", "0"];

  function inMarkedDomain(z, tolerance = 1e-12) {
    return finiteTau(z) && z.re >= -tolerance && z.re <= 1 + tolerance &&
      z.im * z.im - z.re * (1 - z.re) >= -tolerance;
  }

  // Connect nearby numerical samples in tau coordinates, then compress height.
  // An exact quadratic minimum checks the whole chord, not just its endpoints.
  function interpolateDomainSegment(a, b, subdivisions = 12) {
    if (!inMarkedDomain(a) || !inMarkedDomain(b)) return [];
    const dx = b.re - a.re, dy = b.im - a.im;
    const A = dx * dx + dy * dy;
    const B = 2 * a.im * dy - dx * (1 - 2 * a.re);
    const C = a.im * a.im - a.re * (1 - a.re);
    const t = A ? Math.max(0, Math.min(1, -B / (2 * A))) : 0;
    if (A * t * t + B * t + C < -1e-12) return [];
    return Array.from({length: subdivisions + 1}, (_, i) => ({re: a.re + dx * i / subdivisions, im: a.im + dy * i / subdivisions}));
  }

  function cuspInfo(path, counts = [0, 0, 0]) {
    if (path?.kind === "edge" && Number.isInteger(path.index) && path.index >= 0 && path.index <= 2) {
      const cusp = cuspForEdge[path.index];
      return {kind: "cusp", cusp, label: `Approaching cusp ${cusp}`,
        detail: `The open source edge tends to cusp ${cusp} in this marked Γ(2) chart, for the base and the selected attachment. The endpoint itself is not computed.`};
    }
    return {kind: "boundary", cusp: null, label: "Source vertex: direction-dependent boundary limit",
      detail: "A source vertex can have a finite boundary limit depending on the direction and attachment. This trace does not assume a cusp endpoint."};
  }

  function drawArrow(ctx, from, to, color, options = {}) {
    if (![...from, ...to].every(Number.isFinite)) return;
    const dx = to[0] - from[0], dy = to[1] - from[1], length = Math.hypot(dx, dy);
    if (length < 2) return;
    const ux = dx / length, uy = dy / length, head = Math.min(9, length * .35);
    ctx.save(); ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = options.width ?? 1.8;
    if (options.dashed) ctx.setLineDash([5, 5]);
    ctx.beginPath(); ctx.moveTo(...from); ctx.lineTo(...to); ctx.stroke();
    ctx.setLineDash([]); ctx.beginPath(); ctx.moveTo(...to);
    ctx.lineTo(to[0] - head * ux + head * .48 * uy, to[1] - head * uy - head * .48 * ux);
    ctx.lineTo(to[0] - head * ux - head * .48 * uy, to[1] - head * uy + head * .48 * ux);
    ctx.closePath(); ctx.fill(); ctx.restore();
  }

  function drawDomain(canvas, options = {}) {
    const width = canvas.getBoundingClientRect().width;
    if (!(width > 0)) return null;
    const height = options.height ?? 580, ratio = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
    canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d"); ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    const chartWidth = Math.max(60, Math.min(width - 82, 240));
    const box = {x: (width - chartWidth) / 2 + 10, y: 45, w: chartWidth, h: height - 99};
    let realMin=0,realMax=1;
    if (options.markedStrip) {
      const values=[];
      for(const x of options.realRange??[])if(Number.isFinite(x))values.push(x);
      const retain=z=>{if(finiteTau(z))values.push(z.re);};
      for(const p of options.points??[])retain(p.tau);
      for(const p of options.paths??[])for(const z of p.points??[])retain(z);
      for(const b of options.bands??[])for(const z of [...(b.left??[]),...(b.right??[])])retain(z);
      for(const layer of options.meshLayers??[])for(const sample of layer.mesh?.samples??[])retain(sample?.[layer.field??"result"]?.tau);
      retain(options.arrow?.from);retain(options.arrow?.to);
      for(const ref of options.references??[])if(Number.isFinite(ref.re))values.push(ref.re);
      if(values.length){const lo=Math.min(...values),hi=Math.max(...values);if(lo< -1e-10)realMin=Math.floor((lo-.06)*4)/4;if(hi>1+1e-10)realMax=Math.ceil((hi+.06)*4)/4;}
    }
    const realSpan=realMax-realMin;
    const map = z => [box.x + (z.re-realMin)/realSpan * box.w, box.y + box.h / (1 + z.im)];
    canvas.dataset.realMin=String(realMin);canvas.dataset.realMax=String(realMax);
    // The new seeds have different marked images. Keep their raw coordinates
    // visible rather than clipping them to the classical n=1 reference cell.
    const valid = options.markedStrip ? z => finiteTau(z) && z.re >= realMin-1e-12 && z.re <= realMax+1e-12 : inMarkedDomain;
    const segment = (a,b) => options.markedStrip
      ? (valid(a) && valid(b) ? Array.from({length:13},(_,i)=>({re:a.re+(b.re-a.re)*i/12,im:a.im+(b.im-a.im)*i/12})) : [])
      : interpolateDomainSegment(a,b);
    const bottom = box.y + box.h, px = x => box.x + (x-realMin)/realSpan * box.w;
    const line = (a, b) => {ctx.moveTo(...a); ctx.lineTo(...b);};
    function domainPath() {
      if (options.markedStrip) {ctx.beginPath();ctx.rect(box.x,box.y,box.w,box.h);return;}
      ctx.beginPath(); ctx.moveTo(px(0), box.y);
      for (let k = 0; k <= 160; k++) {
        const x = k / 160; ctx.lineTo(...map({re: x, im: Math.sqrt(x * (1 - x))}));
      }
      ctx.lineTo(px(1), box.y); ctx.closePath();
    }
    domainPath(); ctx.fillStyle = "#f0f2ed"; ctx.fill();
    ctx.save(); domainPath(); ctx.clip();
    ctx.strokeStyle = "#dce2da"; ctx.lineWidth = .7; ctx.beginPath();
    for (const y of [.5, 1, 2, 4, 9]) line(map({re: realMin, im: y}), map({re: realMax, im: y}));
    for (let x=Math.ceil(realMin*4)/4;x<realMax;x+=.25) line([px(x), box.y], [px(x), bottom]);
    ctx.stroke();

    for (const reference of options.references ?? []) {
      if (!Number.isFinite(reference.re) || reference.re < realMin || reference.re > realMax) continue;
      ctx.save(); ctx.strokeStyle = reference.color ?? "#7c7567"; ctx.lineWidth = 1.2; ctx.setLineDash([3, 5]);
      ctx.beginPath(); ctx.moveTo(px(reference.re), box.y); ctx.lineTo(px(reference.re), bottom); ctx.stroke(); ctx.restore();
    }

    for (const band of options.bands ?? []) {
      const left = band.left ?? [], right = band.right ?? [];
      ctx.fillStyle = (band.color ?? "#285f65") + "24";
      for (let i = 1; i < Math.min(left.length, right.length); i++) {
        const quad = [left[i - 1], right[i - 1], right[i], left[i]];
        for (const tri of [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]]) {
          if (!tri.every(valid) || !tri.every((p, k) => segment(p, tri[(k + 1) % 3]).length)) continue;
          ctx.beginPath(); ctx.moveTo(...map(tri[0])); ctx.lineTo(...map(tri[1])); ctx.lineTo(...map(tri[2])); ctx.closePath(); ctx.fill();
        }
      }
    }

    for (const layer of options.meshLayers ?? []) {
      const mesh = layer.mesh; if (!mesh) continue;
      const color = layer.color ?? "#285f65", field = layer.field ?? "result";
      const data = mesh.samples.map(sample => {
        const branch = sample?.[field];
        return branch?.status === "ok" && finiteTau(branch.tau) ? branch.tau : null;
      });
      const perimeter = new Map();
      ctx.fillStyle = color + "1a"; ctx.strokeStyle = color + "35"; ctx.lineWidth = .45;
      for (const face of mesh.triangles) {
        for (let j = 0; j < 3; j++) {
          const a = face[j], b = face[(j + 1) % 3], key = [Math.min(a, b), Math.max(a, b)].join(",");
          perimeter.set(key, (perimeter.get(key) ?? 0) + 1);
        }
        const vertices = face.map(i => data[i]); if (vertices.some(p => !p)) continue;
        ctx.beginPath(); ctx.moveTo(...map(vertices[0])); ctx.lineTo(...map(vertices[1]));
        ctx.lineTo(...map(vertices[2])); ctx.closePath(); ctx.fill(); ctx.stroke();
      }
      ctx.strokeStyle = color; ctx.lineWidth = 1.3; ctx.setLineDash([5, 4]); ctx.beginPath();
      for (const [key, count] of perimeter) {
        if (count !== 1) continue;
        const [a, b] = key.split(",").map(Number); if (data[a] && data[b]) line(map(data[a]), map(data[b]));
      }
      ctx.stroke(); ctx.setLineDash([]);

    }

    for (const path of options.paths ?? []) {
      ctx.strokeStyle = path.color ?? "#b55e3c"; ctx.lineWidth = path.width ?? 2; ctx.beginPath();
      (path.points ?? []).forEach((z, i, all) => {
        if (!valid(z)) return;
        const joined = i ? segment(all[i - 1], z) : [];
        if (joined.length) joined.slice(1).forEach(p => ctx.lineTo(...map(p)));
        else ctx.moveTo(...map(z));
      });
      ctx.stroke();
      (path.points ?? []).forEach((z, i) => {
        if (!valid(z)) return;
        ctx.fillStyle = path.color ?? "#b55e3c"; ctx.strokeStyle = "white"; ctx.lineWidth = 1.2;
        ctx.beginPath(); ctx.arc(...map(z), i === path.selected ? 5 : (path.radius ?? 2.5), 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
      });
    }

    if (options.arrow && segment(options.arrow.from, options.arrow.to).length) {
      const a = map(options.arrow.from), b = map(options.arrow.to), len = Math.hypot(b[0] - a[0], b[1] - a[1]);
      if (len > 12) {
        const direction = [(b[0] - a[0]) / len, (b[1] - a[1]) / len];
        drawArrow(ctx, [a[0] + 6 * direction[0], a[1] + 6 * direction[1]],
          [b[0] - 7 * direction[0], b[1] - 7 * direction[1]], "#253b33", {width: 2});
      }
    }
    for (const point of options.points ?? []) {
      if (!valid(point.tau)) continue;
      const xy = map(point.tau), color = point.color ?? "#285f65";
      ctx.fillStyle = color + "28"; ctx.beginPath(); ctx.arc(...xy, 11, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = color; ctx.strokeStyle = "white"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(...xy, 4.7, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    }
    ctx.restore();
    if (options.markedStrip) {
      ctx.strokeStyle="#a6b2a4";ctx.lineWidth=.8;
      for (const [lo,hi] of [[0,1],[0,.5],[.5,1]]) {
        ctx.beginPath();for(let k=0;k<=120;k++){
          const x=lo+(hi-lo)*k/120,p=map({re:x,im:Math.sqrt(Math.max(0,(x-lo)*(hi-x)))});
          if(k)ctx.lineTo(...p);else ctx.moveTo(...p);
        }ctx.stroke();
      }
    }
    ctx.beginPath(); ctx.moveTo(px(0), box.y);
    for (let k = 0; k <= 160; k++) {
      const x = k / 160; ctx.lineTo(...map({re: x, im: Math.sqrt(x * (1 - x))}));
    }
    ctx.lineTo(px(1), box.y); ctx.strokeStyle = "#a6b2a4"; ctx.lineWidth = 1; ctx.stroke();
    // The dashed top represents the ideal infinity end of the compactified chart.
    ctx.setLineDash([4, 4]); ctx.beginPath(); line([px(realMin), box.y], [px(realMax), box.y]); ctx.stroke(); ctx.setLineDash([]);
    ctx.font = "11px system-ui"; ctx.fillStyle = "#596e62"; ctx.textAlign = "center";
    ctx.fillText("cusp ∞ · open upper end", box.x+box.w/2, box.y - 17);
    let referenceLabelRow = 0;
    for (const reference of options.references ?? []) {
      if (!Number.isFinite(reference.re) || reference.re < realMin || reference.re > realMax) continue;
      const x = px(reference.re), text = reference.label ?? String(reference.re);
      ctx.save(); ctx.font = "11px ui-monospace, monospace";
      const fraction=(reference.re-realMin)/realSpan;
      ctx.textAlign = fraction > .72 ? "right" : fraction < .28 ? "left" : "center";
      const padding = 5, measured = ctx.measureText(text).width;
      const left = ctx.textAlign === "right" ? x - measured - padding : ctx.textAlign === "left" ? x - padding : x - measured / 2 - padding;
      const labelY = box.y + 19 + 19 * referenceLabelRow++;
      ctx.fillStyle = "#fcfcf7ed"; ctx.fillRect(left, labelY - 13, measured + padding * 2, 18);
      ctx.fillStyle = reference.color ?? "#675e4f"; ctx.fillText(text, x, labelY); ctx.restore();
    }
    for (const cusp of (options.markedStrip ? [0,.5,1] : [0, 1])) {
      ctx.fillStyle = "#fcfcf7"; ctx.strokeStyle = "#748777"; ctx.lineWidth = 1.3;
      ctx.beginPath(); ctx.arc(px(cusp), bottom, 4, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
      ctx.fillStyle = "#596e62"; ctx.fillText(`cusp ${cusp===.5?"½":cusp}`, px(cusp), bottom + 20);
    }
    ctx.fillStyle = "#596e62"; ctx.textAlign = "right";
    for (const y of [.5, 1, 2, 4, 9]) ctx.fillText(String(y), box.x - 8, map({re: 0, im: y})[1] + 4);
    if(options.markedStrip && (realMin<0 || realMax>1)){
      ctx.textAlign="left";ctx.fillText(String(realMin),box.x,bottom+36);
      ctx.textAlign="right";ctx.fillText(String(realMax),box.x+box.w,bottom+36);
    }
    ctx.textAlign = "center"; ctx.fillText("Re τ", box.x+box.w/2, height - 7);
    ctx.save(); ctx.translate(11, box.y + box.h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("Im τ · compressed height", 0, 0); ctx.restore();

    const corners = options.zoomWindow?.corners?.filter(finiteTau);
    if (corners?.length > 1) {
      const mapped = corners.map(map), xs = mapped.map(p => p[0]), ys = mapped.map(p => p[1]);
      const raw = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
      const bound = (n, lo, hi) => Math.max(lo, Math.min(hi, n));
      const x0 = bound(raw[0], box.x, box.x + box.w), x1 = bound(raw[1], box.x, box.x + box.w);
      const y0 = bound(raw[2], box.y, bottom), y1 = bound(raw[3], box.y, bottom);
      const tiny = x1 - x0 < 12 || y1 - y0 < 12;
      const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
      const rw = Math.max(12, x1 - x0), rh = Math.max(12, y1 - y0);
      ctx.save(); ctx.strokeStyle = "#253b33"; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
      ctx.strokeRect(cx - rw / 2, cy - rh / 2, rw, rh); ctx.setLineDash([]);
      const label = tiny ? "zoom location*" : "zoom window";
      const labelY = Math.max(box.y + 15, Math.min(bottom - 12, cy - rh / 2 - 7));
      ctx.font = "11px system-ui"; ctx.textAlign = "center"; ctx.fillStyle = "#fcfcf7e8";
      ctx.fillRect(cx - 53, labelY - 12, 106, 17); ctx.fillStyle = "#253b33"; ctx.fillText(label, cx, labelY);
      ctx.restore();
    }
    return {map, box, ctx,realMin,realMax};
  }

  return {sourceXY, sourceRS, sourceVertices, sourceAtlasViewBox, sourceAtlasXY, sourceAtlasRS, monodromyMap, monodromyInverse, monodromyVertices, cuspRational,
    finiteTau, inMarkedDomain, interpolateDomainSegment, cuspForEdge, cuspInfo, drawArrow, drawDomain};
})();

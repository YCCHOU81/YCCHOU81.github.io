"use strict";

/* Ordinary-character coordinates retain the signs lost by the projective
   character.  All translations below are exact half-integer parity data. */
globalThis.GraphShift = (() => {
  const base = "#315d83", attached = "#b6403b", grey = "#6f777a", faint = "#e5e9e5";
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const mod2 = n => ((Math.round(n) % 2) + 2) % 2;
  const half = n => n === .5 ? "½" : n === 1.5 ? "1½" : String(n);
  const pair = p => `(${p.map(half).join(", ")})`;
  const label = e => e.id.replace("l1+", "ℓ₁").replace("l2+", "ℓ₂").replace("l3+", "ℓ₃").replace("l+", "ℓ").replace("m+", "m");
  const fmt = n => Number(n.toFixed(3));
  const math = body => `<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">${body}</math>`;
  const number = value => `<mn>${escape(value)}</mn>`;
  const coordinate = value => value === .5 ? '<mfrac><mn>1</mn><mn>2</mn></mfrac>' : number(value);
  const mathPair = values => `<mrow><mo stretchy="false">(</mo>${values.map(coordinate).join('<mo>,</mo>')}<mo stretchy="false">)</mo></mrow>`;
  const modulo = '<mspace width="0.5em"/><mo stretchy="false">(</mo><mi mathvariant="normal">mod</mi><mspace width="0.25em"/><msup><mi>ℤ</mi><mn>2</mn></msup><mo stretchy="false">)</mo>';
  const character = (period, decorated = false) => `<mrow><msub><mi>χ</mi>${decorated ? '<mtext>cur</mtext>' : '<mn>0</mn>'}</msub><mo stretchy="false">(</mo>${period === "1" ? '<mn>1</mn>' : '<mi>τ</mi>'}<mo stretchy="false">)</mo></mrow>`;
  let instance = 0;

  function create({ onSource } = {}) {
    const get = id => document.getElementById(id);
    const svg = get("shift-triangle");
    if (!svg) throw new Error("The monodromy-shift drawing is missing.");
    const marker = `shift-arrow-${++instance}`;
    let current = null, coordinates = null, dragging = null;
    svg.style.touchAction = "none";
    svg.setAttribute("viewBox", "0 0 640 390");
    svg.setAttribute("tabindex", "0");
    svg.setAttribute("role", "group");
    svg.setAttribute("aria-label", "Undecorated seed and current total decoration: ordinary monodromy triangles. Drag either point, or use the arrow keys to change the source point.");

    function triangle(seed) {
      const x = seed === "bouquet" ? .5 : 1, y = seed === "tetra" ? 1 : .5;
      return { x, y, vertices: [[x, 0], [0, y], [x, y]] };
    }
    function fallbackProject(p) {
      const { x, y, vertices } = triangle(current.config.seed), delta = .06;
      if (p[0] <= x * (1-delta) && p[1] <= y * (1-delta) && p[0]/x+p[1]/y >= 1+delta) return p;
      const inset = vertices.map(v => [(1-3*delta)*v[0]+2*delta*x, (1-3*delta)*v[1]+2*delta*y]);
      let nearest = inset[0], distance = Infinity;
      inset.forEach((a, i) => {
        const b = inset[(i+1)%3], dx = b[0]-a[0], dy = b[1]-a[1];
        const t = Math.max(0, Math.min(1, ((p[0]-a[0])*dx+(p[1]-a[1])*dy)/(dx*dx+dy*dy)));
        const q = [a[0]+t*dx, a[1]+t*dy], d = (q[0]-p[0])**2+(q[1]-p[1])**2;
        if (d < distance) { distance = d; nearest = q; }
      });
      return nearest;
    }
    function changeSource(p) {
      if (typeof onSource === "function") onSource(p);
      else { current.source = fallbackProject(p); render(); }
    }
    function drawTriangle(vertices, shift, decorated, inset) {
      const t = triangle(current.config.seed), delta = inset ? .06 : 0;
      const pts = vertices.map(v => coordinates.xy([(1-3*delta)*v[0]+2*delta*t.x+shift[0], (1-3*delta)*v[1]+2*delta*t.y+shift[1]]));
      const color = decorated ? attached : base;
      const attrs = inset
        ? `fill="${color}" fill-opacity=".045" stroke="${color}" stroke-opacity=".6" stroke-width="1"`
        : `fill="none" stroke="${color}" stroke-width="${decorated ? 1.8 : 2.8}"${decorated ? "" : ' stroke-dasharray="7 5"'}`;
      return `<polygon data-shift-triangle="${decorated ? "decorated" : "base"}" points="${pts.map(p=>p.map(fmt).join(",")).join(" ")}" ${attrs}/>`;
    }
    function render() {
      if (!current) return;
      const {config,source} = current, {x,y,vertices} = triangle(config.seed);
      const shift = config.shift.map(n => mod2(2*n)/2), zero = shift.every(n => n === 0);
      const totalPairs = config.graph.stable_edges.reduce((sum, edge) => sum + Number(edge.count || 0), 0);
      const pairCount = `${totalPairs} pair${totalPairs === 1 ? "" : "s"}`;
      // Keep one coordinate frame for all four possible half-shifts. Counts
      // change only the translated triangle, never the picture's scale.
      const maxX = x + .5, maxY = y + .5;
      const scale = Math.min(500/(maxX+.13), 278/(maxY+.08));
      const left = 72 + (500-scale*(maxX+.13))/2, bottom = 329;
      coordinates = {scale,left,bottom,xy:p=>[left+scale*p[0],bottom-scale*p[1]]};
      const xy = coordinates.xy, basePoint = xy(source), decoratedSource = source.map((n,i)=>n+shift[i]), decoratedPoint = xy(decoratedSource);
      const grid = get("shift-grid")?.checked !== false;
      const ticksX = Array.from({length:Math.floor(maxX*2)+1},(_,i)=>i/2);
      const ticksY = Array.from({length:Math.floor(maxY*2)},(_,i)=>(i+1)/2);
      const parts = [`<defs><marker id="${marker}" markerWidth="8" markerHeight="8" refX="6.8" refY="4" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0 L8 4 L0 8" fill="none" stroke="${attached}" stroke-width="1.4"/></marker></defs>`];
      if (grid) {
        for (const n of ticksX.filter(n=>n>0)) parts.push(`<path d="M${fmt(xy([n,0])[0])} ${bottom} V${fmt(xy([0,maxY+.03])[1])}" stroke="${faint}" stroke-dasharray="2 5" fill="none"/>`);
        for (const n of ticksY) parts.push(`<path d="M${fmt(left)} ${fmt(xy([0,n])[1])} H${fmt(xy([maxX+.03,0])[0])}" stroke="${faint}" stroke-dasharray="2 5" fill="none"/>`);
      }
      parts.push(`<path d="M${fmt(left)} ${fmt(xy([0,maxY+.04])[1])} V${bottom} H${fmt(xy([maxX+.06,0])[0])}" fill="none" stroke="#bac2bb" stroke-width="1.2"/>`);
      for (const n of ticksX) parts.push(`<text x="${fmt(xy([n,0])[0])}" y="350" text-anchor="middle" font-size="12" fill="#66716c">${half(n)}</text>`);
      for (const n of ticksY) parts.push(`<text x="${fmt(left-12)}" y="${fmt(xy([0,n])[1]+4)}" text-anchor="end" font-size="12" fill="#66716c">${half(n)}</text>`);
      parts.push(`<text x="${fmt(xy([maxX+.06,0])[0]+16)}" y="333" font-size="14" fill="#46534e">R</text><text x="${fmt(left-3)}" y="${fmt(xy([0,maxY+.04])[1]-12)}" text-anchor="middle" font-size="14" fill="#46534e">S</text>`);
      // Filled compact set first; true boundary and base dashes remain visible.
      parts.push(drawTriangle(vertices,shift,true,true),drawTriangle(vertices,[0,0],false,true),drawTriangle(vertices,[0,0],false,false),drawTriangle(vertices,shift,true,false));
      if (!zero) {
        const dx = decoratedPoint[0]-basePoint[0],dy = decoratedPoint[1]-basePoint[1],d = Math.hypot(dx,dy);
        const start = [basePoint[0]+dx/d*10,basePoint[1]+dy/d*10],end = [decoratedPoint[0]-dx/d*10,decoratedPoint[1]-dy/d*10];
        parts.push(`<path d="M${start.map(fmt).join(",")} L${end.map(fmt).join(",")}" fill="none" stroke="white" stroke-width="5"/><path d="M${start.map(fmt).join(",")} L${end.map(fmt).join(",")}" fill="none" stroke="${attached}" stroke-width="1.9" marker-end="url(#${marker})"/>`);
      }
      parts.push(`<g data-shift-point="base" style="cursor:grab"><title>Undecorated seed: (${source.map(n=>n.toFixed(4)).join(", ")}); drag to change the lengths.</title><circle cx="${fmt(basePoint[0])}" cy="${fmt(basePoint[1])}" r="15" fill="transparent"/><circle cx="${fmt(basePoint[0])}" cy="${fmt(basePoint[1])}" r="7" fill="white" stroke="${base}" stroke-width="2.2"/></g>`);
      parts.push(`<g data-shift-point="decorated" style="cursor:grab"><title>Current decoration, ${pairCount}: (${decoratedSource.map(n=>n.toFixed(4)).join(", ")}); the same length parameters.</title><circle cx="${fmt(decoratedPoint[0])}" cy="${fmt(decoratedPoint[1])}" r="${zero ? 9 : 15}" fill="transparent"/><circle cx="${fmt(decoratedPoint[0])}" cy="${fmt(decoratedPoint[1])}" r="4.2" fill="${attached}" stroke="white" stroke-width="1"/></g>`);
      parts.push(`<text x="320" y="378" text-anchor="middle" font-size="12" fill="${grey}">Current result: ${pairCount} · ${zero ? "no net shift; triangles coincide" : `shift ${pair(shift)} modulo ℤ²`}</text>`);
      svg.innerHTML = parts.join("");
      get("shift-legend") && (get("shift-legend").innerHTML = `<strong class="current-decoration-summary">Current decoration: ${pairCount} applied${Number.isFinite(config.hemispheres) ? ` · ${config.hemispheres} hemispheres in total` : ""}</strong><div class="current-result-legend"><span style="color:${base}"><span aria-hidden="true" style="display:inline-block;width:24px;border-top:3px dashed ${base};vertical-align:middle;margin-right:6px"></span>Undecorated seed</span> <span style="color:${attached}"><span aria-hidden="true" style="display:inline-block;width:24px;border-top:2px solid ${attached};vertical-align:middle;margin-right:6px"></span>Current decoration</span></div>`);
      get("shift-point-readout") && (get("shift-point-readout").innerHTML = `<span style="color:${base}">Seed${math(`<mrow><mo stretchy="false">(</mo><mi>R</mi><mo>,</mo><mi>S</mi><mo stretchy="false">)</mo></mrow><mo>≈</mo>${mathPair(source.map(n=>n.toFixed(4)))}`)}</span><span style="color:${attached}">Current decoration${math(`<mrow><mo stretchy="false">(</mo><msub><mi>R</mi><mtext>cur</mtext></msub><mo>,</mo><msub><mi>S</mi><mtext>cur</mtext></msub><mo stretchy="false">)</mo></mrow><mo>≈</mo>${mathPair(decoratedSource.map(n=>n.toFixed(4)))}`)}</span><span>Current total shift${math(`<mrow><mo stretchy="false">(</mo><mi>Δ</mi><mi>R</mi><mo>,</mo><mi>Δ</mi><mi>S</mi><mo stretchy="false">)</mo></mrow><mo>≡</mo>${mathPair(shift)}${modulo}`)}</span>`);
      renderBreakdown(shift,zero);
    }
    function renderBreakdown(shift,zero) {
      const out = get("shift-breakdown");
      if (!out) return;
      const chains = current.config.graph.stable_edges, sum = [0,0];
      const rows = chains.map((e,i) => {
        const q = e.count || 0, perPair = e.shift_twice_per_pair.map(n=>mod2(n)/2), contribution = e.shift_twice_per_pair.map(n=>mod2(q*n)/2);
        sum[0] += contribution[0];sum[1] += contribution[1];
        return `<tr${i===current.chainIndex ? ' class="is-selected"' : ""}><td>${escape(label(e))} <small>F${e.faces[0]} / F${e.faces[1]}</small></td><td>${q}</td><td>${pair(perPair)}</td><td>${pair(contribution)}</td></tr>`;
      });
      const coreShift = shift.map((n,i)=>mod2(2*(n-sum[i]))/2), extra = coreShift.some(n=>n!==0) ? `<p class="quiet">Core contribution: ${pair(coreShift)}.</p>` : "";
      const selected = chains[current.chainIndex];
      let note = "";
      if (selected) {
        const next = shift.map((n,i)=>mod2(2*n+selected.shift_twice_per_pair[i])/2);
        const changed = selected.shift_twice_per_pair.map(n=>mod2(n));
        const signNote = changed.every(n=>n === 0) ? "Neither ordinary sign would change." : changed.every(n=>n === 1) ? "Both ordinary signs would flip." : `The ordinary sign on period ${changed[0] ? "1" : "τ"} would flip.`;
        note = `<div class="next-pair-note"><strong>If you apply one more pair on ${escape(label(selected))}</strong><p>Current shift ${pair(shift)} → proposed shift ${pair(next)}, modulo ℤ².</p><p>${signNote} This pair has not been applied.</p></div>`;
      }
      const signedCharacter = (period, index) => `${character(period, true)}<mo>=</mo>${shift[index] ? '<mo>−</mo>' : ''}${character(period)}`;
      out.innerHTML = `<p class="shift-total"><strong>Current total shift: ${pair(shift)}, modulo ℤ².</strong></p><p class="quiet"><var>χ</var><sub>0</sub>: undecorated seed; <var>χ</var><sub>cur</sub>: current decoration.</p>${math(signedCharacter("1", 0))}${math(signedCharacter("τ", 1))}<p class="quiet">The current decoration ${zero ? "preserves both ordinary signs" : "changes the ordinary signs shown above"}; the spherical rotations are unchanged.</p>${extra}<table class="shift-contribution-table"><thead><tr><th>Chain / faces</th><th>Applied pairs</th><th>Per pair</th><th>Current shift</th></tr></thead><tbody>${rows.join("")}</tbody></table>${note}`;
    }
    function position(event) {
      const matrix = svg.getScreenCTM();
      if (!matrix || !coordinates) return null;
      const p = svg.createSVGPoint();p.x=event.clientX;p.y=event.clientY;
      const q = p.matrixTransform(matrix.inverse());
      return [(q.x-coordinates.left)/coordinates.scale,(coordinates.bottom-q.y)/coordinates.scale];
    }
    svg.addEventListener("pointerdown", event => {
      if (!current || event.button !== 0) return;
      const target = event.target.closest("[data-shift-point]");
      if (!target) return;
      event.preventDefault();svg.focus({preventScroll:true});dragging = target.dataset.shiftPoint;
      svg.setPointerCapture(event.pointerId);
    });
    svg.addEventListener("pointermove", event => {
      if (!dragging) return;
      const p = position(event);
      if (p) changeSource(p.map((n,i)=>n-(dragging==="decorated" ? current.config.shift[i] : 0)));
    });
    const end = () => { dragging=null; };
    svg.addEventListener("pointerup", end);svg.addEventListener("pointercancel", end);svg.addEventListener("lostpointercapture", end);
    svg.addEventListener("keydown", event => {
      if (!current) return;
      const direction = {ArrowLeft:[-1,0],ArrowRight:[1,0],ArrowUp:[0,1],ArrowDown:[0,-1]}[event.key];
      if (!direction) return;
      event.preventDefault();const step=event.shiftKey?.04:.01;
      changeSource(current.source.map((n,i)=>n+step*direction[i]));
    });
    get("shift-grid")?.addEventListener("change", render);
    return {update(options) {
      if (!options?.config || !Array.isArray(options.source) || !options.source.every(Number.isFinite)) throw new Error("A metric source and decorated graph are needed.");
      current = {...options,source:options.source.slice()};render();
    }};
  }
  return {create};
})();

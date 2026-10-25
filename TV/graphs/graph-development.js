"use strict";

/* A common developing chart, rather than a planar model of the glued surface.
   A spherical half-turn with equatorial axis Q_j restricts to
   exp(it) -> exp(i(2 theta_j-t)) on this equator. */
globalThis.GraphDevelopment = (() => {
  const PI = Math.PI, TAU = 2 * PI;
  const ink = "#25353d", muted = "#6f777a", line = "#d8deda", accent = "#b6403b", reference = "#315d83";
  const subs = ["₀", "₁", "₂", "₃"];
  const esc = x => String(x).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll('"', "&quot;");
  const number = (x, digits = 5) => Math.abs(x) < .5 * 10 ** -digits ? "0" : Number(x.toFixed(digits)).toString().replace("-", "−");
  const principal = x => Math.atan2(Math.sin(x), Math.cos(x));
  const complex = z => `${number(z.re)} ${z.im < 0 ? "−" : "+"} ${number(Math.abs(z.im))}i`;
  const math = body => `<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">${body}</math>`;
  const sub = (symbol, index) => `<msub><mi>${symbol}</mi><mn>${index}</mn></msub>`;
  const num = x => `<mn>${number(x)}</mn>`;
  const phaseNumber = x => `<mrow>${num(x)}<mi>π</mi></mrow>`;
  const exponential = body => `<mrow><mi mathvariant="normal">exp</mi><mo stretchy="false">(</mo>${body}<mo stretchy="false">)</mo></mrow>`;
  const complexMath = z => `<mrow>${num(z.re)}<mo>${z.im < 0 ? "−" : "+"}</mo>${num(Math.abs(z.im))}<mi>i</mi></mrow>`;
  const scientificMath = x => {
    const [coefficient, exponent] = x.toExponential(2).split("e");
    return `<mrow><mn>${coefficient}</mn><mo>×</mo><msup><mn>10</mn><mn>${Number(exponent)}</mn></msup></mrow>`;
  };
  const ordinaryPair = i => `<mrow>${sub("γ", i)}${sub("γ", 3)}</mrow>`;
  const wordMath = word => word.terms.length ? word.terms.map(term => {
    const pair = `<mrow><mo stretchy="false">(</mo>${ordinaryPair(term.i)}<mo stretchy="false">)</mo></mrow>`;
    return term.power === 1 ? pair : `<msup>${pair}${num(term.power)}</msup>`;
  }).join("") : '<mi mathvariant="normal">I</mi>';
  function lengthSeparation(seed, first, second) {
    const [a,b,c] = first.affinePi.map((value, j) => value - second.affinePi[j]);
    const coefficients = seed === "bouquet" ? [a+c/2,a+b/2,a+(b+c)/2] : seed === "tetra" ? [(a+c)/2,(a+b)/2,(a+b+c)/2] : [a/2+c/4,a/2+b/2+c/4,a+b];
    const names = seed === "double" ? ["x","y","L"] : ["x","y","z"];
    const terms = [], labels = [];
    coefficients.forEach((coefficient, j) => {
      if (!coefficient) return;
      const value = Math.abs(coefficient);
      let numerator = Math.round(value*4), denominator = 4;
      while (denominator > 1 && numerator % 2 === 0) {numerator /= 2; denominator /= 2;}
      const factor = value === 1 ? "" : denominator === 1 ? `<mn>${numerator}</mn>` : `<mfrac><mn>${numerator}</mn><mn>${denominator}</mn></mfrac>`;
      const sign = coefficient < 0 ? "−" : terms.length ? "+" : "";
      terms.push(`${sign ? `<mo>${sign}</mo>` : ""}${factor}<mi>${names[j]}</mi>`);
      labels.push(`${sign}${value === 1 ? "" : denominator === 1 ? numerator : `(${numerator}/${denominator})`}${names[j]}`);
    });
    return {math:terms.join("") || '<mn>0</mn>', label:labels.join(" ") || "0"};
  }
  function periodLengthAndSign(seed, holonomy, word) {
    const affinePi = [0,0,0], weights = [0,0,0,0];
    for (const term of word.terms) {
      for (let j=0; j<3; j++) affinePi[j] += term.power*(holonomy.axes[term.i].affinePi[j]-holonomy.axes[3].affinePi[j]);
      weights[term.i] += term.power;
      weights[3] += term.power;
    }
    const length = lengthSeparation(seed,{affinePi},{affinePi:[0,0,0]});
    const terms = [];
    weights.forEach((coefficient,i) => {
      if (Math.abs(coefficient) % 2 === 0) return;
      terms.push(`${terms.length ? '<mo>+</mo>' : ''}${sub("n",i)}`);
    });
    const sign = terms.length ? `<msup><mrow><mo stretchy="false">(</mo><mo>−</mo><mn>1</mn><mo stretchy="false">)</mo></mrow><mrow>${terms.join("")}</mrow></msup>` : '';
    const negative = length.math.startsWith('<mo>−</mo>') && !length.math.includes('<mo>+</mo>');
    const argument = negative ? length.math.slice('<mo>−</mo>'.length).replaceAll('<mo>−</mo>','<mo>+</mo>') : length.math;
    const grouped = /^<mi>[xyzL]<\/mi>$/.test(argument) ? argument : `<mo stretchy="false">(</mo>${argument}<mo stretchy="false">)</mo>`;
    return `${sign}<msup><mi>e</mi><mrow>${negative ? '<mo>−</mo>' : ''}<mi>i</mi>${grouped}</mrow></msup>`;
  }
  const Q = i => `Q${subs[i]}`;
  const gamma = i => `γ${subs[i]}`;
  const point = (a, r, cx = 220, cy = 170) => [cx + r * Math.cos(a), cy - r * Math.sin(a)];
  const xy = p => p.map(x => Number(x.toFixed(3))).join(" ");
  const text = (x, y, content, extra = "") => `<text x="${x}" y="${y}" ${/\bfill=/.test(extra) ? "" : `fill="${ink}"`} font-family="system-ui,sans-serif" ${/\bfont-size=/.test(extra) ? "" : 'font-size="12"'} ${extra}>${esc(content)}</text>`;
  const svgCircle = (a, r, radius, fill, extra = "", cx = 220, cy = 170) => {
    const [x, y] = point(a, r, cx, cy);
    return `<circle cx="${x}" cy="${y}" r="${radius}" fill="${fill}" ${extra}/>`;
  };
  function arcPath(from, change, radius, cx = 220, cy = 170, principalArc = true) {
    if (Math.abs(change) < 1e-10) return "";
    const d = principalArc ? principal(change) : change, a = point(from, radius, cx, cy), b = point(from + d, radius, cx, cy);
    return `M ${xy(a)} A ${radius} ${radius} 0 ${Math.abs(d) > PI ? 1 : 0} ${d > 0 ? 0 : 1} ${xy(b)}`;
  }
  function trace(holonomy, cycle = "a") {
    const key = cycle === "b" ? "tau" : "one", word = holonomy.selectedwords[key];
    const steps = [];
    // The rightmost matrix acts first. Each positive γ_i γ_3 factor
    // consequently applies h_3 and then h_i to the developed value.
    for (const term of word.terms.slice().reverse()) {
      for (let k = 0; k < Math.abs(term.power); k++) {
        const poles = term.power > 0 ? [3, term.i] : [term.i, 3];
        for (const j of poles) steps.push({pole: j, term, copy: k});
      }
    }
    const initial = .23 * PI, states = [initial];
    for (const step of steps) states.push(2 * holonomy.axes[step.pole].phase - states.at(-1));
    const phase = states.at(-1) - initial;
    const endpoint = {re: Math.cos(phase), im: Math.sin(phase)};
    const expected = holonomy.projective[key];
    const error = Math.hypot(endpoint.re - expected.re, endpoint.im - expected.im);
    return {key, word, steps, states, initial, phase, endpoint, expected, error};
  }
  function labelPositions(axes) {
    const placed = [];
    return axes.map(axis => {
      const label = Q(axis.i), width = 76, height = 34;
      let best = null;
      for (const r of [150, 170, 190]) for (const shift of [0, .15, -.15, .3, -.3, .45, -.45, .6, -.6]) {
        const [px, py] = point(axis.phase + shift, r);
        const cx = Math.max(width / 2 + 6, Math.min(434 - width / 2, px));
        const cy = Math.max(34, Math.min(270, py));
        const box = {left: cx - width / 2, right: cx + width / 2, top: cy - height / 2, bottom: cy + height / 2};
        const overlap = placed.reduce((sum, p) => sum + Math.max(0, Math.min(box.right, p.right + 6) - Math.max(box.left, p.left - 6)) * Math.max(0, Math.min(box.bottom, p.bottom + 4) - Math.max(box.top, p.top - 4)), 0);
        const nearestX = Math.max(box.left,Math.min(box.right,220)), nearestY = Math.max(box.top,Math.min(box.bottom,170));
        const circleIntrusion = Math.max(0,128-Math.hypot(nearestX-220,nearestY-170));
        const score = overlap * 100 + circleIntrusion**2 * 200 + Math.abs(shift) * 70 + Math.abs(r - 150);
        if (!best || score < best.score) best = {cx, cy, box, score, label};
      }
      placed.push(best.box);
      return best;
    });
  }
  function create({onPole = () => {}, onArc = () => {}, onViewChange = () => {}} = {}) {
    const get = id => document.getElementById(id);
    const equator = get("developed-equator"), sequence = get("halfturn-sequence"), slider = get("development-step");
    const prev = get("development-prev"), next = get("development-next"), start = get("development-start");
    let current = null, step = 0, caseKey = "", selectedArc = null;
    const requestId = request => !request ? "" : request.kind === "period" ? `period:${request.cycle}` : request.kind === "axis-location" ? `axis-location:${request.i}` : `${request.kind}:${request.i}:${request.reversed ? 1 : 0}`;
    const requestFromId = id => {
      const [kind,value,reversed] = id.split(":");
      return kind === "period" ? {kind,cycle:value} : kind === "axis-location" ? {kind,i:Number(value)} : {kind,i:Number(value),reversed:reversed === "1"};
    };
    const arcHosts = () => [equator,sequence,get("development-readout"),get("halfturn-readout")].filter(Boolean);
    const requestLabel = request => request.kind === "period" ? `Period ${request.cycle === "b" ? "τ" : "1"} rotation` : request.kind === "axis-location" ? `Axis h${subs[request.i]}` : `${request.reversed ? "−" : ""}${request.kind === "pair" ? "2" : ""}L${subs[request.i]}`;
    function visibleArcRequests() {
      const found = new Map();
      for (const host of arcHosts()) for (const element of host.querySelectorAll("[data-development-arc]")) {
        const id = element.dataset.developmentArc;
        if (!found.has(id)) {
          const request = requestFromId(id);
          found.set(id,{id,request,label:requestLabel(request)});
        }
      }
      return [...found.values()];
    }
    function highlightEdgeMatches(matches = [],edgeLabel = "") {
      const coefficients = new Map(matches.map(match => [match.id,match.coefficient]));
      for (const host of arcHosts()) for (const element of host.querySelectorAll("[data-development-arc]")) {
        element.classList.remove("is-edge-match");
        delete element.dataset.edgeCoefficient;
        element.querySelectorAll(".edge-arc-coefficient").forEach(badge => badge.remove());
        element.removeAttribute("aria-description");
        const coefficient = coefficients.get(element.dataset.developmentArc);
        if (!Number.isFinite(coefficient) || coefficient === 0) continue;
        element.classList.add("is-edge-match");
        element.dataset.edgeCoefficient = String(coefficient);
        const signed = `${coefficient < 0 ? "−" : "+"}${Math.abs(coefficient)}`;
        element.setAttribute("aria-description",`${edgeLabel} contributes ${signed} times its length to this ${element.dataset.developmentArc.startsWith("axis-location:") ? "axis location" : "arc"}.`);
        if (element.tagName.toLowerCase() === "button") {
          const badge = document.createElement("span");
          badge.className = "edge-arc-coefficient";
          badge.textContent = `${signed} × ${edgeLabel}`;
          badge.setAttribute("aria-hidden","true");
          element.append(badge);
        }
      }
    }
    function arcAttributes(request,label,svg = true) {
      const selected = requestId(selectedArc) === requestId(request);
      return `class="development-arc-button${selected ? " is-selected" : ""}" data-development-arc="${requestId(request)}" ${svg ? 'role="button" tabindex="0"' : 'type="button"'} aria-label="${esc(label)}" aria-pressed="${selected}"`;
    }
    function arcMarkup(path,request,label,{marker,width = 3,extra = ""} = {}) {
      const selected = requestId(selectedArc) === requestId(request);
      return `<g ${arcAttributes(request,label)} style="cursor:pointer"><title>${esc(label)}</title>${path ? `<path class="development-arc-hit" d="${path}" fill="none" stroke="transparent" stroke-width="19" pointer-events="stroke"/><path class="development-arc-line" d="${path}" fill="none" stroke="${accent}" stroke-width="${selected ? width+1.5 : width}" ${marker ? `marker-end="url(#${marker})"` : ''} pointer-events="none"/>` : ''}${extra}</g>`;
    }
    const arcLink = (request,label) => `<button ${arcAttributes(request,label,false)}>${label}</button>`;
    function renderGuide(path) {
      const guide = get("halfturn-guide");
      if (!guide) return;
      const reduced = [], cancelled = new Set(), partners = new Map(), cancellations = [];
      path.steps.forEach((item,index) => {
        const previous = reduced.at(-1);
        if (previous !== undefined && path.steps[previous].pole === item.pole) {
          reduced.pop();cancelled.add(previous+1);cancelled.add(index+1);partners.set(previous+1,index+1);partners.set(index+1,previous+1);cancellations.push([previous,index]);
        } else reduced.push(index);
      });
      const reducedWord = reduced.length ? reduced.map(index => `h${subs[path.steps[index].pole]}`).join(" → ") : "the identity";
      const cancellationCopy = cancellations.length === 1 && cancellations[0][1] === cancellations[0][0]+1 ? `The two consecutive h${subs[path.steps[cancellations[0][0]].pole]} half-turns cancel; this period is equivalent to <span class="halfturn-reduced-word">${reducedWord}</span>.` : `After adjacent cancellations, this period is represented by <span class="halfturn-reduced-word">${reducedWord}</span>.`;
      let copy = guide.querySelector(".halfturn-guide-copy");
      if (!copy) {
        copy = document.createElement("div");copy.className = "halfturn-guide-copy";guide.prepend(copy);
      }
      copy.innerHTML = `<p><strong>Locate the marked points.</strong> Place one hemisphere on the sphere and continue its chart across the glued edges. The boundary lengths locate <var>Q</var><sub>i</sub>, the image in this continued chart of each marked vertex <var>p</var><sub>i</sub> of the quotient sphere.</p><p><strong>Draw its half-turn axis.</strong> A loop around <var>p</var><sub>i</sub> gives <var>h</var><sub>i</sub>: a 180° rotation about the diameter through <var>Q</var><sub>i</sub> and its antipode. On the equator this acts as reflection across that diameter.</p><p><strong>Follow the selected period.</strong> Each pair uses <var>h</var><sub>3</sub> and one <var>h</var><sub>i</sub>, in the indicated order. This case has <strong>${path.steps.length/2} ${path.steps.length===2?"pair":"pairs"}, or ${path.steps.length} displayed half-turns</strong>. Read the chips from left to right.</p>${cancellations.length ? `<p class="halfturn-cancel-note">${cancellationCopy} We show the expanded pairs so their length contributions remain visible.</p>` : '<p>Two half-turns give a rotation; composing all the pairs gives the period rotation. Choose Start to inspect the local moves.</p>'}`;
      let steps = get("halfturn-steps");
      if (!steps) {steps = document.createElement("div");steps.id="halfturn-steps";guide.append(steps);}
      steps.setAttribute("aria-label","Half-turns in their order of application");
      const chip = (index,label,complete = false) => `<button type="button" class="halfturn-step-chip ${index === step ? 'is-active' : index < step ? 'is-completed' : 'is-pending'}${index === 0 ? ' is-start' : ''}${complete ? ' is-final' : ''}${cancelled.has(index) ? ' is-cancelled' : ''}" data-development-step="${index}" aria-current="${index === step ? 'step' : 'false'}" aria-label="${index === 0 ? 'Start, before any half-turn' : `Apply through half-turn ${index} of ${path.steps.length}: h${path.steps[index-1].pole}${complete ? ', complete period' : ''}${cancelled.has(index) ? `; cancels with half-turn ${partners.get(index)} in the expanded word` : ''}`}" aria-pressed="${index === step}"><span>${label}</span>${complete ? '<small>Complete</small>' : ''}${cancelled.has(index) ? '<small class="halfturn-cancel-badge">cancel</small>' : ''}</button>`;
      steps.innerHTML = chip(0,"Start") + path.steps.map((item,index) => `<span class="halfturn-step-arrow" aria-hidden="true">→</span>${chip(index+1,`h${subs[item.pole]}`,index === path.steps.length-1)}`).join("");
    }
    function renderEquator(active) {
      const {holonomy: h, path} = current, axes = h.axes;
      const pair = active?.term ?? path.word.terms[0] ?? {i: 0, power: 0};
      const i = pair.i, theta3 = axes[3].phase, diff = axes[i].phase - theta3;
      const length = lengthSeparation(current.config.seed, axes[i], axes[3]);
      const reversed = pair.power < 0;
      const arcStart = reversed ? axes[i].phase : theta3, arcChange = reversed ? -diff : diff;
      const arcName = `${reversed ? "−" : ""}L${subs[i]}`;
      const arcLength = reversed ? lengthSeparation(current.config.seed, axes[3], axes[i]) : length;
      const separationMath = `${reversed ? '<mo>−</mo>' : ''}${sub("L",i)}<mo>=</mo>${arcLength.math}`;
      const labels = labelPositions(axes);
      let html = `<title>The developed equator and the axes of the local monodromies</title><desc>The points Q zero through Q three are obtained by continuing hemisphere charts through the chosen gluing tree. Each diameter is the axis of a spherical half-turn. The oriented arc represents the actual signed angular separation from Q three to the selected point, including a major arc when required by the edge lengths.</desc><defs><marker id="development-arc-arrow" markerUnits="userSpaceOnUse" markerWidth="9" markerHeight="9" refX="7.3" refY="4.5" orient="auto"><path d="M 0 0 L 9 4.5 L 0 9 Z" fill="${accent}"/></marker></defs>`;
      html += `<circle cx="220" cy="170" r="100" fill="none" stroke="${line}" stroke-width="1.5"/><circle cx="220" cy="170" r="2" fill="${muted}"/>`;
      html += text(229, 184, "O", `font-size="11" fill="${muted}"`);
      axes.forEach(axis => {
        const relevant = axis.i === i || axis.i === 3, activeAxis = step > 0 && axis.i === active?.pole, a = point(axis.phase, 100), b = point(axis.phase + PI, 100);
        html += `<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${relevant ? reference : muted}" stroke-width="${activeAxis ? 2.6 : relevant ? 1.8 : 1.2}" ${relevant ? "" : 'stroke-dasharray="3 4" opacity=".6"'}/>`;
        html += svgCircle(axis.phase + PI, 100, 2.6, "white", `stroke="${relevant ? reference : muted}"`);
        if (relevant) {
          const label = point(axis.phase + PI, axis.i === 3 ? 64 : 43);
          html += text(label[0], label[1], `h${subs[axis.i]}`, `text-anchor="middle" font-size="13" fill="${reference}" style="paint-order:stroke;stroke:white;stroke-width:4px"`);
        }
      });
      const arc = arcPath(arcStart, arcChange, 116,220,170,false);
      const lengthLabel = point(arcStart + arcChange / 2, 77);
      const axisRequest = {kind:"axis",i,reversed};
      html += arcMarkup(arc,axisRequest,`Highlight the edges contributing to ${reversed ? 'minus ' : ''}L${i}`,{marker:"development-arc-arrow",extra:text(lengthLabel[0], lengthLabel[1], arcName, `text-anchor="middle" font-size="16" fill="${accent}" style="paint-order:stroke;stroke:white;stroke-width:5px"`)});
      axes.forEach((axis, index) => {
        const label = labels[index], p = point(axis.phase, 100), relevant = axis.i === i || axis.i === 3;
        html += `<g data-development-pole="${axis.i}" role="button" tabindex="0" aria-label="Select marked point p${axis.i}, with weight n${axis.i} equal to ${axis.weight}" style="cursor:pointer">`;
        html += `<line x1="${p[0]}" y1="${p[1]}" x2="${label.cx}" y2="${label.cy}" stroke="${line}" stroke-width="1"/>`;
        html += svgCircle(axis.phase, 100, 13, "transparent");
        html += svgCircle(axis.phase, 100, relevant ? 5.3 : 4.2, relevant ? reference : muted, 'stroke="white" stroke-width="1.5"');
        html += `<rect x="${label.box.left}" y="${label.box.top}" width="76" height="34" rx="3" fill="white" fill-opacity=".95"/>`;
        html += text(label.cx, label.cy - 2, label.label, `text-anchor="middle" font-size="14" fill="${relevant ? reference : muted}"`);
        html += text(label.cx, label.cy + 13, `h${subs[axis.i]}`, `text-anchor="middle" font-size="12" fill="${relevant ? reference : muted}"`);
        html += "</g>";
      });
      html += text(220, 317, `${reversed ? `h${subs[i]} then h${subs[3]}` : `h${subs[3]} then h${subs[i]}`}: rotation by ${reversed ? "−" : ""}2L${subs[i]}`, `text-anchor="middle" font-size="15" fill="${accent}"`);
      equator.innerHTML = html;
      get("axis-length-formula").innerHTML = math(`<mrow>${separationMath}</mrow>`);
      const signFormula = `<mrow>${sub("L",i)}<mo>=</mo>${sub("θ",i)}<mo>−</mo>${sub("θ",3)}<mo>,</mo><mspace width="1em"/>${sub("λ",i)}<mo>=</mo><msup><mrow><mo stretchy="false">(</mo><mo>−</mo><mn>1</mn><mo stretchy="false">)</mo></mrow><mrow>${sub("n",i)}<mo>+</mo>${sub("n",3)}</mrow></msup>${exponential(`<mi>i</mi>${sub("L",i)}`)}<mo>.</mo></mrow>`;
      const chartDetails = `<details class="calculation-details"><summary>Axes and ordinary signs</summary><p>Gluing transports the boundary angle <var>α</var><sub>i</sub> by a sign <var>ε</var><sub>i</sub> and offset <var>β</var><sub>i</sub>; the ordinary lift also retains the weight sign.</p>${math(`<mrow><msub><mi>θ</mi><mi>i</mi></msub><mo>=</mo><msub><mi>ε</mi><mi>i</mi></msub><msub><mi>α</mi><mi>i</mi></msub><mo>+</mo><msub><mi>β</mi><mi>i</mi></msub><mo>,</mo><mspace width="1em"/><msub><mi>Q</mi><mi>i</mi></msub><mo>=</mo>${exponential('<mi>i</mi><msub><mi>θ</mi><mi>i</mi></msub>')}<mo>.</mo></mrow>`)}${math(signFormula)}</details>`;
      get("development-readout").innerHTML = `<p>The red arc gives the signed separation. Two half-turns double its angle.</p><div class="development-arc-links">${arcLink(axisRequest,arc ? "Highlight the arc’s contributing edges" : "Zero separation: highlight the contributing edges")}</div>${chartDetails}`;
    }
    function renderSequence() {
      if (!current || !sequence || !equator) return;
      const {holonomy: h, path, cycle} = current, N = path.steps.length;
      step = Math.max(0, Math.min(N, step));
      const active = path.steps[Math.max(0, step - 1)] ?? path.steps[0];
      renderGuide(path);
      const detailsOpen=get("development-readout")?.querySelector("details")?.open;
      renderEquator(active);
      if(detailsOpen){const details=get("development-readout").querySelector("details");if(details)details.open=true;}
      const period = cycle === "b" ? "τ" : "1", finished = step === N;
      let html = `<title>Successive developed images under the selected period word</title><desc>The open marker denotes the initial developed value. The filled marker denotes its image under the half-turns applied so far. At the final step it is the developed value after continuation along the selected marked period. The spherical multiplier is the square of the ordinary period character.</desc><defs><marker id="development-probe-arrow" markerUnits="userSpaceOnUse" markerWidth="9" markerHeight="9" refX="7.3" refY="4.5" orient="auto"><path d="M 0 0 L 9 4.5 L 0 9 Z" fill="${accent}"/></marker></defs><circle cx="220" cy="170" r="102" fill="none" stroke="${line}" stroke-width="1.5"/>`;
      if (active) {
        const theta = h.axes[active.pole].phase, a = point(theta, 102), b = point(theta + PI, 102);
        html += `<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${reference}" stroke-width="1.8"/>`;
        html += svgCircle(theta, 102, 4, reference);
        const label = point(theta, 67);
        html += text(label[0], label[1], `h${subs[active.pole]}`, `text-anchor="middle" font-size="14" fill="${reference}" style="paint-order:stroke;stroke:white;stroke-width:4px"`);
      }
      for (let j = 0; j < path.states.length; j++) {
        if (j !== step) html += svgCircle(path.states[j], 102, 2.2, muted, `opacity="${j < step ? ".34" : ".13"}"`);
      }
      const p0 = point(path.initial, 102), cur = point(path.states[step], 102);
      html += `<circle cx="${p0[0]}" cy="${p0[1]}" r="6.5" fill="white" stroke="${reference}" stroke-width="2"/>`;
      if (step) {
        const before = path.states[step - 1], arc = arcPath(before, path.states[step] - before, 110);
        if (arc) html += !finished && step % 2 === 1 ? arcMarkup(arc,{kind:"axis-location",i:active.pole},`Highlight the edge lengths locating the h${active.pole} axis`,{marker:"development-probe-arrow",width:2.5}) : `<path d="${arc}" fill="none" stroke="${line}" stroke-width="1.4"/>`;
        const p = point(before, 102);
        html += `<circle cx="${p[0]}" cy="${p[1]}" r="4.5" fill="white" stroke="${muted}" stroke-width="1.3"/>`;
      }
      html += `<line x1="220" y1="170" x2="${cur[0]}" y2="${cur[1]}" stroke="${accent}" stroke-width="1" stroke-dasharray="3 4"/><circle cx="${cur[0]}" cy="${cur[1]}" r="5.3" fill="${accent}" stroke="white" stroke-width="1.5"/>`;
      if (step > 0 && step % 2 === 0 && !finished) {
        const change = path.states[step] - path.states[step - 2];
        const pairArc = arcPath(path.states[step - 2], change, 121);
        const pairLabel = point(path.states[step - 2] + principal(change)/2, 76);
        html += arcMarkup(pairArc,{kind:"pair",i:active.term.i,reversed:active.term.power < 0},`Highlight the edges contributing to this pair rotation`,{marker:"development-probe-arrow",extra:text(pairLabel[0],pairLabel[1],`${active.term.power < 0 ? "−" : ""}2L${subs[active.term.i]}`,`text-anchor="middle" font-size="15" fill="${accent}" style="paint-order:stroke;stroke:white;stroke-width:5px"`)});
        const start = point(path.states[step - 2], 102);
        html += `<circle cx="${start[0]}" cy="${start[1]}" r="4.5" fill="white" stroke="${muted}" stroke-width="1.3"/>`;
      }
      if (finished) {
        const netArc = arcPath(path.initial, path.phase, 122);
        const netLabel = point(path.initial + principal(path.phase)/2, 73);
        html += arcMarkup(netArc,{kind:"period",cycle},`Highlight the edges contributing to marked period ${period}`,{marker:"development-probe-arrow",extra:text(netLabel[0],netLabel[1],`χ(${period})²`,`text-anchor="middle" font-size="15" fill="${accent}" style="paint-order:stroke;stroke:white;stroke-width:5px"`)});
      }
      const occupied = [];
      function probeLabel(angles, label, color) {
        const width = Math.max(47, label.length * 6.7 + 14), height = 23;
        const central = Math.atan2(angles.reduce((s,a)=>s+Math.sin(a),0),angles.reduce((s,a)=>s+Math.cos(a),0));
        let best;
        for (const r of [151, 168, 137]) for (const offset of [0,.2,-.2,.38,-.38,.6,-.6]) {
          const p = point(central + offset, r);
          const x = Math.max(width/2+8,Math.min(432-width/2,p[0])), y = Math.max(27,Math.min(285,p[1]));
          const box = {left:x-width/2,right:x+width/2,top:y-height/2,bottom:y+height/2};
          const overlap = occupied.reduce((sum,b)=>sum+Math.max(0,Math.min(box.right,b.right+8)-Math.max(box.left,b.left-8))*Math.max(0,Math.min(box.bottom,b.bottom+5)-Math.max(box.top,b.top-5)),0);
          const score=overlap*100+Math.abs(offset)*80+Math.abs(r-151);
          if (!best||score<best.score) best={x,y,box,score};
        }
        occupied.push(best.box);
        for(const a of angles){const p=point(a,102);html+=`<line x1="${p[0]}" y1="${p[1]}" x2="${best.x}" y2="${best.y}" stroke="${line}" stroke-width="1"/>`;}
        html += `<rect x="${best.box.left}" y="${best.box.top}" width="${width}" height="${height}" rx="3" fill="white" fill-opacity=".96"/>`;
        html += text(best.x,best.y+4,label,`text-anchor="middle" font-size="12" fill="${color}"`);
      }
      const subStep = String(step).replace(/\d/g,d=>"₀₁₂₃₄₅₆₇₈₉"[Number(d)]), currentLabel=finished?`f(z+${period})`:`P${subStep}`;
      const proximity=Math.hypot(cur[0]-p0[0],cur[1]-p0[1]);
      if(proximity<15) probeLabel([path.initial,path.states[step]],`f(z) ${proximity<1e-7?"=":";"} ${currentLabel}`,accent);
      else {probeLabel([path.initial],"f(z)",reference);probeLabel([path.states[step]],currentLabel,accent);}
      if(finished){
        html += text(220,312,`Period ${period}: f(z) → f(z+${period})`, `text-anchor="middle" font-size="15" fill="${accent}"`);
        html += text(220,337,`Total spherical multiplier χ(${period})²`, `text-anchor="middle" font-size="12" fill="${muted}"`);
      }else{
        const pairCaption = active ? `${active.term.power < 0 ? "−" : ""}2L${subs[active.term.i]}` : "";
        html += text(220, 312, step ? step % 2 === 0 ? `Two half-turns → rotation by ${pairCaption}` : `h${subs[active.pole]}: first half-turn of the pair` : "Start at the developed value f(z)", `text-anchor="middle" font-size="14" fill="${step && step % 2 === 0 ? accent : ink}"`);
        html += text(220, 337, "Follow the chips from left to right", `text-anchor="middle" font-size="11" fill="${muted}"`);
      }
      sequence.innerHTML = html;
      if (slider) {
        slider.max = String(N); slider.value = String(step);
        slider.setAttribute("aria-valuetext", step ? `Step ${step} of ${N}: spherical half-turn h${active.pole}` : `Initial value, before ${N} spherical half-turns`);
      }
      if (prev) prev.disabled = step === 0;
      if (next) next.disabled = step === N;
      if (start) start.disabled = step === 0;
      get("development-step-label").textContent = finished ? `Complete · ${N} / ${N}` : step ? `Step ${step} / ${N} · h${subs[active.pole]}` : `Start · 0 / ${N}`;
      const periodMath = cycle === "b" ? '<mi>τ</mi>' : '<mn>1</mn>';
      const character = `<mrow><mi>χ</mi><mo stretchy="false">(</mo>${periodMath}<mo stretchy="false">)</mo></mrow>`;
      const eigenvalueWord = path.word.terms.length ? path.word.terms.map(term => term.power === 1 ? sub("λ",term.i) : `<msup>${sub("λ",term.i)}${num(term.power)}</msup>`).join("") : '<mn>1</mn>';
      const lengthAndSign = periodLengthAndSign(current.config.seed,h,path.word);
      const caption = finished ? `The red point is the final developed value <var>f</var>(<var>z</var> + ${cycle === "b" ? '<var>τ</var>' : '1'}). The outer arrow shows the complete period’s principal angular increment.` : step ? `The red point <var>P</var><sub>${step}</sub> is an intermediate image of the blue initial value <var>f</var>(<var>z</var>). The chips show the moves in application order.` : 'The blue and red points coincide at the initial value <var>f</var>(<var>z</var>). Follow the chips from left to right.';
      const rightRequest = finished ? {kind:"period",cycle} : step > 0 && step % 2 === 0 ? {kind:"pair",i:active.term.i,reversed:active.term.power < 0} : {kind:"axis-location",i:active.pole};
      const rightLink = arcLink(rightRequest,step === 0 ? "Highlight the next half-turn’s axis lengths" : finished ? "Highlight the period’s contributing edges" : step % 2 === 0 ? "Highlight this pair’s contributing edges" : "Highlight this half-turn’s axis lengths");
      get("period-length-formula").innerHTML = math(`<mrow>${character}<mo>=</mo>${lengthAndSign}</mrow>`);
      get("halfturn-readout").innerHTML = `<p>${caption}</p><div class="development-arc-links">${rightLink}</div><details class="calculation-details"><summary>The selected ordinary monodromy word</summary>${math(`<mrow><msub><mi>W</mi>${periodMath}</msub><mo>=</mo>${wordMath(path.word)}<mo>,</mo></mrow>`)}${math(`<mrow>${character}<mo>=</mo>${eigenvalueWord}<mo>.</mo></mrow>`)}</details>`;
      onViewChange();
    }
    const selectPole = event => {
      if (event.type === "keydown" && !["Enter", " "].includes(event.key)) return;
      const target = event.target.closest?.("[data-development-pole]");
      if (!target) return;
      event.preventDefault();
      onPole(Number(target.dataset.developmentPole));
    };
    const selectArc = event => {
      if (event.type === "keydown" && !["Enter"," "].includes(event.key)) return;
      const target = event.target.closest?.("[data-development-arc]");
      if (!target) return;
      event.preventDefault();event.stopPropagation();
      const id = target.dataset.developmentArc, host = event.currentTarget;
      selectedArc = requestFromId(id);
      renderSequence();
      onArc({...selectedArc});
      if (event.type === "keydown") host.querySelector(`[data-development-arc="${id}"]`)?.focus({preventScroll:true});
    };
    function changeStep(index) {
      const hadSelection = selectedArc !== null;
      step = index;selectedArc = null;renderSequence();
      if (hadSelection) onArc(null);
    }
    const selectStep = event => {
      const target = event.target.closest?.("[data-development-step]");
      if (!target || !current) return;
      event.preventDefault();
      const index = Number(target.dataset.developmentStep);
      changeStep(index);
      get("halfturn-steps")?.querySelector(`[data-development-step="${index}"]`)?.focus({preventScroll:true});
    };
    equator?.addEventListener("click", selectPole);
    equator?.addEventListener("keydown", selectPole);
    for (const element of [equator,sequence,get("development-readout"),get("halfturn-readout")]) {
      element?.addEventListener("click",selectArc);
      element?.addEventListener("keydown",selectArc);
    }
    get("halfturn-guide")?.addEventListener("click",selectStep);
    get("halfturn-steps")?.addEventListener("click",selectStep);
    slider?.addEventListener("input", () => changeStep(Number(slider.value)));
    prev?.addEventListener("click", () => changeStep(step-1));
    next?.addEventListener("click", () => changeStep(step+1));
    start?.addEventListener("click", () => changeStep(0));
    return {
      visibleArcRequests,
      highlightEdgeMatches,
      update({config, source, holonomy, cycle = "a", arcRequest = undefined}) {
        if (!holonomy || !equator || !sequence) return;
        const path = trace(holonomy, cycle);
        // Pair counts alter lift signs, while the transported spherical axes
        // and word stay fixed. Keep the current development step for comparison.
        const key = `${config.seed}|${cycle}|${path.word.exponents.join(",")}`;
        if (key !== caseKey) {step = path.steps.length; caseKey = key;selectedArc = null;}
        if (arcRequest !== undefined) selectedArc = arcRequest;
        current = {config, source, holonomy, cycle, path};
        renderSequence();
      }
    };
  }
  return {create, trace};
})();

"use strict";

/* Metric readout for an exact labelled ribbon graph. Screen embeddings remain
   topological; only this equator diagram uses the actual spherical lengths. */
globalThis.GraphMetric = (() => {
  const TWO_PI = 2 * Math.PI;
  const accent = "#b6403b", blue = "#315d83", secondary = "#6f777a";
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const finite = x => Number.isFinite(x);
  const decimal = (x, digits = 4) => {
    if (!finite(x)) return "unavailable";
    const threshold = .5 * 10 ** -digits;
    return (Math.abs(x) < threshold ? 0 : x).toFixed(digits);
  };
  const half = x => x === .5 ? "½" : x === 0 ? "0" : decimal(x, 2);
  const complex = angle => {
    const re = Math.cos(angle), im = Math.sin(angle);
    return `${decimal(re)} ${im < -.00005 ? "−" : "+"} ${decimal(Math.abs(im))}i`;
  };
  const vectorLength = (coefficients, variables) => coefficients.reduce((sum, coefficient, i) => sum + Number(coefficient) * variables[i], 0);
  const math = body => `<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline">${body}</math>`;
  const symbol = name => `<mi>${escape(name)}</mi>`;
  const number = value => `<mn>${escape(value).replace(/^-/u, "−")}</mn>`;
  const row = body => `<mrow>${body}</mrow>`;
  const equals = '<mo>=</mo>';
  const pi = '<mi>π</mi>';
  const frac = (top, bottom) => `<mfrac>${row(top)}${row(bottom)}</mfrac>`;
  const character = period => `<mi>χ</mi><mo stretchy="false">(</mo>${period === "1" ? number(1) : symbol("τ")}<mo stretchy="false">)</mo>`;
  const delta = variable => `<mi>Δ</mi>${symbol(variable)}`;
  function metricEquations(seed) {
    const R = symbol("R"), S = symbol("S"), one = number(1), two = number(2);
    const factor = body => `<mo stretchy="false">(</mo>${body}<mo stretchy="false">)</mo>`;
    let formulas, perimeter;
    if (seed === "double") {
      formulas = [["x", `${two}${pi}${factor(`${one}<mo>−</mo>${R}`)}`], ["y", `${two}${pi}${factor(`${R}<mo>+</mo>${two}${S}<mo>−</mo>${one}`)}`], ["L", `${pi}${factor(`${one}<mo>−</mo>${two}${S}`)}`]];
      perimeter = `${symbol("x")}<mo>+</mo>${symbol("y")}<mo>+</mo>${two}${symbol("L")}${equals}${two}${pi}`;
    } else if (seed === "tetra") {
      formulas = [["x", `${two}${pi}${factor(`${one}<mo>−</mo>${R}`)}`], ["y", `${two}${pi}${factor(`${one}<mo>−</mo>${S}`)}`], ["z", `${two}${pi}${factor(`${R}<mo>+</mo>${S}<mo>−</mo>${one}`)}`]];
      perimeter = `${symbol("x")}<mo>+</mo>${symbol("y")}<mo>+</mo>${symbol("z")}${equals}${two}${pi}`;
    } else {
      formulas = [["x", `${pi}${factor(`${one}<mo>−</mo>${two}${R}`)}`], ["y", `${pi}${factor(`${one}<mo>−</mo>${two}${S}`)}`], ["z", `${pi}${factor(`${two}${R}<mo>+</mo>${two}${S}<mo>−</mo>${one}`)}`]];
      perimeter = `${symbol("x")}<mo>+</mo>${symbol("y")}<mo>+</mo>${symbol("z")}${equals}${pi}`;
    }
    return `<details class="metric-formula-details"><summary>Length formulas</summary>${formulas.map(([name, rhs]) => `<div>${math(`${symbol(name)}${equals}${rhs}`)}</div>`).join("")}<div class="quiet">${math(perimeter)}</div></details>`;
  }
  function lengthLaw(seed, cycle) {
    const plus = '<mo>+</mo>', two = number(2), four = number(4);
    let numerator, denominator;
    if (seed === "double") {
      numerator = cycle === "a" ? `${symbol("y")}${plus}${two}${symbol("L")}` : `${symbol("x")}${plus}${symbol("y")}`;
      denominator = `${cycle === "a" ? two : four}${pi}`;
    } else {
      numerator = `${symbol(cycle === "a" ? "y" : "x")}${plus}${symbol("z")}`;
      denominator = `${two}${pi}`;
    }
    return math(`${symbol(cycle === "a" ? "R" : "S")}${equals}${frac(numerator, denominator)}`);
  }
  function expressionMath(coefficients, names) {
    const terms = [];
    coefficients.forEach((coefficient, i) => {
      const n = Number(coefficient);
      if (!n) return;
      const absolute = Math.abs(n);
      const coefficientMath = absolute === 1 ? "" : absolute === .5 ? frac(number(1), number(2)) : number(absolute);
      terms.push(`${n < 0 ? '<mo>−</mo>' : terms.length ? '<mo>+</mo>' : ""}${coefficientMath}${symbol(names[i])}`);
    });
    return terms.join("") || number(0);
  }
  function complexMath(angle) {
    return `${number(decimal(Math.cos(angle)))}<mo>${Math.sin(angle) < -.00005 ? "−" : "+"}</mo>${number(decimal(Math.abs(Math.sin(angle))))}${symbol("i")}`;
  }
  function metric(seed, R, S) {
    if (seed === "double") return {
      names:["x", "y", "L"],
      normalized:[2 * (1 - R), 2 * (R + 2 * S - 1), 1 - 2 * S],
      formulas:["x = 2π(1 − R)", "y = 2π(R + 2S − 1)", "L = π(1 − 2S)"],
      perimeter:"x + y + 2L = 2π",
      laws:{a:"R = (y + 2L)/(2π)", b:"S = (x + y)/(4π)"}
    };
    if (seed === "tetra") return {
      names:["x", "y", "z"],
      normalized:[2 * (1 - R), 2 * (1 - S), 2 * (R + S - 1)],
      formulas:["x = 2π(1 − R)", "y = 2π(1 − S)", "z = 2π(R + S − 1)"],
      perimeter:"x + y + z = 2π",
      laws:{a:"R = (y + z)/(2π)", b:"S = (x + z)/(2π)"}
    };
    return {
      names:["x", "y", "z"],
      normalized:[1 - 2 * R, 1 - 2 * S, 2 * R + 2 * S - 1],
      formulas:["x = π(1 − 2R)", "y = π(1 − 2S)", "z = π(2R + 2S − 1)"],
      perimeter:"x + y + z = π",
      laws:{a:"R = (y + z)/(2π)", b:"S = (x + z)/(2π)"}
    };
  }
  function expression(coefficients, names) {
    const terms = [];
    coefficients.forEach((coefficient, i) => {
      const n = Number(coefficient);
      if (!n) return;
      const abs = Math.abs(n), value = abs === 1 ? "" : abs === .5 ? "½" : String(abs);
      terms.push(`${n < 0 ? "−" : terms.length ? "+" : ""} ${value}${names[i]}`.trim());
    });
    return terms.join(" ") || "0";
  }
  function lengthExpression(coefficients, names, perimeter) {
    const direct = expression(coefficients, names);
    if (!Array.isArray(perimeter)) return direct;
    const difference = coefficients.map((value, i) => Number(value) - Number(perimeter[i]));
    const nonzero = difference.map((value, i) => [value, i]).filter(([value]) => value !== 0);
    if (!nonzero.length) return "2π";
    if (nonzero.length !== 1 || nonzero[0][0] !== -1) return direct;
    const complement = `2π − ${names[nonzero[0][1]]}`;
    return complement.length <= direct.length ? complement : direct;
  }
  const point = (cx, cy, radius, angle) => [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)];
  const xy = p => p.map(x => decimal(x, 3)).join(" ");
  function arc(cx, cy, radius, start, extent) {
    const from = point(cx, cy, radius, start), to = point(cx, cy, radius, start + extent);
    if (extent >= TWO_PI - 1e-9) {
      const middle = point(cx, cy, radius, start + Math.PI);
      return `M${xy(from)} A${radius} ${radius} 0 1 1 ${xy(middle)} A${radius} ${radius} 0 1 1 ${xy(from)}`;
    }
    return `M${xy(from)} A${radius} ${radius} 0 ${extent > Math.PI ? 1 : 0} 1 ${xy(to)}`;
  }
  function create({onEdge = () => {}, onHemisphere = () => {}, onCycle = () => {}} = {}) {
    const get = id => document.getElementById(id);
    let current = null, activeCycle = "a", hemisphereIndex = 0, previousSelection = "";
    const choice = get("hemisphere-choice"), boundary = get("hemisphere-boundary");
    function selectedEdgeIDs() {
      if (!current) return new Set();
      const {scene, selected, mode, chainIndex} = current;
      if (mode === "edge" && selected) return new Set([selected]);
      if (selected || !Number.isInteger(chainIndex)) return new Set();
      return new Set(scene.metadata.edges.filter(e => e.chainIndex === chainIndex).map(e => e.id));
    }
    function drawBoundary(data) {
      const {config, scene} = current, graph = config.graph;
      const hemispheres = scene.metadata.hemispheres;
      if (!hemispheres.length) return;
      hemisphereIndex = Math.max(0, Math.min(hemispheres.length - 1, hemisphereIndex));
      const hemisphere = hemispheres[hemisphereIndex], darts = graph.vertices[hemisphere.index ?? hemisphereIndex];
      const byDart = new Map(scene.metadata.edges.flatMap(e => (e.darts || []).map(d => [d, e])));
      const active = selectedEdgeIDs(), variables = data.normalized.map(x => Math.PI * x);
      const pieces = darts.map((dart, side) => {
        const edge = byDart.get(dart);
        const coefficients = graph.length_coefficients?.[dart];
        return {dart, side, edge, coefficients, length:coefficients ? vectorLength(coefficients, variables) : NaN};
      });
      const perimeter = pieces.reduce((sum, p) => sum + p.length, 0);
      const check = finite(perimeter) && Math.abs(perimeter - TWO_PI) <= 1e-9 && pieces.every(p => finite(p.length) && p.length > 0);
      if (choice) {
        const signature = hemispheres.map(h => `${h.id}:${h.degree}`).join(";");
        if (choice.dataset.metricOptions !== signature) {
          choice.innerHTML = hemispheres.map((h, i) => `<option value="${i}">${escape(h.label)} · ${h.degree} boundary side${h.degree === 1 ? "" : "s"}</option>`).join("");
          choice.dataset.metricOptions = signature;
        }
        choice.value = String(hemisphereIndex);
      }
      if (!check) {
        boundary.innerHTML = '<text x="130" y="115" text-anchor="middle" font-size="12" fill="#394444">Boundary metric is unresolved.</text>';
        get("perimeter-readout").textContent = `The boundary data do not satisfy the perimeter condition: ${finite(perimeter) ? decimal(perimeter / Math.PI, 8) + "π" : "edge data are missing"}.`;
        return;
      }
      const repeated = new Map();
      pieces.forEach(p => repeated.set(p.edge?.id, (repeated.get(p.edge?.id) || 0) + 1));
      const occurrences = new Map();
      let angle = -Math.PI / 2, markup = '<title>Spherical boundary of the selected hemisphere</title><desc>Arc lengths equal the spherical edge lengths. Each label gives an edge and its length expression. Repeated edges occur on both sides of a bridge.</desc><circle cx="130" cy="117" r="78" fill="none" stroke="#e1e5e2" stroke-width="1"/>';
      for (const piece of pieces) {
        const edgeID = piece.edge?.id || "", edgeLabel = piece.edge?.label || `Side ${piece.side + 1}`;
        const occurrence = (occurrences.get(edgeID) || 0) + 1;
        occurrences.set(edgeID, occurrence);
        const label = `${edgeLabel}${repeated.get(edgeID) > 1 ? `·${occurrence}` : ""}`;
        const mid = angle + piece.length / 2, p = point(130, 117, 100, mid), startPoint = point(130, 117, 78, angle);
        const lengthLabel = lengthExpression(piece.coefficients, data.names, graph.P);
        const labelX = p[0] < 62 ? 10 : p[0] > 198 ? 250 : p[0], anchor = p[0] < 62 ? "start" : p[0] > 198 ? "end" : "middle";
        const isSelected = active.has(edgeID), title = `${label}: ${lengthLabel} = ${decimal(piece.length / Math.PI)}π`;
        markup += `<g data-metric-edge="${escape(edgeID)}" role="button" tabindex="0" aria-label="${escape(title + ". Locate this edge in both graphs.")}"><title>${escape(title)}</title><path d="${arc(130, 117, 78, angle, piece.length)}" fill="none" stroke="transparent" stroke-width="18"/><path class="boundary-arc${isSelected ? " is-selected" : ""}" d="${arc(130, 117, 78, angle, piece.length)}" style="stroke:${isSelected ? accent : secondary}"/><text class="boundary-label" x="${decimal(labelX, 3)}" y="${decimal(p[1] - 3, 3)}" text-anchor="${anchor}" style="fill:${isSelected ? accent : secondary}"><tspan x="${decimal(labelX, 3)}" font-size="11.5" font-weight="600">${escape(label)}</tspan><tspan x="${decimal(labelX, 3)}" dy="13" font-size="10.5">${escape(lengthLabel)}</tspan></text></g><circle cx="${decimal(startPoint[0], 3)}" cy="${decimal(startPoint[1], 3)}" r="2.8" fill="white" stroke="${secondary}" stroke-width="1"/>`;
        angle += piece.length;
      }
      markup += `<text x="130" y="111" text-anchor="middle" font-size="18" fill="#394444">${escape(hemisphere.label)}</text><text x="130" y="132" text-anchor="middle" font-size="14" fill="${secondary}">2π</text>`;
      boundary.innerHTML = markup;
      const selectedPiece = pieces.find(p => active.has(p.edge?.id));
      const repeatedNote = [...repeated.values()].some(n => n > 1) ? " Repeated labels are the two sides of one edge." : "";
      get("perimeter-readout").innerHTML = `<strong>Boundary length 2π.</strong>${selectedPiece ? `<div class="selected-length" style="color:${accent}">${escape(selectedPiece.edge.label)}: <var>${escape(lengthExpression(selectedPiece.coefficients, data.names, graph.P))}</var> (${decimal(selectedPiece.length / Math.PI, 3)}π)</div>` : ""}<small>${escape(repeatedNote)}</small>`;
    }
    function drawCharacter(data) {
      const {config, source} = current, shift = config.shift || [0, 0];
      const coordinate = source[activeCycle === "a" ? 0 : 1] + shift[activeCycle === "a" ? 0 : 1];
      const angle = (activeCycle === "a" ? -1 : 1) * TWO_PI * coordinate;
      const ordinary = [148 + 72 * Math.cos(angle), 91 - 72 * Math.sin(angle)];
      const rotation = [148 + 72 * Math.cos(2 * angle), 91 - 72 * Math.sin(2 * angle)];
      document.querySelectorAll("[data-cycle]").forEach(button => button.setAttribute("aria-pressed", String(button.dataset.cycle === activeCycle)));
      const variable = activeCycle === "a" ? "R" : "S", period = activeCycle === "a" ? "1" : "τ";
      const exponent = `${activeCycle === "a" ? '<mo>−</mo>' : ""}${number(2)}${pi}${symbol("i")}<mo stretchy="false">(</mo>${symbol(variable)}<mo>+</mo>${delta(variable)}<mo stretchy="false">)</mo>`;
      get("length-character-law").innerHTML = `<div>${lengthLaw(config.seed, activeCycle)}</div><div>${math(`${character(period)}${equals}<msup><mi mathvariant="normal">e</mi>${row(exponent)}</msup>`)}</div>`;
      const rotationAngle = ((2 * angle + Math.PI) % TWO_PI + TWO_PI) % TWO_PI - Math.PI;
      const arrowEnd = [148 + 59 * Math.cos(rotationAngle), 91 - 59 * Math.sin(rotationAngle)];
      const angularArrow = Math.abs(rotationAngle) > .04 ? `<path d="M207 91 A59 59 0 0 ${rotationAngle > 0 ? 0 : 1} ${xy(arrowEnd)}" fill="none" stroke="${accent}" stroke-width="1.5" marker-end="url(#metric-rotation-arrow)"/>` : "";
      get("period-character").innerHTML = `<title>Monodromy on the period ${period}</title><desc>The blue point is the ordinary eigencharacter χ. The red hollow point is its square, the spherical rotation multiplier. The red arc indicates that rotation angle.</desc><defs><marker id="metric-rotation-arrow" markerWidth="6" markerHeight="6" refX="4.8" refY="3" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0L6 3L0 6" fill="none" stroke="${accent}" stroke-width="1.3"/></marker></defs><path d="M53 91H244 M148 8V173" fill="none" stroke="#d5ddd7" stroke-width="1"/><circle cx="148" cy="91" r="72" fill="none" stroke="${secondary}" stroke-width="1.4"/><text x="248" y="95" font-size="11" fill="${secondary}">Re</text><text x="153" y="12" font-size="11" fill="${secondary}">Im</text><text x="225" y="107" font-size="10" fill="${secondary}">1</text><text x="67" y="107" font-size="10" fill="${secondary}">−1</text>${angularArrow}<path d="M148 91L${xy(rotation)}" fill="none" stroke="${accent}" stroke-dasharray="3 3" stroke-width="1.2"/><path d="M148 91L${xy(ordinary)}" fill="none" stroke="${blue}" stroke-width="1.5"/><circle cx="${decimal(rotation[0], 3)}" cy="${decimal(rotation[1], 3)}" r="7.5" fill="white" stroke="${accent}" stroke-width="1.8"/><circle cx="${decimal(ordinary[0], 3)}" cy="${decimal(ordinary[1], 3)}" r="4.5" fill="${blue}"/><circle cx="66" cy="191" r="4" fill="${blue}"/><text x="77" y="195" font-size="11" fill="#394444">χ · ordinary character</text><circle cx="66" cy="207" r="4.5" fill="white" stroke="${accent}" stroke-width="1.3"/><text x="77" y="211" font-size="11" fill="#394444">χ² · spherical rotation</text>`;
      get("multiplier-readout").textContent = "Red arc: the spherical rotation χ².";
      const signA = Math.round(2 * shift[0]) % 2 ? "−1" : "+1", signB = Math.round(2 * shift[1]) % 2 ? "−1" : "+1";
      get("shift-readout").textContent = `Decoration shift (ΔR, ΔS) = (${half(shift[0])}, ${half(shift[1])}).`;
      get("character-readout").textContent = `Ordinary signs on (1, τ): (${signA}, ${signB}); spherical rotations are unchanged.`;
    }
    function repaint() {
      if (!current) return;
      const {config, source} = current, data = metric(config.seed, ...source);
      const formulaHost = get("metric-formulas");
      if (formulaHost.dataset.metricSeed !== config.seed) {
        formulaHost.innerHTML = metricEquations(config.seed);
        formulaHost.dataset.metricSeed = config.seed;
      }
      get("metric-values").innerHTML = data.names.map((name, i) => `<span class="length-chip"><var>${escape(name)}</var>/π = <strong>${decimal(data.normalized[i]).replace(/\.?0+$/, "")}</strong></span>`).join("");
      drawBoundary(data);drawCharacter(data);
    }
    function selectBoundaryEdge(event) {
      const node = event.target.closest?.("[data-metric-edge]");
      if (!node || !node.dataset.metricEdge) return;
      event.preventDefault();onEdge(node.dataset.metricEdge);
    }
    boundary?.addEventListener("click", selectBoundaryEdge);
    boundary?.addEventListener("keydown", event => {if (event.key === "Enter" || event.key === " ") selectBoundaryEdge(event);});
    choice?.addEventListener("change", () => {
      hemisphereIndex = Number(choice.value) || 0;
      repaint();
      const hemisphere = current?.scene.metadata.hemispheres[hemisphereIndex];
      if (hemisphere) onHemisphere(hemisphere.id);
    });
    document.querySelectorAll("[data-cycle]").forEach(button => button.addEventListener("click", () => {
      activeCycle = button.dataset.cycle === "b" ? "b" : "a";repaint();onCycle(activeCycle);
    }));
    return {update(next) {
      if (!next?.config?.graph || !next?.scene?.metadata || !Array.isArray(next.source) || !next.source.every(finite)) return;
      current = next;
      const signature = `${next.config.id}|${next.mode}|${next.selected || ""}`;
      if (signature !== previousSelection) {
        if (next.mode === "hemisphere" && next.selected) {
          const index = next.scene.metadata.hemispheres.findIndex(h => h.id === next.selected);
          if (index >= 0) hemisphereIndex = index;
        } else if (next.mode === "edge" && next.selected) {
          const edge = next.scene.metadata.edges.find(e => e.id === next.selected);
          if (edge && !edge.hemispheres.includes(hemisphereIndex)) hemisphereIndex = edge.hemispheres[0];
        }
        previousSelection = signature;
      }
      repaint();
    }};
  }
  return {create, metric};
})();

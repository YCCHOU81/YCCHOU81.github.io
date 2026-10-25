"use strict";

/* Concrete edge witnesses for developed angles.  Distinct edges remain
   distinct even when the metric assigns them equal lengths.  No perimeter
   relation is used to cancel a concrete edge contribution. */
globalThis.GraphArcEdges = (() => {
  const subscripts = "₀₁₂₃₄₅₆₇₈₉";
  const subscript = n => String(n).replace(/\d/g, d => subscripts[Number(d)]);
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const add = (target, source, scale = 1) => {
    for (const [id, coefficient] of source) {
      const n = (target.get(id) || 0) + scale * coefficient;
      if (n) target.set(id, n); else target.delete(id);
    }
    return target;
  };
  const scaled = (source, scale) => add(new Map(), source, scale);
  const textExpression = (coefficients, names) => {
    const terms = [];
    coefficients.forEach((n, i) => {
      if (!n) return;
      terms.push(`${n < 0 ? "−" : terms.length ? "+" : ""}${Math.abs(n) === 1 ? "" : Math.abs(n)}${names[i]}`);
    });
    return terms.join(" ") || "0";
  };
  const mathExpression = (coefficients, symbols) => {
    const terms = [];
    coefficients.forEach((n, i) => {
      if (!n) return;
      terms.push(`${n < 0 ? '<mo>−</mo>' : terms.length ? '<mo>+</mo>' : ''}${Math.abs(n) === 1 ? '' : `<mn>${Math.abs(n)}</mn>`}${symbols[i]}`);
    });
    return `<mrow>${terms.join("") || '<mn>0</mn>'}</mrow>`;
  };

  function resolve({config, holonomy, scene, request} = {}) {
    const graph = config?.graph, metadata = scene?.metadata;
    if (!graph?.vertices || !Array.isArray(graph.pairs) || !metadata?.edges || !holonomy?.axes || !request)
      throw new Error("A ribbon graph, its developed charts, its drawn edges, and an arc request are required.");

    const edgeById = new Map(), edgeByDart = new Map();
    for (const edge of metadata.edges) {
      if (!Array.isArray(edge.darts) || edge.darts.length !== 2) throw new Error("Each drawn edge must identify both ribbon darts.");
      edgeById.set(edge.id, edge);
      for (const dart of edge.darts) {
        if (edgeByDart.has(dart)) throw new Error("A ribbon dart belongs to more than one drawn edge.");
        edgeByDart.set(dart, edge);
      }
    }

    // α_d is the boundary prefix before d in its own hemisphere chart.
    const prefixes = new Map();
    graph.vertices.forEach(darts => {
      const prefix = new Map();
      for (const dart of darts) {
        const edge = edgeByDart.get(dart);
        if (!edge) throw new Error(`The drawing does not contain ribbon dart ${dart}.`);
        prefixes.set(dart, new Map(prefix));
        add(prefix, new Map([[edge.id, 1]]));
      }
    });

    const transitionVectors = holonomy.transitions.map(transition => {
      const first = prefixes.get(transition.d), second = prefixes.get(transition.mate);
      const edge = edgeByDart.get(transition.d);
      if (!first || !second || !edge || edgeByDart.get(transition.mate)?.id !== edge.id)
        throw new Error("A hemisphere transition does not match a drawn edge.");
      // c_d = α_d + α_d̄ + ℓ_d, including repeated edge incidences.
      return add(add(new Map(first), second), new Map([[edge.id, 1]]));
    });

    const chartOffsets = [], pending = new Set();
    function offset(vertex) {
      if (chartOffsets[vertex]) return chartOffsets[vertex];
      const chart = holonomy.charts[vertex];
      if (!chart || pending.has(vertex)) throw new Error("The developed chart tree is invalid.");
      pending.add(vertex);
      let result;
      if (chart.parent === null) result = new Map();
      else {
        const parent = holonomy.charts[chart.parent], transition = transitionVectors[chart.viaTransition];
        if (!parent || !transition) throw new Error("A developed chart transition is missing.");
        // β_child = β_parent + ε_parent c_d.
        result = add(new Map(offset(chart.parent)), transition, parent.sign);
      }
      pending.delete(vertex);
      chartOffsets[vertex] = result;
      return result;
    }

    const axisVectors = holonomy.axes.map(axis => {
      const chart = holonomy.charts[axis.vertex], prefix = prefixes.get(axis.d);
      if (!chart || !prefix) throw new Error("A developed marked point has no boundary prefix.");
      // θ_i = ε_i α_i + β_i.
      return add(new Map(offset(axis.vertex)), prefix, chart.sign);
    });

    const kind = request.kind || request.type;
    const reverse = request.reversed ? -1 : 1;
    let witness = new Map(), title, symbolText, symbolMath, anglePi = 0, axisContributions = [];
    let ordinaryPhasePi = null, ordinaryWeightCorrectionPi = 0;
    const angleSymbol = i => `<msub><mi>θ</mi><mn>${i}</mn></msub>`;
    const separationMath = i => `<mrow><mo>(</mo>${angleSymbol(i)}<mo>−</mo>${angleSymbol(3)}<mo>)</mo></mrow>`;
    const axisIndex = Number(request.i);
    if (["axis", "pair", "axis-location"].includes(kind)) {
      if (!Number.isInteger(axisIndex) || axisIndex < 0 || axisIndex >= axisVectors.length)
        throw new Error("A developed marked-point index between zero and three is required.");
      const i = axisIndex, factor = kind === "pair" ? 2 * reverse : kind === "axis" ? reverse : 1;
      const location = kind === "axis-location";
      witness = scaled(axisVectors[i], factor);
      if (!location) add(witness, axisVectors[3], -factor);
      anglePi = factor * (holonomy.axes[i].phasePi - (location ? 0 : holonomy.axes[3].phasePi));
      axisContributions = location ? [{i, coefficient: 1}] : [{i, coefficient: factor}, {i: 3, coefficient: -factor}];
      title = location ? `Locate Q${subscript(i)}` : kind === "pair"
        ? `${request.reversed ? "Reverse pair" : "Pair"} γ${subscript(i)}γ₃`
        : `Axis separation ${request.reversed ? `Q${subscript(i)} → Q₃` : `Q₃ → Q${subscript(i)}`}`;
      symbolText = location ? `θ${subscript(i)}` : `${factor === 1 ? "" : factor === -1 ? "−" : factor}(${`θ${subscript(i)} − θ₃`})`;
      symbolMath = location ? angleSymbol(i) : `${factor < 0 ? '<mo>−</mo>' : ''}${Math.abs(factor) === 1 ? '' : `<mn>${Math.abs(factor)}</mn>`}${separationMath(i)}`;
    } else if (kind === "period") {
      const key = request.cycle === "b" ? "tau" : "one", word = holonomy.selectedwords[key];
      if (!word?.terms) throw new Error("The selected marked period has no local monodromy word.");
      const period = key === "tau" ? "τ" : "1", symbolic = [];
      for (const term of word.terms) {
        const coefficient = 2 * term.power;
        add(witness, axisVectors[term.i], coefficient);
        add(witness, axisVectors[3], -coefficient);
        anglePi += coefficient * (holonomy.axes[term.i].phasePi - holonomy.axes[3].phasePi);
        ordinaryWeightCorrectionPi += term.power * (holonomy.axes[term.i].weight + holonomy.axes[3].weight);
        axisContributions.push({i: term.i, coefficient}, {i: 3, coefficient: -coefficient});
        symbolic.push({coefficient, i: term.i});
      }
      ordinaryPhasePi = word.phasePi;
      title = `Period ${period}`;
      symbolText = symbolic.map((term, j) => `${term.coefficient < 0 ? "−" : j ? "+" : ""}${Math.abs(term.coefficient)}(θ${subscript(term.i)} − θ₃)`).join(" ") || "0";
      symbolMath = mathExpression(symbolic.map(t => t.coefficient), symbolic.map(t => separationMath(t.i)));
    } else throw new Error(`Unknown developed-arc request: ${kind}.`);

    const terms = [...witness].map(([edgeId, coefficient]) => {
      const edge = edgeById.get(edgeId), dart = edge.darts[0], lengthPi = holonomy.lengthsPi[dart];
      if (!Number.isFinite(lengthPi) || !Number.isSafeInteger(coefficient)) throw new Error("A developed-edge coefficient or length is invalid.");
      return {edgeId, label: edge.label, index: edge.index, coefficient, lengthPi, lengthLabel: `ℓ${subscript(edge.index + 1)}`,
        chain: edge.chain, chainIndex: edge.chainIndex, darts: edge.darts.slice()};
    }).sort((a, b) => a.index - b.index);
    const dotProductPi = terms.reduce((sum, term) => sum + term.coefficient * term.lengthPi, 0);
    const errorPi = Math.abs(dotProductPi - anglePi);
    if (errorPi > 1e-8 * Math.max(1, Math.abs(anglePi))) throw new Error("The concrete edge sum does not reproduce the developed angle.");
    const metricCoefficients = [0, 0, 0];
    for (const term of terms) {
      const coefficients = graph.length_coefficients[term.darts[0]];
      coefficients.forEach((n, i) => {metricCoefficients[i] += term.coefficient * n;});
    }
    const metricNames = config.seed === "double" ? ["x", "y", "L"] : ["x", "y", "z"];
    const principalAnglePi = Math.atan2(Math.sin(Math.PI * anglePi), Math.cos(Math.PI * anglePi)) / Math.PI;
    const fullTurns = Math.round((anglePi - principalAnglePi) / 2);
    const note = kind === "axis-location"
      ? "These edges locate the half-turn axis. The current probe displacement also depends on its incoming point."
      : "Signed sum along the chosen chart continuation. Negative coefficients reverse direction; repeated coefficients count an edge again.";
    return {
      kind, title, note, terms, edgeIds: terms.map(term => term.edgeId), axisContributions,
      anglePi, dotProductPi, errorPi, principalAnglePi, fullTurns,
      symbolText, symbolMath,
      expressionText: textExpression(terms.map(t => t.coefficient), terms.map(t => t.lengthLabel)),
      expressionMath: mathExpression(terms.map(t => t.coefficient), terms.map(t => `<msub><mi>ℓ</mi><mn>${t.index + 1}</mn></msub>`)),
      metricNames, metricCoefficients,
      metricExpressionText: textExpression(metricCoefficients, metricNames),
      metricExpressionMath: mathExpression(metricCoefficients, metricNames.map(name => `<mi>${escape(name)}</mi>`)),
      ordinaryPhasePi, ordinaryWeightCorrectionPi,
      // This difference is a known whole-turn contribution, not an edge-length
      // simplification. The concrete edge witness itself is never shortened by
      // silently applying a hemisphere perimeter relation.
      ordinaryToProjectiveFullTurns: ordinaryWeightCorrectionPi,
      principalAngleNote: fullTurns ? `The visible principal displacement omits ${fullTurns} signed full turn${Math.abs(fullTurns) === 1 ? "" : "s"}.` : "The displayed principal displacement needs no full-turn adjustment.",
      perimeterRelationsUsed: [],
      convention: "Each ℓⱼ is the length of the concrete drawn edge Eⱼ; equal lengths do not identify distinct edges."
    };
  }
  return {resolve};
})();

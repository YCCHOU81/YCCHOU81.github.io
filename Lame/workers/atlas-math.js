/* Shared double-precision elliptic evaluation and compact mesh geometry. */
const atlasMath = (() => {
/* Elliptic monodromy atlas. No external dependencies.
 * Numerical roots, not interval certificates. Parameter continuation preserves
 * the chosen F0 chart; modular transformations are used ONLY inside evaluation.
 */

const PI = Math.PI;
const TAU = 2 * PI;
const SOURCE_MARGIN = 0.004;
const ROOT_TOL = 1e-13;
const MAX_Y = 60;
const c = (re = 0, im = 0) => ({ re, im });
const add = (a, b) => c(a.re + b.re, a.im + b.im);
const sub = (a, b) => c(a.re - b.re, a.im - b.im);
const mul = (a, b) => c(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const scale = (a, k) => c(a.re * k, a.im * k);
const abs = a => Math.hypot(a.re, a.im);
const div = (a, b) => {
  const d = b.re * b.re + b.im * b.im;
  return c((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
};
const finite = z => Number.isFinite(z.re) && Number.isFinite(z.im);
const polar = (m, a) => c(m * Math.cos(a), m * Math.sin(a));
const asComplex = z => Array.isArray(z) ? c(z[0], z[1]) : c(z.re, z.im);
const centered = x => x - Math.round(x);

/** Standard SL2 reduction for evaluation, with its covariant character action.
 * If t' = gamma t, then r'=a r-b s, s'=-c r+d s, and Z'=factor Z.
 * This function NEVER changes a solver's root or displayed marking.
 */
function evaluationChart(tau, r, s) {
  let t = asComplex(tau), rr = r, ss = s, factor = c(1, 0);
  if (!finite(t) || !(t.im > 0)) throw new RangeError('The modulus must lie in the upper half-plane.');
  for (let k = 0; k < 48; k++) {
    const integer = Math.round(t.re);
    if (integer) {
      t.re -= integer;
      rr += integer * ss;
    }
    rr = centered(rr);
    ss = centered(ss);
    if (abs(t) >= 1 - 2e-15) return { tau: t, r: rr, s: ss, factor, transformations: k };
    factor = mul(t, factor);
    t = div(c(-1, 0), t);
    [rr, ss] = [ss, -rr];
  }
  throw new RangeError('Elliptic evaluation could not find a stable chart.');
}

/** Fourier evaluation after reduction. Each exponential has nonpositive real
 * part, including near cusps, so no large theta terms are formed and cancelled.
 */
function reducedElliptic(t, r, s) {
  const ar = r + s * t.re, ai = s * t.im;
  const sign = ai >= 0 ? 1 : -1;
  const w = polar(Math.exp(-TAU * Math.abs(ai)), sign * TAU * ar);
  const denominator = sub(w, c(1, 0));
  if (abs(denominator) < 2e-15) throw new RangeError('The elliptic argument is too close to a lattice pole.');
  const quotient = div(add(w, c(1, 0)), denominator);
  const cot = c(-sign * quotient.im, sign * quotient.re);
  const csc2 = scale(div(w, mul(denominator, denominator)), -4);
  let sinSum = c(), pSum = c(), ppSum = c(), eisen2 = c(), eisen4 = c(), eisen6 = c();
  let last = Infinity, terms = 0;
  for (let k = 1; k <= 32; k++) {
    const qk = polar(Math.exp(-TAU * k * t.im), TAU * k * t.re);
    const den = sub(c(1, 0), qk);
    const coeff = div(qk, den);
    const plus = div(polar(Math.exp(-TAU * k * t.im * (1 + s)),
      TAU * k * (r + (1 + s) * t.re)), den);
    const minus = div(polar(Math.exp(-TAU * k * t.im * (1 - s)),
      TAU * k * (-r + (1 - s) * t.re)), den);
    const difference = sub(plus, minus);
    const sine = c(difference.im / 2, -difference.re / 2);
    const cosine = scale(add(plus, minus), .5);
    sinSum = add(sinSum, sine);
    pSum = add(pSum, scale(sub(coeff, cosine), k));
    ppSum = add(ppSum, scale(sine, k * k));
    eisen2 = add(eisen2, scale(coeff, k));
    eisen4 = add(eisen4, scale(coeff, k ** 3));
    eisen6 = add(eisen6, scale(coeff, k ** 5));
    last = Math.max(abs(plus), abs(minus), abs(coeff)) * k ** 5;
    terms = k;
    if (k >= 5 && last < 3e-18) break;
  }
  if (!(last < 1e-14)) throw new RangeError('The elliptic Fourier series did not converge.');
  const z = add(add(scale(cot, PI), scale(sinSum, 4 * PI)), c(0, TAU * s));
  const wp = add(scale(sub(csc2, c(1 / 3, 0)), PI * PI), scale(pSum, 8 * PI * PI));
  const wpPrime = add(scale(mul(cot, csc2), -2 * PI ** 3), scale(ppSum, 16 * PI ** 3));
  const result = { z, wp, wpPrime, terms };
  const argument = c(ar, ai);
  if (abs(argument) < .22) {
    // Remove the universal Laurent pole BEFORE forming Z^3-3 wp Z-wp'.
    // Direct polynomial evaluation loses digits when both character entries
    // approach integers, despite a misleadingly small relative residual.
    const eta1 = scale(sub(c(1, 0), scale(eisen2, 24)), PI * PI / 3);
    const g2 = scale(add(c(1, 0), scale(eisen4, 240)), 4 * PI ** 4 / 3);
    const g3 = scale(sub(c(1, 0), scale(eisen6, 504)), 8 * PI ** 6 / 27);
    const coefficients = [c(), scale(g2, 1 / 20), scale(g3, 1 / 28)];
    for (let m = 3; m <= 14; m++) {
      let sum = c();
      for (let i = 1; i <= m - 2; i++) sum = add(sum, mul(coefficients[i], coefficients[m - 1 - i]));
      coefficients.push(scale(sum, 3 / ((m - 2) * (2 * m + 3))));
    }
    let H = sub(c(0, TAU * s), mul(eta1, argument)), P = c(), Pp = c();
    const a2 = mul(argument, argument);
    let evenPower = a2, oddPower = argument;
    for (let m = 1; m <= 14; m++) {
      P = add(P, mul(coefficients[m], evenPower));
      Pp = add(Pp, scale(mul(coefficients[m], oddPower), 2 * m));
      H = sub(H, scale(mul(coefficients[m], mul(evenPower, argument)), 1 / (2 * m + 1)));
      evenPower = mul(evenPower, a2);
      oddPower = mul(oddPower, a2);
    }
    const H2 = mul(H, H);
    const stableTerms = [scale(div(H2, argument), 3), mul(H2, H),
      scale(div(P, argument), -3), scale(mul(P, H), -3), scale(Pp, -1)];
    result.stableZ2 = stableTerms.reduce(add, c());
    result.stableNormalization = Math.max(1, stableTerms.reduce((sum, term) => sum + abs(term), 0));
  }
  return result;
}

/** Return Z or Z^(2), together with elliptic values and a cancellation-aware
 * relative residual. tau can be {re,im} or [re,im].
 */
function evaluate(tau, r, s, n = 1) {
  if (n !== 1 && n !== 2) throw new RangeError('This engine implements indices one and two.');
  const chart = evaluationChart(tau, r, s);
  const e = reducedElliptic(chart.tau, chart.r, chart.s);
  const z3 = mul(mul(e.z, e.z), e.z), pz = scale(mul(e.wp, e.z), 3);
  const fReduced = n === 1 ? e.z : (e.stableZ2 ?? sub(sub(z3, pz), e.wpPrime));
  const normalization = n === 1 ? Math.max(1, abs(e.z), TAU) :
    (e.stableNormalization ?? Math.max(1, abs(z3) + abs(pz) + abs(e.wpPrime)));
  const f2 = mul(chart.factor, chart.factor), f3 = mul(f2, chart.factor);
  const value = div(fReduced, n === 1 ? chart.factor : f3);
  const result = { value, z: div(e.z, chart.factor), wp: div(e.wp, f2),
    wpPrime: div(e.wpPrime, f3), residual: abs(value),
    relativeResidual: abs(fReduced) / normalization,
    terms: e.terms, transformations: chart.transformations };
  if (!finite(value)) throw new RangeError('Elliptic evaluation lost numerical precision.');
  return result;
}

const thetaZ = (tau, r, s) => evaluate(tau, r, s, 1).value;
const Z2 = (tau, r, s) => evaluate(tau, r, s, 2).value;
const cauchyZ = thetaZ;

function inF0(tau, strict = true) {
  const t = asComplex(tau), epsilon = strict ? 2e-9 : -2e-12;
  return finite(t) && t.re > epsilon && t.re < 1 - epsilon &&
    t.im > 2e-6 && t.im < MAX_Y &&
    (t.re - .5) ** 2 + t.im ** 2 > .25 + epsilon;
}

function sourceStatus(r, s, margin = SOURCE_MARGIN) {
  if (![r, s].every(Number.isFinite)) return 'invalid_parameters';
  if (!(r > 0 && s > 0 && r < .5 && s < .5 && r + s > .5)) return 'outside_source_triangle';
  if (r < margin || s < margin || .5 - r < margin || .5 - s < margin || r + s - .5 < margin)
    return 'near_boundary';
  return 'ok';
}


function buildSourceMesh(resolution = 14, inset = .025) {
  const n = Math.max(2, Math.min(80, Math.round(resolution)));
  if (!(inset >= .009 && inset < 1 / 3)) throw new RangeError('Use a barycentric inset between .009 and 1/3.');
  const points = [], triangles = [], lookup = new Map();
  for (let i = 0; i <= n; i++) for (let j = 0; j <= n - i; j++) {
    const a = inset + (1 - 3 * inset) * i / n;
    const b = inset + (1 - 3 * inset) * j / n;
    const d = 1 - a - b;
    lookup.set(`${i},${j}`, points.length);
    points.push({ r: .5 * (a + b), s: .5 * (b + d), barycentric: [a, b, d] });
  }
  const at = (i, j) => lookup.get(`${i},${j}`);
  for (let i = 0; i < n; i++) for (let j = 0; j < n - i; j++) {
    triangles.push([at(i, j), at(i + 1, j), at(i, j + 1)]);
    if (i + j < n - 1) triangles.push([at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
  }
  return { points, triangles, resolution: n, inset };
}


return {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh};
})();

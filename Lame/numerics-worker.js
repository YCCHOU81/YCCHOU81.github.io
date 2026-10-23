/* Hemisphere Atlas: self-contained classic numerical worker.
 * Bundled without changing the numerical equations, continuation, or branch
 * labels in the checked atlas engine. Geometric anchors select the sheet.
 * Supported order: n = 1 + a + b + c <= 21. Interior samples only.
 * Interactive point requests stay in the compact source triangle with all
 * barycentric weights >= 0.04; the monodromy boundary is never sampled.
 * Results are numerical approximations, not interval-certified values.
 */
"use strict";

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


const classicalPrecision = (()=>{
/* Hemisphere Atlas: self-contained classic asymptotic worker.
 * Preserves the checked atlas numerical equations, geometric branch anchors,
 * adaptive continuation, and independent 70/100 digit comparison.
 * Results are numerical estimates, not interval-certified limits.
 */
"use strict";
/*
The MIT Licence.

Copyright (c) 2025 Michael Mclaughlin

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


*/
const Decimal = (() => {
/*!
 *  decimal.js v10.6.0
 *  An arbitrary-precision Decimal type for JavaScript.
 *  https://github.com/MikeMcl/decimal.js
 *  Copyright (c) 2025 Michael Mclaughlin <M8ch88l@gmail.com>
 *  MIT Licence
 */


// -----------------------------------  EDITABLE DEFAULTS  ------------------------------------ //


  // The maximum exponent magnitude.
  // The limit on the value of `toExpNeg`, `toExpPos`, `minE` and `maxE`.
var EXP_LIMIT = 9e15,                      // 0 to 9e15

  // The limit on the value of `precision`, and on the value of the first argument to
  // `toDecimalPlaces`, `toExponential`, `toFixed`, `toPrecision` and `toSignificantDigits`.
  MAX_DIGITS = 1e9,                        // 0 to 1e9

  // Base conversion alphabet.
  NUMERALS = '0123456789abcdef',

  // The natural logarithm of 10 (1025 digits).
  LN10 = '2.3025850929940456840179914546843642076011014886287729760333279009675726096773524802359972050895982983419677840422862486334095254650828067566662873690987816894829072083255546808437998948262331985283935053089653777326288461633662222876982198867465436674744042432743651550489343149393914796194044002221051017141748003688084012647080685567743216228355220114804663715659121373450747856947683463616792101806445070648000277502684916746550586856935673420670581136429224554405758925724208241314695689016758940256776311356919292033376587141660230105703089634572075440370847469940168269282808481184289314848524948644871927809676271275775397027668605952496716674183485704422507197965004714951050492214776567636938662976979522110718264549734772662425709429322582798502585509785265383207606726317164309505995087807523710333101197857547331541421808427543863591778117054309827482385045648019095610299291824318237525357709750539565187697510374970888692180205189339507238539205144634197265287286965110862571492198849978748873771345686209167058',

  // Pi (1025 digits).
  PI = '3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989380952572010654858632789',


  // The initial configuration properties of the Decimal constructor.
  DEFAULTS = {

    // These values must be integers within the stated ranges (inclusive).
    // Most of these values can be changed at run-time using the `Decimal.config` method.

    // The maximum number of significant digits of the result of a calculation or base conversion.
    // E.g. `Decimal.config({ precision: 20 });`
    precision: 20,                         // 1 to MAX_DIGITS

    // The rounding mode used when rounding to `precision`.
    //
    // ROUND_UP         0 Away from zero.
    // ROUND_DOWN       1 Towards zero.
    // ROUND_CEIL       2 Towards +Infinity.
    // ROUND_FLOOR      3 Towards -Infinity.
    // ROUND_HALF_UP    4 Towards nearest neighbour. If equidistant, up.
    // ROUND_HALF_DOWN  5 Towards nearest neighbour. If equidistant, down.
    // ROUND_HALF_EVEN  6 Towards nearest neighbour. If equidistant, towards even neighbour.
    // ROUND_HALF_CEIL  7 Towards nearest neighbour. If equidistant, towards +Infinity.
    // ROUND_HALF_FLOOR 8 Towards nearest neighbour. If equidistant, towards -Infinity.
    //
    // E.g.
    // `Decimal.rounding = 4;`
    // `Decimal.rounding = Decimal.ROUND_HALF_UP;`
    rounding: 4,                           // 0 to 8

    // The modulo mode used when calculating the modulus: a mod n.
    // The quotient (q = a / n) is calculated according to the corresponding rounding mode.
    // The remainder (r) is calculated as: r = a - n * q.
    //
    // UP         0 The remainder is positive if the dividend is negative, else is negative.
    // DOWN       1 The remainder has the same sign as the dividend (JavaScript %).
    // FLOOR      3 The remainder has the same sign as the divisor (Python %).
    // HALF_EVEN  6 The IEEE 754 remainder function.
    // EUCLID     9 Euclidian division. q = sign(n) * floor(a / abs(n)). Always positive.
    //
    // Truncated division (1), floored division (3), the IEEE 754 remainder (6), and Euclidian
    // division (9) are commonly used for the modulus operation. The other rounding modes can also
    // be used, but they may not give useful results.
    modulo: 1,                             // 0 to 9

    // The exponent value at and beneath which `toString` returns exponential notation.
    // JavaScript numbers: -7
    toExpNeg: -7,                          // 0 to -EXP_LIMIT

    // The exponent value at and above which `toString` returns exponential notation.
    // JavaScript numbers: 21
    toExpPos:  21,                         // 0 to EXP_LIMIT

    // The minimum exponent value, beneath which underflow to zero occurs.
    // JavaScript numbers: -324  (5e-324)
    minE: -EXP_LIMIT,                      // -1 to -EXP_LIMIT

    // The maximum exponent value, above which overflow to Infinity occurs.
    // JavaScript numbers: 308  (1.7976931348623157e+308)
    maxE: EXP_LIMIT,                       // 1 to EXP_LIMIT

    // Whether to use cryptographically-secure random number generation, if available.
    crypto: false                          // true/false
  },


// ----------------------------------- END OF EDITABLE DEFAULTS ------------------------------- //


  inexact, quadrant,
  external = true,

  decimalError = '[DecimalError] ',
  invalidArgument = decimalError + 'Invalid argument: ',
  precisionLimitExceeded = decimalError + 'Precision limit exceeded',
  cryptoUnavailable = decimalError + 'crypto unavailable',
  tag = '[object Decimal]',

  mathfloor = Math.floor,
  mathpow = Math.pow,

  isBinary = /^0b([01]+(\.[01]*)?|\.[01]+)(p[+-]?\d+)?$/i,
  isHex = /^0x([0-9a-f]+(\.[0-9a-f]*)?|\.[0-9a-f]+)(p[+-]?\d+)?$/i,
  isOctal = /^0o([0-7]+(\.[0-7]*)?|\.[0-7]+)(p[+-]?\d+)?$/i,
  isDecimal = /^(\d+(\.\d*)?|\.\d+)(e[+-]?\d+)?$/i,

  BASE = 1e7,
  LOG_BASE = 7,
  MAX_SAFE_INTEGER = 9007199254740991,

  LN10_PRECISION = LN10.length - 1,
  PI_PRECISION = PI.length - 1,

  // Decimal.prototype object
  P = { toStringTag: tag };


// Decimal prototype methods


/*
 *  absoluteValue             abs
 *  ceil
 *  clampedTo                 clamp
 *  comparedTo                cmp
 *  cosine                    cos
 *  cubeRoot                  cbrt
 *  decimalPlaces             dp
 *  dividedBy                 div
 *  dividedToIntegerBy        divToInt
 *  equals                    eq
 *  floor
 *  greaterThan               gt
 *  greaterThanOrEqualTo      gte
 *  hyperbolicCosine          cosh
 *  hyperbolicSine            sinh
 *  hyperbolicTangent         tanh
 *  inverseCosine             acos
 *  inverseHyperbolicCosine   acosh
 *  inverseHyperbolicSine     asinh
 *  inverseHyperbolicTangent  atanh
 *  inverseSine               asin
 *  inverseTangent            atan
 *  isFinite
 *  isInteger                 isInt
 *  isNaN
 *  isNegative                isNeg
 *  isPositive                isPos
 *  isZero
 *  lessThan                  lt
 *  lessThanOrEqualTo         lte
 *  logarithm                 log
 *  [maximum]                 [max]
 *  [minimum]                 [min]
 *  minus                     sub
 *  modulo                    mod
 *  naturalExponential        exp
 *  naturalLogarithm          ln
 *  negated                   neg
 *  plus                      add
 *  precision                 sd
 *  round
 *  sine                      sin
 *  squareRoot                sqrt
 *  tangent                   tan
 *  times                     mul
 *  toBinary
 *  toDecimalPlaces           toDP
 *  toExponential
 *  toFixed
 *  toFraction
 *  toHexadecimal             toHex
 *  toNearest
 *  toNumber
 *  toOctal
 *  toPower                   pow
 *  toPrecision
 *  toSignificantDigits       toSD
 *  toString
 *  truncated                 trunc
 *  valueOf                   toJSON
 */


/*
 * Return a new Decimal whose value is the absolute value of this Decimal.
 *
 */
P.absoluteValue = P.abs = function () {
  var x = new this.constructor(this);
  if (x.s < 0) x.s = 1;
  return finalise(x);
};


/*
 * Return a new Decimal whose value is the value of this Decimal rounded to a whole number in the
 * direction of positive Infinity.
 *
 */
P.ceil = function () {
  return finalise(new this.constructor(this), this.e + 1, 2);
};


/*
 * Return a new Decimal whose value is the value of this Decimal clamped to the range
 * delineated by `min` and `max`.
 *
 * min {number|string|bigint|Decimal}
 * max {number|string|bigint|Decimal}
 *
 */
P.clampedTo = P.clamp = function (min, max) {
  var k,
    x = this,
    Ctor = x.constructor;
  min = new Ctor(min);
  max = new Ctor(max);
  if (!min.s || !max.s) return new Ctor(NaN);
  if (min.gt(max)) throw Error(invalidArgument + max);
  k = x.cmp(min);
  return k < 0 ? min : x.cmp(max) > 0 ? max : new Ctor(x);
};


/*
 * Return
 *   1    if the value of this Decimal is greater than the value of `y`,
 *  -1    if the value of this Decimal is less than the value of `y`,
 *   0    if they have the same value,
 *   NaN  if the value of either Decimal is NaN.
 *
 */
P.comparedTo = P.cmp = function (y) {
  var i, j, xdL, ydL,
    x = this,
    xd = x.d,
    yd = (y = new x.constructor(y)).d,
    xs = x.s,
    ys = y.s;

  // Either NaN or ±Infinity?
  if (!xd || !yd) {
    return !xs || !ys ? NaN : xs !== ys ? xs : xd === yd ? 0 : !xd ^ xs < 0 ? 1 : -1;
  }

  // Either zero?
  if (!xd[0] || !yd[0]) return xd[0] ? xs : yd[0] ? -ys : 0;

  // Signs differ?
  if (xs !== ys) return xs;

  // Compare exponents.
  if (x.e !== y.e) return x.e > y.e ^ xs < 0 ? 1 : -1;

  xdL = xd.length;
  ydL = yd.length;

  // Compare digit by digit.
  for (i = 0, j = xdL < ydL ? xdL : ydL; i < j; ++i) {
    if (xd[i] !== yd[i]) return xd[i] > yd[i] ^ xs < 0 ? 1 : -1;
  }

  // Compare lengths.
  return xdL === ydL ? 0 : xdL > ydL ^ xs < 0 ? 1 : -1;
};


/*
 * Return a new Decimal whose value is the cosine of the value in radians of this Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-1, 1]
 *
 * cos(0)         = 1
 * cos(-0)        = 1
 * cos(Infinity)  = NaN
 * cos(-Infinity) = NaN
 * cos(NaN)       = NaN
 *
 */
P.cosine = P.cos = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (!x.d) return new Ctor(NaN);

  // cos(0) = cos(-0) = 1
  if (!x.d[0]) return new Ctor(1);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + Math.max(x.e, x.sd()) + LOG_BASE;
  Ctor.rounding = 1;

  x = cosine(Ctor, toLessThanHalfPi(Ctor, x));

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return finalise(quadrant == 2 || quadrant == 3 ? x.neg() : x, pr, rm, true);
};


/*
 *
 * Return a new Decimal whose value is the cube root of the value of this Decimal, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 *  cbrt(0)  =  0
 *  cbrt(-0) = -0
 *  cbrt(1)  =  1
 *  cbrt(-1) = -1
 *  cbrt(N)  =  N
 *  cbrt(-I) = -I
 *  cbrt(I)  =  I
 *
 * Math.cbrt(x) = (x < 0 ? -Math.pow(-x, 1/3) : Math.pow(x, 1/3))
 *
 */
P.cubeRoot = P.cbrt = function () {
  var e, m, n, r, rep, s, sd, t, t3, t3plusx,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite() || x.isZero()) return new Ctor(x);
  external = false;

  // Initial estimate.
  s = x.s * mathpow(x.s * x, 1 / 3);

   // Math.cbrt underflow/overflow?
   // Pass x to Math.pow as integer, then adjust the exponent of the result.
  if (!s || Math.abs(s) == 1 / 0) {
    n = digitsToString(x.d);
    e = x.e;

    // Adjust n exponent so it is a multiple of 3 away from x exponent.
    if (s = (e - n.length + 1) % 3) n += (s == 1 || s == -2 ? '0' : '00');
    s = mathpow(n, 1 / 3);

    // Rarely, e may be one less than the result exponent value.
    e = mathfloor((e + 1) / 3) - (e % 3 == (e < 0 ? -1 : 2));

    if (s == 1 / 0) {
      n = '5e' + e;
    } else {
      n = s.toExponential();
      n = n.slice(0, n.indexOf('e') + 1) + e;
    }

    r = new Ctor(n);
    r.s = x.s;
  } else {
    r = new Ctor(s.toString());
  }

  sd = (e = Ctor.precision) + 3;

  // Halley's method.
  // TODO? Compare Newton's method.
  for (;;) {
    t = r;
    t3 = t.times(t).times(t);
    t3plusx = t3.plus(x);
    r = divide(t3plusx.plus(x).times(t), t3plusx.plus(t3), sd + 2, 1);

    // TODO? Replace with for-loop and checkRoundingDigits.
    if (digitsToString(t.d).slice(0, sd) === (n = digitsToString(r.d)).slice(0, sd)) {
      n = n.slice(sd - 3, sd + 1);

      // The 4th rounding digit may be in error by -1 so if the 4 rounding digits are 9999 or 4999
      // , i.e. approaching a rounding boundary, continue the iteration.
      if (n == '9999' || !rep && n == '4999') {

        // On the first iteration only, check to see if rounding up gives the exact result as the
        // nines may infinitely repeat.
        if (!rep) {
          finalise(t, e + 1, 0);

          if (t.times(t).times(t).eq(x)) {
            r = t;
            break;
          }
        }

        sd += 4;
        rep = 1;
      } else {

        // If the rounding digits are null, 0{0,4} or 50{0,3}, check for an exact result.
        // If not, then there are further digits and m will be truthy.
        if (!+n || !+n.slice(1) && n.charAt(0) == '5') {

          // Truncate to the first rounding digit.
          finalise(r, e + 1, 1);
          m = !r.times(r).times(r).eq(x);
        }

        break;
      }
    }
  }

  external = true;

  return finalise(r, e, Ctor.rounding, m);
};


/*
 * Return the number of decimal places of the value of this Decimal.
 *
 */
P.decimalPlaces = P.dp = function () {
  var w,
    d = this.d,
    n = NaN;

  if (d) {
    w = d.length - 1;
    n = (w - mathfloor(this.e / LOG_BASE)) * LOG_BASE;

    // Subtract the number of trailing zeros of the last word.
    w = d[w];
    if (w) for (; w % 10 == 0; w /= 10) n--;
    if (n < 0) n = 0;
  }

  return n;
};


/*
 *  n / 0 = I
 *  n / N = N
 *  n / I = 0
 *  0 / n = 0
 *  0 / 0 = N
 *  0 / N = N
 *  0 / I = 0
 *  N / n = N
 *  N / 0 = N
 *  N / N = N
 *  N / I = N
 *  I / n = I
 *  I / 0 = I
 *  I / N = N
 *  I / I = N
 *
 * Return a new Decimal whose value is the value of this Decimal divided by `y`, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 */
P.dividedBy = P.div = function (y) {
  return divide(this, new this.constructor(y));
};


/*
 * Return a new Decimal whose value is the integer part of dividing the value of this Decimal
 * by the value of `y`, rounded to `precision` significant digits using rounding mode `rounding`.
 *
 */
P.dividedToIntegerBy = P.divToInt = function (y) {
  var x = this,
    Ctor = x.constructor;
  return finalise(divide(x, new Ctor(y), 0, 1, 1), Ctor.precision, Ctor.rounding);
};


/*
 * Return true if the value of this Decimal is equal to the value of `y`, otherwise return false.
 *
 */
P.equals = P.eq = function (y) {
  return this.cmp(y) === 0;
};


/*
 * Return a new Decimal whose value is the value of this Decimal rounded to a whole number in the
 * direction of negative Infinity.
 *
 */
P.floor = function () {
  return finalise(new this.constructor(this), this.e + 1, 3);
};


/*
 * Return true if the value of this Decimal is greater than the value of `y`, otherwise return
 * false.
 *
 */
P.greaterThan = P.gt = function (y) {
  return this.cmp(y) > 0;
};


/*
 * Return true if the value of this Decimal is greater than or equal to the value of `y`,
 * otherwise return false.
 *
 */
P.greaterThanOrEqualTo = P.gte = function (y) {
  var k = this.cmp(y);
  return k == 1 || k === 0;
};


/*
 * Return a new Decimal whose value is the hyperbolic cosine of the value in radians of this
 * Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [1, Infinity]
 *
 * cosh(x) = 1 + x^2/2! + x^4/4! + x^6/6! + ...
 *
 * cosh(0)         = 1
 * cosh(-0)        = 1
 * cosh(Infinity)  = Infinity
 * cosh(-Infinity) = Infinity
 * cosh(NaN)       = NaN
 *
 *  x        time taken (ms)   result
 * 1000      9                 9.8503555700852349694e+433
 * 10000     25                4.4034091128314607936e+4342
 * 100000    171               1.4033316802130615897e+43429
 * 1000000   3817              1.5166076984010437725e+434294
 * 10000000  abandoned after 2 minute wait
 *
 * TODO? Compare performance of cosh(x) = 0.5 * (exp(x) + exp(-x))
 *
 */
P.hyperbolicCosine = P.cosh = function () {
  var k, n, pr, rm, len,
    x = this,
    Ctor = x.constructor,
    one = new Ctor(1);

  if (!x.isFinite()) return new Ctor(x.s ? 1 / 0 : NaN);
  if (x.isZero()) return one;

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + Math.max(x.e, x.sd()) + 4;
  Ctor.rounding = 1;
  len = x.d.length;

  // Argument reduction: cos(4x) = 1 - 8cos^2(x) + 8cos^4(x) + 1
  // i.e. cos(x) = 1 - cos^2(x/4)(8 - 8cos^2(x/4))

  // Estimate the optimum number of times to use the argument reduction.
  // TODO? Estimation reused from cosine() and may not be optimal here.
  if (len < 32) {
    k = Math.ceil(len / 3);
    n = (1 / tinyPow(4, k)).toString();
  } else {
    k = 16;
    n = '2.3283064365386962890625e-10';
  }

  x = taylorSeries(Ctor, 1, x.times(n), new Ctor(1), true);

  // Reverse argument reduction
  var cosh2_x,
    i = k,
    d8 = new Ctor(8);
  for (; i--;) {
    cosh2_x = x.times(x);
    x = one.minus(cosh2_x.times(d8.minus(cosh2_x.times(d8))));
  }

  return finalise(x, Ctor.precision = pr, Ctor.rounding = rm, true);
};


/*
 * Return a new Decimal whose value is the hyperbolic sine of the value in radians of this
 * Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-Infinity, Infinity]
 *
 * sinh(x) = x + x^3/3! + x^5/5! + x^7/7! + ...
 *
 * sinh(0)         = 0
 * sinh(-0)        = -0
 * sinh(Infinity)  = Infinity
 * sinh(-Infinity) = -Infinity
 * sinh(NaN)       = NaN
 *
 * x        time taken (ms)
 * 10       2 ms
 * 100      5 ms
 * 1000     14 ms
 * 10000    82 ms
 * 100000   886 ms            1.4033316802130615897e+43429
 * 200000   2613 ms
 * 300000   5407 ms
 * 400000   8824 ms
 * 500000   13026 ms          8.7080643612718084129e+217146
 * 1000000  48543 ms
 *
 * TODO? Compare performance of sinh(x) = 0.5 * (exp(x) - exp(-x))
 *
 */
P.hyperbolicSine = P.sinh = function () {
  var k, pr, rm, len,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite() || x.isZero()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + Math.max(x.e, x.sd()) + 4;
  Ctor.rounding = 1;
  len = x.d.length;

  if (len < 3) {
    x = taylorSeries(Ctor, 2, x, x, true);
  } else {

    // Alternative argument reduction: sinh(3x) = sinh(x)(3 + 4sinh^2(x))
    // i.e. sinh(x) = sinh(x/3)(3 + 4sinh^2(x/3))
    // 3 multiplications and 1 addition

    // Argument reduction: sinh(5x) = sinh(x)(5 + sinh^2(x)(20 + 16sinh^2(x)))
    // i.e. sinh(x) = sinh(x/5)(5 + sinh^2(x/5)(20 + 16sinh^2(x/5)))
    // 4 multiplications and 2 additions

    // Estimate the optimum number of times to use the argument reduction.
    k = 1.4 * Math.sqrt(len);
    k = k > 16 ? 16 : k | 0;

    x = x.times(1 / tinyPow(5, k));
    x = taylorSeries(Ctor, 2, x, x, true);

    // Reverse argument reduction
    var sinh2_x,
      d5 = new Ctor(5),
      d16 = new Ctor(16),
      d20 = new Ctor(20);
    for (; k--;) {
      sinh2_x = x.times(x);
      x = x.times(d5.plus(sinh2_x.times(d16.times(sinh2_x).plus(d20))));
    }
  }

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return finalise(x, pr, rm, true);
};


/*
 * Return a new Decimal whose value is the hyperbolic tangent of the value in radians of this
 * Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-1, 1]
 *
 * tanh(x) = sinh(x) / cosh(x)
 *
 * tanh(0)         = 0
 * tanh(-0)        = -0
 * tanh(Infinity)  = 1
 * tanh(-Infinity) = -1
 * tanh(NaN)       = NaN
 *
 */
P.hyperbolicTangent = P.tanh = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite()) return new Ctor(x.s);
  if (x.isZero()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + 7;
  Ctor.rounding = 1;

  return divide(x.sinh(), x.cosh(), Ctor.precision = pr, Ctor.rounding = rm);
};


/*
 * Return a new Decimal whose value is the arccosine (inverse cosine) in radians of the value of
 * this Decimal.
 *
 * Domain: [-1, 1]
 * Range: [0, pi]
 *
 * acos(x) = pi/2 - asin(x)
 *
 * acos(0)       = pi/2
 * acos(-0)      = pi/2
 * acos(1)       = 0
 * acos(-1)      = pi
 * acos(1/2)     = pi/3
 * acos(-1/2)    = 2*pi/3
 * acos(|x| > 1) = NaN
 * acos(NaN)     = NaN
 *
 */
P.inverseCosine = P.acos = function () {
  var x = this,
    Ctor = x.constructor,
    k = x.abs().cmp(1),
    pr = Ctor.precision,
    rm = Ctor.rounding;

  if (k !== -1) {
    return k === 0
      // |x| is 1
      ? x.isNeg() ? getPi(Ctor, pr, rm) : new Ctor(0)
      // |x| > 1 or x is NaN
      : new Ctor(NaN);
  }

  if (x.isZero()) return getPi(Ctor, pr + 4, rm).times(0.5);

  // TODO? Special case acos(0.5) = pi/3 and acos(-0.5) = 2*pi/3

  Ctor.precision = pr + 6;
  Ctor.rounding = 1;

  // See https://github.com/MikeMcl/decimal.js/pull/217
  x = new Ctor(1).minus(x).div(x.plus(1)).sqrt().atan();

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return x.times(2);
};


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic cosine in radians of the
 * value of this Decimal.
 *
 * Domain: [1, Infinity]
 * Range: [0, Infinity]
 *
 * acosh(x) = ln(x + sqrt(x^2 - 1))
 *
 * acosh(x < 1)     = NaN
 * acosh(NaN)       = NaN
 * acosh(Infinity)  = Infinity
 * acosh(-Infinity) = NaN
 * acosh(0)         = NaN
 * acosh(-0)        = NaN
 * acosh(1)         = 0
 * acosh(-1)        = NaN
 *
 */
P.inverseHyperbolicCosine = P.acosh = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (x.lte(1)) return new Ctor(x.eq(1) ? 0 : NaN);
  if (!x.isFinite()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + Math.max(Math.abs(x.e), x.sd()) + 4;
  Ctor.rounding = 1;
  external = false;

  x = x.times(x).minus(1).sqrt().plus(x);

  external = true;
  Ctor.precision = pr;
  Ctor.rounding = rm;

  return x.ln();
};


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic sine in radians of the value
 * of this Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-Infinity, Infinity]
 *
 * asinh(x) = ln(x + sqrt(x^2 + 1))
 *
 * asinh(NaN)       = NaN
 * asinh(Infinity)  = Infinity
 * asinh(-Infinity) = -Infinity
 * asinh(0)         = 0
 * asinh(-0)        = -0
 *
 */
P.inverseHyperbolicSine = P.asinh = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite() || x.isZero()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + 2 * Math.max(Math.abs(x.e), x.sd()) + 6;
  Ctor.rounding = 1;
  external = false;

  x = x.times(x).plus(1).sqrt().plus(x);

  external = true;
  Ctor.precision = pr;
  Ctor.rounding = rm;

  return x.ln();
};


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic tangent in radians of the
 * value of this Decimal.
 *
 * Domain: [-1, 1]
 * Range: [-Infinity, Infinity]
 *
 * atanh(x) = 0.5 * ln((1 + x) / (1 - x))
 *
 * atanh(|x| > 1)   = NaN
 * atanh(NaN)       = NaN
 * atanh(Infinity)  = NaN
 * atanh(-Infinity) = NaN
 * atanh(0)         = 0
 * atanh(-0)        = -0
 * atanh(1)         = Infinity
 * atanh(-1)        = -Infinity
 *
 */
P.inverseHyperbolicTangent = P.atanh = function () {
  var pr, rm, wpr, xsd,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite()) return new Ctor(NaN);
  if (x.e >= 0) return new Ctor(x.abs().eq(1) ? x.s / 0 : x.isZero() ? x : NaN);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  xsd = x.sd();

  if (Math.max(xsd, pr) < 2 * -x.e - 1) return finalise(new Ctor(x), pr, rm, true);

  Ctor.precision = wpr = xsd - x.e;

  x = divide(x.plus(1), new Ctor(1).minus(x), wpr + pr, 1);

  Ctor.precision = pr + 4;
  Ctor.rounding = 1;

  x = x.ln();

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return x.times(0.5);
};


/*
 * Return a new Decimal whose value is the arcsine (inverse sine) in radians of the value of this
 * Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-pi/2, pi/2]
 *
 * asin(x) = 2*atan(x/(1 + sqrt(1 - x^2)))
 *
 * asin(0)       = 0
 * asin(-0)      = -0
 * asin(1/2)     = pi/6
 * asin(-1/2)    = -pi/6
 * asin(1)       = pi/2
 * asin(-1)      = -pi/2
 * asin(|x| > 1) = NaN
 * asin(NaN)     = NaN
 *
 * TODO? Compare performance of Taylor series.
 *
 */
P.inverseSine = P.asin = function () {
  var halfPi, k,
    pr, rm,
    x = this,
    Ctor = x.constructor;

  if (x.isZero()) return new Ctor(x);

  k = x.abs().cmp(1);
  pr = Ctor.precision;
  rm = Ctor.rounding;

  if (k !== -1) {

    // |x| is 1
    if (k === 0) {
      halfPi = getPi(Ctor, pr + 4, rm).times(0.5);
      halfPi.s = x.s;
      return halfPi;
    }

    // |x| > 1 or x is NaN
    return new Ctor(NaN);
  }

  // TODO? Special case asin(1/2) = pi/6 and asin(-1/2) = -pi/6

  Ctor.precision = pr + 6;
  Ctor.rounding = 1;

  x = x.div(new Ctor(1).minus(x.times(x)).sqrt().plus(1)).atan();

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return x.times(2);
};


/*
 * Return a new Decimal whose value is the arctangent (inverse tangent) in radians of the value
 * of this Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-pi/2, pi/2]
 *
 * atan(x) = x - x^3/3 + x^5/5 - x^7/7 + ...
 *
 * atan(0)         = 0
 * atan(-0)        = -0
 * atan(1)         = pi/4
 * atan(-1)        = -pi/4
 * atan(Infinity)  = pi/2
 * atan(-Infinity) = -pi/2
 * atan(NaN)       = NaN
 *
 */
P.inverseTangent = P.atan = function () {
  var i, j, k, n, px, t, r, wpr, x2,
    x = this,
    Ctor = x.constructor,
    pr = Ctor.precision,
    rm = Ctor.rounding;

  if (!x.isFinite()) {
    if (!x.s) return new Ctor(NaN);
    if (pr + 4 <= PI_PRECISION) {
      r = getPi(Ctor, pr + 4, rm).times(0.5);
      r.s = x.s;
      return r;
    }
  } else if (x.isZero()) {
    return new Ctor(x);
  } else if (x.abs().eq(1) && pr + 4 <= PI_PRECISION) {
    r = getPi(Ctor, pr + 4, rm).times(0.25);
    r.s = x.s;
    return r;
  }

  Ctor.precision = wpr = pr + 10;
  Ctor.rounding = 1;

  // TODO? if (x >= 1 && pr <= PI_PRECISION) atan(x) = halfPi * x.s - atan(1 / x);

  // Argument reduction
  // Ensure |x| < 0.42
  // atan(x) = 2 * atan(x / (1 + sqrt(1 + x^2)))

  k = Math.min(28, wpr / LOG_BASE + 2 | 0);

  for (i = k; i; --i) x = x.div(x.times(x).plus(1).sqrt().plus(1));

  external = false;

  j = Math.ceil(wpr / LOG_BASE);
  n = 1;
  x2 = x.times(x);
  r = new Ctor(x);
  px = x;

  // atan(x) = x - x^3/3 + x^5/5 - x^7/7 + ...
  for (; i !== -1;) {
    px = px.times(x2);
    t = r.minus(px.div(n += 2));

    px = px.times(x2);
    r = t.plus(px.div(n += 2));

    if (r.d[j] !== void 0) for (i = j; r.d[i] === t.d[i] && i--;);
  }

  if (k) r = r.times(2 << (k - 1));

  external = true;

  return finalise(r, Ctor.precision = pr, Ctor.rounding = rm, true);
};


/*
 * Return true if the value of this Decimal is a finite number, otherwise return false.
 *
 */
P.isFinite = function () {
  return !!this.d;
};


/*
 * Return true if the value of this Decimal is an integer, otherwise return false.
 *
 */
P.isInteger = P.isInt = function () {
  return !!this.d && mathfloor(this.e / LOG_BASE) > this.d.length - 2;
};


/*
 * Return true if the value of this Decimal is NaN, otherwise return false.
 *
 */
P.isNaN = function () {
  return !this.s;
};


/*
 * Return true if the value of this Decimal is negative, otherwise return false.
 *
 */
P.isNegative = P.isNeg = function () {
  return this.s < 0;
};


/*
 * Return true if the value of this Decimal is positive, otherwise return false.
 *
 */
P.isPositive = P.isPos = function () {
  return this.s > 0;
};


/*
 * Return true if the value of this Decimal is 0 or -0, otherwise return false.
 *
 */
P.isZero = function () {
  return !!this.d && this.d[0] === 0;
};


/*
 * Return true if the value of this Decimal is less than `y`, otherwise return false.
 *
 */
P.lessThan = P.lt = function (y) {
  return this.cmp(y) < 0;
};


/*
 * Return true if the value of this Decimal is less than or equal to `y`, otherwise return false.
 *
 */
P.lessThanOrEqualTo = P.lte = function (y) {
  return this.cmp(y) < 1;
};


/*
 * Return the logarithm of the value of this Decimal to the specified base, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * If no base is specified, return log[10](arg).
 *
 * log[base](arg) = ln(arg) / ln(base)
 *
 * The result will always be correctly rounded if the base of the log is 10, and 'almost always'
 * otherwise:
 *
 * Depending on the rounding mode, the result may be incorrectly rounded if the first fifteen
 * rounding digits are [49]99999999999999 or [50]00000000000000. In that case, the maximum error
 * between the result and the correctly rounded result will be one ulp (unit in the last place).
 *
 * log[-b](a)       = NaN
 * log[0](a)        = NaN
 * log[1](a)        = NaN
 * log[NaN](a)      = NaN
 * log[Infinity](a) = NaN
 * log[b](0)        = -Infinity
 * log[b](-0)       = -Infinity
 * log[b](-a)       = NaN
 * log[b](1)        = 0
 * log[b](Infinity) = Infinity
 * log[b](NaN)      = NaN
 *
 * [base] {number|string|bigint|Decimal} The base of the logarithm.
 *
 */
P.logarithm = P.log = function (base) {
  var isBase10, d, denominator, k, inf, num, sd, r,
    arg = this,
    Ctor = arg.constructor,
    pr = Ctor.precision,
    rm = Ctor.rounding,
    guard = 5;

  // Default base is 10.
  if (base == null) {
    base = new Ctor(10);
    isBase10 = true;
  } else {
    base = new Ctor(base);
    d = base.d;

    // Return NaN if base is negative, or non-finite, or is 0 or 1.
    if (base.s < 0 || !d || !d[0] || base.eq(1)) return new Ctor(NaN);

    isBase10 = base.eq(10);
  }

  d = arg.d;

  // Is arg negative, non-finite, 0 or 1?
  if (arg.s < 0 || !d || !d[0] || arg.eq(1)) {
    return new Ctor(d && !d[0] ? -1 / 0 : arg.s != 1 ? NaN : d ? 0 : 1 / 0);
  }

  // The result will have a non-terminating decimal expansion if base is 10 and arg is not an
  // integer power of 10.
  if (isBase10) {
    if (d.length > 1) {
      inf = true;
    } else {
      for (k = d[0]; k % 10 === 0;) k /= 10;
      inf = k !== 1;
    }
  }

  external = false;
  sd = pr + guard;
  num = naturalLogarithm(arg, sd);
  denominator = isBase10 ? getLn10(Ctor, sd + 10) : naturalLogarithm(base, sd);

  // The result will have 5 rounding digits.
  r = divide(num, denominator, sd, 1);

  // If at a rounding boundary, i.e. the result's rounding digits are [49]9999 or [50]0000,
  // calculate 10 further digits.
  //
  // If the result is known to have an infinite decimal expansion, repeat this until it is clear
  // that the result is above or below the boundary. Otherwise, if after calculating the 10
  // further digits, the last 14 are nines, round up and assume the result is exact.
  // Also assume the result is exact if the last 14 are zero.
  //
  // Example of a result that will be incorrectly rounded:
  // log[1048576](4503599627370502) = 2.60000000000000009610279511444746...
  // The above result correctly rounded using ROUND_CEIL to 1 decimal place should be 2.7, but it
  // will be given as 2.6 as there are 15 zeros immediately after the requested decimal place, so
  // the exact result would be assumed to be 2.6, which rounded using ROUND_CEIL to 1 decimal
  // place is still 2.6.
  if (checkRoundingDigits(r.d, k = pr, rm)) {

    do {
      sd += 10;
      num = naturalLogarithm(arg, sd);
      denominator = isBase10 ? getLn10(Ctor, sd + 10) : naturalLogarithm(base, sd);
      r = divide(num, denominator, sd, 1);

      if (!inf) {

        // Check for 14 nines from the 2nd rounding digit, as the first may be 4.
        if (+digitsToString(r.d).slice(k + 1, k + 15) + 1 == 1e14) {
          r = finalise(r, pr + 1, 0);
        }

        break;
      }
    } while (checkRoundingDigits(r.d, k += 10, rm));
  }

  external = true;

  return finalise(r, pr, rm);
};


/*
 * Return a new Decimal whose value is the maximum of the arguments and the value of this Decimal.
 *
 * arguments {number|string|bigint|Decimal}
 *
P.max = function () {
  Array.prototype.push.call(arguments, this);
  return maxOrMin(this.constructor, arguments, -1);
};
 */


/*
 * Return a new Decimal whose value is the minimum of the arguments and the value of this Decimal.
 *
 * arguments {number|string|bigint|Decimal}
 *
P.min = function () {
  Array.prototype.push.call(arguments, this);
  return maxOrMin(this.constructor, arguments, 1);
};
 */


/*
 *  n - 0 = n
 *  n - N = N
 *  n - I = -I
 *  0 - n = -n
 *  0 - 0 = 0
 *  0 - N = N
 *  0 - I = -I
 *  N - n = N
 *  N - 0 = N
 *  N - N = N
 *  N - I = N
 *  I - n = I
 *  I - 0 = I
 *  I - N = N
 *  I - I = N
 *
 * Return a new Decimal whose value is the value of this Decimal minus `y`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 */
P.minus = P.sub = function (y) {
  var d, e, i, j, k, len, pr, rm, xd, xe, xLTy, yd,
    x = this,
    Ctor = x.constructor;

  y = new Ctor(y);

  // If either is not finite...
  if (!x.d || !y.d) {

    // Return NaN if either is NaN.
    if (!x.s || !y.s) y = new Ctor(NaN);

    // Return y negated if x is finite and y is ±Infinity.
    else if (x.d) y.s = -y.s;

    // Return x if y is finite and x is ±Infinity.
    // Return x if both are ±Infinity with different signs.
    // Return NaN if both are ±Infinity with the same sign.
    else y = new Ctor(y.d || x.s !== y.s ? x : NaN);

    return y;
  }

  // If signs differ...
  if (x.s != y.s) {
    y.s = -y.s;
    return x.plus(y);
  }

  xd = x.d;
  yd = y.d;
  pr = Ctor.precision;
  rm = Ctor.rounding;

  // If either is zero...
  if (!xd[0] || !yd[0]) {

    // Return y negated if x is zero and y is non-zero.
    if (yd[0]) y.s = -y.s;

    // Return x if y is zero and x is non-zero.
    else if (xd[0]) y = new Ctor(x);

    // Return zero if both are zero.
    // From IEEE 754 (2008) 6.3: 0 - 0 = -0 - -0 = -0 when rounding to -Infinity.
    else return new Ctor(rm === 3 ? -0 : 0);

    return external ? finalise(y, pr, rm) : y;
  }

  // x and y are finite, non-zero numbers with the same sign.

  // Calculate base 1e7 exponents.
  e = mathfloor(y.e / LOG_BASE);
  xe = mathfloor(x.e / LOG_BASE);

  xd = xd.slice();
  k = xe - e;

  // If base 1e7 exponents differ...
  if (k) {
    xLTy = k < 0;

    if (xLTy) {
      d = xd;
      k = -k;
      len = yd.length;
    } else {
      d = yd;
      e = xe;
      len = xd.length;
    }

    // Numbers with massively different exponents would result in a very high number of
    // zeros needing to be prepended, but this can be avoided while still ensuring correct
    // rounding by limiting the number of zeros to `Math.ceil(pr / LOG_BASE) + 2`.
    i = Math.max(Math.ceil(pr / LOG_BASE), len) + 2;

    if (k > i) {
      k = i;
      d.length = 1;
    }

    // Prepend zeros to equalise exponents.
    d.reverse();
    for (i = k; i--;) d.push(0);
    d.reverse();

  // Base 1e7 exponents equal.
  } else {

    // Check digits to determine which is the bigger number.

    i = xd.length;
    len = yd.length;
    xLTy = i < len;
    if (xLTy) len = i;

    for (i = 0; i < len; i++) {
      if (xd[i] != yd[i]) {
        xLTy = xd[i] < yd[i];
        break;
      }
    }

    k = 0;
  }

  if (xLTy) {
    d = xd;
    xd = yd;
    yd = d;
    y.s = -y.s;
  }

  len = xd.length;

  // Append zeros to `xd` if shorter.
  // Don't add zeros to `yd` if shorter as subtraction only needs to start at `yd` length.
  for (i = yd.length - len; i > 0; --i) xd[len++] = 0;

  // Subtract yd from xd.
  for (i = yd.length; i > k;) {

    if (xd[--i] < yd[i]) {
      for (j = i; j && xd[--j] === 0;) xd[j] = BASE - 1;
      --xd[j];
      xd[i] += BASE;
    }

    xd[i] -= yd[i];
  }

  // Remove trailing zeros.
  for (; xd[--len] === 0;) xd.pop();

  // Remove leading zeros and adjust exponent accordingly.
  for (; xd[0] === 0; xd.shift()) --e;

  // Zero?
  if (!xd[0]) return new Ctor(rm === 3 ? -0 : 0);

  y.d = xd;
  y.e = getBase10Exponent(xd, e);

  return external ? finalise(y, pr, rm) : y;
};


/*
 *   n % 0 =  N
 *   n % N =  N
 *   n % I =  n
 *   0 % n =  0
 *  -0 % n = -0
 *   0 % 0 =  N
 *   0 % N =  N
 *   0 % I =  0
 *   N % n =  N
 *   N % 0 =  N
 *   N % N =  N
 *   N % I =  N
 *   I % n =  N
 *   I % 0 =  N
 *   I % N =  N
 *   I % I =  N
 *
 * Return a new Decimal whose value is the value of this Decimal modulo `y`, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 * The result depends on the modulo mode.
 *
 */
P.modulo = P.mod = function (y) {
  var q,
    x = this,
    Ctor = x.constructor;

  y = new Ctor(y);

  // Return NaN if x is ±Infinity or NaN, or y is NaN or ±0.
  if (!x.d || !y.s || y.d && !y.d[0]) return new Ctor(NaN);

  // Return x if y is ±Infinity or x is ±0.
  if (!y.d || x.d && !x.d[0]) {
    return finalise(new Ctor(x), Ctor.precision, Ctor.rounding);
  }

  // Prevent rounding of intermediate calculations.
  external = false;

  if (Ctor.modulo == 9) {

    // Euclidian division: q = sign(y) * floor(x / abs(y))
    // result = x - q * y    where  0 <= result < abs(y)
    q = divide(x, y.abs(), 0, 3, 1);
    q.s *= y.s;
  } else {
    q = divide(x, y, 0, Ctor.modulo, 1);
  }

  q = q.times(y);

  external = true;

  return x.minus(q);
};


/*
 * Return a new Decimal whose value is the natural exponential of the value of this Decimal,
 * i.e. the base e raised to the power the value of this Decimal, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 */
P.naturalExponential = P.exp = function () {
  return naturalExponential(this);
};


/*
 * Return a new Decimal whose value is the natural logarithm of the value of this Decimal,
 * rounded to `precision` significant digits using rounding mode `rounding`.
 *
 */
P.naturalLogarithm = P.ln = function () {
  return naturalLogarithm(this);
};


/*
 * Return a new Decimal whose value is the value of this Decimal negated, i.e. as if multiplied by
 * -1.
 *
 */
P.negated = P.neg = function () {
  var x = new this.constructor(this);
  x.s = -x.s;
  return finalise(x);
};


/*
 *  n + 0 = n
 *  n + N = N
 *  n + I = I
 *  0 + n = n
 *  0 + 0 = 0
 *  0 + N = N
 *  0 + I = I
 *  N + n = N
 *  N + 0 = N
 *  N + N = N
 *  N + I = N
 *  I + n = I
 *  I + 0 = I
 *  I + N = N
 *  I + I = I
 *
 * Return a new Decimal whose value is the value of this Decimal plus `y`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 */
P.plus = P.add = function (y) {
  var carry, d, e, i, k, len, pr, rm, xd, yd,
    x = this,
    Ctor = x.constructor;

  y = new Ctor(y);

  // If either is not finite...
  if (!x.d || !y.d) {

    // Return NaN if either is NaN.
    if (!x.s || !y.s) y = new Ctor(NaN);

    // Return x if y is finite and x is ±Infinity.
    // Return x if both are ±Infinity with the same sign.
    // Return NaN if both are ±Infinity with different signs.
    // Return y if x is finite and y is ±Infinity.
    else if (!x.d) y = new Ctor(y.d || x.s === y.s ? x : NaN);

    return y;
  }

   // If signs differ...
  if (x.s != y.s) {
    y.s = -y.s;
    return x.minus(y);
  }

  xd = x.d;
  yd = y.d;
  pr = Ctor.precision;
  rm = Ctor.rounding;

  // If either is zero...
  if (!xd[0] || !yd[0]) {

    // Return x if y is zero.
    // Return y if y is non-zero.
    if (!yd[0]) y = new Ctor(x);

    return external ? finalise(y, pr, rm) : y;
  }

  // x and y are finite, non-zero numbers with the same sign.

  // Calculate base 1e7 exponents.
  k = mathfloor(x.e / LOG_BASE);
  e = mathfloor(y.e / LOG_BASE);

  xd = xd.slice();
  i = k - e;

  // If base 1e7 exponents differ...
  if (i) {

    if (i < 0) {
      d = xd;
      i = -i;
      len = yd.length;
    } else {
      d = yd;
      e = k;
      len = xd.length;
    }

    // Limit number of zeros prepended to max(ceil(pr / LOG_BASE), len) + 1.
    k = Math.ceil(pr / LOG_BASE);
    len = k > len ? k + 1 : len + 1;

    if (i > len) {
      i = len;
      d.length = 1;
    }

    // Prepend zeros to equalise exponents. Note: Faster to use reverse then do unshifts.
    d.reverse();
    for (; i--;) d.push(0);
    d.reverse();
  }

  len = xd.length;
  i = yd.length;

  // If yd is longer than xd, swap xd and yd so xd points to the longer array.
  if (len - i < 0) {
    i = len;
    d = yd;
    yd = xd;
    xd = d;
  }

  // Only start adding at yd.length - 1 as the further digits of xd can be left as they are.
  for (carry = 0; i;) {
    carry = (xd[--i] = xd[i] + yd[i] + carry) / BASE | 0;
    xd[i] %= BASE;
  }

  if (carry) {
    xd.unshift(carry);
    ++e;
  }

  // Remove trailing zeros.
  // No need to check for zero, as +x + +y != 0 && -x + -y != 0
  for (len = xd.length; xd[--len] == 0;) xd.pop();

  y.d = xd;
  y.e = getBase10Exponent(xd, e);

  return external ? finalise(y, pr, rm) : y;
};


/*
 * Return the number of significant digits of the value of this Decimal.
 *
 * [z] {boolean|number} Whether to count integer-part trailing zeros: true, false, 1 or 0.
 *
 */
P.precision = P.sd = function (z) {
  var k,
    x = this;

  if (z !== void 0 && z !== !!z && z !== 1 && z !== 0) throw Error(invalidArgument + z);

  if (x.d) {
    k = getPrecision(x.d);
    if (z && x.e + 1 > k) k = x.e + 1;
  } else {
    k = NaN;
  }

  return k;
};


/*
 * Return a new Decimal whose value is the value of this Decimal rounded to a whole number using
 * rounding mode `rounding`.
 *
 */
P.round = function () {
  var x = this,
    Ctor = x.constructor;

  return finalise(new Ctor(x), x.e + 1, Ctor.rounding);
};


/*
 * Return a new Decimal whose value is the sine of the value in radians of this Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-1, 1]
 *
 * sin(x) = x - x^3/3! + x^5/5! - ...
 *
 * sin(0)         = 0
 * sin(-0)        = -0
 * sin(Infinity)  = NaN
 * sin(-Infinity) = NaN
 * sin(NaN)       = NaN
 *
 */
P.sine = P.sin = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite()) return new Ctor(NaN);
  if (x.isZero()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + Math.max(x.e, x.sd()) + LOG_BASE;
  Ctor.rounding = 1;

  x = sine(Ctor, toLessThanHalfPi(Ctor, x));

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return finalise(quadrant > 2 ? x.neg() : x, pr, rm, true);
};


/*
 * Return a new Decimal whose value is the square root of this Decimal, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 *  sqrt(-n) =  N
 *  sqrt(N)  =  N
 *  sqrt(-I) =  N
 *  sqrt(I)  =  I
 *  sqrt(0)  =  0
 *  sqrt(-0) = -0
 *
 */
P.squareRoot = P.sqrt = function () {
  var m, n, sd, r, rep, t,
    x = this,
    d = x.d,
    e = x.e,
    s = x.s,
    Ctor = x.constructor;

  // Negative/NaN/Infinity/zero?
  if (s !== 1 || !d || !d[0]) {
    return new Ctor(!s || s < 0 && (!d || d[0]) ? NaN : d ? x : 1 / 0);
  }

  external = false;

  // Initial estimate.
  s = Math.sqrt(+x);

  // Math.sqrt underflow/overflow?
  // Pass x to Math.sqrt as integer, then adjust the exponent of the result.
  if (s == 0 || s == 1 / 0) {
    n = digitsToString(d);

    if ((n.length + e) % 2 == 0) n += '0';
    s = Math.sqrt(n);
    e = mathfloor((e + 1) / 2) - (e < 0 || e % 2);

    if (s == 1 / 0) {
      n = '5e' + e;
    } else {
      n = s.toExponential();
      n = n.slice(0, n.indexOf('e') + 1) + e;
    }

    r = new Ctor(n);
  } else {
    r = new Ctor(s.toString());
  }

  sd = (e = Ctor.precision) + 3;

  // Newton-Raphson iteration.
  for (;;) {
    t = r;
    r = t.plus(divide(x, t, sd + 2, 1)).times(0.5);

    // TODO? Replace with for-loop and checkRoundingDigits.
    if (digitsToString(t.d).slice(0, sd) === (n = digitsToString(r.d)).slice(0, sd)) {
      n = n.slice(sd - 3, sd + 1);

      // The 4th rounding digit may be in error by -1 so if the 4 rounding digits are 9999 or
      // 4999, i.e. approaching a rounding boundary, continue the iteration.
      if (n == '9999' || !rep && n == '4999') {

        // On the first iteration only, check to see if rounding up gives the exact result as the
        // nines may infinitely repeat.
        if (!rep) {
          finalise(t, e + 1, 0);

          if (t.times(t).eq(x)) {
            r = t;
            break;
          }
        }

        sd += 4;
        rep = 1;
      } else {

        // If the rounding digits are null, 0{0,4} or 50{0,3}, check for an exact result.
        // If not, then there are further digits and m will be truthy.
        if (!+n || !+n.slice(1) && n.charAt(0) == '5') {

          // Truncate to the first rounding digit.
          finalise(r, e + 1, 1);
          m = !r.times(r).eq(x);
        }

        break;
      }
    }
  }

  external = true;

  return finalise(r, e, Ctor.rounding, m);
};


/*
 * Return a new Decimal whose value is the tangent of the value in radians of this Decimal.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-Infinity, Infinity]
 *
 * tan(0)         = 0
 * tan(-0)        = -0
 * tan(Infinity)  = NaN
 * tan(-Infinity) = NaN
 * tan(NaN)       = NaN
 *
 */
P.tangent = P.tan = function () {
  var pr, rm,
    x = this,
    Ctor = x.constructor;

  if (!x.isFinite()) return new Ctor(NaN);
  if (x.isZero()) return new Ctor(x);

  pr = Ctor.precision;
  rm = Ctor.rounding;
  Ctor.precision = pr + 10;
  Ctor.rounding = 1;

  x = x.sin();
  x.s = 1;
  x = divide(x, new Ctor(1).minus(x.times(x)).sqrt(), pr + 10, 0);

  Ctor.precision = pr;
  Ctor.rounding = rm;

  return finalise(quadrant == 2 || quadrant == 4 ? x.neg() : x, pr, rm, true);
};


/*
 *  n * 0 = 0
 *  n * N = N
 *  n * I = I
 *  0 * n = 0
 *  0 * 0 = 0
 *  0 * N = N
 *  0 * I = N
 *  N * n = N
 *  N * 0 = N
 *  N * N = N
 *  N * I = N
 *  I * n = I
 *  I * 0 = N
 *  I * N = N
 *  I * I = I
 *
 * Return a new Decimal whose value is this Decimal times `y`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 */
P.times = P.mul = function (y) {
  var carry, e, i, k, r, rL, t, xdL, ydL,
    x = this,
    Ctor = x.constructor,
    xd = x.d,
    yd = (y = new Ctor(y)).d;

  y.s *= x.s;

   // If either is NaN, ±Infinity or ±0...
  if (!xd || !xd[0] || !yd || !yd[0]) {

    return new Ctor(!y.s || xd && !xd[0] && !yd || yd && !yd[0] && !xd

      // Return NaN if either is NaN.
      // Return NaN if x is ±0 and y is ±Infinity, or y is ±0 and x is ±Infinity.
      ? NaN

      // Return ±Infinity if either is ±Infinity.
      // Return ±0 if either is ±0.
      : !xd || !yd ? y.s / 0 : y.s * 0);
  }

  e = mathfloor(x.e / LOG_BASE) + mathfloor(y.e / LOG_BASE);
  xdL = xd.length;
  ydL = yd.length;

  // Ensure xd points to the longer array.
  if (xdL < ydL) {
    r = xd;
    xd = yd;
    yd = r;
    rL = xdL;
    xdL = ydL;
    ydL = rL;
  }

  // Initialise the result array with zeros.
  r = [];
  rL = xdL + ydL;
  for (i = rL; i--;) r.push(0);

  // Multiply!
  for (i = ydL; --i >= 0;) {
    carry = 0;
    for (k = xdL + i; k > i;) {
      t = r[k] + yd[i] * xd[k - i - 1] + carry;
      r[k--] = t % BASE | 0;
      carry = t / BASE | 0;
    }

    r[k] = (r[k] + carry) % BASE | 0;
  }

  // Remove trailing zeros.
  for (; !r[--rL];) r.pop();

  if (carry) ++e;
  else r.shift();

  y.d = r;
  y.e = getBase10Exponent(r, e);

  return external ? finalise(y, Ctor.precision, Ctor.rounding) : y;
};


/*
 * Return a string representing the value of this Decimal in base 2, round to `sd` significant
 * digits using rounding mode `rm`.
 *
 * If the optional `sd` argument is present then return binary exponential notation.
 *
 * [sd] {number} Significant digits. Integer, 1 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toBinary = function (sd, rm) {
  return toStringBinary(this, 2, sd, rm);
};


/*
 * Return a new Decimal whose value is the value of this Decimal rounded to a maximum of `dp`
 * decimal places using rounding mode `rm` or `rounding` if `rm` is omitted.
 *
 * If `dp` is omitted, return a new Decimal whose value is the value of this Decimal.
 *
 * [dp] {number} Decimal places. Integer, 0 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toDecimalPlaces = P.toDP = function (dp, rm) {
  var x = this,
    Ctor = x.constructor;

  x = new Ctor(x);
  if (dp === void 0) return x;

  checkInt32(dp, 0, MAX_DIGITS);

  if (rm === void 0) rm = Ctor.rounding;
  else checkInt32(rm, 0, 8);

  return finalise(x, dp + x.e + 1, rm);
};


/*
 * Return a string representing the value of this Decimal in exponential notation rounded to
 * `dp` fixed decimal places using rounding mode `rounding`.
 *
 * [dp] {number} Decimal places. Integer, 0 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toExponential = function (dp, rm) {
  var str,
    x = this,
    Ctor = x.constructor;

  if (dp === void 0) {
    str = finiteToString(x, true);
  } else {
    checkInt32(dp, 0, MAX_DIGITS);

    if (rm === void 0) rm = Ctor.rounding;
    else checkInt32(rm, 0, 8);

    x = finalise(new Ctor(x), dp + 1, rm);
    str = finiteToString(x, true, dp + 1);
  }

  return x.isNeg() && !x.isZero() ? '-' + str : str;
};


/*
 * Return a string representing the value of this Decimal in normal (fixed-point) notation to
 * `dp` fixed decimal places and rounded using rounding mode `rm` or `rounding` if `rm` is
 * omitted.
 *
 * As with JavaScript numbers, (-0).toFixed(0) is '0', but e.g. (-0.00001).toFixed(0) is '-0'.
 *
 * [dp] {number} Decimal places. Integer, 0 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 * (-0).toFixed(0) is '0', but (-0.1).toFixed(0) is '-0'.
 * (-0).toFixed(1) is '0.0', but (-0.01).toFixed(1) is '-0.0'.
 * (-0).toFixed(3) is '0.000'.
 * (-0.5).toFixed(0) is '-0'.
 *
 */
P.toFixed = function (dp, rm) {
  var str, y,
    x = this,
    Ctor = x.constructor;

  if (dp === void 0) {
    str = finiteToString(x);
  } else {
    checkInt32(dp, 0, MAX_DIGITS);

    if (rm === void 0) rm = Ctor.rounding;
    else checkInt32(rm, 0, 8);

    y = finalise(new Ctor(x), dp + x.e + 1, rm);
    str = finiteToString(y, false, dp + y.e + 1);
  }

  // To determine whether to add the minus sign look at the value before it was rounded,
  // i.e. look at `x` rather than `y`.
  return x.isNeg() && !x.isZero() ? '-' + str : str;
};


/*
 * Return an array representing the value of this Decimal as a simple fraction with an integer
 * numerator and an integer denominator.
 *
 * The denominator will be a positive non-zero value less than or equal to the specified maximum
 * denominator. If a maximum denominator is not specified, the denominator will be the lowest
 * value necessary to represent the number exactly.
 *
 * [maxD] {number|string|bigint|Decimal} Maximum denominator. Integer >= 1 and < Infinity.
 *
 */
P.toFraction = function (maxD) {
  var d, d0, d1, d2, e, k, n, n0, n1, pr, q, r,
    x = this,
    xd = x.d,
    Ctor = x.constructor;

  if (!xd) return new Ctor(x);

  n1 = d0 = new Ctor(1);
  d1 = n0 = new Ctor(0);

  d = new Ctor(d1);
  e = d.e = getPrecision(xd) - x.e - 1;
  k = e % LOG_BASE;
  d.d[0] = mathpow(10, k < 0 ? LOG_BASE + k : k);

  if (maxD == null) {

    // d is 10**e, the minimum max-denominator needed.
    maxD = e > 0 ? d : n1;
  } else {
    n = new Ctor(maxD);
    if (!n.isInt() || n.lt(n1)) throw Error(invalidArgument + n);
    maxD = n.gt(d) ? (e > 0 ? d : n1) : n;
  }

  external = false;
  n = new Ctor(digitsToString(xd));
  pr = Ctor.precision;
  Ctor.precision = e = xd.length * LOG_BASE * 2;

  for (;;)  {
    q = divide(n, d, 0, 1, 1);
    d2 = d0.plus(q.times(d1));
    if (d2.cmp(maxD) == 1) break;
    d0 = d1;
    d1 = d2;
    d2 = n1;
    n1 = n0.plus(q.times(d2));
    n0 = d2;
    d2 = d;
    d = n.minus(q.times(d2));
    n = d2;
  }

  d2 = divide(maxD.minus(d0), d1, 0, 1, 1);
  n0 = n0.plus(d2.times(n1));
  d0 = d0.plus(d2.times(d1));
  n0.s = n1.s = x.s;

  // Determine which fraction is closer to x, n0/d0 or n1/d1?
  r = divide(n1, d1, e, 1).minus(x).abs().cmp(divide(n0, d0, e, 1).minus(x).abs()) < 1
      ? [n1, d1] : [n0, d0];

  Ctor.precision = pr;
  external = true;

  return r;
};


/*
 * Return a string representing the value of this Decimal in base 16, round to `sd` significant
 * digits using rounding mode `rm`.
 *
 * If the optional `sd` argument is present then return binary exponential notation.
 *
 * [sd] {number} Significant digits. Integer, 1 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toHexadecimal = P.toHex = function (sd, rm) {
  return toStringBinary(this, 16, sd, rm);
};


/*
 * Returns a new Decimal whose value is the nearest multiple of `y` in the direction of rounding
 * mode `rm`, or `Decimal.rounding` if `rm` is omitted, to the value of this Decimal.
 *
 * The return value will always have the same sign as this Decimal, unless either this Decimal
 * or `y` is NaN, in which case the return value will be also be NaN.
 *
 * The return value is not affected by the value of `precision`.
 *
 * y {number|string|bigint|Decimal} The magnitude to round to a multiple of.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 * 'toNearest() rounding mode not an integer: {rm}'
 * 'toNearest() rounding mode out of range: {rm}'
 *
 */
P.toNearest = function (y, rm) {
  var x = this,
    Ctor = x.constructor;

  x = new Ctor(x);

  if (y == null) {

    // If x is not finite, return x.
    if (!x.d) return x;

    y = new Ctor(1);
    rm = Ctor.rounding;
  } else {
    y = new Ctor(y);
    if (rm === void 0) {
      rm = Ctor.rounding;
    } else {
      checkInt32(rm, 0, 8);
    }

    // If x is not finite, return x if y is not NaN, else NaN.
    if (!x.d) return y.s ? x : y;

    // If y is not finite, return Infinity with the sign of x if y is Infinity, else NaN.
    if (!y.d) {
      if (y.s) y.s = x.s;
      return y;
    }
  }

  // If y is not zero, calculate the nearest multiple of y to x.
  if (y.d[0]) {
    external = false;
    x = divide(x, y, 0, rm, 1).times(y);
    external = true;
    finalise(x);

  // If y is zero, return zero with the sign of x.
  } else {
    y.s = x.s;
    x = y;
  }

  return x;
};


/*
 * Return the value of this Decimal converted to a number primitive.
 * Zero keeps its sign.
 *
 */
P.toNumber = function () {
  return +this;
};


/*
 * Return a string representing the value of this Decimal in base 8, round to `sd` significant
 * digits using rounding mode `rm`.
 *
 * If the optional `sd` argument is present then return binary exponential notation.
 *
 * [sd] {number} Significant digits. Integer, 1 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toOctal = function (sd, rm) {
  return toStringBinary(this, 8, sd, rm);
};


/*
 * Return a new Decimal whose value is the value of this Decimal raised to the power `y`, rounded
 * to `precision` significant digits using rounding mode `rounding`.
 *
 * ECMAScript compliant.
 *
 *   pow(x, NaN)                           = NaN
 *   pow(x, ±0)                            = 1

 *   pow(NaN, non-zero)                    = NaN
 *   pow(abs(x) > 1, +Infinity)            = +Infinity
 *   pow(abs(x) > 1, -Infinity)            = +0
 *   pow(abs(x) == 1, ±Infinity)           = NaN
 *   pow(abs(x) < 1, +Infinity)            = +0
 *   pow(abs(x) < 1, -Infinity)            = +Infinity
 *   pow(+Infinity, y > 0)                 = +Infinity
 *   pow(+Infinity, y < 0)                 = +0
 *   pow(-Infinity, odd integer > 0)       = -Infinity
 *   pow(-Infinity, even integer > 0)      = +Infinity
 *   pow(-Infinity, odd integer < 0)       = -0
 *   pow(-Infinity, even integer < 0)      = +0
 *   pow(+0, y > 0)                        = +0
 *   pow(+0, y < 0)                        = +Infinity
 *   pow(-0, odd integer > 0)              = -0
 *   pow(-0, even integer > 0)             = +0
 *   pow(-0, odd integer < 0)              = -Infinity
 *   pow(-0, even integer < 0)             = +Infinity
 *   pow(finite x < 0, finite non-integer) = NaN
 *
 * For non-integer or very large exponents pow(x, y) is calculated using
 *
 *   x^y = exp(y*ln(x))
 *
 * Assuming the first 15 rounding digits are each equally likely to be any digit 0-9, the
 * probability of an incorrectly rounded result
 * P([49]9{14} | [50]0{14}) = 2 * 0.2 * 10^-14 = 4e-15 = 1/2.5e+14
 * i.e. 1 in 250,000,000,000,000
 *
 * If a result is incorrectly rounded the maximum error will be 1 ulp (unit in last place).
 *
 * y {number|string|bigint|Decimal} The power to which to raise this Decimal.
 *
 */
P.toPower = P.pow = function (y) {
  var e, k, pr, r, rm, s,
    x = this,
    Ctor = x.constructor,
    yn = +(y = new Ctor(y));

  // Either ±Infinity, NaN or ±0?
  if (!x.d || !y.d || !x.d[0] || !y.d[0]) return new Ctor(mathpow(+x, yn));

  x = new Ctor(x);

  if (x.eq(1)) return x;

  pr = Ctor.precision;
  rm = Ctor.rounding;

  if (y.eq(1)) return finalise(x, pr, rm);

  // y exponent
  e = mathfloor(y.e / LOG_BASE);

  // If y is a small integer use the 'exponentiation by squaring' algorithm.
  if (e >= y.d.length - 1 && (k = yn < 0 ? -yn : yn) <= MAX_SAFE_INTEGER) {
    r = intPow(Ctor, x, k, pr);
    return y.s < 0 ? new Ctor(1).div(r) : finalise(r, pr, rm);
  }

  s = x.s;

  // if x is negative
  if (s < 0) {

    // if y is not an integer
    if (e < y.d.length - 1) return new Ctor(NaN);

    // Result is positive if x is negative and the last digit of integer y is even.
    if ((y.d[e] & 1) == 0) s = 1;

    // if x.eq(-1)
    if (x.e == 0 && x.d[0] == 1 && x.d.length == 1) {
      x.s = s;
      return x;
    }
  }

  // Estimate result exponent.
  // x^y = 10^e,  where e = y * log10(x)
  // log10(x) = log10(x_significand) + x_exponent
  // log10(x_significand) = ln(x_significand) / ln(10)
  k = mathpow(+x, yn);
  e = k == 0 || !isFinite(k)
    ? mathfloor(yn * (Math.log('0.' + digitsToString(x.d)) / Math.LN10 + x.e + 1))
    : new Ctor(k + '').e;

  // Exponent estimate may be incorrect e.g. x: 0.999999999999999999, y: 2.29, e: 0, r.e: -1.

  // Overflow/underflow?
  if (e > Ctor.maxE + 1 || e < Ctor.minE - 1) return new Ctor(e > 0 ? s / 0 : 0);

  external = false;
  Ctor.rounding = x.s = 1;

  // Estimate the extra guard digits needed to ensure five correct rounding digits from
  // naturalLogarithm(x). Example of failure without these extra digits (precision: 10):
  // new Decimal(2.32456).pow('2087987436534566.46411')
  // should be 1.162377823e+764914905173815, but is 1.162355823e+764914905173815
  k = Math.min(12, (e + '').length);

  // r = x^y = exp(y*ln(x))
  r = naturalExponential(y.times(naturalLogarithm(x, pr + k)), pr);

  // r may be Infinity, e.g. (0.9999999999999999).pow(-1e+40)
  if (r.d) {

    // Truncate to the required precision plus five rounding digits.
    r = finalise(r, pr + 5, 1);

    // If the rounding digits are [49]9999 or [50]0000 increase the precision by 10 and recalculate
    // the result.
    if (checkRoundingDigits(r.d, pr, rm)) {
      e = pr + 10;

      // Truncate to the increased precision plus five rounding digits.
      r = finalise(naturalExponential(y.times(naturalLogarithm(x, e + k)), e), e + 5, 1);

      // Check for 14 nines from the 2nd rounding digit (the first rounding digit may be 4 or 9).
      if (+digitsToString(r.d).slice(pr + 1, pr + 15) + 1 == 1e14) {
        r = finalise(r, pr + 1, 0);
      }
    }
  }

  r.s = s;
  external = true;
  Ctor.rounding = rm;

  return finalise(r, pr, rm);
};


/*
 * Return a string representing the value of this Decimal rounded to `sd` significant digits
 * using rounding mode `rounding`.
 *
 * Return exponential notation if `sd` is less than the number of digits necessary to represent
 * the integer part of the value in normal notation.
 *
 * [sd] {number} Significant digits. Integer, 1 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 */
P.toPrecision = function (sd, rm) {
  var str,
    x = this,
    Ctor = x.constructor;

  if (sd === void 0) {
    str = finiteToString(x, x.e <= Ctor.toExpNeg || x.e >= Ctor.toExpPos);
  } else {
    checkInt32(sd, 1, MAX_DIGITS);

    if (rm === void 0) rm = Ctor.rounding;
    else checkInt32(rm, 0, 8);

    x = finalise(new Ctor(x), sd, rm);
    str = finiteToString(x, sd <= x.e || x.e <= Ctor.toExpNeg, sd);
  }

  return x.isNeg() && !x.isZero() ? '-' + str : str;
};


/*
 * Return a new Decimal whose value is the value of this Decimal rounded to a maximum of `sd`
 * significant digits using rounding mode `rm`, or to `precision` and `rounding` respectively if
 * omitted.
 *
 * [sd] {number} Significant digits. Integer, 1 to MAX_DIGITS inclusive.
 * [rm] {number} Rounding mode. Integer, 0 to 8 inclusive.
 *
 * 'toSD() digits out of range: {sd}'
 * 'toSD() digits not an integer: {sd}'
 * 'toSD() rounding mode not an integer: {rm}'
 * 'toSD() rounding mode out of range: {rm}'
 *
 */
P.toSignificantDigits = P.toSD = function (sd, rm) {
  var x = this,
    Ctor = x.constructor;

  if (sd === void 0) {
    sd = Ctor.precision;
    rm = Ctor.rounding;
  } else {
    checkInt32(sd, 1, MAX_DIGITS);

    if (rm === void 0) rm = Ctor.rounding;
    else checkInt32(rm, 0, 8);
  }

  return finalise(new Ctor(x), sd, rm);
};


/*
 * Return a string representing the value of this Decimal.
 *
 * Return exponential notation if this Decimal has a positive exponent equal to or greater than
 * `toExpPos`, or a negative exponent equal to or less than `toExpNeg`.
 *
 */
P.toString = function () {
  var x = this,
    Ctor = x.constructor,
    str = finiteToString(x, x.e <= Ctor.toExpNeg || x.e >= Ctor.toExpPos);

  return x.isNeg() && !x.isZero() ? '-' + str : str;
};


/*
 * Return a new Decimal whose value is the value of this Decimal truncated to a whole number.
 *
 */
P.truncated = P.trunc = function () {
  return finalise(new this.constructor(this), this.e + 1, 1);
};


/*
 * Return a string representing the value of this Decimal.
 * Unlike `toString`, negative zero will include the minus sign.
 *
 */
P.valueOf = P.toJSON = function () {
  var x = this,
    Ctor = x.constructor,
    str = finiteToString(x, x.e <= Ctor.toExpNeg || x.e >= Ctor.toExpPos);

  return x.isNeg() ? '-' + str : str;
};


// Helper functions for Decimal.prototype (P) and/or Decimal methods, and their callers.


/*
 *  digitsToString           P.cubeRoot, P.logarithm, P.squareRoot, P.toFraction, P.toPower,
 *                           finiteToString, naturalExponential, naturalLogarithm
 *  checkInt32               P.toDecimalPlaces, P.toExponential, P.toFixed, P.toNearest,
 *                           P.toPrecision, P.toSignificantDigits, toStringBinary, random
 *  checkRoundingDigits      P.logarithm, P.toPower, naturalExponential, naturalLogarithm
 *  convertBase              toStringBinary, parseOther
 *  cos                      P.cos
 *  divide                   P.atanh, P.cubeRoot, P.dividedBy, P.dividedToIntegerBy,
 *                           P.logarithm, P.modulo, P.squareRoot, P.tan, P.tanh, P.toFraction,
 *                           P.toNearest, toStringBinary, naturalExponential, naturalLogarithm,
 *                           taylorSeries, atan2, parseOther
 *  finalise                 P.absoluteValue, P.atan, P.atanh, P.ceil, P.cos, P.cosh,
 *                           P.cubeRoot, P.dividedToIntegerBy, P.floor, P.logarithm, P.minus,
 *                           P.modulo, P.negated, P.plus, P.round, P.sin, P.sinh, P.squareRoot,
 *                           P.tan, P.times, P.toDecimalPlaces, P.toExponential, P.toFixed,
 *                           P.toNearest, P.toPower, P.toPrecision, P.toSignificantDigits,
 *                           P.truncated, divide, getLn10, getPi, naturalExponential,
 *                           naturalLogarithm, ceil, floor, round, trunc
 *  finiteToString           P.toExponential, P.toFixed, P.toPrecision, P.toString, P.valueOf,
 *                           toStringBinary
 *  getBase10Exponent        P.minus, P.plus, P.times, parseOther
 *  getLn10                  P.logarithm, naturalLogarithm
 *  getPi                    P.acos, P.asin, P.atan, toLessThanHalfPi, atan2
 *  getPrecision             P.precision, P.toFraction
 *  getZeroString            digitsToString, finiteToString
 *  intPow                   P.toPower, parseOther
 *  isOdd                    toLessThanHalfPi
 *  maxOrMin                 max, min
 *  naturalExponential       P.naturalExponential, P.toPower
 *  naturalLogarithm         P.acosh, P.asinh, P.atanh, P.logarithm, P.naturalLogarithm,
 *                           P.toPower, naturalExponential
 *  nonFiniteToString        finiteToString, toStringBinary
 *  parseDecimal             Decimal
 *  parseOther               Decimal
 *  sin                      P.sin
 *  taylorSeries             P.cosh, P.sinh, cos, sin
 *  toLessThanHalfPi         P.cos, P.sin
 *  toStringBinary           P.toBinary, P.toHexadecimal, P.toOctal
 *  truncate                 intPow
 *
 *  Throws:                  P.logarithm, P.precision, P.toFraction, checkInt32, getLn10, getPi,
 *                           naturalLogarithm, config, parseOther, random, Decimal
 */


function digitsToString(d) {
  var i, k, ws,
    indexOfLastWord = d.length - 1,
    str = '',
    w = d[0];

  if (indexOfLastWord > 0) {
    str += w;
    for (i = 1; i < indexOfLastWord; i++) {
      ws = d[i] + '';
      k = LOG_BASE - ws.length;
      if (k) str += getZeroString(k);
      str += ws;
    }

    w = d[i];
    ws = w + '';
    k = LOG_BASE - ws.length;
    if (k) str += getZeroString(k);
  } else if (w === 0) {
    return '0';
  }

  // Remove trailing zeros of last w.
  for (; w % 10 === 0;) w /= 10;

  return str + w;
}


function checkInt32(i, min, max) {
  if (i !== ~~i || i < min || i > max) {
    throw Error(invalidArgument + i);
  }
}


/*
 * Check 5 rounding digits if `repeating` is null, 4 otherwise.
 * `repeating == null` if caller is `log` or `pow`,
 * `repeating != null` if caller is `naturalLogarithm` or `naturalExponential`.
 */
function checkRoundingDigits(d, i, rm, repeating) {
  var di, k, r, rd;

  // Get the length of the first word of the array d.
  for (k = d[0]; k >= 10; k /= 10) --i;

  // Is the rounding digit in the first word of d?
  if (--i < 0) {
    i += LOG_BASE;
    di = 0;
  } else {
    di = Math.ceil((i + 1) / LOG_BASE);
    i %= LOG_BASE;
  }

  // i is the index (0 - 6) of the rounding digit.
  // E.g. if within the word 3487563 the first rounding digit is 5,
  // then i = 4, k = 1000, rd = 3487563 % 1000 = 563
  k = mathpow(10, LOG_BASE - i);
  rd = d[di] % k | 0;

  if (repeating == null) {
    if (i < 3) {
      if (i == 0) rd = rd / 100 | 0;
      else if (i == 1) rd = rd / 10 | 0;
      r = rm < 4 && rd == 99999 || rm > 3 && rd == 49999 || rd == 50000 || rd == 0;
    } else {
      r = (rm < 4 && rd + 1 == k || rm > 3 && rd + 1 == k / 2) &&
        (d[di + 1] / k / 100 | 0) == mathpow(10, i - 2) - 1 ||
          (rd == k / 2 || rd == 0) && (d[di + 1] / k / 100 | 0) == 0;
    }
  } else {
    if (i < 4) {
      if (i == 0) rd = rd / 1000 | 0;
      else if (i == 1) rd = rd / 100 | 0;
      else if (i == 2) rd = rd / 10 | 0;
      r = (repeating || rm < 4) && rd == 9999 || !repeating && rm > 3 && rd == 4999;
    } else {
      r = ((repeating || rm < 4) && rd + 1 == k ||
      (!repeating && rm > 3) && rd + 1 == k / 2) &&
        (d[di + 1] / k / 1000 | 0) == mathpow(10, i - 3) - 1;
    }
  }

  return r;
}


// Convert string of `baseIn` to an array of numbers of `baseOut`.
// Eg. convertBase('255', 10, 16) returns [15, 15].
// Eg. convertBase('ff', 16, 10) returns [2, 5, 5].
function convertBase(str, baseIn, baseOut) {
  var j,
    arr = [0],
    arrL,
    i = 0,
    strL = str.length;

  for (; i < strL;) {
    for (arrL = arr.length; arrL--;) arr[arrL] *= baseIn;
    arr[0] += NUMERALS.indexOf(str.charAt(i++));
    for (j = 0; j < arr.length; j++) {
      if (arr[j] > baseOut - 1) {
        if (arr[j + 1] === void 0) arr[j + 1] = 0;
        arr[j + 1] += arr[j] / baseOut | 0;
        arr[j] %= baseOut;
      }
    }
  }

  return arr.reverse();
}


/*
 * cos(x) = 1 - x^2/2! + x^4/4! - ...
 * |x| < pi/2
 *
 */
function cosine(Ctor, x) {
  var k, len, y;

  if (x.isZero()) return x;

  // Argument reduction: cos(4x) = 8*(cos^4(x) - cos^2(x)) + 1
  // i.e. cos(x) = 8*(cos^4(x/4) - cos^2(x/4)) + 1

  // Estimate the optimum number of times to use the argument reduction.
  len = x.d.length;
  if (len < 32) {
    k = Math.ceil(len / 3);
    y = (1 / tinyPow(4, k)).toString();
  } else {
    k = 16;
    y = '2.3283064365386962890625e-10';
  }

  Ctor.precision += k;

  x = taylorSeries(Ctor, 1, x.times(y), new Ctor(1));

  // Reverse argument reduction
  for (var i = k; i--;) {
    var cos2x = x.times(x);
    x = cos2x.times(cos2x).minus(cos2x).times(8).plus(1);
  }

  Ctor.precision -= k;

  return x;
}


/*
 * Perform division in the specified base.
 */
var divide = (function () {

  // Assumes non-zero x and k, and hence non-zero result.
  function multiplyInteger(x, k, base) {
    var temp,
      carry = 0,
      i = x.length;

    for (x = x.slice(); i--;) {
      temp = x[i] * k + carry;
      x[i] = temp % base | 0;
      carry = temp / base | 0;
    }

    if (carry) x.unshift(carry);

    return x;
  }

  function compare(a, b, aL, bL) {
    var i, r;

    if (aL != bL) {
      r = aL > bL ? 1 : -1;
    } else {
      for (i = r = 0; i < aL; i++) {
        if (a[i] != b[i]) {
          r = a[i] > b[i] ? 1 : -1;
          break;
        }
      }
    }

    return r;
  }

  function subtract(a, b, aL, base) {
    var i = 0;

    // Subtract b from a.
    for (; aL--;) {
      a[aL] -= i;
      i = a[aL] < b[aL] ? 1 : 0;
      a[aL] = i * base + a[aL] - b[aL];
    }

    // Remove leading zeros.
    for (; !a[0] && a.length > 1;) a.shift();
  }

  return function (x, y, pr, rm, dp, base) {
    var cmp, e, i, k, logBase, more, prod, prodL, q, qd, rem, remL, rem0, sd, t, xi, xL, yd0,
      yL, yz,
      Ctor = x.constructor,
      sign = x.s == y.s ? 1 : -1,
      xd = x.d,
      yd = y.d;

    // Either NaN, Infinity or 0?
    if (!xd || !xd[0] || !yd || !yd[0]) {

      return new Ctor(// Return NaN if either NaN, or both Infinity or 0.
        !x.s || !y.s || (xd ? yd && xd[0] == yd[0] : !yd) ? NaN :

        // Return ±0 if x is 0 or y is ±Infinity, or return ±Infinity as y is 0.
        xd && xd[0] == 0 || !yd ? sign * 0 : sign / 0);
    }

    if (base) {
      logBase = 1;
      e = x.e - y.e;
    } else {
      base = BASE;
      logBase = LOG_BASE;
      e = mathfloor(x.e / logBase) - mathfloor(y.e / logBase);
    }

    yL = yd.length;
    xL = xd.length;
    q = new Ctor(sign);
    qd = q.d = [];

    // Result exponent may be one less than e.
    // The digit array of a Decimal from toStringBinary may have trailing zeros.
    for (i = 0; yd[i] == (xd[i] || 0); i++);

    if (yd[i] > (xd[i] || 0)) e--;

    if (pr == null) {
      sd = pr = Ctor.precision;
      rm = Ctor.rounding;
    } else if (dp) {
      sd = pr + (x.e - y.e) + 1;
    } else {
      sd = pr;
    }

    if (sd < 0) {
      qd.push(1);
      more = true;
    } else {

      // Convert precision in number of base 10 digits to base 1e7 digits.
      sd = sd / logBase + 2 | 0;
      i = 0;

      // divisor < 1e7
      if (yL == 1) {
        k = 0;
        yd = yd[0];
        sd++;

        // k is the carry.
        for (; (i < xL || k) && sd--; i++) {
          t = k * base + (xd[i] || 0);
          qd[i] = t / yd | 0;
          k = t % yd | 0;
        }

        more = k || i < xL;

      // divisor >= 1e7
      } else {

        // Normalise xd and yd so highest order digit of yd is >= base/2
        k = base / (yd[0] + 1) | 0;

        if (k > 1) {
          yd = multiplyInteger(yd, k, base);
          xd = multiplyInteger(xd, k, base);
          yL = yd.length;
          xL = xd.length;
        }

        xi = yL;
        rem = xd.slice(0, yL);
        remL = rem.length;

        // Add zeros to make remainder as long as divisor.
        for (; remL < yL;) rem[remL++] = 0;

        yz = yd.slice();
        yz.unshift(0);
        yd0 = yd[0];

        if (yd[1] >= base / 2) ++yd0;

        do {
          k = 0;

          // Compare divisor and remainder.
          cmp = compare(yd, rem, yL, remL);

          // If divisor < remainder.
          if (cmp < 0) {

            // Calculate trial digit, k.
            rem0 = rem[0];
            if (yL != remL) rem0 = rem0 * base + (rem[1] || 0);

            // k will be how many times the divisor goes into the current remainder.
            k = rem0 / yd0 | 0;

            //  Algorithm:
            //  1. product = divisor * trial digit (k)
            //  2. if product > remainder: product -= divisor, k--
            //  3. remainder -= product
            //  4. if product was < remainder at 2:
            //    5. compare new remainder and divisor
            //    6. If remainder > divisor: remainder -= divisor, k++

            if (k > 1) {
              if (k >= base) k = base - 1;

              // product = divisor * trial digit.
              prod = multiplyInteger(yd, k, base);
              prodL = prod.length;
              remL = rem.length;

              // Compare product and remainder.
              cmp = compare(prod, rem, prodL, remL);

              // product > remainder.
              if (cmp == 1) {
                k--;

                // Subtract divisor from product.
                subtract(prod, yL < prodL ? yz : yd, prodL, base);
              }
            } else {

              // cmp is -1.
              // If k is 0, there is no need to compare yd and rem again below, so change cmp to 1
              // to avoid it. If k is 1 there is a need to compare yd and rem again below.
              if (k == 0) cmp = k = 1;
              prod = yd.slice();
            }

            prodL = prod.length;
            if (prodL < remL) prod.unshift(0);

            // Subtract product from remainder.
            subtract(rem, prod, remL, base);

            // If product was < previous remainder.
            if (cmp == -1) {
              remL = rem.length;

              // Compare divisor and new remainder.
              cmp = compare(yd, rem, yL, remL);

              // If divisor < new remainder, subtract divisor from remainder.
              if (cmp < 1) {
                k++;

                // Subtract divisor from remainder.
                subtract(rem, yL < remL ? yz : yd, remL, base);
              }
            }

            remL = rem.length;
          } else if (cmp === 0) {
            k++;
            rem = [0];
          }    // if cmp === 1, k will be 0

          // Add the next digit, k, to the result array.
          qd[i++] = k;

          // Update the remainder.
          if (cmp && rem[0]) {
            rem[remL++] = xd[xi] || 0;
          } else {
            rem = [xd[xi]];
            remL = 1;
          }

        } while ((xi++ < xL || rem[0] !== void 0) && sd--);

        more = rem[0] !== void 0;
      }

      // Leading zero?
      if (!qd[0]) qd.shift();
    }

    // logBase is 1 when divide is being used for base conversion.
    if (logBase == 1) {
      q.e = e;
      inexact = more;
    } else {

      // To calculate q.e, first get the number of digits of qd[0].
      for (i = 1, k = qd[0]; k >= 10; k /= 10) i++;
      q.e = i + e * logBase - 1;

      finalise(q, dp ? pr + q.e + 1 : pr, rm, more);
    }

    return q;
  };
})();


/*
 * Round `x` to `sd` significant digits using rounding mode `rm`.
 * Check for over/under-flow.
 */
 function finalise(x, sd, rm, isTruncated) {
  var digits, i, j, k, rd, roundUp, w, xd, xdi,
    Ctor = x.constructor;

  // Don't round if sd is null or undefined.
  out: if (sd != null) {
    xd = x.d;

    // Infinity/NaN.
    if (!xd) return x;

    // rd: the rounding digit, i.e. the digit after the digit that may be rounded up.
    // w: the word of xd containing rd, a base 1e7 number.
    // xdi: the index of w within xd.
    // digits: the number of digits of w.
    // i: what would be the index of rd within w if all the numbers were 7 digits long (i.e. if
    // they had leading zeros)
    // j: if > 0, the actual index of rd within w (if < 0, rd is a leading zero).

    // Get the length of the first word of the digits array xd.
    for (digits = 1, k = xd[0]; k >= 10; k /= 10) digits++;
    i = sd - digits;

    // Is the rounding digit in the first word of xd?
    if (i < 0) {
      i += LOG_BASE;
      j = sd;
      w = xd[xdi = 0];

      // Get the rounding digit at index j of w.
      rd = w / mathpow(10, digits - j - 1) % 10 | 0;
    } else {
      xdi = Math.ceil((i + 1) / LOG_BASE);
      k = xd.length;
      if (xdi >= k) {
        if (isTruncated) {

          // Needed by `naturalExponential`, `naturalLogarithm` and `squareRoot`.
          for (; k++ <= xdi;) xd.push(0);
          w = rd = 0;
          digits = 1;
          i %= LOG_BASE;
          j = i - LOG_BASE + 1;
        } else {
          break out;
        }
      } else {
        w = k = xd[xdi];

        // Get the number of digits of w.
        for (digits = 1; k >= 10; k /= 10) digits++;

        // Get the index of rd within w.
        i %= LOG_BASE;

        // Get the index of rd within w, adjusted for leading zeros.
        // The number of leading zeros of w is given by LOG_BASE - digits.
        j = i - LOG_BASE + digits;

        // Get the rounding digit at index j of w.
        rd = j < 0 ? 0 : w / mathpow(10, digits - j - 1) % 10 | 0;
      }
    }

    // Are there any non-zero digits after the rounding digit?
    isTruncated = isTruncated || sd < 0 ||
      xd[xdi + 1] !== void 0 || (j < 0 ? w : w % mathpow(10, digits - j - 1));

    // The expression `w % mathpow(10, digits - j - 1)` returns all the digits of w to the right
    // of the digit at (left-to-right) index j, e.g. if w is 908714 and j is 2, the expression
    // will give 714.

    roundUp = rm < 4
      ? (rd || isTruncated) && (rm == 0 || rm == (x.s < 0 ? 3 : 2))
      : rd > 5 || rd == 5 && (rm == 4 || isTruncated || rm == 6 &&

        // Check whether the digit to the left of the rounding digit is odd.
        ((i > 0 ? j > 0 ? w / mathpow(10, digits - j) : 0 : xd[xdi - 1]) % 10) & 1 ||
          rm == (x.s < 0 ? 8 : 7));

    if (sd < 1 || !xd[0]) {
      xd.length = 0;
      if (roundUp) {

        // Convert sd to decimal places.
        sd -= x.e + 1;

        // 1, 0.1, 0.01, 0.001, 0.0001 etc.
        xd[0] = mathpow(10, (LOG_BASE - sd % LOG_BASE) % LOG_BASE);
        x.e = -sd || 0;
      } else {

        // Zero.
        xd[0] = x.e = 0;
      }

      return x;
    }

    // Remove excess digits.
    if (i == 0) {
      xd.length = xdi;
      k = 1;
      xdi--;
    } else {
      xd.length = xdi + 1;
      k = mathpow(10, LOG_BASE - i);

      // E.g. 56700 becomes 56000 if 7 is the rounding digit.
      // j > 0 means i > number of leading zeros of w.
      xd[xdi] = j > 0 ? (w / mathpow(10, digits - j) % mathpow(10, j) | 0) * k : 0;
    }

    if (roundUp) {
      for (;;) {

        // Is the digit to be rounded up in the first word of xd?
        if (xdi == 0) {

          // i will be the length of xd[0] before k is added.
          for (i = 1, j = xd[0]; j >= 10; j /= 10) i++;
          j = xd[0] += k;
          for (k = 1; j >= 10; j /= 10) k++;

          // if i != k the length has increased.
          if (i != k) {
            x.e++;
            if (xd[0] == BASE) xd[0] = 1;
          }

          break;
        } else {
          xd[xdi] += k;
          if (xd[xdi] != BASE) break;
          xd[xdi--] = 0;
          k = 1;
        }
      }
    }

    // Remove trailing zeros.
    for (i = xd.length; xd[--i] === 0;) xd.pop();
  }

  if (external) {

    // Overflow?
    if (x.e > Ctor.maxE) {

      // Infinity.
      x.d = null;
      x.e = NaN;

    // Underflow?
    } else if (x.e < Ctor.minE) {

      // Zero.
      x.e = 0;
      x.d = [0];
      // Ctor.underflow = true;
    } // else Ctor.underflow = false;
  }

  return x;
}


function finiteToString(x, isExp, sd) {
  if (!x.isFinite()) return nonFiniteToString(x);
  var k,
    e = x.e,
    str = digitsToString(x.d),
    len = str.length;

  if (isExp) {
    if (sd && (k = sd - len) > 0) {
      str = str.charAt(0) + '.' + str.slice(1) + getZeroString(k);
    } else if (len > 1) {
      str = str.charAt(0) + '.' + str.slice(1);
    }

    str = str + (x.e < 0 ? 'e' : 'e+') + x.e;
  } else if (e < 0) {
    str = '0.' + getZeroString(-e - 1) + str;
    if (sd && (k = sd - len) > 0) str += getZeroString(k);
  } else if (e >= len) {
    str += getZeroString(e + 1 - len);
    if (sd && (k = sd - e - 1) > 0) str = str + '.' + getZeroString(k);
  } else {
    if ((k = e + 1) < len) str = str.slice(0, k) + '.' + str.slice(k);
    if (sd && (k = sd - len) > 0) {
      if (e + 1 === len) str += '.';
      str += getZeroString(k);
    }
  }

  return str;
}


// Calculate the base 10 exponent from the base 1e7 exponent.
function getBase10Exponent(digits, e) {
  var w = digits[0];

  // Add the number of digits of the first word of the digits array.
  for ( e *= LOG_BASE; w >= 10; w /= 10) e++;
  return e;
}


function getLn10(Ctor, sd, pr) {
  if (sd > LN10_PRECISION) {

    // Reset global state in case the exception is caught.
    external = true;
    if (pr) Ctor.precision = pr;
    throw Error(precisionLimitExceeded);
  }
  return finalise(new Ctor(LN10), sd, 1, true);
}


function getPi(Ctor, sd, rm) {
  if (sd > PI_PRECISION) throw Error(precisionLimitExceeded);
  return finalise(new Ctor(PI), sd, rm, true);
}


function getPrecision(digits) {
  var w = digits.length - 1,
    len = w * LOG_BASE + 1;

  w = digits[w];

  // If non-zero...
  if (w) {

    // Subtract the number of trailing zeros of the last word.
    for (; w % 10 == 0; w /= 10) len--;

    // Add the number of digits of the first word.
    for (w = digits[0]; w >= 10; w /= 10) len++;
  }

  return len;
}


function getZeroString(k) {
  var zs = '';
  for (; k--;) zs += '0';
  return zs;
}


/*
 * Return a new Decimal whose value is the value of Decimal `x` to the power `n`, where `n` is an
 * integer of type number.
 *
 * Implements 'exponentiation by squaring'. Called by `pow` and `parseOther`.
 *
 */
function intPow(Ctor, x, n, pr) {
  var isTruncated,
    r = new Ctor(1),

    // Max n of 9007199254740991 takes 53 loop iterations.
    // Maximum digits array length; leaves [28, 34] guard digits.
    k = Math.ceil(pr / LOG_BASE + 4);

  external = false;

  for (;;) {
    if (n % 2) {
      r = r.times(x);
      if (truncate(r.d, k)) isTruncated = true;
    }

    n = mathfloor(n / 2);
    if (n === 0) {

      // To ensure correct rounding when r.d is truncated, increment the last word if it is zero.
      n = r.d.length - 1;
      if (isTruncated && r.d[n] === 0) ++r.d[n];
      break;
    }

    x = x.times(x);
    truncate(x.d, k);
  }

  external = true;

  return r;
}


function isOdd(n) {
  return n.d[n.d.length - 1] & 1;
}


/*
 * Handle `max` (`n` is -1) and `min` (`n` is 1).
 */
function maxOrMin(Ctor, args, n) {
  var k, y,
    x = new Ctor(args[0]),
    i = 0;

  for (; ++i < args.length;) {
    y = new Ctor(args[i]);

    // NaN?
    if (!y.s) {
      x = y;
      break;
    }

    k = x.cmp(y);

    if (k === n || k === 0 && x.s === n) {
      x = y;
    }
  }

  return x;
}


/*
 * Return a new Decimal whose value is the natural exponential of `x` rounded to `sd` significant
 * digits.
 *
 * Taylor/Maclaurin series.
 *
 * exp(x) = x^0/0! + x^1/1! + x^2/2! + x^3/3! + ...
 *
 * Argument reduction:
 *   Repeat x = x / 32, k += 5, until |x| < 0.1
 *   exp(x) = exp(x / 2^k)^(2^k)
 *
 * Previously, the argument was initially reduced by
 * exp(x) = exp(r) * 10^k  where r = x - k * ln10, k = floor(x / ln10)
 * to first put r in the range [0, ln10], before dividing by 32 until |x| < 0.1, but this was
 * found to be slower than just dividing repeatedly by 32 as above.
 *
 * Max integer argument: exp('20723265836946413') = 6.3e+9000000000000000
 * Min integer argument: exp('-20723265836946411') = 1.2e-9000000000000000
 * (Math object integer min/max: Math.exp(709) = 8.2e+307, Math.exp(-745) = 5e-324)
 *
 *  exp(Infinity)  = Infinity
 *  exp(-Infinity) = 0
 *  exp(NaN)       = NaN
 *  exp(±0)        = 1
 *
 *  exp(x) is non-terminating for any finite, non-zero x.
 *
 *  The result will always be correctly rounded.
 *
 */
function naturalExponential(x, sd) {
  var denominator, guard, j, pow, sum, t, wpr,
    rep = 0,
    i = 0,
    k = 0,
    Ctor = x.constructor,
    rm = Ctor.rounding,
    pr = Ctor.precision;

  // 0/NaN/Infinity?
  if (!x.d || !x.d[0] || x.e > 17) {

    return new Ctor(x.d
      ? !x.d[0] ? 1 : x.s < 0 ? 0 : 1 / 0
      : x.s ? x.s < 0 ? 0 : x : 0 / 0);
  }

  if (sd == null) {
    external = false;
    wpr = pr;
  } else {
    wpr = sd;
  }

  t = new Ctor(0.03125);

  // while abs(x) >= 0.1
  while (x.e > -2) {

    // x = x / 2^5
    x = x.times(t);
    k += 5;
  }

  // Use 2 * log10(2^k) + 5 (empirically derived) to estimate the increase in precision
  // necessary to ensure the first 4 rounding digits are correct.
  guard = Math.log(mathpow(2, k)) / Math.LN10 * 2 + 5 | 0;
  wpr += guard;
  denominator = pow = sum = new Ctor(1);
  Ctor.precision = wpr;

  for (;;) {
    pow = finalise(pow.times(x), wpr, 1);
    denominator = denominator.times(++i);
    t = sum.plus(divide(pow, denominator, wpr, 1));

    if (digitsToString(t.d).slice(0, wpr) === digitsToString(sum.d).slice(0, wpr)) {
      j = k;
      while (j--) sum = finalise(sum.times(sum), wpr, 1);

      // Check to see if the first 4 rounding digits are [49]999.
      // If so, repeat the summation with a higher precision, otherwise
      // e.g. with precision: 18, rounding: 1
      // exp(18.404272462595034083567793919843761) = 98372560.1229999999 (should be 98372560.123)
      // `wpr - guard` is the index of first rounding digit.
      if (sd == null) {

        if (rep < 3 && checkRoundingDigits(sum.d, wpr - guard, rm, rep)) {
          Ctor.precision = wpr += 10;
          denominator = pow = t = new Ctor(1);
          i = 0;
          rep++;
        } else {
          return finalise(sum, Ctor.precision = pr, rm, external = true);
        }
      } else {
        Ctor.precision = pr;
        return sum;
      }
    }

    sum = t;
  }
}


/*
 * Return a new Decimal whose value is the natural logarithm of `x` rounded to `sd` significant
 * digits.
 *
 *  ln(-n)        = NaN
 *  ln(0)         = -Infinity
 *  ln(-0)        = -Infinity
 *  ln(1)         = 0
 *  ln(Infinity)  = Infinity
 *  ln(-Infinity) = NaN
 *  ln(NaN)       = NaN
 *
 *  ln(n) (n != 1) is non-terminating.
 *
 */
function naturalLogarithm(y, sd) {
  var c, c0, denominator, e, numerator, rep, sum, t, wpr, x1, x2,
    n = 1,
    guard = 10,
    x = y,
    xd = x.d,
    Ctor = x.constructor,
    rm = Ctor.rounding,
    pr = Ctor.precision;

  // Is x negative or Infinity, NaN, 0 or 1?
  if (x.s < 0 || !xd || !xd[0] || !x.e && xd[0] == 1 && xd.length == 1) {
    return new Ctor(xd && !xd[0] ? -1 / 0 : x.s != 1 ? NaN : xd ? 0 : x);
  }

  if (sd == null) {
    external = false;
    wpr = pr;
  } else {
    wpr = sd;
  }

  Ctor.precision = wpr += guard;
  c = digitsToString(xd);
  c0 = c.charAt(0);

  if (Math.abs(e = x.e) < 1.5e15) {

    // Argument reduction.
    // The series converges faster the closer the argument is to 1, so using
    // ln(a^b) = b * ln(a),   ln(a) = ln(a^b) / b
    // multiply the argument by itself until the leading digits of the significand are 7, 8, 9,
    // 10, 11, 12 or 13, recording the number of multiplications so the sum of the series can
    // later be divided by this number, then separate out the power of 10 using
    // ln(a*10^b) = ln(a) + b*ln(10).

    // max n is 21 (gives 0.9, 1.0 or 1.1) (9e15 / 21 = 4.2e14).
    //while (c0 < 9 && c0 != 1 || c0 == 1 && c.charAt(1) > 1) {
    // max n is 6 (gives 0.7 - 1.3)
    while (c0 < 7 && c0 != 1 || c0 == 1 && c.charAt(1) > 3) {
      x = x.times(y);
      c = digitsToString(x.d);
      c0 = c.charAt(0);
      n++;
    }

    e = x.e;

    if (c0 > 1) {
      x = new Ctor('0.' + c);
      e++;
    } else {
      x = new Ctor(c0 + '.' + c.slice(1));
    }
  } else {

    // The argument reduction method above may result in overflow if the argument y is a massive
    // number with exponent >= 1500000000000000 (9e15 / 6 = 1.5e15), so instead recall this
    // function using ln(x*10^e) = ln(x) + e*ln(10).
    t = getLn10(Ctor, wpr + 2, pr).times(e + '');
    x = naturalLogarithm(new Ctor(c0 + '.' + c.slice(1)), wpr - guard).plus(t);
    Ctor.precision = pr;

    return sd == null ? finalise(x, pr, rm, external = true) : x;
  }

  // x1 is x reduced to a value near 1.
  x1 = x;

  // Taylor series.
  // ln(y) = ln((1 + x)/(1 - x)) = 2(x + x^3/3 + x^5/5 + x^7/7 + ...)
  // where x = (y - 1)/(y + 1)    (|x| < 1)
  sum = numerator = x = divide(x.minus(1), x.plus(1), wpr, 1);
  x2 = finalise(x.times(x), wpr, 1);
  denominator = 3;

  for (;;) {
    numerator = finalise(numerator.times(x2), wpr, 1);
    t = sum.plus(divide(numerator, new Ctor(denominator), wpr, 1));

    if (digitsToString(t.d).slice(0, wpr) === digitsToString(sum.d).slice(0, wpr)) {
      sum = sum.times(2);

      // Reverse the argument reduction. Check that e is not 0 because, besides preventing an
      // unnecessary calculation, -0 + 0 = +0 and to ensure correct rounding -0 needs to stay -0.
      if (e !== 0) sum = sum.plus(getLn10(Ctor, wpr + 2, pr).times(e + ''));
      sum = divide(sum, new Ctor(n), wpr, 1);

      // Is rm > 3 and the first 4 rounding digits 4999, or rm < 4 (or the summation has
      // been repeated previously) and the first 4 rounding digits 9999?
      // If so, restart the summation with a higher precision, otherwise
      // e.g. with precision: 12, rounding: 1
      // ln(135520028.6126091714265381533) = 18.7246299999 when it should be 18.72463.
      // `wpr - guard` is the index of first rounding digit.
      if (sd == null) {
        if (checkRoundingDigits(sum.d, wpr - guard, rm, rep)) {
          Ctor.precision = wpr += guard;
          t = numerator = x = divide(x1.minus(1), x1.plus(1), wpr, 1);
          x2 = finalise(x.times(x), wpr, 1);
          denominator = rep = 1;
        } else {
          return finalise(sum, Ctor.precision = pr, rm, external = true);
        }
      } else {
        Ctor.precision = pr;
        return sum;
      }
    }

    sum = t;
    denominator += 2;
  }
}


// ±Infinity, NaN.
function nonFiniteToString(x) {
  // Unsigned.
  return String(x.s * x.s / 0);
}


/*
 * Parse the value of a new Decimal `x` from string `str`.
 */
function parseDecimal(x, str) {
  var e, i, len;

  // TODO BigInt str: no need to check for decimal point, exponential form or leading zeros.
  // Decimal point?
  if ((e = str.indexOf('.')) > -1) str = str.replace('.', '');

  // Exponential form?
  if ((i = str.search(/e/i)) > 0) {

    // Determine exponent.
    if (e < 0) e = i;
    e += +str.slice(i + 1);
    str = str.substring(0, i);
  } else if (e < 0) {

    // Integer.
    e = str.length;
  }

  // Determine leading zeros.
  for (i = 0; str.charCodeAt(i) === 48; i++);

  // Determine trailing zeros.
  for (len = str.length; str.charCodeAt(len - 1) === 48; --len);
  str = str.slice(i, len);

  if (str) {
    len -= i;
    x.e = e = e - i - 1;
    x.d = [];

    // Transform base

    // e is the base 10 exponent.
    // i is where to slice str to get the first word of the digits array.
    i = (e + 1) % LOG_BASE;
    if (e < 0) i += LOG_BASE;

    if (i < len) {
      if (i) x.d.push(+str.slice(0, i));
      for (len -= LOG_BASE; i < len;) x.d.push(+str.slice(i, i += LOG_BASE));
      str = str.slice(i);
      i = LOG_BASE - str.length;
    } else {
      i -= len;
    }

    for (; i--;) str += '0';
    x.d.push(+str);

    if (external) {

      // Overflow?
      if (x.e > x.constructor.maxE) {

        // Infinity.
        x.d = null;
        x.e = NaN;

      // Underflow?
      } else if (x.e < x.constructor.minE) {

        // Zero.
        x.e = 0;
        x.d = [0];
        // x.constructor.underflow = true;
      } // else x.constructor.underflow = false;
    }
  } else {

    // Zero.
    x.e = 0;
    x.d = [0];
  }

  return x;
}


/*
 * Parse the value of a new Decimal `x` from a string `str`, which is not a decimal value.
 */
function parseOther(x, str) {
  var base, Ctor, divisor, i, isFloat, len, p, xd, xe;

  if (str.indexOf('_') > -1) {
    str = str.replace(/(\d)_(?=\d)/g, '$1');
    if (isDecimal.test(str)) return parseDecimal(x, str);
  } else if (str === 'Infinity' || str === 'NaN') {
    if (!+str) x.s = NaN;
    x.e = NaN;
    x.d = null;
    return x;
  }

  if (isHex.test(str))  {
    base = 16;
    str = str.toLowerCase();
  } else if (isBinary.test(str))  {
    base = 2;
  } else if (isOctal.test(str))  {
    base = 8;
  } else {
    throw Error(invalidArgument + str);
  }

  // Is there a binary exponent part?
  i = str.search(/p/i);

  if (i > 0) {
    p = +str.slice(i + 1);
    str = str.substring(2, i);
  } else {
    str = str.slice(2);
  }

  // Convert `str` as an integer then divide the result by `base` raised to a power such that the
  // fraction part will be restored.
  i = str.indexOf('.');
  isFloat = i >= 0;
  Ctor = x.constructor;

  if (isFloat) {
    str = str.replace('.', '');
    len = str.length;
    i = len - i;

    // log[10](16) = 1.2041... , log[10](88) = 1.9444....
    divisor = intPow(Ctor, new Ctor(base), i, i * 2);
  }

  xd = convertBase(str, base, BASE);
  xe = xd.length - 1;

  // Remove trailing zeros.
  for (i = xe; xd[i] === 0; --i) xd.pop();
  if (i < 0) return new Ctor(x.s * 0);
  x.e = getBase10Exponent(xd, xe);
  x.d = xd;
  external = false;

  // At what precision to perform the division to ensure exact conversion?
  // maxDecimalIntegerPartDigitCount = ceil(log[10](b) * otherBaseIntegerPartDigitCount)
  // log[10](2) = 0.30103, log[10](8) = 0.90309, log[10](16) = 1.20412
  // E.g. ceil(1.2 * 3) = 4, so up to 4 decimal digits are needed to represent 3 hex int digits.
  // maxDecimalFractionPartDigitCount = {Hex:4|Oct:3|Bin:1} * otherBaseFractionPartDigitCount
  // Therefore using 4 * the number of digits of str will always be enough.
  if (isFloat) x = divide(x, divisor, len * 4);

  // Multiply by the binary exponent part if present.
  if (p) x = x.times(Math.abs(p) < 54 ? mathpow(2, p) : Decimal.pow(2, p));
  external = true;

  return x;
}


/*
 * sin(x) = x - x^3/3! + x^5/5! - ...
 * |x| < pi/2
 *
 */
function sine(Ctor, x) {
  var k,
    len = x.d.length;

  if (len < 3) {
    return x.isZero() ? x : taylorSeries(Ctor, 2, x, x);
  }

  // Argument reduction: sin(5x) = 16*sin^5(x) - 20*sin^3(x) + 5*sin(x)
  // i.e. sin(x) = 16*sin^5(x/5) - 20*sin^3(x/5) + 5*sin(x/5)
  // and  sin(x) = sin(x/5)(5 + sin^2(x/5)(16sin^2(x/5) - 20))

  // Estimate the optimum number of times to use the argument reduction.
  k = 1.4 * Math.sqrt(len);
  k = k > 16 ? 16 : k | 0;

  x = x.times(1 / tinyPow(5, k));
  x = taylorSeries(Ctor, 2, x, x);

  // Reverse argument reduction
  var sin2_x,
    d5 = new Ctor(5),
    d16 = new Ctor(16),
    d20 = new Ctor(20);
  for (; k--;) {
    sin2_x = x.times(x);
    x = x.times(d5.plus(sin2_x.times(d16.times(sin2_x).minus(d20))));
  }

  return x;
}


// Calculate Taylor series for `cos`, `cosh`, `sin` and `sinh`.
function taylorSeries(Ctor, n, x, y, isHyperbolic) {
  var j, t, u, x2,
    i = 1,
    pr = Ctor.precision,
    k = Math.ceil(pr / LOG_BASE);

  external = false;
  x2 = x.times(x);
  u = new Ctor(y);

  for (;;) {
    t = divide(u.times(x2), new Ctor(n++ * n++), pr, 1);
    u = isHyperbolic ? y.plus(t) : y.minus(t);
    y = divide(t.times(x2), new Ctor(n++ * n++), pr, 1);
    t = u.plus(y);

    if (t.d[k] !== void 0) {
      for (j = k; t.d[j] === u.d[j] && j--;);
      if (j == -1) break;
    }

    j = u;
    u = y;
    y = t;
    t = j;
    i++;
  }

  external = true;
  t.d.length = k + 1;

  return t;
}


// Exponent e must be positive and non-zero.
function tinyPow(b, e) {
  var n = b;
  while (--e) n *= b;
  return n;
}


// Return the absolute value of `x` reduced to less than or equal to half pi.
function toLessThanHalfPi(Ctor, x) {
  var t,
    isNeg = x.s < 0,
    pi = getPi(Ctor, Ctor.precision, 1),
    halfPi = pi.times(0.5);

  x = x.abs();

  if (x.lte(halfPi)) {
    quadrant = isNeg ? 4 : 1;
    return x;
  }

  t = x.divToInt(pi);

  if (t.isZero()) {
    quadrant = isNeg ? 3 : 2;
  } else {
    x = x.minus(t.times(pi));

    // 0 <= x < pi
    if (x.lte(halfPi)) {
      quadrant = isOdd(t) ? (isNeg ? 2 : 3) : (isNeg ? 4 : 1);
      return x;
    }

    quadrant = isOdd(t) ? (isNeg ? 1 : 4) : (isNeg ? 3 : 2);
  }

  return x.minus(pi).abs();
}


/*
 * Return the value of Decimal `x` as a string in base `baseOut`.
 *
 * If the optional `sd` argument is present include a binary exponent suffix.
 */
function toStringBinary(x, baseOut, sd, rm) {
  var base, e, i, k, len, roundUp, str, xd, y,
    Ctor = x.constructor,
    isExp = sd !== void 0;

  if (isExp) {
    checkInt32(sd, 1, MAX_DIGITS);
    if (rm === void 0) rm = Ctor.rounding;
    else checkInt32(rm, 0, 8);
  } else {
    sd = Ctor.precision;
    rm = Ctor.rounding;
  }

  if (!x.isFinite()) {
    str = nonFiniteToString(x);
  } else {
    str = finiteToString(x);
    i = str.indexOf('.');

    // Use exponential notation according to `toExpPos` and `toExpNeg`? No, but if required:
    // maxBinaryExponent = floor((decimalExponent + 1) * log[2](10))
    // minBinaryExponent = floor(decimalExponent * log[2](10))
    // log[2](10) = 3.321928094887362347870319429489390175864

    if (isExp) {
      base = 2;
      if (baseOut == 16) {
        sd = sd * 4 - 3;
      } else if (baseOut == 8) {
        sd = sd * 3 - 2;
      }
    } else {
      base = baseOut;
    }

    // Convert the number as an integer then divide the result by its base raised to a power such
    // that the fraction part will be restored.

    // Non-integer.
    if (i >= 0) {
      str = str.replace('.', '');
      y = new Ctor(1);
      y.e = str.length - i;
      y.d = convertBase(finiteToString(y), 10, base);
      y.e = y.d.length;
    }

    xd = convertBase(str, 10, base);
    e = len = xd.length;

    // Remove trailing zeros.
    for (; xd[--len] == 0;) xd.pop();

    if (!xd[0]) {
      str = isExp ? '0p+0' : '0';
    } else {
      if (i < 0) {
        e--;
      } else {
        x = new Ctor(x);
        x.d = xd;
        x.e = e;
        x = divide(x, y, sd, rm, 0, base);
        xd = x.d;
        e = x.e;
        roundUp = inexact;
      }

      // The rounding digit, i.e. the digit after the digit that may be rounded up.
      i = xd[sd];
      k = base / 2;
      roundUp = roundUp || xd[sd + 1] !== void 0;

      roundUp = rm < 4
        ? (i !== void 0 || roundUp) && (rm === 0 || rm === (x.s < 0 ? 3 : 2))
        : i > k || i === k && (rm === 4 || roundUp || rm === 6 && xd[sd - 1] & 1 ||
          rm === (x.s < 0 ? 8 : 7));

      xd.length = sd;

      if (roundUp) {

        // Rounding up may mean the previous digit has to be rounded up and so on.
        for (; ++xd[--sd] > base - 1;) {
          xd[sd] = 0;
          if (!sd) {
            ++e;
            xd.unshift(1);
          }
        }
      }

      // Determine trailing zeros.
      for (len = xd.length; !xd[len - 1]; --len);

      // E.g. [4, 11, 15] becomes 4bf.
      for (i = 0, str = ''; i < len; i++) str += NUMERALS.charAt(xd[i]);

      // Add binary exponent suffix?
      if (isExp) {
        if (len > 1) {
          if (baseOut == 16 || baseOut == 8) {
            i = baseOut == 16 ? 4 : 3;
            for (--len; len % i; len++) str += '0';
            xd = convertBase(str, base, baseOut);
            for (len = xd.length; !xd[len - 1]; --len);

            // xd[0] will always be be 1
            for (i = 1, str = '1.'; i < len; i++) str += NUMERALS.charAt(xd[i]);
          } else {
            str = str.charAt(0) + '.' + str.slice(1);
          }
        }

        str =  str + (e < 0 ? 'p' : 'p+') + e;
      } else if (e < 0) {
        for (; ++e;) str = '0' + str;
        str = '0.' + str;
      } else {
        if (++e > len) for (e -= len; e-- ;) str += '0';
        else if (e < len) str = str.slice(0, e) + '.' + str.slice(e);
      }
    }

    str = (baseOut == 16 ? '0x' : baseOut == 2 ? '0b' : baseOut == 8 ? '0o' : '') + str;
  }

  return x.s < 0 ? '-' + str : str;
}


// Does not strip trailing zeros.
function truncate(arr, len) {
  if (arr.length > len) {
    arr.length = len;
    return true;
  }
}


// Decimal methods


/*
 *  abs
 *  acos
 *  acosh
 *  add
 *  asin
 *  asinh
 *  atan
 *  atanh
 *  atan2
 *  cbrt
 *  ceil
 *  clamp
 *  clone
 *  config
 *  cos
 *  cosh
 *  div
 *  exp
 *  floor
 *  hypot
 *  ln
 *  log
 *  log2
 *  log10
 *  max
 *  min
 *  mod
 *  mul
 *  pow
 *  random
 *  round
 *  set
 *  sign
 *  sin
 *  sinh
 *  sqrt
 *  sub
 *  sum
 *  tan
 *  tanh
 *  trunc
 */


/*
 * Return a new Decimal whose value is the absolute value of `x`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function abs(x) {
  return new this(x).abs();
}


/*
 * Return a new Decimal whose value is the arccosine in radians of `x`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function acos(x) {
  return new this(x).acos();
}


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic cosine of `x`, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function acosh(x) {
  return new this(x).acosh();
}


/*
 * Return a new Decimal whose value is the sum of `x` and `y`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 * y {number|string|bigint|Decimal}
 *
 */
function add(x, y) {
  return new this(x).plus(y);
}


/*
 * Return a new Decimal whose value is the arcsine in radians of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function asin(x) {
  return new this(x).asin();
}


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic sine of `x`, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function asinh(x) {
  return new this(x).asinh();
}


/*
 * Return a new Decimal whose value is the arctangent in radians of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function atan(x) {
  return new this(x).atan();
}


/*
 * Return a new Decimal whose value is the inverse of the hyperbolic tangent of `x`, rounded to
 * `precision` significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function atanh(x) {
  return new this(x).atanh();
}


/*
 * Return a new Decimal whose value is the arctangent in radians of `y/x` in the range -pi to pi
 * (inclusive), rounded to `precision` significant digits using rounding mode `rounding`.
 *
 * Domain: [-Infinity, Infinity]
 * Range: [-pi, pi]
 *
 * y {number|string|bigint|Decimal} The y-coordinate.
 * x {number|string|bigint|Decimal} The x-coordinate.
 *
 * atan2(±0, -0)               = ±pi
 * atan2(±0, +0)               = ±0
 * atan2(±0, -x)               = ±pi for x > 0
 * atan2(±0, x)                = ±0 for x > 0
 * atan2(-y, ±0)               = -pi/2 for y > 0
 * atan2(y, ±0)                = pi/2 for y > 0
 * atan2(±y, -Infinity)        = ±pi for finite y > 0
 * atan2(±y, +Infinity)        = ±0 for finite y > 0
 * atan2(±Infinity, x)         = ±pi/2 for finite x
 * atan2(±Infinity, -Infinity) = ±3*pi/4
 * atan2(±Infinity, +Infinity) = ±pi/4
 * atan2(NaN, x) = NaN
 * atan2(y, NaN) = NaN
 *
 */
function atan2(y, x) {
  y = new this(y);
  x = new this(x);
  var r,
    pr = this.precision,
    rm = this.rounding,
    wpr = pr + 4;

  // Either NaN
  if (!y.s || !x.s) {
    r = new this(NaN);

  // Both ±Infinity
  } else if (!y.d && !x.d) {
    r = getPi(this, wpr, 1).times(x.s > 0 ? 0.25 : 0.75);
    r.s = y.s;

  // x is ±Infinity or y is ±0
  } else if (!x.d || y.isZero()) {
    r = x.s < 0 ? getPi(this, pr, rm) : new this(0);
    r.s = y.s;

  // y is ±Infinity or x is ±0
  } else if (!y.d || x.isZero()) {
    r = getPi(this, wpr, 1).times(0.5);
    r.s = y.s;

  // Both non-zero and finite
  } else if (x.s < 0) {
    this.precision = wpr;
    this.rounding = 1;
    r = this.atan(divide(y, x, wpr, 1));
    x = getPi(this, wpr, 1);
    this.precision = pr;
    this.rounding = rm;
    r = y.s < 0 ? r.minus(x) : r.plus(x);
  } else {
    r = this.atan(divide(y, x, wpr, 1));
  }

  return r;
}


/*
 * Return a new Decimal whose value is the cube root of `x`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function cbrt(x) {
  return new this(x).cbrt();
}


/*
 * Return a new Decimal whose value is `x` rounded to an integer using `ROUND_CEIL`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function ceil(x) {
  return finalise(x = new this(x), x.e + 1, 2);
}


/*
 * Return a new Decimal whose value is `x` clamped to the range delineated by `min` and `max`.
 *
 * x {number|string|bigint|Decimal}
 * min {number|string|bigint|Decimal}
 * max {number|string|bigint|Decimal}
 *
 */
function clamp(x, min, max) {
  return new this(x).clamp(min, max);
}


/*
 * Configure global settings for a Decimal constructor.
 *
 * `obj` is an object with one or more of the following properties,
 *
 *   precision  {number}
 *   rounding   {number}
 *   toExpNeg   {number}
 *   toExpPos   {number}
 *   maxE       {number}
 *   minE       {number}
 *   modulo     {number}
 *   crypto     {boolean|number}
 *   defaults   {true}
 *
 * E.g. Decimal.config({ precision: 20, rounding: 4 })
 *
 */
function config(obj) {
  if (!obj || typeof obj !== 'object') throw Error(decimalError + 'Object expected');
  var i, p, v,
    useDefaults = obj.defaults === true,
    ps = [
      'precision', 1, MAX_DIGITS,
      'rounding', 0, 8,
      'toExpNeg', -EXP_LIMIT, 0,
      'toExpPos', 0, EXP_LIMIT,
      'maxE', 0, EXP_LIMIT,
      'minE', -EXP_LIMIT, 0,
      'modulo', 0, 9
    ];

  for (i = 0; i < ps.length; i += 3) {
    if (p = ps[i], useDefaults) this[p] = DEFAULTS[p];
    if ((v = obj[p]) !== void 0) {
      if (mathfloor(v) === v && v >= ps[i + 1] && v <= ps[i + 2]) this[p] = v;
      else throw Error(invalidArgument + p + ': ' + v);
    }
  }

  if (p = 'crypto', useDefaults) this[p] = DEFAULTS[p];
  if ((v = obj[p]) !== void 0) {
    if (v === true || v === false || v === 0 || v === 1) {
      if (v) {
        if (typeof crypto != 'undefined' && crypto &&
          (crypto.getRandomValues || crypto.randomBytes)) {
          this[p] = true;
        } else {
          throw Error(cryptoUnavailable);
        }
      } else {
        this[p] = false;
      }
    } else {
      throw Error(invalidArgument + p + ': ' + v);
    }
  }

  return this;
}


/*
 * Return a new Decimal whose value is the cosine of `x`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function cos(x) {
  return new this(x).cos();
}


/*
 * Return a new Decimal whose value is the hyperbolic cosine of `x`, rounded to precision
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function cosh(x) {
  return new this(x).cosh();
}


/*
 * Create and return a Decimal constructor with the same configuration properties as this Decimal
 * constructor.
 *
 */
function clone(obj) {
  var i, p, ps;

  /*
   * The Decimal constructor and exported function.
   * Return a new Decimal instance.
   *
   * v {number|string|bigint|Decimal} A numeric value.
   *
   */
  function Decimal(v) {
    var e, i, t,
      x = this;

    // Decimal called without new.
    if (!(x instanceof Decimal)) return new Decimal(v);

    // Retain a reference to this Decimal constructor, and shadow Decimal.prototype.constructor
    // which points to Object.
    x.constructor = Decimal;

    if (isDecimalInstance(v)) {
      x.s = v.s;

      if (external) {
        if (!v.d || v.e > Decimal.maxE) {

          // Infinity.
          x.e = NaN;
          x.d = null;
        } else if (v.e < Decimal.minE) {

          // Zero.
          x.e = 0;
          x.d = [0];
        } else {
          x.e = v.e;
          x.d = v.d.slice();
        }
      } else {
        x.e = v.e;
        x.d = v.d ? v.d.slice() : v.d;
      }

      return;
    }

    t = typeof v;

    if (t === 'number') {
      if (v === 0) {
        x.s = 1 / v < 0 ? -1 : 1;
        x.e = 0;
        x.d = [0];
        return;
      }

      if (v < 0) {
        v = -v;
        x.s = -1;
      } else {
        x.s = 1;
      }

      // Fast path for small integers.
      if (v === ~~v && v < 1e7) {
        for (e = 0, i = v; i >= 10; i /= 10) e++;

        if (external) {
          if (e > Decimal.maxE) {
            x.e = NaN;
            x.d = null;
          } else if (e < Decimal.minE) {
            x.e = 0;
            x.d = [0];
          } else {
            x.e = e;
            x.d = [v];
          }
        } else {
          x.e = e;
          x.d = [v];
        }

        return;
      }

      // Infinity or NaN?
      if (v * 0 !== 0) {
        if (!v) x.s = NaN;
        x.e = NaN;
        x.d = null;
        return;
      }

      return parseDecimal(x, v.toString());
    }

    if (t === 'string') {
      if ((i = v.charCodeAt(0)) === 45) {  // minus sign
        v = v.slice(1);
        x.s = -1;
      } else {
        if (i === 43) v = v.slice(1);  // plus sign
        x.s = 1;
      }

      return isDecimal.test(v) ? parseDecimal(x, v) : parseOther(x, v);
    }

    if (t === 'bigint') {
      if (v < 0) {
        v = -v;
        x.s = -1;
      } else {
        x.s = 1;
      }

      return parseDecimal(x, v.toString());
    }

    throw Error(invalidArgument + v);
  }

  Decimal.prototype = P;

  Decimal.ROUND_UP = 0;
  Decimal.ROUND_DOWN = 1;
  Decimal.ROUND_CEIL = 2;
  Decimal.ROUND_FLOOR = 3;
  Decimal.ROUND_HALF_UP = 4;
  Decimal.ROUND_HALF_DOWN = 5;
  Decimal.ROUND_HALF_EVEN = 6;
  Decimal.ROUND_HALF_CEIL = 7;
  Decimal.ROUND_HALF_FLOOR = 8;
  Decimal.EUCLID = 9;

  Decimal.config = Decimal.set = config;
  Decimal.clone = clone;
  Decimal.isDecimal = isDecimalInstance;

  Decimal.abs = abs;
  Decimal.acos = acos;
  Decimal.acosh = acosh;        // ES6
  Decimal.add = add;
  Decimal.asin = asin;
  Decimal.asinh = asinh;        // ES6
  Decimal.atan = atan;
  Decimal.atanh = atanh;        // ES6
  Decimal.atan2 = atan2;
  Decimal.cbrt = cbrt;          // ES6
  Decimal.ceil = ceil;
  Decimal.clamp = clamp;
  Decimal.cos = cos;
  Decimal.cosh = cosh;          // ES6
  Decimal.div = div;
  Decimal.exp = exp;
  Decimal.floor = floor;
  Decimal.hypot = hypot;        // ES6
  Decimal.ln = ln;
  Decimal.log = log;
  Decimal.log10 = log10;        // ES6
  Decimal.log2 = log2;          // ES6
  Decimal.max = max;
  Decimal.min = min;
  Decimal.mod = mod;
  Decimal.mul = mul;
  Decimal.pow = pow;
  Decimal.random = random;
  Decimal.round = round;
  Decimal.sign = sign;          // ES6
  Decimal.sin = sin;
  Decimal.sinh = sinh;          // ES6
  Decimal.sqrt = sqrt;
  Decimal.sub = sub;
  Decimal.sum = sum;
  Decimal.tan = tan;
  Decimal.tanh = tanh;          // ES6
  Decimal.trunc = trunc;        // ES6

  if (obj === void 0) obj = {};
  if (obj) {
    if (obj.defaults !== true) {
      ps = ['precision', 'rounding', 'toExpNeg', 'toExpPos', 'maxE', 'minE', 'modulo', 'crypto'];
      for (i = 0; i < ps.length;) if (!obj.hasOwnProperty(p = ps[i++])) obj[p] = this[p];
    }
  }

  Decimal.config(obj);

  return Decimal;
}


/*
 * Return a new Decimal whose value is `x` divided by `y`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 * y {number|string|bigint|Decimal}
 *
 */
function div(x, y) {
  return new this(x).div(y);
}


/*
 * Return a new Decimal whose value is the natural exponential of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} The power to which to raise the base of the natural log.
 *
 */
function exp(x) {
  return new this(x).exp();
}


/*
 * Return a new Decimal whose value is `x` round to an integer using `ROUND_FLOOR`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function floor(x) {
  return finalise(x = new this(x), x.e + 1, 3);
}


/*
 * Return a new Decimal whose value is the square root of the sum of the squares of the arguments,
 * rounded to `precision` significant digits using rounding mode `rounding`.
 *
 * hypot(a, b, ...) = sqrt(a^2 + b^2 + ...)
 *
 * arguments {number|string|bigint|Decimal}
 *
 */
function hypot() {
  var i, n,
    t = new this(0);

  external = false;

  for (i = 0; i < arguments.length;) {
    n = new this(arguments[i++]);
    if (!n.d) {
      if (n.s) {
        external = true;
        return new this(1 / 0);
      }
      t = n;
    } else if (t.d) {
      t = t.plus(n.times(n));
    }
  }

  external = true;

  return t.sqrt();
}


/*
 * Return true if object is a Decimal instance (where Decimal is any Decimal constructor),
 * otherwise return false.
 *
 */
function isDecimalInstance(obj) {
  return obj instanceof Decimal || obj && obj.toStringTag === tag || false;
}


/*
 * Return a new Decimal whose value is the natural logarithm of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function ln(x) {
  return new this(x).ln();
}


/*
 * Return a new Decimal whose value is the log of `x` to the base `y`, or to base 10 if no base
 * is specified, rounded to `precision` significant digits using rounding mode `rounding`.
 *
 * log[y](x)
 *
 * x {number|string|bigint|Decimal} The argument of the logarithm.
 * y {number|string|bigint|Decimal} The base of the logarithm.
 *
 */
function log(x, y) {
  return new this(x).log(y);
}


/*
 * Return a new Decimal whose value is the base 2 logarithm of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function log2(x) {
  return new this(x).log(2);
}


/*
 * Return a new Decimal whose value is the base 10 logarithm of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function log10(x) {
  return new this(x).log(10);
}


/*
 * Return a new Decimal whose value is the maximum of the arguments.
 *
 * arguments {number|string|bigint|Decimal}
 *
 */
function max() {
  return maxOrMin(this, arguments, -1);
}


/*
 * Return a new Decimal whose value is the minimum of the arguments.
 *
 * arguments {number|string|bigint|Decimal}
 *
 */
function min() {
  return maxOrMin(this, arguments, 1);
}


/*
 * Return a new Decimal whose value is `x` modulo `y`, rounded to `precision` significant digits
 * using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 * y {number|string|bigint|Decimal}
 *
 */
function mod(x, y) {
  return new this(x).mod(y);
}


/*
 * Return a new Decimal whose value is `x` multiplied by `y`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 * y {number|string|bigint|Decimal}
 *
 */
function mul(x, y) {
  return new this(x).mul(y);
}


/*
 * Return a new Decimal whose value is `x` raised to the power `y`, rounded to precision
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} The base.
 * y {number|string|bigint|Decimal} The exponent.
 *
 */
function pow(x, y) {
  return new this(x).pow(y);
}


/*
 * Returns a new Decimal with a random value equal to or greater than 0 and less than 1, and with
 * `sd`, or `Decimal.precision` if `sd` is omitted, significant digits (or less if trailing zeros
 * are produced).
 *
 * [sd] {number} Significant digits. Integer, 0 to MAX_DIGITS inclusive.
 *
 */
function random(sd) {
  var d, e, k, n,
    i = 0,
    r = new this(1),
    rd = [];

  if (sd === void 0) sd = this.precision;
  else checkInt32(sd, 1, MAX_DIGITS);

  k = Math.ceil(sd / LOG_BASE);

  if (!this.crypto) {
    for (; i < k;) rd[i++] = Math.random() * 1e7 | 0;

  // Browsers supporting crypto.getRandomValues.
  } else if (crypto.getRandomValues) {
    d = crypto.getRandomValues(new Uint32Array(k));

    for (; i < k;) {
      n = d[i];

      // 0 <= n < 4294967296
      // Probability n >= 4.29e9, is 4967296 / 4294967296 = 0.00116 (1 in 865).
      if (n >= 4.29e9) {
        d[i] = crypto.getRandomValues(new Uint32Array(1))[0];
      } else {

        // 0 <= n <= 4289999999
        // 0 <= (n % 1e7) <= 9999999
        rd[i++] = n % 1e7;
      }
    }

  // Node.js supporting crypto.randomBytes.
  } else if (crypto.randomBytes) {

    // buffer
    d = crypto.randomBytes(k *= 4);

    for (; i < k;) {

      // 0 <= n < 2147483648
      n = d[i] + (d[i + 1] << 8) + (d[i + 2] << 16) + ((d[i + 3] & 0x7f) << 24);

      // Probability n >= 2.14e9, is 7483648 / 2147483648 = 0.0035 (1 in 286).
      if (n >= 2.14e9) {
        crypto.randomBytes(4).copy(d, i);
      } else {

        // 0 <= n <= 2139999999
        // 0 <= (n % 1e7) <= 9999999
        rd.push(n % 1e7);
        i += 4;
      }
    }

    i = k / 4;
  } else {
    throw Error(cryptoUnavailable);
  }

  k = rd[--i];
  sd %= LOG_BASE;

  // Convert trailing digits to zeros according to sd.
  if (k && sd) {
    n = mathpow(10, LOG_BASE - sd);
    rd[i] = (k / n | 0) * n;
  }

  // Remove trailing words which are zero.
  for (; rd[i] === 0; i--) rd.pop();

  // Zero?
  if (i < 0) {
    e = 0;
    rd = [0];
  } else {
    e = -1;

    // Remove leading words which are zero and adjust exponent accordingly.
    for (; rd[0] === 0; e -= LOG_BASE) rd.shift();

    // Count the digits of the first word of rd to determine leading zeros.
    for (k = 1, n = rd[0]; n >= 10; n /= 10) k++;

    // Adjust the exponent for leading zeros of the first word of rd.
    if (k < LOG_BASE) e -= LOG_BASE - k;
  }

  r.e = e;
  r.d = rd;

  return r;
}


/*
 * Return a new Decimal whose value is `x` rounded to an integer using rounding mode `rounding`.
 *
 * To emulate `Math.round`, set rounding to 7 (ROUND_HALF_CEIL).
 *
 * x {number|string|bigint|Decimal}
 *
 */
function round(x) {
  return finalise(x = new this(x), x.e + 1, this.rounding);
}


/*
 * Return
 *   1    if x > 0,
 *  -1    if x < 0,
 *   0    if x is 0,
 *  -0    if x is -0,
 *   NaN  otherwise
 *
 * x {number|string|bigint|Decimal}
 *
 */
function sign(x) {
  x = new this(x);
  return x.d ? (x.d[0] ? x.s : 0 * x.s) : x.s || NaN;
}


/*
 * Return a new Decimal whose value is the sine of `x`, rounded to `precision` significant digits
 * using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function sin(x) {
  return new this(x).sin();
}


/*
 * Return a new Decimal whose value is the hyperbolic sine of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function sinh(x) {
  return new this(x).sinh();
}


/*
 * Return a new Decimal whose value is the square root of `x`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function sqrt(x) {
  return new this(x).sqrt();
}


/*
 * Return a new Decimal whose value is `x` minus `y`, rounded to `precision` significant digits
 * using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal}
 * y {number|string|bigint|Decimal}
 *
 */
function sub(x, y) {
  return new this(x).sub(y);
}


/*
 * Return a new Decimal whose value is the sum of the arguments, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * Only the result is rounded, not the intermediate calculations.
 *
 * arguments {number|string|bigint|Decimal}
 *
 */
function sum() {
  var i = 0,
    args = arguments,
    x = new this(args[i]);

  external = false;
  for (; x.s && ++i < args.length;) x = x.plus(args[i]);
  external = true;

  return finalise(x, this.precision, this.rounding);
}


/*
 * Return a new Decimal whose value is the tangent of `x`, rounded to `precision` significant
 * digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function tan(x) {
  return new this(x).tan();
}


/*
 * Return a new Decimal whose value is the hyperbolic tangent of `x`, rounded to `precision`
 * significant digits using rounding mode `rounding`.
 *
 * x {number|string|bigint|Decimal} A value in radians.
 *
 */
function tanh(x) {
  return new this(x).tanh();
}


/*
 * Return a new Decimal whose value is `x` truncated to an integer.
 *
 * x {number|string|bigint|Decimal}
 *
 */
function trunc(x) {
  return finalise(x = new this(x), x.e + 1, 1);
}


P[Symbol.for('nodejs.util.inspect.custom')] = P.toString;
P[Symbol.toStringTag] = 'Decimal';

// Create and configure initial Decimal constructor.
var Decimal = P.constructor = clone(DEFAULTS);

// Create the internal constants from their string values.
LN10 = new Decimal(LN10);
PI = new Decimal(PI);



return Decimal;
})();

const precisionEngineCache=new Map();

function createPrecisionEngine(precision=75){
 if(precisionEngineCache.has(precision))return precisionEngineCache.get(precision);
 const D=Decimal.clone({precision,rounding:Decimal.ROUND_HALF_EVEN}),d=x=>new D(x);
 const C=(x=0,y=0)=>[d(x),d(y)],A=(a,b)=>[a[0].plus(b[0]),a[1].plus(b[1])],S=(a,b)=>[a[0].minus(b[0]),a[1].minus(b[1])];
 const M=(a,b)=>[a[0].times(b[0]).minus(a[1].times(b[1])),a[0].times(b[1]).plus(a[1].times(b[0]))];
 const K=(a,k)=>[a[0].times(k),a[1].times(k)],norm=a=>a[0].times(a[0]).plus(a[1].times(a[1])),abs=a=>norm(a).sqrt();
 const V=(a,b)=>K([a[0].times(b[0]).plus(a[1].times(b[1])),a[1].times(b[0]).minus(a[0].times(b[1]))],d(1).div(norm(b)));
 const P=(z,n)=>{if(n<0)return V(C(1),P(z,-n));let o=C(1),a=z;while(n){if(n%2)o=M(o,a);a=M(a,a);n=Math.floor(n/2);}return o;};
 const pi=D.acos(-1),twoPi=pi.times(2),one=C(1),threshold=d(10).pow(-precision+12),thresholdSquared=threshold.times(threshold);
 const exp=z=>{const v=z[0].exp();return[v.times(z[1].cos()),v.times(z[1].sin())];};
 const imaginaryExp=(real,imag)=>exp([d(imag).times(twoPi).neg(),d(real).times(twoPi)]);
 const asComplex=z=>Array.isArray(z)?C(z[0],z[1]):C(z.re,z.im);
 const inChart=t=>t[0].gt(0)&&t[0].lt(1)&&t[1].gt(0)&&norm(S(t,C('.5'))).gt('.25');
 function reduce(tau,r,s){
  let t=asComplex(tau),rr=d(r),ss=d(s),factor=one,ma=d(1),mb=d(0),mc=d(0),md=d(1);
  for(let k=0;k<100;k++){
   const integer=t[0].toNearest(1);t=[t[0].minus(integer),t[1]];rr=rr.plus(integer.times(ss));
   ma=ma.minus(integer.times(mc));mb=mb.minus(integer.times(md));
   rr=rr.minus(rr.toNearest(1));ss=ss.minus(ss.toNearest(1));
   if(norm(t).gte(d(1).minus(threshold)))return{t,r:rr,s:ss,factor,matrixC:mc};
   factor=M(t,factor);t=V(C(-1),t);[rr,ss]=[ss,rr.neg()];
   [ma,mb,mc,md]=[mc.neg(),md.neg(),ma,mb];
  }
  throw Error('Could not reduce the elliptic evaluation chart.');
 }
 function basis(tau,r,s,n=8,withDerivative=false){
  const chart=reduce(tau,r,s),{t}=chart; r=chart.r;s=chart.s;
  const ar=r.plus(s.times(t[0])),ai=s.times(t[1]),sign=ai.gte(0)?1:-1;
  const w=imaginaryExp(ar.times(sign),ai.abs()),den=S(w,one),quot=V(A(w,one),den);
  const cot=K([quot[1].neg(),quot[0]],sign),csc=K(V(w,M(den,den)),-4);
  let sin=C(),psum=C(),ppsum=C(),eisen2=C();
  let dsin=C(),dpsum=C(),dppsum=C();
  const timesI=(z,k)=>K([z[1].neg(),z[0]],k);
  const q=imaginaryExp(t[0],t[1]);
  const up=sign===1?M(q,w):V(q,w),down=sign===1?V(q,w):M(q,w);
  let qk=q,uk=up,vk=down;
  for(let k=1;k<400;k++){
   const denominator=S(one,qk),coeff=V(qk,denominator),plus=V(uk,denominator),minus=V(vk,denominator);
   const dif=S(plus,minus),sine=K([dif[1],dif[0].neg()],'.5'),cos=K(A(plus,minus),'.5');
   sin=A(sin,sine);psum=A(psum,K(S(coeff,cos),k));ppsum=A(ppsum,K(sine,k*k));
   eisen2=A(eisen2,K(coeff,k));
   if(withDerivative){
    // Differentiate the convergent Fourier series in the reduced chart.
    // r and s are fixed there; q' = 2 pi i q and exp(2 pi i z)' = 2 pi i s exp(2 pi i z).
    const dcoeff=timesI(V(coeff,denominator),twoPi.times(k));
    const dplus=timesI(M(plus,A(C(d(1).plus(s)),coeff)),twoPi.times(k));
    const dminus=timesI(M(minus,A(C(d(1).minus(s)),coeff)),twoPi.times(k));
    const dsine=K([S(dplus,dminus)[1],S(dplus,dminus)[0].neg()],'.5');
    const dcos=K(A(dplus,dminus),'.5');
    dsin=A(dsin,dsine);dpsum=A(dpsum,K(S(dcoeff,dcos),k));dppsum=A(dppsum,K(dsine,k*k));
   }
   // The equivalent squared tail test avoids three arbitrary-precision
   // square roots per Fourier term without weakening the stopping threshold.
   if(k>5&&D.max(norm(plus),norm(minus),norm(coeff)).times(k**10).lt(thresholdSquared))break;
   if(k===399)throw Error('High-precision Fourier series did not converge.');
   qk=M(qk,q);uk=M(uk,up);vk=M(vk,down);
  }
  const z=A(A(K(cot,pi),K(sin,pi.times(4))),C(0,twoPi.times(s)));
  const wp=A(K(S(csc,C(d(1).div(3))),pi.pow(2)),K(psum,pi.pow(2).times(8)));
  const pp=A(K(M(cot,csc),pi.pow(3).times(-2)),K(ppsum,pi.pow(3).times(16)));
  const eta1=K(S(one,K(eisen2,24)),pi.pow(2).div(3));
  let derivative;
  if(withDerivative){
   const dcot=K(csc,pi.times(s).neg()),dcsc=K(M(cot,csc),twoPi.times(s).neg());
   derivative={z:A(K(dcot,pi),K(dsin,pi.times(4))),
    wp:A(K(dcsc,pi.pow(2)),K(dpsum,pi.pow(2).times(8))),
    pp:A(K(A(M(dcot,csc),M(cot,dcsc)),pi.pow(3).times(-2)),K(dppsum,pi.pow(3).times(16)))};
  }
  if(n<=2)return{z,wp,pp,eta1,factor:chart.factor,matrixC:chart.matrixC,derivative};
  let t3=one,t4=one,t2sum=one;
  let dt3=C(),dt4=C(),dt2sum=C();
  const sqrt=z=>{const radius=abs(z);return[radius.plus(z[0]).div(2).sqrt(),radius.minus(z[0]).div(2).sqrt().times(z[1].isNeg()?-1:1)];};
  const halfq=sqrt(q),eighthq=sqrt(sqrt(halfq));
  for(let k=1;k<40;k++){
   const term=P(halfq,k*k),term2=P(halfq,k*(k+1));t3=A(t3,K(term,2));t4=A(t4,K(term,k%2?-2:2));
   t2sum=A(t2sum,term2);
   if(withDerivative){
    const dt=timesI(term,pi.times(k*k));dt3=A(dt3,K(dt,2));dt4=A(dt4,K(dt,k%2?-2:2));
    dt2sum=A(dt2sum,timesI(term2,pi.times(k*(k+1))));
   }
   if(k>3&&norm(term).lt(thresholdSquared))break;
  }
  const t2=K(M(eighthq,t2sum),2);
  const e1=K(A(P(t3,4),P(t4,4)),pi.pow(2).div(3));
  const e2=K(A(P(t2,4),P(t3,4)),pi.pow(2).div(-3));
  const e3=K(S(P(t2,4),P(t4,4)),pi.pow(2).div(3));
  if(withDerivative){
   const dt2=K(A(M(timesI(eighthq,pi.div(4)),t2sum),M(eighthq,dt2sum)),2);
   const p2=K(M(P(t2,3),dt2),4),p3=K(M(P(t3,3),dt3),4),p4=K(M(P(t4,3),dt4),4);
   derivative.e1=K(A(p3,p4),pi.pow(2).div(3));
   derivative.e2=K(A(p2,p3),pi.pow(2).div(-3));
   derivative.e3=K(S(p2,p4),pi.pow(2).div(3));derivative.D=S(derivative.e2,derivative.e1);
  }
  return{z,wp,pp,e1,e2,e3,D:S(e2,e1),factor:chart.factor,matrixC:chart.matrixC,derivative};
 }
 function fromBasis(b,n){
  if(n===1)return b.z;
  if(n===2)return S(S(P(b.z,3),K(M(b.wp,b.z),3)),b.pp);
  const T=V(S(b.e3,b.e1),b.D),Q={'-2':one,'-1':one,0:b.z};
  let G=V(b.pp,b.D),R=A(M(V(S(b.wp,b.e1),b.D),b.z),K(G,'.5'));
  for(let j=1;j<n;j++){
   const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1],phi=A(M(R,q2),K(G,d(j-1).div(2))),Y=M(M(q3,q2),q1),u=S(phi,Y),v=S(phi,M(T,Y));
   const H=M(M(phi,u),v),Hp=A(A(M(phi,u),M(phi,v)),M(u,v));
   const nextG=V(H,M(G,P(q3,3))),nextQ=S(V(M(q3,nextG),G),K(V(Hp,M(G,P(q3,2))),d(j).div(2)));
   Q[j]=nextQ;
   if(j<n-1)R=V(A(M(phi,nextQ),K(M(q3,nextG),d(j+1).div(2))),M(q3,q1));
   G=nextG;
  }
  const j=n-1,k=j%2?(j+1)**2/4:j*(j+2)/4,lead=K(P(b.D,-k),d(2).pow(-2*k));
  return V(Q[n-1],lead);
 }
 function value(tau,r,s,n){const b=basis(tau,r,s,n);return n>12?fromBasis(b,n):V(fromBasis(b,n),P(b.factor,n*(n+1)/2));}
 function valueAndDerivative(tau,r,s,n){
  if(n===1)return baseValueAndDerivative(tau,r,s);
  const b=basis(tau,r,s,n,true),zero=C();
  // First-order complex arithmetic propagates the exact series derivative
  // through the same premodule recurrence used for the equation itself.
  const constant=z=>[z,zero],jet=k=>[b[k],b.derivative[k]],
   ja=(a,b)=>[A(a[0],b[0]),A(a[1],b[1])],js=(a,b)=>[S(a[0],b[0]),S(a[1],b[1])],
   jm=(a,b)=>[M(a[0],b[0]),A(M(a[1],b[0]),M(a[0],b[1]))],
   jk=(a,k)=>[K(a[0],k),K(a[1],k)],
   jv=(a,b)=>{const q=V(a[0],b[0]);return[q,V(S(a[1],M(q,b[1])),b[0])];},
   jp=(a,n)=>[P(a[0],n),K(M(P(a[0],n-1),a[1]),n)];
  const z=jet('z'),wp=jet('wp'),pp=jet('pp');let f;
  if(n===2)f=js(js(jp(z,3),jk(jm(wp,z),3)),pp);
  else{
   const e1=jet('e1'),e3=jet('e3'),delta=jet('D'),T=jv(js(e3,e1),delta),Q={'-2':constant(one),'-1':constant(one),0:z};
   let G=jv(pp,delta),R=ja(jm(jv(js(wp,e1),delta),z),jk(G,'.5'));
   for(let j=1;j<n;j++){
    const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1],phi=ja(jm(R,q2),jk(G,d(j-1).div(2))),Y=jm(jm(q3,q2),q1),u=js(phi,Y),v=js(phi,jm(T,Y));
    const H=jm(jm(phi,u),v),Hp=ja(ja(jm(phi,u),jm(phi,v)),jm(u,v));
    const nextG=jv(H,jm(G,jp(q3,3))),nextQ=js(jv(jm(q3,nextG),G),jk(jv(Hp,jm(G,jp(q3,2))),d(j).div(2)));
    Q[j]=nextQ;
    if(j<n-1)R=jv(ja(jm(phi,nextQ),jk(jm(q3,nextG),d(j+1).div(2))),jm(q3,q1));
    G=nextG;
   }
   const j=n-1,k=j%2?(j+1)**2/4:j*(j+2)/4,lead=jk(jp(delta,-k),d(2).pow(-2*k));
   f=jv(Q[n-1],lead);
  }
  const weight=n*(n+1)/2;
  if(n>12)return{value:f[0],derivative:V(f[1],P(b.factor,2))};
  return{value:V(f[0],P(b.factor,weight)),
   derivative:S(V(f[1],P(b.factor,weight+2)),K(V(f[0],P(b.factor,weight+1)),d(weight).times(b.matrixC)))};
 }
 function baseValueAndDerivative(tau,r,s){
  const b=basis(tau,r,s,1),sum=A(b.pp,K(M(A(b.wp,b.eta1),b.z),2));
  const reducedDerivative=K([sum[1].neg(),sum[0]],d(1).div(pi.times(4)));
  const derivative=S(V(reducedDerivative,P(b.factor,3)),K(V(b.z,P(b.factor,2)),b.matrixC));
  return{value:V(b.z,b.factor),derivative};
 }
 function correct(seed,r,s,n,{maxIterations=18,trustRadius=.3,tolerance='1e-26'}={}){
  let tau=asComplex(seed);const start=tau,tol=d(tolerance);
  for(let iteration=0;iteration<maxIterations;iteration++){
   let f,df;
   {const exact=valueAndDerivative(tau,r,s,n);f=exact.value;df=exact.derivative;}
   const step=V(f,df),error=abs(step);
   if(!error.isFinite())return{status:'precision_limited',iterations:iteration};
   if(error.lt(tol))return{status:'ok',tau,iterations:iteration,correction:error};
   let accepted=false;
   for(let back=0;back<18;back++){
    const trial=S(tau,K(step,d(2).pow(-back)));
    if(!inChart(trial)||abs(S(trial,start)).gt(trustRadius))continue;
    if(abs(value(trial,r,s,n)).lt(abs(f))){tau=trial;accepted=true;break;}
   }
   if(!accepted)return{status:'correction_failed',iterations:iteration,tau};
  }
  return{status:'iteration_limit',tau};
 }
 const engine={D,d,C,A,S,M,V,K,P,abs,inChart,basis,value,valueAndDerivative,baseValueAndDerivative,correct,precision,asComplex};
 precisionEngineCache.set(precision,engine);return engine;
}

const schurPredictor=(()=>{// Standalone Bloch-jet determinant. No website imports or mutations.
const C=x=>typeof x==='number'?[x,0]:x;
const add=(...xs)=>xs.reduce((a,b)=>{b=C(b);return[a[0]+b[0],a[1]+b[1]]},[0,0]);
const neg=a=>[-a[0],-a[1]];
const sub=(a,b)=>add(a,neg(C(b)));
const mul=(a,b)=>{a=C(a);b=C(b);return[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]]};
const div=(a,b)=>{a=C(a);b=C(b);const d=b[0]**2+b[1]**2;return[(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d]};
const abs=a=>Math.hypot(...a);
const exp=a=>{const r=Math.exp(a[0]);return[r*Math.cos(a[1]),r*Math.sin(a[1])]};

function reduceModular(weights,tau,r,s){
 weights=[...weights];tau=[...tau];let matrix=[1,0,0,1];
 for(let steps=0;steps<60;steps++){
  const m=Math.round(tau[0]);
  if(m){
   tau[0]-=m;r+=m*s;
   if(Math.abs(m)%2)weights=[weights[0],weights[1],weights[3],weights[2]];
   matrix=[matrix[0]-m*matrix[2],matrix[1]-m*matrix[3],matrix[2],matrix[3]];
  }
  if(abs(tau)>=1-1e-13)break;
  tau=div(-1,tau);[r,s]=[s,-r];weights=[weights[0],weights[2],weights[1],weights[3]];
  matrix=[-matrix[2],-matrix[3],matrix[0],matrix[1]];
 }
 r-=Math.round(r);s-=Math.round(s);
 return {weights,tau,r,s,matrix};
}

// Taylor coefficients in v=2*pi*z; theta common q^(1/4) omitted.
function thetaJets(z,tau,N){
 z=C(z);const count=Math.max(4,Math.ceil(Math.sqrt(42/(Math.PI*tau[1]))+Math.abs(z[1])/tau[1])+3);
 const out=Array.from({length:N+1},()=>[0,0]);
 for(let n=0;n<count;n++){
  const k=n+.5,q=exp(mul([0,Math.PI*n*(n+1)],tau));
  const ep=exp(mul([0,2*Math.PI*k],z)),em=div(1,ep),factor=mul([0,-((-1)**n)],q);
  let fp=[1,0],fm=[1,0];
  for(let d=0;d<=N;d++){
   out[d]=add(out[d],mul(factor,sub(mul(fp,ep),mul(fm,em))));
   fp=mul(fp,[0,k/(d+1)]);fm=mul(fm,[0,-k/(d+1)]);
  }
 }
 return out;
}
function polyDiv(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=a[k]??[0,0];
  for(let j=1;j<=Math.min(k,b.length-1);j++)v=sub(v,mul(b[j],out[k-j]));
  out.push(div(v,b[0]));
 }
 return out;
}
function polyMul(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=[0,0];for(let j=Math.max(0,k-b.length+1);j<=Math.min(k,a.length-1);j++)v=add(v,mul(a[j],b[k-j]));
  out.push(v);
 }
 return out;
}
function choose(n,k){let v=1;for(let i=1;i<=k;i++)v=v*(n-k+i)/i;return v}

function tvMatrix(weights,tau,r,s,{reduce=true}={}){
 if(!(tau[1]>0))throw Error('Upper half-plane required');
 if(reduce)({weights,tau,r,s}=reduceModular(weights,tau,r,s));
 else {r-=Math.round(r);s-=Math.round(s)}
 const sites=[];
 for(let i=0;i<4;i++)for(let h=weights[i]-1;h>=0;h-=2)sites.push([i,h]);
 const degree=2*Math.max(0,...sites.map(x=>x[1])),N=degree+1;
 const alpha=add(r,mul(s,tau)),zero=thetaJets(0,tau,N+1),num=thetaJets(alpha,tau,N),den=zero.slice(1);
 const expjet=[[1,0]];
 for(let d=1;d<=N;d++)expjet.push(mul(expjet[d-1],[0,s/d]));
 const factor=div(zero[1],num[0]),diag=polyMul(expjet,polyDiv(num,den,N),N).map(x=>mul(x,factor));
 const bits=[[0,0],[1,0],[0,1],[1,1]],points=[[0,0],[.5,0],mul(.5,tau),mul(.5,add(1,tau))],jets=new Map();
 const M=sites.map(([i,l])=>sites.map(([j,h])=>{
  const d=l+h;
  if(i===j)return mul(choose(d,l),diag[d+1]);
  const du=bits[i][0]-bits[j][0],dv=bits[i][1]-bits[j][1],hu=(du+2)%2,hv=(dv+2)%2;
  const key=hu+2*hv,pointIndex=key,periodU=(du-hu)/2,periodV=(dv-hv)/2;
  if(!jets.has(key)){
   const z=points[pointIndex],scalar=mul(exp(mul([0,2*Math.PI*s],z)),factor);
   const jet=polyMul(expjet,polyDiv(thetaJets(add(z,alpha),tau,degree),thetaJets(z,tau,degree),degree),degree).map(x=>mul(x,scalar));
   jets.set(key,jet);
  }
  const periodFactor=exp([0,2*Math.PI*(s*periodU-r*periodV)]);
  return mul(choose(d,l),mul(periodFactor,jets.get(key)[d]));
 }));
 return {M,sites,weights,tau,r,s};
}
function determinant(weights,tau,r,s,options={}){
 // This shortcut selects a square-root character lift. It is NOT the full
 // zero locus of the larger-period equation. Enable only after the graph
 // branch and its lift have been calibrated; full generic is the default.
 if(options.isogeny && weights[0]===weights[1] && weights[2]===weights[3]){
  const [u,v]=options.isogenyLift??[0,0];
  if(weights[0]===weights[2])return determinant([weights[0],0,0,0],tau,(r+u)/2,(s+v)/2,{...options,isogeny:false});
  if(weights[2]===0)return determinant([weights[0],0,0,0],mul(2,tau),r,(s+v)/2,{...options,isogeny:false});
  return determinant([weights[0],weights[2],0,0],div(-1,mul(2,tau)),(s+v)/2,-r,{...options,isogeny:false});
 }
 const {M}=tvMatrix(weights,tau,r,s,options),N=M.length;
 let product=[1,0],minRelativePivot=1;
 // Row normalization keeps intermediate products in the double range.
 for(let i=0;i<N;i++){
  // A 1x1 row is the residual itself: dividing by its absolute value would
  // erase its zero and make secant iteration impossible.
  const scale=N===1?1:Math.max(...M[i].map(abs));
  if(!scale||!Number.isFinite(scale))return {value:[NaN,NaN],minRelativePivot:0,dimension:N};
  for(let j=0;j<N;j++)M[i][j]=mul(1/scale,M[i][j]);
 }
 for(let j=0;j<N;j++){
  let p=j;for(let i=j+1;i<N;i++)if(abs(M[i][j])>abs(M[p][j]))p=i;
  const pivot=abs(M[p][j]);minRelativePivot=Math.min(minRelativePivot,pivot);
  if(pivot<1e-300)return {value:[0,0],minRelativePivot:0,dimension:N};
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];product=neg(product)}
  product=mul(product,M[j][j]);
  for(let i=j+1;i<N;i++){
   const ratio=div(M[i][j],M[j][j]);
   for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]));
  }
 }
 return {value:product,minRelativePivot,dimension:N};
}

function solve(weights,r,s,guess,{tol=1e-11,maxsteps=32,isogeny=false,isogenyLift=[0,0]}={}){
 let a=guess,b=add(guess,[1e-6,2e-6]),fa=determinant(weights,a,r,s,{isogeny,isogenyLift}).value,fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 for(let i=0;i<maxsteps;i++){
  const denominator=sub(fb,fa);if(abs(denominator)===0)break;
  const step=div(mul(fb,sub(b,a)),denominator);
  if(abs(step)<tol)return {tau:b,iterations:i+2,residual:abs(fb)};
  a=b;fa=fb;b=sub(b,step);
  if(!Number.isFinite(abs(b))||b[1]<=.005)throw Error('Unsafe iterate');
  fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 }
 throw Error('No convergence');
}

// Fixed row/column equilibration and complete-pivot Schur residual.
// Unlike renormalizing a determinant at every iterate, this preserves a
// holomorphic residual and resolves the formerly unstable count-32 anchor.
function fixedSchur(E,get,initial){
 const abs=E.abs,div=E.div,mul=E.mul,sub=E.sub;
 const ref=get(initial),N=ref.length;if(!N)throw Error('Empty elliptic matching matrix');
 const greater=(a,b)=>a?.gt?a.gt(b):a>b;
 const maximum=a=>a.reduce((m,x)=>greater(abs(x),m)?abs(x):m,abs(a[0]));
 const scales=ref.map(maximum),M=ref.map((row,i)=>row.map(x=>div(x,scales[i])));
 const cs=Array.from({length:N},(_,j)=>maximum(M.map(row=>row[j])));
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=div(M[i][j],cs[j]);
 const rows=Array.from({length:N},(_,i)=>i),cols=[...rows];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=abs(M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const v=abs(M[i][k]);if(greater(v,best)){p=i;q=k;best=v}}
  if(!greater(best,0))throw Error('The Schur complement has a singular fixed pivot');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]]}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]]}
  for(let i=j+1;i<N;i++){const ratio=div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]))}
 }
 return{dimension:N,evaluate:z=>{
  const raw=get(z),A=rows.map(i=>cols.map(j=>div(div(raw[i][j],scales[i]),cs[j])));
  for(let j=0;j<N-1;j++)for(let i=j+1;i<N;i++){const ratio=div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=sub(A[i][k],mul(ratio,A[j][k]))}
  return A[N-1][N-1];
 }};
}

function predict(initial,n,r,s){
 try{
  const seed=[...initial],reduced=reduceModular([n,0,0,0],seed,r,s),[a,b,c,d]=reduced.matrix;
  const get=z=>tvMatrix([n,0,0,0],div(add(mul(a,z),b),add(mul(c,z),d)),reduced.r,reduced.s,{reduce:false}).M;
  const schur=fixedSchur({abs,div,mul,sub},get,seed);
  let tau=seed;
  for(let iteration=0;iteration<40;iteration++){
   const value=schur.evaluate(tau),h=1e-5;
   const derivative=div(sub(schur.evaluate(add(tau,[h,0])),schur.evaluate(sub(tau,[h,0]))),[2*h,0]);
   let step=div(value,derivative);if(!Number.isFinite(abs(step)))return null;
   if(abs(step)<3e-10)return tau;
   if(abs(step)>.06)step=mul(step,.06/abs(step));
   let accepted=false;
   for(let back=0;back<16;back++){
    const trial=sub(tau,mul(step,2**-back));
    if(!(trial[0]>0&&trial[0]<1&&trial[1]>0&&abs(sub(trial,.5))>.5)||abs(sub(trial,seed))>.24)continue;
    if(abs(schur.evaluate(trial))<abs(value)){tau=trial;accepted=true;break;}
   }
   if(!accepted)return null;
  }
  return null;
 }catch{return null;}
}

return{predict};})();
function transport(r,s,counts,previous,data,id){
 const started=performance.now(),n=1+counts.reduce((sum,x)=>sum+x,0),key=counts.join(',');
 const entry=data.anchors[key];
 if(!entry)return{id,n,status:'missing_anchor',message:'This labelled sheet has no geometric anchor.'};
 const p=createPrecisionEngine(45),v=createPrecisionEngine(70);
 const warm=previous?.status==='ok'&&previous?.result?.status==='ok'&&previous.counts?.every((x,i)=>x===counts[i]);
 const originR=warm?previous.r:data.r,originS=warm?previous.s:data.s;
 let tau=p.asComplex(warm?(previous.result.tauExact??previous.result.tau):(entry.highPrecisionTau??entry.tau));
 const dr=p.d(counts[0]+counts[2]).div(2),ds=p.d(counts[1]+counts[2]).div(2);
 const distance=Math.hypot(r-originR,s-originS);let fraction=0,increment=distance?Math.min(1,.02/distance):1,steps=0,iterations=0,last,priorTau=null,priorFraction=null;
 for(let attempts=0;fraction<1&&attempts<200;attempts++){
  if(performance.now()-started>4000)return{id,n,status:'time_budget',continuationSteps:steps,message:'The point budget was reached. The same labelled sheet can be retried.'};
  const next=Math.min(1,fraction+increment),rr=p.d(String(originR)).plus(p.d(String(r)).minus(String(originR)).times(String(next))).plus(dr),ss=p.d(String(originS)).plus(p.d(String(s)).minus(String(originS)).times(String(next))).plus(ds);
  let seed=tau;
  if(priorTau&&fraction>priorFraction){
   const ratio=(next-fraction)/(fraction-priorFraction);let movement=p.K(p.S(tau,priorTau),p.d(String(ratio)));
   if(ratio>0&&ratio<=3&&p.abs(movement).lt('.12'))seed=p.A(tau,movement);
  }
  if(distance*(next-fraction)>.0005){
   const proposed=schurPredictor.predict(seed.map(x=>x.toNumber()),n,rr.toNumber(),ss.toNumber());
   if(proposed&&p.abs(p.S(p.asComplex(proposed),tau)).lt('.16'))seed=p.asComplex(proposed);
  }
  last=p.correct(seed,rr,ss,n,{tolerance:'1e-18',trustRadius:.23,maxIterations:20});
  iterations+=last.iterations??0;
  if(last.status==='ok'){
   const movement=p.abs(p.S(last.tau,tau)).toNumber();
   if(movement>.16&&increment>1e-5){increment/=2;continue;}
   priorTau=tau;priorFraction=fraction;tau=last.tau;fraction=next;steps++;increment=Math.min(1-fraction,increment*(last.iterations>7?.8:1.4));
  }else{increment/=2;if(increment<1e-6)break;}
 }
 if(fraction!==1||last?.status!=='ok')return{id,n,status:last?.status??'unresolved',continuationSteps:steps,iterations,message:'This labelled branch did not pass the correction checks; no other root was substituted.'};
 const vr=v.d(String(r)).plus(dr.toString()),vs=v.d(String(s)).plus(ds.toString());
 const check=v.correct(tau.map(x=>x.toString()),vr,vs,n,{tolerance:'1e-28',trustRadius:.025,maxIterations:16});
 if(check.status!=='ok')return{id,n,status:'precision_limited',message:'The independent higher-precision point correction did not converge.'};
 const difference=v.abs(v.S(check.tau,v.asComplex(tau.map(x=>x.toString()))));
 if(difference.gt('2e-9'))return{id,n,status:'precision_limited',comparisonTauError:difference.toNumber(),message:'The two working precisions did not agree.'};
 const estimate=v.D.max(difference,check.correction,v.d('1e-30'));
 return{id,n,status:'ok',counts:[...counts],sheetKey:key,target:{r:r+dr.toNumber(),s:s+ds.toNumber()},tau:{re:check.tau[0].toNumber(),im:check.tau[1].toNumber()},tauExact:check.tau.map(x=>x.toPrecision(55)),correctionEstimate:estimate.toNumber(),comparisonTauError:difference.toNumber(),numericalErrorEstimate:estimate.toString(),iterations:iterations+check.iterations,continuationSteps:steps,method:'marked_continuation/scaled_multiprecision_45_70',validation:'two_precisions_not_interval_certified'};
}

return{transport};
})();

const atlasAttachment = (() => {
const {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh} = atlasMath;
/* Validated finite-range attachment engine. Geometry chooses the sheet;
 * numerical continuation transports it. No global root reseeding is used. */

// Smith division avoids squaring the tiny recurrence denominators. Such
// squares can underflow long before the quotient itself loses precision.
function div(a,b){
 if(Math.abs(b.re)>=Math.abs(b.im)){
  const q=b.im/b.re,d=b.re+b.im*q;
  return c((a.re+a.im*q)/d,(a.im-a.re*q)/d);
 }
 const q=b.re/b.im,d=b.im+b.re*q;
 return c((a.re*q+a.im)/d,(a.im*q-a.re)/d);
}

const MAX_N = 21;
const MAX_ATTACHMENTS = MAX_N - 1;
const DEFAULT_SOURCE={r:.31,s:.29};
const BASE_TAU=c(.5904953761444916,.7963943351098978);
let seedData=null;
const anchorCache=new Map();
const finite=z=>z&&Number.isFinite(z.re)&&Number.isFinite(z.im);
const zcopy=z=>Array.isArray(z)?c(z[0],z[1]):c(z.re,z.im);
const keyOf=counts=>counts.join(',');
const sameCounts=(a,b)=>Array.isArray(a)&&Array.isArray(b)&&a.length===3&&a.every((v,i)=>v===b[i]);
const pow=(z,n)=>{if(n<0)return div(c(1),pow(z,-n));let out=c(1),base=z;while(n){if(n&1)out=mul(out,base);base=mul(base,base);n=Math.floor(n/2);}return out;};

function validateCounts(counts,anchorData=seedData){
 if(!Array.isArray(counts)||counts.length!==3||!counts.every(v=>Number.isInteger(v)&&v>=0))
  return {status:'invalid_counts',message:'Enter three nonnegative whole numbers.'};
 const n=1+counts.reduce((a,b)=>a+b,0);
 const available=Math.min(MAX_N,anchorData?.maxN??anchorData?.max_n_completed??MAX_N);
 if(n>available)return {status:'unsupported_index',n,maxN:available,
  message:`This release validates at most ${available-1} attachments in total.`};
 return {status:'ok',n,maxN:available};
}

function setAnchorData(data){
 if(!data?.anchors||!Number.isFinite(data.r)||!Number.isFinite(data.s))
  throw new TypeError('Attachment anchors must include r, s and an anchors dictionary.');
 seedData=data;anchorCache.clear();
 return {status:'ok',maxN:Math.min(MAX_N,data.maxN??data.max_n_completed??MAX_N),
  anchorCount:Object.keys(data.anchors).length};
}

function ellipticBasis(tau,r,s){
 const e=evaluate(tau,r,s,1);
 let e1=evaluate(tau,.5,0,1).wp,e2=evaluate(tau,0,.5,1).wp,e3=evaluate(tau,.5,.5,1).wp;
 // The premodular polynomial is symmetric in the three half-period values.
 // Choose the largest Legendre denominator; a nearly coalescing pair makes
 // the rational recurrence grow although the underlying polynomial is finite.
 const d12=abs(sub(e2,e1)),d13=abs(sub(e3,e1)),d23=abs(sub(e3,e2));
 if(d13>=d12&&d13>=d23)[e2,e3]=[e3,e2];
 else if(d23>d12)[e1,e2,e3]=[e2,e3,e1];
 const D=sub(e2,e1),t=div(sub(e3,e1),D);
 const natural=Math.max(1,abs(e.z),Math.sqrt(abs(e.wp)),Math.cbrt(abs(e.wpPrime)),Math.sqrt(abs(D)));
 // The equation has weight n(n+1)/2. Normalize its weight-one scale BEFORE
 // the recurrence: large intermediate H and squared division denominators
 // otherwise overflow even while the final modulus is well conditioned.
 // This multiplies the equation by a nonzero real factor, preserving roots.
 const normalizedSize=1,inv=normalizedSize/natural,inv2=inv*inv;
 return {...e,z:scale(e.z,inv),wp:scale(e.wp,inv2),wpPrime:scale(e.wpPrime,inv2*inv),
  e1:scale(e1,inv2),e2:scale(e2,inv2),e3:scale(e3,inv2),D:scale(D,inv2),t,natural:normalizedSize,
  weightScale:natural/normalizedSize};
}
function leading(j,D){const k=j%2===0?j*(j+2)/4:(j+1)**2/4;return scale(pow(D,-k),2**(-2*k));}

function recurrence(b,n){
 const Q={'-2':c(1),'-1':c(1),0:b.z};
 let G=div(b.wpPrime,b.D),R=add(mul(div(sub(b.wp,b.e1),b.D),b.z),scale(G,.5));
 let suspect=false,minRelativePrior=Math.min(1,abs(b.z)/b.natural);
 for(let j=1;j<n;j++){
  const q3=Q[j-3],q2=Q[j-2],q1=Q[j-1];
  const phi=add(mul(R,q2),scale(G,(j-1)/2));
  const product=mul(mul(q3,q2),q1),u=sub(phi,product),v=sub(phi,mul(b.t,product));
  const H=mul(mul(phi,u),v),Hp=add(add(mul(phi,u),mul(phi,v)),mul(u,v));
  const nextG=div(H,mul(G,pow(q3,3)));
  const nextQ=sub(div(mul(q3,nextG),G),scale(div(Hp,mul(G,pow(q3,2))),j/2));
  Q[j]=nextQ;
  if(j<n-1){
   const first=mul(phi,nextQ),second=scale(mul(q3,nextG),(j+1)/2),numerator=add(first,second);
   const cancellation=abs(numerator)/Math.max(1e-300,abs(first)+abs(second));
   if(cancellation<1e-10&&minRelativePrior<1e-7)suspect=true;
   R=div(numerator,mul(q3,q1));
   const degree=(j+1)*(j+2)/2;
   const typical=abs(leading(j,b.D))*b.natural**degree;
   minRelativePrior=Math.min(minRelativePrior,abs(nextQ)/Math.max(1e-300,typical));
  }
  G=nextG;
 }
 if(n>=3&&minRelativePrior<1e-11)suspect=true;
 const value=div(Q[n-1],leading(n-1,b.D));
 return {value,suspect: suspect||!finite(value),minRelativePrior};
}

function explicitTerms(b,n){
 const Z=b.z,p=b.wp,pp=b.wpPrime,p2=mul(p,p),pp2=mul(pp,pp);
 const g2=scale(add(add(mul(b.e1,b.e2),mul(b.e1,b.e3)),mul(b.e2,b.e3)),-4);
 const g3=scale(mul(mul(b.e1,b.e2),b.e3),4);
 if(n===3)return[pow(Z,6),scale(mul(p,pow(Z,4)),-15),scale(mul(pp,pow(Z,3)),-20),
  mul(sub(scale(g2,27/4),scale(p2,45)),pow(Z,2)),scale(mul(mul(p,pp),Z),-12),scale(pp2,-5/4)];
 if(n===4)return[pow(Z,10),scale(mul(p,pow(Z,8)),-45),scale(mul(pp,pow(Z,7)),-120),
  mul(sub(scale(g2,399/4),scale(p2,630)),pow(Z,6)),scale(mul(mul(p,pp),pow(Z,5)),-504),
  scale(mul(sub(sub(scale(pow(p,3),280),scale(mul(g2,p),49)),scale(g3,115)),pow(Z,4)),-15/4),
  scale(mul(mul(sub(scale(g2,11),scale(p2,24)),pp),pow(Z,3)),15),
  scale(mul(add(add(sub(scale(pow(p,4),140),scale(mul(g2,p2),245)),scale(mul(g3,p),190)),scale(mul(g2,g2),21)),pow(Z,2)),-9/4),
  scale(mul(mul(add(sub(scale(pow(p,3),40),scale(mul(g2,p),163)),scale(g3,125)),pp),Z),-1),
  scale(mul(sub(scale(g2,25),scale(p2,3)),pp2),3/4)];
 return null;
}
function safelyAwayFromLattice(tau,r,s){
 const rr=r-Math.round(r),ss=s-Math.round(s);
 let shortest=Math.min(1,abs(tau),abs(sub(tau,c(1)))),distance=Infinity;
 for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)distance=Math.min(distance,abs(add(c(rr-i),scale(tau,ss-j))));
 return distance>.2*shortest;
}

/** Polynomial projection cancels removable recurrence singularities. It is a
 * fallback, not a universal replacement for the pointwise recurrence. */
function project(b,n,radiusRatio=.025){
 const degree=n*(n+1)/2,M=2**Math.ceil(Math.log2(degree+1)),radius=radiusRatio*b.natural;
 let sr=0,si=0,cr=0,ci=0,maxNode=0;
 for(let k=0;k<M;k++){
  const angle=2*Math.PI*(k+.1732050807568877)/M;
  const x=add(b.z,c(radius*Math.cos(angle),radius*Math.sin(angle)));
  const value=recurrence({...b,z:x},n).value;
  if(!finite(value))throw new RangeError('A polynomial-projection node lost precision.');
  maxNode=Math.max(maxNode,abs(value));
  let y=value.re-cr,t=sr+y;cr=(t-sr)-y;sr=t;
  y=value.im-ci;t=si+y;ci=(t-si)-y;si=t;
 }
 return {value:c(sr/M,si/M),maxNode,nodes:M,radius};
}

function evaluateAttachmentEquation(tau,r,s,n,options={}){
 if(!Number.isInteger(n)||n<1||n>MAX_N)throw new RangeError(`Supported indices are one through ${MAX_N}.`);
 tau=zcopy(tau);
 if(n<=2){const e=evaluate(tau,r,s,n);return {...e,method:n===1?'Hecke':'pole_subtracted_explicit',normalizer:e.residual/Math.max(e.relativeResidual,1e-300)};}
 const b=ellipticBasis(tau,r,s),raw=recurrence(b,n),degree=n*(n+1)/2;
 let value=raw.value,method='CKL_recurrence',normalizer=Math.max(1,b.natural**degree),projection=null;
 const explicit=explicitTerms(b,n),safe=safelyAwayFromLattice(tau,r,s);
 if(options.mode==='projection'){
  projection=project(b,n,options.radiusRatio??.025);value=projection.value;method='polynomial_projection';
 }else if(raw.suspect){
  if(explicit&&safe){value=explicit.reduce(add,c());normalizer=Math.max(1,explicit.reduce((a,v)=>a+abs(v),0));method='explicit_removed_division';}
  else{projection=project(b,n);value=projection.value;method='polynomial_projection';}
 }
 if(!finite(value))throw new RangeError('The equation could not be evaluated reliably.');
 return {value,residual:abs(value),relativeResidual:abs(value)/normalizer,normalizer,
  method,projectionNodes:projection?.nodes??0,recurrenceCondition:raw.minRelativePrior};
}

function differential(tau,r,s,n){
 const h=1e-5*Math.max(1,abs(tau));
 const plus=evaluateAttachmentEquation(c(tau.re+h,tau.im),r,s,n).value;
 const minus=evaluateAttachmentEquation(c(tau.re-h,tau.im),r,s,n).value;
 return scale(sub(plus,minus),1/(2*h));
}
function finalCheck(tau,r,s,n,e,derivative){
 if(n<=2)return{status:'ok',comparisonTauError:0};
 const b=ellipticBasis(tau,r,s),terms=explicitTerms(b,n);
 let difference;
 if(terms&&safelyAwayFromLattice(tau,r,s))difference=abs(sub(e.value,terms.reduce(add,c())));
 else{
  const a=project(b,n,.0125),d=project(b,n,.025);
  difference=Math.max(abs(sub(e.value,a.value)),abs(sub(a.value,d.value)));
 }
 const comparisonTauError=difference/Math.max(1e-300,abs(derivative));
 return {status:comparisonTauError<2e-8?'ok':'precision_limited',comparisonTauError};
}

function correct(tau,r,s,n,trustRadius=.2){
 tau=zcopy(tau);const seed=zcopy(tau);let evaluations=0;
 try{
  if(!inF0(tau))return{status:'outside_chart',tau,iterations:0};
  let e=evaluateAttachmentEquation(tau,r,s,n);evaluations++;
  for(let iteration=0;iteration<28;iteration++){
   const derivative=differential(tau,r,s,n);evaluations+=2;
   if(!finite(derivative)||abs(derivative)<1e-250)return{status:'small_derivative',tau,iterations:iteration};
   const step=div(e.value,derivative),correctionEstimate=abs(step);
   if(correctionEstimate<2e-10){
    const checked=finalCheck(tau,r,s,n,e,derivative);
    return {...e,...checked,tau,iterations:iteration,evaluations,correctionEstimate};
   }
   let accepted=false;
   for(let back=0;back<22;back++){
    const trial=sub(tau,scale(step,2**-back));
    if(!inF0(trial)||abs(sub(trial,seed))>trustRadius)continue;
    const next=evaluateAttachmentEquation(trial,r,s,n);evaluations++;
    if(next.residual<e.residual*(1-.00001*2**-back)){
     tau=trial;e=next;accepted=true;break;
    }
   }
   if(!accepted)return{status:'correction_failed',tau,iterations:iteration+1,evaluations};
  }
  return{status:'iteration_limit',tau,iterations:28,evaluations};
 }catch(error){return{status:'precision_limited',tau,iterations:0,evaluations,message:error.message};}
}

function initialAnchor(counts,data){
 const key=keyOf(counts);
 if(anchorCache.has(key))return anchorCache.get(key);
 if(counts.every(x=>x===0))return{r:DEFAULT_SOURCE.r,s:DEFAULT_SOURCE.s,tau:BASE_TAU,status:'validated_anchor'};
 const entry=data?.anchors?.[key];
 if(!entry)return null;
 const tau=entry.tau??entry.analyticTau??entry.fine?.tau;
 if(!tau)return null;
 const answer={r:data.r,s:data.s,tau:zcopy(tau),status:entry.status??'geometric_anchor'};
 anchorCache.set(key,answer);return answer;
}

function transport(r,s,counts,previous,data,id){
 const transportStarted=globalThis.performance?.now?.()??Date.now();
 const n=1+counts.reduce((a,b)=>a+b,0),dr=(counts[0]+counts[2])/2,ds=(counts[1]+counts[2])/2;
 if(n>12)return classicalPrecision.transport(r,s,counts,previous,data,id);
 const prior=id==='base'?previous?.base:(sameCounts(previous?.counts,counts)?previous?.result:null);
 const anchor=prior?.status==='ok'&&sourceStatus(previous.r,previous.s)==='ok'?
  {r:previous.r,s:previous.s,tau:prior.tau}:initialAnchor(counts,data);
 const target={r:r+dr,s:s+ds};
 if(!anchor)return{id,n,target,status:'missing_anchor',message:'The geometric sheet anchor is unavailable.'};
 let tau=zcopy(anchor.tau),fraction=0,steps=0,iterations=0,evaluations=0,last=null;
 const distance=Math.hypot(r-anchor.r,s-anchor.s);
 let increment=distance?Math.min(1,.022/distance):1;
 for(let attempts=0;fraction<1&&attempts<220;attempts++){
  if((globalThis.performance?.now?.()??Date.now())-transportStarted>350){
   last={status:'precision_limited'};break;
  }
  const next=Math.min(1,fraction+increment),rr=anchor.r+next*(r-anchor.r),ss=anchor.s+next*(s-anchor.s);
  last=correct(tau,rr+dr,ss+ds,n,Math.min(.23,Math.max(.035,.8*abs(tau))));
  iterations+=last.iterations??0;evaluations+=last.evaluations??0;
  if(last.status==='ok'){
   const movement=abs(sub(last.tau,tau));
   if(movement>.16&&increment>3e-5){increment*=.5;continue;}
   tau=last.tau;fraction=next;steps++;
   increment=Math.min(1-fraction,increment*(last.iterations>7?.8:1.35));
  }else{increment*=.5;if(increment<1e-5)break;}
 }
 if(fraction!==1||last?.status!=='ok')return{id,n,target,status:last?.status??'unresolved',
  iterations,evaluations,continuationSteps:steps,message:'This sheet could not be transported reliably; no other root was substituted.'};
 return {...last,id,n,target,tau,iterations,evaluations,continuationSteps:steps,
  counts:[...counts],sheetKey:keyOf(counts),method:`marked_continuation/${last.method}`};
}

/** Input counts refer to CWW's clockwise edge labels L1,L2,L3. */
function solveAttachment(r,s,counts,previousState=null,anchorData=null){
 const started=globalThis.performance?.now?.()??Date.now();
 const elapsed=()=> (globalThis.performance?.now?.()??Date.now())-started;
 if(anchorData&&anchorData!==seedData)setAnchorData(anchorData);
 const checked=validateCounts(counts,seedData),source=sourceStatus(r,s);
 if(checked.status!=='ok'||source!=='ok')return{r,s,counts:Array.isArray(counts)?[...counts]:counts,
  n:checked.n,status:checked.status!=='ok'?checked.status:source,branches:[],elapsedMs:elapsed(),
  message:checked.message??'Choose a point slightly inside the source triangle.'};
 let base=transport(r,s,[0,0,0],previousState,seedData,'base');
 if(base.status!=='ok'&&previousState?.base?.status==='ok')
  base={...transport(r,s,[0,0,0],null,seedData,'base'),restartedFromAnchor:true};
 let result=counts.every(x=>x===0)?{...base,id:'attached',counts:[0,0,0],sheetKey:'0,0,0'}:
  transport(r,s,counts,previousState,seedData,'attached');
 // A failed long warm path may restart only from this same sheet's fixed
 // geometric anchor. It must never borrow another sheet's analytic root.
 if(result.status!=='ok'&&sameCounts(previousState?.counts,counts)&&previousState?.result?.status==='ok')
  result={...transport(r,s,counts,null,seedData,'attached'),restartedFromAnchor:true};
 return{r,s,counts:[...counts],n:checked.n,status:base.status==='ok'&&result.status==='ok'?'ok':'partial',
  base,result,branches:[base,result],maxN:checked.maxN,elapsedMs:elapsed(),
  validation:{kind:'numerical_not_interval_certified',sourceMargin:.004,
   equationDegree:checked.n*(checked.n+1)/2,geometricAnchorKey:keyOf(counts)}};
}

async function sampleAttachment(options={}){
 if(options.anchorData)setAnchorData(options.anchorData);
 const counts=options.counts??[0,0,0],checked=validateCounts(counts),meshStarted=performance.now();
 if(checked.status!=='ok')return{status:checked.status,message:checked.message,counts,points:[],triangles:[],samples:[]};
 const mesh=buildSourceMesh(options.resolution??16,options.inset??.025),samples=new Array(mesh.points.length);
 // Alternate rows so continuation never jumps across the source triangle.
 const order=[];let cursor=0;
 for(let row=0;row<=mesh.resolution;row++){
  const indices=Array.from({length:mesh.resolution-row+1},()=>cursor++);
  order.push(...(row%2?indices.reverse():indices));
 }
 let previous=options.previousState??null,completed=0;
 for(const index of order){
  if(performance.now()-meshStarted>45000)break;
  if(options.signal?.aborted)return{...mesh,counts:[...counts],n:checked.n,samples,status:'aborted',completed};
  const p=mesh.points[index],sample=solveAttachment(p.r,p.s,counts,previous);
  samples[index]=sample;completed++;
  if(sample.status==='ok')previous=sample;
  if(options.onSample)await options.onSample(sample,index,completed,mesh.points.length);
  // Every point yields to queued interactive jobs in a worker.
  await new Promise(resolve=>setTimeout(resolve,0));
 }
 return{...mesh,counts:[...counts],n:checked.n,samples,completed,
  status:completed===mesh.points.length&&samples.every(s=>s?.status==='ok')?'ok':'partial',
  stopReason:completed<mesh.points.length?'time_budget':null,
  message:completed<mesh.points.length?'The 45-second mesh budget was reached. Checked points are retained; choose the quick mesh for a complete preview.':null,
  metadata:{sampledInteriorOnly:true,sourceCoordinateSideGap:mesh.inset/2,
   maximumRelativeResidual:Math.max(...samples.flatMap(s=>s.branches.filter(b=>b.status==='ok'&&Number.isFinite(b.relativeResidual)).map(b=>b.relativeResidual)))}};
}

return {MAX_N,validateCounts,setAnchorData,solveAttachment,sampleAttachment};
})();

const attachmentAnchors = {"r":0.31,"s":0.29,"maxN":21,"anchors":{"0,0,0":{"tau":[0.5904953761444915,0.7963943351098978],"highPrecisionTau":["0.59049537614449157416284192577349226908012778620242941899071673872","0.79639433510989779706922649871403435841001942722694363801613000114"]},"0,0,1":{"tau":[0.7515095138734595,0.6524433021803057],"highPrecisionTau":["0.75150951387345949645384410891502963838365256442799368836445356092","0.65244330218030567273735558173074768316920664224957470615138134323"]},"0,1,0":{"tau":[0.35448685090813264,0.6879997451009816],"highPrecisionTau":["0.35448685090813261961279501711880500203196199222685656403084028255","0.6879997451009816549969185915679019939473208000296978983723724225"]},"1,0,0":{"tau":[0.5823785379137855,1.1341826530853316],"highPrecisionTau":["0.58237853791378546366575728619209822223639032030552707907671234284","1.1341826530853315962403696846543904016643758273288518054099512232"]},"0,0,2":{"tau":[0.8103860604950234,0.5797764171757713],"highPrecisionTau":["0.8103860604950234175731630008717259854202578066353195524099023101","0.57977641717577128277353126631812503383246140420934427464566736476"]},"0,1,1":{"tau":[0.5400832084805506,0.646138229333322],"highPrecisionTau":["0.54008320848055053686221181766791611682008256560754518625148399504","0.64613822933332195609918022468290574950003947134280503981826771478"]},"0,2,0":{"tau":[0.2661530927667447,0.6170110529155463],"highPrecisionTau":["0.26615309276674471073268400652330747298467210193293483517054268986","0.61701105291554622938162772230971657851016530026401569330520689838"]},"1,0,1":{"tau":[0.7554851567438122,0.9044222083961427],"highPrecisionTau":["0.75548515674381221284180199628711477404486436987017071537132345048","0.90442220839614275355237867598541062254703419291474191429440448789"]},"1,1,0":{"tau":[0.33641698617555044,0.9361692214053354],"highPrecisionTau":["0.33641698617555046406636948191940451238830501503669948991252830308","0.93616922140533537502706857312216666531397406482493637466547587488"]},"2,0,0":{"tau":[0.5802309037096132,1.3519173892371166],"highPrecisionTau":["0.58023090370961317901706581078627485590837001710168383595215299674","1.3519173892371164981016024007805815496337166724563258239745718894"]},"0,0,3":{"tau":[0.841728094023958,0.5344069082972436],"highPrecisionTau":["0.84172809402395794437798963984200296198278994039341481977615430925","0.53440690829724356011926677948157995589544993666502620691599366763"]},"0,1,2":{"tau":[0.6299910352020642,0.6065667355063957],"highPrecisionTau":["0.62999103520206424666514873154487720771148631318073206502307391102","0.6065667355063956826424102709933489384973360426935705179468808005"]},"0,2,1":{"tau":[0.4336528965362644,0.6170485958998836],"highPrecisionTau":["0.43365289653626434940007447921084977524877221510063418025648297513","0.61704859589988364076591207303555004368890260561499544131519013268"]},"0,3,0":{"tau":[0.2194990271876052,0.5699249827074794],"highPrecisionTau":["0.21949902718760520093099855513898729088449236729422662883266496726","0.56992498270747942086610151461424252140108037271255773478539646094"]},"1,0,2":{"tau":[0.8190763243004178,0.788631148844644],"highPrecisionTau":["0.81907632430041772892328519180239142067540837434400890714722488429","0.78863114884464407776975156091075380608999434055926951769016416544"]},"1,1,1":{"tau":[0.5345911839088004,0.8413672802455601],"highPrecisionTau":["0.53459118390880048696444534995708777020751403305227293502385883776","0.84136728024556007641866915000768572689992943651461139960498303591"]},"1,2,0":{"tau":[0.24508363041652606,0.820541713936849],"highPrecisionTau":["0.24508363041652606668447818507255165119622906630612429387480143497","0.82054171393684906578724201430520828588042921667795591595711869001"]},"2,0,1":{"tau":[0.7566877827649808,1.0780300107149037],"highPrecisionTau":["0.75668778276498081344247991255363812083684005946317176155617978348","1.0780300107149035645178205761694741260449757878805112991332941373"]},"2,1,0":{"tau":[0.330728071930688,1.1089889089128258],"highPrecisionTau":["0.33072807193068796164601677265754057598583464328590510218069029164","1.1089889089128256834164626285410938482130621636056091194820439256"]},"3,0,0":{"tau":[0.5793596671311855,1.5133680043418907],"highPrecisionTau":["0.57935966713118549646785059370440192075542422181054615870192472471","1.5133680043418905691687554543054482224506433570565979892916896397"]},"0,0,4":{"tau":[0.8615254401911367,0.5026248506350188],"highPrecisionTau":["0.86152544019113679176238261036090362791599906106155313418191988599","0.50262485063501883243903202792683037354349133615644838533614111682"]},"0,1,3":{"tau":[0.6842316944717046,0.575021459683328],"highPrecisionTau":["0.68423169447170460977163152039873638600324844760931760601219117439","0.57502145968332792830429631206393700577703847815327820997756480831"]},"0,2,2":{"tau":[0.5267688263467021,0.6000381495242691],"highPrecisionTau":["0.52676882634670203761321702168536353822548443379150011373399405036","0.60003814952426916341894021800570857749777995150537650109461929878"]},"0,3,1":{"tau":[0.36881306450740986,0.5894886050882322],"highPrecisionTau":["0.36881306450740986600019944766522646292726181838047444235391206642","0.58948860508823214720729626181194754655378803182106548012591917938"]},"0,4,0":{"tau":[0.19030803642972227,0.5361316608784794],"highPrecisionTau":["0.19030803642972226195923039906606751094429657648058034784256067124","0.53613166087847942481143944943590013202057526501703482121631460924"]},"1,0,3":{"tau":[0.8526909974897778,0.7164178524340136],"highPrecisionTau":["0.85269099748977779302243528157438245855071283373839515788272474449","0.71641785243401359520687806770305032026740422441408457852575423185"]},"1,1,2":{"tau":[0.6314478387366457,0.7723654650396221],"highPrecisionTau":["0.63144783873664569897096624501800919017731117667232515026955020756","0.7723654650396220668988983506195457141834741745512073129968633332"]},"1,2,1":{"tau":[0.42284498422152705,0.7812601291251633],"highPrecisionTau":["0.42284498422152703173937940212605588027360784149563578338332795252","0.78126012912516329735344916647031757865170924841893912207017006811"]},"1,3,0":{"tau":[0.1973984260185703,0.7458993031656952],"highPrecisionTau":["0.19739842601857030394926260489959434252067563703764074971139778595","0.74589930316569522719944247157262638098458789659112178159557837239"]},"2,0,2":{"tau":[0.8220118840367037,0.9366554849150815],"highPrecisionTau":["0.82201188403670369127734153917229160521180002468082310926184341498","0.93665548491508158237731357882655183193651114764700543756265658965"]},"2,1,1":{"tau":[0.5326043104373224,0.9845389328712051],"highPrecisionTau":["0.5326043104373223345175149757944280437126504808483973725883114957","0.98453893287120513983018885997343556056975419330537722734883394418"]},"2,2,0":{"tau":[0.23787172685968883,0.9673218781297337],"highPrecisionTau":["0.23787172685968881552678510133220799921056458967044165201797200653","0.96732187812973369267972860672174032570951911869253887776187155073"]},"3,0,1":{"tau":[0.7572248859548885,1.2123248746321915],"highPrecisionTau":["0.75722488595488845319906528414650769917798664631152593431883823719","1.2123248746321916028050942702024680304975404206377622726359490042"]},"3,1,0":{"tau":[0.32814204362632093,1.2429897592333141],"highPrecisionTau":["0.32814204362632094045321879972854872701078605734163594623952594241","1.2429897592333140928882380021799093100987073416475871014131038247"]},"4,0,0":{"tau":[0.5789215292300838,1.64185870386297],"highPrecisionTau":["0.57892152923008382347305331105318138443897892382885769774757269725","1.6418587038629699722056853929650416131181086950589412496293698848"]},"0,0,5":{"tau":[0.87532736448948,0.47873158903313423],"highPrecisionTau":["0.87532736448947998989392931131756188132882542758413995587904279248","0.47873158903313425512431564199566686010994468294081736618363131118"]},"0,1,4":{"tau":[0.7209441947002367,0.5497539960286885],"highPrecisionTau":["0.72094419470023678379441655762213503456221642166883993470250605756","0.5497539960286884913533071351943236046348998525070588908654107854"]},"0,2,3":{"tau":[0.587659408841123,0.5813159221898846],"highPrecisionTau":["0.58765940884112298984260853113034659798112206512704164358904334387","0.58131592218988460152940071150247813259238731326111584413191558218"]},"0,3,2":{"tau":[0.45853298654130986,0.5864119195544156],"highPrecisionTau":["0.45853298654130985413521440762088757942486287528760971471361864912","0.58641191955441559775222197743266137078719362107241107908411946493"]},"0,4,1":{"tau":[0.32481896292121665,0.5658832421781945],"highPrecisionTau":["0.32481896292121666360821188628315927742123723995897777835699154832","0.56588324217819451154632017772738250071382839776731011019832529447"]},"0,5,0":{"tau":[0.17012729074846442,0.5104168701229547],"highPrecisionTau":["0.17012729074846442719830977101281023773357590692605419970522340727","0.51041687012295462879895022007481363011995661392148045537708790572"]},"1,0,4":{"tau":[0.8737146117761868,0.6659979038828044],"highPrecisionTau":["0.87371461177618683686792676045109496013823980693340644678345352783","0.66599790388280449506727724738527214440301055670187071001000178766"]},"1,1,3":{"tau":[0.6897768161485427,0.7213488596924934],"highPrecisionTau":["0.68977681614854273670238952302105090526885629109728594642841292682","0.72134885969249345991885241542859087965666040248495668524284454564"]},"1,2,2":{"tau":[0.5228133583527507,0.7416723275535352],"highPrecisionTau":["0.52281335835275071335372516715059654417468751928851430680670848517","0.74167232755353515308656362678206612559167239852240288279847073706"]},"1,3,1":{"tau":[0.3550011865750481,0.7333908261815947],"highPrecisionTau":["0.35500118657504811105054704307691854997202539677053553183165233342","0.73339082618159470258805910817218290282473827443673705600063241579"]},"1,4,0":{"tau":[0.16791151078300529,0.6931050520647739],"highPrecisionTau":["0.16791151078300528008349169871694909977370254364645526411680998573","0.69310505206477390360081980771087257096288705206867300177110263357"]},"2,0,3":{"tau":[0.856600950061261,0.8473150080377062],"highPrecisionTau":["0.85660095006126107603085759331448131591926957388724166238879057861","0.84731500803770616682614732436409075459340519026671566753033440705"]},"2,1,2":{"tau":[0.632157377106861,0.8969694583022713],"highPrecisionTau":["0.63215737710686109457554655443471075462615074021461477502987149979","0.89696945830227129525532275557527030842890861675509371485537132268"]},"2,2,1":{"tau":[0.41854260474439214,0.905436867652474],"highPrecisionTau":["0.41854260474439212256403006588805169214482855376088308305384923866","0.9054368676524740437391142374833323370833201913735984967688947774"]},"2,3,0":{"tau":[0.18945832036306007,0.875279890371322],"highPrecisionTau":["0.18945832036306005460419402300552539794039105226516813325044573107","0.87527989037132206200438493392664737788090466683572871410143389135"]},"3,0,2":{"tau":[0.8234166474721952,1.0536888952271188],"highPrecisionTau":["0.82341664747219527044209338600618231731810406559927273990913608647","1.0536888952271188137532989429423864912752877788734776048425626928"]},"3,1,1":{"tau":[0.5316145465016989,1.0991956202308],"highPrecisionTau":["0.53161454650169893714610633865980810574296630645367251874081216952","1.0991956202308000481404997115596830952381833549336505136119099664"]},"3,2,0":{"tau":[0.23439333176450308,1.0838565533760522],"highPrecisionTau":["0.23439333176450307736411955181909579165510046821841070533424084074","1.0838565533760521676096512597819470030794662973820065995755979604"]},"4,0,1":{"tau":[0.7575144695144952,1.3223411987468239],"highPrecisionTau":["0.75751446951449524024547708864000577406543475566745071038382009349","1.3223411987468237660502454926151864447380280705288792306802344496"]},"4,1,0":{"tau":[0.3267340104740559,1.3528628092637536],"highPrecisionTau":["0.32673401047405592949129418813456998877711332448204661597705588541","1.3528628092637535868084849863346567031313176593542605399277377056"]},"5,0,0":{"tau":[0.5786704237611745,1.7486344886564826],"highPrecisionTau":["0.57867042376117452490104457756525778178840133017895355965498943366","1.7486344886564825628276727408487155710443017132425650661551689344"]},"0,0,6":{"tau":[0.8855883749676,0.459893337978774],"highPrecisionTau":["0.88558837496760001689701039401265396872705556184812745947278918631","0.45989333797877403982308818181416148616913797206418153515844414766"]},"0,1,5":{"tau":[0.7476547390629024,0.5290741282185445],"highPrecisionTau":["0.74765473906290235728704380927307720678498484041861837770577659765","0.52907412821854450758465974785199194628435627206421227162885795137"]},"0,2,4":{"tau":[0.6310714701953587,0.5638367265401152],"highPrecisionTau":["0.63107147019535872176236397873568064523026928063100004154683640219","0.56383672654011517453895440387700648674074333888088718032256425625"]},"0,3,3":{"tau":[0.5204211943977096,0.5770842731764217],"highPrecisionTau":["0.52042119439770957192456156023878881487418254556011419283712549052","0.57708427317642170967626542290450102697494494132529032617857613792"]},"0,4,2":{"tau":[0.40965039211323434,0.5716922698326926],"highPrecisionTau":["0.40965039211323433337323945301435064843358891821475307217089056163","0.57169226983269256245385693415069817794105394707859884186279865714"]},"0,5,1":{"tau":[0.29281381542373114,0.5458775986357257],"highPrecisionTau":["0.2928138154237311125115774733454388486464267988743328631177230274","0.54587759863572566722113689732651636836904989182044236774066676159"]},"0,6,0":{"tau":[0.1552308493390186,0.49000353533730145],"highPrecisionTau":["0.15523084933901860067231531816440399926991551480064015099819265645","0.49000353533730144852416362697322401358068339166426027403075242072"]},"1,0,5":{"tau":[0.8882196634076015,0.6282589367010324],"highPrecisionTau":["0.88821966340760146444024861666919190735157592365813410524143767743","0.62825893670103244082411810756964209488335829231643078625987129953"]},"1,1,4":{"tau":[0.7290769007922019,0.6820042681833915],"highPrecisionTau":["0.72907690079220191578437785402494521803752743126677680463838969438","0.68200426818339153710338327159037902508138632118593445197064685026"]},"1,2,3":{"tau":[0.5881532493915438,0.7075651684523052],"highPrecisionTau":["0.58815324939154380106576707309851677133732149528119124291864145041","0.70756516845230517135711282291658655520618390872770731883257701078"]},"1,3,2":{"tau":[0.4510576512137605,0.7117896923823168],"highPrecisionTau":["0.45105765121376054069923593773251282600855024608173282334640725286","0.71178969238231680684338984336637680933874415010923337820968652717"]},"1,4,1":{"tau":[0.3092118506224589,0.695216006203829],"highPrecisionTau":["0.30921185062245885994015969729226740965302035887882430105729306513","0.69521600620382894610330745140008141569485121516791446485127257225"]},"1,5,0":{"tau":[0.1477519141514894,0.6533558004801285],"highPrecisionTau":["0.14775191415148940214730895203560951569170286495389917395531618832","0.65335580048012852092837195860983012506322195276232498978091689156"]},"2,0,4":{"tau":[0.8782133522035309,0.7844578175952746],"highPrecisionTau":["0.87821335220353094042726294393664672553156314471002550120013074007","0.78445781759527455599119007985660224061405900219129400475043842081"]},"2,1,3":{"tau":[0.6922991189465529,0.8329147799736228],"highPrecisionTau":["0.69229911894655296676861149736159614062421720434183465611635121365","0.83291477997362284981357378329660109760235747885753783710889870637"]},"2,2,2":{"tau":[0.5212230759002172,0.8511681141462283],"highPrecisionTau":["0.5212230759002172614306987826350007892286443672814369554169548737","0.85116811414622827757140159012376847577678375596443688605848424793"]},"2,3,1":{"tau":[0.34917986947047064,0.8442598735593233],"highPrecisionTau":["0.34917986947047062637633555303720313451862408185632408162946329534","0.84425987355932335021435375726319731171474916828709971677314908112"]},"2,4,0":{"tau":[0.15960719886127595,0.8098782385880396],"highPrecisionTau":["0.1596071988612759526338024768523755539977552356625024048095187703","0.80987823858803967823550265368919752342431644661615530130239384758"]},"3,0,3":{"tau":[0.8585464147633195,0.9521664394926173],"highPrecisionTau":["0.85854641476331944673754273460655262786745359767357990349315470301","0.95216643949261724572267354663318922313081007923077376617509340288"]},"3,1,2":{"tau":[0.6325417195306986,0.9985719839741891],"highPrecisionTau":["0.6325417195306986182747090801291883629494670865767526355407651868","0.99857198397418915262142482918692510875520880015575339466021101146"]},"3,2,1":{"tau":[0.4162831369003377,1.0068608466821334],"highPrecisionTau":["0.41628313690033770660882007957803817225646066870984727250026943391","1.006860846682133340413689058595090881435203678968653282849719714"]},"3,3,0":{"tau":[0.18548711432422435,0.979494178808416],"highPrecisionTau":["0.18548711432422435111335643795778531309896435686984842787575145588","0.9794941788084160989640805981540994458511397355180929070281671434"]},"4,0,2":{"tau":[0.8242124660453741,1.151180299071085],"highPrecisionTau":["0.82421246604537416282967221432402634619356765568250688524567187383","1.1511802990710850586240856728966466896535347463643848840578202466"]},"4,1,1":{"tau":[0.5310399425071837,1.1953234039601657],"highPrecisionTau":["0.53103994250718363680378579454853962668277413214754692536867270272","1.1953234039601656143849404359854521159840082767998228545372767591"]},"4,2,0":{"tau":[0.2324122394543271,1.1810929728116042],"highPrecisionTau":["0.23241223945432709546839307043893110858611149687303660229889116872","1.1810929728116041574948426309424751714771273342904699809303568356"]},"5,0,1":{"tau":[0.7576893917453332,1.415707662218239],"highPrecisionTau":["0.757689391745333247528234198525341799396657791035095473555384135","1.4157076622182388621098222179203795658677349003691051237424452893"]},"5,1,0":{"tau":[0.32587834922456627,1.4461479825119008],"highPrecisionTau":["0.32587834922456628166152581375478867700790983474920457069435268965","1.4461479825119008242495525393020983405971588814043943620480785326"]},"6,0,0":{"tau":[0.578513211933155,1.8400026346529397],"highPrecisionTau":["0.57851321193315493012856392864328060551480893868003335606331070938","1.8400026346529395740690358325599088985257271554023491653570396556"]},"0,0,7":{"tau":[0.8935695170631885,0.4445238550546686],"highPrecisionTau":["0.89356951706318849075508769891489852966763893799898237342261918614","0.44452385505466860120741290521925276103626487376772335258913648391"]},"0,1,6":{"tau":[0.7680816619446168,0.5117859055568753],"highPrecisionTau":["0.76808166194461680712201859559477977622690524208049922407970224875","0.51178590555687528760780020317802243616545967509102477315345214022"]},"0,2,5":{"tau":[0.6638087929376337,0.5481312045758051],"highPrecisionTau":["0.66380879293763374423740883452960200491604101268363617298727852572","0.54813120457580514691695397091290671917269913644764008484439164416"]},"0,3,4":{"tau":[0.5662502583558813,0.5660674216966612],"highPrecisionTau":["0.56625025835588122226192078988337018277252002008060330619115071576","0.56606742169666121806160811835941575812501913284662678530877988483"]},"0,4,3":{"tau":[0.47038431580611356,0.5691191152884885],"highPrecisionTau":["0.47038431580611354786360872914299885857360822710116304628900296598","0.56911911528848842356759777786567364167955659435964057048289114733"]},"0,5,2":{"tau":[0.37270851191302595,0.5575773347273169],"highPrecisionTau":["0.37270851191302597714723449691379776878768339656970990418084405079","0.55757733472731683212656249973699847580085809461235626080737651275"]},"0,6,1":{"tau":[0.2683647278635283,0.52879748719793],"highPrecisionTau":["0.26836472786352825943652763604182126894324183389907044806611236917","0.52879748719793000028373366987019175669553487877778934863617084826"]},"0,7,0":{"tau":[0.14371472207625965,0.4732798229549187],"highPrecisionTau":["0.14371472207625965452142122890917411932368895878540360538890123715","0.47327982295491869472361722979403071729528778884555460327658048054"]},"1,0,6":{"tau":[0.8988946785558809,0.5986474112765273],"highPrecisionTau":["0.89889467855588088476482507030608875844906881159172819300197611728","0.59864741127652726744888854782277050263463337549008763442081842092"]},"1,1,5":{"tau":[0.757515289092806,0.6505645672211018],"highPrecisionTau":["0.75751528909280597565043295356054253809747916470560228477010955912","0.65056456722110181225943745550001991687361517005313930485920498805"]},"1,2,4":{"tau":[0.6345895968123779,0.6786954967391343],"highPrecisionTau":["0.63458959681237787124775009447253622250312407612340239694020994685","0.67869549673913432174000036648635359401118017554777036927023171463"]},"1,3,3":{"tau":[0.5173362701255009,0.6896700070199766],"highPrecisionTau":["0.5173362701255008789636984007982471794580046293471055590681803961","0.68967000701997655653283838374336209148302550892942842229928132056"]},"1,4,2":{"tau":[0.39982624715111464,0.6851433733448178],"highPrecisionTau":["0.39982624715111466252482764020428321786024951551736380761213301904","0.68514337334481773526812707435432375182190702429812505495888157146"]},"1,5,1":{"tau":[0.276090469234877,0.6641428959307979],"highPrecisionTau":["0.27609046923487702449296912278986955822912752301445163465783154524","0.66414289593079784394398336748013410864438610897279802540444625927"]},"1,6,0":{"tau":[0.13302346902452356,0.6220753749174971],"highPrecisionTau":["0.1330234690245235593649361691312010034161408487946490611721650217","0.62207537491749704169103662037455469005599701762019976070917570568"]},"2,0,5":{"tau":[0.8930918423482282,0.7371942141027298],"highPrecisionTau":["0.89309184234822820384087079989037466992686501861297942938348765661","0.73719421410272978828101180246342420590830176787854697483933198902"]},"2,1,4":{"tau":[0.7328497447335313,0.7837854959610272],"highPrecisionTau":["0.73284974473353128358311006140309555124211936797982910951047604633","0.78378549596102716748266199763819421380300956444152877203360952113"]},"2,2,3":{"tau":[0.5885851008708306,0.8064731758072106],"highPrecisionTau":["0.58858510087083060333941607072529884486673386606847125504384227563","0.80647317580721053851314565755508170731956816583397503777328240864"]},"2,3,2":{"tau":[0.447737001043542,0.8104257477661121],"highPrecisionTau":["0.44773700104354204935935945504578860061428127216814639730175223981","0.81042574776611212830034045503909405378152042702940772164159040232"]},"2,4,1":{"tau":[0.3023745099048237,0.7961194599484742],"highPrecisionTau":["0.30237450990482365453672080171501620847786451082741108020321925467","0.79611945994847424258042172955922167626857329996595125285483857958"]},"2,5,0":{"tau":[0.13926952084402813,0.7604944279863987],"highPrecisionTau":["0.1392695208440281220643756277114309463142828721138326913463173153","0.76049442798639864973560588342119045920129937127135056002140654301"]},"3,0,4":{"tau":[0.8805111459535212,0.88014891120271],"highPrecisionTau":["0.88051114595352122095271303893219459268701620479921293821688917677","0.88014891120270998576037563733841175149821439615730837210067285841"]},"3,1,3":{"tau":[0.6936929117459214,0.9249420837481607],"highPrecisionTau":["0.69369291174592142370358703435042566341084249106166355422090702999","0.92494208374816068509230694985157704983128795812717153649895360852"]},"3,2,2":{"tau":[0.5203604410850786,0.9420472610477201],"highPrecisionTau":["0.52036044108507859103368682415975442196995487464893298896143931498","0.94204726104772003961591098857670758612623804713004590225271851619"]},"3,3,1":{"tau":[0.3460120862317194,0.9359838200637809],"highPrecisionTau":["0.34601208623171941889594885436728463319572829310694634678387610005","0.93598382006378090292257396590979919261184353144775735791927728069"]},"3,4,0":{"tau":[0.1553506593979241,0.9048369846987677],"highPrecisionTau":["0.15535065939792410344823686283246145703605557170120958895039749881","0.90483698469876767247081923853466180748106878882758301381509290058"]},"4,0,3":{"tau":[0.8596823845481332,1.040458892557427],"highPrecisionTau":["0.85968238454813314310993681891020962542934931004995104211142327842","1.040458892557426918942268055991296300689459012675318158110228267"]},"4,1,2":{"tau":[0.6327740811462541,1.0849285948795253],"highPrecisionTau":["0.63277408114625406112482043900783098950217704768036632484737604631","1.0849285948795252130086891797011699341213425384076514870570773624"]},"4,2,1":{"tau":[0.41492179738690005,1.0931249845327629],"highPrecisionTau":["0.41492179738690006817033915385071066717175334835301144470402042493","1.0931249845327628327127570471298683493023793817481643128156505858"]},"4,3,0":{"tau":[0.18315877131009872,1.0674498166149522],"highPrecisionTau":["0.1831587713100987090648297757547107449872057702355942301767512777","1.0674498166149522614176492249613412522448428854832789315950362611"]},"5,0,2":{"tau":[0.8247117287312125,1.2350094503071851],"highPrecisionTau":["0.82471172873121243404698884179436968293657667625424122986196646102","1.2350094503071850913706625684913210728359604519277619270930447363"]},"5,1,1":{"tau":[0.5306736884551937,1.2782904040218237],"highPrecisionTau":["0.53067368845519374423756622828343644977136066014540875220950512178","1.2782904040218237009963416085671048823809411191519655999477783855"]},"5,2,0":{"tau":[0.23116478113594885,1.264772313747665],"highPrecisionTau":["0.23116478113594885143927173333528984842476679352740649306009519574","1.2647723137476648760834420769044670699633605902752187034493342103"]},"6,0,1":{"tau":[0.7578034945743788,1.4968886398046868],"highPrecisionTau":["0.75780349457437882579292464108728269126231236265293503927381062953","1.4968886398046867075347062265150151550308777611294659576584592251"]},"6,1,0":{"tau":[0.3253179805035468,1.5272780873967173],"highPrecisionTau":["0.32531798050354681179920532379330847929799973615650668238463649851","1.5272780873967174303502841641138728735512841220412202572335486014"]},"7,0,0":{"tau":[0.5784082889971616,1.9198626778550154],"highPrecisionTau":["0.57840828899716160506181526297230049942770937340638262617192138181","1.9198626778550155068119237124633891221957876946091110865794277259"]},"0,0,8":{"tau":[0.8999886735985638,0.43165749057710767],"highPrecisionTau":["0.8999886735985638164935848730960737450247202396415483724266427115768604283764402816197","0.4316574905771076544906715328770268110639221134798922035891507255165463507430055272952"]},"0,1,7":{"tau":[0.7842843731758058,0.4970654360110489],"highPrecisionTau":["0.7842843731758058397309777888971647632761796975825942988830995551360814827084533385268","0.4970654360110489107820601546828928928563490804352351244219463002985153939018937886154"]},"0,2,6":{"tau":[0.6895018910795871,0.5341374020413163],"highPrecisionTau":["0.6895018910795871753760376912251855668524255078764198502779889149271066593023441955548","0.5341374020413163427617500712603758865632680914192855574974659526749740424487595705793"]},"0,3,5":{"tau":[0.6017919435613742,0.5549293341084078],"highPrecisionTau":["0.6017919435613742434020853736825321015854385747331274463028246695118022883579667339876","0.5549293341084078254887743623922555341929209879360764449086098403692562753231974765941"]},"0,4,4":{"tau":[0.5166510006237209,0.5631702795678296],"highPrecisionTau":["0.516651000623720883116475998713225599500741782219765151749194361921473673608631603524","0.5631702795678295673957416145579126596000363549551665845620286538133578210173621235154"]},"0,5,3":{"tau":[0.43147276876014107,0.5599209902678184],"highPrecisionTau":["0.4314727687601410501564668413641215256018386829718760328072355605671341499577179538788","0.5599209902678183605394205106091168376941611002634460400772134261481596244303357158151"]},"0,6,2":{"tau":[0.34369136170007636,0.5445324732829234],"highPrecisionTau":["0.3436913617000763447889912360811007432833242900553914487043378409227138580589179569249","0.544532473282923393526382737406269702929944984532267787048661614332148386939806384443"]},"0,7,1":{"tau":[0.24900061019025993,0.5140511474271879],"highPrecisionTau":["0.2490006101902599309323841379787908470778136960958963763044606387550853240204880717713","0.5140511474271879565954863614195891279316968821113849215669962668096775906539479356999"]},"0,8,0":{"tau":[0.1345010488513675,0.4592428537479833],"highPrecisionTau":["0.1345010488513674747064786752784890636372475293797493130126073449753058890730886863367","0.4592428537479832671404288337348455128546409697258791351046367885587198212544396588371"]},"1,0,7":{"tau":[0.9071185345899363,0.5746069889852019],"highPrecisionTau":["0.9071185345899363278372097646264172857443786689534955120386637716828501262636686437669","0.5746069889852018900102985434864218834269072382808042934493071311347777359092138986597"]},"1,1,6":{"tau":[0.7791416217273306,0.6247271194210997],"highPrecisionTau":["0.7791416217273306223622639562432323489698525995328403137073075629181197203253756367587","0.62472711942109970928824380957934603279079541275540804549001984800031191345361659353"]},"1,2,5":{"tau":[0.6694686003680538,0.6541175224491951],"highPrecisionTau":["0.6694686003680537390665232793364906056666339918438252431210586334030552931605840825071","0.6541175224491951145640216557642679196714271957653673390161109912225795352751794158835"]},"1,3,4":{"tau":[0.566291131542106,0.6690524526333079],"highPrecisionTau":["0.5662911315421059204394568817858508875384426412826816837973955712903484445083638077717","0.6690524526333079566909013323927705271924132196031198802999164839146771740778318689242"]},"1,4,3":{"tau":[0.464760493774037,0.6715543759066558],"highPrecisionTau":["0.4647604937740370227584414812700978795910449557712499174293962044279250807578501219229","0.6715543759066557788772843701403794451902709308749971146448698837122285558235126842945"]},"1,5,2":{"tau":[0.36126522625152746,0.6618033415718381],"highPrecisionTau":["0.3612652262515274668268158481593664780905506693431912865663737729555101013166657697171","0.6618033415718380212886265194553006401416139769100329039092307002449810382162505803362"]},"1,6,1":{"tau":[0.25093182248719836,0.6383140187936306],"highPrecisionTau":["0.2509318224871983618895031803697801432482922624221840048068551628779876080225521168318","0.6383140187936306350227618950471425062991285846748631712272145209047033796041581332055"]},"1,7,0":{"tau":[0.1217443668562202,0.5966418365607582],"highPrecisionTau":["0.1217443668562201924831472807841784691460039856086937993125849550515639824414951148891","0.5966418365607581937428214995963749352026852257197013559116158636284099306154131443286"]},"2,0,6":{"tau":[0.9040110611158669,0.7000080185294845],"highPrecisionTau":["0.9040110611158669192212899767339517092859912569672296711930418130194009771378797468059","0.7000080185294845207709644371088323823310194951350510459115530919567347543198532930885"]},"2,1,5":{"tau":[0.762178227429426,0.7446696257600844],"highPrecisionTau":["0.762178227429426008334650355692060459044584055421072031145214488255220027521091390815","0.744669625760084408492485931938331135369286863075035741130579751277945534181254931658"]},"2,2,4":{"tau":[0.6365164099219819,0.7694924688840932],"highPrecisionTau":["0.6365164099219819233115329832823904349680124416273358021939568759436547560854796837576","0.7694924688840931250272809559336781813396628888497276728559472137995345283082370502598"]},"2,3,3":{"tau":[0.5160314298136841,0.7793499753523696],"highPrecisionTau":["0.5160314298136840773287430455102130483058568777325736402203904188413316682973284364644","0.7793499753523696120236716195764260907048640834395330090223451491697699167631740106525"]},"2,4,2":{"tau":[0.3952342014435136,0.7754820407072519],"highPrecisionTau":["0.3952342014435136105768065107510190187285360266268170407099615843419028627271402523893","0.7754820407072519053336140151413559192355718141497807549592089822457884095823398926408"]},"2,5,1":{"tau":[0.2685551858462416,0.7572438938373648],"highPrecisionTau":["0.2685551858462416166091166510379829149774321566118218860283119901997337751780057639078","0.7572438938373648147238626332796820571442177928818211650736515791135103186416589479297"]},"2,6,0":{"tau":[0.12446556844497446,0.7215668307205353],"highPrecisionTau":["0.1244655684449744567335898344554386933336740338018613940627148813701859649233727997165","0.7215668307205353291108006898233156799601461323776097474733411934294358772204737073322"]},"3,0,5":{"tau":[0.8956278040791181,0.8256929702307054],"highPrecisionTau":["0.8956278040791181766362410671082806693878148925799627692802619835614552793284561300064","0.8256929702307054407364489756982487154259369643879653123157769864986685012231286177811"]},"3,1,4":{"tau":[0.734983197185718,0.8684124065360601],"highPrecisionTau":["0.734983197185717977415741206254667155314129814427745195094410575275568613134058207791","0.8684124065360600896404854753483447537996260339470407575185347424660983095523652841342"]},"3,2,3":{"tau":[0.5888624815053458,0.8894315983125628],"highPrecisionTau":["0.5888624815053457812234529859582162121033141947974769120714694195576542531956609343246","0.8894315983125628284818622787466645824845646961184169636893083724284428517030951880484"]},"3,3,2":{"tau":[0.4458530141184615,0.8932607935389507],"highPrecisionTau":["0.4458530141184614948909730378671069293379104423876682264704248273395439047754475125425","0.893260793538950667945861595574269455057507511931454750488622554087194283182796446784"]},"3,4,1":{"tau":[0.2985581051723347,0.8803522880710194],"highPrecisionTau":["0.298558105172334699699348692838903794703226393210214394727204514924600933580655644289","0.880352288071019408639009384849866734341486675579347529845509799340532820017056554082"]},"3,5,0":{"tau":[0.13484413406375953,0.8481938901081237],"highPrecisionTau":["0.1348441340637595294487988027615264452300053552253394968450590649594523588427045319184","0.8481938901081236866354956739558237308734144436659394048455927806863835978804137603003"]},"4,0,4":{"tau":[0.8818814555448774,0.9613304545714055],"highPrecisionTau":["0.8818814555448774303777166962619584899489866777458802989156559383799143135179583120558","0.9613304545714055660604808543295819702836673929048036204792860925256891658809841513811"]},"4,1,3":{"tau":[0.694557975951264,1.0038902962802296],"highPrecisionTau":["0.6945579759512639533249365412249316318480445828129750872483750005482214033421956674357","1.003890296280229667582218742599672925347417775283102450071853425599964355055305625045"]},"4,2,2":{"tau":[0.5198248176004105,1.0202785499707214],"highPrecisionTau":["0.519824817600410494676432329601875250765541030744704431707671613941568899641873656566","1.020278549970721335552817636423128537116831120156337869872637584961899666910906207613"]},"4,3,1":{"tau":[0.3440525661399386,1.0147704677679885],"highPrecisionTau":["0.3440525661399386179504371460083575879998307370204140775545541480487990090823671661438","1.014770467767988566836792152646826632321566933477638846470435210795667768780335317458"]},"4,4,0":{"tau":[0.15280393431861738,0.9856218132307234],"highPrecisionTau":["0.1528039343186173880912209608097533022058734750092344218853444296829475871139170608958","0.9856218132307233698969403434951282830923849871186270759267338520071259966171781428124"]},"5,0,3":{"tau":[0.8604126867429136,1.1170602611367253],"highPrecisionTau":["0.8604126867429135954201254710346728968823570240960828252707248194540525927398456111579","1.117060261136725297497179952935359299579514849995811442334319078374838728592265070999"]},"5,1,2":{"tau":[0.632926040787829,1.1602697414617],"highPrecisionTau":["0.6329260407878289625512173544495867123425845453657077828799821023217010887650452831027","1.160269741461699995802906186250831789887795939748422852135522420300440656647329447624"]},"5,2,1":{"tau":[0.41402886114936516,1.1684116115104175],"highPrecisionTau":["0.4140288611493651403721220929424816174758680347367932385744845296641902121203350427462","1.168411611510417498322966742775506754480472624882525594052032778721043914882319927641"]},"5,3,0":{"tau":[0.18165721310513627,1.1438482089549935],"highPrecisionTau":["0.1816572131051362693410927088861038310729204404284640215573081050890424425080507494477","1.143848208954993511769399155950813139119190596055360096171194247232730072972953191538"]},"6,0,2":{"tau":[0.8250473662854825,1.308670545276684],"highPrecisionTau":["0.8250473662854824877592757110301806364515649933252744358139354774599049200077009604375","1.308670545276683954951880192856140977382015735287680324880743042005027259195551576678"]},"6,1,1":{"tau":[0.5304247761233722,1.351368916988763],"highPrecisionTau":["0.5304247761233722929035515273117735734691063121038204917591072161672333160619649915843","1.351368916988762970770294745889279163084223251591872489570127655742303613668144595193"]},"6,2,0":{"tau":[0.23032393881576027,1.338337097595142],"highPrecisionTau":["0.2303239388157602807394385777262604174230204012561773770413000536344732021757366853968","1.338337097595142137458159423329153353790495744771964112051984353575422073599371549918"]},"7,0,1":{"tau":[0.7578821950635308,1.568741781822955],"highPrecisionTau":["0.757882195063530816243557963662140429799946864116175307844386971540554317879462319071","1.568741781822954938435694040696023428476506280080328571933415288907741482138038118804"]},"7,1,0":{"tau":[0.3249304151620226,1.599097138180921],"highPrecisionTau":["0.3249304151620226062317938303992132620307821655213726322458405249021832229235787045152","1.599097138180920899455858452199476350642370889496044188709018979357701815330376975143"]},"8,0,0":{"tau":[0.5783347884054936,1.9907964524404396],"highPrecisionTau":["0.5783347884054935449942169978418890387880981656433728541099263446594450869351617530501","1.990796452440439632300880746142202742461386479495388196960284304927881257137275945998"]},"0,0,9":{"tau":[0.9052863407274091,0.42066837982167543],"highPrecisionTau":["0.9052863407274090933574655407508129107675384959949349928159072893343066456665583657345","0.4206683798216754501742518183810390914028967073829707442709249865630389854379704428759"]},"0,1,8":{"tau":[0.7975003961237308,0.4843360211540446],"highPrecisionTau":["0.7975003961237307925185014346728138893363248495699354904329837911932262889758660018863","0.4843360211540446346988434748030080160755093550880718164236123623353545999493363737415"]},"0,2,7":{"tau":[0.7102810744526495,0.5216545187363899],"highPrecisionTau":["0.7102810744526495547836604462144826759260183481047176257615361329946656200108277121228","0.5216545187363899425828053261378549346801346650131571586022134506060360766383946306552"]},"0,3,6":{"tau":[0.6302881075104035,0.5442438663912446],"highPrecisionTau":["0.6302881075104035903229456678402202537373441938150574962530570038080310451433110931012","0.54424386639124453316982290351997029623721989185919928106236717211215366332431450678"]},"0,4,5":{"tau":[0.5533301666697178,0.5558690621223308],"highPrecisionTau":["0.5533301666697178126144911742925652820883795753852642588814807433782622324798520716104","0.5558690621223307634721158603176359462800825367786263051443129160328631551311961403577"]},"0,5,4":{"tau":[0.47722162225118536,0.5579177435088031],"highPrecisionTau":["0.4772216222511853600858835269279170748359638770973882762194610560980397222435732472427","0.5579177435088031317513698140040679174416197672779382340718145864876419620510972934248"]},"0,6,3":{"tau":[0.4002268745047575,0.5505193291289712],"highPrecisionTau":["0.400226874504757525678208186415780277659221503894824154505797929223420448510859259029","0.5505193291289711242101637705322417783379089512842106339708839947776748937096196802685"]},"0,7,2":{"tau":[0.32022070355933313,0.5326224962894985],"highPrecisionTau":["0.3202207035593331111754600935860983837288505790236501903824677154138703400711322442626","0.5326224962894985510170968657169428341472627844882339418116028707818250276314916765483"]},"0,8,1":{"tau":[0.2332316901885675,0.5011748764210534],"highPrecisionTau":["0.2332316901885674762306188725043922003629316132064274016113584946891834876912515886629","0.5011748764210534281316651678349490925535081635727841955527104148974785157704082629047"]},"0,9,0":{"tau":[0.12693194518398745,0.44723328196147555],"highPrecisionTau":["0.1269319451839874673464908159514107961815052180809058674960209823010786530691501634155","0.4472332819614755266932946235641348228192144149762099758844047448197555924034847402909"]},"1,0,8":{"tau":[0.9136737045620708,0.5545792734711206],"highPrecisionTau":["0.9136737045620708090302697229127568758004394669925635161537679997550649005312924402306","0.5545792734711205896967319504248535500487531465710751478477666849992396305070789808501"]},"1,1,7":{"tau":[0.7962008484973406,0.6030158039280128],"highPrecisionTau":["0.7962008484973406178721941447930601291602255808323300774585965539749463184125226623231","0.603015803928012846324211040985473753242084682657123462525188774150043598689918961324"]},"1,2,6":{"tau":[0.6967284862987426,0.6329664899109504],"highPrecisionTau":["0.6967284862987426641866793036719047699719434886639437814026095067528690999987261328901","0.6329664899109503432732136999176182388449811432040839603155515564190443427464516323853"]},"1,3,5":{"tau":[0.6041306363884229,0.6503436645774036],"highPrecisionTau":["0.6041306363884228557502006920331734344724603529353731706179607306800144795339118769892","0.6503436645774036081017912365743501542354299246502485333885710582198931470298102580479"]},"1,4,4":{"tau":[0.5141192369559593,0.6572773685404725],"highPrecisionTau":["0.5141192369559592792709500114077572264859329963900439733200169500581353076623811810229","0.6572773685404725191778986828762327668829139286897701984468574000011194892605643916367"]},"1,5,3":{"tau":[0.42400586746770635,0.6544141023426473],"highPrecisionTau":["0.4240058674677063535778671070524931226949066990464206820370368509730643787020103804959","0.6544141023426472490498000251083046156077314153838757226308415551784632887905331333022"]},"1,6,2":{"tau":[0.331101689698656,0.6413589740543262],"highPrecisionTau":["0.3311016896986560142877321597604885927890313077809931787475044536816301506431682881646","0.6413589740543261928013104360404355912574671049591463290180222742320911755078117563938"]},"1,7,1":{"tau":[0.23111429379408077,0.6164446751774848],"highPrecisionTau":["0.231114293794080774734709916435274105640496373845170191112186784511365692178417056024","0.6164446751774848557927995217538536054054945620681890807514058669178615449037385721614"]},"1,8,0":{"tau":[0.1127983490862128,0.5754377482669324],"highPrecisionTau":["0.1127983490862128026723207129704993286534318581960549906864635796149972458618263272367","0.5754377482669324021247511003045992902305942659582499942333007252203119200267442514942"]},"2,0,7":{"tau":[0.912397230479871,0.6697708427064374],"highPrecisionTau":["0.9123972304798710233551449392158518420623014388482547438593901095763382398682881544749","0.6697708427064373722550906970451528586060565160155593271002276014466950828663989680567"]},"2,1,6":{"tau":[0.7844552780175824,0.712614200135365],"highPrecisionTau":["0.7844552780175824559523677328156541431641699376798310741052114611927945425980838945711","0.7126142001353650784933042507499629538189504197944506071473225801889392567889974182011"]},"2,2,5":{"tau":[0.6725187345711167,0.7384467602617243],"highPrecisionTau":["0.672518734571116699509648248567620869590052697600849927689976492589033743122088277113","0.7384467602617242638895158750170916804449339309447690274974213314192192288412516703405"]},"2,3,4":{"tau":[0.5665518369109653,0.7517816841966645],"highPrecisionTau":["0.566551836910965376811215515499981989396160740102495819024262764287302550402437209547","0.7517816841966644987920304654669251654429645428328287775084946792219433562280960512687"]},"2,4,3":{"tau":[0.462099501317581,0.7540950926399174],"highPrecisionTau":["0.4620995013175809751365881589646015186172964949128080456782672694646345233863322063767","0.7540950926399174051575433734801414490575100550333929504127231535815505275945543894688"]},"2,5,2":{"tau":[0.35572927723620884,0.7455407913806367],"highPrecisionTau":["0.3557292772362088599462014360184847363117094202723312641418986815804674191840850848558","0.7455407913806366890729194851334020783531036014306584619029509221923032574648917307516"]},"2,6,1":{"tau":[0.2429065728888235,0.7251063477057245],"highPrecisionTau":["0.2429065728888234914520767490887756869121668542065150079358598021691238029292635037643","0.7251063477057244233653477535706407303515192708196040651153239964134703973409902566328"]},"2,7,0":{"tau":[0.11317042601141826,0.6898883177182566],"highPrecisionTau":["0.1131704260114182603773193399040641504065796279324918980588070242665956990051494346158","0.6898883177182566998889399435842244800485389539686358983124274398680990397603885024203"]},"3,0,6":{"tau":[0.9067124091344732,0.7826787101347856],"highPrecisionTau":["0.9067124091344731688970616042450424253422087573943836975189503163502097036248656387906","0.7826787101347856328659255790076509721280010198900544646213013365081847336184223401277"]},"3,1,5":{"tau":[0.7648628516050452,0.8233650732281796],"highPrecisionTau":["0.7648628516050451901649737417484720131456302364043250320530011234460666149009172628365","0.8233650732281796908786256721945353599698313939964480950998495671957493211159273081736"]},"3,2,4":{"tau":[0.6376872271367137,0.8462169074874232],"highPrecisionTau":["0.637687227136713692389101076206900374077254958710176828284816497222860470409771500164","0.8462169074874231820291405633032274111181907960952365772135115521883230929141456037989"]},"3,3,3":{"tau":[0.5152900929628823,0.8554069855901278],"highPrecisionTau":["0.5152900929628822836528581388475037732007921843235275932162599754167869923539633709309","0.8554069855901278287832422052870640178468365049630618050236835470987623251728795540705"]},"3,4,2":{"tau":[0.39255422302925636,0.8519925774317841],"highPrecisionTau":["0.3925542230292563710671258291958081693541369434624626137500580087337046524167407801448","0.8519925774317840853410334584823676540853431354663232449464487151427740029581124192859"]},"3,5,1":{"tau":[0.2642680884070581,0.8354776252924652],"highPrecisionTau":["0.2642680884070581451919624813383072907130963964892052924114005973018074609245077641383","0.8354776252924651517505784017489025830115114615024064362949949022959021922828390982745"]},"3,6,0":{"tau":[0.11994084773381948,0.8033902657430243],"highPrecisionTau":["0.1199408477338194822076144059766371824840521494910512711288580447690026245484636724763","0.8033902657430243138076330885258966277013214562154481744079041095114580665843360007086"]},"4,0,5":{"tau":[0.8971639809366251,0.9011761716805341],"highPrecisionTau":["0.8971639809366250672743042601927497449066134340298483546784207548116755588474590333828","0.9011761716805341006341387697503417771662709037822146382054344951347691544755590878896"]},"4,1,4":{"tau":[0.7363335146569225,0.9414993163681471],"highPrecisionTau":["0.7363335146569225896056095986897527868704162765378335035840805068846183501390171958811","0.9414993163681471045605509404793011004426765444060786528979429861753359784130129819063"]},"4,2,3":{"tau":[0.589046567364497,0.9614434107222732],"highPrecisionTau":["0.5890465673644970253339346805175838287103213932787464968388493550602444130469152502118","0.9614434107222732830256634631536037337458410261825001830071642686004350475668468823742"]},"4,3,2":{"tau":[0.4446503039921923,0.9652054073582298],"highPrecisionTau":["0.4446503039921922773753010806485444721337889831379920696748438854533250544046536426695","0.9652054073582298162383827405868421132685020117824362794574607050477241085062372695125"]},"4,4,1":{"tau":[0.2961510986030259,0.9532245778049533],"highPrecisionTau":["0.2961510986030258968323712103383015971453260966887442406457633320456087390199002912863","0.9532245778049533073255361540309515960503288153590936935387863246660665481194948722813"]},"4,5,0":{"tau":[0.1321565057912849,0.9232359840162704],"highPrecisionTau":["0.1321565057912849159063484985339535710937959438495906816609912167782658319337586904144","0.9232359840162703935026412858810281078165234233917570458407097405231346892818578265263"]},"5,0,4":{"tau":[0.8827780402701756,1.0322177292958412],"highPrecisionTau":["0.8827780402701756350777133174304177959051939955519362855959347478554030785119848518481","1.032217729295841307202897115612417457515307690342044670481650791943150937219642832706"]},"5,1,3":{"tau":[0.6951374123938552,1.073296201030138],"highPrecisionTau":["0.6951374123938551278050364276625081983985782982238312778866440698487107086006684426653","1.07329620103013794432353602281106860970730282782262277430450120688996975906743252629"]},"5,2,2":{"tau":[0.5194643875961141,1.0892017355913224],"highPrecisionTau":["0.519464387596114085212663421340600326707612922598126042320064738646790716376657901353","1.089201735591322339158645736454516754174905175669510787091707552668557524950068709975"]},"5,3,1":{"tau":[0.3427400369916056,1.0840791608266542],"highPrecisionTau":["0.3427400369916055953361965128442098826897607337436190218394253286071881262306928307647","1.084079160826654239184596381051429318681089404693444784370670120726591839633086721442"]},"5,4,0":{"tau":[0.1511330994806695,1.0562653963191864],"highPrecisionTau":["0.1511330994806695059453339788005540858805268067218644916558261871392963397063835393901","1.056265396319186446466899157028464104543469282505696853785902643266826837405147040975"]},"6,0,3":{"tau":[0.8609137234808509,1.1848760155011235],"highPrecisionTau":["0.8609137234808508804462395713111362996433954975936203226683278384115035837740848475941","1.184876015501123553264030995614024506893891229057355377403327855230051540120354729435"]},"6,1,2":{"tau":[0.6330312822893928,1.2272136268889113],"highPrecisionTau":["0.6330312822893927931063097990881707581122506709052680204369328790861868759015100731396","1.227213626888911297497891711506906898884561583949231655321987259356226450167081898486"]},"6,2,1":{"tau":[0.4134077837747312,1.2353204330404501],"highPrecisionTau":["0.4134077837747312055298928357946782724953315515742731611173272239220235684394359034402","1.235320433040450184674122879853663872085784308935157944300098695796774451355206071167"]},"6,3,0":{"tau":[0.18062458321973887,1.2115308098835826],"highPrecisionTau":["0.1806245832197388628016613301944947864290494186961581561406386285046303316419963568514","1.211530809883582598521283222319243504451398101136800193851900794565956289818332369982"]},"7,0,2":{"tau":[0.8252846473115993,1.3744349481415081],"highPrecisionTau":["0.8252846473115993000583257520761830566102880039392885966401551946888580915744941425337","1.374434948141508142007558931006539798467975867422874670639601394769408060263885480287"]},"7,1,1":{"tau":[0.5302474343356233,1.4167199794724403],"highPrecisionTau":["0.530247434335623291370265435945443505852681546543328215391861699974266702319666524666","1.416719979472440375979215389697480104012510454645450885047454196998656734733217679439"]},"7,2,0":{"tau":[0.22972834561982114,1.404035548910081],"highPrecisionTau":["0.2297283456198211473373909361349059383789173447292782244090900027264500810104787817743","1.404035548910081081333490725876124823816072384770694760343358163296941925058929077583"]},"8,0,1":{"tau":[0.7579388342521736,1.6332140920183909],"highPrecisionTau":["0.7579388342521736071410984610155705822045679910995072560583750497771945291323161058004","1.633214092018390830262229809864757523547625202964731672135995441047809910767320523534"]},"8,1,0":{"tau":[0.3246509413617829,1.6635454197877955],"highPrecisionTau":["0.3246509413617829319741417580776544701746355968294362578770034702582368764437007533204","1.663545419787795430809478641798195533619487705548094193136762971460344688364459721837"]},"9,0,0":{"tau":[0.5782813053022888,2.054602837847386],"highPrecisionTau":["0.5782813053022888666224578892797046972378670816647624564493195497899442557309325870068","2.054602837847386278392741945434590330852437132605779492967572371930483245461452510423"]},"0,0,10":{"tau":[0.9097487527704582,0.4111306827603229],"highPrecisionTau":["0.9097487527704582311829954035110662226462386637572822797406447692212361802232315758738","0.4111306827603229112076606665239937209981640300906934810190325052363006172169388998746"]},"0,1,9":{"tau":[0.8085205030604694,0.4731841000501837],"highPrecisionTau":["0.8085205030604694461325300212907446310711489163497119382313925288517186144290727204435","0.4731841000501836812793093361737933746741014432665465208931020638496840365821414556122"]},"0,2,8":{"tau":[0.727484913591236,0.5104695998726397],"highPrecisionTau":["0.7274849135912359935395427827278660431160290699415985677085502471519823022808833463531","0.5104695998726397307961867603762266663607305513313859679442355736331953669971344523805"]},"0,3,7":{"tau":[0.6537225575156654,0.5342070175203936],"highPrecisionTau":["0.6537225575156654015733107991827132789060879075169381315911251374713589078510508778357","0.5342070175203936190463111646672820642209455189936083210614496995384313437659833730767"]},"0,4,6":{"tau":[0.5832557827619747,0.5481280598591305],"highPrecisionTau":["0.5832557827619746525809213563212223583141684552349144256131657007918514102492268615786","0.548128059859130543893974458794355599356452913421144200562674016359558380508406435414"]},"0,5,5":{"tau":[0.5141321200052162,0.5537657245792036],"highPrecisionTau":["0.5141321200052162539130890920804087891171774889777019904260412534715358255103508080571","0.5537657245792036411719277629055726605103724618333836811983801444200981167317507660221"]},"0,6,4":{"tau":[0.4449961846412119,0.551605558803055],"highPrecisionTau":["0.4449961846412118680214941970304425110577902587156706714825489809457540713632857984239","0.5516055588030549979192325176288147562056255097815412790737486782235211234190840691204"]},"0,7,3":{"tau":[0.3745100535182204,0.5413557437284982],"highPrecisionTau":["0.3745100535182204283520743143631119892186922743516572439401310895923405283394766872266","0.5413557437284981652834254377242361187235346496078510938741787571912318084191642101086"]},"0,8,2":{"tau":[0.3007934566085034,0.5217792175341178],"highPrecisionTau":["0.3007934566085033701396994595313023709638607462420587085137884836541699047048467153618","0.5217792175341177899617030766835046022844365339677729222142071619310117329969797166622"]},"0,9,1":{"tau":[0.22010452946094078,0.4898137606886881],"highPrecisionTau":["0.2201045294609407835751005158402713186693271698348910471300566678598983735334896432693","0.4898137606886881105952215977268522842130944006575020736278451175042405721308286921784"]},"0,10,0":{"tau":[0.12058200288812793,0.4367980258823749],"highPrecisionTau":["0.120582002888127924364146423289292650954714253320446667603822165196076530125928182836","0.4367980258823748597960042388667750158583041418653129033955093339879957825050510735667"]},"1,0,9":{"tau":[0.9190384211101494,0.5375538951023546],"highPrecisionTau":["0.9190384211101494734414068157153184941830995069960151844715559069973313292736108285195","0.5375538951023546051355336994158317210659610343969470049771658910507894633127937016922"]},"1,1,8":{"tau":[0.8100412482120526,0.5844412215281697],"highPrecisionTau":["0.8100412482120525647923264364386024363933817660384897794440716628263333836182740328968","0.5844412215281697052887081226217038735876655266253464959166425205695879671939935788096"]},"1,2,7":{"tau":[0.7186835859387825,0.6145586210627302],"highPrecisionTau":["0.7186835859387825008478073229841623513711627998302965995382667425814616059843835336601","0.6145586210627301849283773345826778611031793580435261196197978751268513870922364041437"]},"1,3,6":{"tau":[0.6343628118808657,0.6334855626949449],"highPrecisionTau":["0.6343628118808657171002008285489444703199905391913219967473203546160192391994569092755","0.6334855626949449364723893762855386002688843998754264474544404896866235869189131962729"]},"1,4,5":{"tau":[0.5531333364818894,0.6433577571540114],"highPrecisionTau":["0.5531333364818894526888977167787875419595852753030661198563342677178280107223810770044","0.6433577571540114968711169143304900486208760167469335533524581659919433770517220488689"]},"1,5,4":{"tau":[0.4727572835445658,0.645028601098364],"highPrecisionTau":["0.4727572835445657897711160583205570614008491925779982896904681554112685775000780754485","0.6450286010983640168787039079495773279267674572764228743627301796170928524208020880567"]},"1,6,3":{"tau":[0.3913913054166512,0.6385772047483168],"highPrecisionTau":["0.3913913054166511952285889006601062657019995261046405197283409805910274780007300042907","0.638577204748316862200419338637901220761559077739381215657451965198834989381961572511"]},"1,7,2":{"tau":[0.3068027654273528,0.6233524874653934],"highPrecisionTau":["0.3068027654273527801354760197334907454378698255229963593347090798372196687991845151513","0.623352487465393433118506052669513049986195884241602511886819611303105885739359089508"]},"1,8,1":{"tau":[0.21505973900580513,0.5976350819356401],"highPrecisionTau":["0.2150597390058051379754469266165959426576872160433563899853342716296071882471572163058","0.5976350819356400397629158313504199220273680331012108034666905625106353839353682466911"]},"1,9,0":{"tau":[0.1055076289568248,0.557406650482888],"highPrecisionTau":["0.1055076289568248023481595975531730096860213820736295322227858294291885041810006943752","0.5574066504828880025758448129572641317054210702959268357067886038946615196839498946032"]},"2,0,8":{"tau":[0.9190606076711629,0.644560676162795],"highPrecisionTau":["0.9190606076711629274439796838945805462146440830850808157776831872000570389849940982253","0.6445606761627950106064310819017378035457942738271093642026916446636153861985829935983"]},"2,1,7":{"tau":[0.8020005852137946,0.6857419783728032],"highPrecisionTau":["0.8020005852137945376806602375315539194663056368284454984861325546528423225117853208063","0.6857419783728031833905837069669880707954548569491661238578471078660157550249000674937"]},"2,2,6":{"tau":[0.7006387670582919,0.7119871111940945],"highPrecisionTau":["0.7006387670582918332387796608590107182883483268733914559694400511330150745910301771306","0.7119871111940945587144826743549621709184614929093861202204839700971223586507505702888"]},"2,3,5":{"tau":[0.6056123636276993,0.7274684775197592],"highPrecisionTau":["0.6056123636276993634973917783410737630624409921247238639963484480079692939253212005269","0.7274684775197591397685519673784359392739974716684421876771786701687398925876687246496"]},"2,4,4":{"tau":[0.5130170313391375,0.7337170748264283],"highPrecisionTau":["0.5130170313391374209379406845900595514276440388821862875501239432627502587866916772142","0.7337170748264283347082754095583634994920907111821322568407822311753797436915978352951"]},"2,5,3":{"tau":[0.4202900175632987,0.7312134305812743],"highPrecisionTau":["0.4202900175632987047160961264515075506090430336500358024595057108820127965409537992383","0.7312134305812743003521169461276234399628337575368879709044542445185909295130946274858"]},"2,6,2":{"tau":[0.3248540459501834,0.7196891904709856],"highPrecisionTau":["0.3248540459501834309459843068053545938134694218663553733026307358606862067187822183142","0.7196891904709855471911517553488473634095238095420442251182240258796868392991067839124"]},"2,7,1":{"tau":[0.22273975015311986,0.6980096115053703],"highPrecisionTau":["0.2227397501531198651433313320988675161196668945109576118145108883241770219817143294766","0.6980096115053702959954238328694307944978066457308445010510781240203225356737588234232"]},"2,8,0":{"tau":[0.10424399065456832,0.6634698710594106],"highPrecisionTau":["0.1042439906545683291718902873849519542274907995590766198004319737674953111976726567331","0.6634698710594105606397609789980444569207715013460381069084018044272690332220807364215"]},"3,0,7":{"tau":[0.9152155571535485,0.7476036605229398],"highPrecisionTau":["0.915215557153548468126933152087089212048919177872991104561428040462611983880806250385","0.747603660522939788804215305558641288059563507552877357195140566637251186545585148456"]},"3,1,6":{"tau":[0.7875582513635886,0.7864254368915393],"highPrecisionTau":["0.7875582513635885825426314228104352165205668181688691612984989469762724327962179888028","0.7864254368915393128042702666085320044317383731117380156665638558676087627062384123145"]},"3,2,5":{"tau":[0.674390219373184,0.8101019886409881],"highPrecisionTau":["0.6743902193731839819010930133783194490048500157888232985715857536813346639692518209884","0.8101019886409880915423626376911826816118959019214989170082962195948229277093482708294"]},"3,3,4":{"tau":[0.5667527982850816,0.8224381075576594],"highPrecisionTau":["0.5667527982850816648839931090648323464432090014620580335562483638099708287623314671345","0.8224381075576593991495193455965995652255251101665397296639905228868299429728167606083"]},"3,4,3":{"tau":[0.4605137349559376,0.8246607204324721],"highPrecisionTau":["0.4605137349559376121837480416003691666454035518483808104097466460820842198677965307948","0.8246607204324720256693087191946098115664984552979044745227484783267919009620108539364"]},"3,5,2":{"tau":[0.3524306517251559,0.8169135208837309],"highPrecisionTau":["0.3524306517251559136326060199716910287116877149620732967528125040174791219301267793335","0.8169135208837308901605600787708190100444585703242345740811855986046857886212957903361"]},"3,6,1":{"tau":[0.23827210887070233,0.7984044507900027],"highPrecisionTau":["0.2382721088707023285657601886176418563283061581885434224176641925618732914189152037137","0.7984044507900026490289665757238428598578639553713252111365817273141494737996838400709"]},"3,7,0":{"tau":[0.10858981240248783,0.7668390401184514],"highPrecisionTau":["0.1085898124024878289095360127603484420903187604005999931922523136044288814151233560111","0.7668390401184514229588083289927586075974034227757814445188459228232170669839402104658"]},"4,0,6":{"tau":[0.9083686362080391,0.853471128910256],"highPrecisionTau":["0.9083686362080391457666791940200675722152313361562062249082895557972276924118284631071","0.853471128910255949936534148664034624515691297799938283900417061917484306378069706874"]},"4,1,5":{"tau":[0.7665877408192001,0.8916706712153726],"highPrecisionTau":["0.7665877408192001329905185106113827771956664157203096030944492303042890073970703898582","0.8916706712153725732472175471566519217860650701134100080232073449920759010124546147843"]},"4,2,4":{"tau":[0.6384606435127549,0.9132281048057832],"highPrecisionTau":["0.6384606435127549352801982705834649560645076611794193652814094787908874199039396926855","0.9132281048057831927041283092170043265766826670672219226464690133765713338128212167624"]},"4,3,3":{"tau":[0.5148113446942287,0.9219781522126247],"highPrecisionTau":["0.5148113446942287652480235115864151282683380867359402676207511631921227755297359205891","0.9219781522126247401418852981007161464767758836182221092137687001648563067180422423846"]},"4,4,2":{"tau":[0.3908091918996164,0.918884894011992],"highPrecisionTau":["0.3908091918996163836458659262274802317775619516588918214775814044704967815186201478188","0.9188848940119919800287743156866409467396296417159187139679171533834508250295825794801"]},"4,5,1":{"tau":[0.26152358969543366,0.9035272110690888],"highPrecisionTau":["0.2615235896954336429776629463742982070762532839152471091340575580685031369792294227192","0.9035272110690887969832850014441621721656072796709315364741178719496096650128120569616"]},"4,6,0":{"tau":[0.11716117880608924,0.8737082193732995],"highPrecisionTau":["0.1171611788060892451857549314140909689415255113781764288659893007546743732263561833823","0.873708219373299492559018781740758144569065771819293798520362242239111949563523937017"]},"5,0,5":{"tau":[0.8981825574879845,0.9674058441093298],"highPrecisionTau":["0.8981825574879844968008303508699159729493082402153648136797114476296636247072736240272","0.9674058441093298608671873605620390148584536349497528289507791184403153614642434777689"]},"5,1,4":{"tau":[0.7372529540960431,1.006118644389196],"highPrecisionTau":["0.7372529540960431118379173417386997017693667787705213431904437142440565016271391053633","1.006118644389196154714145333114457346482127526953547296292074479679154083568279297935"]},"5,2,3":{"tau":[0.5891746214659662,1.0253215576534365],"highPrecisionTau":["0.5891746214659662540269123127360659631125538044747944560314951278854532165731587993903","1.02532155765343641257042384978835300428377858234390881499639696812439624915913642686"]},"5,3,2":{"tau":[0.4438246155489609,1.0290425985006992],"highPrecisionTau":["0.4438246155489609235405631714983246615783681388064192360472679579110397061724478507396","1.029042598500699192827885142119258984667613643194400667502949400694164298975229725723"]},"5,4,1":{"tau":[0.29451315034278686,1.0177125624401901],"highPrecisionTau":["0.2945131503427868776148001721857462099366908671280518446110750507678732539068311898454","1.017712562440190132856550601915848999195062467463777237329074818318824177224000077518"]},"5,5,0":{"tau":[0.13037027880383903,0.989190627497009],"highPrecisionTau":["0.1303702788038390392375269262955746047929111977664795859567259612370579644788552600754","0.9891906274970090102436547282805374741391632654015093294324427640920943712698190506554"]},"6,0,4":{"tau":[0.883402440379279,1.0953251469013674],"highPrecisionTau":["0.8834024403792789959274272999066395263440085710054821217378132088959223602246557017956","1.095325146901367339503445097605014389740876032639153409556472191418512794380560134073"]},"6,1,3":{"tau":[0.6955470786137932,1.1353624094097583],"highPrecisionTau":["0.6955470786137932186484695754440777118446226735144491279950848159844484478751774071276","1.13536240940975824371812634597277727453954446849563882286301424200000001944233963053"]},"6,2,2":{"tau":[0.5192082355101879,1.1509252351605825],"highPrecisionTau":["0.5192082355101878586987562359564658629213431705675248341361258774251117101408530897706","1.150925235160582623979969560400804818889613891330726103926416332383358085565399627666"]},"6,3,1":{"tau":[0.34181112754341164,1.1460818676471523],"highPrecisionTau":["0.3418111275434116429031379723650755776158122146642467238990110130692397616301956528584","1.146081867647152239808267507603993618993477696656240909349007403047927359806213269666"]},"6,4,0":{"tau":[0.14996694760821253,1.1192104404698613],"highPrecisionTau":["0.1499669476082125125658653848884694861319420582282715304973718039290242583759344480003","1.119210440469861389612073201726076277602979073268243535774180966814943695342483581426"]},"7,0,3":{"tau":[0.8612740853110933,1.2458070785396733],"highPrecisionTau":["0.8612740853110933211804932394673768494157801333943319613603535929576795724269746790612","1.245807078539673288035104164996755037698840444366572842188697050903851575749902876496"]},"7,1,2":{"tau":[0.6331074021030786,1.2875138975974518],"highPrecisionTau":["0.6331074021030785710283419891071225812777509542159872519948644404427526703392432852283","1.287513897597451828712165164451478867456389610023731712742623509621866226395591879642"]},"7,2,1":{"tau":[0.4129566717215887,1.2955967017886751],"highPrecisionTau":["0.4129566717215886590107519143645261951590381832713242953429149710234336264142508331324","1.295596701788675059412625564017237143508295043317531832958960714583366077844275595403"]},"7,3,0":{"tau":[0.17988051687742626,1.2723691650393238],"highPrecisionTau":["0.1798805168774262682696285910572521078397276674812932367075854438862054030446352100603","1.272369165039323829413312405120135008250489204223627371718507850988920287205181474767"]},"8,0,2":{"tau":[0.8254589686702685,1.4338726835189464],"highPrecisionTau":["0.8254589686702684530890600983682051093715907895346290630049311433994735816854119749681","1.433872683518946311896827847023965149288501880569966830612317177863291293901961316042"]},"8,1,1":{"tau":[0.5301163985559172,1.4758532559781143],"highPrecisionTau":["0.5301163985559171730436915979856548030398756299364931963909465181080504029114021596286","1.47585325597811434419078961162265158574456721548040651872987703420597664938281028863"]},"8,2,0":{"tau":[0.22929014552917062,1.4634259603517479],"highPrecisionTau":["0.2292901455291706287898755210663644004394053033984805609364849069833568307067577726391","1.46342596035174785503170908040767280450708732751666609644746365378068628360584781516"]},"9,0,1":{"tau":[0.7579809867273274,1.6916950098860004],"highPrecisionTau":["0.7579809867273274340394245737183159245765860180673385166831415555930107829854021216083","1.691695009886000421107935827069443226085815593055247133634483862812335516931393463699"]},"9,1,0":{"tau":[0.32444264395471506,1.7220087308402454],"highPrecisionTau":["0.3244426439547150418335005243447246150096517966583535753764461516245644223022983886183","1.722008730840245388050020524575257807829901432857733440353060252871934701972630824612"]},"10,0,0":{"tau":[0.578241177463598,2.112585796401073],"highPrecisionTau":["0.5782411774635980251332596435413468642018407820538830672579350519913895905666842583884","2.112585796401072819440156958687628931747672980445556593406961587102019295346519233618"]},"0,0,11":{"tau":[0.9135705198512508,0.40274319084995835],"highPrecisionTau":["0.9135705198512508562467599637767572126303768773539714967606643564737534529199501803088","0.4027431908499583799739497243638518639527090875514927017746256062012694588057514798861"]},"0,1,10":{"tau":[0.8178748914036359,0.4633057764560679],"highPrecisionTau":["0.8178748914036359330092484774787552186013220338633505644636773410181721796082837790445","0.4633057764560679163945783701572315116570122859502497957648261008274761921686029932926"]},"0,2,9":{"tau":[0.7419995615345809,0.5003921913931301],"highPrecisionTau":["0.7419995615345809607337706927779890021035540272722804221086610593180162286931281295759","0.5003921913931301079879780165454001523106736869464091347533425313685401314329906588707"]},"0,3,8":{"tau":[0.6733854903201331,0.5248610171571945],"highPrecisionTau":["0.6733854903201330830647098645012613107325080345853528708303151699694364994244030372238","0.524861017157194505522759918294313434010238866869863828407647708961446291235610751469"]},"0,4,7":{"tau":[0.6082149762714262,0.5403819371013867],"highPrecisionTau":["0.6082149762714261919572433561154190972248505472803556382347788467162908053644981295088","0.5403819371013866664308420901061839104598229849475971034264565423782795442616787303288"]},"0,5,6":{"tau":[0.5446816420698254,0.5485497483853404],"highPrecisionTau":["0.5446816420698254002648353292563199420447806096737747389266421237310057593862858126666","0.548549748385340430887979951425264355742901157502387647892862556943679236353527525421"]},"0,6,5":{"tau":[0.4816333320402924,0.5500281308875302],"highPrecisionTau":["0.4816333320402923783139300154037572748650955990957191979509067933286525864316144463136","0.5500281308875302720116404501786520424173878835068226151855829525676565878821136979835"]},"0,7,4":{"tau":[0.4180892754139823,0.544884048450245],"highPrecisionTau":["0.4180892754139822880432662587876217148713962526809297802163036529310828775384026233337","0.544884048450244936957198173870475992132206981917353671039463452967242552613265747856"]},"0,8,3":{"tau":[0.3529239050607367,0.5326145495488299],"highPrecisionTau":["0.3529239050607366655454066099317006323253306972639389023586768412831034236155953611566","0.5326145495488299203852863647359907738304041120920926111598882625057846486571514217447"]},"0,9,2":{"tau":[0.2844107691180769,0.5118959676022632],"highPrecisionTau":["0.2844107691180768892572073966021934596492498012139616150078080480046475700445396314516","0.5118959676022632057725436009030211592419540009734332206161802664685760032039710926685"]},"0,10,1":{"tau":[0.2089795013447162,0.4796957634301651],"highPrecisionTau":["0.2089795013447162036851714954156573216162778675977840078032491571880934273857624564103","0.4796957634301650971968214709792949534522231574215224352623510794378931038034110545238"]},"0,11,0":{"tau":[0.11516327836048905,0.42761431079597856],"highPrecisionTau":["0.1151632783604890533432331821087081478349047040535333043184212405714026241269235430836","0.4276143107959785683501567588329211986750411941954409290434918200923864812322377450748"]},"1,0,10":{"tau":[0.923522083580442,0.5228436920813326],"highPrecisionTau":["0.9235220835804419729577901991715419987250549366151670161340742648752549472330305850573","0.5228436920813325850993154657244999545263440604383660804820052713429442048377723484474"]},"1,1,9":{"tau":[0.8215231642780433,0.5683143809096617],"highPrecisionTau":["0.8215231642780433163580964151508197060617665715690054629989497703111854809781069165924","0.5683143809096617399450528390543231460022710504087478773373460556851723344594939528437"]},"1,2,8":{"tau":[0.7367880463090198,0.59836942505751],"highPrecisionTau":["0.7367880463090198132721744920254927279897686533349195982415736543507545514329209996966","0.5983694250575099878466405617122802907415133315223020647848599322052468157456215326461"]},"1,3,7":{"tau":[0.6591386056231198,0.6182916746489376],"highPrecisionTau":["0.6591386056231197271435140457983655938322555943486422364954586633700828542286786548335","0.6182916746489376934191796386591662131170555386956679966406374860962732272471142317643"]},"1,4,6":{"tau":[0.5848640159391224,0.6301805299640165],"highPrecisionTau":["0.5848640159391224500202942814738757445345565111037668203820102124086932036139173125612","0.6301805299640165664700033413926517568730686010027946935039420256944468510877846810768"]},"1,5,5":{"tau":[0.5119820863111095,0.6349826256776173],"highPrecisionTau":["0.511982086311109493507537708981070304293983501675918253805120518826341679194922287459","0.634982625677617368933222150397826033005457464953172667965156308115136494190607712075"]},"1,6,4":{"tau":[0.43905263991071825,0.6330082494050312],"highPrecisionTau":["0.4390526399107182211685785980424288875456818676329292572606016986009538480026142471781","0.6330082494050312081414004218894735169002607728985982453499232710112852169840014532377"]},"1,7,3":{"tau":[0.36463837528057746,0.6240648804722684],"highPrecisionTau":["0.3646383752805774757056101600163852569046239707737847813070260240295458678045085279706","0.6240648804722683448542266968189657173543669460794700737676483979638659290423613522375"]},"1,8,2":{"tau":[0.28676838049023473,0.6073809424113679],"highPrecisionTau":["0.2867683804902347196048770132181734408377724815427078859471807379372273847519093888711","0.607380942411367894129323381781689047775785734647530377864089215402518663188723312206"]},"1,9,1":{"tau":[0.2017602409055539,0.5812406655903369],"highPrecisionTau":["0.2017602409055539024237780505374384754377283992962996719077680066461514799036460578812","0.5812406655903369127449874816711285234110524492044475708270495517964386733075509926213"]},"1,10,0":{"tau":[0.09943620012791543,0.5418265137711081],"highPrecisionTau":["0.09943620012791542854090753386088448607309338262289711771613578994760203294436265508793","0.5418265137711081193891078239328210091146579330393318042662484261089917850757301491147"]},"2,0,9":{"tau":[0.9244965945055537,0.6231238419640506],"highPrecisionTau":["0.9244965945055537033866823112469756256011078600394006562629326277668886820519369298468","0.6231238419640505857391914396910870295710390379396611990125659594416028634580588418946"]},"2,1,8":{"tau":[0.816210312738608,0.6628008359309334],"highPrecisionTau":["0.8162103127386079787121032444016782168771948896389734779048405187939297714940239719734","0.6628008359309334041428724950659488356060212683078733243359898863105974991188074412735"]},"2,2,7":{"tau":[0.7232639602197353,0.6891250629829302],"highPrecisionTau":["0.7232639602197352254594581415315386831025598502118060684430781336504283050946623488299","0.6891250629829301832698850544321089553722294994043497262523444491149964340155231611045"]},"2,3,6":{"tau":[0.6368087234473495,0.7059664572287357],"highPrecisionTau":["0.6368087234473494638359635669485754181674330806923235345901781168427281678922737541681","0.7059664572287357333740441603834797216648977568595418084074630706603682411036702325727"]},"2,4,5":{"tau":[0.5532820844409287,0.7148409867471551],"highPrecisionTau":["0.5532820844409287833357437879776394838444857435115188940921776552241306533686179426903","0.714840986747155085843323734383676438236257483028694181966811019171667380418960605625"]},"2,5,4":{"tau":[0.4705583079566005,0.7163735393510486],"highPrecisionTau":["0.470558307956600474507400414620380593782529567954615469077282222716876704060368922447","0.7163735393510486857681860643708441185343843942072203464329650772129084309593403998324"]},"2,6,3":{"tau":[0.38684796399943494,0.7106298176971176],"highPrecisionTau":["0.3868479639994349430225714719899321473338348149817255007728172931900630467559838116734","0.7106298176971176784786254972163817223954278795314535989028442801686490149947490554444"]},"2,7,2":{"tau":[0.300010590519598,0.6971539737455481],"highPrecisionTau":["0.300010590519597999647105645893867578735816833285816684458976981468653224201451326359","0.6971539737455481075217255949168136349537271229129299777667027749599942349319687208537"]},"2,8,1":{"tau":[0.2064339173398105,0.6747839892977446],"highPrecisionTau":["0.2064339173398105128226331189422485377291326658964226835679961980105183182403894913696","0.6747839892977446495462629976967210131987897901569154818185291436122593194174486047625"]},"2,9,0":{"tau":[0.09699457093004263,0.6410066821773549],"highPrecisionTau":["0.09699457093004263119384513204303825470175925112677859259127493842740811552845596814454","0.6410066821773549290063426871108965278667776990600113582278140323452289289917342030057"]},"3,0,8":{"tau":[0.9219625862618696,0.718300170447209],"highPrecisionTau":["0.9219625862618696061758688595737421905240876791857212140805567502537527809082236774776","0.7183001704472089347864251584625502352725791620628955207762600282659421751650401793775"]},"3,1,7":{"tau":[0.8054263691316562,0.7554462829596336],"highPrecisionTau":["0.8054263691316562442319122204624785224895826152683646303563315767336616436374168645858","0.7554462829596335624871084168398299413593267332432904986218220490608763299456931131483"]},"3,2,6":{"tau":[0.7030651717204842,0.7794188051572261],"highPrecisionTau":["0.7030651717204842244282132142943324863367355149652892781798256180198037492291777877131","0.7794188051572260179160896286115659931530288578172995012257030027097997007607704905724"]},"3,3,5":{"tau":[0.6065802114192296,0.793688038869172],"highPrecisionTau":["0.6065802114192295741928984953593369887253824383555008067435292972915409424512204255707","0.7936880388691720546704724054798711150377555932522159322470428498729507578499856056018"]},"3,4,4":{"tau":[0.5123726481422266,0.7995062578759802],"highPrecisionTau":["0.5123726481422265599740854126566306903558057254710494262822308803553720402195893492211","0.7995062578759801552128594227458599781736672533371420798185933127730759762653207654099"]},"3,5,3":{"tau":[0.4180178046003677,0.7972728654869734],"highPrecisionTau":["0.4180178046003676656672558986699584577526784127736206810543667586321849056077347280767","0.7972728654869733647550574697703155877700121476831940458820577605284989641581452164532"]},"3,6,2":{"tau":[0.3210707964935976,0.7867852240665781],"highPrecisionTau":["0.3210707964935976112350320523427308997437626026791105650375553431048167081070466800616","0.7867852240665781772616183315830587848102322068013026533684236450886326865187558994586"]},"3,7,1":{"tau":[0.21784528541401077,0.7671610731254896],"highPrecisionTau":["0.217845285414010761513443469105781761603094063108797437106675305710648210187395728737","0.7671610731254896729605482785285286417431246238828229964875914080279358042378041904809"]},"3,8,0":{"tau":[0.09963559115192062,0.7363017420446022],"highPrecisionTau":["0.09963559115192061956899548107391425144481966370200223286897439919284118394594493401995","0.7363017420446021480470998456497556090176589653007422197580369205508593401627484394277"]},"4,0,7":{"tau":[0.916960169040881,0.8144531131952817],"highPrecisionTau":["0.9169601690408809186246636828724219343185859652237106020939546142606121149980169840437","0.8144531131952817299749420385298760925556821143594853871666711333635953177386089857408"]},"4,1,6":{"tau":[0.7895758107902655,0.8507391923521915],"highPrecisionTau":["0.7895758107902655512486732280829055876718546644080402701156934234957882013813389147809","0.8507391923521914810130231776890859870621086206062339584841885362154923116651306697272"]},"4,2,5":{"tau":[0.6756412810453973,0.8729810884061157],"highPrecisionTau":["0.6756412810453973904055663821081764168174224819735103757989300866233262020124105056617","0.8729810884061156195313865115240512359438503425280966671173991774064624538167007606631"]},"4,3,4":{"tau":[0.5668984611129122,0.884639177811809],"highPrecisionTau":["0.5668984611129121994589108955069267205636449388407635843461345256383200046194607289302","0.8846391778118089982934913390838681225022673136458098671405745611633748223500631312644"]},"4,4,3":{"tau":[0.45946160702678335,0.8868103076086719],"highPrecisionTau":["0.4594616070267833741943909767192444741196871359463386593822496857776718624092806038259","0.8868103076086718204060716556428265437134748332249803252414822823836682594667293516968"]},"4,5,2":{"tau":[0.35024985048842067,0.8796329618867217],"highPrecisionTau":["0.3502498504884206435063172396714807211849384803266306496700798520331481697482089513944","0.8796329618867216644413114364845278139381497450162649475601183453383681095935270726631"]},"4,6,1":{"tau":[0.23526991624878854,0.8624294725194341],"highPrecisionTau":["0.2352699162487885284366729164465139443476285328176830084648017389138659937554696848992","0.8624294725194341708860914159579891683441220997592890794821697286906786086199163173775"]},"4,7,0":{"tau":[0.10575014763779851,0.8331886558415479],"highPrecisionTau":["0.105750147637798504363794784310411560877099474467030513040625825984552168658528493866","0.8331886558415479054521844678553935594189932959953897031016212910555971498483236217661"]},"5,0,6":{"tau":[0.909478336130215,0.9158150216727714],"highPrecisionTau":["0.9094783361302150557108391902872134878151302511125261183905395470129389570610486454509","0.9158150216727714061180676088322787470162951439676538495345755813981965715875651490709"]},"5,1,5":{"tau":[0.7677772064188196,0.952327367648329],"highPrecisionTau":["0.7677772064188196649284631197180791841321858366345961937191071888776850123902300524408","0.9523273676483290353445104806012995746565741833906869227211291770674063497632723065479"]},"5,2,4":{"tau":[0.6390028604706922,0.9729788345009572],"highPrecisionTau":["0.6390028604706921686131628287473855067271967539930432855073123454006013604195334238048","0.9729788345009571892015887753426938460860487074758568881528461516129652308485569437618"]},"5,3,3":{"tau":[0.5144782751688971,0.9814199027712702],"highPrecisionTau":["0.5144782751688970918713312926290899896154236098105216374807884654669169209476358732632","0.9814199027712701980263237095685257542441618212941891126421480657187929375164556881108"]},"5,4,2":{"tau":[0.38959268539149444,0.9785620709554678],"highPrecisionTau":["0.389592685391494444626287364037887056147695881806248057170551690733193720477209607699","0.9785620709554677761187422520986295121086858261294875040506384487197412080870911810467"]},"5,5,1":{"tau":[0.25963277795822404,0.9640250154428804],"highPrecisionTau":["0.2596327779582240347786605170179396712824517430939205156014848372751262484317789116044","0.9640250154428804291726894752281089905547195297946662794988161739487143233126190082015"]},"5,6,0":{"tau":[0.11529505721141707,0.935752735570142],"highPrecisionTau":["0.1152950572114170703893530136458360072725746142689572847685918920209894903142273661161","0.935752735570142029910630096606450211988774383154740880639740882367062441333513030863"]},"6,0,5":{"tau":[0.898900147925036,1.0266196244194499],"highPrecisionTau":["0.8989001479250360826020908296251193327980945899320215236940797700370075731198803956087","1.026619624419449899129963787214252374200126543433037019167616344311070475509497444695"]},"6,1,4":{"tau":[0.7379121304054201,1.0641875888152226],"highPrecisionTau":["0.7379121304054201072914374133142229745104293681082260911944907498704607030831524366738","1.064187588815222678016022646216693715098806258784928697832668016710333267443837216241"]},"6,2,3":{"tau":[0.5892673990558824,1.0828543465832914],"highPrecisionTau":["0.5892673990558824334081194886208872474692048165176690789166317775794190349145836274304","1.082854346583291306319814091792350454094987147847469199206105812456279483213549876207"]},"6,3,2":{"tau":[0.44322851898129456,1.0865484268975043],"highPrecisionTau":["0.4432285189812945556654389124176920227244013331260448749302817511867544285587546095924","1.086548426897504291271466771835886332481478380282823432225603503940668201792125567091"]},"6,4,1":{"tau":[0.29333838206637347,1.0756945523849286],"highPrecisionTau":["0.2933383820663734876143170777222620996630025901507038580475372079530507079181399650256","1.075694552384928558999578233253509634807386921547107769770771502313579670247796692677"]},"6,5,0":{"tau":[0.12910937368416608,1.0482188756958317],"highPrecisionTau":["0.1291093736841660978415936903634978882903302453792273876973833400043405925700048757768","1.048218875695831706997893658057598133796981413306640672636755721399383301736979324091"]},"7,0,4":{"tau":[0.8838573865685705,1.1523009263231543],"highPrecisionTau":["0.8838573865685704500944383792515242195836393120995967023832957374767500028505003361808","1.152300926323154222214222354062904767781832200894453505682903058107558516692613909228"]},"7,1,3":{"tau":[0.6958486702907646,1.1915745597379706],"highPrecisionTau":["0.6958486702907645589399707376838871480383138477220002282574834491371340761573506991421","1.191574559737970556016219188405528746475727333367392837091317087432202191982400679437"]},"7,2,2":{"tau":[0.5190187547853987,1.206884256897688],"highPrecisionTau":["0.5190187547853986731772602043277033549724554993818830449918504254568060454985503344797","1.206884256897687982650573772914494550978849532151797469442396240524396622418158710142"]},"7,3,1":{"tau":[0.3411264129888459,1.2022500080874907],"highPrecisionTau":["0.3411264129888458942294885644772922262796535119242351116675528934567107442822996262415","1.202250008087490795494728593931030177706114438407810760330002294102075733932591626784"]},"7,4,0":{"tau":[0.14911578112600785,1.1760716818208807],"highPrecisionTau":["0.1491157811260078541039496021823864906414465677516271271250007556208241244717948208489","1.176071681820880760877765400286112572008613113636319176406919723986683905479600152155"]},"8,0,3":{"tau":[0.8615427773268887,1.3011769323276752],"highPrecisionTau":["0.8615427773268886139879252450336681552126440034252649047498804971539866609268446597465","1.301176932327675126633026282892309352708089881706277570105579063232174367412963012843"]},"8,1,2":{"tau":[0.6331643609690645,1.3424113761332421],"highPrecisionTau":["0.6331643609690645133728193714605564506875411490191975167557432316176208189389026460193","1.34241137613324217470805706768958059916785599283596215216046255481264973233131841837"]},"8,2,1":{"tau":[0.4126178490351031,1.350476962940132],"highPrecisionTau":["0.4126178490351031033567431401002192881778207226364897695081527823189070343974142380196","1.350476962940132143661314446148842438090532344296801321282229274506424738406548009729"]},"8,3,0":{"tau":[0.17932493075464095,1.3276715574532627],"highPrecisionTau":["0.1793249307546409464908047704803337719381577931882926663647123159969628568842405782946","1.327671557453262642868664507926703143137033164427689571613915762144484120737301505696"]},"9,0,2":{"tau":[0.8255909930696255,1.4881190103136215],"highPrecisionTau":["0.8255909930696254385697971254395621258546023922022563765842967135839609326471939993589","1.488119010313621415543923385067656097928792988852525057929594667330087782494390661077"]},"9,1,1":{"tau":[0.5300167236508954,1.5298685448852112],"highPrecisionTau":["0.5300167236508953312960493092019995772294675595928494076986599394220193201498810906052","1.529868544885211157035627998543097612071083948731358961551148166444533027239211506622"]},"9,2,0":{"tau":[0.22895789502527408,1.5176370820909544],"highPrecisionTau":["0.2289578950252740884004843355837457432670495654230058020510583388769893642561626099664","1.517637082090954479274640733603402303340614603532317153326775269385697102043976096036"]},"10,0,1":{"tau":[0.7580132224444402,1.745212327495684],"highPrecisionTau":["0.7580132224444402370842221420850855951868034863207474220892669411656637373529934248402","1.745212327495683822505253994103933648187178247744828935507184861821561378726431900062"]},"10,1,0":{"tau":[0.32428317208656704,1.7755127431036775],"highPrecisionTau":["0.3242831720865670354123909088632218890267988839170007772557912146083859335600299058348","1.775512743103677402323272501950754971075152678273899048412310728558364569289836384503"]},"11,0,0":{"tau":[0.578210300862536,2.1657209292992188],"highPrecisionTau":["0.5782103008625360088420900191486876606838751668915027312211958223157391736243088166203","2.165720929299218616175686199131973625932322574539811368854081776190694379553616923172"]},"0,2,10":{"tau":[0.7544364288904193,0.4912609259848916],"highPrecisionTau":["0.75443642889041929142891105548579757470698568032001150861046500890344677212727603399302911729039965449014677389229913319","0.491260925984891613021666299931098658702281301665198687940900502883524008907790259896172972306810647137433407892328510508"]},"0,0,12":{"tau":[0.9168889024623954,0.39528598502874823],"highPrecisionTau":["0.916888902462395446063926092199181997804197190656871981466553789129342409087844621831777983137138519184532878678243621449","0.395285985028748252478110848995115375646303831272396180422395290918030346553297417617120194564763468741634104605127874768"]},"0,5,7":{"tau":[0.5704658865227374,0.5428437214463359],"highPrecisionTau":["0.57046588652273741229899160813508555296033848929495176546385864227115492891200934927497244695962788653629785567538641761","0.542843721446335809224624071001665006102825887743030625881408542006605768556528706862513249485445321154770274723051984479"]},"0,1,11":{"tau":[0.8259330855733517,0.45447266431652494],"highPrecisionTau":["0.825933085573351677391226445872271388334203840414105455627885860762864773194810359326995371663989661236097825407897101833","0.454472664316524918852150070413115618111037105175483193400755488592433604461088822786898155022840522465812149901723861895"]},"0,6,6":{"tau":[0.512320365083433,0.5469512669334818],"highPrecisionTau":["0.512320365083433030930835806782361205350810515413491228396463038304883972952725659179863027507220234534278330098268946271","0.546951266933481862570788465514290864068879470830681619724175370211125254546555092948608827882372259348927802445398929515"]},"0,4,8":{"tau":[0.6294008579845138,0.5328404533379976],"highPrecisionTau":["0.629400857984513834204721403186844221352872426006099067681273181193094448601283130233319615934332098027605844265648245062","0.532840453337997550275045780621288562889609859448673784169953060775424939273930062403842957803480369106548296060068330845"]},"0,7,5":{"tau":[0.45417140318787264,0.5454176752006793],"highPrecisionTau":["0.454171403187872625913849253740235194617508094075441721347912240236401107174828760781323310858574553833474217056416766745","0.545417675200679270860472455345847350668915275172164649761392804594818286200439717315760002274890845801739281502163787602"]},"0,3,9":{"tau":[0.6901557606245385,0.5161845803933376],"highPrecisionTau":["0.690155760624538538509340437370589151671486215458656024161572009394684598071782423519598890585567785393803052664726787893","0.516184580393337538171567183932246686956880219357848215667173374031645689732111963143712409139746289714564843696563572248"]},"0,8,4":{"tau":[0.3952346169131572,0.5380927070925525],"highPrecisionTau":["0.395234616913157203954289025823832390152530969089237389881009652138426458050826463045903614930933392123422008526845192988","0.538092707092552507375828653514589136236208823222576270630924679429140941224958567802597531572231067061480892793919996272"]},"0,12,0":{"tau":[0.11047350993933404,0.41944521458261735],"highPrecisionTau":["0.110473509939334041922415246835750143872916595810967332297653932039938410571577228763625199822889669980888885346795091939","0.419445214582617326245481835612061937234348650917142279412642932765469778740256762737575292261288346189363510154110992799"]},"0,9,3":{"tau":[0.33451156183459035,0.5243611982095907],"highPrecisionTau":["0.33451156183459032697597954438537656853250371414296827368255790739695165616724921944830525451397843430896870286701619581","0.524361198209590682948095132021897773401973489569816044054863522557154983634178471225191925795320550859869103290333306374"]},"1,2,9":{"tau":[0.7520036997418469,0.5839968907710933],"highPrecisionTau":["0.752003699741846859609520941431511099164694414315214042206386726907702709250490183576002168736166730952830230549074957012","0.583996890771093296701115510083401708091678398330194015018561965979046341792388858827717405893316531544770599896359716911"]},"0,10,2":{"tau":[0.2703817336308169,0.5028622115694782],"highPrecisionTau":["0.270381733630816922129097859039113807427309926125243907992531113509946588020788025072669229719295381637884756389071547298","0.502862211569478159697823188543607592494827011077093516638291792702963754247369909860337128151769663354618462051989198545"]},"1,1,10":{"tau":[0.8312223409603716,0.5541397463511615],"highPrecisionTau":["0.831222340960371576098081734634656238853936257102157771671809313647815448316141769495335228769332049919223170224450286513","0.554139746351161428138542263323443759099765765809166606337106641760717904122413250420983557340025982135687472573172221687"]},"0,11,1":{"tau":[0.19941082102566102,0.4706103459943639],"highPrecisionTau":["0.199410821025661025333263038808868698188732580980114845738866484136237531291546781349224920876870231916740054444407033757","0.470610345994363928246506858841717433758240462121480257187060388034038679269722222812151920324151772464212546103271570709"]},"1,0,11":{"tau":[0.9273340917788601,0.509963083239301],"highPrecisionTau":["0.92733409177886001039408106699869968672842412480739447308140352004203524987700995625240515663399566871067780639671265613","0.509963083239300901650641330768661935134746563566651690388183242215058302663984858239185132906576473059606507148817445377"]},"1,5,6":{"tau":[0.5443523900104076,0.6249012215109597],"highPrecisionTau":["0.544352390010407594838028008468773244292959389739415046817292475587988962763431390875808694945554700658607666378917971514","0.62490122151095976171071728533514639196993094130320583626600622238869343500999468699189853693586010552206059586063501002"]},"1,3,8":{"tau":[0.6798571825047544,0.6045554644681472],"highPrecisionTau":["0.67985718250475438609409726747680176386163958312668243189994832005882897172035931018223380913163486327304054139878604694","0.604555464468147235645851311583561323301622791475613771037447619230034162745749510870787500973308931600208851080189442993"]},"1,4,7":{"tau":[0.6112467139478244,0.6178617993077462],"highPrecisionTau":["0.611246713947824362387882066415744945899627872390616723274201626661267145192249920054329485615782658928772754805176040642","0.617861799307746235992611943859989350977148658223812071310706738237198762034448278419145636260742760226565514455659083068"]},"1,8,3":{"tau":[0.34225549746841977,0.61079019952884],"highPrecisionTau":["0.34225549746841977993627069497038559870863646023934290370956252961887310468089818681856048211338683405594200185909971918","0.610790199528839956798365997154797642189549492792937389455030940224139566424463894635598992321585745404494861404975734175"]},"1,9,2":{"tau":[0.26993649978495465,0.5931116851134267],"highPrecisionTau":["0.269936499784954655022901280373828426974792513924535985060834941391688800600730052633128596902890671206317487052654571537","0.593111685113426677363534094399795889205081306365230976281361201972648963451174879802568833773935787065683679648933114432"]},"1,6,5":{"tau":[0.47795527594499804,0.6261042808102064],"highPrecisionTau":["0.477955275944998026972479036863439683096411571907267358042141781726829166756417625785095109364420439797618246767010288504","0.626104280810206453434146641374291075191342769050337528103130158702682557051492029688635384616852068709764940406716883877"]},"1,7,4":{"tau":[0.4109939840629066,0.62151187654227],"highPrecisionTau":["0.410993984062906585838391834314468249070047330755734494982195012246190859979400877824569705758941247411367370133043947534","0.621511876542270003571646102501431700204806249069679234514082081523987638898735641742772099171967156174311316187669939491"]},"1,10,1":{"tau":[0.19054122518628272,0.5667886872702113],"highPrecisionTau":["0.190541225186282730230808873793639749609154156884533852873755800163035268027951107650487243623810091731054084169617699815","0.566788687270211273960782153120996798295531796417669334287353394817163309513658672517957826638742847904426381996748485399"]},"1,11,0":{"tau":[0.09429044867058786,0.5281854129680246],"highPrecisionTau":["0.0942904486705878621023620030705736814449025186908530159879670608396917656317856307212667857447708497937748144099319897484","0.528185412968024515564982538609874598303893568242217392028236852255570522444372611496789960783725177102547803512594171751"]},"2,0,10":{"tau":[0.9290256507160454,0.604603625008783],"highPrecisionTau":["0.929025650716045385801213223428921121797895439043006649022482310046161562695191056834384207211476473522701671257432209866","0.604603625008782915312114714091263683466566651715321772005953813677454950706263267585958120737762523517296389308947757367"]},"2,2,8":{"tau":[0.7418982880000281,0.6691329974537845],"highPrecisionTau":["0.741898288000028116998184406425036179577403246680070006332389780697411149748172476047366381155907151999306727003913847119","0.66913299745378452922982326945675995743785384533041425756592493368766285261501256073192649639377301664894037836203456649"]},"2,4,6":{"tau":[0.5860237643127753,0.6975285694573233],"highPrecisionTau":["0.586023764312775286900502171964785904777091849146883099497008252302565899039309777539103413180155283961130337018097914335","0.697528569457323315477341238934393725084552135501255370477334787346881903909799159504137455542781384001152431410826385396"]},"2,3,7":{"tau":[0.6623561153474223,0.6868441923852341],"highPrecisionTau":["0.662356115347422359456241514652929873999223374980368280653689702185873274157971692381920244848880826152195135879229459165","0.686844192385234125885172883524651335390615241532607252879823963695439056497608942548540624049123610036873014883331223091"]},"2,6,4":{"tau":[0.4359689730349024,0.7001131108987093],"highPrecisionTau":["0.435968973034902405093395741593304198182222257400111706554479892554489298668431652159180638756558308950237184233074416134","0.700113110898709268845106064768285241770633532273441302211398542420280116777288193073118339197831276821954866389440779797"]},"2,5,5":{"tau":[0.5110288526677386,0.7018741660461616],"highPrecisionTau":["0.511028852667738532943937466267381303889935028938728976256864618953295007339082563184567156846767915950402677370382665246","0.701874166046161619303876907780140276958262886687599463500704098486793563884526873705092191631427435875671596581018235449"]},"2,1,9":{"tau":[0.8279765767054595,0.6429215888056666],"highPrecisionTau":["0.827976576705459552352308683589321917184052696341590586172064041415958809676225082909426926930219673927088217779362985393","0.642921588805666595138554704602605193018078356875310737795434834074547295398266663104139234039189415073910474381277311344"]},"2,9,1":{"tau":[0.19295317393448094,0.6545999890300367],"highPrecisionTau":["0.192953173934480929000968349763522998627140134850568176993599218082527016563633182222709560022832349471213845276673285388","0.654599989030036709676161683799940387714553581388588864912177682157744220420111343802265277844971347070429983178583011309"]},"2,7,3":{"tau":[0.35943873606023946,0.6921100903969195],"highPrecisionTau":["0.35943873606023945455396466425373504195807340051652616710812156471130845101082299972903671944585695524249370518381191556","0.692110090396919473642765079548124229460711611117974445921077032856536090275945553080173192170076541235977762903868756958"]},"2,10,0":{"tau":[0.0909777157670327,0.621603861407134],"highPrecisionTau":["0.0909777157670327042893388613874265610966678666204211485429169537160701666730605737419439486607592597658249766027426189677","0.621603861407134014910395463410433723554728780374074066671677842264023500650623539760365981166910947557245989620939116152"]},"3,0,9":{"tau":[0.9274586254026601,0.6933458420450004],"highPrecisionTau":["0.927458625402660158631096772200587833784048522647408050779361152881577456240530784443051983237948045999509801442043476243","0.693345842045000414078754179653711634929350506268647226335448939133837377543390251600639384367546977882072903354056556276"]},"2,8,2":{"tau":[0.2795543429862328,0.677322203402726],"highPrecisionTau":["0.27955434298623275078594579397644830170569535390668937397628135348921027561522316580177175845871952098171438473953237857","0.67732220340272593157464757518860714846270741840875799077620070252725497906019426223522763196388397650109531085319795928"]},"3,2,7":{"tau":[0.7261349409172455,0.752970560521211],"highPrecisionTau":["0.726134940917245525048074991276167503782860681744376795903725047570232665632649261762813446165574643389319734975317793629","0.752970560521210998706134257556853554819309109955631808121745659444949364526247495792523129906571241630630533922170800275"]},"3,1,8":{"tau":[0.8198886845084018,0.728993442826741],"highPrecisionTau":["0.81988868450840182260712571029124423958049219282206777223019756531654259172293528795633210766477715662518264957773667118","0.728993442826740929219195743273552079751320683629192169280419628312008102563456194575321753754972943476123322454372512147"]},"3,3,6":{"tau":[0.6384023022751889,0.7684570129377951],"highPrecisionTau":["0.638402302275188884599932186236406818484770934206827700483769246125791958044902974403046668371682272911831488080255748233","0.768457012937795152008852167053516598477623247153473999006132153228282832808067398294565762811025341166379961951203518152"]},"3,6,3":{"tau":[0.38401793848871507,0.7729059826121494],"highPrecisionTau":["0.384017938488715084360308166560970991486076752322194712772198047842505513481771552932581207668665289039691703146155370825","0.772905982612149394800266123589751940824766405566970666193443619321021404939862072751198299739188985910258882842106685533"]},"3,9,0":{"tau":[0.09237704662228559,0.7103022872547312],"highPrecisionTau":["0.0923770466222855927468146610377207795318727263242924114280104416297565036590481927782299118242480982695571802556970418747","0.710302287254731225605235820021972156692635872661854655829691068670063970111675545620771635668398524219650737844695362815"]},"3,4,5":{"tau":[0.5534269290461918,0.7766771348721835],"highPrecisionTau":["0.553426929046191893534710484199618338782893332857577031743969672299338918452864135037564511369259618391997471197828592756","0.776677134872183460339073154146404459335218330435287881912422749528208060163868100160007736836871016790758786606475559093"]},"4,1,7":{"tau":[0.8076754772201781,0.816365852427291],"highPrecisionTau":["0.807675477220178084200471852478008388776226441358725456306100301514726754891086821624432635145849403773476843863275590801","0.816365852427291027867102321453876401095346990879829745369178333184063475660537796645453735676639624177014312096483708293"]},"3,5,4":{"tau":[0.46920249164962297,0.7781400105642429],"highPrecisionTau":["0.469202491649622971122376363987982954441434120168033456349637372117131614043409029365739379148826888277647142216314983683","0.77814001056424291655979883155693983805519024103353381182916615149510414084747173042857155644369865195409975444889795798"]},"3,7,2":{"tau":[0.29584386545179764,0.760621894720252],"highPrecisionTau":["0.295843865451797646473911676071645199622974093150354410337650209218221518995236305302265353752989893717851878895693094524","0.760621894720251997193702512232223678810131017748908829816253359042497185682217369907712260160994209798123277583506099402"]},"4,0,8":{"tau":[0.9237730721314442,0.7817792061082146],"highPrecisionTau":["0.923773072131444180349524687577818861167501127452308202328532416813716066879464603131285012349749654525829488000558720511","0.781779206108214662390621452840041108588868463711862264736632377637311656671394244958839288467704528592964010950936243034"]},"3,8,1":{"tau":[0.20134273358064064,0.7403926968587028],"highPrecisionTau":["0.201342733580640641348502954349724088122547249898829365160278266177225978898022004672589095622142581836450792445371180263","0.740392696858702799431659248240516964589923464936283144580436764073544485674507415385906446355997069017181312431228005384"]},"4,3,5":{"tau":[0.6072509094223427,0.8522435875817199],"highPrecisionTau":["0.607250909422342750339819842384881763287744436916504407628606043335143432062751649497579217564340040510822382344479999227","0.852243587581719942083356883899707724458384481716899128642500828777610625112680423185566460874429231681018081293449577802"]},"4,2,6":{"tau":[0.7047046261815285,0.8388114411121474],"highPrecisionTau":["0.704704626181528520839176867068784092497253676184546769749842534984760004055421295104732235481539290767632686468905938411","0.838811441112147410647715164278392822603124342943702462156556566777894246884524645271620670815576506153195358434030232221"]},"4,4,4":{"tau":[0.5119457148967927,0.8577668897456762],"highPrecisionTau":["0.511945714896792742731798748189323850249988441008197291983268581362741975121191580176284802651466883614340765200039926663","0.85776688974567613716881456620587129449999718829926232059567656159703403136022274722634093023486981002456795940671830107"]},"4,5,3":{"tau":[0.4164842720017269,0.8557361824671413],"highPrecisionTau":["0.416484272001726883225853072170571969687918186407077940601980428353267933213532793255555204751801000175886445107778331922","0.855736182467141351369870783768870684286393505437129513115291681782700036783727803257199815108137870695727261189155579616"]},"4,8,0":{"tau":[0.09675761593064516,0.7992614555813463],"highPrecisionTau":["0.0967576159306451609218846817122110334539784154883319916430727487189582579134338401952863714464718096847421522422158468565","0.799261455581346339185945513239269254520131379203910299721271756821101563362673874621897716553955879854668500745862877944"]},"4,6,2":{"tau":[0.3185389943807543,0.845985479878159],"highPrecisionTau":["0.318538994380754256763943798289680693371819278385117467076254972967426936871319563166952798323207715217979441238520168364","0.845985479878158995945995959626070368774838560101771856070746933946638168192809591590411281210093162130563787447977152781"]},"4,7,1":{"tau":[0.21464399104586068,0.8277657433354307],"highPrecisionTau":["0.21464399104586067465215168781914263704303683767349211531399030058628267504177329110051404570622577353174848022740199864","0.827765743335430784325621682262572183829768007211156071911795392473481104253885800223803436957690298182263267935464657449"]},"5,0,7":{"tau":[0.9181389606664684,0.8734950882199283],"highPrecisionTau":["0.918138960666468351811767744571160616248787190547878749482105346651210293541863164089940091853350060883625207146927107372","0.873495088219928312825649716160887620289396562656137703402030571974395470635967349207008817930534625878462617585220501381"]},"5,2,5":{"tau":[0.6765282082632454,0.9292751636040629],"highPrecisionTau":["0.676528208263245393438704013790687977716817800882111625514917270602431638311314772639260167132335346472686192417451433148","0.929275163604062968739181284259533577225135346039559162419681557717449858380766395176711811162153393672309873620034518268"]},"5,1,6":{"tau":[0.7909812804834715,0.9080483718333872],"highPrecisionTau":["0.790981280483471442510629664498206733547688925426230601135041024778466144956819495741850863990758906547174372652730697515","0.908048371833387182846969854162726535377888079467189346706198816537210250458701630229950345649881605142087192566526760214"]},"5,4,3":{"tau":[0.45871633900106784,0.9425856335304564],"highPrecisionTau":["0.45871633900106784602146827700726993712003526755987734434700565418053015533785345447219297807153046354134016768634808701","0.942585633530456388101331183746967125567303204789028671796203350515984421719454440646209858630890304616490980229023278023"]},"5,5,2":{"tau":[0.34871099524793453,0.9358269666824759],"highPrecisionTau":["0.348710995247934542338592966143561534744992587938472243336753517745113041540007623228179351137107658389770620016626714663","0.935826966682475903114240082431624771242283857611878411667820404021888555171722767009677692358942881089242909835186593604"]},"5,3,4":{"tau":[0.5670055965859698,0.9404467831979073],"highPrecisionTau":["0.567005596585969809937704780509815151576261240144727536751793033136805764289208127995021307414913968328269498734484480775","0.940446783197907320369109745775918638555383215217776057425057803763785265062655339628160826157496928749839547816586028707"]},"5,7,0":{"tau":[0.103828286013891,0.8919121526790303],"highPrecisionTau":["0.103828286013891007006121341982231998497544394429836996310550089431263613478526833264493346535618122470153821368478442767","0.891912152679030332609050749240483422621510823029474370100691604499437432450583997423076801016276758581726159580416832835"]},"5,6,1":{"tau":[0.23318093853812538,0.9195573356598994],"highPrecisionTau":["0.233180938538125390540850050631741211472388535399838098101692481736196060289741609468917997511794211452894369710565239834","0.919557335659899350999110583975408138935955858652434536019687392114787832623556159008484261554796183042098328922298654039"]},"6,0,6":{"tau":[0.9102673293722553,0.9717426501850467],"highPrecisionTau":["0.910267329372255266380551030483338049182641904188716267818491360376381564137178214976836635621129099595727449799516830238","0.971742650185046651246767931904217996001305172889983347051588726748895108234916703305656837273779085200838199720385554444"]},"6,1,5":{"tau":[0.768639298560527,1.0070452210760192],"highPrecisionTau":["0.768639298560526995061994171488309920816375053157899833660177478219429756554043387539528236859810602784248034608873064462","1.00704522107601929709006848299033634943698906286903198239637898415684770725548150547896960272225591119497939532492886085"]},"6,3,3":{"tau":[0.5142347154201922,1.035247234103994],"highPrecisionTau":["0.514234715420192168930665039587576945865381861235722738103058863004919044073350787743110096312104651525713499702260322831","1.03524723410399388605907180467368021963062373669048415405763719390058493099999553954761011023034344627019655339720509269"]},"6,5,1":{"tau":[0.2582622585682527,1.0186364748020347],"highPrecisionTau":["0.258262258568252694475253146777378751767752990557230476091355049036699099814550078121273885715780072514962448658920858407","1.01863647480203469524216220252475778506776619907228219128817978726513278313990720015869331270846959407631482376562777303"]},"6,2,4":{"tau":[0.6394000799407816,1.0270330607660314],"highPrecisionTau":["0.639400079940781611233972507977568675728822828325742221353543860625789780362486253219528277076811365031229092472180099293","1.02703306076603129302137747972501297090681074127140294694378008452748287268566150613191246947182685994982068573994321899"]},"6,4,2":{"tau":[0.3887033506187352,1.0325673132171065],"highPrecisionTau":["0.388703350618735168919312666942576971396515019819498653272268690837340680941380736400338662063427493427783221790341538618","1.03256731321710647666097962994123472647234485227036632048956354419642353942131523761148592172375334889883696458174680031"]},"7,1,4":{"tau":[0.7384032713687888,1.1170034163666092],"highPrecisionTau":["0.738403271368788805410935008962406154645023629818116757058366706932682550148007523361969378202216966854191113867205442098","1.11700341636660907588909516343829816522136390295273929620094552413058220244409083430656988512149599149813245449413116506"]},"6,6,0":{"tau":[0.11396588674782296,0.9914763713774831],"highPrecisionTau":["0.113965886747822957763169066188505925750709486678627786777863586094021580977465001787287033037335035768040328162088544196","0.991476371377483155981318247035710117354374563933622806169652391573266216174294756638592790622822024760762862599581764496"]},"7,0,5":{"tau":[0.899428322101563,1.080283576640981],"highPrecisionTau":["0.899428322101563032635972838040878298768274480321923219538641694172260281416319780452997962877339854418895881645656791898","1.0802835766409810745995007585451215847700573112510809440853141125463804516504170391041379969282926991025496889870607587"]},"7,4,1":{"tau":[0.29246243639143316,1.1284494547239543],"highPrecisionTau":["0.292462436391433152985428421815492988011043362027741330675451365454392449214758883962001640024533013473797102564195760793","1.12844945472395423913826578358890398623171896065497352863493198125038285220627711989551311528266853342391715882464474678"]},"7,5,0":{"tau":[0.12817976132532832,1.101750571486296],"highPrecisionTau":["0.128179761325328312693713282214761019060770142510179563973467129306469399957548994098156948631929721960188583818161920144","1.10175057148629617029525293027131394861881435796715450073606314431512037254247374426919421405485214104769960487569428897"]},"8,2,2":{"tau":[0.5188741959385742,1.2581081665614666],"highPrecisionTau":["0.518874195938574169212904100495495144183266192037420494273489498511645886244831263043538562332553457848810817901616945875","1.25810816656146658356872768009645460317251416055017194713318770914224012917807972339264210133856959349916815256457548108"]},"8,1,3":{"tau":[0.6960778051448823,1.2429912873399964],"highPrecisionTau":["0.696077805144882299816150003313025011537441764935196529894019137476202727211113923842591275788167646006201461158696687157","1.24299128733999631048921294205979212088818810081184120756830032687929059861015684233553920531622459401312412608524757625"]},"7,3,2":{"tau":[0.4427818248642249,1.138943320204907],"highPrecisionTau":["0.442781824864224894622500903174468733914238068912981041313137145349249371855471214092275600737704041705319807550180242501","1.13894332020490704731233758095381288736757898624074523813872230116580400890165337966651542606160724230212819677830623522"]},"8,0,4":{"tau":[0.8842004733615165,1.2042963284376886],"highPrecisionTau":["0.884200473361516554933801007552249601342347010785452242503832297678736650022365084832262316400109221858701931864525479529","1.20429632843768856360732884788867469160676836677785252470982678960200156782683356716169535802432246573160855805905211155"]},"8,3,1":{"tau":[0.34060555405446236,1.2536348429194102],"highPrecisionTau":["0.340605554054462340755391261700097376453962424268235180080281578965448319217225151435421110908949199723313198278751534945","1.25363484291941013903578280255807469111749289724508621989423474724573021626854031653062109152128398225761547753931538774"]},"7,2,3":{"tau":[0.5893368984217073,1.135268015313458],"highPrecisionTau":["0.589336898421707323004618564795810488272228079531725893228995637962121416404418999710964383262821302020211649667312853261","1.13526801531345800086813680691982577501218392110874538933986021759851403270305833748034178549567010483590531864198686336"]},"8,4,0":{"tau":[0.14847298442375348,1.2279828005686262],"highPrecisionTau":["0.148472984423753491890994439509946187476414680175808454066175494798537283369270854418225628677225304114501164449923563828","1.22798280056862618602146086743898004070291835448996419692985524698771743760739258460812660777874445274868729834258719153"]},"9,2,1":{"tau":[0.4123564557928534,1.400873228402031],"highPrecisionTau":["0.412356455792853378104246628540922890932849683572010512739263726079096824545965665595914556922985907193769619875700431263","1.40087322840203103133153953129467044232772855540534142866753720773359458282053832812518564277573560483471236679671749215"]},"10,0,2":{"tau":[0.8256934887583117,1.5380237527888825],"highPrecisionTau":["0.825693488758311693488369974383152407657273056045829467151732701376341490630365504234548588974838746695442232735716163599","1.53802375278888254176138914814604021776600249408346214109892286124616949307772597114197153364978235938204431860996174553"]},"9,0,3":{"tau":[0.8617489121425846,1.3519495966088786],"highPrecisionTau":["0.861748912142584617165328507955124688497165265193276912841624778156569666290951254663226693885260582646134759964489471913","1.35194959660887851223416061655798849398330229403835812341812857003750398103531358594698256339279157876169060584383491593"]},"10,2,0":{"tau":[0.22869972733237606,1.5675148718751424],"highPrecisionTau":["0.228699727332376059870790766170131626993320872121184202407790752520262757533205146503243629519907511626614433275148673986","1.56751487187514242548836999156795488109710485319942976167408974326466538368974164908160953244292150230320656660879155852"]},"10,1,1":{"tau":[0.5299390795730483,1.5795936525615895],"highPrecisionTau":["0.529939079573048349334669476579011862229761382482138776046947765899676962834929448434866165907469127854335579046483686661","1.57959365256158943970880454784514518951593471146011589319382440695163196921245208556117466402472000739953293692019851226"]},"9,1,2":{"tau":[0.6332081622850778,1.3928204500789263],"highPrecisionTau":["0.633208162285077824194060251671028133884175223677500119783882215459812172831169980349784217938446723804882947663409184494","1.39282045007892628831237403142542194829882244525549518319977544224150285155763537005639748646018821092474732551976385137"]},"9,3,0":{"tau":[0.1788982106282843,1.3783934288232884],"highPrecisionTau":["0.178898210628284297871716985897372372555351928342500920144209628503860356089116061184229630322498047574784388033621510269","1.37839342882328848182897304030557532561976736994535989101593373331696307424992665121917113441277320679920131698378960038"]},"11,1,0":{"tau":[0.32415833173286235,1.8248385905393214],"highPrecisionTau":["0.324158331732862330395768715727168691998311973553820958019648971200411102153986877459139128047334730622623864054227493231","1.82483859053932148960087302588912086960508209875542928041576051304018391763864528361250830530759698510792703947100913512"]},"11,0,1":{"tau":[0.7580384357058384,1.794548485328673],"highPrecisionTau":["0.758038435705838432542904880561936084468261124010884078955375910513770133346875741874155425197342697726973807523121475132","1.79454848532867298126416354409305866944687164399314474942363189038449479780635299657737633337673583736424298318807073336"]},"0,1,12":{"tau":[0.8329609711818248,0.4465096376676281],"highPrecisionTau":["0.832960971181824802647565707228348813803503403809431426740951733110767099066997365892571776383603007757522745885436511938","0.44650963766762806887591187324981856097831324140082516374128711669371437063065492671547045219213813217904950901261103155"]},"0,4,9":{"tau":[0.6476447475101138,0.5256011961827077],"highPrecisionTau":["0.647644747510113799285506731895623499104163412431874500021015041126052681310626500759195520398914938125001614552788968756","0.525601196182707752900078543976291667139706714642034004299442908516629654780626195736577097767283489862302856522906340646"]},"12,0,0":{"tau":[0.5781860351724207,2.2147573220972965],"highPrecisionTau":["0.578186035172420762243915168836182079092625414823323805806634165621325417082986850148497337237328511364582603138998148872","2.21475732209729667736511319616002468183060146041741695838422430965946781320131225033460966869176820061834133601942828205"]},"0,0,13":{"tau":[0.9198036951384831,0.38859420020386826],"highPrecisionTau":["0.919803695138483130685835370707397228157664500125216968518103525936986221002922798774301959871583389367895936321591426217","0.388594200203868266117678747127143702231479230570292819982137408025956230308690545305102321026050721981741166248130891232"]},"0,2,11":{"tau":[0.765231808139672,0.48294160786474505],"highPrecisionTau":["0.765231808139672008251467994915467464360470826318501147584037957941472981518343159189500829277630430418526255376482738533","0.482941607864745075498177451207541527697787116679314492618221340608840004692158429025854661341585189195270576480106164303"]},"0,6,7":{"tau":[0.5384833370976886,0.5430269807890544],"highPrecisionTau":["0.538483337097688568129691387351372723845933709239062818256845315648310140212903920622211700416384988332191149489975688872","0.543026980789054383368043661387524071552682468169591549285855643743992135526499090592327906220482478600144139177326222152"]},"0,5,8":{"tau":[0.592571357971151,0.536964313980208],"highPrecisionTau":["0.592571357971151034045108858862113611838048155233408385680523147835165443882322629805167454535522467311368421879151138211","0.536964313980208084309867196881703235715951113472490935856121509165179483339834273884223870589272794828774074839904158764"]},"0,3,10":{"tau":[0.7046544185986553,0.5081316561666004],"highPrecisionTau":["0.704654418598655321500375620846316380434330324408071858883018002561963315857647410972386289896392318555508856783141828605","0.508131656166600429141953383188674004495357862839641827146969025461028021222584128861722201819618892346365340730608830774"]},"0,7,6":{"tau":[0.48469779304702043,0.5441483763149451],"highPrecisionTau":["0.48469779304702041150255477438515783066916352932899734988489698603119012121274509560086427992039229719926360288106787895","0.544148376314945106533159120601177450292830829766948520224106051351982410744565552394481032568580627308352418994813280178"]},"0,10,3":{"tau":[0.31859463224142,0.5166048967342561],"highPrecisionTau":["0.318594632241419968937207024717265217118177084776252501570989826943032098811433210169283060272025812894493001384185645227","0.516604896734256141150544548721716091036360826625759939455528391870415874755320741457110143639715553159180058337091816722"]},"0,8,5":{"tau":[0.43060853778846203,0.5403668820904648],"highPrecisionTau":["0.43060853778846203512814436341092600754444984517594447362485595221165811389079693111652070808865086891667099854251792617","0.540366882090464751521257839216770007591021087678542599939635110715827233295981093323768190264113242838329629703607894971"]},"0,9,4":{"tau":[0.3755459210507888,0.5314111722797763],"highPrecisionTau":["0.375545921050788831390788623225778026303523290419564752394784273589547724331600856936899509175719646228124875553694034304","0.531411172279776264611846229948875512362532302736661820267117969179551838818222618213349341128778327801835222323638127512"]},"0,11,2":{"tau":[0.25821250432335996,0.49457578527290263],"highPrecisionTau":["0.258212504323359944649675426149833277429773757380090044807692430557656566198681107739442013503916969792616454984522831672","0.494575785272902616551965250756094742984364095392167106585477744340268888964985535801581308902072046429968694971313268162"]},"0,13,0":{"tau":[0.10636619295018897,0.41211240299877916],"highPrecisionTau":["0.106366192950188967907547769444482545986654992164712507481065013994135562121343687973498076448587300566067813406590518306","0.412112402998779184754665556766155153614792295241081592439010753933797399460563732373836200131963630255517844192956577002"]},"0,12,1":{"tau":[0.19107785082052353,0.4623925538783791],"highPrecisionTau":["0.191077850820523518741335333950066800556162767054254188163168868261627837702599978798179927721183381011244404285876606665","0.462392553878379096917707699391714210467874128277708138566963946062605751117514024501388624779517018460647015968008807535"]},"1,1,11":{"tau":[0.8395391568589752,0.5415510366202874],"highPrecisionTau":["0.8395391568589751678694506621305891989234671384110887088785389317936100960201512057175757803074885452540845115127021924","0.541551036620287333240881377754406006422542190338926431690193895660639097004727121032412155697901459161149591240790916718"]},"1,2,10":{"tau":[0.7649932800736748,0.5711301529848141],"highPrecisionTau":["0.764993280073674854380888795502575654325130840576043696109158787349022568560479847383193182313237551753507432504192972844","0.571130152984814042512632972561723967496884340294975977214472050614857631385598061597886176230204682303386761433488073081"]},"1,3,9":{"tau":[0.6974709499288771,0.592085167646258],"highPrecisionTau":["0.69747094992887717536855988795144613609551522503469904279065897215998491854675790884747510027563809112964505963187114112","0.592085167646257999777952947532230750848091210039579486912019279883493315044420366952376211301527626313242909372381583785"]},"1,0,12":{"tau":[0.9306214528799485,0.49855800794515936],"highPrecisionTau":["0.930621452879948431319246150511016294023806183239084738678396890808748806340569072325387011533147986243887034196427381975","0.49855800794515937454177166529553353003967521105321687471724163888433500237900699354968378387811762554513129577531931948"]},"1,4,8":{"tau":[0.6335740836081337,0.606403879456759],"highPrecisionTau":["0.63357408360813366585943103046975651752017228257142797888545109578328037516810260229416778797975255625105668363188387753","0.606403879456759031428569193353687932031966742417184115015523937171981840694102186537740208302406602749542267415774615578"]},"1,5,7":{"tau":[0.5715950500450319,0.6150825713144787],"highPrecisionTau":["0.571595050045031852059603520902915043381188842266886478274218414624923762012585338637854787791233193320962241899504116498","0.615082571314478790547283039274137158573081190377513456733027512072915061188451174838306924844245833652822939220242735195"]},"1,6,6":{"tau":[0.5104496559705697,0.6186161769463604],"highPrecisionTau":["0.510449655970569680494319347603915449371209653406261343954560682558730950130740324291285947885053653050495889397068849469","0.618616176946360434340539251479610488690556413090099467753594890703996240546815597980952786711462782196144480181408914214"]},"1,8,4":{"tau":[0.38722928121100264,0.6106499513154057],"highPrecisionTau":["0.387229281211002662326948474948764777314310222634542459131210822090955006200602544112120972946173819638713695960022866759","0.610649951315405654905321965727149159870117212476672245141333889709645071064373492822045388624386162230744872282586369414"]},"1,9,3":{"tau":[0.3232228092836815,0.5986339983018435],"highPrecisionTau":["0.323222809283681526100613375055066161497712403236228137881167752145996897651141578750422888245882886658708081495189946129","0.59863399830184350148739494026707567396044489847218023102019957903906127642956890791354777380351282369908467431971834607"]},"1,7,5":{"tau":[0.44927969119942723,0.6171737930812051],"highPrecisionTau":["0.44927969119942723538723248214080718249792766868983188408245434509177903306989154594335201141395229364288054502980264099","0.617173793081205125350667050454805451551556610917426600924962204081676024608458732158872772640150924344783124519195881306"]},"1,11,1":{"tau":[0.18093374220121114,0.5539247270598293],"highPrecisionTau":["0.180933742201211127979051497876010246523305623579602775175817608432543232980387192215348993419956425240167888997996464303","0.553924727059829331257526344792492316939888476901883128719566202749637236381536070006658499800014512813181301165757111252"]},"1,12,0":{"tau":[0.08986515291281698,0.5161092071550749],"highPrecisionTau":["0.0898651529128169768277789858657281811213316060606184562631024881591967025076812312371365308227391354154413290898709285443","0.516109207155074887909189872278764514842972396210170393137731855518397661373702969892850191675618589148531873138948924491"]},"1,10,2":{"tau":[0.25557368129679947,0.5802753015425673],"highPrecisionTau":["0.255573681296799459966029422370018549964362587474961943937664348551322368800220719574675960889212675845363333424942932885","0.58027530154256738978263457361518156835456410130295966736325334054801087653196306484133968097815632547721541692197249487"]},"2,0,11":{"tau":[0.9328645482697911,0.5883923257027699],"highPrecisionTau":["0.932864548269791112432015284904060616567784466621533375978092392933643727296019861435237980184774143548698292135011945969","0.588392325702769864950851337423823717516265909237074663798439666532325029129364300184176806344925286876714467588739128921"]},"2,1,10":{"tau":[0.8378968649192508,0.6254805064568096],"highPrecisionTau":["0.837896864919250757788563723486778216012436611954376917902952672344006282232193758254632631499013508028909364070165261475","0.625480506456809602065879744355534869831074332090256433263704995098394394170328466362145872365434766526944880531825051509"]},"2,2,9":{"tau":[0.7575382442776027,0.6514673883651478],"highPrecisionTau":["0.757538244277602792774287893772794249964568129612332206131497323042076785198212339921409500992791751931501733331025528486","0.651467388365147835959640593970925218365827084272562967279154198877541904909529413820270156682886289909780634271384877392"]},"2,4,7":{"tau":[0.6132317977090082,0.6816902887688328],"highPrecisionTau":["0.613231797709008103995280718736664417575848329308212103817076630101580753840317361697072756363512656258370386732106736593","0.681690288768832739999809708288216665465635635585051647134890490407426768455611405132758702130894550144293799091620107859"]},"2,8,3":{"tau":[0.33652973402243547,0.6753958537445383],"highPrecisionTau":["0.336529734022435453775672778322401562651889278210810614276587251025165996752287740965361858792221563280240308367799156456","0.675395853744538237770720573292564311948525257894541657517523293210803596557158372014304036762949486818957135123587862554"]},"2,6,5":{"tau":[0.47609374329106136,0.6891565147028521],"highPrecisionTau":["0.476093743291061362314696003124959031648346541274316259078178842594312087764874203360637817777496756203702475668759553702","0.689156514702852083597802623033966188796457090645344744595526737668277877821832822708774638995745742345862533073039347367"]},"2,3,8":{"tau":[0.6836998929881162,0.6697303165306909],"highPrecisionTau":["0.683699892988116225760127868505824751463797731275428274449010545501278781215673727913771555851589469500128154298204513914","0.669730316530690969830191643008292254232280302409415842127452470771537886147913166413690173497029754727071949541307219417"]},"2,7,4":{"tau":[0.40719132848659084,0.6850156522744049],"highPrecisionTau":["0.407191328486590850767595243671867600228333529881011780317584352619077111427125975265529027481449718847142829639381230963","0.685015652274404888370458267281604609885606289857010848643968415490046026958353525881550485133012685892305134842409770612"]},"2,5,6":{"tau":[0.544424791209615,0.6880591507712184],"highPrecisionTau":["0.544424791209614925865137255878886754550817581914764878146558028123119060288998663537051979210207135280534866031604361453","0.688059150771218357220260076141688297523836047722837538026075173919859540898396276837388046653920305710823654898376443228"]},"2,10,1":{"tau":[0.1816040184653351,0.6368534634008758],"highPrecisionTau":["0.181604018465335123266393859068711731091210447744211586943944243852409658822282797747450080782099526221008964464940424459","0.636853463400875813299097600856656172849231246265731876497261685581986583384581011352960594015176014761878847059334836715"]},"2,9,2":{"tau":[0.26239231791131307,0.6597151426550364],"highPrecisionTau":["0.262392317911313087609251875537282770341591028514625389891922853111058197533752310142667511000969161617770680608275627896","0.65971514265503633922677960858420788544827087095423967037236294574017522178798420724324180309409922436461936564206865031"]},"3,0,10":{"tau":[0.9320306311005777,0.6717638250368093],"highPrecisionTau":["0.932030631100577744178481946559444557881438639712164913806857364291718012487204055879259349358034080550987317021448830057","0.671763825036809345608661333661886967862052390830890410149943888539667989434894971888999031241631664336430112424996676167"]},"3,1,9":{"tau":[0.8318548767504345,0.7060699124207054],"highPrecisionTau":["0.831854876750434473300566135207720044634641372833725140978806492885621923393865149346065515683139766381615315050933896191","0.706069912420705349486766088002657474429466862839029762262738918059062753534429875998644714254676713753916720199944549973"]},"2,11,0":{"tau":[0.08589451914703533,0.6046251300568509],"highPrecisionTau":["0.0858945191470353240631012554273297275140171595108374403290554058260205358889698562530431757676362381360878329516444537224","0.604625130056850937240936956182815257103720125826284601471359898911598222917846732646124813111664201658161068084727279149"]},"3,3,7":{"tau":[0.6644639113650598,0.7461447326606251],"highPrecisionTau":["0.664463911365059723692362007393944259774320648163209381481294400838015115616049373214902243050341675766643120905541435402","0.746144732660625072181333565524685191882379798867602447562022378639003061031494353106368434830045820800282874710117461266"]},"3,2,8":{"tau":[0.7451293975708193,0.7298867653926417],"highPrecisionTau":["0.745129397570819300786410352147584214251128373016618878157043689332677440315921419683927788258081046272154099255311458417","0.729886765392641731929767395802731671817090997121997285480743042059359879542628837540496544067097855607654853252540379828"]},"3,4,6":{"tau":[0.5868277711127584,0.7560207259804486],"highPrecisionTau":["0.586827771112758427796998291549799251026590096889390327412745759245079077053678947681500791286739354225436749216265635877","0.756020725980448576829367556924096497401922081624341311128699300682780788427798748930597639539082034917532048631179751693"]},"3,8,2":{"tau":[0.27508120680441,0.7376648861610304],"highPrecisionTau":["0.275081206804410013885274722358652396704177580980495289919725377236716899904046200511451874093439050631677917543041935678","0.737664886161030450999985755617849506657595148163456042809977184659813871333821783595502686689450244566486433291645178621"]},"3,6,4":{"tau":[0.43401905695895615,0.7584807014619963],"highPrecisionTau":["0.434019056958956126986493207669639610888838503524855324106713023773774794734627446320126645910934493682289783237385089757","0.758480701461996360180231684948752412069950968885304184240841541494948230534980279280096025374172674428553302205042951601"]},"3,5,5":{"tau":[0.5104608024582875,0.7600690958426131],"highPrecisionTau":["0.510460802458287546317606751095385861796795567226092339124852287836320183182551130352097461068368790281784498465510031314","0.76006909584261302941378806947462664989919676789420362330985465863261876612839169136560095493026285652909523910333172413"]},"3,7,3":{"tau":[0.35615296250986356,0.7511529758816169],"highPrecisionTau":["0.356152962509863534188231638375153075739592715520556041739804102451693864688433755084400416383768025672243254191220535865","0.751152975881616859986112112437311163830640043060269048353293421428205921909972379252223761161120863414875190734937486151"]},"3,10,0":{"tau":[0.08636368702472443,0.6878238923846131],"highPrecisionTau":["0.0863636870247244265789536245661538275169279278850149803073366781141718925270851069504915942134087020079994718392913062571","0.687823892384613048111240262815379167882515544938386990749775412229727938188600693701327123607280326108487994960884247403"]},"3,9,1":{"tau":[0.18771204578327447,0.7171398804255601],"highPrecisionTau":["0.18771204578327447602962037574659667999023299271894911726494523123293858661543160438095025879183810139385399276848369543","0.717139880425560047082524545031406150334010354190724363161191153233437026440383377247326629436058263367293109580965227591"]},"4,1,8":{"tau":[0.8223230746485642,0.7869839430850978],"highPrecisionTau":["0.822323074648564189942268526008206855942756164337093422174974908267220282466533768075720997277239069261986758635926171584","0.786983943085097832557462034006024996840661599807125773156305400112091456852946660601411211990291826456623403396757506514"]},"4,0,9":{"tau":[0.9293185934879381,0.7539039249246032],"highPrecisionTau":["0.929318593487938135997159447666389160342532328364854848336001612295047682237258484526953363690437342701401824143256674707","0.753903924924603238821894717063500148648407928476882774084575429954536616156518110011924126871005135943010778237177284303"]},"4,8,1":{"tau":[0.1979860551557774,0.798047166630608],"highPrecisionTau":["0.197986055155777407927226292666320316781656763136192407181899142604892426074404734796599878922110644518646596863500532608","0.798047166630608088129608138695006845214813239517831445995921234210790567842953816607546457484911639126448858174006443764"]},"4,2,7":{"tau":[0.7280924963804316,0.8093724225581056],"highPrecisionTau":["0.728092496380431555527474008963714572070546533659022350409816753732722794841720361585216808204283803442402379973227246523","0.809372422558105600091830932267938076509935798482851577563441497738340512264486754164492510603128571376074817768546967656"]},"4,3,6":{"tau":[0.6395124987102401,0.8239118577525792],"highPrecisionTau":["0.639512498710240038603808426570816926512611771934965677728867148300647051449839067730276276795645788197244427031883282145","0.823911857752579277753776045869290691641003683098730733598385763331668049639190493521753689934857300028528344981590920737"]},"4,5,4":{"tau":[0.4682775613123293,0.8330926989492609],"highPrecisionTau":["0.46827756131232927246897875696589872686625364101525394916077370136642127461891856432218808223293901256087687533891537223","0.833092698949260841521756926787122415356035086572114979784632910821350918755880397585026465063187989343435572515450695515"]},"4,7,2":{"tau":[0.29302749339260387,0.8168047756240928],"highPrecisionTau":["0.293027493392603864400537266929664110647874167546801790342959496645081910213648627962752458504214024160541529799538168131","0.816804775624092855049646376352697288832856979016190687971335436890218707069332734747549516932760550182118403944043018605"]},"4,4,5":{"tau":[0.5535419604520415,0.8316706280357204],"highPrecisionTau":["0.553541960452041515077931567079471858000126543034661924007439500510382710214384122081698966198209173409151088679930608544","0.831670628035720423958253582469205228420595780930868162024015728306743222502239530269454870496386368556568078764892935064"]},"4,9,0":{"tau":[0.0894758080484876,0.770325978698775],"highPrecisionTau":["0.0894758080484876018141778607832361166723160333405115859886404741941171303417963947954858572175998356546537032551413324884","0.770325978698774930953923187498314942457426324050155503930232398409883649851602221855595080521325380016936322028416169893"]},"4,6,3":{"tau":[0.38208285341924186,0.8282356884374248],"highPrecisionTau":["0.382082853419241846418848345167100128159740859350637611655478198290946555531657316489783980757552126060080750018375046742","0.828235688437424810814368772924861684114806702922719025459646608549519516935437997314560912790161289767193390586377644454"]},"5,0,8":{"tau":[0.9250048698214043,0.8379733009222077],"highPrecisionTau":["0.925004869821404326125806234262240344713828960647003485539925654096189550668881180033493998885339509898432160284177051593","0.837973300922207703011260076923904648934167580124974035979967171994303560473121653784109396281066868612933898575827392241"]},"5,4,4":{"tau":[0.5116419201270405,0.9102906338951099],"highPrecisionTau":["0.511641920127040442970651497788367229271337704340476171193466840957681294716001125845299521704963447729473609968566127389","0.910290633895109862261506956494097602442792640142525909841143869336568175779221044769705411822528984499634851335694051521"]},"5,3,5":{"tau":[0.6077380296423583,0.9049808837445056],"highPrecisionTau":["0.607738029642358302873176611027713117400732267563982573774202640934300249922264559929850412254320910343380508658271106781","0.904980883744505611689564936786155440449208098423195826701329328131933514109606587370147393764704906699626309654423380324"]},"5,2,6":{"tau":[0.7058779441111815,0.892157411817112],"highPrecisionTau":["0.705877944111181532456510717711603936594944924424245069798623088620077023662110540184717913194422574890609952712755348295","0.89215741181711191850960900944815162547113912298764645068375514843629593744502899762349886961243947912277142315102253723"]},"5,1,7":{"tau":[0.8092553374992725,0.8708011191689307],"highPrecisionTau":["0.809255337499272517886511730215513619147691186378372903856119537339907159408907119607088030227585879875711995472858128976","0.870801119168930657783898321827349605612234563079491391514968462024221010929468701015395471011970542584057009385377404609"]},"5,5,3":{"tau":[0.4153840554750778,0.9084149509442729],"highPrecisionTau":["0.415384055475077809534967147010032937048932284949584073335038775322952594224568857347758514443482346511326799800167403191","0.908414950944272973549579249423653697180399797250796796307238396760132774785956857181304387329212425839767614381766678102"]},"5,6,2":{"tau":[0.31673473723510703,0.8992095246179124],"highPrecisionTau":["0.31673473723510703553069350571460708786278759761808882417885233183734721190515057597630399746549280268363055992992484372","0.899209524617912352282962903097028362586240188209786276957267154303328079556599877961519115736090595621302058816129097177"]},"5,8,0":{"tau":[0.09479690615658962,0.8551219559111504],"highPrecisionTau":["0.0947969061565896104652044891504658998557231437127788926383072991363469059829102717694270554004737323138342689516032343633","0.855121955911150325151245928603804754557074721072939542042400338683707245737385901016071734590640788366717651233557520723"]},"5,7,1":{"tau":[0.21239822891844024,0.8820012923501439],"highPrecisionTau":["0.212398228918440228191515684018515773436229483772328201132915342111582904882880513999425633981419687128853063563709776443","0.882001292350143903552181241619054715037667037318733998161735479925440554250245096953146904875247561396166445841885690208"]},"6,0,7":{"tau":[0.918983366059116,0.9266037514019184],"highPrecisionTau":["0.918983366059116060087879807093273755511157641796323587202264334924058488708919964865050133363644431694605223183100924692","0.926603751401918469523113917825929748079615060002280913598931464319467254327422306872218368575370637582210385640404666493"]},"6,1,6":{"tau":[0.7920089240685594,0.9599063268589241],"highPrecisionTau":["0.792008924068559350768413094620881537448599659239929015154085755901652067245894998264077357389124647971873935031525324302","0.959906326858924071722010053864723453728179165385258047411959009411079903099160113356810547079716864438247479465557048971"]},"6,3,4":{"tau":[0.5670863598517215,0.9911909296457426],"highPrecisionTau":["0.567086359851721482085122840168151125643517697398740196902694370186422549824995985908785653372101035735269860472735157827","0.991190929645742576290851316439180735220152512965304826820839620207699019671744121961092792103263844770082304992464400714"]},"6,2,5":{"tau":[0.6771844902243679,0.9803826190863811],"highPrecisionTau":["0.677184490224367915716980601424676911014110880854256265469835065687466558775413106877673283919144356036646466726585480633","0.980382619086381088304358823533709612517874373894095941840823893427994864259033499616285376217940701822690628708032680561"]},"6,5,2":{"tau":[0.34757455210792376,0.9868670873787456],"highPrecisionTau":["0.34757455210792374402962660221621723225438872569555704022401428659434755847329915197774077968138038306501429974124653393","0.98686708737874565393536912154539306673062840741660417115509451450782114134738409051996199465183318342853635950703498751"]},"6,4,3":{"tau":[0.4581640681695595,0.9933080860317534],"highPrecisionTau":["0.458164068169559525513775879522036380620320510543473335469741945847049242750390452017633444324643226641005566766862815268","0.993308086031753325972154917700405038808894502130413193902569037039838924764981530019426811960470218234117891446821404763"]},"6,6,1":{"tau":[0.23165373714813736,0.9712927043004851],"highPrecisionTau":["0.231653737148137354866993663116932166361709610560142168062298670519249829171480983470817453125139643441965331523577605658","0.971292704300485139628823104943870399657978375274508116924784632793318049516970830722322999716671348788989658667590721331"]},"6,7,0":{"tau":[0.10244943042988142,0.9448019983083755],"highPrecisionTau":["0.102449430429881411099022043840545149292255383394144163751310156475295293147174192338012886156136692121037674335136233641","0.944801998308375519127138490123223428026024640300773764540520937062747297571366896264788617074089422626826559132253207655"]},"7,2,4":{"tau":[0.6397010767375967,1.07646728720198],"highPrecisionTau":["0.639701076737596668909763748292827964562740162787053492911202919884402555112030193037222078179648421246541436614072602175","1.07646728720198002136061390048319129023842167538551309709611139682130108215607467605425404704110104214852914375614641333"]},"7,5,1":{"tau":[0.2572308978909206,1.0684985689722353],"highPrecisionTau":["0.257230897890920568341659014135382202429048557340296038497348443451412760634077731505333754554080954018923334893677069567","1.06849856897223521316660806591110765131101545795839976715060491932426876310126581236189556810675849544846444160172429722"]},"7,6,0":{"tau":[0.11297803148654446,1.0421702176562773],"highPrecisionTau":["0.112978031486544463418945862792511851990052481131351899585522530072651858719544521558166831338904166936061402806151595936","1.04217021765627720667198024002536047208098810026144756973205930713653850600235540131954763605844952739775745897186218993"]},"7,4,2":{"tau":[0.3880297928810439,1.0819671105594861],"highPrecisionTau":["0.388029792881043914303660384166087067098653289565570379087259281663620335908562101946419745466871351320864157881977318696","1.08196711055948615687899630145284171694138507664948902875152469173770522017319782330774098652832434251321147684104052729"]},"7,0,6":{"tau":[0.910852816984025,1.0225831728802406],"highPrecisionTau":["0.910852816984025008361071821690569103824103255007947904571455418170979217812506417962547654129912289809431072241316837438","1.02258317288024066052208173703177176002793328252112802248425831169915015648116764715341748264470656640812120483651514889"]},"7,3,3":{"tau":[0.5140500392122411,1.084509113119114],"highPrecisionTau":["0.514050039212241073342137860374557724105822332457356954382596530648159116684959991777211586926925996745272829500402762823","1.08450911311911406820656009189666962638805908062528325117455060154958026527651692928505121013456950836018584295204204998"]},"7,1,5":{"tau":[0.7692877456230733,1.0569825661020904],"highPrecisionTau":["0.76928774562307328604150640709855460350608373298563912922105745277642308770973587847087169923086024346488953711791184332","1.05698256610209050516595705021039633970371648009798236745011248463031070403712671664530951081147335765334076342722441761"]},"8,0,5":{"tau":[0.8998302333743319,1.1294226909475802],"highPrecisionTau":["0.899830233374331845352329167909830662011430402952649729859695729769340011176906171049137439910392882187026502301546942562","1.12942269094758018878976963506884951503319147837293933833163786687861873715381727966576519192931741935260860599361321458"]},"8,1,4":{"tau":[0.7387803469487564,1.165494034977333],"highPrecisionTau":["0.73878034694875637357019626335680872485019749875580533572351597636204280512542232309259695628145933100550896625236861465","1.16549403497733296923791060975221175015694046280607615406102384048507540292402817296131698891829526314745941290739858973"]},"8,3,2":{"tau":[0.442437256286368,1.1871099421760405],"highPrecisionTau":["0.442437256286368009143739138528826964465845446345187292266070751056752472662181754343100945939638247628546240770468000252","1.18710994217604039809667808605487788836955670817208053067863307535258984664626340982934943205043541964751949425955482892"]},"8,5,0":{"tau":[0.12747140473958643,1.1507914742417054],"highPrecisionTau":["0.127471404739586430915388007703823368690938399365740119907036604886873606048845127919089889374885464118928530905521277654","1.15079147424170544162206170806832325390132109044864097900188731819449169211706831114768986522054646993109434798194708393"]},"8,4,1":{"tau":[0.2917893755720934,1.1768955181059753],"highPrecisionTau":["0.291789375572093405582008345013773015587211608181012137602139388152070454612176005773553306278896688021711931958450594704","1.17689551810597529155609969809218115800892927567770355536483684489842264464391535845277188005901960198322133233158579872"]},"8,2,3":{"tau":[0.5893904025804179,1.1834482865633247],"highPrecisionTau":["0.589390402580417892090628297995206629901594138181763573614512761270134559582509277201760001076428710611190557424616811331","1.18344828656332462886426063303607264767330589232335569828751080479935562567583198415194164010854900075909226264540126071"]},"9,2,2":{"tau":[0.5187611525015304,1.3053643899458511],"highPrecisionTau":["0.518761152501530400849907268508149218668317105575878151804166259024285507376758094265622263685425268251530584211604970326","1.30536438994585111728639474208677883521556380228765478781089521210425951949061589047268820408678476604662261910379220373"]},"9,3,1":{"tau":[0.34019923134213154,1.3010176988254518],"highPrecisionTau":["0.340199231342131553366343875533100401457500102878591238206988331132898096971957220652501498692944388873789807332304352195","1.30101769882545176868660564518465800011250124489508342726144398016091847670157420863333010854305389796763619429104487895"]},"9,0,4":{"tau":[0.8844663361603385,1.2521532029119602],"highPrecisionTau":["0.884466336160338587989528292433128547116155312357896865622532522680316589045177618175033333245056197026090538437474564595","1.2521532029119601571247393814987453089050354695693366605595788502604675681284260614244296219304705184944652032130437021"]},"9,1,3":{"tau":[0.6962563510836425,1.2903980737965617],"highPrecisionTau":["0.696256351083642466380517577380963326788037548605185923471585545921999349547915790939652293496309042788133958615889077626","1.2903980737965616471763621893141051485234898924978850308440187788982393893873215506564350025132676139027961498244513281"]},"9,4,0":{"tau":[0.14797429740283394,1.2757755965903872],"highPrecisionTau":["0.147974297402833929830960852082617735469367330341243146283781280002057587710158314008588211330152936771966022102666019862","1.27577559659038721858795035803871368188634867732824220966215143890656762594739910979572294719424517720479632029422738064"]},"10,0,3":{"tau":[0.8619107637740004,1.3988518670956993],"highPrecisionTau":["0.861910763774000452528383011073285979889883310609056000355944034615598616794768233311399622271384162498217941583401409899","1.3988518670956993473042950450484736714842318092867966950042899770792468295186910946769980020692397866684051182600894007"]},"10,1,2":{"tau":[0.633242610341318,1.4394365080219085],"highPrecisionTau":["0.63324261034131803452439795359178696062207361121778889103728046572182092025116846022649974933034878607052984158338064486","1.43943650802190837128391833419140856486287679008331875071178454356114772255661708691026316136923135044343923631314523071"]},"11,0,2":{"tau":[0.8257747128529203,1.584240849911943],"highPrecisionTau":["0.825774712852920257473297573687404492087572483915761614142165446976373112483973111648583319153019157221457841667942732287","1.58424084991194284022061416367396835045257765050606682796449305355710865907645218123424881892417939583534013806437891588"]},"10,3,0":{"tau":[0.17856285345892262,1.4252564030051202],"highPrecisionTau":["0.178562853458922616423093874764285600675467451463286197693309698444685634575584925430777129187917400775007002901222576054","1.42525640300512026508114119925383911718470053305696899809391311854984445185805634026191618902647644887975154487694464008"]},"10,2,1":{"tau":[0.4121503177019055,1.4474794764952492],"highPrecisionTau":["0.4121503177019055433006925106317310550694683605566856677911377617770511458089606937170163205746792179450987859302734806","1.4474794764952491599038779038067265229528887136361141172555118438895064178961957446305991669962247268664196817720246561"]},"12,0,1":{"tau":[0.7580585338981912,1.8403133497557587],"highPrecisionTau":["0.758058533898191270622712426209608345656942552850880836818009161433655742103723039804054952516634705721471099136574861326","1.84031334975575864818728266217204652300870930162508073210228440265815214074790908480025366040755147804127057011217202788"]},"12,1,0":{"tau":[0.32405874889385744,1.8705952969768187],"highPrecisionTau":["0.324058748893857414163424977968663392911562059735996042853100306292520787134777927685789283412342876614087218428658477473","1.87059529697681872672501196542483292482048162564298811253126005885396015056316859414382105072823166955259058174931497945"]},"11,1,1":{"tau":[0.5298773840192985,1.6256682262583737],"highPrecisionTau":["0.52987738401929854778238870179826720705227122534957569530964317414204199598670434775774754543485334146236368091369272823","1.62566822625837363083157588787285563155422063801106331381414578179407444902945208238459313215673513604488031642810770339"]},"13,0,0":{"tau":[0.5781666191366013,2.2602827047410177],"highPrecisionTau":["0.578166619136601329401649790549791030539851236846342380312852917588515785497406177363130882796687487291393972920337393149","2.26028270474101763651234278392280600945975291406013579838850254406698957566414345610254641034502516297522155956563152889"]},"11,2,0":{"tau":[0.22849499383738303,1.6137108449735116],"highPrecisionTau":["0.228494993837383022659506949022627204612702968343419601851653580726088929538655624955443465914014070034812284183690177877","1.61371084497351166214871921913024241868783575936958070290573445880364820150496076758232375213367890751472972325419428096"]},"0,0,14":{"tau":[0.9223893287870664,0.3825414515698568],"highPrecisionTau":["0.922389328787066393091823076833314742269015414474254455626270218438489835547521825918340495025449742063607964078231648329","0.382541451569856821433823838823616938605912524402766566627998327381822692035648953334019068733541411779713096612323317276"]},"0,1,13":{"tau":[0.8391550425654452,0.43927999078430247],"highPrecisionTau":["0.839155042565445197651583568104612309805686956594838956803205139188144378861557847927645801360456518633638989379447019117","0.439279990784302479200169241984995474580827788238105522263909568587349864620118792629324276557053810275528536651112530543"]},"0,2,12":{"tau":[0.7747058688568497,0.4753231578374572],"highPrecisionTau":["0.774705868856849697860025068464783241113111333556308401506588469229692269019859776779765971025617358124815394863871529789","0.475323157837457237234128689236852430152051146522015653882874831447820691632839683999689251540522941760949987904801819593"]},"0,3,11":{"tau":[0.7173337045611946,0.5006486337674284],"highPrecisionTau":["0.717333704561194577735814441091392674452995666412670783888290674896316962407165617147859510494679873939694034860364517148","0.500648633767428414279171655048824135226548977244631522117452213620083410179474629418311073247725262049614168829162161184"]},"0,4,10":{"tau":[0.6635457183303154,0.5187037485325875],"highPrecisionTau":["0.663545718330315342272120790266050950991033258432644953735032521115484907550121682112703573010768771850257745314488743216","0.518703748532587448920128251208279574101120356208936969102692853858001170471480712361710367328652385960078400688548422397"]},"0,8,6":{"tau":[0.4607911423987462,0.540627432657807],"highPrecisionTau":["0.46079114239874618006623113735691667191591748370279007513143187217171844864835481332176548488966065962131207312057640998","0.540627432657806961675667376556604159023142716922935675803153570385942171768869710055670290469191423649530666138139286046"]},"0,6,8":{"tau":[0.5611075555794764,0.5386385701699462],"highPrecisionTau":["0.561107555579476447287939911963862298097707348380587532843070818128460752817083240065277400617636638847073554725390279875","0.538638570169946209263878168812675109111249235639860596995125364221816340352707602325099211243368896355366252316382363201"]},"0,10,4":{"tau":[0.35838229941636246,0.5249330103052348],"highPrecisionTau":["0.358382299416362468557719485835939970694985107569914616824611248864565611270091854800668921355262272383312086959275445517","0.524933010305234814924067008567043952910845120620355544198073617925861042385307600876057570660488042041368357329244871857"]},"0,9,5":{"tau":[0.4101343776765566,0.5351295125497875],"highPrecisionTau":["0.410134377676556573101490467115283995960640883673640983725628791622244685366957590416123082127741164829507521127912546023","0.535129512549787440599352036027489570266128453143427345031414434290208927250186577411170697084408780060179453585868650944"]},"0,5,9":{"tau":[0.6117686894067903,0.5310901508562409],"highPrecisionTau":["0.611768689406790249570116239174994063243772183985511321957121967242999603967687885047575999676663261797175335345987518175","0.531090150856240935877468845493820707643872959503761645380157485293465373933856410327796414555415885406591594495408674016"]},"0,7,7":{"tau":[0.510949378125259,0.5417687670253204],"highPrecisionTau":["0.510949378125258935995424084371148946190291613970916714814496444445326913460187568135897964288687974153628719032358953528","0.541768767025320438429746552997373724999204820974754659339263856893718242384643590305287532931365205326626115010160332382"]},"0,11,3":{"tau":[0.3046779338502163,0.5093284861852875],"highPrecisionTau":["0.304677933850216273292903941489669330941711287129615986842596891913462740454302516767721178094851948889758180927707948871","0.509328486185287518227466299531604841950845375482203744180692626928822043124865820805438161732953632815600734544216889777"]},"0,12,2":{"tau":[0.24754027937121056,0.4869464278458852],"highPrecisionTau":["0.247540279371210572854217305754608463325744051228467942337790451102291462745306204222432596742560284922526350008450518045","0.486946427845885225473630386617982004809214669300960316664102885375738178139478453379576299005114065348975813876737899611"]},"0,13,1":{"tau":[0.18374377460961577,0.4549114940052373],"highPrecisionTau":["0.183743774609615770315940921277825115117366443911586486587737331804584610813812464785918100810180469882066442369201294994","0.454911494005237274695674498307361645497882696409776940911784206738280537488191017411968895325922587747728406285993056484"]},"1,0,13":{"tau":[0.933490553150132,0.488363463296833],"highPrecisionTau":["0.933490553150132004269046097026479413027664561082574088890286407280083749696153419555881185757084348940944614365536826707","0.488363463296832997413323770812127880169696209177822746376726492333988513834994625293920448506046396293215638355602280879"]},"1,1,12":{"tau":[0.8467610051409991,0.530271112724664],"highPrecisionTau":["0.846761005140999020248958060996254010949402617803015657637580476447997398256447320141045346823517099222106237314699311668","0.530271112724664029661925549468431744558032006780876184337114298665484286319209118115543719015587903920071351209629493524"]},"0,14,0":{"tau":[0.10273243236640553,0.40547873329466727],"highPrecisionTau":["0.102732432366405534248890304517661539731224128262596110783138582443343326492905312041603181529478714189396835708211013686","0.405478733294667281641667150443983342403315829458393505052705604744256537310918445375090753492649409244139887317904747599"]},"1,4,9":{"tau":[0.6527462375098245,0.5957617217211458],"highPrecisionTau":["0.652746237509824547518839866519140082147685458910556483784103810486456942260618731412272387002967949338327684255251788851","0.595761721721145809631510640533233600234648056473071456662364248175297420112333440764799498662314498480680578544264987834"]},"1,2,11":{"tau":[0.7762289325422875,0.5595259498620627],"highPrecisionTau":["0.776228932542287583839052637680865944108857441496391869488947790552908001091372885704364486036427666593600214939557593069","0.559525949862062769027546631602988821160517669477834523114389943827461383660498805305595164127845862872285391843581310713"]},"1,5,8":{"tau":[0.5948866722001753,0.6056643390492901],"highPrecisionTau":["0.594886672200175304162594215982699396928085842574799928142755383770106503520933909382374299383577557385625377985970544047","0.605664339049290045761188390594839511543791823700296762439642821731961026912032237179957484527231455155607436523566542522"]},"1,3,10":{"tau":[0.7126521353490264,0.5807131515259145],"highPrecisionTau":["0.712652135349026348829658327814592946959917120144745807070497823098482554644951400693080439375318337989746782754340992674","0.58071315152591456344746392475603400026756251349799300209901080659402368261810465489016222106000171923937239444302539058"]},"1,6,7":{"tau":[0.5380782861160247,0.6109501710378051],"highPrecisionTau":["0.538078286116024761202798031383030477163456588694949513737042489665967774166931476641949918785975358304213804463269711221","0.610950171037805154591543975815426019036584510180446626387680067013675551602153416843745339977131955666130312446210167149"]},"1,7,6":{"tau":[0.4815840391391055,0.6118622049267902],"highPrecisionTau":["0.481584039139105519205820619953610284765966309979819978061774429897602559095684917680084914538083923453290980172081388056","0.611862204926790239334045134242920650927028858209809845084466126612523437799820240711098986746047418114046573007159305565"]},"1,9,4":{"tau":[0.36681259044382547,0.6004458070775532],"highPrecisionTau":["0.366812590443825451976548340912366091317665093680507972621147513604695085442879509212467776201507700576629928848650546866","0.600445807077553276032758093592214107401414100815219383885147955420896401911512103506603312252981375761038223834878484873"]},"1,8,5":{"tau":[0.4247397670683974,0.608423995962291],"highPrecisionTau":["0.424739767068397414075863974844156328247447411691103099587190985773576246506210165804737578336123938711938947863375571503","0.608423995962291019340385849521729071491156295797995173191283258129425625187463561195904706736310635331386363264141972468"]},"1,10,3":{"tau":[0.30681815822350594,0.587474597545092],"highPrecisionTau":["0.306818158223505957641857000607356812558898562935455302581727953473997615144613707379543239418825649998356696883569973597","0.587474597545092043086975890832414540142173921872257687342071131558586040227108384737472967406942722902040399961870077141"]},"1,11,2":{"tau":[0.24315662556397058,0.5686542359103528],"highPrecisionTau":["0.243156625563970567096297206090832144390527771753853426047488769467687419706009165760481677899896955572461742222530014539","0.568654235910352745663182890841551342537651356984956409268229505525068771633114676321192823456660923849860832472081503126"]},"1,13,0":{"tau":[0.08601237830721474,0.5053173994821852],"highPrecisionTau":["0.0860123783072147381728496869846792344491122851390760710648454749414886368174017936553234779772538939637333141188923484228","0.505317399482185189683634022760060305969068719139006819506288096919256112370174225265409441480656335677347248467554864674"]},"2,0,12":{"tau":[0.9361653493365252,0.5740454898864948],"highPrecisionTau":["0.936165349336525198912619029712655236986206742066822955294126421530619130700912121242947593045183762477483535706898453033","0.574045489886494801490066176705559774048775698873499716051137169248730736632573887798702932395066554908828069300883721487"]},"1,12,1":{"tau":[0.17260137570975198,0.5423777434023371],"highPrecisionTau":["0.172601375709751970992088111474529418958697896651698785504680799603689035896469284513293525021996420710610473634189338678","0.542377743402337077994117861975585148588480062977073466741508858533199827368740326252677897961344513945489329765243496005"]},"2,1,11":{"tau":[0.8463868411018071,0.6100173401104462],"highPrecisionTau":["0.846386841101807033096965350194135394231037217478798404267222072868519523579792785911505044588638366109583242466962564224","0.610017340110446175193290108164485584950608127556887368426825063434370386083539960792015145064064431892324410806129841469"]},"2,2,10":{"tau":[0.7708712102473464,0.6357152874710401],"highPrecisionTau":["0.770871210247346402271132496038506713876878992636260000949858547550516907835906973412816753229838509202229848202850145118","0.635715287471040067255595979879990755969701279374427781915024385915828248865727280696336251410904040140114611305898912208"]},"2,3,9":{"tau":[0.7018258777680003,0.654317180507618],"highPrecisionTau":["0.701825877768000353311107982976547467816832581477646883872467365027552175885073962255174211571378769425069022283795551258","0.654317180507618011188987423151619658840513187035413447355751999479563990890050789026950911319712595785220970093767227516"]},"2,4,8":{"tau":[0.6362398818713875,0.6671895876876294],"highPrecisionTau":["0.636239881871387432276786224000616937357107975862709842266400048415432666733637575873865272845318923536585338680206186929","0.667189587687629414704465063753229078858504493272442610202779397595302631655713586229979756510407033321705758017760241162"]},"2,6,6":{"tau":[0.5096098981883405,0.6782592296312535],"highPrecisionTau":["0.509609898188340443917287345633150842348860071893026133652963267211301572740057172089982234975443529598771184381701137776","0.678259229631253468802678943403047786599118770756821236982604532489052220690035883739669036028801158639032878974452793244"]},"2,7,5":{"tau":[0.4466661835321044,0.6769508299974706],"highPrecisionTau":["0.446666183532104411116971612296402127721627041544252050896344815350644840995372436084858211703821969213029446060557674144","0.676950829997470539795141291409935587278510081477736399851390575330534753128612550710562030739343595788706595970631962557"]},"2,5,7":{"tau":[0.572518051050392,0.6750481363244598],"highPrecisionTau":["0.572518051050392023205159730284042768806613482457802994314884092098021474500613135502722969832334747110636652123799640535","0.675048136324459852793768734382255244913646911709249458992371136609558861198896144567369175048946203410881579291685265963"]},"2,8,4":{"tau":[0.3828372528287884,0.6710463429432526],"highPrecisionTau":["0.382837252828788382206310699116996017524928257272806359084607829682031398267086668190772095225538415325377933303958681023","0.671046342943252564135300351702483031978395882924640352956453285175269336278621466306930414911757079082494847447681624577"]},"2,11,1":{"tau":[0.17190435495254627,0.6210937774663854],"highPrecisionTau":["0.17190435495254626036343057090816849132950023122082026533073343586601004356228354110805805426852521860850657347778668809","0.621093777466385336098248236032682741907597059065762578514962280059131028977415161300872577553281049933818134577817420209"]},"2,10,2":{"tau":[0.24776911278772046,0.6439581501875191],"highPrecisionTau":["0.247769112787720466290249301180518093303805827215054209107097208455258219872506057070000139211409479226203152016790339242","0.643958150187519139702316339705275744400771794407947990074501044313912454968617309144751310326501567457856769048553395507"]},"2,9,3":{"tau":[0.3170711985172672,0.6602476273067223],"highPrecisionTau":["0.317071198517267222790289594005776134883224139397273118857794752900733145304670290106618915928011520369361966675946359121","0.660247627306722227619258321916442990924196878122939636721520468352271188020526751672081453937821812882522541995439827293"]},"3,0,11":{"tau":[0.93589984417963,0.6528583170828339],"highPrecisionTau":["0.935899844179629986644240906141502062883563705118537602817795348427572204035577317467801785378185477638226668631568520754","0.652858317082833961339830557328640250900797587885580901437798008027106583065886003536897244416509292040531634727753396902"]},"2,12,0":{"tau":[0.08153636280486469,0.5896044042626103],"highPrecisionTau":["0.081536362804864683192382463571220765963097871745698287989905000610532642538854702242060193794710355254991588168347430378","0.589604404262610328840393420327542249699966909807097402345455557951731772981124293843701433309805192090011713617511313461"]},"3,1,10":{"tau":[0.8419348919130294,0.6859592435855795],"highPrecisionTau":["0.841934891913029418310277389471082450524535272467022703450288795892656336628908317012463261352803099268117223982168654357","0.685959243585579437519208956287719832269224076140782147234888612661900854260045546305317811242356519949500568087625115646"]},"3,2,9":{"tau":[0.7610639969557145,0.7095220949889205],"highPrecisionTau":["0.761063996955714523880713517928226700433356955636517692849269069335347960277324253268158106804738572512124590581420942554","0.709522094988920539075878469727526314098181636210182773495692432130283084214973546266198052907704404540217532834025469277"]},"3,3,8":{"tau":[0.686233761303362,0.7262628596655234],"highPrecisionTau":["0.68623376130336199558867138728803244494327949983407241597004721095449002495749483790288443436496155984417774629945806465","0.726262859665523406258629448953498521987599430174026695412357876672328859175323434829400455335846920158495154516093154374"]},"3,4,7":{"tau":[0.6145881887317797,0.7373053158713555],"highPrecisionTau":["0.614588188731779659669961517787315105226867234777469110454132655106829260878720728845759530915723439557145868270638793601","0.737305315871355456655292763911147069316643400054142110994934346220968362438345779089845721941151309993940996283848569526"]},"3,5,6":{"tau":[0.5445275088282272,0.7432181866271056],"highPrecisionTau":["0.544527508828227171258478936331873795150169254633797950671921053305398050716067119947823762348353567215818848848303074045","0.743218186627105608003877726106883167262040671039147523602104398612370089020599485570446783753821899640156454814875321621"]},"3,6,5":{"tau":[0.4749171060055317,0.7442603975088518],"highPrecisionTau":["0.474917106005531694069933099090511727164313582143089573484811617886916474042433384677194003066339573716360905270505304367","0.744260397508851795351387780399524099178244395461052460282365995953604108961854433684078109361033287203288333959758716201"]},"3,10,1":{"tau":[0.1762481419800606,0.6967042387735546],"highPrecisionTau":["0.176248141980060624536358230215613438713399679248029113987064345644260658096567123329867513540005244771458961443508285591","0.696704238773554580993988586297132640215936670979810116455308534116981287268260268693465838277310014271035746316766966671"]},"3,7,4":{"tau":[0.40474469067732316,0.7404620292602487],"highPrecisionTau":["0.404744690677323142408340495865455557110217311344192750723875620527336385293185401255961924103902238848538154646691943711","0.740462029260248775839583348277560766745296919313820042218407438331615072131775110760863911931319915956011159397973341859"]},"3,11,0":{"tau":[0.08129262339280045,0.6681409362851168],"highPrecisionTau":["0.081292623392800449945348534309769241359762495544508446445030953048175765639748868211467873385216939585708641466222062885","0.66814093628511678859264402668712899013942435863727599613819899657797093003653676857055775771755309969198400305846298003"]},"3,8,3":{"tau":[0.33286896997498544,0.731635741310402],"highPrecisionTau":["0.332868969974985458058345892805913484027242401993206319345314778188604535529964364142321914791529496579444094403317485109","0.731635741310401941596566084966185083814060498675375264787036313343345457325421885877088054348242221077742468797667326768"]},"3,9,2":{"tau":[0.2576722189848771,0.7173316743092328],"highPrecisionTau":["0.257672218984877102406413568926058387094709122190341508780518196621626040627700902359897926538511031698643366543301197971","0.717331674309232788768388799833976787986403145886729532529426032092782496147783105987947408479844205775364289819564165326"]},"4,0,10":{"tau":[0.9339278998715803,0.7297611436294077],"highPrecisionTau":["0.933927899871580355886374305622309635921023063224463654894547750814782672798620943498501975721900785520908291361598734134","0.729761143629407673319733070332504058565981275946216509825809688848517017101028753338732996332462575467853900169624477909"]},"4,1,9":{"tau":[0.8344389956685113,0.7615015185051726],"highPrecisionTau":["0.834438995668511363765393273072407472892244320474559519761933402739942935837031413475478836116167446844231040817526382928","0.761501518505172609660966077500383527354554298747164410353088211834010287519364185387420162964885554048819379008398303437"]},"4,2,8":{"tau":[0.7473495958881731,0.7836882660697881],"highPrecisionTau":["0.747349595888173052742973949734594575594266688641315629587860173109502100917855798736479923071937795525244099061565376384","0.783688266069788135228294701371455565728596822128080033456296458745179636512889609605069981158595663551370595519649623412"]},"4,4,6":{"tau":[0.5874072147672713,0.8082182028275158],"highPrecisionTau":["0.587407214767271355227993802870631375193733305625542755013517875175862015956591967850241185992497274160708229793844022649","0.808218202827515768497204861233294695710968752390581404931815641209103816740113284561207868286513472094891919082487004012"]},"4,5,5":{"tau":[0.5100776589420866,0.812056909282993],"highPrecisionTau":["0.51007765894208651498834897854481756810099105259624392047456604399018899204538562713729124531560773526048287068533111401","0.81205690928299297258923732543877453326490905424457179477548660729779337230132974688904203733391987831648252831125230423"]},"4,3,7":{"tau":[0.6659429535510664,0.7989213403939668],"highPrecisionTau":["0.665942953551066348955224884911607304301114905040821172089709053855375119337779050839784790750107195947587673810752345568","0.798921340393966798454837666835083196359239725874751675022727732536084174062921825999359699805406102270245021828160152571"]},"4,6,4":{"tau":[0.43266736415643564,0.8106043004124315],"highPrecisionTau":["0.432667364156435625533504554768761244861331253172496629123279832762194243684388434704089157031764663382893476252990888748","0.810604300412431454407930606610304637646446871347569417193993828258731316510488932539449864650899397244100420964997014207"]},"4,7,3":{"tau":[0.35388246422648717,0.8037774767618773],"highPrecisionTau":["0.353882464226487155229560560163025133426307256880499421817981815674403359785741329825733979971482140109905378405094227973","0.803777476761877325625642421382264216831651843226642458559322729940449995355004196489522097985490831716528568149389785027"]},"4,8,2":{"tau":[0.2720324065608931,0.7912255473108809],"highPrecisionTau":["0.272032406560893092115139193630086047358863963341553270062948664506452743165998578810656500146095219692873869467459899424","0.791225547310880907821444998211235953619839273259543516130825561526224462437570886184994304140752893982878240580011820928"]},"4,9,1":{"tau":[0.18423313200767433,0.7722182459851897],"highPrecisionTau":["0.184233132007674337913625831371781453432745887695303207446846609195243823187995374973356760100812582960886908153186344884","0.772218245985189745820768108487581626406572192017276755609252391812792121651659841959214613556813662101277060039180345681"]},"5,1,8":{"tau":[0.8240450211031198,0.8389190477603129],"highPrecisionTau":["0.82404502110311978423923867735952453412649579158285505481809031678013542145108746279535722972986210681995796068609525547","0.838919047760312952346796343421592310855862496145457528967778638942303147993041418249786660828488603880920936053109630329"]},"4,10,0":{"tau":[0.08344983549141918,0.7452751156882867],"highPrecisionTau":["0.0834498354914191877330378952101541776400401315943159681462575580874722800420048230368936154322954667412967374906095773361","0.745275115688286759525433397731166543297331688281846948383036934766922893445581434367141865963501908254575367862348599347"]},"5,0,9":{"tau":[0.9305914123031541,0.8076112380501875],"highPrecisionTau":["0.930591412303154077341209773362382181329147969657943285787573386106894646780655408189860822009008527461788995582292277478","0.807611238050187515562813906765831493128902338283091743290016849800147091395164758821265685740108555081174041573864661767"]},"5,2,7":{"tau":[0.7295046610508585,0.8601669497943283],"highPrecisionTau":["0.729504661050858583115810534363789560114962324866687989541955453868661437198483057827923942691591517350127895713585146277","0.860166949794328247815796802219515697615365565086699446528896586780466170060511177827657631235669119875714291201738375407"]},"5,4,5":{"tau":[0.5536312593467271,0.8814287789575078],"highPrecisionTau":["0.553631259346727081304460452979919067431151292244816942278643069165101670331480750014944127727517833718456535650521246778","0.881428778957507831503771791538087453625142887701674841390546176936488067394762645752387802547881497884593761383106932042"]},"5,3,6":{"tau":[0.6403243988160838,0.8740110694775064],"highPrecisionTau":["0.640324398816083841540217557316221592042437370034785762603082723169383634746156555853042565394985906793485047623534908954","0.874011069477506412995203876008948270133663406013053549659843523360144366160554500353576499044509072503663092164421297636"]},"5,5,4":{"tau":[0.46760706523483453,0.882824646867966],"highPrecisionTau":["0.467607065234834553328207460328968773481578055175011663208202727082162689358063692902544740488243943047591810569628475431","0.882824646867965973255654797801012846388254207221005869330993365765648770967861269470762503119648702326977916404049652453"]},"5,6,3":{"tau":[0.38068036780464204,0.8782546544324845],"highPrecisionTau":["0.380680367804642041485786937656695327244580599080182024687070686368805525766137695264650430646637884871133622679872729207","0.878254654432484518438639517084826346592022787531023900098530849363734900691197655602594103698879734210839049689127095646"]},"5,7,2":{"tau":[0.2910039141828122,0.8674591842631038],"highPrecisionTau":["0.291003914182812188447367661193084965541218586453239846100028913180249869446791526812088440650955317503189856224332249066","0.867459184263103794868464188265748769629053951506253917394671917553798583232876508229328913437922530272453993581570870349"]},"6,0,8":{"tau":[0.9258927401109209,0.888631289177905],"highPrecisionTau":["0.925892740110920945815147547781812180060057358119012209332130833753280533640831335572433517295029315058223900423721004852","0.888631289177905025071900645879399634091936852328456219910695560282589740628535358280019423737158774517933609069519696752"]},"5,9,0":{"tau":[0.0874883351077125,0.8236876488581545],"highPrecisionTau":["0.087488335107712505787535396017799982646479111886149476689139812464284128140947223959107809864596063742693636629188254393","0.823687648858154514117998213635605430455885690319861360508167960551031387444782976766348910803034773061687256956329457923"]},"5,8,1":{"tau":[0.1956151534094621,0.8497673392847127],"highPrecisionTau":["0.195615153409462079138574025236306900248781844002039892925357144296420193455513696350467121924160562215053114764023821103","0.849767339284712749923425607448295439726221945262674538048261619646354471826035402349291875638607167430923399042727229155"]},"6,1,7":{"tau":[0.8104189295491027,0.9201829648159305],"highPrecisionTau":["0.810418929549102607791944748888613913750707484033145011601600539442025991600802257882741327457015997848553715228989898994","0.920182964815930560791922961544979672971059917867443683424636815073405608273044899320010488094342583734502723786982731199"]},"6,2,6":{"tau":[0.7067533873778477,0.940728120126954],"highPrecisionTau":["0.70675338737784766709278795552792923625820186243812755395165869313303824436399724674388387041964516627693096823760323266","0.940728120126953931703213677201419153225558058704569442351886197169296513416147633851257702094458620577484009123157236595"]},"6,5,3":{"tau":[0.41456038961295893,0.9564863814184388],"highPrecisionTau":["0.414560389612958928449031484372748191799374850509968866240055831252848178406147507455068663715121687462685401432101105497","0.956486381418438784501967258290835020123720143648939880871707077471639082714408697279190380690620000923591452452469082778"]},"6,3,5":{"tau":[0.6081049183991147,0.953092022270483],"highPrecisionTau":["0.608104918399114713712168202927899015844595609034061988481672753438158421534358007531551637520865770294439691682188380793","0.95309202227048291444660856641307194764639011740910650137621144532793130113481599291474473245061056820539464389548736119"]},"6,4,4":{"tau":[0.5114153093144821,0.9582409437676357],"highPrecisionTau":["0.511415309314482047116723111780316608705836266591623156755700900712378635939271992646736499839311786559203330811289108972","0.958240943767635725156597585173003338236874518956618376272877809802199346362858254618503761529571625683400358774280980715"]},"6,8,0":{"tau":[0.09338170805448154,0.9055491342191027],"highPrecisionTau":["0.0933817080544815319343715670795194007725970040946327886120401129739566849111814278368876327380427035824354539537100581082","0.905549134219102630928680774130490108752710254811740733921537118189815550192175721636524061570931317763957583669954201699"]},"6,7,1":{"tau":[0.21074468372366348,0.9312469907433482],"highPrecisionTau":["0.210744683723663487415768172918059052642975971788778234466971635488505542694960082593233809213351127971362246489777242294","0.931246990743348168322262883456342317535996805765918310554996997079926498557989553770055517212253955612769863273037042011"]},"6,6,2":{"tau":[0.3153911252757734,0.9476974474707535],"highPrecisionTau":["0.31539112527577339274525651951641265603307265153884515813921673976785155040094151604328810971773785308922769512293110918","0.947697447470753424246241747546985597863849528676136866861035515501664159957175575396516250096790185075691670998078755126"]},"7,0,7":{"tau":[0.9196141938186791,0.9750019229729491],"highPrecisionTau":["0.91961419381867909484102066637770039941150709967562113796508203597383983027750837579302859701376066102621052494442254787","0.975001922972949064513534946546940200580565324659512302265148881205107081825330534313469640789972763426289231635689976609"]},"7,1,6":{"tau":[0.7927879005413166,1.0073650739104931],"highPrecisionTau":["0.792787900541316693668096949546696717503304881313558495301086834664062947450303715281146233057211589725600221788108321384","1.00736507391049316326504932155054262945398470306721945317958479855348924281599447629009238415898992337569168329270056673"]},"7,3,4":{"tau":[0.5671487172535868,1.0377962971051948],"highPrecisionTau":["0.567148717253586818422790453636895879066563602300601191366208285733512255847910665635095153134317889899817331267789512069","1.03779629710519473726314327296640928214991038851799580305115986824040467688475536120091738690589942750944171461847456216"]},"7,6,1":{"tau":[0.23049574913064075,1.0186640417440136],"highPrecisionTau":["0.230495749130640756701820212394506367446390198453671391703706038618381984872841329124156404624384125019423740281490984433","1.01866404174401351098280252571954173966618620343543929690925321033130616239398091845645894132418764022427585351502133516"]},"7,2,5":{"tau":[0.6776862364272683,1.0272677664661283],"highPrecisionTau":["0.677686236427268347408708673572821312581127414439116734372869962685204362618199869284587252270534037909801638657252280739","1.02726776646612836784562186660524305990921659950425129203368866698657278473753188556703798203764653921992654837112656965"]},"7,4,3":{"tau":[0.45774089882505925,1.039898107393273],"highPrecisionTau":["0.457740898825059273251416992004039008797767895887405016575441403396471649741894936503259562294602940739666069052799738929","1.03989810739327301007460432886364602040136839055096272455864384628821642941920579685831551504533443236279500325917110225"]},"7,7,0":{"tau":[0.10141787809401678,0.9930416950277594],"highPrecisionTau":["0.101417878094016785708495803609618243146734928370446575418825098353789382259942196977007339081047273427394537056500049798","0.993041695027759424396322145420491472515517085961320565072843292253123383622906347113193978087347746050421079588435933471"]},"7,5,2":{"tau":[0.34670631979453,1.0337045077349707],"highPrecisionTau":["0.346706319794529993592259413926264095847700210468121517162655314440665521772465203097945827860152976783672415972112752643","1.03370450773497076960616307561478685559523177935033506798578694806346179362759146157437442048895138680627644837884142069"]},"8,0,6":{"tau":[0.9113016149775144,1.0692654846435605],"highPrecisionTau":["0.911301614977514360774657584271820181014700174706066178638231955981730073221974123112414423559092536901365974266461874809","1.06926548464356051200427578940127141208429507331754446032470137807956114446855385203377653153285201098083219382453997731"]},"8,2,4":{"tau":[0.6399353545836429,1.1220622796881707],"highPrecisionTau":["0.639935354583642883489147888942886127033327691768989813083486294152055671904030642859649290828456612519526891922557131224","1.12206227968817075312958487297062346298891019285793418503561045756892645367337811521451505808840465562861826727634128655"]},"8,1,5":{"tau":[0.7697897747957696,1.1029694353315909],"highPrecisionTau":["0.769789774795769581498981887208322109739112882188136657633250533138346567379270870554374557947184099978602078467602923425","1.10296943533159083726301023764283303664200878542553931029076904118565498989934546195588970343247844453108932523181346379"]},"8,6,0":{"tau":[0.11221980890603077,1.0887427839632644],"highPrecisionTau":["0.112219808906030769052368249914435469562569882791535878999405749134965656452136264638514363790638751847648001988580277454","1.0887427839632644615591157585265416248489198393631464537822754582317468269513764393531118674307126569314427691488093694"]},"8,4,2":{"tau":[0.38750542281261113,1.127536884602037],"highPrecisionTau":["0.387505422812611114035321772590822456083939044301089904692517171380526705681106678277663906195539541685001010023088615005","1.12753688460203714691407922211645604192584926939060553672518663731862830385586245153087793713787680815105179113133467508"]},"8,3,3":{"tau":[0.513906062764676,1.1299696643561783],"highPrecisionTau":["0.513906062764675949939589257772798235198317218100123659385114121394323446208305416418188545637124775976606433541551454279","1.12996966435617837367233354563393237096323974173304268507964255385410594478417579957005329215920418095153593600171160493"]},"9,0,5":{"tau":[0.9001441992649009,1.1747878349815115],"highPrecisionTau":["0.900144199264900903440834272882239868579427911945680594865083516005527847729972507872481877955435781981352468538689334442","1.17478783498151157343446361827167131813149328981980272595035131392468835001485517892015804261877342579027588913160715809"]},"8,5,1":{"tau":[0.2564319332140167,1.1144299658538943],"highPrecisionTau":["0.256431933214016711962861253714671223914329657031567441776979237972849720424793975656242876324066199792080051264833750003","1.11442996585389427057892505182953429050461728986510724553621011546918955195625932644640630609358881750749952241023905427"]},"9,1,4":{"tau":[0.7390768986909234,1.2103507080811866],"highPrecisionTau":["0.73907689869092337512109953586747716627304095444736905780323414847373916164032414423034758345367964165026847534518629502","1.21035070808118660108043001024527459097653947019901217457383299419262712925403809636581867193964251376121617319771890924"]},"9,2,3":{"tau":[0.5894325357815896,1.2280599004974975],"highPrecisionTau":["0.58943253578158961588843543290836259292563929953104424459504538406550024532734680227114641891892322225476076098040884892","1.22805990049749754707949547257034333848509185325689847526633344250237219671440814325449699962742261170700999042459495066"]},"9,3,2":{"tau":[0.442165216211689,1.231711291750181],"highPrecisionTau":["0.442165216211689031559688785375500739893365967440171009609813383437835908945124528691937004682727028613510263678937653928","1.23171129175018101035599864381636827591914271814152922258337558014580206157408173142088877175691513527098627358767546132"]},"9,4,1":{"tau":[0.2912596251986874,1.2217184983917744],"highPrecisionTau":["0.291259625198687449031589267865320724633922038043376709264478015584948438879723258264856864324453357797129491133117372925","1.22171849839177440294803164341110694646670832218672781288456435871957700038791103067896700777982138050661239881409170867"]},"10,1,3":{"tau":[0.6963984019538398,1.3343961655483831],"highPrecisionTau":["0.696398401953839789454036647985844210200976277721082395576042985731879738848935730023954322157120589033023477488034960217","1.33439616554838320464265002257192301977575266077257636076123417236006956411111874719428709458584786406560576519295527527"]},"10,2,2":{"tau":[0.5186709440485201,1.3492424889953707],"highPrecisionTau":["0.518670944048520069479927540579862102332640124294065433585287513615205381188112873995758259223603062048548812118273102573","1.34924248899537074295994225923432953948017297951961169519979318074410166668892307605588202958594310670998387190597224772"]},"9,5,0":{"tau":[0.12691741521545255,1.196081350537873],"highPrecisionTau":["0.126917415215452543067254985191810982824163502960621363760380910608650055132421666168093592288769338666079441573887106053","1.19608135053787299900198528687758684964376402073983670475146609037925446825189562925801878100151344317472209338482639371"]},"10,3,1":{"tau":[0.33987563682215755,1.3449973263843118],"highPrecisionTau":["0.339875636822157551714194449964649944666081731990886783719937771929597440601326813540717808114918278186148802485267212799","1.34499732638431188159736235952596207433798620789658389011798258223075703968456688193655829434209284553909115527736756582"]},"10,0,4":{"tau":[0.884676959501032,1.2965089018086715],"highPrecisionTau":["0.884676959501031989621917735821667754364367662530930455084174473611752682266548227573437890968605292998081273754872029224","1.29650890180867157834585513194995534205573960249262300442357526729946271145378270651573474976643405308165825245152040185"]},"11,0,3":{"tau":[0.8620403177434205,1.4424468169562077],"highPrecisionTau":["0.862040317743420491947341567095992827132346123512452799628736915541759493716208070512603092207636546665166847208760782428","1.44244681695620779146349139519003452476523858138803639735871564953139253391723916030932217579738704580173452116481828594"]},"10,4,0":{"tau":[0.14757884989708808,1.3200812942415847],"highPrecisionTau":["0.147578849897088070064174747585894210622953135610006853393904290340315496337958455140709181300434575412203835281211962697","1.32008129424158470620355522896728866206168162928597619557578673019892819380895552510109177147828364280515429269737947496"]},"11,3,0":{"tau":[0.178294215989465,1.4688202687762588],"highPrecisionTau":["0.178294215989464983072725089781620155477552876887625288371683985219465016336066196044850844220416033675217463401614907402","1.4688202687762587136212219094926643070418254669847019148115236071610148998519326056286180301360734795032717492872606493"]},"11,2,1":{"tau":[0.4119847381542441,1.4908371683483794],"highPrecisionTau":["0.411984738154244090047453370388631283989329917819571737285126780503349202439489639823086224888313791928469881756757437195","1.49083716834837939113361321860374417906671889441555252062239152369874814060448993090481311198562557836239303969253839895"]},"12,0,2":{"tau":[0.8258402078901242,1.6272849710709034],"highPrecisionTau":["0.825840207890124297095806393998863119061832969363003163998526442436977175564210242590270425522139527923134086109871160048","1.62728497107090333133613590758904821429422730169017976536459915839741518223378514494135015252154471508528129577810794065"]},"11,1,2":{"tau":[0.633270216568141,1.4828018938318066],"highPrecisionTau":["0.633270216568140947485374636224942696776617341122211709215441830376966546342695051555541068739208955536253629135557778817","1.48280189383180654402825673621109114645987073713565690214260627636986377026528854589791471155440378161666019870880456178"]},"12,1,1":{"tau":[0.5298275280043219,1.6685973140147004],"highPrecisionTau":["0.529827528004321839073809450306899837871081088573238234693368250880512322218931226289035511053060387220246810806354247032","1.66859731401470030962961043284209934289852001081228466178900744462004620878145879780881122906315354416499362672778267853"]},"13,0,1":{"tau":[0.7580748164486868,1.8829917403937926],"highPrecisionTau":["0.758074816448686804843524937364656509780286919310379634761635803688771599997961917202218488949770762140336126487971884924","1.88299174039379251555738402007236155357610406771487965383510013205263655241905118331802991333446762583136714924764891288"]},"13,1,0":{"tau":[0.3239780265760303,1.913267117962907],"highPrecisionTau":["0.323978026576030290476232650160476445434453296599650709502605957731666531479439348156822600670762474969034835260618128678","1.91326711796290699113946567019856108939236999208304121960648385529816478494157015216796238067826989758463066071764235743"]},"12,2,0":{"tau":[0.22832981204295763,1.656738084167294],"highPrecisionTau":["0.228329812042957624334138566698715125542859783538147035919275827470745855699344165446943358261697094357854431557084842909","1.65673808416729403650286469296461785514923137924707602491537413354042432707928901483432708879110251672755870385451203894"]},"14,0,0":{"tau":[0.5781508413290928,2.302766732009113],"highPrecisionTau":["0.578150841329092805509984107276979831464844945449302920214851988355188719412905324041414621116902920683030032642470589467","2.30276673200911297290414051151721260650984673038862484624241079451548635720576634213269637960243569376215294313479324672"]},"0,2,13":{"tau":[0.7830991929871853,0.468313489599611],"highPrecisionTau":["0.783099192987185281819717747730331168209558719489382493809957146919851017963002553652553952072104859558768850925568971355","0.468313489599611021087113060570101366080209559412570329356770193689921017922672060937020945064922143801832518870166179991"]},"0,1,14":{"tau":[0.8446638499671623,0.4326753058003434],"highPrecisionTau":["0.844663849967162221997143185989858886234684552934645492369540452875335665797140176401261675475897478124690223641746589908","0.432675305800343392996355277612437293855780774483528619366452891921331605266974554016826405278244883708003591467007728503"]},"0,0,15":{"tau":[0.9247025387002802,0.37702898116058414],"highPrecisionTau":["0.924702538700280204837473380705723534514495705140129565956473528697971560439708945483804184408986256204601589310638488156","0.377028981160584164787706607394390244460259610853819403845922817205915192885908457233330313949095376242974480306735529288"]},"0,3,12":{"tau":[0.7285312939411894,0.4936820001321328],"highPrecisionTau":["0.728531293941189420265691866999916378082495612161215672472580974482882844528202632049202596566833122818844622242350741449","0.49368200013213279071124836060325419805005169332574637040277334660158061484945100268480082016459842372446043888316626419"]},"0,5,10":{"tau":[0.6286221249164963,0.5253223061266369],"highPrecisionTau":["0.628622124916496252496517747776950690066615997151628214686401580687267001044470294657005847785406422029914088232532671573","0.525322306126636922135660598742253726766956338920941891679886675938130092787815991633640605229216730578049737774793390948"]},"0,4,11":{"tau":[0.677547493836691,0.5121571824814954],"highPrecisionTau":["0.677547493836691046903647293072812409742368280655016946404735051015697201109281271481625523437630607655582882479709503817","0.512157182481495479633787287789644287067435535597769478841506859428665133388266893775800164448644741585639671759878111136"]},"0,6,9":{"tau":[0.5809016497397761,0.5340188557597028],"highPrecisionTau":["0.580901649739776054289617149093275487207443715201739854623752765583213387474625029784235246679445518403048934822751032271","0.534018855759702804268531775859866866885960513799489416574028166519458320784408131478087056098729028623086051982277870489"]},"0,7,8":{"tau":[0.533820682337016,0.5387022771544175],"highPrecisionTau":["0.533820682337015938017493357925016675095589000604090581328466537150506005470217254677241384583300512609169588117483592047","0.538702277154417525927056773487718159108533031535172202713524497167185472796824049582311946269367432325072686272234840474"]},"0,9,6":{"tau":[0.43986220961190897,0.5366893156926101],"highPrecisionTau":["0.439862209611908984828994461823127493503192083982976996652533309301505254138175023709155815847519012539834606830641276958","0.536689315692610142036562310962032484676540988220523195592079209361223594677938320835099554008657915315813133552163491587"]},"0,10,5":{"tau":[0.3921536273717254,0.529857241178694],"highPrecisionTau":["0.392153627371725399626715023885403407347171426645221380825539719537901198126958422928928326113313224022653154667327569776","0.52985724117869401367324697627173220393152240961565511560597941037941841759568390571374281599540093223284614086941270837"]},"0,8,7":{"tau":[0.48694084117886965,0.539584533603839],"highPrecisionTau":["0.486940841178869675506059502798331283639009115681304905948490870767854492153141561445544434526970952451398140058576760693","0.539584533603838920215226904571135891651025385303497109624680133585282214328324148372050916618240833178850403711444361834"]},"0,11,4":{"tau":[0.3432677934518387,0.5187039743999083],"highPrecisionTau":["0.343267793451838695999155091178490989147978769271728352881718116143552317915363491966933129384678697385546647449821108908","0.518703974399908352751835496664136482908769164062860857712038208985741092141133235765259247054036303318754609024998300433"]},"0,13,2":{"tau":[0.23809225361437555,0.47989593234699796],"highPrecisionTau":["0.238092253614375551902606018741565341829485516826117201505931868733569344202431374777987200654680616807147282930761469693","0.479895932346997934650396424891389635950694971041340579937390641330890290546983605826757021838841186468685427636812047241"]},"0,12,3":{"tau":[0.2923910697053928,0.5025032066406875],"highPrecisionTau":["0.292391069705392797877935986853812208315714863852721969326472867562774942842062522782047078837714878012969832293966979358","0.502503206640687482892850326666589251064770585104282533882955744739359537979011364671194314344004946076628865465578115667"]},"0,14,1":{"tau":[0.17722968965291586,0.44806201713539157],"highPrecisionTau":["0.177229689652915844299954572143866966502867868741897415302487975954242942120515155740657605235843805744443827579378861656","0.448062017135391546234192188999646896457032925995471148542802583263824057888443340829952786522903196835164037128546280717"]},"1,3,11":{"tau":[0.72588953393176,0.5702963153623049],"highPrecisionTau":["0.725889533931760063885140748210944242022641099612410357764319546739025316459810794415097691069813968787554371231418313668","0.570296315362304899377526212458139913888098524579246150872859669486895170215374291637331292483203169191420742168243870863"]},"1,1,13":{"tau":[0.8530997408239182,0.5200860210144161],"highPrecisionTau":["0.853099740823918187474807405289282916049770748344319679423578289849886341692739718001111990203356188648199471298808043965","0.520086021014416171358781552825230804770591890246747803284478570908882127353507996168904791466587863127485375882711702358"]},"1,0,14":{"tau":[0.9360203773244657,0.4791766596638121],"highPrecisionTau":["0.936020377324465776716000303742092704770263115560193063031050687984651381818184077411018477359090891968194445343919013307","0.479176659663812104168251144906889653285379505732917125085434096704494141022358612903129560176872820935661415152479652835"]},"1,2,12":{"tau":[0.7860565134189939,0.548991520342759],"highPrecisionTau":["0.78605651341899386288836524620257874975017651070300305693837093843685278504222910740277408791973493660628241332995323896","0.54899152034275900673109529964634413767748786364986627717256237132125864847748078446654384666042400690359520572858691757"]},"0,15,0":{"tau":[0.0994894843846173,0.39943677330367994],"highPrecisionTau":["0.0994894843846173027061754868070224389648150291285112401864516568440012260497864801030926718290727897636933138457392200661","0.399436773303679946438749534408504557649726303437326580805032842578595756376958388349527506840090702519876156942114313648"]},"1,4,10":{"tau":[0.6694110786390894,0.5858727286987461],"highPrecisionTau":["0.669411078639089399173289361395062851342407866700711529605879552648169248575421818072312774487193058515031346152645638595","0.585872728698746047517993415261085301607395762909651920887326806512877250455427916888580787595047135731303494124958157073"]},"1,7,7":{"tau":[0.5092920691085489,0.6060456479318813],"highPrecisionTau":["0.509292069108548833215117840672846272448856461935549050104451221033413241550684817523564460820337637363298828372471834032","0.606045647931881291816703291186362176737981170490643947024208292642415341886844925497087349887141953338052814172065433584"]},"1,6,8":{"tau":[0.5619080069175524,0.6033308510638451],"highPrecisionTau":["0.561908006917552348875529423297171601140117940594634509502495656660806630901868198104205213493930359052573659438325993411","0.603330851063845122255227837050331562738018311759699533799922621375745269541279670455301588855568064954168777257927527551"]},"1,5,9":{"tau":[0.6150613397886892,0.5967021352808437],"highPrecisionTau":["0.615061339788689192051766786150410373702377599587097531665697844839287401032913915707392438287842516377697006619053643181","0.596702135280843705164110594081899898020075494154185706544223195322545950804362588730834641925073659136036166970094089901"]},"1,8,6":{"tau":[0.45666247918286623,0.604947097475667],"highPrecisionTau":["0.456662479182866245178717912725606214768735827831639784176881850161880489491509431035012418046147708569152573143239401039","0.604947097475667027921873189329075189713523485035068297414115009384520121791480783610288063280160762403951542723843423134"]},"1,11,3":{"tau":[0.2925152147483873,0.577198996547863],"highPrecisionTau":["0.292515214748387336425397847511154897292112867169557138948312245061313941398294647502653494886162362938236664062127471082","0.577198996547862997385514639543368527399602519525628821509886834583257606275471970337200355526009813754580425589840900764"]},"1,9,5":{"tau":[0.40346978959956736,0.599972479630595],"highPrecisionTau":["0.403469789599567367621499705171241816633633096767238056558641194517239648113897138810113779040251284853678956804461772058","0.599972479630595010198068397910829319966108733812483129471538441730385098778849356967522272839138983664466510897481457234"]},"1,10,4":{"tau":[0.349060588617138,0.5908817623809324],"highPrecisionTau":["0.34906058861713799817603544575024526358511567109872600816509763380977866832956843680073147374007370332955253402678005767","0.590881762380932347918043073175558490788580448548414932813724894241795884610091979108895607466646318611107711090495967258"]},"1,14,0":{"tau":[0.08262267706868062,0.4955950885169834],"highPrecisionTau":["0.0826226770686806280084176251348793643863707700744499446733065537874018386738516339853395624414257051861262814058541361264","0.495595088516983419423863040323973151940405954956935074708316727522878793026914724359397098042623491852475437508799298084"]},"2,0,13":{"tau":[0.9390379898501805,0.561229641883284],"highPrecisionTau":["0.939037989850180455420293098857457897459699167413390732631875382032748730680120838280433781632319665071284635151453122711","0.561229641883283907483878889206832888644662073409640024610276331756039938910412103161788819886486493559752137419510868361"]},"1,13,1":{"tau":[0.16529631384430996,0.5319367525905705],"highPrecisionTau":["0.165296313844309964471269148850368134614342256930670459226091767517302524608132019581793686363131384858391588617017505016","0.531936752590570446568094003784488427585764568778138287644652555008775179096291689281871515706987284652142937585530911372"]},"1,12,2":{"tau":[0.23230163405726711,0.5580721512709124],"highPrecisionTau":["0.232301634057267111740846126580888357736965624519714203717752130868893480062362811126698866353686970721226270139366037913","0.558072151270912326056165645288025020428962229439349905650068564692335248566817049374656282541842968002872857484326313557"]},"2,2,11":{"tau":[0.7823872261203892,0.6215576792077651],"highPrecisionTau":["0.782387226120389191688161812107294556566810302322186274233675056852081300956479955852512941031328532893487441552947842844","0.621557679207765061918999393531686546106858644652965285622619516575728304886284922000682966322561477790344945403521804884"]},"2,1,12":{"tau":[0.8537448968559918,0.596184335935649],"highPrecisionTau":["0.853744896855991751387445383451282658393640157772425374062050997387485230284306370280121349827790103649634937898015334747","0.596184335935649018323710013843589778354705980313445728361697908093864388951384442082154551408183005192890765979134648773"]},"2,3,10":{"tau":[0.7174308017927045,0.6403528267995242],"highPrecisionTau":["0.71743080179270450063960205635526800642213819741765256760204154896960167683074561377899202422539188320203647591980158323","0.640352826799524183823343402245998107834377698674316860283686116919536500651803359828749256610868020717419088941519623784"]},"2,4,9":{"tau":[0.6559788771969745,0.653883350026478],"highPrecisionTau":["0.655978877196974503674147029298175460048426337374758321600086498526494804747649104906133612530084031101640974689039403733","0.653883350026478006568887924546455092074878359126088652139445542578529952610932137052709836923525900185919359295132764094"]},"2,5,8":{"tau":[0.5965210663392805,0.6628582350258831],"highPrecisionTau":["0.596521066339280445785356763863142886175450704223173244233092472468539082564499464324486834988402626420390018958366175753","0.662858235025883129001794520048158448743726615323336552933980989150658969843762449396371188038219576667675145909282043749"]},"2,6,7":{"tau":[0.5380970098279463,0.6676683573900454],"highPrecisionTau":["0.538097009827946247898535708568143998818270598054402581747906773129380157639740271566145227829443631043299108202351601631","0.66766835739004535214156560804579272457481287975708501045586174143938911024770923511556157098961759263018134257156906462"]},"2,9,4":{"tau":[0.3619333765016562,0.6581261498587079],"highPrecisionTau":["0.361933376501656175165121565369612663076199878529627008789046289069872552486422239269500759116655877660226529292739140934","0.65812614985870782151569919707842576999924774769994195923984739577310517317321649932637799060455437160589753573191831813"]},"2,7,6":{"tau":[0.4799777699137084,0.6684970437085425],"highPrecisionTau":["0.479977769913708411564193305958136225727240415706044468644636784419542187432561561539939082765492091764658353945889599932","0.668497043708542449318725424970974396451386222222117652126181665461507964015034180587963512505233882330447866892148755844"]},"2,8,5":{"tau":[0.42149953800282086,0.6653631787860715],"highPrecisionTau":["0.421499538002820860150486479867410261463445719598032210524283029293692959092117124964482462653648759384545662777776935063","0.665363178786071582280082537978468430597331918939827326259274428675978137985592667635163980646456662530895376584579957408"]},"2,11,2":{"tau":[0.2351455387679288,0.6297556896278045],"highPrecisionTau":["0.235145538767928793556797262518982438848715942529734634251805253671826211683318144158371130358318454785663307054948033161","0.629755689627804461523173151837373311673619569355847920862856156190778591432375334627996805709423842725783444370768543597"]},"2,12,1":{"tau":[0.1635084235422822,0.6069776651315492],"highPrecisionTau":["0.163508423542282203815990145796572681479872741478588473785486974953522473538590143894892911239925751573442395211816272114","0.606977665131549213399114331584774730528276778869344205773284266449079638626761117120177448256373708886017308733750551181"]},"2,10,3":{"tau":[0.30031897885282605,0.6464563372141714],"highPrecisionTau":["0.300318978852826039385039639720183765241377462591103277865613158818783854064892945602064078978941688923532968180990804422","0.646456337214171418685646610968171159821697601578967299180562242965004391667297428453886977266237866985780120550668790526"]},"3,1,11":{"tau":[0.8505534170266602,0.6681321028597156],"highPrecisionTau":["0.850553417026660169180153010998914054623216055929647549533027505315301512578651996259601036062353714555194340863925181662","0.668132102859715602801175428267484225197635036886515739593932715376469865019244759989602798657989712340447695805134433753"]},"2,13,0":{"tau":[0.07775309367630712,0.5761916169396546],"highPrecisionTau":["0.0777530936763071212282694158099764869707992890804337156820979275800996815183252759251689904114482525897084472051062999391","0.576191616939654558424722863367963241692619553316862998550558248102685287771849981708529522267591334213447512801674729388"]},"3,2,10":{"tau":[0.7746401483939952,0.691388800792102],"highPrecisionTau":["0.774640148393995248910075690730414696374415046546728280592063590060025775724706278488930929006615861504651898840047084215","0.691388800792101940059876600833169644954949251316903029571060331983360143695684394926478302334079907827400613122624041665"]},"3,0,12":{"tau":[0.9392214786277894,0.6361186385564389],"highPrecisionTau":["0.939221478627789410053111811613889320749466595176675221867249463178886266815362497260017944766734741369864204219666058862","0.636118638556438923851261872095335893037120795174568076314600107500085353552426355735296702375687058318031920904134088055"]},"3,4,8":{"tau":[0.6380618539686047,0.7202953385181686],"highPrecisionTau":["0.638061853968604736028383912431870936001105076990227540434915667294405783202166053420774393618754569268860147004882546385","0.720295338518168645628815174642208682659538213585127218430312529094225974980529651269280381373779353027114400469799955142"]},"3,3,9":{"tau":[0.7047156575241047,0.7084196201617546],"highPrecisionTau":["0.704715657524104682481703803216035398569367446160785287407956379716229605905684058364564989575366103404929319421251858088","0.708419620161754644754979351095322224889214631968869983061271294056038184390646131746286192811274298542948312006763969232"]},"3,6,6":{"tau":[0.5091026627475932,0.7305789831195664],"highPrecisionTau":["0.509102662747593192608922474186961596326837825049901023929043382248773043059470370940574836459590297242428685122854076751","0.73057898311956641885977987403543329861210096774350578709162112680440210824288454074471570491823020951609316967607029064"]},"3,5,7":{"tau":[0.5731919241400807,0.727583552775755],"highPrecisionTau":["0.573191924140080709530053918607988417951448775738672467598405174367439054213872271479600796569167200594708429291095995347","0.72758355277575509093893449112350024909237805958305903497105694865715688590899935138173755822718268228240291468104902247"]},"3,7,5":{"tau":[0.44497134110717323,0.7293865572259005],"highPrecisionTau":["0.44497134110717326058430202840467028860832123980930122550772921367815672975209508480768843674083296484816991337805972631","0.729386557225900467201339690118124272505032682071813628618009501898131498338332783474041114658722872290604746514831809841"]},"3,8,4":{"tau":[0.3799734667640775,0.7239474966866714],"highPrecisionTau":["0.379973466764077504312055700263406108891759924411992251140060665735500695977410788407438421727246452082559141514759306742","0.723947496686671395732036974608136374377964144552039162271033659917343788629316923931331048669387017880946060530537143028"]},"3,11,1":{"tau":[0.16646062881828264,0.6785649440651692],"highPrecisionTau":["0.166460628818282641958724782883711819399914041772254606739288640659593515156731491108049242980260305430587584920200047637","0.678564944065169260869913669169753326794897807501039797385889878537966539811649001877419791494600228262810499744587780173"]},"3,10,2":{"tau":[0.2428484965374024,0.6991713857616352],"highPrecisionTau":["0.242848496537402405826701779327200920293397096968294277613823504765476937515629653500095225918959159487284108070098737443","0.699171385761635164952971110163947571186064282724410269285301180547380335701375859917493058774730874203656346571390427686"]},"3,9,3":{"tau":[0.31309967481103046,0.7140287572855422],"highPrecisionTau":["0.31309967481103047346879316221824343798604534218919408223790309335968573233553810325964868087448209292189292328992281315","0.714028757285542104443875664223620526035386931566564485417216168054651491050945616898852248250822960984978776450702840859"]},"3,12,0":{"tau":[0.07695254449550046,0.650720217086467],"highPrecisionTau":["0.0769525444955004535166969426888519635924348460918365678151891745320945346387586925435355554567706818153280334418009601373","0.650720217086467041126738773676515278879419072985272940130359940182522839245610743238322118290226767400338866845400289045"]},"4,0,11":{"tau":[0.9378252030207109,0.7085885540311195],"highPrecisionTau":["0.937825203020710954404753128828174743003057974494643916445905786816453354121629507768503059958372906717842990858717396573","0.708588554031119524715217996007842390551046093930760215144592061601267445677718414120499209951061684448897218543740938715"]},"4,1,10":{"tau":[0.8446410242074076,0.7391320988536378],"highPrecisionTau":["0.844641024207407645820448075858067142128763469781013413227711137864241274444752812152997921360771362120880174618626536124","0.739132098853637881938422437421300555448788078764341944536962792351833510777641528339727322628505460982341773299732808807"]},"4,2,9":{"tau":[0.7635027584323566,0.7610371847719761],"highPrecisionTau":["0.763502758432356626902719418607974545858596175945586695215395838880288357752469482658646365035721185348014776119509441032","0.76103718477197612628333740192663885956354372202553349827625424179021437255312744857213736059128844744421004468382638843"]},"4,7,4":{"tau":[0.40302847092294025,0.7901302292603876],"highPrecisionTau":["0.40302847092294025551757138958493170338262618108659022661451689935960147352705221959662344420997846745325870680503462894","0.790130229260387567558601195067768758801579019436773753562836850215789807327935699959959247605852081084797617556432648963"]},"4,5,6":{"tau":[0.5446179507896601,0.7926566339417336],"highPrecisionTau":["0.544617950789660084534995718637225445365375653659325721168037087417010503985927614935940662089119863649162326364012932673","0.792656633941733533330552737978315838773540761418410819282914981613658230138079182717700442406013512868529992708959929478"]},"4,4,7":{"tau":[0.6155649717259968,0.7870746312713428],"highPrecisionTau":["0.615564971725996754794708219765522603002210310651938321339022638157360997205988805704398823756942157738504960531201020361","0.787074631271342783784058491494394151153818934849924256043469174760206324333309263364278898859052168806352586437181790139"]},"4,6,5":{"tau":[0.47409733391301406,0.793665685409713],"highPrecisionTau":["0.47409733391301408644924071561683305224705373868063128669085883618771706811365657949771132159383034627656744507064292391","0.793665685409713010789353079962939972276209404729203202761304444293449617479647291099159402710073394532360366499175397215"]},"4,3,8":{"tau":[0.6880239303601602,0.7766968959091713],"highPrecisionTau":["0.688023930360160205938681151800946207331461593357099998623268316180013903334562874823566195411205877324710382691182231606","0.776696895909171335729132873514466586411802663005250362713900206585869729517379758615274046646895500021698100154478825497"]},"4,9,2":{"tau":[0.2544322156063508,0.7685878003335446],"highPrecisionTau":["0.254432215606350815085500846028375566007481958113949363978933606777436826373716990338323696886004951899237126880911421834","0.768587800333544623962413711846458641702094267770219983728981173650065062816639318844681085413586100928772949449579423637"]},"4,10,1":{"tau":[0.17267250605681517,0.7495097375651063],"highPrecisionTau":["0.172672506056815164898187271334359785561488400375450226548738422338336127938783072573496699842175003725229487348922480056","0.749509737565106334758206719343180171914621548794731711183550640768421302306615550908084714911354017008701003249669228394"]},"4,8,3":{"tau":[0.3303171051938833,0.7818953114392704],"highPrecisionTau":["0.330317105193883320747243512695454421990656035122855642179858565875656567549777849525067946498453110677947411985220495309","0.781895311439270404340019009070925791407966888069185471728751016977797666382165598874537691383680089216960005524500686429"]},"4,11,0":{"tau":[0.07837381518135739,0.7233159753550892],"highPrecisionTau":["0.0783738151813573918558831297750403551403834908040033146567803004700386925993908095492939729767041001048780275315336692805","0.72331597535508914508979140999170666287169580667340745018014947850494344515428156114082590502204941270219950919853521314"]},"5,2,8":{"tau":[0.748962091793748,0.8322478619817727],"highPrecisionTau":["0.748962091793747911462810903797625515516421048388841462518800907379400193805982367668895789757815705432457004674480715763","0.832247861981772750947900502018244497399123414772881605176041148386781506627513261276517528492560421366635088947615046544"]},"5,0,10":{"tau":[0.9352326502641265,0.7812742667233942],"highPrecisionTau":["0.9352326502641264418846066747548749950975390725048558444622327477057083068237042744165950078145521165476745898308777707","0.781274266723394117932299646236948161865647504448349507451327279837964923801982581112718617988040648889144842098048150362"]},"5,1,9":{"tau":[0.8362776797881908,0.8112374317078765],"highPrecisionTau":["0.83627767978819083762401959290918365041996294879615443229900596452682116952254682316134902387322057115498842392188730182","0.811237431707876554644069272673634900765437953017395986900747265227959390048945137575973188086796613749214786302743918407"]},"5,4,6":{"tau":[0.58784045709184,0.8555862365457835],"highPrecisionTau":["0.587840457091840009569869820227961706371162406285163613488266463862253224350780547374582224899139008055265730435176881867","0.855586236545783442084825023158895145365257943141296731863800802959258903560798654818056506118347461289687813128408126936"]},"5,5,5":{"tau":[0.5098005680822356,0.8592694766612341],"highPrecisionTau":["0.509800568082235608910185446429636005466474119876604140819343268783873336072770572311862006821625373780989951070638503099","0.859269476661234130925216198038786552693209240376950532320686677097094215001154025182478916119837310533152200758308991302"]},"5,6,4":{"tau":[0.4316761062611138,0.8579244195685505],"highPrecisionTau":["0.431676106261113823780432059240307320891110872245735041926810554456130831530112277260019337082113882761573296789234196956","0.85792441956855053069379451630150735543178444678847401030698153291623009562454868887701693738579469522121335657813609159"]},"5,3,7":{"tau":[0.6670321782724893,0.8467228186759682],"highPrecisionTau":["0.667032178272489334173774173558603682717357021131329097440037698865357131900831278837328333644787474600013857687579920253","0.846722818675968190222719071680175714887416849634203507594933116824162905123219082521724901822246026863112645935740025153"]},"5,8,2":{"tau":[0.26982659360672234,0.8396297548869927],"highPrecisionTau":["0.269826593606722322829993381769043196345329054114518008479607089119581054817097262893028255744911660651796029451891094752","0.839629754886992660503033280606920173733764244065905395864888581256895709259472897731382260504879650791842583961269898352"]},"5,7,3":{"tau":[0.3522230067210753,0.8514805214914426],"highPrecisionTau":["0.352223006721075317182668302663466339395033202889498403697781941598522085525088635488120942493510565275984586071357232036","0.851480521491442545449736437164982591131171915027100709767888404629354353972587959940563496105154814415578349966679812865"]},"5,10,0":{"tau":[0.08144441009020648,0.7964329392658757],"highPrecisionTau":["0.0814444100902064777811932113408335259461777810596545357510331494379418652638636130898876323202167855156357505857658851501","0.79643293926587566345391390504209706038641730814189388693666210509935685623817580139627750506296229100955353676977511514"]},"5,9,1":{"tau":[0.1817615933523272,0.8217265619928039],"highPrecisionTau":["0.181761593352327183650689760521144819065760783168181790224500891733944362603365516206913144098489927624471146274616850101","0.821726561992803964675453911516667335124614960385452875511213534262692965204170319438347455389299940615292677423110614415"]},"6,0,9":{"tau":[0.9315136688735636,0.8561148380510791],"highPrecisionTau":["0.931513668873563541632178823832794209939678680779481908754129770430387261269785974955930340361218513822565992448090189355","0.856114838051079180529280991132913424128379651462555343467771435274049728681734987340114858324529554305732391729043604135"]},"6,1,8":{"tau":[0.8253210630919122,0.8861313857660946],"highPrecisionTau":["0.825321063091912227490132787780417952877786395752456087702664014021652464865988479937534934324249124770897129111028470885","0.886131385766094653658708056401960401573676707700718458470578031662832727556763761293366799310313147685407945308386079525"]},"6,4,5":{"tau":[0.5537011311223006,0.9269966656400465],"highPrecisionTau":["0.553701131122300615617800625327372268921772150039556399102298585303081900384773116426975281737704385605982691096104806626","0.926996665640046502853983940245528506810177793952032843835981060095512058201207993593754565059651855894362011412810892275"]},"6,3,6":{"tau":[0.6409401629502428,0.9198399712588273],"highPrecisionTau":["0.640940162950242781520563591804943371588907026056162842046203633469230406645759507911736984289193185650086863681222924841","0.919839971258827344193342395830521666958649507619961760357086816517133044559865297506593016190894112546414647842100231979"]},"6,2,7":{"tau":[0.7305657043321823,0.9065251303388827],"highPrecisionTau":["0.730565704332182247850400544637427707409433465589812745497618293156000279394105968829313367680750154320726108294765533034","0.906525130338882686862249299009174055405219162664940033140708518022508613798721877465740197037363983745735404444977934832"]},"6,6,3":{"tau":[0.3796215978942538,0.9240286513383941],"highPrecisionTau":["0.379621597894253796519368929805963669958990129766110879118518588633431765095016185199160887792823711613444050048576719585","0.924028651338394093252819657286351310232018225162128273625888256665154462996794155156266601153051753615641606674121062997"]},"6,7,2":{"tau":[0.2894863816499319,0.9137214715466462],"highPrecisionTau":["0.289486381649931908146622858941633812197385235283853693114405023202661460610719785308814805046087625521578321937991403325","0.913721471546646125770511111817151111291173758034816734172739701444725652002197210609464957135964437858806863469280015066"]},"6,8,1":{"tau":[0.19385896371355252,0.8968321096910856],"highPrecisionTau":["0.193858963713552525620202484678976835443535689446148007496830282265379724948815484038807318697122945126602791391341165898","0.896832109691085612525853591719946532176323762124304875965658949855884810953954636763852942047007864523351674765486525996"]},"6,5,4":{"tau":[0.467100354906823,0.9283746005111562],"highPrecisionTau":["0.467100354906823023303494978694391800256248991903840687962658706562907687235348668031025877457152871839664584840903110234","0.928374600511156192811549532354911839353920081894531210871790190540450362219465575184210665728602923961313565940603707143"]},"7,1,7":{"tau":[0.8113066721543174,0.9654796296109466],"highPrecisionTau":["0.811306672154317400177237303731971645061351719616511778231227218621378060647181889930488440873168786233054373126614796194","0.965479629610946610703705217799959451568216929893835225433117273196086631746279449308578727021434402177709521900607580354"]},"6,9,0":{"tau":[0.08604657190331183,0.871950674074973],"highPrecisionTau":["0.0860465719033118310256557437971255312546273109164808160059023085275836843990554006353432259372531223295100764965053602977","0.871950674074972995976020625658988012223069144217327136401807297853026635938868460840225534289718680110514035051297123706"]},"7,0,8":{"tau":[0.9265597778316539,0.9348911668624053],"highPrecisionTau":["0.926559777831653865448925002576395985337908563966827810616492535015374541652479115830530009753796129806937136353026996489","0.934891166862405294436216262846716838640301169476429154206768749132102753285408390795536118874625234991917181863726231692"]},"7,2,6":{"tau":[0.707427625181293,0.9854008876195175],"highPrecisionTau":["0.707427625181292990881209040573967506036025409630300970068226818376871398828315166274121938477956907978612881665595115108","0.985400887619517456754285669738621391289979497883132720622100508102989621502518299494302659194609933038612422948037868443"]},"7,3,5":{"tau":[0.6083892978595518,0.9974076041389508],"highPrecisionTau":["0.6083892978595517329934543168546468159239167668183420474092632844423238510509269125668568100036565285195048321554161633","0.99740760413895084664543526345475565763785881686112205414184860942244482191757036637487955039573919360334033243248791847"]},"7,4,4":{"tau":[0.5112404316164261,1.0024317397033984],"highPrecisionTau":["0.511240431616426075362145922331625177718742736310083840245456839320658094284455695092832201740822440883375384075270207255","1.0024317397033982937088243224913384667646959392390669461354030661587274977364753512003101881546797403447983926801463513"]},"7,7,1":{"tau":[0.2094829500132592,0.9764458922528222],"highPrecisionTau":["0.209482950013259189317076332640833708377785450117751020637717948700762181098511885139685797158223588181476866276878396279","0.976445892252822170105623853998544368222199448690943220432536766594143285248027387113139519983196395823567738785223051146"]},"7,6,2":{"tau":[0.3143571535853885,0.9923110942909021],"highPrecisionTau":["0.314357153585388505649857507762616995537165910340721720376571704657788053362610324441188079364565563452199241693490389127","0.992311094290902195262847350977529906935084858555540095547943778402611206509266424656292688420050942108342308744826243105"]},"7,5,3":{"tau":[0.41392386608488524,1.0007736597637804],"highPrecisionTau":["0.413923866084885218636981115128579895067709203240043706761013831373029855619825583688119069791208290141326708666613854876","1.00077365976378030910450864622408215806359836425528147115278615264782906646953515375414583694616980723377229806035921942"]},"7,8,0":{"tau":[0.09231714157859527,0.9516408259881424],"highPrecisionTau":["0.0923171415785952688570325674911244067595308777738295791441840034360825486879200200073401933182846350124431683061926892995","0.95164082598814233257999575378238637720375913026777628366106699520316705823017112374463534496822408834613824992171845565"]},"8,2,5":{"tau":[0.678079871991247,1.0706316139932608],"highPrecisionTau":["0.678079871991246949576722481299505277351860204726906232977540912129077419554788474365935959009626401053124099226247066993","1.07063161399326087860667177240748713894487120148912065522037546198164513923122734944147307678141425201595912362889086136"]},"8,0,7":{"tau":[0.9201006942713453,1.0195431642974655],"highPrecisionTau":["0.920100694271345330924527988531544799327671441227108145508867660193303346257618735172528907437721141589866392558639723616","1.01954316429746555741821950569677121594183042559490147777269701606986180951773170783011345227527044018146772390305070651"]},"8,3,4":{"tau":[0.5671978963270895,1.0809393785972574],"highPrecisionTau":["0.567197896327089504438780780798436792615160068019373790406451746081442167147348473708956665071819752752350536907356304536","1.08093937859725731408573949678581023480658933138661252719097212363544598179998539374116049368114460108862809645121087844"]},"8,1,6":{"tau":[0.7933951559821849,1.0511789055564684],"highPrecisionTau":["0.793395155982184957487886615732929823495329805222171776730503058191748380001397867761773015137663901712080857684780677233","1.0511789055564683232502528538421287233407335641302941559523833510013333433024845802391239468184928656431434016804971272"]},"8,6,1":{"tau":[0.2295926172845578,1.0624129887181417],"highPrecisionTau":["0.229592617284557783850292813297115586610169181291956503544535956045853851540744021416308596194707770210261278004731872043","1.0624129887181416460958955748578984135842238718696938343658748424007712998440191569234736611833557766068781435515487787"]},"8,7,0":{"tau":[0.10062135055507786,1.0374636058220685],"highPrecisionTau":["0.100621350555077852079501791264944879538072991312817413078001564266465928607651630153019150568613646293522286076064167644","1.03746360582206846083767215244910540040513046719829912564424729853331691479906921084890755138708306102845343246043497166"]},"8,5,2":{"tau":[0.3460252235159376,1.07703315057545],"highPrecisionTau":["0.346025223515937609067032194086127004515355181235011717989196841161124475820175479045500833098880276398658143370823834683","1.0770331505754500889472892584472799430644009939719649414260449610047032107521810287651630411444941991548082455637892936"]},"9,0,6":{"tau":[0.911654535950164,1.1124703911943978],"highPrecisionTau":["0.911654535950164032358854463930990764331626187767303069448315749963924368636853440182714233759423932544418835924576550792","1.11247039119439788028843247391923186587129052289330000109066056273204915884248133553449627241074005520669318662352969951"]},"8,4,3":{"tau":[0.45740811475805754,1.0830298959387457],"highPrecisionTau":["0.457408114758057560557879825208974353596633029612772429396009817258096731487178702650077489387252025790022015126045062203","1.08302989593874557835209147128360778082924831307970573897872625478478981121597766425460476120598888475622427821318860799"]},"9,1,5":{"tau":[0.7701875507980696,1.1456255711794847],"highPrecisionTau":["0.770187550798069601334753141434923501683810574189313175344020218303686116296493076159665773722111193362316733762923151219","1.14562557117948472923663166926964171136332193409738427847857798258394803888720840340279326004705222253301136695340684303"]},"9,2,4":{"tau":[0.6401217215893285,1.1644064399048693],"highPrecisionTau":["0.640121721589328429198563748791783055347110427185395076543934978266970007471508807258227699666146921000572804744166096531","1.16440643990486928096430758435782262959211247470059318851364045652486605655170653267462691262316401186982363565512022795"]},"9,4,2":{"tau":[0.3870880611237631,1.1698619606860758],"highPrecisionTau":["0.387088061123763111140503580807929658995609090153492488119966174292848040857849600716931397432474947248973840715183877339","1.16986196068607576252079039219463359290909032950422739050987007363519547370604008213063749498825059487218137188839480945"]},"9,6,0":{"tau":[0.11162289930970945,1.1318630179892937],"highPrecisionTau":["0.111622899309709447957283604437787600173006402220702867379650214532489381398877802833273592231512356120726248539754362234","1.13186301798929376226209984112854040656978190140538362057744581372517672028686450420848209354356948991573646283614086241"]},"9,5,1":{"tau":[0.25579848702560465,1.1570438638857115],"highPrecisionTau":["0.255798487025604624412655310551545994316148755430308004954707333903980146512553350116987992074927262589418212631590858969","1.15704386388571155073850565041834454753929634153421277038385050932031880899490624990658383744993459361552093325711806495"]},"9,3,3":{"tau":[0.5137913009078899,1.172206668726127],"highPrecisionTau":["0.513791300907889923351589733214045396440203972138290551439922541973610450515141310332315278668889125812178849770556042765","1.17220666872612698119612434254381721311603064156248042627997335761062449235059917719855631819732929330399362973422625508"]},"10,0,5":{"tau":[0.9003947433945816,1.2169487947302262],"highPrecisionTau":["0.900394743394581603816342677801698059549726087992732568613494362233398554823503093507053393246315604249449999838239286343","1.2169487947302262473426463764178156315657758705472247458933227878773655704720967849900735615260876786046999900973029364"]},"10,2,3":{"tau":[0.5894663508721522,1.269616704664962],"highPrecisionTau":["0.589466350872152195502912169386587212286270785275142281429678648553957354657123890916069366370922203413069927616784206717","1.26961670466496209366365586463565740929360007366781912437850870678004362283488029514713683193602750100172158824997373593"]},"10,3,2":{"tau":[0.4419462894619868,1.2732601630763927],"highPrecisionTau":["0.441946289461986796758959830042530178091639148376838507521554191324736100646363256529496528475476798101861544356310296603","1.27326016307639259183334143148436637335146653078271473335578979957942973488895486630485687291674617299903478326627513116"]},"10,5,0":{"tau":[0.12647490830891742,1.2381831583236318],"highPrecisionTau":["0.126474908308917408868490760447177491379719637697386454277587155798013590453035146808121567398689899356007708383491321929","1.23818315832363174220498930944287980060331782143107760328801689413023773783752727381874593283674089643391658301229137666"]},"10,4,1":{"tau":[0.2908343646446471,1.263446344585875],"highPrecisionTau":["0.29083436464464714674186503336535503671095559786607247223552271489831575001596181148396330835655479151838391692609022287","1.26344634458587512447951365646094514128657861067823623348192289240984309030186096016060917421072447944691198049677562374"]},"10,1,4":{"tau":[0.7393147831542055,1.2521047269045147],"highPrecisionTau":["0.73931478315420550128982673158739722483249777559047073252892590122096970522429525184015614371678060339292221097973630993","1.25210472690451474504380423972950428369185785715709957644426327142599130241634426347521478618625055646802023022431133521"]},"11,0,4":{"tau":[0.8848469109483587,1.3378593122853588],"highPrecisionTau":["0.88484691094835865299417025817583285545221911192022155244390963760531058333950344103354093444055141642672038335916887574","1.33785931228535875398443367946213114885048882148014436653359419983890611903662375582403602320067221642368316459127058846"]},"11,2,2":{"tau":[0.5185977248059518,1.3902063954719277],"highPrecisionTau":["0.518597724805951775925830240443998651316195307861319336624347291123844165847190015775288149937369567364617002985020071576","1.39020639547192774729894519779188605103643380201757934171012380561578714504631712903615081796432421472152083820302753177"]},"11,3,1":{"tau":[0.33961342619206547,1.3860439393295256],"highPrecisionTau":["0.339613426192065490018951061617363579039641297872758365501947743269092740623591350419837027343363055223364851405368383983","1.38604393932952553847273846927400896318233115477866813804032619475501839590779221929722114393394218547021192492706759434"]},"11,1,3":{"tau":[0.696513405543701,1.3754573542960915],"highPrecisionTau":["0.696513405543700978118144965223343281520914603924429116108756869090534700660089584361928521344363027661286768870519703996","1.37545735429609153543428028051767544857954977467107070060293138588256340188783834784711918889575227636847862873733680424"]},"11,4,0":{"tau":[0.14725951290108097,1.361391850078212],"highPrecisionTau":["0.147259512901080976058674048740872707022607795169539731698551006606955207016117487209840350799872669704172694709427979079","1.36139185007821199154485498659268540971328131703260825854169584955395345405307400723978208788491406125902036197291691508"]},"12,0,3":{"tau":[0.8621457225759497,1.4831804395676644],"highPrecisionTau":["0.862145722575949635221128011598642176320406146806188924102241372661840877823319424212293038775949269877277919367681897206","1.48318043956766444028453797759843312608552869235265579468305606383031324715922797431083382665882931079713494033911793431"]},"12,2,1":{"tau":[0.41184964123519446,1.5313775566281667],"highPrecisionTau":["0.411849641235194469055781505724389469504191275939636221664586638651002043358391046468807144664225245931313262501486988002","1.53137755662816671369654570113151659134348448579970214678718680526099781673470522835714856620653123722184802134941267727"]},"12,1,2":{"tau":[0.6332926962736524,1.5233484371092858],"highPrecisionTau":["0.633292696273652419758508578517319485859827463446353063358609367925589376559891188551618291546896727556714084941058194538","1.52334843710928582882206642098665991093715563762803922989353672761630111468798551974786783754966391105964984586032802982"]},"13,0,2":{"tau":[0.8258938119441659,1.6675688592881914],"highPrecisionTau":["0.82589381194416594272282047866263254271102046295518055235483079836295346128353119985937546734021400033758747381590480612","1.66756885928819144765943284669296346657459915420405497315180062753824494234154477319519059812591640982319163832315130665"]},"12,3,0":{"tau":[0.17807551787122455,1.5095288444672899],"highPrecisionTau":["0.178075517871224542147046119305059752362430337310118073362131131639367348214949481673459921682556845053877252275575814068","1.50952884446728975774428623619932796345666992992111921293385389149677118958095622321236091050644581380580257430855661546"]},"13,1,1":{"tau":[0.5297866513065048,1.7087869806127551],"highPrecisionTau":["0.529786651306504782639462382519789606319414815613125554162041724217924405019251253053950002391721248908046127272078143221","1.70878698061275511660742885656939106253549790825290683970362647972022432830580142654159982556806213747671302431948073086"]},"14,0,1":{"tau":[0.7580881939033295,1.9229755802635837],"highPrecisionTau":["0.75808819390332946630199254434724823236659013278158510713985481418956134515509037634931455669260752579657733914461186958","1.92297558026358379013071526483731759411551943782922021335253960184828359695818455881932337063738171891897375220216166328"]},"14,1,0":{"tau":[0.3239116760464443,1.9532455868368472],"highPrecisionTau":["0.323911676046444247290443380675046830230106225443353955679741318364258560174773437372144989428319398168640742290453579582","1.95324558683684721200986373556817492786974309144689437803067458444429797225733507069456190304561249571140486792854608598"]},"13,2,0":{"tau":[0.2281945559519048,1.6970082560458013],"highPrecisionTau":["0.228194555951904803270059983291196649991286453441600060086551344621474521491251057182149209304442197405693459467352832063","1.69700825604580133323051370177106728894562827668316226447344822353876003249112885593721094798291770366536717948209493983"]},"0,1,15":{"tau":[0.8496019165586416,0.42660838007694074],"highPrecisionTau":["0.849601916558641551720773496071696130100595602752900157360905225995692187209569543907900967510986826716681870384339143653","0.426608380076940716356269152779621491137008267485851736166213283398147252082864835865252701487437042317162538155805526118"]},"15,0,0":{"tau":[0.5781378463077104,2.342590652742551],"highPrecisionTau":["0.578137846307710356715604860695673265121758562655160846126352097398029195879442251710619296716299775937990557039406398761","2.34259065274255091830428356536346479573591670273837116261334099878819163176668367138854381704632503503729817068427269542"]},"0,0,16":{"tau":[0.9267873915942256,0.3719783286753756],"highPrecisionTau":["0.926787391594225610699046801029893508886300304667639669726359832131224875356228560696827107650254664607197370297155310066","0.371978328675375628167470030019272373439751457004827300791862610023890994728619116538404652681638415322929977041532781361"]},"0,2,14":{"tau":[0.7905962650496281,0.4618359238228509],"highPrecisionTau":["0.790596265049628103333862983816532694098097059377455089996039309914195091328091828970889202936608380984858567703369556978","0.461835923822850856575790217185380694198321337550502099375220598124784333765408015506751817194843944799939032501664811322"]},"0,3,13":{"tau":[0.7385047487908512,0.4871815984001738],"highPrecisionTau":["0.738504748790851170851441495570246397003955063694338390332160489734741411107084991148521014783678570771939095664385163051","0.487181598400173751173197970989150208277388106532220387117154243131988130553007824897634634864736501770104620961683511169"]},"0,4,12":{"tau":[0.6899863538203479,0.5059545842241427],"highPrecisionTau":["0.689986353820347947338677194731078996366419797043222654697427175832813267276565769882205961583300958861574415616081519734","0.505954584224142727152929145142011816578097620700370277263819143999989797012924481964309784935402771724499649264639765834"]},"0,5,11":{"tau":[0.6435555560584743,0.5197168229277463],"highPrecisionTau":["0.643555556058474199245895329404324479817898953435307596024251333384382874965002046756061020751956847000893257806805355244","0.519716822927746313637010698291443141957310863294487409638426813568633467773337276994002170915006630048020259024166090804"]},"0,7,9":{"tau":[0.5539619790144574,0.535217242023182],"highPrecisionTau":["0.553961979014457401484570991000629381262010344668454429963914148422924124937503928438178380367696147589111416515957077448","0.535217242023182016167736803665526164187142217079979492700910728855631622955319961104793415554265932992292298402467336553"]},"0,6,10":{"tau":[0.5983911740938758,0.5293122600049306],"highPrecisionTau":["0.598391174093875804803441021338045208086899882826692385004906811140707614018239679390427391371279617601938632179544217725","0.52931226000493062142414787461119549216379574960875490884618274796745888638391182584227241575264399486892749108630403555"]},"0,8,8":{"tau":[0.5098727163913618,0.5376843408427727],"highPrecisionTau":["0.509872716391361756892119243559869666128722948477304009785881156707289915232523742476099239566370227883363031936891399386","0.537684340842772730671775935352273998680830713516624276685870113608041655560077992564751633713182536181011932080490665139"]},"0,9,7":{"tau":[0.4657847041191043,0.5368041799331726],"highPrecisionTau":["0.465784704119104298197619661497152376577579574527978466277654518327159039806947523073876682661947560610265997668806001568","0.53680417993317259451440933502335477719530351589653111245277373024963177153612088385731712691268391085161506414195773668"]},"0,12,4":{"tau":[0.32984094572509065,0.5127427551052125],"highPrecisionTau":["0.329840945725090647110194209948107603473187920514964056356708998390611883549021554568329683106063680462911984188762419612","0.512742755105212492027211085211967459166959198083674197909031628045968630084630783617810163926791231661764989496338468282"]},"0,10,6":{"tau":[0.421361868146264,0.5325250842312076],"highPrecisionTau":["0.421361868146263953413294741326496786585617534354672447021492459024532670820205711872702721985046497439599587460851177785","0.532525084231207589593105682577261112976019199483472983757594137242835327649982285792214268747172623220325835477173807296"]},"0,11,5":{"tau":[0.37621787893028574,0.5246412444958528],"highPrecisionTau":["0.376217878930285747732279189843683055186706543183219953956045787316143192421578709232661305558096218582389998882259378246","0.524641244495852821341268671336019802743365802756923777400880485346221978509893494413607669105522775321919336096523996518"]},"0,13,3":{"tau":[0.28145114568162616,0.4960961268042735],"highPrecisionTau":["0.281451145681626158188704153988150655746097788241048980909511074483084883860489656104666237605249555307979827172068609443","0.496096126804273501978909835621011581308559140975557656668103674980933290057894963354598999402448194685303998465575120959"]},"0,15,1":{"tau":[0.17139778493093405,0.44175867458981644],"highPrecisionTau":["0.171397784930934056730331279410023849103303214026966003348364894336606097886288196559590581620688217722845912676993436244","0.441758674589816463479790837152336881750032204798602416153390033935199996785515151135608764146981055519508826264745888113"]},"0,16,0":{"tau":[0.09657326880699513,0.3939010006318137],"highPrecisionTau":["0.0965732688069951334479397971966698899005253010874964040777576661934700340359634335556155416703344967980327928916237926165","0.393901000631813716936247564773116641387518261910352262334607505419162688132744421737027145831641892627173995406507206771"]},"1,2,13":{"tau":[0.7947354127728681,0.5393722453045587],"highPrecisionTau":["0.794735412772868030650186283276858242231193507961038094568073019211262048375335233160628930586924551750504753099643873836","0.539372245304558738604939324108300361503453201856686675271716977316581770108133485597874381079980677439147029202536044284"]},"0,14,2":{"tau":[0.2296591487003129,0.4733570694958038],"highPrecisionTau":["0.229659148700312910897207577230656222006579160551397842886837516677036618460134043552128170275008214191545596826686334036","0.473357069495803821123341105393497649733295230750698406671716603342381661958034662926806887011613366019747232895815321842"]},"1,3,12":{"tau":[0.7375475309211614,0.5607135605418765],"highPrecisionTau":["0.737547530921161395076011750962724803426444632195973958742312407674807862485320078802571181694407863386118257138923697921","0.560713560541876434405585020617710179193425861422965063537050988953828229534582930927037783965059393127379315662252991597"]},"1,0,15":{"tau":[0.9382708620710256,0.47083943619788005],"highPrecisionTau":["0.938270862071025643633525282624273601826794735737435410856823527849281300247056616249839778513514878835552222989801726978","0.470839436197880062117902992678930951113598686739084867190295607849788975707408319779082954570572292963070071522690207003"]},"1,5,10":{"tau":[0.6327290579265455,0.5882088137724734],"highPrecisionTau":["0.632729057926545502278100087121086413240598571894490291744432261060917637314289738600602423232959179864707570759578373835","0.58820881377247337080784129135320214829059185798207986274897290492585350958858416364462742130980372445537167328064450527"]},"1,1,14":{"tau":[0.8587151156835329,0.5108276731850847],"highPrecisionTau":["0.858715115683532856470327274848958145300296350873126400629666825768236187296352737186334591465443364934008462846160940846","0.510827673185084757583913223953672404953002470826493846705877491776188785304346279199687390579318693131806076388514193604"]},"1,7,8":{"tau":[0.5333723486439484,0.6000037801298679],"highPrecisionTau":["0.533372348643948469932178604229135401089274380246327597663952695887810973402619348794242091144179134656935007238147791102","0.600003780129867847084725664215713334992273658048281489836738404694614695558183612667595243968096689495221572608877404653"]},"1,8,7":{"tau":[0.48424962918587566,0.6007215763468017],"highPrecisionTau":["0.484249629185875668809949086011729518794361407204901770020711522145168040151840036887918052809661037040362232546341600092","0.600721576346801715967666010961774199252180937806712140575349081235770921113699146245614961576654338174933152893393171781"]},"1,11,4":{"tau":[0.3334664873390728,0.5819217406882892],"highPrecisionTau":["0.333466487339072777381420357397100490116016726762090219687971191186650534682769253817892867220004488499063272990460072573","0.581921740688289223711686230434358255653971802380881462379503962123720965808481147880009328253858055592772293729411002128"]},"1,12,3":{"tau":[0.27992093200227075,0.5677063829920862],"highPrecisionTau":["0.27992093200227072827287718301882588570579037844563079719889694310768033385750248565214878705565447349943580123630083725","0.567706382992086243947311694180707708526802656152024792389976454582267505126061330510791226497043405619813189915140656918"]},"1,9,6":{"tau":[0.4348957839538801,0.5980503342279124],"highPrecisionTau":["0.434895783953880111007278670345883660926037454247372328042114766603918708831530340015278339206304634013517417148744760314","0.598050334227912433373131927319992206235464162793234945946315937000565794267484183744663966306075458565858289772297854011"]},"1,6,9":{"tau":[0.5827057461227122,0.5958822955361184],"highPrecisionTau":["0.582705746122712147742173559633080886980335644086156689621432660823274222795780930750323177404426214875136629787523439121","0.595882295536118358909629422980349777948467377354160321631579474299805078475389553694812206084458568533178421711418142788"]},"1,10,5":{"tau":[0.3848343306421952,0.591873926080061],"highPrecisionTau":["0.384834330642195180295882841433890119477248018129787055236046579532963372404049910685508499872281082675982757222425851766","0.591873926080060899391019839265256689481563003007895514216228261868042176044113576932815377985847323138978691886066532359"]},"1,14,1":{"tau":[0.15883184204088072,0.5224349092627871],"highPrecisionTau":["0.158831842040880713247254542855688346595378414505306698095189285834434773454495391082166421061802100927967510987755441749","0.522434909262787145460810139641897120618386016216443825750658482405674057068563232343924245184695012038390200468717821111"]},"1,13,2":{"tau":[0.22272076147714764,0.5483851782859462],"highPrecisionTau":["0.222720761477147656612522678676036462092865430868129015084836172552322736120250667092297109532538359204312017757270213978","0.548385178285946243328820782493914247869443084365052636749092962196285152177034904836380356648165271240730261446444437786"]},"2,0,14":{"tau":[0.941563998021656,0.5496890895577903],"highPrecisionTau":["0.941563998021656004988597369525923422833331440056421434621600534249955693650454379253438436005789634040305399560895285422","0.549689089557790323005777569701993282861976543799668486332639457978986056831900932430556515205479259655193788843828098163"]},"1,15,0":{"tau":[0.07961325716528028,0.48677452382241787],"highPrecisionTau":["0.0796132571652802846007996221929994579745451634850723394399255136120510218720834616011603482883657200603590101069300380834","0.486774523822417860778561504198288922525142758950082734902086501144300987006590061895418325942382574907674383833966945342"]},"2,1,13":{"tau":[0.8601909635222494,0.5837133433964489],"highPrecisionTau":["0.86019096352224949286616968895881712368000436117940814651594145643146010449789104111370324836015694406292266112538777889","0.583713343396448867500353484057963610333971003288935145899091839592687833644441322161354704696521603102874155879802148645"]},"2,4,10":{"tau":[0.6731198456825497,0.6416372705075959],"highPrecisionTau":["0.673119845682549704688701276906469573950895331354690985668302547792944354535968537560559082233797466950783886211150622694","0.641637270507595919658434035358631018358545742392912969842839005384674726738344680010136120182968487951373120639402669606"]},"2,2,12":{"tau":[0.7924453897710773,0.6087442472076777],"highPrecisionTau":["0.792445389771077288488961572884217871802580856866691256023993816197250290313638975317820716764050275960438753094549063369","0.608744247207677693035226177350131813546342162343784824849876846470111047521722576718882640709881094944089210242897423543"]},"2,3,11":{"tau":[0.7310217102471352,0.6276310643772127],"highPrecisionTau":["0.731021710247135280517450436827095077906134734336833737420289325221409040958596742258857261347563337591224016744963121014","0.627631064377212712746800121683976893396347117377707050830541675911188699324344543090691080888983928421178536545570313656"]},"2,5,9":{"tau":[0.6172956493259958,0.6514598252779894],"highPrecisionTau":["0.617295649325995768539801226843135290859730189245327589527292550600352366313361733622464215174288872367106699758476972332","0.651459825277989444481727457751256164793648422114930187646851168581125386037342956789802550777153152938595650724091837377"]},"2,7,7":{"tau":[0.5085414040748613,0.6599780423741551],"highPrecisionTau":["0.508541404074861265517191785122127158968737654255989493809991553641398559725693792561157581791687420334599085316490244799","0.659978042374155121807234352804330699677480786577914830618400779001017772448885733305893003731343058848126321999937442207"]},"2,8,6":{"tau":[0.4544088942040294,0.6589671651544582],"highPrecisionTau":["0.454408894204029421344169510387423377186962740281823822632078981888250952966847447946294837941100170146659352707792143727","0.658967165154458150607702837851368335843524077293246330493021336310665619116293909729522420933235334343304078168822628651"]},"2,11,3":{"tau":[0.28573034605867714,0.6338437033511397],"highPrecisionTau":["0.285730346058677148417667015574024120015266321870880712569949832135411654210081790239640643804527252254891144339684627989","0.633843703351139726274182671077977184829938966797888681862478317649635985552692059858832728793436449148245086847192752968"]},"2,9,5":{"tau":[0.3997035079662859,0.6544217501283516],"highPrecisionTau":["0.399703507966285879140529330209995654911606712340238571917482970739552615258158809433622276936498503796444157796261261492","0.654421750128351568730779735506710040031030324807179992355703641536738360443817157574670546139127516366243240664403734024"]},"2,10,4":{"tau":[0.34377555379949676,0.6461626951712005],"highPrecisionTau":["0.343775553799496771065704004982109133977180223721946814953725297006924467394468036655704791859818768192629509444460305267","0.646162695171200480621651796793474326786975670611919093653701294349261255464263294224652076050569543955071885318650377589"]},"2,6,8":{"tau":[0.5626529551361849,0.6575012923030078],"highPrecisionTau":["0.562652955136184863524665576365416749859405924762173475758408514719976842324686268042709412846842668269818590810475097407","0.657501292303007800765876978238712986315073664920815827772734346594386767363316038972104352536072433349947519301154571477"]},"2,13,1":{"tau":[0.15616163749851897,0.594238769074766],"highPrecisionTau":["0.156161637498518969283328925850071471587241694148418101225600213535094457149288327258884864455864703580990678185692401054","0.594238769074766029189800331914042310998606014016920479293536290053873994791477097407784569362922386373279544823411734843"]},"3,0,13":{"tau":[0.9421077415035402,0.6211604788039776],"highPrecisionTau":["0.942107741503540159540968825461777833960960504590054366748882286817868007241442774386360481589146414376553553331376691173","0.621160478803977604649699338673853688117501060216455951832574040714714068830970113264645365106152792681515910276221727745"]},"2,14,0":{"tau":[0.07443379797671193,0.5641181723702157],"highPrecisionTau":["0.0744337979767119269384586349892570451084445466508561232084204485825431967694085666905910345732284162437689732020609723551","0.564118172370215619635440926345702767708921495913572066671711101308262501600611111889751277364854379720051842767895108142"]},"3,1,12":{"tau":[0.8580155041708522,0.652188087087183],"highPrecisionTau":["0.858015504170852118259510331697946774641494724358450369067547479637362332779542741471606913761534608795783998599155353238","0.652188087087182971284637459692629533179155787628966799600676747504183375314515961493577215711503580768256640812844493112"]},"2,12,2":{"tau":[0.22412618991994268,0.6168722534907474],"highPrecisionTau":["0.224126189919942689025309969041240772818090748459825670222376698733434785536066541334419925823219089648660862474876604355","0.616872253490747379777099899814564209833113226316280360662908811876041822486392339801457772231266306712408898044324140921"]},"3,3,10":{"tau":[0.7206201368188154,0.692300708951353],"highPrecisionTau":["0.720620136818815454928092317278611482513999143415987585068911636783973559421786678426668660178211344833848368056006829914","0.69230070895135303164848204224033039071705254064287765434717647691303795986732370740517267026542941530900275927160967834"]},"3,2,11":{"tau":[0.7863584085876046,0.6751114935755776],"highPrecisionTau":["0.786358408587604672482956139514981370193284542292999744559132967959502776487501094664287164999859798039585069910683695642","0.675111493575577652937640656180276812049829346908030165563989761551917183432248955268151546972023305495985128921969889835"]},"3,7,6":{"tau":[0.47894321931010886,0.7183578623524647],"highPrecisionTau":["0.478943219310108864622797983396253415224282247325808803470349911471826923665565471909632051021934598543694542191601343005","0.71835786235246471241673911645143642894693653710355131144553544027135606438242745074196081918717856634230042463979905482"]},"3,5,8":{"tau":[0.5976830269033863,0.713096220486911],"highPrecisionTau":["0.597683026903386219721804663819647169330233446752155382161409009915065836951191666927413463663179472651949342735412356954","0.713096220486910982813720695518007166086684137706775874708270434946448069721315528789867914483764357331839325909581298236"]},"3,4,9":{"tau":[0.658195791704761,0.7047760864267403],"highPrecisionTau":["0.658195791704760959112928945670434516410246326623489422359486388262247340386071275525192762522481519858077051801396247665","0.70477608642674024441604104583105448062310931748042000084560544922269953240231570836184414985089189668581317605473404041"]},"3,6,7":{"tau":[0.5381675190674902,0.7175739010132858],"highPrecisionTau":["0.53816751906749018030481139365683253350431094424127161769478594871971858097729330350438873205205128573699135989823381737","0.717573901013285724846229008767371178877629986472677297745734374313386893930187636703570242878548106188075930686763010207"]},"3,10,3":{"tau":[0.29608819271815456,0.6980592196743013],"highPrecisionTau":["0.296088192718154571903522579751387839465413423343279110087895314715724462849438194540558252920570688958036815983522067217","0.698059219674301273585434574629540131859916429811017360626939277459969044059814302841514397639491569585252451129003816781"]},"3,8,5":{"tau":[0.419362541791453,0.7154649481359829],"highPrecisionTau":["0.419362541791452981754326223422339706729101051826632001929900029211546731446915775446224733010006392078332868760321861978","0.715464948135982975939449250212754288101558638734792241290849337352385001414013854919338148030657591224935821236368717818"]},"3,11,2":{"tau":[0.23006113212835375,0.6828310564029121],"highPrecisionTau":["0.230061132128353742579608680799552976727694558443894417936863650030787907873838849017014251484791263639810758748926753952","0.682831056402912104716115747317622217924228811839167441217153246563025727682310355096125307903714619145690175773361654752"]},"3,9,4":{"tau":[0.35871734682894074,0.7087846147542501],"highPrecisionTau":["0.35871734682894074662112766531494319425796267379140763439955506255178365576174251154923708310595035320816722977746657835","0.708784614754250168059003213220314925544791724751240096028836594964350920859836206253053757040393158920692432256307640638"]},"3,12,1":{"tau":[0.15799768451882154,0.6623252578937677],"highPrecisionTau":["0.157997684518821532165241174984649929335872180478612680090375918988533420987740383657091667864848627504087584116089496259","0.662325257893767729457421520584724120600718153685456976476978601246312036378557560317454258053795414327899987335823339596"]},"4,0,12":{"tau":[0.9411678598891168,0.6898248066447193],"highPrecisionTau":["0.94116785988911681077837067929658006130367944479828521360280513547216630474066418519709480323810124196904072440992276293","0.689824806644719385415494138154110463911398239025055641587463244410547575326340568160306615389423765144149056291714725117"]},"4,1,11":{"tau":[0.8533596710212288,0.7192934194830155],"highPrecisionTau":["0.853359671021228823268763805257510677759282490223883497379217696799812407570257005334519892727416611021645773981167823352","0.719293419483015509258341987411432441640784054691504328436907586414779709589103576558499712344089329480783395746202139343"]},"4,2,10":{"tau":[0.7772621130398515,0.7408743435149233],"highPrecisionTau":["0.777262113039851514799260377637450935176517971562423096572969102560738331461954418630255583287782702409281136940037929861","0.740874343514923303479645843923901226064568382682568161093486120814824183098382966355914541263571912901536562223658814011"]},"4,3,9":{"tau":[0.7067698008679769,0.7567831150126075],"highPrecisionTau":["0.706769800867976920183586072651724330639539552582745533860961657475368073249584458678390074429807068566378664336961312002","0.75678311501260748837445901649399151817222490215572016607646327571970953978776548394947549276591604779449294336644930263"]},"4,4,8":{"tau":[0.6393794095396858,0.7679306165623513],"highPrecisionTau":["0.639379409539685752490595893742611703743752310206147495896319853443803816373527158310451026680829220073648806476206694975","0.767930616562351260505500599418901519144690662785305215160797121729300640600378365116582801018740169134770949406947302891"]},"3,13,0":{"tau":[0.07319143838669522,0.6351601781483206],"highPrecisionTau":["0.0731914383866952184809512939276954785031138706864196866643123440300452695517083541916527342776359012523847466521120763595","0.635160178148320579137298768324621437466669553066694993262114508653946732462453527789526372595624825126657319337911447697"]},"4,6,6":{"tau":[0.5087560422698675,0.777638415648086],"highPrecisionTau":["0.508756042269867441824713092051571940929998685768582105696598316388781010874773699423671963481022720162372929071662623313","0.777638415648086035687933002391380265427140912630310767533356652705749074943906715782671239011393025230517629797882332766"]},"4,5,7":{"tau":[0.5736941341192566,0.774798578970541],"highPrecisionTau":["0.57369413411925660529470755029835968828962937417931786743340483311585929818382176315650707636019123131665499153713873611","0.774798578970541059821138562022602118112306727656305089641162759057464538900000368788481460362306462347595182381612182375"]},"4,7,5":{"tau":[0.4437717670776126,0.776541019250416],"highPrecisionTau":["0.443771767077612647058085685285062471837710588766351560189637731456162922554334421362298273886185206019365200083916226838","0.776541019250415914044750561624832836215927416229066541362470351409733675834276047154595946802444792244488016604907593189"]},"4,8,4":{"tau":[0.3779454574088418,0.771458947862016],"highPrecisionTau":["0.377945457408841753507557215812790888393695224586711926175168430934008937223869519920599904938513532359550322766315677807","0.771458947862016018806744449480951684963375816561061253883802152693803015851115445747599441423596810268988746871239817998"]},"4,9,3":{"tau":[0.3103105519466542,0.7621990685713833],"highPrecisionTau":["0.310310551946654199151152353502341340037094821523569069666534909337368316393135466998779309551140077326025623879377552564","0.762199068571383268778449699042003499498880165570091920796998826156727402920287866496888327671592829939409830112489178942"]},"4,11,1":{"tau":[0.16280817365660716,0.7293472585627945],"highPrecisionTau":["0.162808173656607146029079902671712781383802184002969333151393443887762091424995907544801370747823542665759763828821371817","0.729347258562794475787950843996543486061879325601734592121308313321557572672032752971686401986249076684337816586450144305"]},"4,10,2":{"tau":[0.23945018273070542,0.7483828208848996],"highPrecisionTau":["0.239450182730705437997510926418513251408727767108404789596275104101975118945014702969288994761000217449567151237220387144","0.748382820884899607183744805486937229737947573646440900314773436027860598182166602082472982941193937757523268325071657118"]},"5,1,10":{"tau":[0.846576330039017,0.7869153136649272],"highPrecisionTau":["0.846576330039017044681758366613240606186288767723837615734150044777445239657285279167670319226235758130331752883997662331","0.786915313664927176781644005469765772697701402447904670452470631737585995245156328039872785844547069029689228727377194592"]},"5,2,9":{"tau":[0.7652843120864515,0.8076191797798643],"highPrecisionTau":["0.765284312086451540379520858769704850867083640381274888270871285723449992585109022983130224685888673822835329488346978928","0.807619179779864270195036996661093237384549595820495483626661789195023979815191479457518709789685093695392935494250228753"]},"4,12,0":{"tau":[0.07403434843744677,0.7038641204015469],"highPrecisionTau":["0.0740343484374467666282142122339319492641331164786705901466751978685407124806518713931573763119372759113584252493199963102","0.703864120401546866122588644762859293744855692937216134118935688121911438620892595202904042271593685740466132409516621346"]},"5,4,7":{"tau":[0.616297545849198,0.8323502708043897],"highPrecisionTau":["0.616297545849198000223727245630090527811461478258785490032948428500765279474893180793669937553740487203299872788084884695","0.832350270804389687504374760144077086769656729734367824535907507345251211775754814164309366890250462773258845983963362625"]},"5,5,6":{"tau":[0.5446921452017452,0.8376818696909301],"highPrecisionTau":["0.54469214520174528634264365426731689243431172541797396349285740179632467689970881599368495830827763762094988812572595124","0.837681869690930039588242389881368553384718565261255326513615463394316965523746800115704007302361305632050932951390600944"]},"5,6,5":{"tau":[0.4734923400174731,0.8386691933312667],"highPrecisionTau":["0.473492340017473094195019158399661404707699456018603600525346445884537656484869105512192766854850363911861304229468351356","0.838669193331266683528904451688849541006949291627315842160144683715983002711904208293571566788744712054729881713976041818"]},"5,7,4":{"tau":[0.40175849546070114,0.8353397142299783],"highPrecisionTau":["0.401758495460701135540337692886660263238608435547071105470522174651351473641702502596807427845791440558350097273595521412","0.835339714229978313722603773749178880175850338528387819568312119258477324072989960272088572657062013297575608735743731318"]},"5,3,8":{"tau":[0.6893505695755464,0.8224744910387468],"highPrecisionTau":["0.689350569575546448148720839306688664510347918714831281820802495619191960126565494197165729795976629114234982101776154088","0.822474491038746819128055599000809660899338819310662037892876389940958438710513575599142940332372045340645082600670853474"]},"5,11,0":{"tau":[0.0763569617194207,0.7725126887511905],"highPrecisionTau":["0.0763569617194207114560007610129757179224793262994512038147608181249895375678677783799994643457141857423091940216203951861","0.772512688751190513190276475154638467128703971225280674980103523599139995074793239589572867287480045129909242510462043327"]},"5,8,3":{"tau":[0.32843872982805045,0.8275590971882855],"highPrecisionTau":["0.328438729828050430197073069427997994200363374016577656857255629123948121341562851264361603857803970404285393513295240559","0.82755909718828541174821659203763643140034511276395142074326660162938559877699862068394733844268879773571530142914001083"]},"5,10,1":{"tau":[0.1701195389226169,0.7970547482097198],"highPrecisionTau":["0.170119538922616891829268264627024320511317382497493703203564909945779780638835780835768055219188210306213246443268913268","0.797054748209719756242393627307574307803645761809979708759663947278946975733039485680715640494488942116709042962101948347"]},"6,0,10":{"tau":[0.9361822876144416,0.8278659260053461],"highPrecisionTau":["0.93618228761444151120618055633767111050976635506184342338478097771523801354499743303700313999720572638511598598228578355","0.827865926005346119377321272798080496365042158216740676359629673480791586406012882982593198183265696354047322025606500826"]},"6,2,8":{"tau":[0.7501808964711224,0.8766551928831733],"highPrecisionTau":["0.750180896471122403887047759163486984046379942629775762723781184253963320516712500487316239029841030363124068610306853524","0.876655192883173307389740715232428741047934649863745857698826282286289541200357521119435766346075418403988277771345888109"]},"5,9,2":{"tau":[0.25207410531467406,0.8150014981431876],"highPrecisionTau":["0.252074105314674079852483944466959939249084919581309370790495613075881772916757273821981534848304155603646169705094092962","0.815001498143187564473100568028148531128813459921532486494230400554360437912406355352792072541026532047903518304847869501"]},"6,1,9":{"tau":[0.8376473902285634,0.8565297373366757],"highPrecisionTau":["0.837647390228563369835393276464436586434384881638463331212842030680666971140634784475019054125138207090708526598932846478","0.85652973733667575571798140196816515588277062603787602509790958178124840983176107643127003717871532823352112088266297494"]},"6,4,6":{"tau":[0.5881743130174233,0.8990777453124434],"highPrecisionTau":["0.588174313017423333883164007218750609470125211221278656525994651374284408410931459574779409404516329487613173797447431729","0.899077745312443479786066878813757046325207605244852807640392642925178190128795696550626595009258290365852342570200555296"]},"6,3,7":{"tau":[0.667863586963237,0.8905493700807333],"highPrecisionTau":["0.667863586963237030198508751300885383568184195437042308380293644729349332676523332261654895169639899512040980193111177628","0.890549370080733309574902352408554006622413159926903632770776848247036126560853627995733257188495678814892705636070712315"]},"6,5,5":{"tau":[0.509590833210253,0.9026415553168862],"highPrecisionTau":["0.509590833210253018006516199156440102186020164492705805061000195139438285208687458701166431369507885632424430962941476918","0.902641555316886142405564062844713181987818835307636621997069109331136491182015582515174857443783457863213806770804080846"]},"6,7,3":{"tau":[0.350961417733907,0.8952392159611482],"highPrecisionTau":["0.350961417733906986138182237375033646786614784301988888509515808365655069543939033303976119593048242274634919119510113993","0.895239215961148149337593170594817387651072174705145622460432729060067027733801417627522698871015939565107008539563963495"]},"6,6,4":{"tau":[0.4309202063903544,0.9013828738897268],"highPrecisionTau":["0.430920206390354409351190913806978677851780494722359414484902287727248327155840731839960943682889182688381346211533369701","0.901382873889726868528169739681953439040209725581628111028115221326258452786357694896462094170293553079457123015873500736"]},"6,8,2":{"tau":[0.26816249796693886,0.8839300010996212],"highPrecisionTau":["0.268162497966938837142459510225292436348390060808475898599271353103211116940251992549867429270526181024640614539926328844","0.883930001099621207003330226994044052749326229580949048470465820376907973648367108850766647718243428079752045125979013882"]},"6,10,0":{"tau":[0.0799833716936036,0.8427760285289642],"highPrecisionTau":["0.0799833716936035996104539559021651147619104674127834003379906781478089798539595490082883165119863636497507513750780175009","0.842776028528964110078669993208805434198321064792194871232268255438984954103665421765580766653880956798148337524716192945"]},"6,9,1":{"tau":[0.1799214439969883,0.8668616305447521],"highPrecisionTau":["0.179921443996988305789556499992763654648805441389165525199804999758913622590088227278770635817531887612499710838742875967","0.866861630544752041266545354832520563031180312826396962534496909007419392749278599628478856172358697402012431316400672654"]},"7,0,9":{"tau":[0.932209856982397,0.9004835521633886],"highPrecisionTau":["0.932209856982396983064974025152178888326766360391338412357099790516463919681899268246580517260466811694978526931165278013","0.900483552163388616616779244984884078987283530759517647361452607069272796013680086683477357740170065993737236030988881859"]},"7,1,8":{"tau":[0.8262999454023061,0.9295215906001116],"highPrecisionTau":["0.826299945402306066488295757686228935898580676119433828702688242558764301731024095528895365418800959903597495169364829522","0.929521590600111674957345281829474872078593238137535494007918562218950695718176085907799462825434319035810845564540334961"]},"7,2,7":{"tau":[0.7313879662896802,0.9492549788546624],"highPrecisionTau":["0.731387966289680184993871760995964536801017951516543413847250833035588779678463731644201772878359776567743029252855496653","0.949254978854662413749196585563233095031831599233332655909282392254505488377586108899259372451865303173967179424496645897"]},"7,6,3":{"tau":[0.3787975443212932,0.9663048279010922],"highPrecisionTau":["0.378797544321293212937204741420019312402330410704513486702969042152899808423336468813615148078139610587376033391882635537","0.966304827901092178879888001912997452185333548470970754190174937737612111760047822658781175445648340057830557273073013012"]},"7,5,4":{"tau":[0.4667054568366157,0.9704721568702986],"highPrecisionTau":["0.466705456836615692823515509037821245594339354350720066332216272119760413567093880766101457013945720100286920417849801897","0.970472156870298605957450028451694396285815167537356010574393416765121925086665723919068226911479549512239271985897444044"]},"7,4,5":{"tau":[0.5537566044033116,0.9691070882444152],"highPrecisionTau":["0.553756604403311600785772116756268958902489975524625818923971692905621278835777810070806961039838607252173882184328277665","0.969107088244415199330102430044559116254965016847060512377549359145558261185715066227340074656483706619515901120430224241"]},"7,3,6":{"tau":[0.6414205843612086,0.9621555568149527],"highPrecisionTau":["0.641420584361208634934106524822948311519176160088359153123062493275459977086730815206969179698733621245054365498638245432","0.962155556814952639206382647550715284659447475171130405939734722406764000331595043320745423134462850546912050670878639315"]},"7,8,1":{"tau":[0.19251170103085552,0.9401156999825755],"highPrecisionTau":["0.192511701030855513844343230518232167623307539960749907512336135225098714371112495596915239045498958593320499995048023794","0.940115699982575469915609696639205269585937423979409921592626378692299228869597976045717621080816860767865943568946567775"]},"7,7,2":{"tau":[0.28831138934989414,0.9563823806874152],"highPrecisionTau":["0.28831138934989413412582980861112554737501119338662320120401570690765748370691880639137522045557611773822825064202617546","0.956382380687415157336728146056394821588075649965332080016368632687501304639192215467464275596234098619096671926549842475"]},"7,9,0":{"tau":[0.08495697805690386,0.916143163401304],"highPrecisionTau":["0.0849569780569038559683864280622711748754500203011884645453962017475858922350727558808850290516304848562500815659086927881","0.91614316340130401952616614808587273052690556455115532399407051455197043298013437936471359131793624568182687401746956562"]},"8,2,6":{"tau":[0.7079600618573457,1.0268142325965934],"highPrecisionTau":["0.707960061857345725398145027157682720613108043100022052021282558099292612067969213623916591288405519711926313649418363975","1.02681423259659338915334401044472509662207549910629861191910693041012624948608998667021399510889755071903210023828014195"]},"8,0,8":{"tau":[0.9270768474165864,0.9775459889415128],"highPrecisionTau":["0.927076847416586351674520556255499458519323414479029666609281673662687304804313437660423954753698399004635744206124285604","0.977545988941512736484153276402653476843099332471652934478185635568178890926256555541032630495513722904802301955970385612"]},"8,3,5":{"tau":[0.6086148931994174,1.0385368562552666],"highPrecisionTau":["0.608614893199417459172876343155194303232832903402084354888019880385784053078338885254414294837166738747331953494165285618","1.03853685625526663833211934945114966099267639374532875844327761413614721697157185413780986335613940066670605700829256826"]},"8,1,7":{"tau":[0.8120027373941217,1.0073852206395777],"highPrecisionTau":["0.812002737394121695784453781239853539542581034080310982487584902975161538180833859496090973443439206115904788082959276644","1.00738522063957771891213096464807827744236509192948494466400692469006381439502545750678840099778598535700858078591476987"]},"8,6,2":{"tau":[0.3135408325391782,1.0336805525318662],"highPrecisionTau":["0.313540832539178215668575157293341707647229544865434445427733665763780200614531954878970517178975358938191722990047242373","1.03368055253186623011950475084323369603073249883466960306846534998849447235620545606099972820827705998883553569816608368"]},"8,5,3":{"tau":[0.41341962022453005,1.0418819534066144],"highPrecisionTau":["0.41341962022453005119994990585947406587519143362098433742672805366176161766231251838674410675106473225079503093324998165","1.04188195340661438411762019467131530539746129593008785898646217030003997454628148974949705034390322226195653938138224408"]},"8,4,4":{"tau":[0.5111019351048782,1.0434618622178717],"highPrecisionTau":["0.511101935104878270498081127193189540834803179667187211122124127114028504622063299049298602691844304725673156078491330273","1.04346186221787162746479331435987120358892664488652323910883082428180571548137466032799436157737650583751028909626173262"]},"8,7,1":{"tau":[0.20849329928978416,1.0182784709846617],"highPrecisionTau":["0.208493299289784175512618152344306483662087741568726247674826159202454251983335888160366135358885743719433999513919188371","1.01827847098466170267377793981279006862688593122494813118336804490310734186556910719849880473374089477684736531301006504"]},"8,8,0":{"tau":[0.09149096823644631,0.994168314553956],"highPrecisionTau":["0.0914909682364463067057087242691246960987718795801541096615274742290747275127669796525707083626411682278655352461889012802","0.994168314553956097422857447560308414642252242252985909155467190897931375600294347869642818203192291680469274877219943804"]},"9,1,6":{"tau":[0.7938792958813226,1.091911630347394],"highPrecisionTau":["0.793879295881322589914159464308066309712710568581396309717186358452434691147521070211125470222580523009381182575544696352","1.09191163034739395609208136503412973469630146858098636775167370553518810349875843491052104193915126888896285359669933195"]},"9,0,7":{"tau":[0.9204853873751433,1.0608528564282176],"highPrecisionTau":["0.920485387375143261931839162056977595001673617197535952510219634838335443583937530035427266940533587033330019242015121202","1.06085285642821760921573683801505047443171882443559326388125301031023884519194249334442254960277024357777939618810364722"]},"9,2,5":{"tau":[0.678395237042239,1.1110038001636702],"highPrecisionTau":["0.678395237042239056047109831150238315182189506445981228986771382236704537064782445646074715903075627045422665045921076566","1.11100380016367009130706160913971342300698497473272101421613375462253219119488038483828968497054258809656335394448081805"]},"9,3,4":{"tau":[0.5672374028546672,1.121133868474373],"highPrecisionTau":["0.567237402854667149201315092946836383023145118489517703317090674294175928621503129539166190538599224503331108163613432981","1.1211338684743730002542303402427906152516246818861261592280199996474324239790099696940774945298254697948571610327266037"]},"9,4,3":{"tau":[0.45714088515433465,1.1232158103302436],"highPrecisionTau":["0.457140885154334644373523188806423634375394624414704682684386369275316008905067172047216557277425328287605075572240009031","1.12321581033024370003153403037719432541081924988305666951539474831291993021360326474964803509501559491826434112795016993"]},"9,5,2":{"tau":[0.345479437151131,1.1173785366528126],"highPrecisionTau":["0.34547943715113097765906680568123318365260744154300915802414066340010868770889099621703756133115232899911784442145541333","1.11737853665281252101438026279000888191348952768872749731864708920137447280953453043443464197023351843851293467362297616"]},"9,7,0":{"tau":[0.09999083090085949,1.078680704956053],"highPrecisionTau":["0.0999908309008594814319432841118928426767911709238488749882941895166130183458506846290352189212673798777238470600245485294","1.07868070495605287274461785595671348976702811238280913145820915546682624736941990789389207244718421911553388713906947147"]},"10,0,6":{"tau":[0.9119378659127046,1.1527150712358742],"highPrecisionTau":["0.911937865912704688987459937248762974038026353696612339062021801869060691988948259913689287717492956502211984921996757161","1.15271507123587427256278050337859474301047523469052631882040914001890336078639729260889086312717350976582497163465547788"]},"9,6,1":{"tau":[0.22887221619383238,1.1030960093325055],"highPrecisionTau":["0.228872216193832395375327472974283904127996257265061753307901446745608450166636553190410434933691788540638653826757292414","1.10309600933250559577727609415349018292375271457363825277355954851544843663029666815254644301739766579764802951686221521"]},"10,1,5":{"tau":[0.7705087792761596,1.1854284865041702],"highPrecisionTau":["0.770508779276159591906684779984197883963030159854131590005539678424340678641347960525085224489307641191095637924102813807","1.18542848650417025092522581447348634171174426686592462221917583046561136175419381787613367253123347859696429797648362468"]},"10,2,4":{"tau":[0.6402726848398136,1.20395645536953],"highPrecisionTau":["0.640272684839813585200211966088460820520489419695797399467289641623367976526234803441224731133814673330490243742067625847","1.20395645536953006962295260994189886357252597694215637211819527539770052714168758569229250970527462415386692047614941025"]},"10,4,2":{"tau":[0.3867497545807647,1.2093971477127803],"highPrecisionTau":["0.386749754580764688595869437369520738241182126017846609004053809880939889787116830590031988922524428274528925790071722067","1.20939714771278029202958826588503444119494800595562972290899618275040457795110080776121354572382719131883810291022210016"]},"10,3,3":{"tau":[0.5136981481305009,1.2116697372412273],"highPrecisionTau":["0.513698148130500808820509836914557505582463372204397782557292328797649362344944158129061681539402789692955660291993913752","1.2116697372412273091477388450754023720303432650631939959065444406399888077509496603779504985883926304268470225934127993"]},"10,6,0":{"tau":[0.1111432349695354,1.1720407541682374],"highPrecisionTau":["0.111143234969535405925031391337575615083990808480056502400295072744335946691096618546857578435146702400533417109137992491","1.17204075416823747213195282614857224838517542566594825765182795051227131286105424952560830357865495534190881420919871765"]},"10,5,1":{"tau":[0.25528663209599917,1.1968137707869133],"highPrecisionTau":["0.255286632095999156448121391877214028661981416218483003086817240877318932770429271438394933264067585607278053291333519072","1.19681377078691323284336073119230371683128657874017298931552618161807228787532582735042625260838038636777692768810062799"]},"11,1,4":{"tau":[0.7395087988185464,1.291174672160155],"highPrecisionTau":["0.739508798818546416167797404795367666077817927369792615515607420367779197306038187521604122863532389756961081251699359167","1.29117467216015486424059512580733203489225070577797018656986373638865025115007446660771608111646228566286911035699684949"]},"11,3,2":{"tau":[0.44176725239466524,1.3121625682565354],"highPrecisionTau":["0.441767252394665228190511169493012199900369287481595003150427417850479061311376473855429034523197842534710996957593843113","1.31216256825653545775574574289111868886594462973093019403348053466480539654118578810413250766397968972677717115290625125"]},"11,0,5":{"tau":[0.9005982414807026,1.256350039862179],"highPrecisionTau":["0.900598241480702523787198086181208493079366485533309999532386269416727746969032534093979592988880120986178389249062481277","1.25635003986217888954421638974441329241904368585688337866140502934472553604399174167134060979620487778713641497959760161"]},"11,2,3":{"tau":[0.589493931733405,1.308525380451205],"highPrecisionTau":["0.589493931733405040660118748242981909375493345002641834955730862322100106353363199061054458645736210105703237905276696994","1.30852538045120497024335737016237149602526225513529994396166209637137439245693273415594692957170389724144840050808084613"]},"11,5,0":{"tau":[0.12611520223152095,1.277536952806415],"highPrecisionTau":["0.126115202231520944136374875781613295302222479512622362478620039865597788019208790167424050931895730757051266039029591021","1.27753695280641512848003066290773258864728358204250707081225587578049561361180412521085415587858979354622418836719945163"]},"11,4,1":{"tau":[0.29048729822107366,1.3024955057272134],"highPrecisionTau":["0.290487298221073637080376452472021296417768101552617212097275436397857580415464487099241920052646819400189962069779573864","1.30249550572721343458228952347002868691804960213726768710921055635290194528882319634636261133037219615508037164979639639"]},"12,0,4":{"tau":[0.8849861883066501,1.3765988704249696],"highPrecisionTau":["0.88498618830665015524914778365866052539250351491467225669864591289050250153713438007345744889252683046648094405690698177","1.37659887042496974889830761871637952946184086711906320758850102264558537223552629000977352915036624029204637494225551282"]},"12,3,1":{"tau":[0.33939780112959966,1.4245343022779933],"highPrecisionTau":["0.339397801129599641193368852581161884058165506912882190838714471578422942002994464892573498365463447171147531906630855017","1.42453430227799328818531236978020841654943789568191127193103365786270837655910832270512955153746077787388931016890348087"]},"12,1,3":{"tau":[0.6966079046875696,1.413959444874052],"highPrecisionTau":["0.696607904687569558508620664247074447546503487355990191410574027931819736383416998343790352617174967127419499683268414324","1.41395944487405192199037651617012218812692113387151001973748007747735601223360246085886421139932722791824842099953467446"]},"12,4,0":{"tau":[0.14699763936594412,1.4000990735847432],"highPrecisionTau":["0.146997639365944133976895554609008137799169615696460284293258403829171753050008082820901080402761739757520428318200601896","1.40009907358474312985295954345702821295759904984711794055162811576931920494724562634786242893423056001359161235978850031"]},"12,2,2":{"tau":[0.5185374288694176,1.4286284548901063],"highPrecisionTau":["0.518537428869417509531966151270240920761627672449407065359604379479894363478996056328560174291046195598886359976708915034","1.42862845489010635180096136180646414327991568832721425474951497970152345369220495096056598235506277300913507527232285141"]},"13,0,3":{"tau":[0.8622326867624782,1.5214125464245447],"highPrecisionTau":["0.862232686762478187987799639895966942050164749004416339118367911879155644341833500172775398661586981087413121923432276196","1.52141254642454473439399917419150029123358808364168166865749880428870802884968470817976369094665912440772721179082977431"]},"13,1,2":{"tau":[0.6333112551559541,1.5614259866186182],"highPrecisionTau":["0.633311255155954109883304144430722701832275778925462428615241433083068932839475705727288941849453846326515186313719341563","1.56142598661861824092058682669705853970573653420326364453660946020293239833325250737077029477271227338388129453733071026"]},"13,2,1":{"tau":[0.4117379197686188,1.5694500992605074],"highPrecisionTau":["0.411737919768618791002315195178514881872710356236162942049592883289660993574461316955152808321580526948629860900514461341","1.5694500992605072833691618202038359246773893069092729116995094430747695121163128261768233131091349221563938424382811577"]},"13,3,0":{"tau":[0.17789498793935168,1.54774045050484],"highPrecisionTau":["0.177894987939351696146160114355905100175212987633301269190290574212850994101615390256865893846363770362124659297079512033","1.54774045050484002164001659597537982858234077033123637242979473569432847466739325713230800083121053213941349769822987953"]},"14,0,2":{"tau":[0.8259382540273967,1.7054288345734199],"highPrecisionTau":["0.825938254027396698851365803262238112283665674166475771911323376668943990599536164815368462385094611066493932902848173986","1.70542883457341985194470560515133744505083609881442137416881960947465902516440224611166609969819707689742119341128320011"]},"14,1,1":{"tau":[0.5297527114832316,1.7465687884704653],"highPrecisionTau":["0.52975271148323153765920108094761090328088730899154892450643882085165443458650492664230904847748659480524084001453875426","1.74656878847046530681852657504308843100078801126577271716454366911345898851735559895160022286547534838629218026314848762"]},"14,2,0":{"tau":[0.2280823736063089,1.7348569280713209],"highPrecisionTau":["0.228082373606308894220877443511677123739558196050913043998916236606779692646719284012407509302869419002424164658950653144","1.7348569280713209352213853980154833001059212529480416863874977734112242859040309231876509563186598920286274743750709492"]},"15,1,0":{"tau":[0.32385647112171834,1.9908518557148498],"highPrecisionTau":["0.323856471121718360878883705959192592518791820226401527947962072512761369840404736669556362095303189214136543683757839163","1.99085185571484986001499175029456593105708534744794140170723025429479330890423031880822771298784566909952253798721398286"]},"15,0,1":{"tau":[0.7580993199583882,1.9605862979761244],"highPrecisionTau":["0.758099319958388174585651075914529265846890726891953954289682037562860743268440832747602834728416586326528872609513240198","1.96058629797612432567709531637075515207638349095438893522452326928446896609718678803109223678948889048871785129907171818"]},"16,0,0":{"tau":[0.5781270160661567,2.380068204250201],"highPrecisionTau":["0.578127016066156707005232425673856806089134526714462406800184016390225289679117162400802793619148849576446200533443861231","2.38006820425020088505351132476027689636354352212204353340323882049023122220184500864717888756331459071926634775657891568"]},"0,0,17":{"tau":[0.9286786794290725,0.36732624907491296],"highPrecisionTau":["0.92867867942907245367306076963134716683639303814873511058070902701240834007168757768963780028865871907026665778876128371","0.367326249074912939998248273709880958705242820211622162028665359148765096336627330385283734204212068971570368547469867264"]},"0,1,16":{"tau":[0.854059048637775,0.4210081880004765],"highPrecisionTau":["0.854059048637775045096255364005865181257662524028931629627342604200988655767705287378758159327652518145304527209212212159","0.421008188000476505344340337156366204850005574201272448042951259491900592292983745332099171134679203537707746405647992891"]},"0,2,15":{"tau":[0.7973410646191416,0.4558262433824152],"highPrecisionTau":["0.797341064619141557940880244965352994253568490118638427347813059013262010797231559564537393674783921438663925958712302772","0.455826243382415152495341096665845119214410322680707196097203367820650001736613365176604902332630142354341474486146968104"]},"0,8,9":{"tau":[0.5301841429320102,0.5352177748130861],"highPrecisionTau":["0.530184142932010262246378529577542871130443098505173001477737810018645688739587309195239096247260468343361762964629222647","0.535217774813086089503703215170181923342221737816356267946822600848737135141279192899849971506776162374779060318497512218"]},"0,3,14":{"tau":[0.7474541638731521,0.48110181973597277],"highPrecisionTau":["0.747454163873152128613532382595613662285172447660344349222433577451998117505731268068293224998230249253259963119513329472","0.481101819735972777047618413476352734119899529737994757511338714830504160789308852231745960586810801240824713174333391333"]},"0,7,10":{"tau":[0.5718605400936111,0.5314881811716506],"highPrecisionTau":["0.571860540093611095904250631234524331987175024946744338611768606858237842527098930131637702194431819252325203539543293876","0.531488181171650683070949640513437934246001836860433408119806998612610594654225296601855729744912462050449002325011459311"]},"0,12,5":{"tau":[0.3619822411459284,0.5195358299018049],"highPrecisionTau":["0.361982241145928438039412679475816902477045150546009576508143296875193479943086201256465815717969696081165829963899813774","0.519535829901804934245135447355507660269550322444805710635513586898421797004107382910280646199325621612894566735056106369"]},"0,13,4":{"tau":[0.31782197227719383,0.5070526306799414],"highPrecisionTau":["0.317821972277193830127184085069754092365822358241550863861991254747122014780937690144832946462680425601241185596728694287","0.507052630679941438605257737439329256886955754316383447095807466685160402881877494830937449681410590383327997901306661822"]},"0,10,7":{"tau":[0.4469753560715659,0.5336449612593405],"highPrecisionTau":["0.446975356071565933048201878710088487537488255596586272448546381390347401688450500522680265397467300949489597986784787123","0.533644961259340477023501084760333711061393092516621698687879462644960776882967388704561688654203472565062530585839046266"]},"0,11,6":{"tau":[0.4048716665849437,0.5282582245333951],"highPrecisionTau":["0.404871666584943747599224010383438718052232018960901446819100347723199800571367155018555519486857437871846871743686198529","0.528258224533395104553923951471974324171017090930104287560979449009835640045706150916833405278447267907280458081224042497"]},"0,9,8":{"tau":[0.4886481686852377,0.535931545135289],"highPrecisionTau":["0.488648168685237709264840333044484310801556492908433789544945422051990132730054325193440374103747805944023337048023146507","0.535931545135288975937398385137146406691271802016873397463198369218910774972528543674009183207202703512976551060515879505"]},"0,16,1":{"tau":[0.16614008517418888,0.4359312814499069],"highPrecisionTau":["0.166140085174188868043365755216778134154191662698739839715780366666915715992660272151187864542674273436823532176218008878","0.435931281449906899166855637770899173998884972346295801452559122488316656068412663467231901387848073286465738984432010238"]},"0,14,3":{"tau":[0.271638171122321,0.49007387830306415],"highPrecisionTau":["0.271638171122321007277508865078170366152633215136393551484923192569229596303525798034534934115796231352110300661058248989","0.490073878303064151750619045678572304633312447509041255711394135885628272390622775809797382462177647860991874559859432997"]},"0,15,2":{"tau":[0.22207760329736564,0.4672721691660805],"highPrecisionTau":["0.22207760329736564237319992607267411086601717971132311040972737409928695890391548260155876426040951840799078887960404055","0.467272169166080514626747097044809567070705357593922105391566986395088932716266749097502677503602125580694679869310220542"]},"0,17,0":{"tau":[0.09393332943849202,0.3888023661997055],"highPrecisionTau":["0.093933329438492014505503607256628691246482863665427942962314188329038965408489220837486040673418924439707855959541116394","0.388802366199705488774919360829244599627647511494009651914891926939092680980460309821538047612294705490674763198860687671"]},"1,1,15":{"tau":[0.8637299702033779,0.5023619745065002],"highPrecisionTau":["0.863729970203377938692581076441727909749555226198419099381141708685807521334312858041747218291357922665356158528341325477","0.502361974506500170104513466419065728962234296648281752464333942664615255526575591478918292245560100977363381798742500238"]},"1,7,9":{"tau":[0.5545287597705815,0.5939068575372601],"highPrecisionTau":["0.554528759770581522149218998449042734135321495471793913153147130501934864058588511965281640791415467468891561361564870377","0.593906857537260180468256363288885939384261745249118078183087504606446440548905913094194080200527933770372761503685762772"]},"1,2,14":{"tau":[0.8024641495585804,0.5305426669204361],"highPrecisionTau":["0.802464149558580325506370614127398541170961377852937710229840943158671742954380076299978660633069288231432966183366390961","0.530542666920436110452135908535787272143621991128437417504054327557632058112560990038145387866480321009241973385183619847"]},"1,0,16":{"tau":[0.9402883569618049,0.4632263869403486],"highPrecisionTau":["0.940288356961804948308598778163863634634909414836345399323432693822998143261577376872922601246519806654415046000497056058","0.463226386940348627966495841501241799190428360277183037400733056855819023871425715173556277936433528665913863686922394678"]},"1,8,8":{"tau":[0.5083838934885024,0.5960611086697406],"highPrecisionTau":["0.50838389348850240156666711139156663695703101841793701873485157949048337950711363379012652894349673384133709745859150906","0.596061108669740619271795319312348731741968152109990952370500166940957889528160928062151000415106742500686971497909497878"]},"1,9,7":{"tau":[0.4622310455173699,0.595197581841303],"highPrecisionTau":["0.462231045517369853538869045431494356592096657283535870030019595417257575012500432622733357280152797238567409538060091139","0.595197581841303052122700944704713061359233250290736868067491944234085349390722180093157009636246906061088590076086472381"]},"1,10,6":{"tau":[0.4156972787777365,0.5912766170605784],"highPrecisionTau":["0.415697278777736485577397594356652298762276020720192599520690943041494755813659708718130713842147481462062436810778858057","0.591276617060578478465188378180665478897014141059836273945514084710413877481447881641466742876450307678426814095582633278"]},"1,12,4":{"tau":[0.319645998000365,0.5735226197235028],"highPrecisionTau":["0.319645998000364947788526435971864862240932613841621156902965333870109400544496482891923378557956899549133195383878636489","0.573522619723502808947514244412597106764281080293847543543563825360718492471175574409494368365802585166543440169213514247"]},"1,11,5":{"tau":[0.36835534512655144,0.5841480387393009],"highPrecisionTau":["0.36835534512655144734905501499488791414185231157143820093335393336842605250088069636179357145025942193347777135586475974","0.584148038739300904475931372300897965970884004808952809476995508087317967658986201992117467230637827586941603790635357167"]},"1,13,3":{"tau":[0.2687356308359957,0.5589084568379655],"highPrecisionTau":["0.26873563083599572713136591495424092166715670501748836154056837696317939910643933644083512402477750120674020805217324016","0.558908456837965475518358342729413760571838118953077896755097394411694573546914105162203230536133038874183465057098082215"]},"1,16,0":{"tau":[0.07692027704897683,0.4787226145125947],"highPrecisionTau":["0.0769202770489768300979691607962298471672778081975757998264954930470181463792362055701100428462572943491792717686592311009","0.47872261451259465911874027303151801854979035728517076405231850978484667782798195257278370539506312205978193739202885524"]},"1,15,1":{"tau":[0.1530645081421755,0.5137384031125072],"highPrecisionTau":["0.153064508142175521020658457451108357794911943652900084703742316445028078685016235015935915601740520631436117784888575015","0.513738403112507134127180148903136855955252640569356027118231797609424902264264977240532510956354099950048479714717894011"]},"1,14,2":{"tau":[0.21419355803325968,0.5394750326138809],"highPrecisionTau":["0.214193558033259668790188371642487520961992926300096018994883545934566820293068248526376308859870254866688296461264142961","0.539475032613880845886511022841375465906729862454997175267965862361958552203358720333693395253779791895677124216341007237"]},"2,0,15":{"tau":[0.9438051647418576,0.5392240942401172],"highPrecisionTau":["0.943805164741857608921228873725959601593951128899738534996638347678221485022959275809132969803000404712749721777091770643","0.539224094240117200885025322889629132163531436566785600701192338015275340468387386498400877360975569830146510812674225669"]},"2,1,14":{"tau":[0.8658908299189761,0.5723939174732112],"highPrecisionTau":["0.865890829918976189636735620786936110437777338983417260209091892342140333811691256498563663556324260293792681029735516939","0.572393917473211187077711472673273226424965312675908536122266133318956185677374115211272766444622937883615283589314445035"]},"2,3,12":{"tau":[0.7429767887014412,0.6159824929298461],"highPrecisionTau":["0.742976788701441220808033454995378046493649644242594921934289369716237320145398480893691482166840663101768203554063533005","0.615982492929846097839988099566270084581418255558025139122895762005015290541514080523108130392324840068166518946628808465"]},"2,2,13":{"tau":[0.8013150586208496,0.5970757627817096],"highPrecisionTau":["0.801315058620849644078714910098475174259419966094978018435535302177444254382767234094940118904165649698632297613909222684","0.597075762781709607059197571427518681732549838452153993848574948539024078846514274468147816798800112935451878482999722365"]},"2,5,10":{"tau":[0.635473197568748,0.6408033095772973],"highPrecisionTau":["0.635473197568747942488117116973818243062544721080772827948837969252085727946415917010945516702061242659758444297333158061","0.640803309577297202569457905342483385349374343435413592407666562913732859188342535896625926379693902283798001218773689428"]},"2,6,9":{"tau":[0.5840694470921156,0.6478063978280714],"highPrecisionTau":["0.584069447092115638498769899295615112891678727697078263880670978832216026895702448857771716774618669548418764753056755691","0.647806397828071394055071791468136900375062231941797297404647709187908288938540101390195565243677845397280676710627147937"]},"2,8,7":{"tau":[0.482842096915679,0.6522267319622157],"highPrecisionTau":["0.482842096915678993312000191049404404286794714393759897815999419481382235252039442444435932389215079378730106028303182667","0.652226731962215733238069615999428550548626541342136813324857764928674027512935720752506900976084681114534920441139315557"]},"2,7,8":{"tau":[0.5333524581880236,0.6515763031954882],"highPrecisionTau":["0.53335245818802362343295166531372888657780108099777504736592475609718838313291252697898386100956680944565809607795453028","0.651576303195488153689026986411480825118389811154268370094236327149518459937263839287350880513663342463124498732819215469"]},"2,9,6":{"tau":[0.43209229003188676,0.6497693132852977],"highPrecisionTau":["0.432092290031886737733609854455842775158406243295969564111231404907430253018085272739885272326163823192944937354459111979","0.649769313285297751400241609969058480494981814599901787622129164367709785749689379744194698561790493663862675781080735569"]},"2,10,5":{"tau":[0.38062346299436206,0.6441158698021008],"highPrecisionTau":["0.380623462994362034395181099270985963759250183513357953150362147813025341190594916375777096225964804649150221509834435","0.644115869802100762887468719219189979989426914427246596670956987114058970321726582810660066718978757463706691298715493939"]},"2,12,3":{"tau":[0.2728998581652332,0.6222590123016707],"highPrecisionTau":["0.272899858165233184719683182180392605728676640635790078741752177205049849824985001356410442329365587435160909712144567494","0.622259012301670656209288142920601554618416780534491089468934615006936791213829541189817365689499762505694117233337219818"]},"2,11,4":{"tau":[0.32784108169534415,0.6350638413304407],"highPrecisionTau":["0.327841081695344175968529081041377433763059788539770504958686643933413748519924679621494565959295326655087926994400976159","0.635063841330440735135422195628024253111487345271350490899424721093893266936957034196885287863116439550212275478439350732"]},"2,14,1":{"tau":[0.14967228362339197,0.5826670116020192],"highPrecisionTau":["0.149672283623391958928927668800565360171288037396324630918564129685067331285929107633679873329702774065317025128390466911","0.582667011602019214591122639999845978110107619445352400784903636820883191071193792458001165883220444583111244415667529888"]},"2,13,2":{"tau":[0.21441435115702806,0.6051181671808262],"highPrecisionTau":["0.214414351157028070408467792736839725163467032315252832351345460546524917947231782555693793031753256905510608374120262383","0.605118167180826232335096713506628159218310619032905796911277611678910840532614226923872025980011642849055625393693333544"]},"2,15,0":{"tau":[0.0714947093740834,0.5531741619481674],"highPrecisionTau":["0.0714947093740834025313308227563729988834276326420819552697054381540470467585175433336732818951324902458398819133224219984","0.553174161948167331397491954238826834629029438298562778286645819664631476950924043953342028717181440085512371698838529364"]},"3,6,8":{"tau":[0.5632229300191659,0.7052860376641538],"highPrecisionTau":["0.563222930019165955319962212212744288682920347161895895785850934520115519340067354389970878295647113506645741864374636341","0.705286037664153814624474917725654152260127238997376934822969738196756734924578719594656684477987479726244024286872498389"]},"3,5,9":{"tau":[0.6188767605409428,0.6996641601750694],"highPrecisionTau":["0.618876760540942744559337776114631943419373085166737018889392542718104788437644450899229475408870604547699908075260300499","0.699664160175069353541226919666984535387361198650700252892929499413474910185928638999771903113589979104984484488629438947"]},"3,8,6":{"tau":[0.45291851437802577,0.7066703374126422],"highPrecisionTau":["0.452918514378025756457579284173258717284075470062752648959340242162523150801947709966432266021387660435726362930834559865","0.706670337412642163948243037563591743729440386326161997130091082996031529065848798903401435648765934238154583026248388101"]},"3,9,5":{"tau":[0.3971877795874781,0.7024595315079066],"highPrecisionTau":["0.397187779587478071966740696499095230808651418524052847905028943141056852452698777782599706656095323107678720972983939465","0.702459531507906638789161083228700492862792825723423697084308244522809557019718781580405059162444131762992100022306557986"]},"3,7,7":{"tau":[0.5080834467413612,0.7076003241353223],"highPrecisionTau":["0.508083446741361216162328322011177375917362249825835317084177994864391729661196517384364395940317684864660997179886510898","0.707600324135322294875868375046278325394884343775316250866636282594715214498241796671796876408458737767260734235841096637"]},"3,12,2":{"tau":[0.21890735949135476,0.6680309869575664],"highPrecisionTau":["0.218907359491354763107748081463188005324212129024513509378253266626993638336398589100433973215201764617369123482746538143","0.668030986957566380130361611686413223213346378307366614498430749486827427258467250710624990142398246150444177555809035242"]},"3,10,4":{"tau":[0.3402604120285516,0.6948263159687679],"highPrecisionTau":["0.340260412028551628628685736619839945073491818837553494809087263618554134180249405799954723658044382698478178978509713465","0.694826315968767901116105515669974167884555290663723177339027325490174000450373701731583472269862577174625986319267186383"]},"3,11,3":{"tau":[0.28128195130961425,0.6835004641226299],"highPrecisionTau":["0.281281951309614245026701643234449829617452708243610503182987391404119076418667970844425044230711056536825075863403926259","0.683500464122629826642107807882171005861087769685779769962135227564824375410901923961687750744063700591227849067265441665"]},"3,14,0":{"tau":[0.06989708610141422,0.6211520044551148],"highPrecisionTau":["0.0698970861014142244072823206065290709115056638144264323561295192670107822604781834610052743400973033130321112906907377656","0.621152004455114788119819109351742576272780416799747777945377013643954875645407130263564256822248398324444820938291844438"]},"3,13,1":{"tau":[0.15060021656737488,0.6476772985887815],"highPrecisionTau":["0.150600216567374871452927970613184929024853737837894681775725595756985045675018660650622229379858333849286489608090949383","0.647677298588781569027846843823223307605128396504856799089298069182356427564569406491364706879851046683215656464139247968"]},"4,3,10":{"tau":[0.7228995624619816,0.7388177779182287],"highPrecisionTau":["0.722899562461981622528857097229898554087310101791267285331030216907611144299382785838192357514724205023524239778663730375","0.738817777918228613388149973948125100907045920104501598611335876572546732496712756623594518710528152485820180351186061448"]},"4,5,8":{"tau":[0.5985432549855003,0.7583495994622269],"highPrecisionTau":["0.598543254985500357934511025561747881746100304198238132372730242411978744747314855475949615196208606276578182034660136559","0.758349599462226925761082209271665942572152734123845265946284178928796846946315356539123484512910401708474513226837334297"]},"4,6,7":{"tau":[0.5382380600200394,0.7625805183981718],"highPrecisionTau":["0.538238060020039388962667495862300883983605344695018872019391737666830198886114615569106121356639630710018064668010279629","0.762580518398171729280108765938444310204918071484283274615284153719814315771231073015524185162564416207139914370322561961"]},"4,7,6":{"tau":[0.47821047455979543,0.7633369992256024],"highPrecisionTau":["0.478210474559795418604090842314319840528507492820294308466989909711498934186955406434995538586208663317816204232762541361","0.763336999225602383513689767671274947471346731366503363226095271859621260057069737488835226982955881962151488802218783397"]},"4,4,9":{"tau":[0.6598067749856117,0.7505170710884933],"highPrecisionTau":["0.659806774985611734907008429823344880215901838155085712556923334536954765807791381086625116399017392429103919533115851533","0.75051707108849329829778648264189009232816605976070303036853379872442639022707104672214110263868119126733672798075293647"]},"4,9,4":{"tau":[0.35642178337034974,0.7543826425340441],"highPrecisionTau":["0.356421783370349746728063692603704059755156213078187451422827052550200456171923250016556010667522969912308603650258529828","0.75438264253404412365765281618823537246939375041563386262128136281326112722145369693559659303950897079463975726889668652"]},"4,8,5":{"tau":[0.4178330308003808,0.7606348183774246],"highPrecisionTau":["0.417833030800380836380312770509680084663867578287270684635614531575464275586148842526131177242249883943063380308356550345","0.760634818377424607698195827840226196090273732849066797378112431481775594682434775274774777482993846298782816509556880339"]},"4,10,3":{"tau":[0.29309794234734987,0.7443673083773925],"highPrecisionTau":["0.293097942347349849960186964939331420547451860126849393875060937543412253650049100672252793732478141332336515055536608438","0.744367308377392478820384776921743561973469901970170922920153735067389547566363530333710462192222015242904018636243140096"]},"4,11,2":{"tau":[0.22653100160983572,0.7302134612165375],"highPrecisionTau":["0.226531001609835718401336852810017888994617461027805857855528172141498202850033757614907453635739932402856357157740770708","0.730213461216537568260783424409965748029391869790005131979237823132128083596130355071653087516071603740545449496453099463"]},"4,13,0":{"tau":[0.0702779385844411,0.6864782157945091],"highPrecisionTau":["0.0702779385844411105861135516436612905144358720352209795331784015804781021925409096116803963249732151410399012425109678748","0.686478215794509114697613021802802059210290973816677503602920953565915972620326966658337871359882397854007312152544880607"]},"4,12,1":{"tau":[0.15428412413916961,0.7112927543233428],"highPrecisionTau":["0.154284124139169612421758029972366905152399139202684393050193689610315309810453351577331127310222058893406642979420243188","0.711292754323342874684212468921779949019699677564094462244170207730001791264232965984045163007809959729425534442763038624"]},"5,3,9":{"tau":[0.7083005013649197,0.8007611177315443],"highPrecisionTau":["0.708300501364919696880546621700762946197532129878294130368396050710589999271302528936345392348574066807490613436366852573","0.800761117731544300125976405965388542558507549380848229470819771808559209399645537022045266675707448851544321477382843226"]},"5,5,7":{"tau":[0.5740790465201467,0.8179010866922533],"highPrecisionTau":["0.574079046520146748246104227037603796892164752651766164445715571666522049695750827953168661755540949711126858025817591173","0.817901086692253288479739283549832330760712767895287236532777817603358864446674946010002531819530635322692615842760416668"]},"5,4,8":{"tau":[0.640372429776983,0.8113549097131436],"highPrecisionTau":["0.640372429776983022520338857277829727211517172828473741856804437346187956316713364320787200771897640952298742713545633479","0.811354909713143586798022173395961911675137987421173181250032610840482934733590016031608992257896933403549222029841661696"]},"5,7,5":{"tau":[0.4428766221066927,0.8196035084128906],"highPrecisionTau":["0.442876622106692712384559810384022419684207999724132077222201999787853352013054643316878560637901004325550428850603345655","0.819603508412890691823193270964725464912099102261107023329883107266712558028095372654630450312086711791528072160779373698"]},"5,9,3":{"tau":[0.30824500501556096,0.8060501610360766],"highPrecisionTau":["0.308245005015560963087424558177223156914269702285073954036506833934171449483390346249607869707352777847247125631319136983","0.8060501610360766615912345646906564837420903781330277730277320548700193505125600612554838331902887453814564633474219681"]},"5,8,4":{"tau":[0.37643354554314046,0.8148016020104509],"highPrecisionTau":["0.376433545543140432745665256551175811327442234350000759689518528753789483797917998642240399392154862307874236056867734342","0.814801602010450942928603164102919796303736358499124911172967463005261086685953376066064559135417986915881029914073298816"]},"5,6,6":{"tau":[0.508502300575453,0.8206233450443328],"highPrecisionTau":["0.508502300575452945249897733951380828727906728913929396521255816848257128598579032007534902714579022016770890008511535537","0.820623345044332785686833174735290472930227071914986001405217243174393013111287871475250846143939077773741254029771076373"]},"5,10,2":{"tau":[0.2369640852236031,0.7930204613730886],"highPrecisionTau":["0.236964085223603101909992889924853072691539495848041095087860995378964343439594734915858575286648213368578140939537135702","0.793020461373088606478109096337693647100712268182070472654115523413626359494208163924299295049072358884171216266655030648"]},"5,12,0":{"tau":[0.0720109603076079,0.7513021827665493],"highPrecisionTau":["0.0720109603076078984932299173969811910502280446548714351256741818060289020669182862866574291795915212303036305163673942454","0.751302182766549212810885763051379537590332133926216455563787785176884067996748601908813862735405642070736210006726581779"]},"5,11,1":{"tau":[0.16018901152494042,0.7751354373013303],"highPrecisionTau":["0.160189011524940423288646430249310191788658955900635655761682761056978982307950120622312629249654020916944406031609449733","0.775135437301330251828987475706553988076946815171180756482066647859579418808801340601494667607081712354081690824672117655"]},"6,4,7":{"tau":[0.6168644816062893,0.8740112756744572],"highPrecisionTau":["0.61686448160628927911643508535154419502186218458231891397517619680170588634388302314144598931615930485972744941242702508","0.874011275674457244286134524960915874886913077751516581474126502307969233671402519168776991936074545273462899959088794432"]},"6,1,10":{"tau":[0.848024568879133,0.8304937415817207],"highPrecisionTau":["0.84802456887913304789037499166615388116015264515402013228753026932066063519103534223320650726863507397645568594898855079","0.830493741581720788416470371584252475645755175270491284191558797920194148056344995894164510951414459634937921424656255153"]},"6,3,8":{"tau":[0.6903689286122873,0.8645263796627695],"highPrecisionTau":["0.690368928612287286019171563868415924007803377814771646180076417647502882531507373612395105572244713762798166450471031271","0.864526379662769496957061305960747362385849440450920261514337827140257359678551185348336972604974219120094735167785187221"]},"6,5,6":{"tau":[0.5447523550649215,0.8791474618160311],"highPrecisionTau":["0.544752355064921500698815226222197387638561259139319883503373243135765054712921106143351075591859171766410713555407612869","0.879147461816031067400110230361455059595540435105202775976648379168146950097666840148338234979119343111636685275947841089"]},"6,2,9":{"tau":[0.766637840723928,0.850290484364859],"highPrecisionTau":["0.766637840723928024152876731459633347030695518121513273669297255362713208530991261090337754715739250372844882325069400079","0.850290484364859021621310721587440160515074059175347128272283944932200576091238460994010793745147406987585197344278776833"]},"6,6,5":{"tau":[0.47302802846103337,0.8801196778137831],"highPrecisionTau":["0.473028028461033386052642436174519968478295873768592373931931843290759154390680277266126260693639493694510970972842089122","0.880119677813783113080949773195928343734001587404326284104617146703616075349000033653398679615359353654075291785337212046"]},"6,7,4":{"tau":[0.40078296389649587,0.876954754128817],"highPrecisionTau":["0.400782963896495882170647702477526107887452438016433781286918445099830123574510065392400041341854920191950265259511316616","0.876954754128817012577979782517440689077996116021764740866125616942857853343442224718799230482003248784830963182003012957"]},"6,8,3":{"tau":[0.3270021216675372,0.8695320038833225],"highPrecisionTau":["0.327002121667537254089718092133188644792660420768132470642159372727564230766147922022529202414942114741128035035621380365","0.869532003883322563888398266787837126807828209838488697808530384205808719607125803197273187621168383057797755818799719219"]},"6,11,0":{"tau":[0.07488216401695037,0.8171385426314481],"highPrecisionTau":["0.0748821640169503722740452142454200448825933380931180783244860126037971916951656945897822736438375072766540868600465089635","0.817138542631448134160462610716165132364467589868540809302379921771233476615642009369857520567284948959661296278109586698"]},"6,9,2":{"tau":[0.25028593147806083,0.8575561084917834],"highPrecisionTau":["0.250285931478060804369135275693322278204615309760092574492884548970284267533687750270273641027340532214646034318103444643","0.857556108491783393511762746424410518505873195966294959212378771102100284597340656759329516643462661050095484543869900442"]},"7,1,9":{"tau":[0.8387030887114131,0.8982234179262738],"highPrecisionTau":["0.838703088711413126226055503567053040255707692493348869500775995154100924693298594912997854502691440397550122511148115609","0.898223417926273803890744215385074568688165378567751909927395779611469757985673586720672763221909676925725928295246591622"]},"6,10,1":{"tau":[0.1682103097155749,0.8404677109420212],"highPrecisionTau":["0.168210309715574907330997296935985029308125880126623001169242024509978531747270992696929887399050357870575792365805599623","0.840467710942021190063239189060522843746937474287829309558718841343282885948056381443696684339240529235485901683444969631"]},"7,2,8":{"tau":[0.7511304433463696,0.9176618185042872],"highPrecisionTau":["0.751130443346369571857148829405081500038681129804886698294291754843745201469354837246146266951100978172636988082482917698","0.917661818504287181329585635487883434543049299341300122952338752921087478637138425756413425495645717475918119795326554808"]},"7,3,7":{"tau":[0.6685160365550262,0.9310988725650912],"highPrecisionTau":["0.668516036555026198463970794424290949397791256050783179217006175379970369599524428793090534625577488730987440857806800284","0.931098872565091156222702494538952971405711119005170390610980386842858238012667602409761966890871677162467534379905239679"]},"7,4,6":{"tau":[0.5884379808410521,0.939361690750273],"highPrecisionTau":["0.588437980841052186797648192196808140513046661530061053602717546812104379910554420118552804271161474937851853199517826253","0.939361690750272972007558011830664290775035074324673603070997381528778256796997313483503764067074460485181923617644518016"]},"7,5,5":{"tau":[0.5094268292828673,0.9428312560151767],"highPrecisionTau":["0.509426829282867263453219572992214258627261324536535540278078176266689582623621021056975197547766783649348410489497250149","0.942831256015176724937684084125281528900168992183313804983675037359723323606281927798969342548792672832171045635327181092"]},"7,9,1":{"tau":[0.1785032465074875,0.9084410984301143],"highPrecisionTau":["0.178503246507487501675690647868696819319380898537818232770238487606916722519358371800720044209405825403628395516944243046","0.908441098430114310458747922704398415551470118103332101038155699088847527581095546079613964847895177793594879693006466652"]},"7,6,4":{"tau":[0.43032674160058154,0.9416429492989836],"highPrecisionTau":["0.43032674160058155478427655926094487421351269750139803750751459951321497059770018053322026744806854349861614976689182744","0.941642949298983662688082733389982009515030319224015328004032715279625540965638153859664770334345447231744047453014319837"]},"7,7,3":{"tau":[0.3499735515186631,0.9357397056495693],"highPrecisionTau":["0.349973551518663082009801586947163366385585909141867203680151038708760759535190506104052922390198480398958581128052677745","0.935739705649569370888009571680769427277212398782275512662710070918860453263365861560207833122455873309875015175179517186"]},"8,2,7":{"tau":[0.7320409090942162,0.9889448343775344],"highPrecisionTau":["0.732040909094216149884679741175159141174733764653259562356871141337320427232506639485439287163151769378904556930502109961","0.988944834377534483778381820726704738208707055140718987294552076710790546181618604004114932801653495898535876046131088897"]},"8,0,9":{"tau":[0.9327518933167438,0.9414608461479848],"highPrecisionTau":["0.932751893316743773141068839189415387716182082018902066933160204341032475904937277563413170986408004692878313498736011929","0.941460846147984785350053649392443975121198602017480961443156011448871130203843946750733038060300786688298120482956633687"]},"7,8,2":{"tau":[0.2668672329888535,0.9248591846418215],"highPrecisionTau":["0.26686723298885349823164061657977088963265010485991948130783108178042880900348498212276511539074920953037565482311954905","0.92485918464182151942025476336243802474922567099113830452214414451686753457014502479425945443103056248052768562003879627"]},"7,10,0":{"tau":[0.0788748155489978,0.8852746811748893],"highPrecisionTau":["0.0788748155489978018956197448973821936821953276910149932805293173442559753572103430301681149753918904569358221661561283058","0.88527468117488932095472048792774077765635045895792716561247212327853427100141740089548747997368685143668813576871575913"]},"8,1,8":{"tau":[0.8270712769373972,0.9697348897758672],"highPrecisionTau":["0.827071276937397281365793534951918559517219354195469820860294639225022532877687747568264797312930383610004544370964765326","0.969734889775867190477356593961112091961469645268163346181741758730273239658729472439734669916285354269281702057929542298"]},"8,5,4":{"tau":[0.466390256684215,1.0096561620400377],"highPrecisionTau":["0.466390256684215029441229965290479470290700384151416805737541380397673511231093404494252369382442113753355824072455306414","1.00965616204003768236620231092353401570727517829150964058787428063497502880270493740844309998142272034608614808465158588"]},"8,6,3":{"tau":[0.37814067906620263,1.0056338027031593],"highPrecisionTau":["0.378140679066202659443290694031414317285001554747193662614334099418258090174197454131885174778098587506017871564951107881","1.00563380270315937959999151676058456589951154709872762483808672678662215911967477455035933723840927711740050338288731731"]},"8,4,5":{"tau":[0.5538013258513007,1.0083006683337234],"highPrecisionTau":["0.553801325851300726277870980213582368285940051086917754980025155343732908979886892918101945336619486578421023334372970547","1.00830066833372338395358847925772475128483111620796110199079498996889818805562454304483786906622300891879852067468055174"]},"8,3,6":{"tau":[0.6418040107205197,1.0015138840828968],"highPrecisionTau":["0.641804010720519707039849750185911818324726215243516787451655415126617630682640913813935456049552693951367950669749750213","1.00151388408289683784872545881439724939701737977824938936000998904811478338200671253485449939579219758784161435622770079"]},"8,9,0":{"tau":[0.0841077458813312,0.9569864091901491],"highPrecisionTau":["0.0841077458813312005440736740984919137691691091172327820654030070280769595812330820308308772055500911323868577164614194285","0.956986409190149093344407688614644093891441929866887226609581462487973810972312997032956355423524957803669435554252875483"]},"8,7,2":{"tau":[0.2873786595796952,0.9960207680740386],"highPrecisionTau":["0.287378659579695209945007875011072998873280647393529759962401572100182048709236809513618945606075176647995553824769870767","0.996020768074038647857852682334155652541476498668078845273721412745612959154004498274677962906867418567652297880985703921"]},"8,8,1":{"tau":[0.19144982559477106,0.9802489505618905],"highPrecisionTau":["0.191449825594771066519442276848157930696440648500431193634202693909210340036668661548944248470380290452829226110060574759","0.980248950561890508684675841366810069440289215321675142736041355916215142436032748603891917728514071100816830038613755067"]},"9,2,6":{"tau":[0.7083891464255129,1.0654512965048926],"highPrecisionTau":["0.708389146425512823274066010567660794768797818666141172133744499151033731713170285572151814221803381986353655410599798455","1.06545129650489262423859305659464138884575818683496982893670498118148030488134648132671819926279375396279753656452018083"]},"9,1,7":{"tau":[0.8125605986270854,1.0464187330622325],"highPrecisionTau":["0.812560598627085354101026670992645555645214017653748550005850177811957534868229908295133759478410418447701231532446042721","1.04641873306223264297932484927436519531545922966321269475981104925090006200573361485189902540332114203325115928021200081"]},"9,0,8":{"tau":[0.9274876394745224,1.0171763077014981],"highPrecisionTau":["0.927487639474522380512069446894679311329682977874147947038261820916487481985090926228196527314011356652658762864827939321","1.01717630770149804256553451810651637012528451497341172034896666479918100193237428317890627183052141861366429550204676779"]},"9,5,3":{"tau":[0.4130120981541309,1.0802727522625821],"highPrecisionTau":["0.413012098154130921331114517225632715672239676266705430263113105810178233168258568488966927649579820133333842689434021693","1.08027275226258222512277018363383314775019209491760127802252175767737589538797636298957083156911172248216492840887971408"]},"9,4,4":{"tau":[0.5109899725200349,1.081788384786487],"highPrecisionTau":["0.510989972520034918040670241683240828057963636686461829162176627986766789225075271363969885681205836701396500675195199375","1.08178838478648717679331256841595582579655444782313971833451468522908948349385774973144465310498930110127550733835654689"]},"9,7,1":{"tau":[0.20769980502473712,1.0572557545160834],"highPrecisionTau":["0.207699805024737103148623924368354731535306875935784907784393370357469798109120778971165352530153915249892347436987704609","1.05725575451608354828015737852488596441249957263450290878460054824260085249380523425523201177748665168476490102294049668"]},"9,3,5":{"tau":[0.6087973086467239,1.0769436534228207],"highPrecisionTau":["0.60879730864672385226436862929848959588369362575641170223852281108636175770040563675829041332628182444765865245080645683","1.07694365342282075139153262627413160600146572018789803688573980263895563197076547748248406226119214079572695361135494678"]},"10,0,7":{"tau":[0.920795794595903,1.0994061990180184],"highPrecisionTau":["0.920795794595902944606958103147113664331633624369358457182887766993475299749061917305184180348000958850862766390399204259","1.0994061990180184931044577763583611998583597214567862228979878486233996063540264426311456008583065673551127949071651174"]},"10,2,5":{"tau":[0.6786523307483878,1.148796491506898],"highPrecisionTau":["0.67865233074838778581516022605090985095148022605328368741768678470234963996897787260210567148075180886049611286618257447","1.14879649150689789064532133714439434094704180418377781700706792835476634247955410326430950979264028156258506460580427524"]},"9,6,2":{"tau":[0.3128829281835048,1.0722840242757687],"highPrecisionTau":["0.312882928183504838963048010941622594000472863651978813098965194238665598025556009744746891246402331431844683788496853708","1.07228402427576866304679170480868659896326995733855232490221963101012886790755279469867172586523317251335864412557341373"]},"9,8,0":{"tau":[0.0908339346328407,1.0336993461544215],"highPrecisionTau":["0.0908339346328407077470695214593398639088046269611904225182701154477757093828917926141405769541776786721858575085402698611","1.03369934615442147872292019079631415092587584310282935430006922204480381813397571879878579186876322443116411539337965826"]},"10,4,3":{"tau":[0.4569225668144257,1.1608564024331984],"highPrecisionTau":["0.456922566814425749507957594321414336294036440617843924131681244587509644469323072196817456741658748461442544151532115717","1.16085640243319837556632558654664061421807784959416277624899455533205583462591010350398910030903831349075605839913125942"]},"10,1,6":{"tau":[0.7942724678828209,1.1299982851763066],"highPrecisionTau":["0.794272467882820935514705226299044650213055574012189063601251602812496715033032911529582750424898122047898871039216466288","1.12999828517630666737379034815090827321683783535923860368394618236414295006571203181299439080485039132645867114486522485"]},"10,5,2":{"tau":[0.34503433556169966,1.1551502972659937],"highPrecisionTau":["0.345034335561699673403327979078910675291139126832727734388141789849388368619364355410022075693889343795923932166347287923","1.15515029726599369846474931274513456404776628278581559278822322796354052388524077744957949029239568914289331483027431615"]},"10,3,4":{"tau":[0.5672696487158922,1.1587811413088127],"highPrecisionTau":["0.56726964871589225070164040595610888157464097479083607743605922908256755912056701563854017455857914537328391263372926417","1.15878114130881267057427010235835055557568202760378441595742779459960009372135776359434985716495632358856726644647632851"]},"10,6,1":{"tau":[0.2282868737203545,1.1411436132417034],"highPrecisionTau":["0.228286873720354512463161362047180213454695329063634299396499299969126412778962870505579886125322744946430139456632447026","1.14114361324170333938251456998737174850256502928631579276219696700240677702305580818100049863662349180231274233280538583"]},"10,7,0":{"tau":[0.09948159297513134,1.1171604648803581],"highPrecisionTau":["0.0994815929751313509133882650571286320632234048127113351843611872544037602375171761762154022771549279672429827465915266833","1.11716046488035823084142655278897894783171410966670251828822080426627491259483848653009093474705262141041213009551351254"]},"11,2,4":{"tau":[0.6403968551694723,1.2410751870212946],"highPrecisionTau":["0.640396855169472281939964052835119777555792775268013935774929372811750231599616260443906522917125819587866750733063184458","1.24107518702129467259320717783147416582010357521621770559189108194400371417007315625046778283537924923826268443725401547"]},"11,1,5":{"tau":[0.770772362834424,1.222755379497638],"highPrecisionTau":["0.770772362834423922951576397847574272014463029846889062581335038595320278688740812691681796410226002949864258144588466333","1.22275537949763800920905075240138168004056781043941067503089787252189616188436771897572734695516633295109141623083871121"]},"11,0,6":{"tau":[0.9121692623228751,1.1904035478032875],"highPrecisionTau":["0.912169262322875118762028308902617767631824526529182252884168084080353337164493084298847145004166201548911333586094012201","1.19040354780328743785243056510320886263366226279428493812054841980578147464416525586072305779008163345944714257061309499"]},"11,3,3":{"tau":[0.5136213754449526,1.248716848008788],"highPrecisionTau":["0.513621375444952604683060595722135606737838595389493158190225765704114159125759523646975747882498801709967054958325357934","1.24871684800878795916055441799124493928366541206877078025189589082283223460561668580263363847302890425214828306664422959"]},"11,6,0":{"tau":[0.11075117193257526,1.2096752372749962],"highPrecisionTau":["0.110751171932575257305549990689627074338156374982406176500601064796865175354136804936919876897721696298133189196706880718","1.20967523727499606321987132296323042517272185531732197837210707954617946129579749954652677365941652068387245837288976368"]},"11,4,2":{"tau":[0.3864712906638143,1.246504103547291],"highPrecisionTau":["0.386471290663814310209531495076954047237738980310017513503538007122664007703931820322732596016216517316275602084259343819","1.24650410354729092332073815210991574536777991576508913467798101125509169411847788082169532456136524082502981119587455527"]},"11,5,1":{"tau":[0.2548663956509313,1.2341143111115036],"highPrecisionTau":["0.254866395650931278396619322994925107421875754832445037976762283527149423713457927970943600063778903424447862606181575816","1.23411431111150356160538542324110325359188195138653092169400106874454707816145685442804416716172506977142446602609186568"]},"12,1,4":{"tau":[0.7396692877450731,1.3278970515273893],"highPrecisionTau":["0.739669287745073068109855987724898556840454195485279634760218665045532503730556927184070722752317100530623691369748418806","1.3278970515273893782488678268722960424975586797462486859074873545980987136384116729415939579082667069167869994581469333"]},"12,2,3":{"tau":[0.5895167421845012,1.3451140683131038],"highPrecisionTau":["0.589516742184501228545408769049913229255798883887526080111758942668135474891813272019638083413093094408246778317258041749","1.34511406831310371476129379373759082586737883441460570477911974913230977984227479439979746106385289413443950540277254472"]},"12,0,5":{"tau":[0.9007660134378412,1.2933460768927483],"highPrecisionTau":["0.900766013437841239505043178758561255690888196906815457730856076496529230900388543341013245517475250654006359524585240733","1.29334607689274837511485234351260075451335655406953534985937067305958816270660672474243367848677053178999458007335478571"]},"12,3,2":{"tau":[0.4416188152194287,1.3487462050250516],"highPrecisionTau":["0.441618815219428721056632340238376646678137393439487559650741083345464108663662894579034869575350614639230873113825610311","1.34874620502505169141788841193359130753135470467075258339485869788643010748483255210246204213221694034614920277305830305"]},"12,4,1":{"tau":[0.2902000373357327,1.3392010734596727],"highPrecisionTau":["0.290200037335732702259639086073481241176658440282011834555689055770965110936267061321551696765605850882989811727183635636","1.33920107345967272591868154577484791048381804903064388129651238353724665085166361286999646739797181285141019643482647081"]},"13,0,4":{"tau":[0.8851018568734509,1.4130470995019084],"highPrecisionTau":["0.885101856873450923905257411208284067539688463987610382139162515200853435993746348345011322773283569934040781174096493329","1.41304709950190835566318737602683208917641915068361383487675809433613497508838050553505000360994265350035963428656436468"]},"13,3,1":{"tau":[0.33921821653213363,1.4607754209043649],"highPrecisionTau":["0.339218216532133648989424062426401287845331759739503843955305129448322256702776648276596567777486747077816059765566510054","1.46077542090436485264072943295632674559896902441346901769915198662478325880409976058239865831598417636575639995659956247"]},"13,1,3":{"tau":[0.6966865564002074,1.450210152948468],"highPrecisionTau":["0.696686556400207377748126309236426351261529184821212575837797805336723205479626010024077554949842653341853634077699718909","1.45021015294846798408692753343438367256185236047862669515031446737081370201861908303336301894033455729733237162238945102"]},"14,2,1":{"tau":[0.4116444375676246,1.6053421692915246],"highPrecisionTau":["0.411644437567624572587515200994719079576901797643398815978676807315439223203352587669921447204420468108991260490459994882","1.60534216929152471390830298414042825100378862026946428205417543009035428690284785448574975718305399269736062430935209028"]},"14,0,3":{"tau":[0.8623053124566538,1.5574379665635063],"highPrecisionTau":["0.862305312456653800106013392390580788918850169413618795824432708270565059262736856377337989483539535467953184624886900777","1.55743796656350627931572260017680414310229379059408793091011262536187641446260436011431955415783630702572447631486612904"]},"14,1,2":{"tau":[0.6333267617998382,1.5973221883070887],"highPrecisionTau":["0.63332676179983820209045578775801886127044531199781203221517187378970861988142632440019667127909589309987182791706376486","1.5973221883070886236772229426094613561553440987967423182652353359020752958392085212518151959753624650894628280641145161"]},"13,2,2":{"tau":[0.5184871507751443,1.4648124873052322],"highPrecisionTau":["0.518487150775144308039827802050186685431756320888968036343399419311262870011127440896507454787743966599684636573873052225","1.46481248730523209576995487104361736400612312550593688579969605214083848459059023231857546426347946362907295718647390016"]},"15,1,1":{"tau":[0.5297242175391552,1.782217100746166],"highPrecisionTau":["0.529724217539155260777235136227917979456171101859389499141622348227868706031192532997905712643812076221713035110810517917","1.78221710074616602778900372488444969744609534192805273507796239489805340873025148873868686284169521706393790703394249854"]},"14,3,0":{"tau":[0.17774415781042416,1.5837488636513262],"highPrecisionTau":["0.177744157810424154149854272611675850198755701945741493609958522608073090555211518120912242683969868469741591466725432463","1.58374886365132615498770747208317542232950453146811543398707124362038945610661366692466204581180843636147965657108349435"]},"15,2,0":{"tau":[0.2279882762066647,1.7705613898974435],"highPrecisionTau":["0.227988276206664695422114192091003611046760950521302850294201103796826038535314398985817944631236717554167318744726532432","1.77056138989744357606265732840152101787637736500243660494779013981709489939209719570163649986890129027537500767521839116"]},"16,1,0":{"tau":[0.32381004411121894,2.0263526679806523],"highPrecisionTau":["0.323810044111218956099344541244641370941389497004241901036118770505664332344596858830336643299293310308555767347368335663","2.02635266798065210238653268171159000531614895870986570955391137857441989853097457193821726091414532778014820893610327433"]},"16,0,1":{"tau":[0.7581086739096295,1.9960908376667024],"highPrecisionTau":["0.758108673909629581793835430719561128806960065895027477876181343651418766238021582521438929754284670383735422983144852243","1.99609083766670237507122599887770784888172323661857501766672325186114539735405686864639506309699251477485556840093582765"]},"15,0,2":{"tau":[0.8259755192585663,1.7411427263091217],"highPrecisionTau":["0.825975519258566318359537337347566071017380429151676552842115656473094004163603095059925163403159585477079318848795587574","1.74114272630912179936955852996314200194202701445914065059397835492211746958762447161975557921419126839396595091035169057"]},"0,1,17":{"tau":[0.8581067308248234,0.4158162248618968],"highPrecisionTau":["0.858106730824823339149430431993586420530776013010993974042055849093316938382087849294574992731695991233336044061608579547","0.415816224861896813801259504398630993125587400158801951868975931154482012335453079675624605296921810622102638518584268153"]},"0,3,15":{"tau":[0.7555374951422722,0.47540182595641955],"highPrecisionTau":["0.755537495142272184949656682099086951282855023759596390633654693586574338808034290887993109310251746344633832804024593185","0.47540182595641954104321034520750894830896913805679052560634845364336752671638066368326228836948356897662654195038705942"]},"17,0,0":{"tau":[0.578117895117709,2.4154606742055083],"highPrecisionTau":["0.578117895117708899386516396705718264442865855386177128970844752623157551862844577075822437513746832157570526015637972651","2.41546067420550816831179296743590661699167526073422674552320481276596966651167835844597668335472098478799023803257863807"]},"0,0,18":{"tau":[0.9304042698638463,0.36302110557897177],"highPrecisionTau":["0.930404269863846310058531929880786935261824426972011605402740891176479632113492668369961070029412609274635280516612358337","0.363021105578971777339838154738904955921704599322635659742808076443869905176437607251964352163663771388425228729185796564"]},"0,8,10":{"tau":[0.5483264390544348,0.5323798754250357],"highPrecisionTau":["0.548326439054434772826408495090671940499891667650045629216814877820353646529052946121720720477305898859093439533311420879","0.532379875425035753888222547377023971804096845669288041357404112649619261308818324782450649970455810939241621120732749649"]},"0,2,16":{"tau":[0.8034477051280103,0.45023033091184944],"highPrecisionTau":["0.803447705128010294145842688835889228824409487291730226263607372383465469577588278218923096420046033802176188159664260369","0.450230330911849465364635032928388121076838017324221892975368724080295259759457444507958578031124571695388633348696932771"]},"0,9,9":{"tau":[0.5090029024985244,0.5343759751163205],"highPrecisionTau":["0.50900290249852439431739500184096225893019232523640126755278661701337557713172587623619080567004644646758494271600929548","0.534375975116320545320666514566741399691302923891017546620235719310462084651526850747457293635721022488910031812136969363"]},"0,10,8":{"tau":[0.46968107337626236,0.5336780817848044],"highPrecisionTau":["0.469681073376262338969927902480710762029729937664865482417011488381693610262049761162901687639135749165405638374924269762","0.533678081784804402892708254621620302173035734259049638376098197085368885515444124091173550268602714977954215528931893308"]},"0,11,7":{"tau":[0.43012385509176465,0.5302530346555953],"highPrecisionTau":["0.430123855091764636828900968464552005097363472362161478272518105971564159272960673795235796885990101569659608409760393558","0.53025303465559520181792229495893810474698132344346567076261278036297572855575770757783986282085689684695650206701733604"]},"0,12,6":{"tau":[0.39006631592197105,0.5239693742808247],"highPrecisionTau":["0.390066315921971044105493558568828075730040234701496159827209172822189920866732586080269132749963677775890575555125323486","0.523969374280824728022685153355247585713622441799036129917108918063204753998160766972227255477353211351051694895263075209"]},"0,16,2":{"tau":[0.21521813746014437,0.46159171463273746],"highPrecisionTau":["0.215218137460144363516980928441132230594344490620124861156118701470151996213518705638063710215825273580803789681838589465","0.461591714632737438236225220065443675802440389710114880381838044626386930401442028624188925930464644267834285114165469506"]},"1,1,16":{"tau":[0.868240384538875,0.49458049152750533],"highPrecisionTau":["0.868240384538874957311621171926812301485042206122385457393762543426769031136137584715790371627114323388412972539317176049","0.494580491527505338016984623101736320961873296121820322633798051262322946289880436816527431648294822095008237597689596138"]},"0,17,1":{"tau":[0.16137071924505175,0.4305216203882836],"highPrecisionTau":["0.16137071924505174673089905957071366716518028799416168727644484905523914233504693200683757416820552700224413860482945194","0.430521620388283631316678807956915556621800780715964157511001431345865254773918923801405441692837324627574387755549869947"]},"0,15,3":{"tau":[0.26277835830337243,0.4844044791623395],"highPrecisionTau":["0.262778358303372419573660154719972863146343963567036170385240667967370197681090850885764320339421033836666804217591999226","0.484404479162339447994286334831943586080338848835706420162884044804600500795616367078282297190775874909486119597712958227"]},"0,18,0":{"tau":[0.09152935378749712,0.38408442019550465],"highPrecisionTau":["0.0915293537874971193555730385062372331235216687913567208722581857883008954483733828639029778150502714390635669350355668888","0.384084420195504676339471943035267743115062168471573180738625264464951009961020870942662156756711688332243549373850919431"]},"1,0,17":{"tau":[0.9421093004878299,0.45623662960321926],"highPrecisionTau":["0.942109300487829929373347010523683451508295177335761638443192389622750306035874431020868153110169297318228449439824343612","0.45623662960321923376963365711613435401421988513173437994051572259407656100867504367367153082478859365067767681927979123"]},"1,8,9":{"tau":[0.5297120930707151,0.5911666899928312],"highPrecisionTau":["0.529712093070715119430338876570871129956722068589198822191095488659526350185223135027543428506171543274906157313296709485","0.591166689992831198776959025320998627009370159438527170198899293008261199018237224971524509412719917640012863284881316223"]},"1,2,15":{"tau":[0.8093973506471155,0.5223998975948725],"highPrecisionTau":["0.809397350647115514260057675686036413520617671847792017471949543187500573818122022630939815113736864581956649948057843229","0.522399897594872463905656091768674835521638913435249437650088029211932324312532906935536254695984582565682690445437006693"]},"1,9,8":{"tau":[0.48628394890531523,0.5917479223720434],"highPrecisionTau":["0.486283948905315234385922274567745899169243089818310217317917447347303634624537333919513623650919735572375247181146669651","0.591747922372043467374916668611380423264332960225199478482737194476021061683162449127639107672737250676218757342260477383"]},"1,10,7":{"tau":[0.4426957169285713,0.5896126847619434],"highPrecisionTau":["0.442695716928571309057831420223271441526117922445362039620652805879317226412261203225725058320892184144426104830523433154","0.589612684761943383717520343302001232866216571126408143184814012048559084606923606323762963167245646811634469068164590187"]},"1,16,1":{"tau":[0.14788220657895748,0.5057385588594627],"highPrecisionTau":["0.147882206578957494732675297031491027208921947733628783197181529319645545967945394801510610958686701780330384013551346841","0.505738558859462781710489772590869234600882041345957568293596758526779340044914455526541002515688865891693896568771964589"]},"1,15,2":{"tau":[0.2065482866019055,0.5312436872220039],"highPrecisionTau":["0.206548286601905498866084803523627589263859096859739121423112671469734661611934425959245950700411507733987143338416455309","0.531243687222003848031048046040592219213310604246710230212771829118859522878896773215249696549829520438822447059554595319"]},"1,17,0":{"tau":[0.07449367671782937,0.4713322500357738],"highPrecisionTau":["0.0744936767178293757895230370902962602903434752451429084721966991149896722603302068906238859936087354366649726204330459376","0.471332250035773756176596676424208477539121849123647465537053811746894817352420351240391459669420911785443530508293417471"]},"2,0,16":{"tau":[0.9458092115020454,0.5296760801204656],"highPrecisionTau":["0.945809211502045479732400504659695405066240940101068255327968084376519063392704904709068597929122544093807835621108414475","0.529676080120465652105225899442424200457572380771351486174423921841292486007106288625943449325985379018918085944194696433"]},"1,11,6":{"tau":[0.39862041063761217,0.5846858046441265],"highPrecisionTau":["0.398620410637612154240775892506684963104277550684774276945200799761570692301213824882977487840534735512618057552233510744","0.584685804644126444241861373260019656383992970980646897078457544831126380446694265299516194279869107578516394993721047089"]},"1,12,5":{"tau":[0.3536656753052071,0.5767950027733381],"highPrecisionTau":["0.353665675305207118474531237757259968540682240864919511815335282810258308887458788828183149314854369920301527884606755104","0.576795002773338105226338511993633933038302171962701639243216388127047693857933610158403070459465145029404127084926745663"]},"1,13,4":{"tau":[0.3073021365880327,0.5656399477881907],"highPrecisionTau":["0.307302136588032687377279910078586282430615962139443356044392229488196248981212517994426221266312814813345303624887398143","0.565639947788190724008980658553261935350734255629672959672083835742789364920439416937437772691089886063150919146493340139"]},"1,14,3":{"tau":[0.25872667253800036,0.5507285031137847],"highPrecisionTau":["0.258726672538000331524601932109405679396960882050337300096413759582010877189480974502119079138929057181801728769884436673","0.550728503113784740485844347914802302918272990765108280077177698210352169595033492595572286485094092789155668585385605017"]},"2,1,15":{"tau":[0.8709719205299945,0.5620583347704555],"highPrecisionTau":["0.870971920529994500435334886986260661316614755387249372903347391246032257050767863022929772505930613333072912087360757132","0.562058334770455537041882787213685460073716753446088745351728816315989761978859048192372464252381509581025986652947262766"]},"2,9,7":{"tau":[0.46025996191187296,0.6445616377422754],"highPrecisionTau":["0.460259961911872963500518374369067387330956772071087733844499546626603817786843630202193429564181645845056804264960008001","0.644561637742275427409213272254898018282243474019923139185573949585427214779464484930019029524823951976209650508447836195"]},"2,7,9":{"tau":[0.5551368930399545,0.6433935979559583],"highPrecisionTau":["0.555136893039954462329174628004120006185629543954243361576594883759049116825026566226023931780216062199893455590730115215","0.643393597955958270782079280591922577909362772530709011242499561543799669133992643194534550311799848113098540118073533404"]},"2,8,8":{"tau":[0.5077049539338357,0.6453660614372363],"highPrecisionTau":["0.507704953933835759558840659721892396856648111573677519959602141243998656938783507697482817545412992730633263327132536723","0.645366061437236332820980514165021830753114635985769705053399695456209100017318166740660060877510922830151333703308008559"]},"2,10,6":{"tau":[0.41242380038981874,0.6409498847956074],"highPrecisionTau":["0.412423800389818751471284494902233173046967247253548705020114660217893817052620221378973440030263658079998161080542561185","0.640949884795607394822930534838807320297041026031456331854426700199187578399834831339183257811515515865967549523970419026"]},"2,13,3":{"tau":[0.2615183413745217,0.6115750672468322],"highPrecisionTau":["0.261518341374521652156690627653392623924775891476526388906343476138917480567550870796683164134777429991822873403305513677","0.61157506724683224082368246399302371627055266035332154204854913722814271354267629130658523139468346055318812888145122725"]},"2,12,4":{"tau":[0.3137333901986152,0.6247434061398599],"highPrecisionTau":["0.313733390198615199005089091666901347495962868365048699521518202046609024944509376351338867127541261202913311425480637579","0.624743406139859896593392468531007587482180752367490452394766409477441262952708968249850227947196534854945436344200264207"]},"2,11,5":{"tau":[0.36376635523194567,0.6344160714796924],"highPrecisionTau":["0.36376635523194568896111956282205506844402058065713263357124907667102127192489860790364314765195006412618278482545902941","0.634416071479692421645640984499896361569486624142437097147059061954120443188899855078059341924663584877267578039198670689"]},"3,4,11":{"tau":[0.69100324140063,0.677485636568661],"highPrecisionTau":["0.691003241400630007737956064648039423348431569665928714167854834515051225601955033236296084747669536627447070429434930865","0.677485636568661014082460448823573922287873679580713880924898794698223206143745029587249506438994364628244495815949875679"]},"2,15,1":{"tau":[0.143893164978368,0.572094291991783],"highPrecisionTau":["0.143893164978367986856815708574590352098987677012388840052161328096588345875884426042501890872839492993351425673245688184","0.572094291991782927145842173386189054161422187227423961321394313894783198954133624633508304277410297919543553965213310568"]},"2,16,0":{"tau":[0.06887133935499416,0.5431928937321099],"highPrecisionTau":["0.0688713393549941581384503799194094917998749012232387089842908641881785969658909996916300618746030951471033261256728524056","0.543192893732109852699709960886676530183816326069469003809338251249275630195748628185455359622717220031936068343181262685"]},"2,14,2":{"tau":[0.20578290474375857,0.5943390725221993],"highPrecisionTau":["0.205782904743758582906661102899575162794079640275579532992885324545725328573261113611933150352011950170290072274880204602","0.594339072522199275153017844666709000143081138179454928212808537182552471618339759757725316052322667491487661224242708852"]},"3,9,6":{"tau":[0.4302069132781998,0.6955588511334477],"highPrecisionTau":["0.430206913278199799282591958617788733577089430695700480357875206913481000520167497316621754836605187295301000372728945685","0.695558851133447731735405370479572003690353213548695844356828507496173347926532138525026031072657938979582758180078996675"]},"3,8,7":{"tau":[0.4819222282721553,0.6978406377890127],"highPrecisionTau":["0.481922228272155290682441071080766601749806257987095378610923805774178335854671605930031216888003712789941988434465090384","0.697840637789012697804813978865467924616217981843620126246356059227352074638138918789864768494077526066701706969138162043"]},"3,7,8":{"tau":[0.5333979609264652,0.6972271921778301],"highPrecisionTau":["0.53339796092646520386179410053706154011395363445268499657753223979718961720877988092784369593652875948552931628100480722","0.697227192177830164531068858811498473106412793151088687452217513019370166509423761444043440985503259774973428673165851425"]},"3,12,3":{"tau":[0.26826784385286256,0.6701642764133235],"highPrecisionTau":["0.26826784385286256389559960129299667460628070047987189458038422767908078233901099205328335766277261240743003580989899973","0.670164276413323517872171473914578492129990391760275678060748609491154520905262457189251061541671980598261639513514153353"]},"3,10,5":{"tau":[0.3777818570031734,0.6903116816371236],"highPrecisionTau":["0.377781857003173423724557481367548615826665988777552562065419889072595133394409034620588859280221059677651316193957289415","0.690311681637123628898600676994080992276556832507218927888368661902517327966237317355759357498346186246705357804466763849"]},"3,11,4":{"tau":[0.32407057625497826,0.6819387314017876],"highPrecisionTau":["0.324070576254978273279044692316836355728054939028697124689924394204903783749060146616041927612563054013928986917258705793","0.681938731401787646119758469665449777350174661535740683035320977332026790107467851944130356369562198070139506018193627821"]},"3,14,1":{"tau":[0.14407311078175086,0.6343782027076906],"highPrecisionTau":["0.144073110781750868627072654003993562559886688676064999836492158149805208935299202159179489749950210659349770769725739092","0.63437820270769064877477633262214917901241080055768942230578740164616203642485282927477939792820038254117673747367734891"]},"3,13,2":{"tau":[0.20908481681707308,0.6545468733305383],"highPrecisionTau":["0.209084816817073070708661592492710685263907776659768933384283792750501998754181631650703133748117775473900583477377753902","0.654546873330538343308534008733930030849939693404833354478370835693840085010775681432979954366728327649679104606129849905"]},"4,5,9":{"tau":[0.6200490006829161,0.7431701673875538],"highPrecisionTau":["0.620049000682916161085096581281224843363926609149563591206502329195088902964987464891139523241395148499176223069599691378","0.743170167387553763324173174473805905684763033484723305823133053403063721620099037423577308370275641851777261627804773392"]},"3,15,0":{"tau":[0.06698479209669943,0.6084538741197194],"highPrecisionTau":["0.0669847920966994325302858940754053558986895624498034463230790102274010751844066865203206908597445030011011385191562545252","0.608453874119719421811820443413731705185906460230876869788557865041810984013626187476890152338575372920111591290446434364"]},"4,6,8":{"tau":[0.5636604538658759,0.7484758944497112],"highPrecisionTau":["0.563660453865875913417226916448670301076055442809869571967787136894955497220454395353062467722267866231119803491667453616","0.748475894449711133804418913750594585919199123978293889084808624904250039138115805794053260522687091656680373935765480843"]},"4,7,7":{"tau":[0.5077673938350309,0.7506708467722549],"highPrecisionTau":["0.507767393835030891594238124767405281813889895719106701772189714080144041133701695747545574445571405141830091750468860266","0.750670846772254979488594783600718221595091466272120209876667651747565628750047179276591959965725970780088574648241783301"]},"4,8,6":{"tau":[0.4518459370109677,0.7498096322903067],"highPrecisionTau":["0.451845937010967707912337645823096366545090451820642949250319091510390179890831511287922137294567764815960819223448075036","0.749809632290306699726974402477301871691661081749857811357280821794350514663609201839199540846058963120158865189070813854"]},"4,9,5":{"tau":[0.39537122542855485,0.7458626170308591],"highPrecisionTau":["0.395371225428554859300957655607016102781424653620929172503431860880001760516615645128449611317170900942309355835673398914","0.745862617030859052121604270655825966755066030436355413031761535726550426763749712128282050785980679538873896393841306661"]},"4,10,4":{"tau":[0.3377343362059571,0.7387123638545804],"highPrecisionTau":["0.33773433620595709506523378747721395127021038508254653997663626401230060949196430755874956798900326037358826493786999013","0.738712363854580337113085151334424150564975517817321373851130320001693297537198762733802477369042625214901034988565544426"]},"4,12,2":{"tau":[0.21526694802382632,0.7137656103430466],"highPrecisionTau":["0.215266948023826313144036181942360157219519739914108445288081441693689487616606957299621531923448066370170951386548054356","0.713765610343046594487114766202606911481382645957564305477315256139113845691791633951143823308819520817244249890956873565"]},"4,11,3":{"tau":[0.27812033941803743,0.7281359308009079],"highPrecisionTau":["0.278120339418037406478396223638800134484892869497061890871956144840254477371870129588202094633106610993852204983492323857","0.728135930800907910402979021653101623428636860717948547625314075240972656903800244214377507759108013163790909633311829837"]},"4,14,0":{"tau":[0.06699129243406791,0.670818046356182],"highPrecisionTau":["0.0669912924340679116323677592406914558028273602209530467366549290752451563118902479682377800762907685446806593274841091196","0.670818046356182033031693173040186763734692007042844314520344987103977591715983323193031655758434166708485941679050772633"]},"5,6,7":{"tau":[0.5382994716337128,0.8037851953730123],"highPrecisionTau":["0.538299471633712881526925530105867704149874259585391847666442997701439977553389817645464936434814260488376347793842610017","0.803785195373012308868368168305488211049037859840118321353983769643419605093461293797079894004217656405098041641146954546"]},"5,5,8":{"tau":[0.5992022282503476,0.7997445996205952],"highPrecisionTau":["0.599202228250347662358487438907221659018177983055716378427117283615601646454044706343564065568988109441856828958438072633","0.799744599620595200842602617388093319947679333376814681431805843793898456906671840125023463174661452504611963182179976401"]},"5,8,5":{"tau":[0.41668214160505107,0.8019742286190437],"highPrecisionTau":["0.416682141605051073792629247455492657216615303549515850641260395698312484469530430414743849524308641224914668728206221594","0.801974228619043679485160784328495128584485603758345543493743656555127986394575681974315390859928402057199590228353920306"]},"5,9,4":{"tau":[0.35469954263238657,0.7960592334705477],"highPrecisionTau":["0.354699542632386542126121605527331307823882703085971399658456367955956583353901161327409405194960995239066361573875018541","0.796059233470547714518863307725404674410406295588637932506762310415971968514146123229474803653968965349927762289310226005"]},"5,10,3":{"tau":[0.29087170608200474,0.7865940293978031],"highPrecisionTau":["0.290871706082004736075801148585848751718450461181982815206042056615042626109901431315753138657386962339644612060319516549","0.78659402939780314536504667490085563832937564108073355074003898422174804670693875610927999506027316825001924044579821854"]},"5,7,6":{"tau":[0.4776619283645782,0.8045233573446753],"highPrecisionTau":["0.477661928364578224149354367465936881943025720243907293216472499650995211280227805768358708366817566984539871069603131665","0.804523357344675303145986614888097953873725849135634339258903294150447711139114947086191746056707477474801342251237273682"]},"5,11,2":{"tau":[0.2239368007164287,0.7732545182471444],"highPrecisionTau":["0.223936800716428723652911478529672600330216633460087067865785516397955423453338675660443220926000964887953984372585673109","0.773254518247144416679581227089161076813859038473605341095992717649421414752620999072434944271863770002098279206305539331"]},"5,12,1":{"tau":[0.15161097690762873,0.755497750451536],"highPrecisionTau":["0.151610976907628734347967132670147925332355236074594746699433851460537592740839465412509400394584547572599019898578568336","0.755497750451536022982905851441086958069097062824237820770448054319338481820244387509138094161827763141567290323938396924"]},"6,5,7":{"tau":[0.5743815100018611,0.8576793453455837],"highPrecisionTau":["0.574381510001861126762363145523492357317214030923500967547231071091435717730198022163821068987381363981169220412565225965","0.857679345345583621260667367476041016469967804479270688565050669339545298646840551313769272336991272103391164735829010043"]},"6,7,5":{"tau":[0.44218381634721937,0.8593537507346175],"highPrecisionTau":["0.442183816347219361066714614754457280247985501312682496047816500229429510688865220920854695694708863795661217249154980152","0.85935375073461747294503479445540384436534665243824048447388362134537211589655323699450019961046833276095666914239783397"]},"6,3,9":{"tau":[0.7094813360206571,0.8412265709238264],"highPrecisionTau":["0.709481336020657118986407507051937461756427952562559312064834606719931268489430171513003348194204403438118521362135518871","0.841226570923826357421309147501813070164622077161156020757163977150283813180896754583491083442133779190411736124045416301"]},"6,9,3":{"tau":[0.3066570239860549,0.8464270150375546],"highPrecisionTau":["0.306657023986054872576303213467324602822364974896889154829479904184705780252807423875910663996040821887601296793542835578","0.846427015037554534210931637911199202586462745747588234563055844122372654152708429451264111170372997743842013859827926703"]},"6,4,8":{"tau":[0.6411446705013957,0.8513866492555786],"highPrecisionTau":["0.64114467050139567229102313269057409974584906533293099816032032622617690030035266214907095653680458248171270467771993471","0.851386649255578626173056058473073724940153241788755355301953151375168470650468605006695405941856692970402486053559152619"]},"6,6,6":{"tau":[0.5083080815462201,0.8603098136314558],"highPrecisionTau":["0.508308081546220049277680297188614411925807297715296310375716676216393949211481828844878448885546189344262001341468252275","0.860309813631455788032081617856694246154786960803753235807035701214193941447287129071063807569128793671934614943199236921"]},"6,10,2":{"tau":[0.2350703744154619,0.8340096626203388],"highPrecisionTau":["0.235070374415461920817977084119720362422826196182165204365859591711170523252195331480621955209720688551916880742106579083","0.834009662620338793720843383453708954469821204262595145452995324971306909702492876770702546604084131889355357819786042339"]},"6,8,4":{"tau":[0.37526500160904913,0.8547762534282708],"highPrecisionTau":["0.375265001609049130150741132890518663166770290088397934860520981675904098775267386005753029262964101624260481689741108769","0.854776253428270762615378242686666338542010440833754045490915999965954695547661397988841918628799445153606343539076742128"]},"6,11,1":{"tau":[0.1582226693957346,0.8170001114891607],"highPrecisionTau":["0.1582226693957345953277957697631114692034435321157904628944108398711062337841591179954805242783272324741683404577534042","0.817000111489160708813448188267628471396585365161574956685978474408204155540083381018058888310440552799499172122167638172"]},"7,3,8":{"tau":[0.6911721556324533,0.9035019433911835],"highPrecisionTau":["0.691172155632453303033352376545826967974478338912011955352904090080585520508121518491769803772368930928927542340286225424","0.903501943391183518129718923956712270255356756022402955735242127918973362384278283323225424168723095133668409619612444427"]},"7,4,7":{"tau":[0.6173142719831752,0.9126750539093664],"highPrecisionTau":["0.617314271983175209652056644565385892767601183864121726784695548506981553155041654621716863068618932847920541907964410673","0.912675053909366399125169919578420769474573716264735425747558361766685414548767287337835483180199685983329007110782574062"]},"7,5,6":{"tau":[0.5448014504978209,0.91765507474021],"highPrecisionTau":["0.544801450497820879660189234963991567274158039583541669604322564637760937896175368289428693721612185356164931488080059004","0.917655074740209998858546869658669310893832196005357237626045839023135512032684804941888966896396428973206966068904485968"]},"7,9,2":{"tau":[0.2488877251769376,0.8969366029273768],"highPrecisionTau":["0.248887725176937602105594350706685877503658677134215012180950947618872532048438483383086476118705426928934277653480307301","0.89693660292737681082898078409235783500614978745500446478035525677211115284366601407186447164625577854242305206216146697"]},"7,6,5":{"tau":[0.4726612674566658,0.9186163137101919],"highPrecisionTau":["0.472661267456665839656339837009690734221884466539212207855259378824526579337145337014311538426479580099800267118744739157","0.918616313710191940951900702185459600053691879658309625277597607650762050295546839963593564735709902668045702151535802416"]},"7,7,4":{"tau":[0.4000123461053018,0.9155851409324114],"highPrecisionTau":["0.400012346105301760676739577222386699750883973672489686456481705089256399496658012841602516758440437690070454167100387073","0.915585140932411394272510273513674848564449907324524052402076067621279230029847751987910434140384304751904532780115620451"]},"7,8,3":{"tau":[0.32587133209902597,0.9084501954598291],"highPrecisionTau":["0.325871332099025987273840615306371665495057286713865421912070372665705329669555732424452961568047397031247656367914246169","0.908450195459829131679929964962732451268612567801726484146049472849261359329221597202685296859260375489269659371595065552"]},"7,11,0":{"tau":[0.07375931442415441,0.8581151398477355],"highPrecisionTau":["0.0737593144241544177262178770664170026733315615270455659638997652194596360646033261867629082639718866870306356443199073729","0.858115139847735509383540006181085127271551176620950022440658723883616883500517862016283133886290824267763947118518640501"]},"7,10,1":{"tau":[0.16673296213132485,0.8805189885972936],"highPrecisionTau":["0.166732962131324861512428949447128655590706366185572067829884313070240402173858517741771100684818527121164298880898616165","0.88051898859729365276606212729770752120005892769879530718443730886290537928270210758134850859335987146174555804120200279"]},"8,1,9":{"tau":[0.8395385129403188,0.9369232712115018],"highPrecisionTau":["0.839538512940318780743614341878048696977828746742066149730729685637013696008889324686608543511324934047640668633870423863","0.936923271211501743360308062807789500810129288776228673430292497619079088295809402665039508263464026621491880185940780852"]},"8,2,8":{"tau":[0.7518880694579414,0.95581495977596],"highPrecisionTau":["0.751888069457941361687282732968786676939391984950116635858693693451099484901850859501473420433046373369311412338565319244","0.955814959775960021751649907499108920734220367329820157585939460407838520473224561677238041224145551512885213505676337772"]},"8,4,6":{"tau":[0.5886504723021578,0.9769324012604373],"highPrecisionTau":["0.588650472302157784751623128435047269522223269588179053903241907423838824674463794289993744480534138003444246329321944223","0.976932401260437331909265810133392344887647720509356168775384213019124226367771792497258101069460887554557040403317475816"]},"8,3,7":{"tau":[0.6690395080213145,0.9688843573768685],"highPrecisionTau":["0.669039508021314517747666775424087967610028058772872750486888006551495580476506527576941732927153930655925677960158639695","0.968884357376868537013971762181993247797251427749984041364484202135894667191716079331103761372909224514536967428865461652"]},"8,5,5":{"tau":[0.5092953828784125,0.9803259950160266],"highPrecisionTau":["0.509295382878412460160874863445273204760199384593009196025583766694586597578735796240569985642693638011028567629553788161","0.980325995016026601660731031087797115692684950146702645703710514995476981599869872297231648799479603358252389797770592254"]},"8,8,2":{"tau":[0.2658341795582519,0.9629542134010324],"highPrecisionTau":["0.265834179558251906184780744174517261130267776356611378241974078003915893808640377880449734960171920395706005005465388361","0.962954213401032350014678224491675817704879662583333195362263409608447729262293915716782125319469275331860489339784363553"]},"8,7,3":{"tau":[0.3491819095207725,0.9734884893913437],"highPrecisionTau":["0.349181909520772483894061450706508699300353583440222469844338988527640787407528160342626324011300628972139201238785568017","0.97348848939134373752225577895605069027337809096182638798692143980629236120755061466833871007749561366099308927735323257"]},"8,6,4":{"tau":[0.42985007814956516,0.9791957998884133],"highPrecisionTau":["0.429850078149565159887574113112832127158646765704092472651404472656711561830028059111965217408106875580807437081244559633","0.979195799888413362268704454977788272782762426318883579806757039695813096177045184542422446956827909421113245084743585855"]},"8,10,0":{"tau":[0.07800760573923333,0.924608569384303],"highPrecisionTau":["0.0780076057392333347810847420146329237860400471300890482341723628990861438108922335062202869252496575576317613575713337923","0.924608569384303038115382179588589345857242749095849790649701849374701200563299567961659868012643647322836980311590419243"]},"8,9,1":{"tau":[0.1773807650596983,0.9470547928983137],"highPrecisionTau":["0.177380765059698310026326643557391732103029335355937767611771593845189500073879485106684580714346598426445996407189164144","0.947054792898313737863270349241747932849273266431847645708582430327868674753038694896065171929972821251225665795458717829"]},"9,0,9":{"tau":[0.9331842649237534,0.979590332525529],"highPrecisionTau":["0.933184264923753397270023914521447535752538346527196027507960641661144351139843905925895943937009339504912978031436645253","0.979590332525528920769751495071710077299300012739252808582589753715075639267576659212186042147503728484517221901069432884"]},"1,3,13":{"tau":[0.7479035941194281,0.551862581693022],"highPrecisionTau":["0.747903594119428108850661714725156720018273118395432890300085816859486496201053217506247880443505360640476612437183956975","0.551862581693022077945491115891748133239500728180073637215342974188936155884351250574732645146835030043688769472603577599"]},"1,4,11":{"tau":[0.6840478176449494,0.5766703831417671],"highPrecisionTau":["0.684047817644949396013768158048674175430751073849205002021067886163558304167677907974702968841722641001785023110251402386","0.57667038314176710027562738339812420864289175899966778118549241769311071316966004579529676488556917040526056017813824328"]},"9,1,8":{"tau":[0.8276922492207691,1.0072536482356895],"highPrecisionTau":["0.827692249220769156905251916710113744177064486793344828945288447120821955105289232821070259380408334551063634713499487368","1.0072536482356895211424868966475471570656096704618529017904035389297740963926889143331883279965302000838518812929025467"]},"0,4,13":{"tau":[0.7011221708287709,0.500080946311254],"highPrecisionTau":["0.701122170828770941901751152954981188201730874688290684659879376923033485555859094830543200539320400799054311535443707255","0.500080946311253966006848521934807908488189686897987104381766537933119699412051944374575116336763905970448485012075239214"]},"0,6,11":{"tau":[0.6139756203817293,0.5246093528992708],"highPrecisionTau":["0.613975620381729206388320152241657577414847786937602623444567907540834980630217372419984797691448542882084468230111465607","0.524609352899270813801358180959834359590212589775249624246450110331602869018096262414135147328662638234961104745547951402"]},"1,5,11":{"tau":[0.6483474411193695,0.5801750153093087],"highPrecisionTau":["0.648347441119369506379869003995699749653315721867577883345147444095189724316913933191371490435882871803183890988098037245","0.580175015309308739031419662219800673850722245998096019814634562091964916029933746257058945020149399018698856826491261381"]},"1,4,12":{"tau":[0.6970190646953709,0.5680903834720774],"highPrecisionTau":["0.697019064695370944184250406822075302003990209066398715723979421089725560947963938473331552025059349324710770981925373455","0.5680903834720773566087209560893350553628794402825262432519727072779048569497756966403644453516373510342541315800304567"]},"5,0,11":{"tau":[0.9391548840197529,0.7581483362783746],"highPrecisionTau":["0.939154884019752849097480861014118829687738107688040490262631246045423125720088644660683240232987396964272857796007718169","0.758148336278374684149620841023681049231825949116191999001783588681668182625866054403805408913829272137469925780502728603"]},"2,4,11":{"tau":[0.6881596223022413,0.6303309114081083],"highPrecisionTau":["0.688159622302241233747665238777100699745979988490833189054782254757950640671792176671332258794174039551805108137563944827","0.630330911408108264719494393217619534090738820769076268918046163519816495784272439548583723828084122373114589509896367055"]},"3,1,13":{"tau":[0.8645460855300499,0.6378181426486293],"highPrecisionTau":["0.864546085530049844886949623636243680951272250989744618371971214025478451318669122347914545420661103565957012773084621223","0.637818142648629369948384829398878457329017035263575859522464886543254559461684772784940344196814813213416259462221999971"]},"3,0,14":{"tau":[0.9446418281366732,0.6076884143324295],"highPrecisionTau":["0.944641828136673203170350631149874087652627307573611203517140593546707079459656741183110278712823006330550235614989479992","0.6076884143324295276694956444936794584398266553332672691381067720409455566327605914476344764222122172405900428383504721"]},"3,3,11":{"tau":[0.7344648993265652,0.6776525723096535],"highPrecisionTau":["0.734464899326565231114165578813678109879274249160871046971025385063539047770972034946452096822262626815861650881231216895","0.677652572309653441903787917144528516821649382218181767179143592896329760009123267363152073451772062307040314236960094948"]},"1,6,10":{"tau":[0.6010396515382281,0.5886710365636824],"highPrecisionTau":["0.601039651538228130503785057843973453738508169972569310781049971040824051617642002425071028825556932746133930384541460645","0.588671036563682394812823129512598411083489462126455825278803298850323626481504489256574777956446893549581436490677729429"]},"0,5,12":{"tau":[0.6568943193216316,0.5143029870426703],"highPrecisionTau":["0.656894319321631550091775802800276377867305374950976297409751387770585981680279491902365408051391678928298422627280464132","0.51430298704267035317015895264570870739726638458539850216500369448726527413038727445317557967709881934219521678881872752"]},"3,2,12":{"tau":[0.7965858678939064,0.6603964560531203],"highPrecisionTau":["0.7965858678939063721991616091874014647041523642964546225929320454970691696060284384989323694743830725823063655755915998","0.660396456053120281763178055587228925741773571586881776436258406445621138335026766802297462944189317873839266320988103239"]},"3,4,10":{"tau":[0.6756739139589181,0.6905599102137322],"highPrecisionTau":["0.675673913958918080926625927423485386578605684930566821805444857648337700983976304573396122080081209257873848120227177709","0.690559910213732172694423386829092271081962038947267537052729344797452690850618218825244180431250400794954512591098838264"]},"4,0,13":{"tau":[0.9440696532223058,0.6730462375664416],"highPrecisionTau":["0.944069653222305811560240057687468760228316111781805852265309720915334524849026789612834021107318043605626759217493997809","0.673046237566441564618427103733688804546490050885906152864783696055003810298764698482583598615335182742877621376327439791"]},"4,1,12":{"tau":[0.8609043812141908,0.7015442379952755],"highPrecisionTau":["0.860904381214190792963848506840507016546687347385424773320601397427901637519792981964816416458061107684621069575726656196","0.701544237995275584314572802994909015472036780909064807976002108838570250628138408063503033307859281744611956646673574408"]},"12,5,0":{"tau":[0.12581844204835965,1.3144942693325665],"highPrecisionTau":["0.125818442048359664634429715377800191739763861218591552427820232281422589780716736199602348061201937308423274156191263013","1.31449426933256641203029995109289894539585723912292509704109819774875574501280762347593310501334633022936830344072644173"]},"4,2,11":{"tau":[0.7891349396393562,0.722780742680429],"highPrecisionTau":["0.789134939639356203880516793694240509679976404906548972017114650057821050707769787701239697654933902025685168146283820431","0.722780742680428990574720910554162901927741948218533262174350938075746379405582381469055969221207398168805191053191152872"]},"5,2,10":{"tau":[0.7791871857384436,0.7856912631515702],"highPrecisionTau":["0.779187185738443583158693672623777795099072805265872446415103505890135971250679378855807636586675290459616162929226228893","0.78569126315157027427427247200976500915521223609014373148261958727088745075676351266261988114141941818158651862776365007"]},"5,1,11":{"tau":[0.8553754602027612,0.7653287324310404],"highPrecisionTau":["0.855375460202761253344429791942398278284067760379750533888602966227559397540071091668308857716197997213755253082225152814","0.765328732431040383452401910002233029616490566891686103970395744724279660004231143034981709475568785150208582996453918681"]},"7,0,10":{"tau":[0.9369020952630526,0.870547492406511],"highPrecisionTau":["0.936902095263052569782964677868155522074651104465082557781366372860865381773895843498414491876436643583935395082090892634","0.870547492406510933147513806434361422555724511760996206208257618577248695312195352153814917226690439663424779696029888273"]},"13,4,0":{"tau":[0.14678003433616835,1.4365206725843267],"highPrecisionTau":["0.146780034336168348826278129842665971070089417089524623547419399199762707770746750906367846282068660657679590548940953754","1.43652067258432664717730573849309890357966047526473457471510912454450070231707200802864110334103566366236830340031996566"]},"0,13,5":{"tau":[0.34917668836827914,0.5145723659523542],"highPrecisionTau":["0.349176688368279147751996901888507108212088638175775920344097694381443801641963299393469302731332473516196280081055406319","0.514572365952354201685532410928291044469207597380441946468466156705280729446499034676165527126574116465520981229401160202"]},"0,7,11":{"tau":[0.587890364909844,0.5276307061839873],"highPrecisionTau":["0.587890364909844023595021593826797226607060864280067751412246719922544339254305630048263870898514918372406305525674574154","0.527630706183987305530701170924692872643204697647985085203474419238330510278114473516550194030206751637508144004508120132"]},"1,3,14":{"tau":[0.7571729124183703,0.5436567942976525],"highPrecisionTau":["0.757172912418370247258234004067836869054186844366209399323344564177939353200092413251216313518430457455495982483711407184","0.543656794297652496906927366014203333430750207923011251806811573388258471385733144656012943251017251168714060649027092743"]},"2,2,14":{"tau":[0.8092023704182623,0.5863915951007022],"highPrecisionTau":["0.809202370418262229340875002652157608953483957706209798967816784892949465175007416542812003884108748645281045804332844504","0.586391595100702171909958886863063846119476281778934482670868343128378003272079858588270560802485204027765655793168935924"]},"0,14,4":{"tau":[0.3069907058592255,0.5016281817944742],"highPrecisionTau":["0.306990705859225479736206630330010480438046560697982251320702215641092995201138946802995860638202854838175529874515305979","0.501628181794474198568706347546915244526415251753044714170152705909943333894008783738326967242375417316289851175629248932"]},"1,7,10":{"tau":[0.5732881589625213,0.5878593365552024],"highPrecisionTau":["0.573288158962521277952599488978601753818943421077522792054590904569083132129241962157207465184300874693718165215154935343","0.58785933655520236253463067366676818502864821555528577566547822267730166365737120429300196992915819320016008979003338613"]},"2,3,13":{"tau":[0.7535840911838629,0.6052671030714027],"highPrecisionTau":["0.753584091183862915571421620227331129494362961911215208359319120279913931690129929524781701610856373845787999279667910775","0.605267103071402758532609517255341361079777427788647414046751123871981437326354518328483406309670680314994848590457735939"]},"2,5,11":{"tau":[0.6515279180621143,0.6308322457280646],"highPrecisionTau":["0.651527918062114265212068838053819095187329325979674071213375773022148222130643430447758635518288255205093731157837542404","0.630832245728064616140317792702786842505314832937584357555278601052723482488392751889392171890850533894016068854585750638"]},"3,0,15":{"tau":[0.9468867730401974,0.5954711621819205],"highPrecisionTau":["0.946886773040197462630576897741594163620307975976787366373652578553248730997155245099828088834877853345769032984459977515","0.59547116218192050667071993808160598422255405502718573432471900990568011563588165389699439611093605193424188162183816373"]},"2,6,10":{"tau":[0.6029341837248035,0.6385947236910482],"highPrecisionTau":["0.602934183724803399766402982069953318607325044854733907580597915982864942135762466648876313664463258211265338545750606822","0.638594723691048224135434311179696823553989652325733808716571476238030751486300608872315387230552301061411559498954198728"]},"3,6,9":{"tau":[0.5850727549830451,0.6937082583491055],"highPrecisionTau":["0.585072754983045095191029764239245772425880432834133060767943430330986708758805367141513686852541836711351817918213669193","0.693708258349105549256657326728889025543270864262375054692701104148415191462911199278455895406020945849084609270795977977"]},"3,5,10":{"tau":[0.6374161558390136,0.687191135854336],"highPrecisionTau":["0.637416155839013581941680030552152993759031030531635121898165076648882802751056343695837483741736055980237474799374883697","0.687191135854336014333790484663973863575289861281338528311888959987943477576253128394870631541037038479727474360052569367"]},"3,1,14":{"tau":[0.8703147624402117,0.6247795146456709],"highPrecisionTau":["0.870314762440211637626250811882676154493827642454944162366334529390864600422346933031645795519351451423586456064262189971","0.624779514645670944362782072259282505215665300684791204866788069770842669047081728277462844254544067173159044432327974982"]},"4,1,13":{"tau":[0.8675034480173437,0.685543429099858],"highPrecisionTau":["0.867503448017343692479691990617108378316515944791667129433081714577316900562166921888697333109872582837967183902052462318","0.685543429099857984326006434464060771399176005392005190364251066431207479478514771177004446839647337681583429895292854343"]},"4,3,11":{"tau":[0.7369376596205486,0.7225104440077347],"highPrecisionTau":["0.736937659620548530883837576350452532315633054934988069113297249552541387300160745671691012475338771923321394680939975282","0.722510444007734709853613678258962768046244688962044579524257022938877427994415472603051709382630519239946032907855714257"]},"4,13,1":{"tau":[0.14683809076330906,0.6950059020454789],"highPrecisionTau":["0.146838090763309058720422094155752536820718799714322718468696570964252172478578928156329876763678403648370122489169924291","0.695005902045478847148492864981646472006762288395069263504001732938573904484151078081016039648619288368715804568151113738"]},"4,4,10":{"tau":[0.6775387078918872,0.7346051655760959],"highPrecisionTau":["0.677538707891887215633848443349119924813823241934677701843338883628044851798564628497639584090826122217658334908233373294","0.734605165576095876258088631913302095361149886627602225913581608382999909881547330027786693926097887301108501950269448596"]},"5,4,9":{"tau":[0.6610268754293382,0.7922885578547302],"highPrecisionTau":["0.661026875429338279410322436160010659425615935562713490529156194095979702937101175421219002703392462163623116323955916344","0.792288557854730229192226795674126680624806095379038560457802669823058537049818059226074804143659724164028335874310253153"]},"5,3,10":{"tau":[0.7246064419895375,0.7811826468905984],"highPrecisionTau":["0.724606441989537511426586509236626453015064683708955425341778008110038513017270354699270721180609452010471415956147353489","0.781182646890598367547042019768780949230953421137850883483061444062832110995092834226321038082223541788166272759462409017"]},"6,1,11":{"tau":[0.8568899141553167,0.8073659053407973],"highPrecisionTau":["0.856889914155316637701433929716594843380100780835715537270658186261239381529872110144739119479851551078833416268341821457","0.807365905340797317072241763732384076253293825998455198748899925894701688298852927313302961819646331661856666785349652594"]},"6,12,0":{"tau":[0.07052662933653298,0.7943811496380393],"highPrecisionTau":["0.0705266293365329718498358848402664852733634524935432750407883655859419207019342237594426234835027267684881478580654267475","0.794381149638039253567653954413965478002344419262785581166104892277788889909098173092234265786421195607451898942784896117"]},"5,13,0":{"tau":[0.06825172384081078,0.7323284828298989],"highPrecisionTau":["0.0682517238408107764409098175182161569805877540829285018336938885583101722190743598742539162854428873629654571815920987728","0.732328482829898965299004804881243627276773321302948322564339910303499416333696690626940014594350118787071805836375516826"]},"6,2,10":{"tau":[0.7806563251926182,0.8268054521702362],"highPrecisionTau":["0.780656325192618134977360463062550064892497379227935332509227763541467331069847493544245311936187037735897536016620535964","0.826805452170236160209380709188301555902411864946108698552387191721050463379954658863133832507783623234349098029146278875"]},"9,2,7":{"tau":[0.732569751126577,1.0260402141037834],"highPrecisionTau":["0.732569751126576911967039992691173043325833955978225752435851182576456896873651656224277281658126449845765218759999915615","1.02604021410378352327196971496460443262476531921785860535866351896288899916170471243708597730745493966494089937301388827"]},"7,1,10":{"tau":[0.849145373450045,0.8706657810081211],"highPrecisionTau":["0.849145373450045029291723153202454744455727435477341012392969129451873225227246006786308787102778570488796192592191038634","0.870665781008121083185373394507548090829096598521817441648161532090109098535161000216824556857460921132724032955929729377"]},"9,3,6":{"tau":[0.6421157739082675,1.0383390694793104],"highPrecisionTau":["0.642115773908267453594084215790506706160821271500307496699883513589511470125194633650021087223051509013452895234463504915","1.03833906947931050498979528591192395722281243008936816457220687916465835967262909922671216090621722310723465909267735379"]},"7,2,9":{"tau":[0.7676972114612337,0.8897558165910235],"highPrecisionTau":["0.767697211461233664858486601985996816926394257927982014780087751206062647051016247513121426597990485577373631090416304291","0.889755816591023542425957736072361130054481545228309310930652204999047253750147360372715098530524725126975554634408131177"]},"8,0,10":{"tau":[0.9374646490366811,0.9100210886737737],"highPrecisionTau":["0.937464649036681169514393981911727679740134095689272635929310270478199605749543371648156379993602185259656112054450919589","0.910021088673773764333858461268677310863564890228739196862925685914805662267148408161196232391534110944452213673941139762"]},"9,5,4":{"tau":[0.466133789954671,1.0463393711929572],"highPrecisionTau":["0.46613378995467097392572583384768088840998675137434059678717005791114419897285320098872402116336523421021328453567993812","1.04633937119295724566483509303572373054097262638721631189382642651976452947213219962082796992928970160632327831308418939"]},"9,8,1":{"tau":[0.1905946444968984,1.017705770211657],"highPrecisionTau":["0.1905946444968983875191723447179448528926956687535374887014051271595196693562163308172953271973433478564756666244239957","1.01770577021165691527747935698210516995819044945922069225547976045537060067701929448455468028422286967156252289546770911"]},"9,7,2":{"tau":[0.2866232363937911,1.0330765571560754],"highPrecisionTau":["0.286623236393791114724982292832903506726000446787355955100126587858072447007603875358251158833693050661672965779115563663","1.03307655715607548839328857028974353856913857579501242100728398012130032048883423638497200060979861978924986671479434797"]},"9,4,5":{"tau":[0.5538379016505133,1.0449912146764724],"highPrecisionTau":["0.553837901650513262279437860307397109897589901421612315561354019604190683234915232859804203337413880730736551827034416031","1.04499121467647234269010782740966516678903186241538425383705498012144041091118857035946093355339878828951350455676435843"]},"9,9,0":{"tau":[0.08342967983783792,0.995010927664116],"highPrecisionTau":["0.0834296798378379099310295314299855561436858563616790424185512201322484583093376116203526258156044230668719678678558060243","0.995010927664116021256403534435840828572442668606461136905813294957956903198375796821469179863959003768948940529218324155"]},"9,6,3":{"tau":[0.3776068880649415,1.042436471202376],"highPrecisionTau":["0.377606888064941512103827450441436488877135213848792464243293730218927805234492027389602856542251234450149247029344986319","1.04243647120237590270602976257671242130883632815481896081023416827457723022674779636352335906507056067898086551495327814"]},"10,1,7":{"tau":[0.8130158099813769,1.0829809042346963],"highPrecisionTau":["0.813015809981376934784561023796982671864720803831758007316332158480922852907391486386795152594125448797849097393910441851","1.08298090423469618769595609455377924730265179996412780112169892899687617117251319037202818413599943953717169569094952492"]},"10,0,8":{"tau":[0.9278205452088593,1.0542231347618938],"highPrecisionTau":["0.927820545208859274009290957818764123071164883699593410873603108797644386384336107946421357746179247776679079237864399599","1.0542231347618937084352273030205041174413867041623183031927233613657965998170348887803978126696118562141714855359657102"]},"10,3,5":{"tau":[0.6089471927876985,1.1129914403440295],"highPrecisionTau":["0.608947192787698536766131745095379568471708822105642941411625242983163206582309158047964868707679945696199074775886361906","1.11299144034402945606793919752600069965476638028738047734371297401450614459678848373539965677592856419925371390526469466"]},"10,2,6":{"tau":[0.7087408148382346,1.1016886921048759],"highPrecisionTau":["0.708740814838234619733620458414802523059623907400119067493200694990952907559958042524480405163800905685724697288569497901","1.10168869210487589800582336174519400201945161816170200971478030007771142384715629386386602744598021726318087333756073254"]},"10,6,2":{"tau":[0.3123436073314346,1.1084950596830052],"highPrecisionTau":["0.312343607331434633941795814277021581967888409352750061994006082688292457552874176801945898762420335216686915720861874742","1.10849505968300516737273616827285126408455521786005145050014108391654161211531220528643547479753186709597824571208648649"]},"10,5,3":{"tau":[0.41267725954145706,1.1163080204684128],"highPrecisionTau":["0.412677259541457057937854199527213802923537763321031410485434475494569868726758730357398885751095598143034372024040481462","1.11630802046841274210760822576076916704230259696674872870606168866526308677627227220699628220438787002592577795212553596"]},"10,4,4":{"tau":[0.5108979260234723,1.1177701179844137],"highPrecisionTau":["0.510897926023472330168971256368263729435167677463481368604313577851461390653781399768373591203360353804263792847943092186","1.11777011798441370587226266269439873907130295779830339498219363976216483676852536107703489454459770223203798074371418245"]},"10,7,1":{"tau":[0.2070520286422243,1.0937735392498382],"highPrecisionTau":["0.207052028642224293916879108355640152262893908615031956409117459030082613364449200001946261799603945519228543031656128538","1.09377353924983812962184060358171960258296717355311095401787542801352321604561038052842734355779164238907639038054331469"]},"10,8,0":{"tau":[0.09030099544116725,1.0706669198749894],"highPrecisionTau":["0.0903009954411672439310671487504302880987472730168465240126173050873566021795676180520850761817781510571675127140087958974","1.07066691987498944394010138174825697831181557898012130627678705854538570630357680359163723236266043707696096926608238575"]},"11,0,7":{"tau":[0.9210504906165637,1.135574673898176],"highPrecisionTau":["0.921050490616563728055843129274470432259689304456786917861824295015497139479235822472982891653407798636302202110052804895","1.13557467389817590961758487417988010610823149278928423418635975440316342966936095621996152101996917752047862627895080255"]},"11,1,6":{"tau":[0.794596738424933,1.1657830804054738],"highPrecisionTau":["0.79459673842493294079612312596813146606650109261268848433781935764660790234213815099184878104855850080929527519536330134","1.16578308040547391419904282868459595945228822895497702015384107924128743542214553065650336285118934049199161606099203422"]},"11,3,4":{"tau":[0.5672963347798877,1.1942020109487477],"highPrecisionTau":["0.567296334779887725147271804786999013957270135778233095030244829364529432906469189597256364018756161500983767031486396107","1.19420201094874772361128712915427087975172904291648065573802522556670289289477398180516620528514881616853607043773145191"]},"11,2,5":{"tau":[0.6788650234315609,1.1843380489812931],"highPrecisionTau":["0.678865023431560903872128520659132133238168249629006994766166306461487067765429392591354750118015923784066514519297866827","1.1843380489812931615152786531011342456618653491106523565767784587612218102476969760751835533982766402917753250308164107"]},"11,4,3":{"tau":[0.45674160265789027,1.1962719556897186],"highPrecisionTau":["0.456741602657890247896302246215658501923579839516751564097038072737644350858920189411923325397222922796275598863080840104","1.1962719556897184783142238111248596361152176883320050054582856403985992765499165295102938182077055214737563898940388804"]},"11,6,1":{"tau":[0.22780387232443094,1.1768970831893688],"highPrecisionTau":["0.227803872324430949616353313634710581829407231505481009170118066096144959397950135991695203279894722372139358531243358208","1.17689708318936881191865845558185309071834957296136229604755144491258523725109025954311494160754390265528319739524973502"]},"11,7,0":{"tau":[0.09906341371529667,1.1532693141306485],"highPrecisionTau":["0.0990634137152966721347318553793374519047398822120903392373648807827181735020735652078066958531108836504530502187210893172","1.15326931413064857762418542866097077503398769667429844219985893295053192819968898526482380711668625803654225935856590266"]},"11,5,2":{"tau":[0.34466594706622933,1.1906751582659951],"highPrecisionTau":["0.344665947066229338172909875874455061065930507915575502288332896275878746308570258639467312276629246926182941306406565072","1.19067515826599521850335591631045502165476053404963379254313538506835302475299595276722904788560927251142500408311270674"]},"12,0,6":{"tau":[0.9123609999982932,1.2258586170350723],"highPrecisionTau":["0.912360999998293258804518572111155324441698492121849463095018913366081627458714932045102040814252424802659870430716923362","1.22585861703507221438739350406585125374929686788203247104461118656095585468297027348583236392742131405998063007642788414"]},"12,1,5":{"tau":[0.7709916050562051,1.257910297664225],"highPrecisionTau":["0.770991605056205131582376178404206487480059932593027736748968981880153472291951507240841347197365403686972994881619846408","1.25791029766422480377237267324499797367194344815174510987808969551874612513784225005062468601388440452686478636337072537"]},"12,2,4":{"tau":[0.6405003356839093,1.2760565169637037],"highPrecisionTau":["0.640500335683909294464237921773644962606424453711365402458680934863174630236762182937565488301621036822050444702723534876","1.27605651696370360008823574421471771417934711000864389442746937384665611757437531304291092811612519604361765367996632531"]},"12,4,2":{"tau":[0.386239061548342,1.2814759090709869],"highPrecisionTau":["0.386239061548341973882722073408686894083980050206448885894605269724748939777598384631777423708791876008363408952427642574","1.28147590907098675016940753309667496551493961636214146301750527567278366867619402182393102704743875646878416148778550607"]},"12,5,1":{"tau":[0.25451667566624775,1.2692478113919023],"highPrecisionTau":["0.254516675666247725686532518182564870338797817555784158392079007092647757295877300326800714533367795658538802187680693758","1.26924781139190244551123106617436710004675069468325394170571923867244520097180521662740481158195805686774950308386521037"]},"12,6,0":{"tau":[0.1104260758825137,1.2450860287274108],"highPrecisionTau":["0.110426075882513703728657729131024799943755430556640719476681560116700862018203526544860562114155853470416110258556985567","1.24508602872741072685230121048212581782247262869971741715144251690215602804245596142871951962127738706552613718082442955"]},"12,3,3":{"tau":[0.5135572759608424,1.283638414610501],"highPrecisionTau":["0.513557275960842412601355069900766131998435981557801731326986641437896468510572205764998341415446356835593155312642359657","1.28363841461050117524896670475879347331105845075862574969397538026916306510521077829109830336302247086944472362420387913"]},"13,2,3":{"tau":[0.5895358362467908,1.3796518450376425],"highPrecisionTau":["0.589535836246790837835435140147615090662226361888628296580219187302867252898808785407431949254505407351494180743597198003","1.37965184503764241348761013312885945170897113484466962613303548331529864724882640346254355661974735249509785453905308353"]},"13,3,2":{"tau":[0.4414942797512834,1.3832798472601187],"highPrecisionTau":["0.441494279751283429498212201827651657745723949758344443601064554412930401911367978641290536005560595720096759351119688014","1.38327984726011866236310679913484244797756027032710543500191090033622605527886549018022012042098668375628129442084896518"]},"13,0,5":{"tau":[0.9009061129213674,1.328224884338827],"highPrecisionTau":["0.900906112921367369435707400029388639236941735083969722517592673151374702159489069237710665094105003736802269150793718113","1.3282248843388268789228457985577515889750211043696095124132952193858546177451522684321090219273155486566551011040037743"]},"13,4,1":{"tau":[0.28995937437061764,1.3738371877543056],"highPrecisionTau":["0.289959374370617633761492390299733146038629840108854253534725081062170076988143371697544093134309926868714044190678759996","1.37383718775430565815495033776890124472099010981566914705224330036465070092554206444086933408053341082864681450515757708"]},"13,1,4":{"tau":[0.7398036728272641,1.362546978047142],"highPrecisionTau":["0.739803672827264098425329901315975612142606861510467294799761336950324049597316501365417937709414064806317663553393136641","1.36254697804714188308262245831423281060948303361327072230333682105730069954040189079179225750048090435862548344043822166"]},"13,5,0":{"tau":[0.12557048382529878,1.3493410192509043],"highPrecisionTau":["0.125570483825298772985712094602694768609075964670861164397302088082047006698483548351999774405449369801808705824422436463","1.34934101925090426195690912806015762442848673092079262282027347795847473469043894372070191937369005090047910237305742653"]},"14,0,4":{"tau":[0.8851990349478351,1.4474668505358033],"highPrecisionTau":["0.885199034947835106738435916368947567396735009765874471158997463646537630253240694354482499421208209454033640389227934469","1.44746685053580323046562707456198248109306006408374105903857039627710269907425783270834197307233032417789479631635572266"]},"14,1,3":{"tau":[0.6967527540007334,1.48446374301556],"highPrecisionTau":["0.696752754000733435866827936161383793215993483545044961933669234284274068022531591027369277255073050727560410235420096362","1.48446374301556005218406604315739096784810844652818455729719188419706538318627364387491586584817947832839830381152932129"]},"14,2,2":{"tau":[0.5184447659152044,1.4990099125451166],"highPrecisionTau":["0.518444765915204434260649407305703327501400605723453524646575405890305636290233222956863089669832635497534886657204578634","1.4990099125451167292251232663311762341525610025219571528729611204588534810882185398027982047196054966289824670701445846"]},"15,1,2":{"tau":[0.6333398553959009,1.631276575543773],"highPrecisionTau":["0.633339855395900921988562874077239540218386323263985554777904158388673953123351580798711133394066999295553193151630127888","1.63127657554377301835954304477933060081412859068412973448589784564191097497562576862931757795111355310636355484196523295"]},"15,0,3":{"tau":[0.8623666127112459,1.5915015239581911],"highPrecisionTau":["0.862366612711245960736313839887752899444080067372538185636836776097620891210610767395742100658709666522602484620575403065","1.5915015239581912041497338682990931709533041273951593526964036928405178202354424207714418142782992657384861870249609211"]},"15,2,1":{"tau":[0.4115654042272635,1.6392931045494745],"highPrecisionTau":["0.411565404227263494376161247165046666704705568631005073762127762340818441415085467462014052752266685023037967187272267082","1.63929310454947456974666568877195911126265806658732666469648587902933752686250458724321344854714567301042405003533747721"]},"14,3,1":{"tau":[0.3390669803303544,1.495021060655508],"highPrecisionTau":["0.339066980330354354791425141407498230807726516338688296911902969309817149654298177335213224969590103072416659483944657675","1.49502106065550816141363961834438859181847138805183456426911848563560022720399511788976994128123921789283872515755296567"]},"16,0,2":{"tau":[0.8260070805513209,1.7749427987758235],"highPrecisionTau":["0.826007080551320891070362789361088014246270056388057941812239090624752622290330792967095783168801648718200372631219352265","1.7749427987758235531049472435656262043875206780076681458658258394814822404913043764213957597995515362376099709621909155"]},"14,4,0":{"tau":[0.1465971280574504,1.4709182072492262],"highPrecisionTau":["0.146597128057450387230092915147563974930616558420772557052803350087681031173328447908954967644060609731657414199362852232","1.47091820724922610009217080363155381885316243959398204659656379140199351941561717737211252704970754929673890939679200369"]},"15,3,0":{"tau":[0.17761680139433272,1.6177981466876377],"highPrecisionTau":["0.177616801394332729371350131238922108135086547169984316237876071209060714763029397511968030177408478896896880501952663435","1.61779814668763755067996855283524301062481616052980924591009074468522175692037681495629896212355242690183603031507511809"]},"16,1,1":{"tau":[0.5297000598753951,1.8159616064482127],"highPrecisionTau":["0.529700059875395146820567767808971183651929923780889606701955818659922018035551558054070646578017325417060006462352221721","1.81596160644821261708863706729498379640672685362916026697861002613839275997587261026833106009803086522786271115353066243"]},"16,2,0":{"tau":[0.2279085592837334,1.8043535100878332],"highPrecisionTau":["0.227908559283733411574695263325971600925728241269855160981927907100718342719159571243679128882289157519757246497264179504","1.80435351008783310866562244167176348389122792541623518702632567874482582443596635735997069634721097864253530900629775359"]},"17,0,1":{"tau":[0.758116613745831,2.0297133538955423],"highPrecisionTau":["0.758116613745830982190387971781361306473131069353615711224146433732985124621002236659571969073296971650941899084723920696","2.02971335389554214855686860229813673279744009504214822301085685201150593996134947388895559111988337017416517464097802495"]},"17,1,0":{"tau":[0.32377062517567284,2.059972029485868],"highPrecisionTau":["0.323770625175672860538169352443640740713563931186886217744950265980920137862583620087274996172833022673919404154731107298","2.05997202948586821201896486308493768536774590834129717560225639990487688709560986490688536461472439021115413274071112941"]},"0,1,18":{"tau":[0.8618026218720678,0.41098380857138184],"highPrecisionTau":["0.861802621872067830921156654203973178998190322677587491088037619599070641067766572244711340260011621240572819644957913086","0.410983808571381835981463389931952808951023027115968810289319618543392023167434819320461407007339056395978054564490385321"]},"18,0,0":{"tau":[0.5781101417761472,2.448987981580475],"highPrecisionTau":["0.578110141776147113392086841102085925343701474024295898408644444551517095754050933824807533029152818692574429706049458022","2.44898798158047506268271127944527260375540035167237941261126809632907320121419687088420579029949454435782838168755333077"]},"0,0,19":{"tau":[0.9319867710015808,0.3590202569846263],"highPrecisionTau":["0.931986771001580728285601951953675791646358366781913594403891551978894797998191298604297937940242433060444103024069706955","0.359020256984626288188496477374917205028617838074162825044452474328771478583686688445724307403323237207504705768530228134"]},"0,2,17":{"tau":[0.8090078690383319,0.44500228984488116],"highPrecisionTau":["0.809007869038331829636294709497007148208205757051968052442637125871199952009681120296551867369283476744027890337544723519","0.445002289844881139240307310854040601304263999610684677441107594003836738504864857527766054579304991296502956521247092825"]},"0,3,16":{"tau":[0.7628812010433695,0.4700453331297952],"highPrecisionTau":["0.762881201043369515836179482254559724086801924602515928943494276512051674067705906292204684268693927135439126777471478935","0.47004533312979515919864676588166576958174534611530720850771063635751116558014051867297665147495321298822926280906322988"]},"0,9,10":{"tau":[0.5272673450733685,0.5323461324120431],"highPrecisionTau":["0.527267345073368492566185002553099385929977695438208441060729796737804073689755754140717556203434713719221459568869379099","0.532346132412043079858098323287946600614500867455208657174669793961834607450509330832523568757047620945561210971621750917"]},"0,12,7":{"tau":[0.41492516182133277,0.526728278604754],"highPrecisionTau":["0.41492516182133279268866361169960052833700446920175985274933436556837925900660844554625380578219725506188590832404267751","0.526728278604754059852577569788811109733294125836910704682600108721387500811915904518316585589724661284037563770003672668"]},"0,11,8":{"tau":[0.4526102415226841,0.5310857438721439],"highPrecisionTau":["0.45261024152268406020702122424818239787447511955686160336803509724006251367492075843248148732721002606643225075889980488","0.531085743872143887342539116916455367223836055647870601807458529550449703025763376544320731527572744678156269803105132842"]},"0,10,9":{"tau":[0.48998779542796955,0.5329364676820303],"highPrecisionTau":["0.489987795427969547863588849247184515760929081806893252178447607460626265831199816262546442294969249718770420963500434936","0.532936467682030315518830346112612749679952196917566061415375371102088576099425585075546582029544328291630752915191540307"]},"0,17,2":{"tau":[0.2089767304222114,0.4562730792254419],"highPrecisionTau":["0.208976730422211406846045569038752254573885747278043119367042925964486585126071294934269772851218560541076637886036731837","0.456273079225441894602482388563980832963796617388071567634067135735385527820029599252773492586956015528796082861651883388"]},"0,13,6":{"tau":[0.3766886514721084,0.5197113431036531],"highPrecisionTau":["0.376688651472108375534248523796340529488700870215970498566630192958427554421495751471640663301068840709705341792557955165","0.519711343103653023181744081217360265746967662190079870418853310130344458600159349043847912441553853215517261629081958842"]},"0,18,1":{"tau":[0.15702048589543108,0.4254809630147997],"highPrecisionTau":["0.157020485895431069067681862496651769608250806715323164654485202026220616168648231250246127461015320467723523011414387888","0.425480963014799717486709950355736153965174143251284628741117483299153797182187299792797928593289595197592156717296718551"]},"1,0,18":{"tau":[0.9437627583813853,0.4497879670397753],"highPrecisionTau":["0.943762758381385293077313992659869685431003024864912008695597382595748259545794440880142363016790776626373596078910373771","0.449787967039775296693327883059645724043705540556605709824557919922591236497458420341426765104505164461511736514108939894"]},"0,15,4":{"tau":[0.2971713723393731,0.4964592351716903],"highPrecisionTau":["0.297171372339373076841771414763079828225493638111835586279949908556172336148707561267239162014689605084244116654496792963","0.496459235171690337510686908552480117515500763536248015588108283653529332558728744101250598759557215361526358828636085074"]},"0,19,0":{"tau":[0.08932871471129729,0.3797004964548042],"highPrecisionTau":["0.089328714711297286994172187001019179275899346658384785717147059299907981349766155376227264580858987803659985792783647922","0.37970049645480419306555678815847238931842324024001699159819339718679912875009230469900706445239410930768275628866309649"]},"0,16,3":{"tau":[0.2547325006976232,0.4790581540363676],"highPrecisionTau":["0.254732500697623193963148412403703744126052737238448633833438800217936737755752230356197499347768963888074118423189314063","0.479058154036367594064589313923400198970495069863841684358649280879010553457556732392573293904632139748561096968740432655"]},"1,6,12":{"tau":[0.6319436728384754,0.575071704670687],"highPrecisionTau":["0.631943672838475412883635081575010826735404451498256593250460349626364497636544312564726831126792677874977769274958276741","0.575071704670687054710937982017846997948123284480345683568843932193649525288190000320214105768652335138337338891960903902"]},"1,2,16":{"tau":[0.8156573260819644,0.5148587262615589],"highPrecisionTau":["0.815657326081964356845744052740272569566665094672260349446851573158245920347623896430465985102997081673959912926784337283","0.514858726261558946157347267307363877626433891230452783319073048577243313611396637701347121171842619203877904497664865915"]},"1,1,17":{"tau":[0.8723226418258045,0.4873944800475971],"highPrecisionTau":["0.87232264182580457769084797686648000637030023504254100421751803701039200188722982914807594161308723622886315060609388055","0.487394480047597131518930140332950448682995127501077963541240576770455207648449289692913592543065776028847924021250383231"]},"1,9,9":{"tau":[0.5076505860332955,0.5879223610897591],"highPrecisionTau":["0.507650586033295506113987899831157013277702016300104479209801740011936751265642940874765386666756072065198775661115874467","0.58792236108975910546683961166996853245964911599880136211418411412293128132581761068415951058412438138630158141995367646"]},"1,11,7":{"tau":[0.4252282209919401,0.584055482592392],"highPrecisionTau":["0.425228220991940053199312747811642826838214641300012709395780408175758304254745033660136449564968704225407054609693788656","0.584055482592392023759891423109253268739607735608791760264431855447436580927972772532270092345769222138912237379166104007"]},"1,10,8":{"tau":[0.46657395748471703,0.5872265075888783],"highPrecisionTau":["0.46657395748471704086790103648232076000678113373498351670078948290880247942026458782674102115499557638402383468987500575","0.587226507588878307642515222422862377885956444076799064741448207447009638546790219507239740193286980133436334371620462797"]},"1,12,6":{"tau":[0.3833185273415113,0.5783103187975565],"highPrecisionTau":["0.383318527341511288057397801122075563239837149274201714064713392588715758298024008434280373167672738108981604757684889693","0.578310318797556492560167887647111311145013646097255835341556015078211439918390180514762697123410495138846352241662698604"]},"1,13,5":{"tau":[0.34047830185940836,0.5698042762468476],"highPrecisionTau":["0.340478301859408342997605492991595363319347511433801980400616572361742354689140255577581187117438439057909592196855291037","0.569804276246847649737111176025772887392573967696944628813710460007911133420031138582313206904368995487752529032800984335"]},"1,8,10":{"tau":[0.548722377699443,0.5861693516503992],"highPrecisionTau":["0.54872237769944294124415179759973253136615272799594725286799765358219607075532014304353692594639939129469822934664816005","0.586169351650399132217056094253526690756901523397862165717255192101419803986020727704355854098043017388974566410025493256"]},"1,17,1":{"tau":[0.14319600666232374,0.4983461116977921],"highPrecisionTau":["0.143196006662323728710824964095633741012567180101535742216263183279895439427479136088869440265064547465484644924948532422","0.498346111697792108498278439795706001213296687110312753716020666600350396037107927035354626436066407775969474743463016892"]},"1,16,2":{"tau":[0.19964909857256055,0.5236092625381856],"highPrecisionTau":["0.199649098572560551661896941287880119656085287136545833999030129515961675844697801500270840520714444908754055683712022006","0.52360926253818566612753345130530931259721051204111396797011446831212543159126317167908722316633146555089784682368135622"]},"1,14,4":{"tau":[0.2962015780122015,0.5582307611278752],"highPrecisionTau":["0.296201578012201475349601895548453971800054820662706761180843000563303948018702815502334941196293653392529834834231431476","0.558230761127875120359573418476878331005853682712861928675254649212961631394646734132368333775767121274780372335158581516"]},"1,15,3":{"tau":[0.24971059665572495,0.5431000579010763],"highPrecisionTau":["0.249710596655724953993365546328626024336649053877385302742951262877171527919873654376910448910512909130651586533026708751","0.543100057901076309878940558940627821147749261877787205826375374546390436358440303292591135718789997673799946879997777559"]},"1,18,0":{"tau":[0.0722936197427982,0.4645161327716073],"highPrecisionTau":["0.0722936197427982075844812015151718468595508962004201291548290167782826282475813720060487705299216500827171523870035416449","0.464516132771607301330470340212831823815995030754252703793722620111291402517731973180641624886877990135131483260522058276"]},"2,0,17":{"tau":[0.9476136048792818,0.5209173511804879],"highPrecisionTau":["0.947613604879281756995000629124089773832529495162396630937494209667846821020708612710588159488185091100562251172069244603","0.520917351180487837366615373539431534998034938802077970923878244820291492335454384986047451446305516843248493415829688455"]},"2,1,16":{"tau":[0.8755338449035339,0.5525710911064837],"highPrecisionTau":["0.875533844903533955016305815203830821399583493131155381245671082903159394870077700056526087103706202860410110298240459663","0.552571091106483697995591737561575322325860726153222994941689337172586683564908786274242824758428464262540658163098877878"]},"2,7,10":{"tau":[0.5744393487881617,0.6354842627971987],"highPrecisionTau":["0.574439348788161763750803953952148558039401963283014748120081683808520808134359792845695261947071397267114162547010550169","0.635484262797198672162081947228684632331223464500600123001408896501216964630728386123523740797144984370145444962954400328"]},"2,8,9":{"tau":[0.5296639431453438,0.638523276061053],"highPrecisionTau":["0.529663943145343791159108929395607960273517317771276505951224519904266703304983787135409465321323343907916118041299618474","0.638523276061052937883912152329277570930755478817738550592827847945713549995277981901953307334121336432835281322618577174"]},"2,17,0":{"tau":[0.06651320220171264,0.5340401192730955],"highPrecisionTau":["0.0665132022017126339345175475620324071872740681194095260069358436958814565740497680315848039664085297985217196472775424694","0.534040119273095455322534495756782554339564246561700506738820002226866786660687551640412444623517044268191970429623205824"]},"2,9,8":{"tau":[0.48503486340611934,0.6390489565122712],"highPrecisionTau":["0.485034863406119359377483642894885512699704764384686641289025423530867725865209906916038932861821917730003799989582006742","0.639048956512271256091369502090353575295388648011123895134372170523726545141319269156760995990983467450868187708758897177"]},"2,11,6":{"tau":[0.3949429358404746,0.6325248427673797],"highPrecisionTau":["0.394942935840474584984043597056778853359150394218237800983530556391591878083683983427442735018043900022984564025244332252","0.632524842767379730542327132774020005033601134984638449693761678252523756065042875541983404123328457553133040545985560896"]},"2,12,5":{"tau":[0.348753239070235,0.625284668946576],"highPrecisionTau":["0.34875323907023499849566161133906443316536749707052214330069931559869413692476905550959393412179832846681627827213985884","0.625284668946575987945591497420328458865593635677212003298159674856634134711206844302130198558971857081159090238891975582"]},"2,16,1":{"tau":[0.13870933266550858,0.5623843572284862],"highPrecisionTau":["0.138709332665508588639971834015894366698549418566181092689103269335933026505174040751075192488177642418901846557256025287","0.562384357228486201608622165928952417415228127505101287228418616927177335421545293197524729054852974772517562124983322537"]},"2,10,7":{"tau":[0.44023858065751753,0.637068915418055],"highPrecisionTau":["0.440238580657517525261651999788091743196081696892696848414379749516038635326210900280257099869004424739354700224804831963","0.637068915418054961178381004690223817637139789758591568777673438588777829174540884454710266788662336222654378815826709141"]},"2,15,2":{"tau":[0.19805497519092993,0.5844081046748366],"highPrecisionTau":["0.198054975190929942749654447412098418660653230159320492767708822291661046955426940391426952684750593496354400182206270412","0.584408104674836570306318940729787465627574211708898366341166745772832496609680673539108431722392862175415131042400250624"]},"2,14,3":{"tau":[0.2513457921468794,0.601684323417455],"highPrecisionTau":["0.251345792146879418157608015825203749945906452769198482005051818111457305469622402670058412706824307792583735495344367153","0.601684323417455008680452824539737604524293296356187991180751380846146335458690243905821295130302451926691701703000681712"]},"3,9,7":{"tau":[0.45893590134021556,0.6883901333649718],"highPrecisionTau":["0.458935901340215551817199229551174230050098992662498152792043148828267879726619783830584818402587314435758651071031502471","0.688390133364971782322412874709037788485279505512045349225639692380501689519100369559607448256813820298908896963449662277"]},"3,11,5":{"tau":[0.3606430528523792,0.6789552323364029],"highPrecisionTau":["0.360643052852379203387986666917126444564775912968045796487866444090340771504501001889848873222968259303397224118408964376","0.678955232336402899868529202544861856680627726417709845450113452707817046877997909789049967676918706341541347949415347113"]},"3,7,9":{"tau":[0.55562303649994,0.6872900572271161],"highPrecisionTau":["0.555623036499940102726268303199752320043589393699863772242821257852413150990730022453780667450535839722360639207737615969","0.687290057227116083722074085728654435217147693651758458513791397725719625988168771145148280528987132626829738251638479431"]},"3,10,6":{"tau":[0.4101949394944247,0.6850269942931958],"highPrecisionTau":["0.410194939494424696034936809650706997137932610829276230225773606735547915780912294988806203394437109192130483059547799085","0.685026994293195777267869084989920709891425068877761923703401901147285387112270294240243418365563048976067195289664174562"]},"3,8,8":{"tau":[0.5072876012234955,0.6891364593123629],"highPrecisionTau":["0.507287601223495432454104383604898697226988876351855879938690402222412234531579249236242434191537278464588171199471137163","0.689136459312362865507562820294011564991881455789937068499416509620971723045201777919659877375689970982824022113778783395"]},"2,13,4":{"tau":[0.30114591785952693,0.6151232156211269],"highPrecisionTau":["0.301145917859526949281416857182660236027793359134141621133414559031631206660484941954920946971430387715578829290526300324","0.615123215621126887495196740652494768365647630003290312106281537130132614447539854588532244799152873864598528155326340461"]},"3,1,15":{"tau":[0.8754519060217141,0.6128785846412682],"highPrecisionTau":["0.875451906021714040540120969351344297494356264333982897340421931254982965728539762772936920971973224742376457905340424071","0.612878584641268210280422162064261626808567100296131447023157465511446539255504224654839395578325494568123437320264689395"]},"3,12,4":{"tau":[0.309743783506423,0.6700031754995566],"highPrecisionTau":["0.309743783506423007000384326164973496049495583003793854117046824820254517975341405644629208032592117842670426658982500367","0.670003175499556610777837132579234965532672329502259819253531488377043559828646152009425558841256543161357850688528177354"]},"3,14,2":{"tau":[0.2003620013167267,0.642196791572032],"highPrecisionTau":["0.200362001316726711805194452810839425131135274851018663465707534652654132716950162814356676177821552599975498223595978653","0.642196791572032004919932217583566862755564717368755068198664557254102610476826441565373553697781978369168160016269572707"]},"3,13,3":{"tau":[0.2567306756501352,0.657893878048059],"highPrecisionTau":["0.256730675650135178204328386108719216675363572398565824713543038750361010404309303720661622722130101126698938408545605071","0.65789387804805896238007278435919647097837219317831835540633543509891009141218677757232819027262845839246852811368133698"]},"3,15,1":{"tau":[0.13826656769720103,0.6222336052106329],"highPrecisionTau":["0.138266567697201036340597987664332306313321336730027868516211201132109635800771683532799152962825003837347095860984443086","0.62223360521063288381258637620568295744857562125994588351671528877161325934064990511170588698568707749493974236590206172"]},"3,16,0":{"tau":[0.06438939895788084,0.5968734205879668],"highPrecisionTau":["0.0643893989578808486400039683619477812243651871623236660541193150891933766215443181415638956938183370598588890575519027046","0.596873420587966804817635008786065782628720322057430202300292074658212730350994442365651751724112676108768137893183215791"]},"4,7,8":{"tau":[0.5334522611085232,0.7386043529863554],"highPrecisionTau":["0.533452261108523192519618888648930369998406037046136293454397057917323093897900366619230200166145166291281769857707548337","0.738604352986355404819630786531912548053784148076532189868856094991304457738910895345043056784047921301279033585005962943"]},"4,5,10":{"tau":[0.6388612258130824,0.729127948798431],"highPrecisionTau":["0.638861225813082385485566670885789063237991524248105210120458516294726850842101645212200637616667428987572093104479164058","0.729127948798430987993279020840368701727188381201208640138618158426753131741588267124640903907720889124704568897430397184"]},"4,8,7":{"tau":[0.4812620596853825,0.7391946570181567],"highPrecisionTau":["0.481262059685382493750101677817777130368969810607687897946240176367918796420803960447816158744596260088372169393867329694","0.739194657018156632596057094368105152658758699094060438428615642163576739045450061879992043631578281328202531788073275901"]},"4,9,6":{"tau":[0.428835294855052,0.7370555884880153],"highPrecisionTau":["0.428835294855052021472387107386957511643477402217772057677080282791934156233425547396389933879537200299030098575719719859","0.737055588488015300464467438013036248364278834984195450185589780081706567355586319776494562426767663382703574915019700861"]},"4,6,9":{"tau":[0.5858333865481948,0.7352751543600502],"highPrecisionTau":["0.585833386548194807417946491605835784437493124346835578899244115856080961058495315286895667323297273619098801153033409583","0.735275154360050183282919771440909734038118577815556369184138329043178321366861854764586227300891541676084329375198796598"]},"4,12,3":{"tau":[0.26495955256684406,0.7132873617357153],"highPrecisionTau":["0.264959552566844071109583755872742464317463447545331435699066149770518589075729248310802892837042280531076318802484504367","0.713287361735715317118251842347622641246014023269131297497949168068972698132198721893105425826742665852201257137163707431"]},"4,10,5":{"tau":[0.37571489164872385,0.7321290035231943],"highPrecisionTau":["0.375714891648723834103370629874383926681965927220237204547140111883721915894533319172413415657373810115809467516475848482","0.732129003523194323839913102653006090731622260105711706283444063280218788913418254979346258433459333816527992938851378663"]},"4,11,4":{"tau":[0.32134507379108673,0.7242818358062783],"highPrecisionTau":["0.321345073791086722794179849951648640023126592367207944694524752640910022951010696989831843919452125484616867425400945813","0.724281835806278207511657899605701743463491479072078108988563865373249180968766658567838446203802768197001424966948108673"]},"4,13,2":{"tau":[0.20535178120318084,0.6987877132371902],"highPrecisionTau":["0.205351781203180855855628478901196931398339311525159993609486747994028282308960557450055043667078577370839208593129986019","0.698787713237190222726719866618534999276624946153101138493342808135178554357464078738783815895538700660565319068980059543"]},"5,0,14":{"tau":[0.9479909558645581,0.7027129451343064],"highPrecisionTau":["0.947990955864558015985580525363837645921076784901546187612806914586240357080248437763235869965196008300102495086318479826","0.70271294513430642601414370905047204285989631438917203765497487984514517168703581162076643910410782895643889571578152693"]},"5,6,8":{"tau":[0.5640031480631195,0.7880948114429046],"highPrecisionTau":["0.564003148063119492220689045366737785922672052830336906161861306345238324563896985841687861022553513532306341228720696995","0.78809481144290462334002815684120257924112799014451760705641985702971039515414988772641877580366029769641765451291569946"]},"5,7,7":{"tau":[0.507533845470928,0.7901981821293996],"highPrecisionTau":["0.507533845470927899349692119776350112317135989427218920741915509044382119186233397438454086699482459884189860721579534041","0.790198182129399616997892237743113303963452043035383205272697698628910710364104813366471570673108672263033820893396571812"]},"5,8,6":{"tau":[0.4510341220398403,0.7893946018358806],"highPrecisionTau":["0.451034122039840293091926602422331143685940694218134709977205361497769293023137211947462254343828038471707555769063284091","0.789394601835880647829600760044903789722657006275509386044697488933001017908387646644262127048533667941347200718456559011"]},"5,9,5":{"tau":[0.39399500974539925,0.7856591688511866],"highPrecisionTau":["0.393995009745399252752321222950656518987717849149501849106644085396580895957793540149878692210351660271535690788441557641","0.785659168851186588973380472577368042959898664278033477041928206573321539942986309654140358435922422096803812515838153291"]},"6,0,13":{"tau":[0.9464371340473069,0.7612127487824928],"highPrecisionTau":["0.946437134047306881912220788410412129724343048318989859903254369607982894683614100793354068873817163075897410144239989097","0.761212748782492801207015594446776428081798573876865439625550443511775192420037537579931475636147323216975469186242589925"]},"5,5,9":{"tau":[0.6209497745406859,0.783035684846519],"highPrecisionTau":["0.620949774540685924620129837285922838202432871209710094171673388116596152648544096466887814017592249456722409360312195696","0.783035684846518954059709727473487291659997123523463254455274802713469344943163699247201794844477583398295473385306897647"]},"6,2,11":{"tau":[0.7927513476644281,0.8057183961078314],"highPrecisionTau":["0.792751347664428129016392453734013241639895884781684864004020461175916524079507023966133996466151794508921846405849363533","0.80571839610783140026508147590197337252000007400305957200882650915674948844974256450298932666476314372373343865846224771"]},"5,12,2":{"tau":[0.21258102868149392,0.7553620264976193],"highPrecisionTau":["0.212581028681493926289328750087235786087265678127427018612028243239823711703106370955693760186848531643581154923571205205","0.755362026497619246650559816138236595607093365076483691981828638466668000328209917253720354951337479479127256636088506129"]},"5,11,3":{"tau":[0.2757555729174594,0.7688968129258307],"highPrecisionTau":["0.275755572917459395180683010108293942051596369224536930200539086249486027241333979858709223888896726975032037571477960027","0.768896812925830675446965252478124687542687910656469382595473921516364948443480415273478887480701919070829012021390081339"]},"5,10,4":{"tau":[0.3358287548413494,0.7788906170302828],"highPrecisionTau":["0.335828754841349399235071851635887537321681334705735630666704797534594063612984274974959212142267185325745655495489462917","0.778890617030282861455871482456952322931766063354339462713707956021610871963337099721232230948888543063979668785646388864"]},"6,5,8":{"tau":[0.5997209873858261,0.8380157297861537],"highPrecisionTau":["0.599720987385826145466865435212485617389990696766713477032282644235193217885186380628886872897593409891361947525719674472","0.83801572978615370791997289629640392733648115424988492000592359976866014320960807358063689138804829039271771729797675326"]},"6,3,10":{"tau":[0.7259289491742442,0.8202194033041901],"highPrecisionTau":["0.725928949174244234636800897392985820525844450489793732467932327529626833934003031215318392439392775840813217902869575168","0.820219403304190129568791923168745739475140333242942624816116269022892624632777966149702943866895629154945252493688756068"]},"6,4,9":{"tau":[0.6619800428805229,0.830858306002961],"highPrecisionTau":["0.661980042880522972665861277283587728849360012535739579193826137586209279997563889822967736645700822080855860120410913743","0.830858306002961009815391811977389888253236067523704879195895471004222247145732122632363017904967579029257171782234945432"]},"6,10,3":{"tau":[0.2891523958420038,0.8255338265041358],"highPrecisionTau":["0.289152395842003832288870742293515334161970221412116053671671157907297030803629993041692823356034155930210214330388228402","0.825533826504135731975986942430311504289501240658739186087604686527825767000885386484707557454557205782992588904390977976"]},"6,9,4":{"tau":[0.3533613627930646,0.8345622410413926],"highPrecisionTau":["0.353361362793064604460681936188593251159449725830457681769802999927914780124402407293025940714155786986916048189880420678","0.834562241041392679398855331582265409538256879661138898245110397038545605439792387381854047290004373690005996370405371949"]},"6,8,5":{"tau":[0.4157854741658521,0.8402062008316801],"highPrecisionTau":["0.41578547416585208059137824332905632100639508755360145485305769353824239628645943496146842440282089239718519766902930219","0.840206200831680037333291534986715957687466837108500143145726310868429995003026994613305493067493261514632250410734032395"]},"6,7,6":{"tau":[0.47723565271910323,0.8426306457234303],"highPrecisionTau":["0.477235652719103214001805638024478742009098426189661057610434998153647428538241747811652370404549465021093283629023523066","0.842630645723430273402559524659126057241936358549105087536325151915539927684888885768768984578838322536168070809812272388"]},"6,6,7":{"tau":[0.5383511885409346,0.8419053982999614],"highPrecisionTau":["0.53835118854093455487298007699672020893781612612062525818308573076574837706749632175889402196723057938348017488792789284","0.841905398299961305504551837538520098276769240516432200673191089499398837978639149788025160495128045134161732867757010705"]},"6,11,2":{"tau":[0.22195292171209294,0.812830262286683],"highPrecisionTau":["0.221952921712092926787578113433736823823535008689970794410045862255557736328735772017119440042614310930758187915735104777","0.81283026228668299144925038187512551830101332092677504364451468817764282876140560802896999542731980213672806863354074639"]},"7,4,8":{"tau":[0.641760159821835,0.8886007320034717],"highPrecisionTau":["0.641760159821834980076426295588963148711930863811248282868822234571087443046018029317473482234418517819370837628527490423","0.888600732003471641511529905973138389821544312972035988825164094635329163245831706072142601087995432056229947387404854046"]},"7,3,9":{"tau":[0.7104168714578707,0.8787883642867614],"highPrecisionTau":["0.710416871457870674893237381861394371455279392036542254522845470889582842627969940058028598064655842435913453580100793738","0.878788364286761452993505676762075824241381680571395475286130541508542855205331306193115877358120627399780369793347535237"]},"7,6,6":{"tau":[0.5081546471275734,0.897246296983203],"highPrecisionTau":["0.508154647127573439696507784767059694011897132229045082554181130069111960159896473014736121468375242270223819607307407567","0.897246296983203011365243619008419411491593387987162762307527962888228840414990255189718465229272855852548226387698902962"]},"7,8,4":{"tau":[0.37433704736118384,0.8919486393310695],"highPrecisionTau":["0.374337047361183857235789176089045950129224973178141154385702992802769442428032994665509851112071209273047790668503737735","0.891948639331069532148149031899367960880042420138577979118700307199521475362423829685871671844208943345736334718832967069"]},"7,7,5":{"tau":[0.44163285197137053,0.8963432078348261],"highPrecisionTau":["0.441632851971370507624086374633920001057304749613803092519031937718251668279406382029347402182785936126849675741373558333","0.896343207834826114577568417073450913179251909501854931867352366332690298171996023223548730777853049959140960776707020653"]},"7,11,1":{"tau":[0.15669576672398067,0.8556715555238423],"highPrecisionTau":["0.15669576672398067129216985366279560812506775936710570182754041519897987998661306535090140680288021790319118086003299997","0.855671555523842305153526585305183950658742445551550041208206088730109521878733081897217465130298753134648427176851193111"]},"7,9,3":{"tau":[0.3054013944468967,0.8839241214277639],"highPrecisionTau":["0.305401394446896671272330059288918735464868121247790330689364207960540711529140308835005337078387099101770496853535038473","0.883924121427763925383080098241761117141969390124622028122707171980967591415799060903233627526985245435293785296833144344"]},"7,5,7":{"tau":[0.574624245710337,0.8946892744385897],"highPrecisionTau":["0.574624245710336990154446126464614235269225154872587453777232967107195017980487789378012333541538807085468666771688908797","0.894689274438589651927702113651684036779977639376580677878685570729825326911051724584870040281185572257639942077534224528"]},"7,10,2":{"tau":[0.23358371553900237,0.8719950656662158],"highPrecisionTau":["0.23358371553900237475696218145135793960682910164368969433577375004339313976980856775910528438989778586703692264956215296","0.871995065666215801940546894688479131128128857666435055504185090280417343725413296322718588614899229697471956764482721082"]},"7,12,0":{"tau":[0.06939312303640456,0.8339810453864935],"highPrecisionTau":["0.0693931230364045600674889034407840573586346685082717997037117063287188138083395985054283779394612449309487439385912529298","0.833981045386493488693475074457485802620214149098393782020784712682757481433090090209417639523155095154404288007791104464"]},"8,2,9":{"tau":[0.7685459852612637,0.9265281264715773],"highPrecisionTau":["0.768545985261263718592839723313763786609576856285145074490803667560742266757340643190883773612107302307971136639202325883","0.926528126471577338160865570141220688711098079275742843066945439200163912653444533405691779285703365522056729487743259351"]},"8,3,8":{"tau":[0.6918195617799308,0.9398787140959912],"highPrecisionTau":["0.6918195617799307114560725724476599518408483657984329371490423606952271390266283063056434270158778248531911636454132509","0.939878714095991221984095705838517546231564382544285734436255890502920724069195970948891167865484390812381256272625619318"]},"8,1,10":{"tau":[0.8500356375134203,0.9080022645480467],"highPrecisionTau":["0.850035637513420271110025392805547835027788738490374200572524737118015620290690327974784421660417469832020734071063025819","0.908002264548046698589083675134786392872178577412333108676858425569716574101400057753321535722597588713506924567908222922"]},"8,4,7":{"tau":[0.61767839243766,0.948798264354949],"highPrecisionTau":["0.617678392437659936883792879968112288633828126996012659351804834368975371072286034538555928482371140618143617052129376676","0.948798264354948929941453319016545651050727216693379890320536847088300424276708650092373592511316291748653921195370075473"]},"8,6,5":{"tau":[0.4723650190882867,0.9546040616976922],"highPrecisionTau":["0.472365019088286699188140876649972756982036989565519026520079711487761625358237786753699523283722523508529848483446972802","0.954604061697692147873971924360680699157701541045268592148148154361410844967228517545718035554306556799696258618261487716"]},"8,5,6":{"tau":[0.5448418595808,0.953651075871003],"highPrecisionTau":["0.544841859580799957889179962092721619573100098485088317952160867754190763320479004607046493466548167485714615576129919214","0.953651075871002912097773667292110038492798604432048978617609806843412469366481732008704894995720922844619856252420647122"]},"8,7,4":{"tau":[0.39939010074735437,0.9516832415576482],"highPrecisionTau":["0.399390100747354357871679063702678295684409541229627573743687518667239757217897559158821091219296468659650819427620908365","0.951683241557648117525708144278552181009279423849480589876609707505882986092829174539919126069377440396836668260759526564"]},"8,9,2":{"tau":[0.24776798499065975,0.9336449573792969],"highPrecisionTau":["0.247767984990659739259486465692125295422063635777339974449685565104268544313842171093106751503325330742661655697318807134","0.933644957379296865264412742613204190662587855938137726234957620841764158915180761887808171065929015480119699913645227522"]},"8,8,3":{"tau":[0.3249609640968331,0.9447837960719521],"highPrecisionTau":["0.324960964096833101319624685176631296802779158012865755204093885619959788993191431150218111437369326385632714400295209012","0.944783796071952009553446105169079326078267322322128104082270487898964502068031479106323828177666283455806658895179016866"]},"8,11,0":{"tau":[0.07287809565394504,0.8960869895662247],"highPrecisionTau":["0.0728780956539450408980125267239631320440469110940287950163337824466081627585898598614004449654454019174055252811771335572","0.896086989566224627102959130317502483611542151030949060417737749942973382852915116183994664645130389136522382163011564568"]},"8,10,1":{"tau":[0.16555937492344952,0.9177639900565017],"highPrecisionTau":["0.165559374923449521271577104180024862283588894200729676570428680297693804774330750175579857789822785184124811046193257842","0.91776399005650169557636472909706564042295667708394401157943032482673754213473793940374447885246395717379826445195838487"]},"9,2,8":{"tau":[0.7525043616635351,0.9915291175667229],"highPrecisionTau":["0.752504361663535106732079115940528897892356812138686021694020816910352117878951989343680845435404522626010390875949260461","0.99152911756672294427536185613588469771167250278203045636729469236853803244743800111347301622183904725297227542812921046"]},"9,0,10":{"tau":[0.9379149663300085,0.9467994016928348],"highPrecisionTau":["0.937914966330008432539075697083708289072725653662576878180949834383130833670354628997510536062592983124689743084131578647","0.946799401692834801377078387164311655226046368212181009537475565092105610108466615263598865436549323065275856697026940242"]},"9,4,6":{"tau":[0.5888246402625387,1.0121689497052944],"highPrecisionTau":["0.588824640262538614117680885353631991095798472610458104653337258496186188924930169172668297949121638143998154967232982221","1.01216894970529444134246227151022728280884374232640974272041158477668110073879716794501217060901882173260881212115402513"]},"9,1,9":{"tau":[0.8402137191687481,0.973081210337101],"highPrecisionTau":["0.840213719168748121139245773469098572248520236196931858855250731717225083015321636379630957204762940957017268953448383448","0.973081210337101017241882970232106177294370896543641020191136793953178025697353055092974576166324015653132754185118011476"]},"9,3,7":{"tau":[0.669467182553563,1.0042975276105233],"highPrecisionTau":["0.669467182553562978815588333303269730763066616363984270314458790125447369000830798838555224149751924531538849453563902063","1.00429752761052331769918951452594249034962367658353522532573893153402523004613578278046218152619957129867850817657532528"]},"9,6,4":{"tau":[0.42946012422118174,1.0144186001613988],"highPrecisionTau":["0.42946012422118172783714107685459528348345866217686542019618361715380065180468845186000335712864287952835269930413000895","1.01441860016139891948397204089690217371715856671422890419836869723026872500821267951481518324003810815765222291717154405"]},"9,5,5":{"tau":[0.5091879581715567,1.0155002268281028],"highPrecisionTau":["0.509187958171556707780933741714314218235627810344588824787560041120887131222033295842035030410951500638676486510921583292","1.01550022682810287365385486204488905480115550326057910830408038911957637913589506956597383975521653332068736757290722267"]},"9,8,2":{"tau":[0.2649939327642938,0.9986234690057386],"highPrecisionTau":["0.264993932764293755215357171655074639208164639100494223255808344331116951664927809583808846217725562838303083488897509147","0.998623469005738566296741101408471726961592326263062964568663269796526555088832977403074375682179896630901919365650118852"]},"9,7,3":{"tau":[0.348535532531593,1.0088733757283106],"highPrecisionTau":["0.348535532531592949597988774268508503694190689199365362978406144538520736873395453825592916750615187344008713694693708397","1.00887337572831069222563565133561942264255243128467433863833327810103914273554258713721806182857366922547581791068698025"]},"9,10,0":{"tau":[0.07731279446405005,0.9612770741387737],"highPrecisionTau":["0.0773127944640500483793142067984937132691648402771497344531632684465177670141008555184231042845298947832585383716559907369","0.961277074138773745263272097393418622384977607635573262428017310991609535070055427073937201706109588562764650674728540486"]},"9,9,1":{"tau":[0.17647329246774576,0.9831457855690638],"highPrecisionTau":["0.176473292467745764923747400387729453898076928740404554893904096817803750407221598731318964166677200576837651145865696961","0.983145785569063794594532098282988637963671802728755621143740633105093186635082882843868248747244972370329587850246477833"]},"10,0,9":{"tau":[0.9335359736301869,1.015284696397039],"highPrecisionTau":["0.933535973630186924728075789873607123632147595498830487820197772316454181982888425445382363889726937162160587085006704116","1.01528469639703889926455292872158225427629477078363917581030105839273794683228081658701987000850337866545672305426529887"]},"10,1,8":{"tau":[0.8282010531441095,1.0424503897654127],"highPrecisionTau":["0.828201053144109532715913535372280187724163643758435962607590053928800679115382514089233047161586688993392136877966053816","1.04245038976541278363015822985279432173097655972671414389033449154923592922593293621298506412906129951129428347024586075"]},"10,4,5":{"tau":[0.5538682077974456,1.0795044937489853],"highPrecisionTau":["0.553868207797445575186489080389136530378895094433962209747919144714136290681654781607686830091449313477035949928480176904","1.07950449374898529189467296365486529881146186493620929311617021364240981485234519647727716708478529154944762781594049856"]},"10,2,7":{"tau":[0.7330051578916865,1.0608887499708195],"highPrecisionTau":["0.733005157891686564427672613522074715009632750053239676182795390903762259787935591110663481631232645198775797050861208672","1.06088874997081949630767784307343125494558613579857336714354853793877264634613523769521669261949223245530392715020043732"]},"10,5,4":{"tau":[0.4659217777112485,1.0808468913380915],"highPrecisionTau":["0.465921777711248505834367592930830023954742367754398380796737775561442640925496984999163703915670485599564826816036887989","1.0808468913380916168697150632717022167013088434756580984301715589191177883082626606351086674687721532174674418546040871"]},"10,6,3":{"tau":[0.37716613374013025,1.0770437237511312],"highPrecisionTau":["0.377166133740130276671618750299966355075482225666215601687852452449827680441184133760531429173179480665744322724319323035","1.07704372375113131395869529341907158769267212086398006877079434770658029251174344247993551207306082726425025254903486207"]},"10,3,6":{"tau":[0.6423732456419513,1.0729640167850305],"highPrecisionTau":["0.64237324564195132364610781454388747034997453564264846172846279544563336523706019928758795168139129642106777618188447507","1.07296401678503057161628204429083518676336726374163142539857606876441736091363990113004744321696271733436845380086374494"]},"10,8,1":{"tau":[0.18989366461811716,1.052853411970201],"highPrecisionTau":["0.189893664618117153705754815243060121415883827652663175107313953404465181318294459682885505154695028803594191347399508145","1.05285341197020118153273126511932411709131436454905851569723249813766886070225595889365323843934057237116105894922727944"]},"10,7,2":{"tau":[0.2860012041958961,1.0678938904448056],"highPrecisionTau":["0.28600120419589607843135748674364377400892032633019447408488908665578316729009447454992950986293003642131484568470571745","1.06789389044480560654199444595618605340489627305484354089555939898647216848446958641406830863737734763091879535939314406"]},"10,9,0":{"tau":[0.0828776386720781,1.0306211784456294],"highPrecisionTau":["0.0828776386720781047866856321990103888880652613478975200644662882761312623553101673495441248342695835931590289880840115506","1.03062117844562930397112758507474046449573805314816578844442573682446233982477148719034628223534746650965394659333530753"]},"11,0,8":{"tau":[0.9280947991818085,1.0890312736873975],"highPrecisionTau":["0.928094799181808513617919258354466107472826122111070670704364521640668555584681025447187648385619715164327829563798878594","1.08903127368739753186417701660930791617512696882207399738993168632632999614725986527106529798407940929477465882550480121"]},"11,2,6":{"tau":[0.7090331629168354,1.135826992040367],"highPrecisionTau":["0.709033162916835393699202243992959057791747441214476840220258337476991090269418070599697024561889551697567955443227854688","1.13582699204036692807143863487321168327025875825402794542799808632076169254701625896909438120976610128055559519624399615"]},"11,5,3":{"tau":[0.4123982831044108,1.1502780821791465],"highPrecisionTau":["0.412398283104410803405703441993601547633668910802000476001993334298445056876447106690719727889892782930697712532554975278","1.15027808217914660991518230745009258559976513784973583181958013186759084204367255350563260760177599098293021900180469372"]},"11,1,7":{"tau":[0.8133928961916208,1.1173891076975633],"highPrecisionTau":["0.813392896191620815717351013968289286149038146377448508261357937578761912902623505770366827877558316890960371226205781114","1.11738910769756323011140838192689449529040704089629202067040494074298511578736576913794515939925605425028115096629022489"]},"11,4,4":{"tau":[0.5108211815976059,1.1516950858454398],"highPrecisionTau":["0.51082118159760590078056519610602398212575548707015726090156998414624962757076483151078268860123277490887734024660710184","1.15169508584543972924803976767618413483239280473830326639342718443631266260394974126030499848669989113744540265845067798"]},"11,6,2":{"tau":[0.31189512919217804,1.1426122438804676],"highPrecisionTau":["0.311895129192178048545511932952318163251382294129405447865203013471602164499367149227092435102492105884063057021652303929","1.14261224388046752506688880158114200873659605274904497931007408365464661423897708561123596083001862111610521432076184703"]},"11,3,5":{"tau":[0.6090720398580755,1.1469715021779754],"highPrecisionTau":["0.609072039858075497939797684570059015586412866858001296254996204818354803934445328602842621770167450271603281797976248298","1.1469715021779753802277838456176529199583244007292558498978140470743900435725180881398352786342035327336711708193774627"]},"11,7,1":{"tau":[0.20651519053869685,1.1281459863027559],"highPrecisionTau":["0.206515190538696851115143898284241522723190554487913585522901781187548352125878214497587040462169133770615622870232553266","1.12814598630275579389606276928181748691469265058754773871785631015143477428524894474184545402733805060060019838424451498"]},"12,0,7":{"tau":[0.9212624459976025,1.1696553667878262],"highPrecisionTau":["0.921262445997602538051211009993443350012523049528349249128315741417813344292824400084634215813617905541600367305220001299","1.16965536678782617843452920769231033803093693970472226687113891684845303618214770407649357295852334136504477206123325784"]},"12,2,5":{"tau":[0.679043213405774,1.2178951159257358],"highPrecisionTau":["0.679043213405774045311103911147317260195140599342931835609786102749102208195278003077589479047289775480249430006265964876","1.2178951159257358276925230117893011958535051439707222412487969339876177364007186636651156291470112850107657528783809269"]},"11,8,0":{"tau":[0.08986160226468287,1.105410577089289],"highPrecisionTau":["0.0898616022646828677459818371361623978368341279305173307214258000795548529096013163901669767266990437626553519346552012634","1.10541057708928901512683196477204553401580414162803562262930166619362017241821615965816448067387814663558269238059961618"]},"12,1,6":{"tau":[0.7948677253892807,1.1995439653140865],"highPrecisionTau":["0.794867725389280729803731556526948264721158285300327377786075088219763722152029596735723554313748786943363398134896081842","1.1995439653140865067717122561055218358883887045265618956471980864788096700769314625734278978667642592060808324786339608"]},"12,3,4":{"tau":[0.5673186880556069,1.2276576968824482],"highPrecisionTau":["0.567318688055606875401096319839298694281337385550248282586515507669407598147680456380773414301414653480385321491049176688","1.22765769688244825856403820153947390501595798610432374777390756041398568023520883540913216858673330565255356008562075373"]},"12,5,2":{"tau":[0.34435718049189806,1.2242186672969189],"highPrecisionTau":["0.344357180491898064289658493173694750268853248129736078669425562183873842218390831564371216292248772676981523942078353142","1.22421866729691878676415391259784025225353229074765928322573172199247451278663293177741537645685538874270680007908405491"]},"12,6,1":{"tau":[0.22740005244857686,1.210632399643944],"highPrecisionTau":["0.227400052448576853848124582580360435920086735351524560114879254569916717585875855050860672495390915638175428756754154754","1.21063239964394401901421616839080422719442507775380651388252773575008819165125945512272906772243058253636007681235136431"]},"12,4,3":{"tau":[0.4565897302741668,1.2297233344881855],"highPrecisionTau":["0.456589730274166819529185831418718440593534818077244903950404925011008245301617486321308943024558129070096192547114170782","1.2297233344881855845832242956638832438838255234847200093731651692999329177216838750482796175486916569544998982285205273"]},"13,0,6":{"tau":[0.9125218585568172,1.2593429834507448],"highPrecisionTau":["0.912521858556817256542630687368431817990277778528570621943512031650497919304303173774913567237146951147682133452858102572","1.25934298345074488192494820589817772994261989591115173447816612839245040034082800349129623070137086036182270277485149127"]},"13,1,5":{"tau":[0.7711761168437719,1.2911424760818695],"highPrecisionTau":["0.771176116843771810223995265443345976464823373164390768317169118837452435749638991273207863435560958637489664659993150417","1.29114247608186959157497585228445145324331078769868508088392538087890192052829548276711811229678009606098437386722055506"]},"13,3,3":{"tau":[0.5135031537932698,1.316673757938982],"highPrecisionTau":["0.513503153793269821393851715840595458898966588363212688066496525078994904448357106264360847617382046331974202145749546404","1.31667375793898202702738539814022436904681961991356301814651284714522401872798294321374165417858039865218394580985839086"]},"13,4,2":{"tau":[0.38604318241844954,1.314553858922319],"highPrecisionTau":["0.386043182418449551967168966701768557058282123685599641191816666282155685262791878415196560196430332409041678688646772496","1.31455385892231881512671387413274278872935919355277561640546704226612678720274496008308077377092416698652633220702317109"]},"13,2,4":{"tau":[0.6405875601495453,1.3091422912445143],"highPrecisionTau":["0.640587560149545336340445428497338570590502213529313282931437343351719977395736971243990397800826602715737959443407353024","1.30914229124451423016345719770975028902094946784222802500834701618626890512425444418072083385076209755087857778808717617"]},"13,5,1":{"tau":[0.25422222014519086,1.302462317077187],"highPrecisionTau":["0.254222220145190830822504228570651973006551011066098701331353328176521064131740704098148947377048812397396198422658070476","1.30246231707718706161245300535898191617095148387025299181944470259310993446016700595751144957057222415981081662834990882"]},"14,2,3":{"tau":[0.5895519893843593,1.412362403864866],"highPrecisionTau":["0.589551989384359298823703049690078570041980528884104087006067619337046753399407564619707930140923202152162976498638481514","1.41236240386486596092738244483265086040078684929845051283434335231810510028176348750857583517482842852821861103243977435"]},"14,1,4":{"tau":[0.7399174044460998,1.3953525962748095],"highPrecisionTau":["0.739917404446099766922792729417095399357820770938655444352278563066325885182890616373699044537499670884270952227543931795","1.39535259627480936257163747665949949177109712601948747958272777940831319744428339633686173363518622040643802287428974969"]},"14,0,5":{"tau":[0.9010244111556281,1.361224022290124],"highPrecisionTau":["0.901024411155628172384547923429583206795041181612706972855386378169470883689946511472676331655241883721546403382911531217","1.36122402229012403864739526011597486825068641139208150095293187743466781583321692178381978323571604941305406277043157572"]},"14,4,1":{"tau":[0.28975560431283015,1.4066313050513464],"highPrecisionTau":["0.289755604312830155688577231176733377939319233349023909761933153456263313602810266593327666838696404578532920437778497804","1.40663130505134643218811431606726934790650881454753039988798465854780924311024766558900187387158328313731051156418762149"]},"14,3,2":{"tau":[0.4413887079963867,1.415986974922707],"highPrecisionTau":["0.441388707996386663973842404736926793197477933949518176634483311344636944205347748190794808881761902222656149622009192225","1.41598697492270691436516412392407580259212995769620290579793258048546851206614644808416439979788802718367711707557888091"]},"15,0,4":{"tau":[0.8852815102227144,1.4800772173072188],"highPrecisionTau":["0.885281510222714370381565682628043255609809304710058682170138985717494250039084978026412945968509866388055411509697290267","1.48007721730721876286970183165759900588821412911048031416504515860240254661685173440952592712605868185717670875260721592"]},"15,1,3":{"tau":[0.6968090206553649,1.5169329342393736],"highPrecisionTau":["0.696809020655364876466766695101718558846754291219028621533356284530697720162747444848912408108542392036090642373544177862","1.51693293423937359673909325244662722068255718284330144359378187091766817653264801062833612793793647995184038387553601553"]},"15,2,2":{"tau":[0.5184086894567791,1.531431330622322],"highPrecisionTau":["0.518408689456779046020626124912823805499128541953133540825034869675736718268078426635453171391768347019197600187715510497","1.53143133062232204064125977425233811486382391253297938556699509409400037269411754145624812219503284458672874221678503607"]},"16,1,2":{"tau":[0.633351015162398,1.6634908448268861],"highPrecisionTau":["0.633351015162398040967431151680742518277669663920994426997622325206283950160301180855188891404580792826429853948742764994","1.66349084482688609341955159460178729367926129007136106249256167255845382938032040231337136941909519607894865277286992127"]},"16,2,1":{"tau":[0.4114979708147442,1.6715044582587142],"highPrecisionTau":["0.411497970814744228015474013734380976448015490074922093508156483167763253689859769736624104923814707832810016339312169967","1.67150445825871425028590589219743210759779143715236477812581643173436565565224357007789182937862631619422980015118294173"]},"16,3,0":{"tau":[0.17750825378335874,1.6500934035160153],"highPrecisionTau":["0.177508253783358747286530344676494609700627236922418470639007395045453634551385787268944502203746858143547186937876834712","1.65009340351601521735519903877371553539962924764978826509787288919560717490099415746413890391614192599188630026681424432"]},"15,3,1":{"tau":[0.3389383670962639,1.5274835809826834],"highPrecisionTau":["0.338938367096263921070369731830347575343772403468622949601990074508716244862758253276860470531155898792949911597165641833","1.52748358098268346109024639910609005631256401182440266777906616616291601727266580151271749305607839653885676126540990114"]},"17,1,1":{"tau":[0.5296793985912197,1.8479965761812465],"highPrecisionTau":["0.529679398591219756050380897059628575743876673534385427002030831703157130049004710765201425254975126341940428821891748056","1.84799657618124648607691391471356570119416503858430472101295819681392373835469428306680101046097610436617898036862797497"]},"17,0,2":{"tau":[0.826034049988733,1.8070252691453603],"highPrecisionTau":["0.826034049988732982073654640351837733901787936448956257512258362802279373297671562636076431714523359611795615764039179342","1.80702526914536040285381253885927370802033485418651617317064644312751566641109659113975257372252460762582248605588339704"]},"16,0,3":{"tau":[0.862418843179909,1.6238088850043548],"highPrecisionTau":["0.86241884317990901208965988129055521841671843593598663885713607824911857511287991834214544546380762672039327428185386748","1.62380888500435491296123786175839111199065793545706296365137439580435457642623616230303663150016084254286808324886892653"]},"17,2,0":{"tau":[0.22784042401377305,1.8364292100596837],"highPrecisionTau":["0.227840424013773048270781549955593566784122548036687036883975941752131846454645353391617879237623249832764210807758595831","1.83642921005968368286137978629899068956324902170203127005182402991647634946514042732634759542708022620832026366570437992"]},"18,0,1":{"tau":[0.7581234112693267,2.0616439195881884],"highPrecisionTau":["0.758123411269326709130139640777917109513596754515846492533738834514425293875880668696422272291738523034284595186948004717","2.06164391958818833765541464036349573558928617273387552163587034165690305136903097750915180768706622123440991124477105763"]},"19,0,0":{"tau":[0.5781034956633583,2.4808369765229408],"highPrecisionTau":["0.578103495663358302519702840691360426131954929157073236061505464128758878280051310574166236363852508344681511036472653561","2.4808369765229405949226941071458622794598506118287701599158859520658083416964264558797380398368607855782700805907806194"]},"18,1,0":{"tau":[0.3237368696672619,2.0918999010402985],"highPrecisionTau":["0.323736869667261893069427041651063972805766302776441908777285699844083794572628772376244916653809317425806741247876945312","2.09189990104029861449334249421302357841148115443164939406813324411207655779939173535614090753877083817141145062456160703"]},"0,3,17":{"tau":[0.7695877950858812,0.4650002155706491],"highPrecisionTau":["0.769587795085881182954471635193200961735716160011726373576688210959646561479086298274360033203054434959701685535953846957","0.465000215570649055857748933930805080082596549152387699703400485925993430204854814909747883744067922216889174779840133733"]},"0,1,19":{"tau":[0.8651937805521878,0.4064700574964499],"highPrecisionTau":["0.86519378055218773687572047294678579729606786231937753131054526359435078502664855833876348731246940283097027776930827691","0.406470057496449875616375076503404246164421395554841056056424538223508758306757946575412248277760276989753677447906835467"]},"0,0,20":{"tau":[0.9334447340901931,0.355288130616513],"highPrecisionTau":["0.933444734090193136344117263419998363263841805580375673771245940602716789605281192547124112675700169927390276197748591555","0.355288130616513008024169902378502213478061826945091724281189239213810334656965074826766227777834968352485955127042273748"]},"0,2,18":{"tau":[0.8140961150214879,0.44010295216451684],"highPrecisionTau":["0.814096115021487923795277356240603745728032622119013494248324590433184279971408891156906019416964277739836008801891294699","0.44010295216451684956628492272703188994885624652405863512499901533104501697122169227327976863107003010333925972439121246"]},"0,6,14":{"tau":[0.6520875592192955,0.5109906209131714],"highPrecisionTau":["0.652087559219295553817008860701413918139116242655372863704595914816267262623559476375378864924028884237679436621637102689","0.510990620913171384357917487231491072653463228275012117892127668784046221684496461384654849130944851810915807765500364624"]},"0,5,15":{"tau":[0.6896363315808847,0.4992980089264604],"highPrecisionTau":["0.689636331580884666076682345891388056920388829857262549102385758019308745133440938135825724615217186342939453430440161597","0.499298008926460406709518662626861739883518194783218877236727980310827609624837874970376582509586750096534786245531470951"]},"0,4,16":{"tau":[0.7285570645695328,0.48424123833729904],"highPrecisionTau":["0.728557064569532794086280010314383329887295598570323782049198852650493852792758725552547976871814449927472772450799206837","0.484241238337299044532846300520119243039293261226467677116714922569428243380076985394737835355652828453234832240196397634"]},"0,7,13":{"tau":[0.6154550520789281,0.5198153612549148],"highPrecisionTau":["0.615455052078928128843805557072490731574875922384204298808338679126451154186102319825097636375800744272891555633118252221","0.519815361254914797019011697291489268152558114087624713452767039557746459660847032428322548758993562370816120970757355178"]},"0,9,11":{"tau":[0.5437673306060249,0.5299881763159936],"highPrecisionTau":["0.543767330606024862181924723960584018509899958066403576621694514899920788145460954033757102824947074157258620621648071747","0.529988176315993548257017413713040793245362707579847647384373765653845986928257549608679474385014115961555171415301218882"]},"0,12,8":{"tau":[0.43715030984426073,0.5282682998773133],"highPrecisionTau":["0.437150309844260709951497159134447172501550704503761679342448645302276333161299596845085964231386873403983293410619585514","0.528268299877313283999435511640401110626377565714330875060762315040307396315125656259853564361751286980103791622315303607"]},"0,10,10":{"tau":[0.508284318831075,0.5316374389127231],"highPrecisionTau":["0.508284318831075005678039148650298645393286805247616608251951139771092171546551779090106409020657289814391406335564301581","0.53163743891272303480001001759789726678178645523007037209593031666348764138381470642615952025858409836751593793312136717"]},"0,11,9":{"tau":[0.4728031079592152,0.5310715568991816],"highPrecisionTau":["0.472803107959215199568172993374274110767580383099142172043598854181685162607545920437508906001608334447781200864820872892","0.531071556899181567360209189282948776125766966086954410093758844252024127342044299593619563835460093686720951692313889072"]},"0,14,6":{"tau":[0.36453239395475256,0.5155185216363228],"highPrecisionTau":["0.364532393954752547655909155415052481406264772823713296700865032596934872899026507329767832259653636034639515246695277357","0.515518521636322824142827449420077655073192743196466891564160563129285984571439716914533850156677782448492682417925837345"]},"0,19,1":{"tau":[0.15303295579235496,0.4207681837195187],"highPrecisionTau":["0.153032955792354955258483017668845952718764142890254149280037272595576085919635682539643614178655417345560980222236749729","0.420768183719518694617939706667216144630541659098194451500937279602625248484840922812789449626738858676318652815113382697"]},"0,18,2":{"tau":[0.20326880304270387,0.4512794394065653],"highPrecisionTau":["0.203268803042703863935371235754430592010165016251802148550609760392023585342933507564034657224862975856314102821675885972","0.451279439406565301645359456825643632508619591380736275283530103723207844797880805116190228690922201248173775184361687434"]},"0,20,0":{"tau":[0.0873046988591041,0.3756116298359268],"highPrecisionTau":["0.0873046988591040978655451168261788204615876848213588591145339675236376278639114362276293704602038558712137995914413832654","0.375611629835926795868100271389718250260377391814880120073854821136321809164004046281500566080520442341173854421386252556"]},"1,0,19":{"tau":[0.9452722160277715,0.4438126612371258],"highPrecisionTau":["0.945272216027771495086030360267074139781504840863387477270278963935554303223185265080285256023489459839489704229498226628","0.443812661237125805924457678963267843774936446717198069434184928119495638963626616539324529667297912531899987882774949698"]},"0,17,3":{"tau":[0.2473877084780122,0.47400762614071124],"highPrecisionTau":["0.247387708478012186235397640713111712616145307197611846189511663735320505910629856606824775194157050355509969595350455596","0.47400762614071125194472077757391077017201527717101830789087541007774327809478307966439724206102175430573138179916779914"]},"1,6,13":{"tau":[0.6451109406297939,0.5686973972156765],"highPrecisionTau":["0.645110940629793955505140515472636614981144289834835727748641515807969450246205927462554616987313868817491441675106759709","0.568697397215676513251488980138590434161843607761487073462146131223149852977071775278052081558611847745521012394234750391"]},"1,5,14":{"tau":[0.6860483650569396,0.5586036544279612],"highPrecisionTau":["0.686048365056939525009839683072167164861005366889008156881620519277549247064799499551405480548453052472237970636246411305","0.558603654427961218948375638280977647146378206517383328146344434346553653496843626407278735478810548239074843213302483277"]},"1,1,18":{"tau":[0.8760381162083136,0.4807305240334679],"highPrecisionTau":["0.876038116208313616125930781133760963686274502981067794246644196640000664286252986358711597942289549590263998731561260282","0.480730524033467900280359900002957176934293215229317863259021276645237542029368286927066073904135890560254112007838904987"]},"1,2,17":{"tau":[0.8213421508456398,0.5078479411217742],"highPrecisionTau":["0.821342150845639755249590898656155669263515629532579473987310772351888056709227948279963869903434779593235613026876464005","0.507847941121774250970268286245432138255781623046534101112202800391871987875691759946138974796731199962516423865616350286"]},"1,4,15":{"tau":[0.7284517843060151,0.5455174120669487],"highPrecisionTau":["0.728451784306015118639184997497005388839383513132618538708135514937989484784099869027209067185218112259569425755783626422","0.545517412066948716070782748712967878478383617210292907680268962359226878017987350548620694284955102632300096235447327372"]},"1,7,12":{"tau":[0.6051428082838736,0.5761417353639434],"highPrecisionTau":["0.605142808283873546998094572091053532005320236330319122233256090656900387520185897660151601621149535084758017622226018449","0.576141735363943365844298502248436367468442110874454230085888498526268723106434737133632459389127353818650874668833715799"]},"1,10,9":{"tau":[0.4878833757414224,0.5843516122695351],"highPrecisionTau":["0.48788337574142239733146515124284326291428200531288506766877252085424251260291793882086692751342489018937199675956200175","0.584351612269535112413836952213697332456688734384310556474691427530577346084373397870891904036908329702292246570458746373"]},"1,9,10":{"tau":[0.5267837542361745,0.5838703219457173],"highPrecisionTau":["0.526783754236174500628834280501682675445848923025196306651859392775679382058579959211216732486833475723078186139401607453","0.583870321945717340926700661533200663192407225870575733793300647643585251817209660375013676513209102830678024569281831809"]},"1,12,7":{"tau":[0.40950318906362465,0.5785822599724885],"highPrecisionTau":["0.409503189063624627642974438274090680072160391374883248964858618730354857093408909082163263867966960897782455684161266379","0.578582259972488482850039816116581181334759293227245977678009551143371974430006633742594301242425533164481767687934874515"]},"1,18,1":{"tau":[0.1389344185803295,0.4914869893676136],"highPrecisionTau":["0.138934418580329500828852757364647709976571130852014269428137798067419970097594711473861811630038181553927934692168280288","0.491486989367613581562747346813195666531567851154648844607240546952676821927489420229327103581324613098239828968906994574"]},"1,13,6":{"tau":[0.369517986239851,0.5721655837523068],"highPrecisionTau":["0.369517986239851019093804146099323706060278420478114586824069907164892211862583534765565300603195047184641234776124917821","0.572165583752306709277731177194619620360758550624227723744265025588155378068994300992803153834008831650103883387571129937"]},"1,11,8":{"tau":[0.4488678327850977,0.5826058767767863],"highPrecisionTau":["0.448867832785097662176496739108649326295945421188634291194503115234527316061526773025728852253894494856457019614386045228","0.582605876776786318603305947140484175322131752785080268214580544739361482337497865720029826637411528637981277787976443363"]},"1,17,2":{"tau":[0.19338707007554284,0.516502848276311],"highPrecisionTau":["0.193387070075542837722489032882337548254888019658423902177201410326614420013289256127852047945079498283045330160193063091","0.51650284827631095068502569341748912917694180806443063315853629431226549508286884896568863188792760735082489050619725948"]},"1,16,3":{"tau":[0.24154070126770247,0.5359655309599943],"highPrecisionTau":["0.24154070126770245497614172452767434670363141025703122916442760511442134700296473664323332626942054153906588749213038888","0.535965530959994382700817110631453500341476172347246759803134314536810683465765345445791543471949033774155422151928737425"]},"1,19,0":{"tau":[0.07028798810304421,0.45820230760935005],"highPrecisionTau":["0.0702879881030442132963253638323711167445176107613714699782231290930060321648848358160685922938817963536914718304124714212","0.458202307609350024136936243776574225386579782865170801574869250835940600797517605837904668481651978907974142392365795035"]},"2,0,18":{"tau":[0.9492481895234294,0.512843778831271],"highPrecisionTau":["0.949248189523429414756116700482980816060863745412766895550525959344660299728411256755788008558291507212854594054989365259","0.512843778831270938992360916304454852983823486013307038987609003528654125071994318974131835730058647262240911447218583616"]},"2,1,17":{"tau":[0.8796556374371811,0.5438213823492419],"highPrecisionTau":["0.879655637437181082273397186983761783827154047920397212680805011296572142504622579446772951567714226485227194181598901409","0.543821382349241952623943146122854040794780733900612340280356082049710567040434202264873380745522007001505176337151192538"]},"2,2,16":{"tau":[0.8226384597200656,0.5674750667889764],"highPrecisionTau":["0.822638459720065600824896878790330332849557610686059922209940762531859181558801733170582723720235858208945771846166422693","0.567475066788976373497675136939384240083302033027521551032492909077970369341099769485240155417163681066793630330354332349"]},"2,6,12":{"tau":[0.6346946191833366,0.6215743906132682],"highPrecisionTau":["0.634694619183336599534615520076300315704397765262948177909848070891161644360213173590973618708835227860740392942990881235","0.621574390613268259222285211663885320848353637396069185710962091215526578531860631346937855295510160307890790217755405123"]},"2,10,8":{"tau":[0.4648293572907084,0.6327386882454753],"highPrecisionTau":["0.464829357290708408530186114252902233699477031782451592784152150655633386317203919935544912174025539668606255838951108311","0.63273868824547530922036887784547561560177519281991950710571767427222324090089984524923088911017028990448927311843812386"]},"2,9,9":{"tau":[0.5070306041408378,0.6333938381986146],"highPrecisionTau":["0.507030604140837814710447673287134153751888653493390249172868544155356371766295615341390743567366513347804130210639370998","0.633393838198614646219834764327602858515595807945545881091567829911093101300264158357129577043647415852382067964729002339"]},"2,4,14":{"tau":[0.7240277655711558,0.601059374415495],"highPrecisionTau":["0.724027765571155812196299983433174817619881581959447629027290228838623351519488523164487577928093245252420748365736345306","0.601059374415494961545032959406600711640034124932211859732255677156694275842750605445972422266932320539366263183106601586"]},"2,11,7":{"tau":[0.4223496821814682,0.6297973838737829],"highPrecisionTau":["0.422349682181468184636386371439409670285694781136767279305115026257414606499212234963469702578124043054093326422505011971","0.629797383873782975556158871936499907493500448567886273194678083362186462086743988558807143213968657814720407080550758492"]},"2,8,10":{"tau":[0.5492233683917567,0.6317834915097977],"highPrecisionTau":["0.549223368391756669651237322697736691995139099387434628546601331546121210978588548608437713520417794277931286166651608476","0.631783491509797767091050996378000014161633660661920918246211898786001508345872627410856715167846295534095974758923770136"]},"2,16,2":{"tau":[0.19109070604342168,0.5752200260834549],"highPrecisionTau":["0.191090706043421694476759139121706927845865499505947659725918409391833380243298107921092220796700227285137231832654386887","0.575220026083454889482975016689847741095859985164727998602766775572621724915996164404298486242973208084654580901335974859"]},"2,12,6":{"tau":[0.37929191658381445,0.6244927583443197],"highPrecisionTau":["0.379291916583814438253421421147131821819561755543329644842836412305270773276767395176856051560343772767110495249719377281","0.624492758344319714945389118921318138955720620427903300533422636771372533826936281008900495218126913842636839894223646262"]},"2,18,0":{"tau":[0.06438018712589687,0.5256063633441977],"highPrecisionTau":["0.0643801871258968705945810777696027343448816385845158590607137956897521857020822173959951894195719565099926367960833569113","0.525606363344197759131079509541386720717942360591832943754185742737213657661747111013251152099114638372093687644823606032"]},"2,17,1":{"tau":[0.134029671414437,0.5534254893934786],"highPrecisionTau":["0.134029671414437014136291925835533264523174930110233445996698631627271928901601453263682270371984347391599035377296623331","0.553425489393478578827704317815586360489865308936614820161324724265349859444750967679867968949506009922618148327173858115"]},"3,0,17":{"tau":[0.9506934119921138,0.5741008045733689],"highPrecisionTau":["0.950693411992113804154204232358017314812587040641582746843416535090403528931253437567713218010575484178402297999741788605","0.574100804573368922155165124927850657517582528108573897651267858177143785461492681167686503884624715538061741174585345428"]},"3,10,7":{"tau":[0.4385599245920038,0.6792961711211726],"highPrecisionTau":["0.438559924592003834479433343233932761637809980592463904238547348127656774753292297362870278179787540091376389993553145767","0.679296171121172533518551960263650749860387794928286846946177431719759164456693182052465200739804898437957431382436366626"]},"3,11,6":{"tau":[0.39241422478589566,0.6750579206017595],"highPrecisionTau":["0.392414224785895640822265242738769535951580765882065427864073642650200605986617936096994328538957701638499258125242071772","0.675057920601759490772709693537419056441700606987707067575886537917199489590397735758249238553715368504103953540821770519"]},"3,8,9":{"tau":[0.5296897773123879,0.6806499153843936],"highPrecisionTau":["0.529689777312387906267517167228331882053728754583013495539855354960942715160748365623825653001738320428751184906359548529","0.680649915384393608014590252750565438359963583211391198603209231326803993365823810568453497189766868250838703542003601583"]},"3,3,14":{"tau":[0.7670720821703318,0.6406492789011242],"highPrecisionTau":["0.767072082170331889335382578990349469302850729449813504709663674987292101849416917056790981747065654048119916504174740752","0.640649278901124189800902802505897437580830749990249344118503767761666052139489977731792102252030805586338122527763831016"]},"3,12,5":{"tau":[0.34538534375008967,0.6683238234807204],"highPrecisionTau":["0.345385343750089644596234069123386379062444277582901470027722814397519164064546485996392534235605940855501810354350443163","0.668323823480720396212992802258773352803675068737478384901762994718503289222850614904947806805547928744885761041999877378"]},"3,17,0":{"tau":[0.06205993536244915,0.5862554856488755],"highPrecisionTau":["0.0620599353624491526813735444118065180101290917205039330816303978596278380829291465014555655358591190810112112502327136631","0.586255485648875508123648841602555719950056565673583173142246113633306797872373788436725222719805551757792101712104919816"]},"3,13,4":{"tau":[0.2969674720366879,0.6589157139515097],"highPrecisionTau":["0.29696747203668786624198476057058485361773147796087256823570706358848893429394714225666898244126242186000566491476985577","0.658915713951509715438195816265772465259888210566438095043321314409522171926131083519979427165984732642751295227242936581"]},"3,9,8":{"tau":[0.4842090667913944,0.6811445172680448],"highPrecisionTau":["0.48420906679139442571208116617188365728253937799681530356559644510825309480009123542073201676095846907342600908863833328","0.681144517268044789901649583665117039909581163533090253171781480728139774214057078004791775573679981295328720059141542943"]},"3,15,2":{"tau":[0.19255859357215308,0.630831630589287],"highPrecisionTau":["0.192558593572153079032551681601504503857991819230082923548685745238271738809243198855818984202789897499045932416140093995","0.630831630589287041220052211509526749707031680597258819567658727615918472731343364370688389139778686386342800031459481447"]},"4,1,15":{"tau":[0.8785136860641116,0.6577676458501202],"highPrecisionTau":["0.878513686064111659224819067253062076220522967180830752111243594367094042630367084981973248918858939922049211889360078799","0.657767645850120153351150838024637476995191342035556428831042231638265911136027204137627538077194756340365183315526546858"]},"4,0,16":{"tau":[0.9508772293837054,0.6316898838865187],"highPrecisionTau":["0.950877229383705455056076789878048215084018388648286769370833659721538926179731379462867412765944046954083400876015867159","0.631689883886518651385473292492548190053334286918246795634035082039964255364797402376256845142334283250184979874660919486"]},"4,9,7":{"tau":[0.4579699461745029,0.7282081644633693],"highPrecisionTau":["0.457969946174502855574447014002544530791412979790989236055572956833284915003427683099055847317287182544201363017230499574","0.72820816446336930505902237494044810612256754815602161435168860808565430557149533343774973623005241137504920879891310739"]},"4,7,9":{"tau":[0.5560063210698681,0.7271509334001041],"highPrecisionTau":["0.556006321069868097537939061326537864604566570184300734427109169958481003024235616184185940460236446794247745440388087081","0.72715093340010404068345951620467450479870201886045895802427244076209970033575180301640425841248966928109230087916750143"]},"4,8,8":{"tau":[0.506997349643245,0.7289033817330131],"highPrecisionTau":["0.506997349643245054492797739354582808050438075788012466240102377234029102894439199727326220124198197974479117103211731787","0.72890338173301307006149453769168632865493072453074425658953884885736281780678542375179750711864622861072964595117099953"]},"4,10,6":{"tau":[0.408559758286704,0.7250456057747774],"highPrecisionTau":["0.408559758286703992141201572864412452261149951484292420441515031754695829101897196920002819647945551823932868047860406228","0.725045605774777438509280240757956385571610047548669557605927725037268620878352016407872395127537603622744832510866459555"]},"4,14,2":{"tau":[0.19655093381658262,0.6850759510014534],"highPrecisionTau":["0.19655093381658261730680430011982783031961665252410035829322742940261821340086846405304131884655072245273572500971700195","0.685075951001453355231772828053546712579149463238234758758272534560013264690932377827440494402936993978951135137994899552"]},"4,11,5":{"tau":[0.35835691895315364,0.7193392496514462],"highPrecisionTau":["0.358356918953153632309131308740887228571893725123229158096831844899037549244873979829779595676460863791170719905682540901","0.719339249651446236750208465987634241124235496091788694500591644476091353630758330205481861223382014885950244804459301625"]},"5,0,15":{"tau":[0.9502526291355056,0.6876820280460794],"highPrecisionTau":["0.950252629135505572923977839912579175361352168559378697767222290024208995739406449720771993078107246549983721042295618981","0.687682028046079312655375603136696788200623709114756024873144443391703173178696354733218501132600082918905746834422806208"]},"5,1,14":{"tau":[0.8755168490177945,0.7127497865854502],"highPrecisionTau":["0.875516849017794519711620019249197911149369084902649753282867820925293288003345534789696042202980503598067977050730193634","0.71274978658545022714277893051491564450139790954362886990551323207568083548952551271663900905779753018461439167808929407"]},"4,12,4":{"tau":[0.3068450459512427,0.7109468300412994],"highPrecisionTau":["0.306845045951242676381134618728122658132154192854776791188650620318642511652383256526205136894411145520450122378674542553","0.710946830041299394864147599625923053075439256252303214976785534300938979406236341167266263774619456191276330289153897002"]},"4,16,0":{"tau":[0.06150509493056121,0.6436612321955112],"highPrecisionTau":["0.0615050949305612061203635103303527182934413440588197564445272305719470237506831840656100660802154887217585956500721691392","0.64366123219551119235984635351004322975074548969554247847827715663007496307072719437758616042758034070672379179657626225"]},"5,2,13":{"tau":[0.8108609030250119,0.7320434766118774],"highPrecisionTau":["0.810860903025011874554337099332733499889692932782222293795209476309154636298609795132810990011095083355482066122706731809","0.732043476611877379798502039272512159736342988752780184392278402818545580857012713852144840191391325902377677442940205231"]},"5,3,12":{"tau":[0.751268719283771,0.7472122702472584],"highPrecisionTau":["0.75126871928377105784190942339114380527657338982224191502896266588183679910318959185097468415925886026621572854562795146","0.747212270247258430462013723571302115163608550226759318379795681532595879110480823270149420290116260514101810642816936145"]},"5,8,7":{"tau":[0.4807620480579557,0.7772249028470218],"highPrecisionTau":["0.480762048057955705202052535612296158333196749273059612063224194949341308262993173196464291670043979523877863479965285443","0.777224902847021844307886651550512505634740329211531997424064261544414538762695635425194157322485091782708000017761037044"]},"5,9,6":{"tau":[0.42778889005995147,0.7752027295471767],"highPrecisionTau":["0.427788890059951455517454938202612937802333228512874201755549690369964442420068418022273353401068494262922812431984637957","0.775202729547176746961833764211121476525917848295979839340174999894381022813227461400648660152488439141005055694447366401"]},"5,7,8":{"tau":[0.5335028315173029,0.7766502546097542],"highPrecisionTau":["0.533502831517302890482698108220269503418167940517713953891179961886880672050851828918254099943587486845959173485286514438","0.776650254609754269885260826946797155913112665821640949918589681149424068843266875366557800701503584394436336369404713339"]},"5,15,0":{"tau":[0.06206489968508734,0.6997080023443177],"highPrecisionTau":["0.0620648996850873396123083500461742061851337957566660677460623333379350504559740950534776530674880275168494366015177756903","0.699708002344317786147722021376291684250411319574976330632971923628386281874844292159332225351179307563092221703779013179"]},"5,6,9":{"tau":[0.5864267974773666,0.7734697017975773],"highPrecisionTau":["0.586426797477366581487820735793850529847160161582991411233537271224890978990792969741414829814424177338962379260446442001","0.773469701797577311774552435137894149561655328997394416437929731674255982932850464728727574612664942153456988746997648604"]},"5,10,5":{"tau":[0.37413988654126384,0.7705338449240868],"highPrecisionTau":["0.374139886541263835281689381196825547897027347760799626432019076165701134260462011251458059043535119198485991047153710827","0.770533844924086855443909684947648457278196434276908008594561003872219088571598022432371938024778934838804439717250378392"]},"5,12,3":{"tau":[0.26247488129638236,0.7527173322726874],"highPrecisionTau":["0.262474881296382382707239440486319710152281364019492714392729767028475764558921314695812048675749144727083228114987330718","0.752717332272687419335193667341637701057130507202357738207779333269378517171079360411367337762939734605399659044014990998"]},"5,11,4":{"tau":[0.31927917250676763,0.7631036013636668],"highPrecisionTau":["0.319279172506767603991133624120772621452351869457712563735075432354301610231543114270448207339959697399673307766627295871","0.763103601363666817665387575006447879498852191695670093365672770180571885466374090415255465866818012157061880530748471949"]},"6,1,13":{"tau":[0.8712616726591564,0.7679461227250056],"highPrecisionTau":["0.871261672659156435996657101509542020352921483210605553774700967320877837555958809359801518813063575436673572591920795948","0.767946122725005594191523273389573306473321766049826470203726418446254821199145731005512042894260495344203434181694292353"]},"6,2,12":{"tau":[0.8033014096247176,0.7866520747196156],"highPrecisionTau":["0.803301409624717625991019125140735875675038124059321029379673441641874848031503006729162951752920319993203924243537473294","0.786652074719615626239701181212742787310559618875731226057226767459771355481552075403114068280956522578780232037454462253"]},"6,0,14":{"tau":[0.949005421967283,0.743380480631949],"highPrecisionTau":["0.949005421967282957319929647537135286771544187076672459441442232406276241178162693126056367699709890530190004226249125538","0.743380480631948927396468819157944643161460940732910563821039572716986184380107659951272260618430097162906501516684986992"]},"6,6,8":{"tau":[0.5642770829010021,0.824811927817701],"highPrecisionTau":["0.564277082901002128954658113873089761254256398707130758445127697128863084798833914222873531611487926842572599969235728855","0.824811927817700945515451711648661380923005017813115136034613809831896362776483796897440986507925672086722912146317467892"]},"6,5,9":{"tau":[0.6216613889126608,0.8199501160630182],"highPrecisionTau":["0.621661388912660775556793864625230391665624773146914295403350557163463199477235466639076061083693759922460640337298981839","0.81995011606301814658498762345456424427740882286447479458412164058902813853081336745207275342161731603264144983543646598"]},"6,8,6":{"tau":[0.45039806981034697,0.8260876430992105],"highPrecisionTau":["0.450398069810346970118993538318813012820758356536283597513969515303194738398030473164729825181777929533953094486548421226","0.826087643099210449628609259193109735905134481159517756540650467475803772868406051374021166277370950354789629343869549527"]},"6,9,5":{"tau":[0.3929168078847036,0.8225247498321271],"highPrecisionTau":["0.392916807884703599846373522093878676405145363405427109044446015583228262149916496862010558566335967650834665216349030886","0.822524749832127085766327861471596680319613483946486789823389150697407302249406788562159404299439433078159264531924563316"]},"6,10,4":{"tau":[0.3343412477072344,0.8160640810544083],"highPrecisionTau":["0.334341247707234363476809509612737819928454032499775193993781713873224628912117761291318158998341062468452907792731494213","0.816064081054408328335784800433203296718083178995129489876693705591495123408392155013626893131420052272674374017104614654"]},"6,7,7":{"tau":[0.5073535058539256,0.8268428425963388],"highPrecisionTau":["0.507353505853925549023905570226357231428478559977197764806880342768381646271813247692294018196191948296105246687474449423","0.826842842596338787498701675427815974060000618058375067500819537530533899627649868692667984921353862292756892635950460289"]},"6,11,3":{"tau":[0.27392193517519225,0.8065343529110806],"highPrecisionTau":["0.273921935175192271387241261876992694962742660532031519244861581318593635949875475612721731377507621694996451762333834181","0.80653435291108051593563509905381365018831471181219479683747665961989736405491519663753528616383178389312034786301741961"]},"7,1,12":{"tau":[0.8657818824511581,0.8241954415213485],"highPrecisionTau":["0.865781882451158018962759619897481603288935007221325108038022377007044504841675453953507949776679163344304739599113117617","0.82419544152134846458079613141607303040618947152155995104342238080864642560687013556822079185495903178899429145303162277"]},"6,12,2":{"tau":[0.21051979622253342,0.7936538002642511],"highPrecisionTau":["0.210519796222533412488610952241580014246025028280600315499481806462969419740640699354346385437672939150421348595450302227","0.793653800264251138247425125052390947780268746896429969930733659938620342160005144893099010312378389735511370241553486913"]},"7,4,9":{"tau":[0.6627429229088004,0.8667655477589068],"highPrecisionTau":["0.662742922908800480355206669704992553675540034553594618326728931833778113068608868126607599911257297272972041599389139943","0.866765547758906791099903843557111187042472479008370141426673866079311048842682027458864655169972551653681315860483060621"]},"7,0,13":{"tau":[0.9472047646577716,0.7997573203196441],"highPrecisionTau":["0.947204764657771606191346704736587740974964062862753690102661686952784219680709910711019335779527347782665273171321326095","0.799757320319644081555755579163427346718562865659062914520994247593563749144556077918326121979350184782579741413057879987"]},"7,6,7":{"tau":[0.5383944826752685,0.8774484344496287],"highPrecisionTau":["0.53839448267526847985783435203052117648944450181261893886537459030760706958814555342966172297160212102725398421334524969","0.877448434449628703034033673597578543424044051810255785373093217058138645284281177694737166616919156521170779250365465074"]},"7,5,8":{"tau":[0.6001384523605778,0.8736810538919624],"highPrecisionTau":["0.600138452360577787771050663210734613110266692205959974698290276399124045903346610000297120327826838833405830696836680385","0.873681053891962353856391297139005552203157079494631240477728210964743469381800197664891265533006241044422907515686446894"]},"7,8,5":{"tau":[0.41506844246528035,0.8758427628238259],"highPrecisionTau":["0.415068442465280345291974750494846539447041691896204125640330398479685721743674675910116522718708028733484772122339818225","0.87584276282382593442789559988219586724457084755463895664545984594065013675564521792287403908663656774063627833293185003"]},"7,9,4":{"tau":[0.35229382468386145,0.8704204959607471],"highPrecisionTau":["0.352293824683861456788170162446417178728239638320686651747383339776375806386673390396419315945376330306947385075230477595","0.870420495960747146061384766502518668145550869280987425129026931522730346225573209230590208616347665543819564897013515564"]},"7,7,6":{"tau":[0.4768952081098166,0.8781641919980808],"highPrecisionTau":["0.476895208109816591814487955474221071651438018676922338896045789612098675150775888457803472402047685782220530461661930569","0.878164191998080809577512411898073563788320864713764183163804741127333225178128499617324271206987020065582035783821903111"]},"7,10,3":{"tau":[0.28778747724203524,0.861745927463801],"highPrecisionTau":["0.287787477242035214232613755073401636166491918645297460977280404776715743313155319475525574107771309153522346842016827374","0.861745927463800982190231052717716686481197802076924137263592033606837102338768034613972148307668611813953212615964636317"]},"7,11,2":{"tau":[0.22038996340972836,0.849550941370944],"highPrecisionTau":["0.22038996340972836387183248400180271185780754012280885413008000932267129848243642786068840778166631444953445798968579134","0.849550941370943990155287096665695199966615551807468992822832917350325289654966370871970908879033190792985309421475805454"]},"8,4,8":{"tau":[0.6422605439149802,0.923422579121353],"highPrecisionTau":["0.642260543914980224197766538573881418988983077709071522768587514612945503238818363991846996307085621199200075065269525616","0.923422579121352924687145743627144755243811722372866382452979867909914950222328013970299670079030868043069628551223487403"]},"8,0,12":{"tau":[0.9448544163374902,0.8576486433653092],"highPrecisionTau":["0.944854416337490277007656824471243116146911246934059707317958266742412124215938065520874591375582550271477421452054290205","0.857648643365309228327832203590366468263622639024907416248083719272298201893709046705873785278179433486183031605522301487"]},"8,3,9":{"tau":[0.7111739631884995,0.9138942617938148],"highPrecisionTau":["0.711173963188499406277454649698466542964571928881496740576174329755311527771220502427396282533731010460851269505107997983","0.913894261793814811328059474568395471806407540820231383685788737560278921373875882160287976435262610213629153179340050295"]},"8,5,7":{"tau":[0.5748225326778329,0.9293437148799205],"highPrecisionTau":["0.574822532677832950618377021024887748697943336014821221206996098550733440089596618229205715017128056438405467197861873427","0.929343714879920473291779443538909412035958977290305859545731521394183544461382833037000015977375876850172866478744058643"]},"8,6,6":{"tau":[0.5080305111623127,0.9318408099809753],"highPrecisionTau":["0.508030511162312734721635089964686994001210935130485858586814919560113179634316675275820864310986230327406399352204980948","0.931840809980975320155317482694054871516039104998725904545427546352233809745183415491016709771951402513500894756176287949"]},"8,7,5":{"tau":[0.4411852972338038,0.9309821829100619],"highPrecisionTau":["0.44118529723380377397853414418286221872956900136102628094683838484617924699466107993867047080211694778562398127278313679","0.930982182910061898401212519855676232200781251710086288952387243425571958493238215297429824667754570823407026180076974565"]},"8,8,4":{"tau":[0.3735843050723766,0.9267389856600805],"highPrecisionTau":["0.373584305072376593206241409841232854863984876287853341968923426720338329661085698470891030173626438566578287607412113649","0.92673898566008053915918188540399899117971515649903660589136225559261041569488339971302878516582389337771203401612094512"]},"8,10,2":{"tau":[0.23238883193294246,0.9074491356117611],"highPrecisionTau":["0.232388831932942471758096199313892224946060323281404822057125131168432478924913754646281874586728626121986309438747346012","0.907449135611761065052627514171285107204569784216439389099865350492579408365507360115987623139468226190345169003809032544"]},"8,9,3":{"tau":[0.3043864226863108,0.9189811290611414],"highPrecisionTau":["0.304386422686310794939072231446558170813408144374482726880986483285927986409037192052152584548805705155568466304593558784","0.918981129061141426393326188856651860114694925140273457230567977754359902625972922434169303838972125423062620088867012048"]},"8,11,1":{"tau":[0.15547889471361806,0.8916759058920652],"highPrecisionTau":["0.155478894713618046319270348626588419829525076441777185787464443383620411829513158572165862433372758920834876543291463837","0.891675905892065216391689104678993148093296444300396034021973965883021218882188303336320674322145789933179527073705261987"]},"8,12,0":{"tau":[0.06850102928746579,0.870716290816529],"highPrecisionTau":["0.0685010292874657921678597200519917278166295555495849917380630513389562589294564902321910767176886357980988874189920121864","0.87071629081652899963851394937131492578778553696685102009487437930087635534259676988206359192079952602775661420158735023"]},"9,1,10":{"tau":[0.8507576390090396,0.9429294207282738],"highPrecisionTau":["0.850757639009039525704638741816394622396046914221288122512940365130564413056099772413117074152645719839953170241255632793","0.942929420728273839048797560125016438699468693727303175528878243162887814907800119735170199797221491063449715550742334009"]},"9,2,9":{"tau":[0.769239041396103,0.9609959403940375],"highPrecisionTau":["0.769239041396103009838345529282748176398482069101566662261953572957842330133433682178751566263215560239811746100926640821","0.96099594039403743914871455989929957149985261156093431766495622228342231583146540006735175052517295308338422104869659516"]},"9,4,7":{"tau":[0.6179781260426472,0.9827315055284488],"highPrecisionTau":["0.617978126042647215117505058119769660370754167360598569433518253305365134917148245442411415036039296550601431095565149893","0.98273150552844875113455960432437770322054148098778478163577556766020213799536781564440785930317399118097173623077854147"]},"9,3,8":{"tau":[0.692350700489226,0.9740214971007513],"highPrecisionTau":["0.69235070048922600606869967542702012323446527485604002338654675245847136488657073514634362482046506092636070880546267201","0.974021497100751235044638421435270511243358019597629041479079093241270981798489687360380792652032385755926873008930761815"]},"9,5,6":{"tau":[0.5448754672556565,0.9874790447088787],"highPrecisionTau":["0.544875467255656494966905572409530925061040430564801092245378825496331739975152651866754054325565948040550736715006308614","0.987479044708878757808423324777903491856406794192561255634848561793240376062763138030424671346620566724150912328507331087"]},"9,6,5":{"tau":[0.4721213945958076,0.9884256526872252],"highPrecisionTau":["0.47212139459580759528904377220181342420792089051387785784117205198216594835445179061499014322859789916644042737692506","0.988425652687225133881169607080749865986737826622345505054003674316290333637788635243146802771863370670914805149882538719"]},"9,8,3":{"tau":[0.3242145598234829,0.978893168791775],"highPrecisionTau":["0.324214559823482889410841873964113493346758785896752006689445542427480695293371039581496988889724060643373416019235659926","0.978893168791775014899385557932647905747255000181274714219915041603830385305187921765402487767456832956610201806715490036"]},"9,10,1":{"tau":[0.16460738160678362,0.9526197969689156],"highPrecisionTau":["0.164607381606783624261251837589600870447834521516785298896722959554358917694057849831820642229712158246921474041062147918","0.952619796968915618466043350265237837160655606417253948664402213124990733538252305988140475265629246562703082725383392931"]},"9,9,2":{"tau":[0.24685382642683565,0.968063153975901],"highPrecisionTau":["0.246853826426835654022004590262195576818231426347022734866620615035444255942587852734247960953139685840610891986333718034","0.968063153975901063599825715394411685629619088641489445950436388943805377477039114531041045958348300428220930004996173115"]},"9,7,4":{"tau":[0.39887866010598255,0.9855970705942478],"highPrecisionTau":["0.398878660105982545203414428745646834656732268989437199588854034721828902362274182339638416719122353536555134872258640474","0.985597070594247826161099485987638596257800002276017912241866068731524691068756646183478726646633746555149897510580978431"]},"9,11,0":{"tau":[0.07216993029329806,0.9315270763035323],"highPrecisionTau":["0.0721699302932980543703095252891486491466527367214262064031364225717732350833465674931918036665401675110382627792134051512","0.931527076303532218069527755575298187701529651809474842175513769227918125366071255460591628003448126103432945062622605331"]},"10,0,10":{"tau":[0.9382824726627461,0.9812712978977031],"highPrecisionTau":["0.938282472662746098226811355908713225200006878466643555874157655680030200181212629132541580239715673898883746397255738095","0.981271297897703088041510308773497530264306898352412630999105934795361931620136996602357447752437336107912890587730761722"]},"10,1,9":{"tau":[0.8407689556224499,1.0070461846339962],"highPrecisionTau":["0.840768955622449936363255332504337643228205912292185103387638327374092765119723095697431633559152116516380085548107867204","1.00704618463399627651326695754297455111950230577411947615481620558291605159419225034326052134949935889322219684537046902"]},"10,3,7":{"tau":[0.669821925756491,1.0376462119677081],"highPrecisionTau":["0.66982192575649096641908814359332427714656839857449441589763447822002768263354388178615810478929049825127351635801279663","1.03764621196770819888963781984828565838567247747942935284321254782512439458855167558598412637783648525729059510314035888"]},"10,2,8":{"tau":[0.7530137795148116,1.0251278916492816],"highPrecisionTau":["0.753013779514811651800290629180447891361583325538272109801375247045431186879422783358575704189457203711904665252335613748","1.0251278916492815971815739366165017321541817529952496402905573848590966637310067378691965465040430068849169557209550785"]},"10,4,6":{"tau":[0.588969456709127,1.0453703353205313],"highPrecisionTau":["0.588969456709127092379602772112473673425772454175130042424119929593108407961176856939976602247617623412045318995763138055","1.04537033532053128996093318728595552295943469265746161962222465168246490541763837513251376710777390616122608159656939655"]},"10,5,5":{"tau":[0.5090987607352357,1.0486497473200374],"highPrecisionTau":["0.509098760735235781365058995488925970164447164748797063950569151379368835931644226272425253430026728991459154179227928984","1.04864974732003744992881083834147688545100585110928934408960885805398880879818144083095816139069965764618934667516086976"]},"10,6,4":{"tau":[0.4291362131994558,1.0476091531083018],"highPrecisionTau":["0.429136213199455797375586093171953068939441774880061293518534002589688724463580530372441667266080649308214853691823081871","1.04760915310830182541734963828056347349673926274926076386440474291798313202064552997067470366094157849734331581986881288"]},"10,8,2":{"tau":[0.26429935912238806,1.0321867139817078],"highPrecisionTau":["0.264299359122388059291841053433123713924575871243633674756358477809584671705266508960844868985654128129601253601398747624","1.03218671398170769807432738469943902319586128630923414793038852056956084814788346183656805621151603456944817597950010611"]},"10,7,3":{"tau":[0.3479995222837529,1.0421997442391353],"highPrecisionTau":["0.347999522283752921497578710718171104108306794687765310323702019556674194677770720080610050754569797807623669431146273265","1.04219974423913533410110166400440090786235941120594927051278516665137258967398494539692000900939048540159271378698982358"]},"10,10,0":{"tau":[0.07674529327873879,0.9956606778220065],"highPrecisionTau":["0.0767452932787387919639518792844186559779358415860221128349907980668401401567606979547510442482382743906183261198592307872","0.995660677822006469140838470017840800890286623276470783860196527345263252943677401348938728616248929854014575481693533491"]},"11,0,9":{"tau":[0.9338267288350928,1.0488665172849838],"highPrecisionTau":["0.93382672883509277154308765164615167403036182283114711841797145540487636876463726991437598243381385903210608383240012958","1.04886651728498381262811266788372971559428536072407427919129276506920730164523299497971619662361489182287403876464674047"]},"10,9,1":{"tau":[0.17572680625854176,1.0170574907186714],"highPrecisionTau":["0.175726806258541761685363282514141424793822263009074166650574495906560465381901704202970526597872408465096590983745514911","1.01705749071867143775923842915718973468292220633268448869785216213654233782457176725092686743937935640235192181638652451"]},"11,3,6":{"tau":[0.6425887131787982,1.1056560377341023],"highPrecisionTau":["0.642588713178798154404480642008931927981887664164696049152195720177490475824596496676568042839302843176287482274109930031","1.10565603773410230644212916071535404566612536178715351150210153242862004792570457739100306742517780305141904939546858238"]},"11,4,5":{"tau":[0.5538936135952062,1.112102725831092],"highPrecisionTau":["0.553893613595206123986299183683554390789007343823272412575718411267596223756346988670564696796807774962629418352232585447","1.11210272583109206143051704339444084988119237336288341727129660948228406080875125212923646698974973442796524131563455978"]},"11,1,8":{"tau":[0.8286241390557273,1.0756202783549904],"highPrecisionTau":["0.828624139055727308911830820968358549488697139017173294554916894435185500732996628866135998772312640048440718039437734479","1.07562027835499038691222603237490525431859805488182879405997849727379893940116349542170543767294574207252208299954520924"]},"11,6,3":{"tau":[0.37679726813510905,1.1097215573123025],"highPrecisionTau":["0.376797268135109061967451095108743469982925413418118022109206180262432779432798653934738791412069341927756790494923857934","1.10972155731230254224656790015283747684854831389341822801977522689882637966255051719276532856002168463155266802900057943"]},"11,2,7":{"tau":[0.7333686313470154,1.0937682059128995],"highPrecisionTau":["0.73336863134701533722171141515723245164271268730710298650956430504635680666413546609240818359592149720519648502495135782","1.09376820591289941626512137207040075004773136398020613415007595304790266844604121349863783766990993585710987727907009317"]},"11,7,2":{"tau":[0.2854818292832064,1.1007482551080328],"highPrecisionTau":["0.285481829283206418907080908867297601814503505275228250319195338808925078329527255483435241756196529994401287277638698751","1.10074825510803285192646721880556528464871946563259929368842331026490803607394032250673112454725182321111250543676834106"]},"11,5,4":{"tau":[0.4657441588662594,1.1134405118212167],"highPrecisionTau":["0.465744158866259430058577194315976080290161360140317639743763148659419079738145704827828468941658126703952384839958605381","1.11344051182121676559029854406432207823358180534391769381709388091165228447575930522331021007008387247258346603709909906"]},"11,8,1":{"tau":[0.18931054885956924,1.0859835971333383],"highPrecisionTau":["0.18931054885956925325216765048117019035073148565653136962000922695288326350336668951976385124839350698586051994322968596","1.08598359713333823833329226614507264396428480823515974774263261398708033670106265880836628848119226047275830053626300726"]},"11,9,0":{"tau":[0.08242091973327052,1.0641343197719606],"highPrecisionTau":["0.082420919733270516242807173332332479728477188024465180732725899009732077831455049778759848707210301678265309542193488701","1.06413431977196049163057707891665982000786141524449524530196888211939204463096209225775297790214339451663802211283000952"]},"12,1,7":{"tau":[0.8137092948956235,1.1498999113183468],"highPrecisionTau":["0.813709294895623460456144661229506338451389136784527634710883780388642497626149687419237622084935632325182375637091411333","1.14989991131834684433644735605269350122103840257686840290422062142263839142492405313813041765600923015455551494167490752"]},"12,2,6":{"tau":[0.7092791764383864,1.1681107244100928],"highPrecisionTau":["0.709279176438386501986758739747612979773137899857538881732685609655504220149636932490756908683406179326071907926801941916","1.16811072441009279785003954187937169798395023237691051689022567997097720209117018985743767472610096274807242714387562793"]},"12,0,8":{"tau":[0.9283238818004315,1.121876597918606],"highPrecisionTau":["0.928323881800431422508372087559565545609198811214919287667906611304382405686251825632597183285043600517922050757349513919","1.12187659791860614340150141828946043644688513369454279455395097256119051615209211291509199418106161078001744248698257295"]},"12,4,4":{"tau":[0.5107564248866314,1.1837987337301703],"highPrecisionTau":["0.510756424886631411830364302053120780904216517169059542736189339236220725675283141954715792183436568249195122010354729941","1.18379873373017019538084041537729537724363675041402802926552017111656553087056762440774298110006709840314627503297425442"]},"12,7,1":{"tau":[0.20606456241038026,1.1606274926557114],"highPrecisionTau":["0.206064562410380270780406970993564281705393348578043363809433936928970496326327366484320531574545156954953534325423208167","1.16062749265571132504105328129557016937697071715941337922120845696420855127979667970613799110264787305251270725482136699"]},"12,6,2":{"tau":[0.3115176079163752,1.1748787657906348],"highPrecisionTau":["0.311517607916375236481599374218981019391023416091650168942227196101270919842574048168525510287676358823175358208616957566","1.17487876579063466222998542558291860242117147004006408229204137518825408625702242627590987401403343119726601835076051657"]},"12,3,5":{"tau":[0.6091772636217091,1.1791216417243613],"highPrecisionTau":["0.609177263621709142114341646109838995653961486317758827836424362069369014878895332731354282286563465427002608145988738416","1.17912164172436137284704981120527539289528008948639054849885929154611859504837869181191308999182873477466456788793059236"]},"12,5,3":{"tau":[0.41216306158035293,1.1824200930480457],"highPrecisionTau":["0.412163061580352920818381680307454614309467663168171535976354783665076609767303208010829647751241878868228196137723017478","1.18242009304804581609793588831212039336384702309277662394707818202780782423868109879631369274094194968635621514688816659"]},"13,0,7":{"tau":[0.9214409778064475,1.2018903409946313],"highPrecisionTau":["0.921440977806447505786957913612722395136540750531318918638924480194665922767909711328363005231232696205492637455114383479","1.20189034099463134122761679029567422725749299598749547505589556937540529359793954128474489200713539927073174872587216454"]},"13,3,4":{"tau":[0.5673376119980529,1.2593642077449851],"highPrecisionTau":["0.567337611998052889832807634304177920410492773578748778302540698550878624162212123541068546208534328148045420282332758633","1.25936420774498504433000137767363311718885627661819348596123898507295866596802092010075539480481231737387516540351942919"]},"13,4,3":{"tau":[0.4564608964663985,1.261426302064711],"highPrecisionTau":["0.45646089646639850416561246045870054075908550391493211420113290225939575209033355804000462154662086838811878576389156324","1.26142630206471085702467085194729640179574293972997702726434842279087871768483677593785999468825168222199945408084674809"]},"14,0,6":{"tau":[0.9126582707654184,1.2910737765372053],"highPrecisionTau":["0.912658270765418383567315446632638410882715900090339722344881508767367498536942361367590964598535079927821976720433669409","1.29107377653720543818938162625975554819491589509014078612895741474903499017580497453929033384983045655790841878953227391"]},"14,1,5":{"tau":[0.7713329980099031,1.3226591389833764],"highPrecisionTau":["0.7713329980099031242381880429897185973868760547243520594159529507668439260995434084339385017580143501880853890361126587","1.32265913898337636064882511952264443068966495119562068994397873381917822768075509068074508753906191085422696226892973227"]},"14,4,2":{"tau":[0.3858763162168333,1.3459393003842848],"highPrecisionTau":["0.385876316216833347885713243629348679640391888937662951680743511406961658861866462901043615274768255227973748021350355637","1.34593930038428485755520627924611512271075534432116533043231534873792508637480971842489634167163285441968408480935589038"]},"13,2,5":{"tau":[0.6791941374830006,1.2496876902740068],"highPrecisionTau":["0.679194137483000572932303360304212980930120629536604904400895793954979320983787865071011222591546550225492660428323455203","1.24968769027400693685351910199896664584970499368645859024671731542893452993988643422977681447685637562787205705471380152"]},"15,1,4":{"tau":[0.7400145642108028,1.4265054288733225],"highPrecisionTau":["0.740014564210802848019054677318761224745997149571145566335207751154204510355345448282287633992945724455333185092951377517","1.42650542887332263534627742335477136094582419095760962003141094138944316787879720560873399874864645448998666174156566861"]},"14,2,4":{"tau":[0.6406618187043187,1.3405342473847826],"highPrecisionTau":["0.640661818704318656346475743780613473119992281569584754784883848988031569831708813804331476325733778393936166633520002663","1.34053424738478247000635415272087463780493816856150856823927974715511191118572578503813490383615358653025601933532605421"]},"15,2,3":{"tau":[0.589565782991508,1.443433921820509],"highPrecisionTau":["0.589565782991508015466284779092673731491608482761295162634111466652544217144463592273018867109363399143473599592393471327","1.44343392182050895280744843450189571709715855887203220223763733386555185124864677990086711466456880523340142415317771433"]},"16,2,2":{"tau":[0.518377718478451,1.5622550296575948],"highPrecisionTau":["0.518377718478451001480532468238572807188696700171180685030465539926224348181565030912712778044007985242314110040885089191","1.5622550296575948072937774607292751939633112553901536475671508580890668622532966620032353780564912813105494126510461717"]},"16,1,3":{"tau":[0.6968572660640867,1.5477976212469635],"highPrecisionTau":["0.69685726606408672528074867992939857906476430525714758422989360306213556182060422064434163286055771393438284389887567425","1.5477976212469634959894597062433756920761041513743838611994102423389528935850205080036813375338666835104954149805259878"]},"14,6,0":{"tau":[0.10992161717088568,1.3102333651357438],"highPrecisionTau":["0.10992161717088568002311143480211726532347479537467961023906471111119094480192301459714303309106379224641257165696738393","1.31023336513574391296567244084430173600479257445293577071952020387130798516562299461818831205391110831127899350975096482"]},"17,0,3":{"tau":[0.8624637210080373,1.6545345842957648],"highPrecisionTau":["0.862463721008037300729546180917544533623569334404104598160963002797706698152530258332252806072466388317007031645879499121","1.65453458429576480434531262582635719921080417178398946759227933223129069381967449917808605327321681929977222672567024457"]},"17,1,2":{"tau":[0.6333606063316902,1.6941365026506623],"highPrecisionTau":["0.63336060633169025757569847567151004979216132056182783584396561879015567067065418765182642334508266030396582001416387933","1.69413650265066233131264323166584455542705382829968378268172859854343566702413264969289119249148232669447126884807232756"]},"17,2,1":{"tau":[0.4114399611648003,1.702147629870708],"highPrecisionTau":["0.411439961164800278133882465000195197867349985821783244995732243876956323122278183337103845392921231628777657866533285844","1.70214762987070785698745985550662766335048576055106814468560168146327102213140136722559455471527435808609242456331978669"]},"18,0,2":{"tau":[0.8260572802927596,1.8375574481974517],"highPrecisionTau":["0.826057280292759561516628834480847473130203378688125252466134679714548938438077012019455901993767791019009825284226558755","1.8375574481974516782825508453375132044348241330610268278251845600702685668519058533937892484083200065818092071971772367"]},"17,3,0":{"tau":[0.17741496118996977,1.680808745272927],"highPrecisionTau":["0.177414961189969768829472751486429148221464066768125505984893851767966348172535104590743821239403478288344605374996881786","1.68080874527292691868317454052655839349353594464112589253324391943255080283734878859722882868182011055055403175790025669"]},"19,1,0":{"tau":[0.32370774103566846,2.1222987823539605],"highPrecisionTau":["0.323707741035668439888088189963959233386610680702834664082160191850142112287441693123598069891912932541624712380240232866","2.12229878235396036688412680611069144266738466652595362358775220027799151421644750652134195244994292102857498739685686453"]},"19,0,1":{"tau":[0.758129275884685,2.0920451203126573],"highPrecisionTau":["0.758129275884685008994757147241346387954915645006152112249043612362316335681485521600827159772862510794415418507120534751","2.09204512031265749655861871126200517376312879547317935889698423647676649758257759438616397582517123113172115198387946202"]},"18,1,1":{"tau":[0.5296615882687118,1.8784878267126017],"highPrecisionTau":["0.529661588268711794788257272065232451211932644265269637276589605032699139410737677228712469095907609616586676129459274911","1.87848782671260167633066329150276548978995686840173376591410491210538039398747793807442051671886048033162852808919002265"]},"18,2,0":{"tau":[0.22778172308595854,1.8669555757837464],"highPrecisionTau":["0.227781723085958531464580792448334721023059017698660954716797198077405725447023078616933069649724871988986978653857654376","1.86695557578374631500877091632256807950979137056153007715937829552866763626430925733921117970711473426727896492775414784"]},"5,0,12":{"tau":[0.9425170062325723,0.7376322020329554],"highPrecisionTau":["0.94251700623257230730558229401512682549060882764896776675076153371784187003981077770396181585100833504384728551315630954","0.737632202032955478133078307228366197692780738307183719424815468614343173046656344676439230280572888500121366249713577551"]},"0,4,14":{"tau":[0.7111591813591145,0.4945175309440681],"highPrecisionTau":["0.711159181359114531769621944916684512419155231239133833706708760070056388695108176034205407154588412080741976336500023793","0.494517530944068077751205526185091509117180317366432074910073283085815368349155666217207923050238837445909022525253168635"]},"20,0,0":{"tau":[0.578097755547965,2.511167757122652],"highPrecisionTau":["0.578097755547965014584610877808462613780712641651030370972062065969221997194864113236035426983884144426801494676900265736","2.51116775712265196126516652500075251463224823739222037524797030032290295301877428441380364431323190756003403051679662331"]},"0,5,13":{"tau":[0.6688926875582869,0.5090939586209665],"highPrecisionTau":["0.668892687558286924916771927393629941079859861388924214539326483528616156931611178817273126457189028999993420709029059733","0.509093958620966540370356447413931920908245867358224389788644999589615088731384153574673313806467868865420515465637135754"]},"6,0,11":{"tau":[0.9401264217988543,0.8030293356828345],"highPrecisionTau":["0.940126421798854272681537336124424411507403541827924310921189873640961226396737073754404273253611447052969593060792959243","0.803029335682834523547168765245735510799954997143249126131473884409241416321907043240116131784460101609531477979589216791"]},"0,6,12":{"tau":[0.6279648943922558,0.5199669649518656],"highPrecisionTau":["0.627964894392255794034003632219903867440768631121677760045671279998701564911362534045128983590746809944243155481550887668","0.519966964951865597796111513974558779181813512173235755888502153169253004400274427007402178594937450134986086681031212703"]},"1,6,11":{"tau":[0.6173409615179278,0.5817298975375367],"highPrecisionTau":["0.617340961517927832237425243276274953267607421325857045466590469113921638459265962339110451796966093384452480669288230832","0.581729897537536660690962068085189237307564406304061307870205303948984124126376305269687258105407870215064940512582344347"]},"1,4,13":{"tau":[0.7086045873536139,0.5600732089621617],"highPrecisionTau":["0.708604587353613857078894966926042148232033973035808237305714365961236652863174076772842630578032419133712370690500100889","0.560073208962161702132333539502703517479030099232108143022787189731096722252294400884598807056718702234161795065711997019"]},"1,5,12":{"tau":[0.6622671584288253,0.5725802331817536],"highPrecisionTau":["0.662267158428825377665517564454229602643241699503970962265719562864191409609858808730753114370028481280261274200323817033","0.57258023318175357721084793426958858438787565842403000114831879034592230487423004765272363130412207042983679544761590907"]},"2,4,12":{"tau":[0.7014743886119195,0.6198584826330498],"highPrecisionTau":["0.701474388611919410934100526977958227090950308474014588982525485942187784044958713052870219451368219796672669641900137878","0.619858482633049717814250921847548870273803579893728549732116905572783351456693746575693111263986232610693447501950892844"]},"3,2,13":{"tau":[0.805598038704657,0.6470104332793236],"highPrecisionTau":["0.805598038704656906979726223934545760644162181530122755764675693225228797263046063586943358494337855071967673194613048673","0.647010433279323591335169143601300572032108730424817267566249216186489548934444748718467879079357141334248774768147972819"]},"4,0,14":{"tau":[0.9466149565782557,0.6579263472558945],"highPrecisionTau":["0.946614956578255706294279478944597763948700334216251641782709912299009011290517990562412622312970182002093328480451966722","0.657926347255894470371645494564482470461171846651801929218231205324395831286629783539589295519448065618496679528335165939"]},"3,3,12":{"tau":[0.7466364253924029,0.6642690791021805],"highPrecisionTau":["0.74663642539240293989571656868552262337429886234082536125369900405987282846884794720681873847767813363975378845246833891","0.664269079102180459031174672741170023332893697379249473587615173222532712751781512219308759323978101285047367175699827792"]},"4,2,12":{"tau":[0.7994935584250712,0.7064286957599494],"highPrecisionTau":["0.799493558425071215824589791444753075651538271694546383711755321415070669570134813414768319869087545795160968495505376949","0.706428695759949338825102473717484181776565696422870007710254508743458867558105653179738947110396277519752276389787887587"]},"5,2,11":{"tau":[0.7911825295943647,0.7660104734139578],"highPrecisionTau":["0.791182529594364755201852906313036210259635648524046260251549251867339926166856584755927330115732800612435484143627440795","0.766010473413957758214927734457708307875322056881809055284946669324045047701288890747926056636421245069394875091225545514"]},"5,0,13":{"tau":[0.9454339482543741,0.7192708895983165],"highPrecisionTau":["0.945433948254374143650138499013139596714312918964178781450101442509089462345664815772605633321208100836599720132366156768","0.719270889598316581205936001311617747593572682682499493643143278346949688589419965173515425762955592897191347426974757929"]},"6,0,12":{"tau":[0.9435061047926451,0.7809717402722857],"highPrecisionTau":["0.943506104792645040950554133421948031755759213252195973791965756172010328226485388251975066672598865136475053745642255545","0.780971740272285645034956029730292145103865484359440664673062994456173025031740703055645663233880750715054347345718806841"]},"0,4,15":{"tau":[0.7202602755977879,0.48924429683721],"highPrecisionTau":["0.720260275597787959100956684988942046797146717050122891462011293026756315985109036892344882390711467330766516181660877262","0.489244296837209966502387196350428786650661658121208914607626126416391502787103523630915897870386695887801912085753976362"]},"7,0,11":{"tau":[0.9408654632153685,0.8441944370861694],"highPrecisionTau":["0.940865463215368408286922857142181303302779374389406934491457600328217738895012195839883824230036964297590961975372082524","0.844194437086169355089100027212189248337093558694987384367742477697632208918268974196424785249101439756231621942635748092"]},"5,1,12":{"tau":[0.8629875741118789,0.7460039204327972],"highPrecisionTau":["0.862987574111878889341473448856412509279070369574978061088166727189193562402015430970957872498109469103328037531345977508","0.746003920432797201833447895264740049154033665469611169458117099236378744958107688153140609879483298404345877608417018628"]},"0,6,13":{"tau":[0.6406035744820531,0.5154203277678647],"highPrecisionTau":["0.640603574482053095595202562712899962822523417470096147078729290243322087224136717120766554318785437602731490190967645684","0.515420327767864701187942526095900503484258334673513718196952264961655971457302957203219656772741195363075693584449555974"]},"0,7,12":{"tau":[0.6023442098757881,0.5237223344933839],"highPrecisionTau":["0.6023442098757880412832696242154872300749855736887085670459821223449592976605655113419067950469033033290380240546498012","0.523722334493383852392256455185633125164964804406949483497668604312831660231595043264303725195507734177408333648054053755"]},"0,5,14":{"tau":[0.6797525319683677,0.504093129842842],"highPrecisionTau":["0.67975253196836767614712196030204585371024057848328306506237393439226481869396021414287259950188789809132645569254124052","0.50409312984284206735616411036819531892569821478628383993707387738803096623881047690609404018914671066449029579072740014"]},"0,8,11":{"tau":[0.5646487734710051,0.5293041800583481],"highPrecisionTau":["0.564648773471005038060927505811596838650728627559755126663208336701675665495387666755818385991715080933622880397510096151","0.529304180058348096878785427029886072779509292782960489086216463368345078448729266348833271276867733181578040333109105521"]},"0,14,5":{"tau":[0.3375865384373533,0.5097677436014357],"highPrecisionTau":["0.337586538437353337124926676903577758718622245893575170529843328408155603027428537571104535007652329420893254810227263413","0.509767743601435729456921122625191201419746560181101033616915853476575281780196050055805134619353613778270654079392435632"]},"1,4,14":{"tau":[0.7190238979859872,0.5525649653300928],"highPrecisionTau":["0.719023897985987202992903943174212073546404265210798403483725890553465416536210143490272620405835212953923866931651171993","0.552564965330092855012943130566219008418519716371805093261264080744753367091443829603800550912930500732440658546601258311"]},"1,5,13":{"tau":[0.6747618382564949,0.5653988831973557],"highPrecisionTau":["0.674761838256494925311075591941689329925391513226952183854013930109693950060522960420550387648972453290567517862989917445","0.565398883197355603737120503256540705841274745891851516610921029030129464539675496779474950538380433703666967967600986588"]},"1,3,15":{"tau":[0.7655250636391404,0.5360226557128306],"highPrecisionTau":["0.765525063639140429414723848591578306826053953583700583700052280205245370813797468050677465143336546068989093878412272474","0.536022655712830603975520143854073302077226358879076969274570539245649889073175908231514615687694707926052003221539175651"]},"1,7,11":{"tau":[0.5900543264006958,0.5819249591724268],"highPrecisionTau":["0.590054326400695720199508261434542340519690797465360875807053357948031439273814011785598951508537740120520112865992623828","0.581924959172426854121657184904910179182348151447094028025702791292784722489201462608017069920473527131927061251186281421"]},"2,3,14":{"tau":[0.7630670354888511,0.595368385830153],"highPrecisionTau":["0.763067035488851159419534579237396668230604951457311094459571730802145677052862305759653049873171173252283697779539551473","0.595368385830152895170109683669603411004414704602712937111126183739799324500147277083613997104582956143213772523047299746"]},"2,2,15":{"tau":[0.8162678564632856,0.5765607079605579],"highPrecisionTau":["0.816267856463285630724605097093167263065088416912085716509733931113663116620344748577851838140576832259864286473336400674","0.5765607079605578860651991000163144748671111634357766950677014179895296139041356131461917122522107936995325523531730091"]},"2,6,11":{"tau":[0.6196938294512346,0.6298572294341847],"highPrecisionTau":["0.619693829451234596596644715342512488390226340171460415400016127753122439460768238482055082273090262325994560162156008542","0.62985722943418466465454484909757994249211857476689926618968623200403689330650989470509208107907248600451024002886229267"]},"2,5,12":{"tau":[0.6658234599354972,0.6214897926853833],"highPrecisionTau":["0.665823459935497149771073774921754296362555082433155290652531555816995356279271928221994282014251622467881762419009443924","0.621489792685383263175835656611985153964082895969024190345925526716915031854865713908457785410447442368133446310096935744"]},"3,2,14":{"tau":[0.8136058577212028,0.6347656954864772],"highPrecisionTau":["0.813605857721202708004922237771986857305397942140432594667400021986088372893553623462804565644622851621672491187064190417","0.634765695486477236284417261656253417646460610381020973282281065525075494118208471631281398891859668580731336076407536782"]},"2,4,13":{"tau":[0.7133544873465288,0.6101279189563364],"highPrecisionTau":["0.71335448734652887859191635417962535159605118229233853041871684187423154261837749061576411773687530829202094856088204068","0.610127918956336467524546614237580118741621171256543080147147657249226070593481960546260221368129594525791933969167086771"]},"3,0,16":{"tau":[0.9488912385726641,0.5843247545468732],"highPrecisionTau":["0.948891238572664096696609310187986170795466545704074893233691189587852180409036995831015985096195665924546183041451220126","0.584324754546873106680984562973031523314797679213198233683452824275248606978727821719384789220153016594956608613152032109"]},"3,3,13":{"tau":[0.7574292888309982,0.6519812285803118],"highPrecisionTau":["0.757429288830998239148953771775656640818953924329102626759479552049770643625298225321602176771924896679644642262155537405","0.651981228580311768927740916953292717695277356790307319950040084435613823756884261991282901518054611445754262769795400632"]},"3,5,11":{"tau":[0.6537848784993011,0.6755849779903084],"highPrecisionTau":["0.653784878499301116114322868476728459903738857148259589104166230900727726865225787228510539184186885830597590895515820119","0.675584977990308375093853855039287748505054030870118092409027202173621139348758205704032715141069749202438889157992604481"]},"3,4,12":{"tau":[0.7045680566086097,0.6654157244673643],"highPrecisionTau":["0.704568056608609697996602679832052421435230387229280554876876005618785023884166976230494435248165159215506685627342371009","0.665415724467364310243064133192775489989434027321409380077079492989506600796207023550843294488158626534680621804556948382"]},"4,3,12":{"tau":[0.7492759608481511,0.7076260469717808],"highPrecisionTau":["0.74927596084815110326189378106341411697849726093105133774845018008802210001140246206126719782896448321362417895910231562","0.707626046971780745686837306612747973735901217924105719962851165699628874275464621278927482157237760461487257241342330219"]},"4,1,14":{"tau":[0.873329078031105,0.6710226619555788],"highPrecisionTau":["0.873329078031105093933083764429101639296004509929366994222287323875030659633918111693076963635877803351036687562516707442","0.67102266195557879921769026364311141238709645626529340475435947079807918431292874905492320068290217688698178685248906369"]},"4,0,15":{"tau":[0.9488676986526617,0.6442089540765913],"highPrecisionTau":["0.948867698652661693499430365212677873754089952116667222949855720105449321987474725581784546238686745332485764129125378053","0.644208954076591277354079548477849023477403494171248547921935881879732051522153375576978836448438150534077073914191683442"]},"3,6,10":{"tau":[0.6043150605882728,0.6828100895958735],"highPrecisionTau":["0.604315060588272768613037882028754288627749599644610523667921520788976489214672541052157949165979969754854280400043354922","0.682810089595873473090962652615422132921488071007152719322158673952282906045265357191186075357256658986016607203866546129"]},"4,4,11":{"tau":[0.6930885930662485,0.7200018533486527],"highPrecisionTau":["0.693088593066248438153940889476926238305631062193534062617905355704058402637805347843660631363159986293722234351016789819","0.720001853348652689174737689481552485333098374039440439388732174033967793408090244691792284779743384854903691953741732668"]},"4,2,13":{"tau":[0.8086175895918946,0.6915580429514468],"highPrecisionTau":["0.808617589591894531002689856420415590568469481475210813137694483297913214592786845020354368051424174808500313569403071408","0.691558042951446900087555345750526516329549502221120592230152795305800576466215831156414282989554433219126509529605656273"]},"4,15,0":{"tau":[0.06408892578470753,0.6566166412699841],"highPrecisionTau":["0.0640889257847075225252974759096627423201736394510369977075367817261265711226733884914987429389349259805032929268531146442","0.656616641269984094928463426572830659430462237103615446064319693150516352987568984504088577952245894348208357302739515799"]},"4,14,1":{"tau":[0.140272521326507,0.6802179728031825],"highPrecisionTau":["0.140272521326507004357489592294231700373937890547559260196116377891674541424440709056533163137453624161192802786131610924","0.680217972803182485151911134840083317270193549360134796580342714176677731067684259656918193016864950035012137214516349576"]},"5,1,13":{"tau":[0.8696433443751258,0.7285739629098119],"highPrecisionTau":["0.869643344375125843938682258540584024126314861258071456069338563263756790723834587718688962246874444238056874995145243898","0.728573962909811873583675676041027881692069846746886999025333318075644207693154982303200603828160552778852104593050447829"]},"5,3,11":{"tau":[0.7387973320422458,0.7634191700679737],"highPrecisionTau":["0.7387973320422458022487902884336735976953964536015727219393862965008536504043375648213300584321113031962607269208198529","0.763419170067973613250429148747883741694060164607258762120459118474226749692462278956528625142064224919122403626338728518"]},"5,4,10":{"tau":[0.6789574087673911,0.7748897453509634],"highPrecisionTau":["0.678957408767391108393850832915210776501854560694113009924480498456370602064175264804173813696917773823069254326270470947","0.774889745350963466767307616210438556211788115150722687040360871521540828993913530574400759925461834685883852062268962746"]},"5,2,12":{"tau":[0.801646268075934,0.748221861432732],"highPrecisionTau":["0.801646268075934009165995913431704316796245886658433657553426382375077470464861300263524163644476603144106657485220366649","0.748221861432731911235014910363047429376628780663832943979983635130473925649152947487353135907968951653971102185474758234"]},"5,13,1":{"tau":[0.1441208393624314,0.7377753386834314],"highPrecisionTau":["0.144120839362431389188700601516632487821320016875469152170201465958475416649017893968710953416218796424182255481288285981","0.737775338683431370523418494185194362444724211599427348893961278553593398740600879866162194221847331862899965192645103259"]},"5,14,0":{"tau":[0.06496508364335254,0.7152261920925459],"highPrecisionTau":["0.0649650836433525464787618926308103466935781201655057792536305796565238605452345794899712697565520562793617138327453499955","0.715226192092545845661770240291698985130210206050223057539207254134229083206102627010458953797797834322626711034755242098"]},"6,1,12":{"tau":[0.8645581421832693,0.7866460463129716],"highPrecisionTau":["0.864558142183269292812581814811081846704333471959336485826983577403042886541959541115947799031410502344465081325628763798","0.786646046312971668511014014096945524751039151562797400912823706444260917635457182892042428823491129108759787088293397345"]},"6,12,1":{"tau":[0.14959723595598146,0.7959613407185725],"highPrecisionTau":["0.149597235955981456909976379574074076428089752375961159694537992240918348524929930083146471768369664658828705944692838361","0.795961340718572547278241841021812840780003784570988836207847909945512835289566746322964675092365681168591976547374301148"]},"6,13,0":{"tau":[0.06676113173025358,0.7740052838018504],"highPrecisionTau":["0.0667611317302535870323401683196062455832473499978730076517667227982277278671881871723686271362363679564347088107897865888","0.774005283801850396381684923588249111832403758185304214671588609539510933834756734811976931918244917273957588169829865979"]},"8,0,11":{"tau":[0.9414449633802877,0.8823107855233865],"highPrecisionTau":["0.941444963380287659222350443358254013627516223151972346740325021546002556170066676721032824834223674793033171831059361729","0.882310785523386528824863968093624430345413123978383189662460801243006730477437291531012480220917875829774637703648944801"]},"7,1,11":{"tau":[0.8580661839300548,0.8461637778086981],"highPrecisionTau":["0.858066183930054821076367065787313659072872415648845308355748094833065637965785914328351881168104803487711207229530538651","0.846163777808698107113213250449480805960889390523680979600478504327223351028568503986387076851160568196022829721624811998"]},"7,2,10":{"tau":[0.7818108240571182,0.8648821522938819],"highPrecisionTau":["0.781810824057118227008282137670671937431784710510000476930261817385657046845816550772475064886574206961737561784153189163","0.864882152293881823696210085413257364828320340169931178585703986789083800781884765420487102209467512224192046251925185599"]},"7,0,12":{"tau":[0.944260864032012,0.8207648295660706],"highPrecisionTau":["0.944260864032011988381448379261239117657908113579522666378627160721342927951592172485681248322655599280289584317550747587","0.820764829566070562509252891216848463473706426350543891190687309263687603759818480619868080205604177405316427824712722119"]},"0,15,5":{"tau":[0.3270388001375981,0.5051296417530291],"highPrecisionTau":["0.32703880013759811422101814396416483172332027602039808201035333789031347343001806307080085868914261955328014610120223217","0.505129641753029095846152772453513200078915594846116341751018114334037754758471449149888034160843415442147140158814961952"]},"0,13,7":{"tau":[0.4011361230441501,0.5231396897941097],"highPrecisionTau":["0.401136123044150103436880759766006918851365604105773101149005503909151991043073199175467287207395083060425350182237154246","0.523139689794109722019709769172075165589010529820429483848542113264985145336762454552272630909969413930996816389293048063"]},"0,16,4":{"tau":[0.28822183664909107,0.49153320724500754],"highPrecisionTau":["0.288221836649091068078240808381147176588762787641544757399948945318052107474773346691296423597034986991515146212927141746","0.491533207245007515204380010922926022597691165516146992604411858971731520988637486002647834752191128588822024688035645629"]},"1,15,4":{"tau":[0.2861583436049288,0.5512548822648211],"highPrecisionTau":["0.286158343604928825547843723367639105616024716421603730985250093323585295913940383363486222702161692073594163519405954314","0.551254882264821082597089989458095357726639545674913268127248462015110813526251515501068545704823589582569062982508064453"]},"0,8,12":{"tau":[0.5794265410977651,0.5260835747826774],"highPrecisionTau":["0.579426541097765141425703566168708327247186703831094362218821822038946406034139122468119440918015532795302676382366710152","0.526083574782677385639846999019818407368050449049021880943070804602827474218755500964782995619349186265399327273086827533"]},"12,7,0":{"tau":[0.09871516396555688,1.187300901205209],"highPrecisionTau":["0.0987151639655568747512951883179059453674215697569775439766084487366120097031919736075544726852910923664372447988645108273","1.18730090120520885991096102439621510123605902139789553924571189188002431373186403039637812267358845081112947089980561746"]},"1,3,16":{"tau":[0.7730955769246913,0.5288974208314254],"highPrecisionTau":["0.773095576924691312503653625748662208504222024055280099717292588822087235397175191519033386793454574237172186413031072258","0.528897420831425395811186663141321965275381243154375610692760035578285131529591572157617062892366903526826061133689904449"]},"1,14,5":{"tau":[0.32856538785355244,0.5631596824368247],"highPrecisionTau":["0.328565387853552453122807565292572074839252844043315660621222139904958039351137710646077813311866903937905321145355036901","0.563159682436824712316753426930739810405777728487228349981905711454730987058511784215610605677213641874073815645443541282"]},"1,8,11":{"tau":[0.5657918406981716,0.5811553571809746],"highPrecisionTau":["0.565791840698171603115733834128105685163140504050364650517387412778729643933884804547817544077506452256923427938513057062","0.581155357180974500952187328902263735572998409198466564731989267567949736824162423222675273380254840838737159180650670263"]},"13,6,0":{"tau":[0.1101531704019854,1.2785335704162577],"highPrecisionTau":["0.110153170401985403802384277113073659733646223843993650684538308623144516848107468180641786933054083728128548229852105052","1.2785335704162577184459640015505058013718328798071998352745203117981630885519128952936166505906241907950672618908104712"]},"15,4,0":{"tau":[0.14644183154375248,1.5035098311190147],"highPrecisionTau":["0.146441831543752470502002935409029573441656852494210412712536522214385592131659105076476754311918578882772285685231249496","1.5035098311190146922952908849336555618750816957614975055966487476280637390854741589951910989520413093112707397573873306"]},"2,14,4":{"tau":[0.2898377924446902,0.606133455757789],"highPrecisionTau":["0.289837792444690222407233727671725240040355090791242673971921487324496794592382365067650984486906997616478617470669502022","0.606133455757788942508355294619335120373117635371185221272552769322051743179492080358427451809358186760757989624496343526"]},"14,5,0":{"tau":[0.1253610047921157,1.3823132838719612],"highPrecisionTau":["0.12536100479211569350640052624799807304358840055999993285522996946579531996328513005018673362454031926177180756235192196","1.38231328387196114024743534507731100467800562104193209954856567257249990977406301864922593429383944755549582949017960258"]},"2,13,5":{"tau":[0.33528776122675835,0.6166813735277872],"highPrecisionTau":["0.335287761226758335560177840183808460226302891151432203525335218835398621274964372802018925565461479750172223220724303526","0.616681373527787201935636453664463905550361036487076113328606306907851744303421947964722600851392355600125312644180062172"]},"2,7,11":{"tau":[0.5916778607567114,0.6278742389130555],"highPrecisionTau":["0.591677860756711355761910587489910931172646898483165089221559628171194940668140557720645838503392260855720845936872055517","0.627874238913055459695795442635401254026778287329833299086260329982589397339834425675737808151365133790771638529506416701"]},"2,15,3":{"tau":[0.24219297044694374,0.5924955244550538],"highPrecisionTau":["0.242192970446943745380280185107463525198079789043039154718120121689464128931298376551525604464847822576153759458262763085","0.592495524455053812231453314095431009177913043519302300596624061262441682333163775730249976927526234369752538638183546207"]},"2,3,15":{"tau":[0.7716016739690534,0.5861887039402627],"highPrecisionTau":["0.771601673969053455003981280381425106891438636369362904409729410433933537187289021392185218453800471885400958402966101837","0.58618870394026272335737541186862877358925014679713964691120064890376077133906366889311723154580229725420176172574325475"]},"3,1,16":{"tau":[0.8800594314334781,0.601958827484376],"highPrecisionTau":["0.880059431433478129323800501997612155111302031542514509574955117200962206627021359161005849880515715969945589891758813016","0.601958827484376041445025437938523013653259761625099083581719762469269304026528996792429350606077433841562725739717583318"]},"3,14,3":{"tau":[0.24642567006040775,0.646558032079764],"highPrecisionTau":["0.246425670060407750063604504560908189880574003246558037341122063562463606037126359987037844356572978405919234702429111385","0.64655803207976402038800839661815763733734508917892590904042824909969156984488162616639656603333556110908396884451684248"]},"3,16,1":{"tau":[0.13306362232229058,0.611085942185133],"highPrecisionTau":["0.133063622322290578175414610093881029267623605347595995076717725867736558285636097541863919949410690259654569700792600993","0.611085942185132982610984902249460725892068241691413579911071759142654902215263267003223102312855493645622250548672856938"]},"3,2,15":{"tau":[0.8207736410443663,0.6235093261573438],"highPrecisionTau":["0.820773641044366353196147161689494277606499383625957020349682752518990053401464999387276139139026251461790824213070818543","0.623509326157343799866803667903044190887203913720159862043393960810440193821166306591138112653010689662516864384809470561"]},"3,7,10":{"tau":[0.5753124667947556,0.6778057140826039],"highPrecisionTau":["0.575312466794755570208067289514958399791851472315068883979200034911309901930305572592541703178934676363793841277590298503","0.677805714082603951535065374580748599523349320145027966075730013213897420724264818167930930059652354410132685913714881802"]},"3,6,11":{"tau":[0.6214051558506005,0.6725506320967659],"highPrecisionTau":["0.621405155850600507995846092794563258702173193236176460258369564193116249790747599975927333297696388637282983994579888493","0.672550632096765885385242140971849411482942076429078045878898262616911606611206626671985883703580300580020563064670379432"]},"4,15,1":{"tau":[0.13443572071810103,0.666713699917096],"highPrecisionTau":["0.134435720718101041814895174011600341449284910582644374530674874668638335888892263719730300771557081104319703010571748131","0.666713699917095933928410660302648310344704738988273591515215540189737220158856063327148353891565060926257867831223305631"]},"3,5,12":{"tau":[0.6683542514248495,0.6647609075879546],"highPrecisionTau":["0.668354251424849457452531294667733127095305610831135037081390672223522149326922707554719880953062682606052175327906890335","0.664760907587954542680235350747414408402222630482020511907745372020931461816504269385161695751547976456452354232713897003"]},"4,13,3":{"tau":[0.2532963074838741,0.6996414232485513],"highPrecisionTau":["0.253296307483874096565364568727632012278518145971132634731479281465103501709189995092788892344920134824440562230368704861","0.699641423248551325034676403537133599699569207642933814213429630593088116465754642349200197574660225716635388574570456916"]},"4,6,10":{"tau":[0.6053608149577492,0.722915647763107],"highPrecisionTau":["0.605360814957749147051287851553849179320959873412514093544295600324225060795883782437025362839291144912019570592142052019","0.722915647763106975183281636794815354958905763305679960732122162870958105770629701248025935041925737449615081989131440042"]},"4,3,13":{"tau":[0.7602133671640198,0.6939725199954984],"highPrecisionTau":["0.760213367164019751671999180846839940102177882776299727778883320991040933793068086446541436377849732653997942816410815502","0.693972519995498318535780485937307208409745808511830952994501686524881317850559476695742035092070172796462862692797119245"]},"4,4,12":{"tau":[0.7068459267683825,0.7065444637861158],"highPrecisionTau":["0.706845926768382573512204854219245743261725054516360647894181001428176397209673976488838195023725345602218793715426300351","0.706544463786115755753109594004717561924587079197701389800377588772392166639440999474924738437352395646103954753567012411"]},"4,5,11":{"tau":[0.6554694076520992,0.7161027344623897],"highPrecisionTau":["0.655469407652099175343060781134142518308133826584098225183553587441534650544983944740846241569078807433380965845067146511","0.716102734462389720632697684936883730557553558177321184855620151558319505141253320147976867716352586202392800688370226624"]},"5,13,2":{"tau":[0.2025877436844616,0.7390694255371073],"highPrecisionTau":["0.202587743684461615270069766961192820082896098091285069251106175331206581554801480505899189967451238239004370562306066974","0.73906942553710731217791958985345424919597621220635764204608252998970862411167790728464716807987314893371347814984646985"]},"5,14,1":{"tau":[0.13751923153617707,0.7216785415057605],"highPrecisionTau":["0.137519231536177074567644293286349597047060702314176329763463283936837460901868236271131808928912661485622387048313226382","0.72167854150576047691992108394392755216361999117543086944843660738055865566465967308250410245332670332320784752896105421"]},"5,5,10":{"tau":[0.6399757253980529,0.7676133315984096],"highPrecisionTau":["0.639975725398052921548534442967162888804555316612159876301624477736204962586697838496800767393649381765903515290424418718","0.767613331598409559735539588676999735759385519356635654738034368733537319475483175713643177827701621227407499768802473571"]},"5,4,11":{"tau":[0.6946815928878074,0.758939862354412],"highPrecisionTau":["0.694681592887807383555578613557797216421374340669099956133338761748295446570924638886079145561449035436112206711324969712","0.758939862354412055358321531465600182082220337307045915493057815979129146200855190830018240820336781782190578739750757081"]},"6,13,1":{"tau":[0.1420676530448785,0.7769637628738744],"highPrecisionTau":["0.142067653044878519383273731411140760190961594688810109311400773105440003315577397652471703772819114921243761389374942673","0.776963762873874316623926618655367591602721255476309605484443981837511030399869783403668397749107545930097524112386304736"]},"6,4,10":{"tau":[0.6800703184805613,0.8121383768380686],"highPrecisionTau":["0.680070318480561355221374105192478427640022057138054405655519502218150199309668608952695104781383608145952698855097187401","0.812138376838068595504647747585135139413259404859712803516310648362705481352348676791530756639159463566430312342247615451"]},"6,14,0":{"tau":[0.06347078638176887,0.7556249555311689],"highPrecisionTau":["0.06347078638176887201642535263766006719816021304689036190293854134791979475484222587567243001255679894408389395502080897","0.755624955531168923006601172027580011483454681215918881514964905089803858680174672105595450089032416616515564786281868423"]},"6,3,11":{"tau":[0.7402438174539494,0.8011609723118959],"highPrecisionTau":["0.740243817453949381148135971131896940574349634933893936695992421649883690753136555581019583434933456838491913648372625626","0.801160972311895916275662027053615333220160773834622683440424004369689303633720586250813877867779683160389121893021089057"]},"7,13,0":{"tau":[0.06561983061957288,0.8123527171139189],"highPrecisionTau":["0.065619830619572881617495410339086737245105049846925713199758417250048375792192935158725829182390739104676638853171905629","0.81235271711391888372484872186877741177906825396336675528565366047333677855578721032389368768608576715275381841049488808"]},"7,12,1":{"tau":[0.14802865365482962,0.8333793867822267],"highPrecisionTau":["0.14802865365482962243715813868741446864260454938238850921464747758038870549027526411909154731115521654708759841923385742","0.833379386782226648357719641982750064539058664774849532184167142152669434389079151316345542835207049238905796873683031108"]},"7,3,10":{"tau":[0.7269808548372917,0.8565026214688495],"highPrecisionTau":["0.726980854837291721859775771065587005092102185060893673620287342502804060876713650152147333290868300143731228395172966397","0.856502621468849502481561520386512242397663621319556006229685931504103878444159834667558185616924489696286933790055387673"]},"7,2,11":{"tau":[0.7939885738362968,0.8425360352459708],"highPrecisionTau":["0.793988573836296855307532669377687576459155121109456913944085294454934760444641757272549822495472127625126273737968189762","0.842536035245970769480760827504221165091641681919751345031814032729389981387173258015583615111721922654649138117500075411"]},"8,1,11":{"tau":[0.8590035847323096,0.8822642628981217],"highPrecisionTau":["0.859003584732309631614342347473079850461320752944891750357325259781980002457266163956156775375740648564202119602650419634","0.882264262898121664521253688187262167907439834083997237935588192690366003263148787284936881700674155356647999218852441717"]},"8,2,10":{"tau":[0.7827392025675338,0.9004053325121216],"highPrecisionTau":["0.782739202567533801834049231114827450982985890464877484102684610430791841235400942629159116990689682517337184264954103715","0.900405332512121619719579194740663200577798730709612363634824321329061814484653256474414247207031770833825050235976267645"]},"12,8,0":{"tau":[0.08949432372228064,1.1382025977363301],"highPrecisionTau":["0.0894943237222806386536330958796761114826096462855303152177144921442828263542679532444493243481154934041622962695614072802","1.13820259773633016032758968180402378805678429532515382743895309584131947287829683758700890739639985730873674600714851605"]},"13,1,6":{"tau":[0.795096770955491,1.231509202187628],"highPrecisionTau":["0.795096770955490998005242213132453444565678450518886482916854118567875565367911176608584026589640927125328386554022117984","1.2315092021876280167533317240426998443716050789238036848434274107674160611171357692897168493220039164649327825467181535"]},"13,5,2":{"tau":[0.34409554521129565,1.256000064114844],"highPrecisionTau":["0.344095545211295663127584986541693297120256514023397544823901548808266808851791223317480935110649151966661104015355672663","1.25600006411484413480849250369415206990563013191858890969725325071778212148835454721592446669830459726698168600559356465"]},"9,0,11":{"tau":[0.941910277677223,0.9178648205820933],"highPrecisionTau":["0.941910277677223012563107240237318798372505541411525567982926389811015128910901095748447536953757932947129202557005439465","0.917864820582093321231762139462869341406918051739221637882707880936246345509641308193899393203896494496807293657973905379"]},"14,3,3":{"tau":[0.5134570053075578,1.348022737989478],"highPrecisionTau":["0.513457005307557827009837590155906909742075272973428707366132410923162671180470513211488372239113079170364736338227655838","1.34802273798947797827270236706667204137978312196149723251206390030053553048484867057879277556646036729948468570264225849"]},"13,7,0":{"tau":[0.09842164950333812,1.2194948744643062],"highPrecisionTau":["0.0984216495033381240662152726249113497735517322295178418009721950893830116728044269268939918961317925916841930166128019676","1.21949487446430612155883176910641055201437454951921475538307559796675662266005717577640543627048349267350145231525685211"]},"14,5,1":{"tau":[0.25397175601001654,1.333964206648046],"highPrecisionTau":["0.253971756010016516143298180185974930884198420942935309827057326705412840208371119321134508002107320217016148444088964084","1.33396420664804591181596830456735134036074323217324400278838291212009686511789121718309198972670657049296544480140653696"]},"13,6,1":{"tau":[0.22705858981560748,1.242576456961025],"highPrecisionTau":["0.227058589815607476010873388675196199153098603961836047456268185312482869630088969282946855693167230855224348946124019487","1.24257645696102489126529440816499814704720002995981043780107690620029143583668821954758673921456678828046398878733704249"]},"15,0,5":{"tau":[0.9011252788313059,1.3925420454067288],"highPrecisionTau":["0.901125278831305882689236797961691819603444621870447612576052814916701364615037998358346822232131246614917042357579333997","1.3925420454067288627006674324033739815121794319657261677234116948012448987137770054328245878529013549422747836279699628"]},"15,4,1":{"tau":[0.2895814541926095,1.4377744493109312],"highPrecisionTau":["0.28958145419260951225392217598505042475910766113757975507813739707947437529723098882382238724591903910913121902746053243","1.43777444931093126885342851938557623624648127852333609332571131728823172608799972539890297570994312348592282279320079135"]},"15,5,0":{"tau":[0.12518231270035904,1.4136085346428897],"highPrecisionTau":["0.125182312700359043747896490980778712556568142756980866554060312110552762840639285686170614565263025745060676272139212097","1.41360853464288967174429838272261427345527776032289175771795439256095952697150332321695272147703591963663906078271887064"]},"16,4,0":{"tau":[0.14630879596372567,1.534479558168463],"highPrecisionTau":["0.146308795963725670450505820419003379948118964451678350430111805418702305306719167771236877715690667285713550975127954747","1.5344795581684630080651967498130039146609960905195034675423194402160073347490641707834630157127991967253833103069385711"]},"16,3,1":{"tau":[0.33882803936870937,1.5583426115290213],"highPrecisionTau":["0.338828039368709355894650979375340903215770133653804804954205925234628294798492331001773636301423939348664195989591878761","1.55834261152902140486729165412714017117959306848547302125974995192322810924340980437958113102241324355382196146686715497"]},"15,3,2":{"tau":[0.4412983890260613,1.447055611117579],"highPrecisionTau":["0.441298389026061285219230346547531687062199843298529330728401284490477167872646858120081245466904689525151327016069384075","1.44705561111757904921908754902403843551466482829806422528960669865262629239048821339985412496934139695005368484516350079"]},"16,0,4":{"tau":[0.8853521381135935,1.5110629130444664],"highPrecisionTau":["0.885352138113593458368608921086254906743655897492186986439727320076493860096941163860398061712338345628147637669841132548","1.51106291304446632316920451779020012589245333008083080381018301685854840736921705162005521748956549659274916347895635906"]},"3,4,13":{"tau":[0.7166651846562491,0.6542329269602049],"highPrecisionTau":["0.71666518465624901967063914078321545160896128866552837625017064122592039226538203620463385321218955504926612338123578786","0.654232926960204896045270868363051588858056002775298014627471806302398344381826548104322988292030227489260059267554498772"]},"2,5,13":{"tau":[0.6786437035901883,0.6127217756926099],"highPrecisionTau":["0.678643703590188318314579843059073313374384476625144426666662058134202873258960812676580368792755360391688606139879206155","0.612721775692609881407138922079627036086704178621192893302094787245874108917990990689772791048660542300762880617764042904"]},"4,2,14":{"tau":[0.8167212279100172,0.6779594334093928],"highPrecisionTau":["0.816721227910017180509944465774177772360925388326743489561739482291818395273794503547499172933436481711476735181594633983","0.677959433409392782411421707620371680076813831574950611179833775949521160078306570161222457053731386454527687270458570742"]}}};

atlasAttachment.setAnchorData(attachmentAnchors);

const POINT_INSET = 0.04;
const POINT_ROUNDOFF = 1e-12;
const pointStateCache = new Map();

function compactSourcePoint(r, s, inset = POINT_INSET) {
  if (![r, s].every(Number.isFinite)) {
    throw new RangeError("The source coordinates must be finite numbers.");
  }
  if (inset !== POINT_INSET) {
    throw new RangeError("Interactive points use the fixed compact inset 0.04.");
  }

  const weights = [1 - 2 * s, 2 * (r + s) - 1, 1 - 2 * r];

  if (weights.some(weight => weight < POINT_INSET - POINT_ROUNDOFF)) {
    throw new RangeError(
      "Choose a point in the compact source triangle: " +
      "all three barycentric weights must be at least 0.04."
    );
  }

  // Pointer projection can undershoot a compact edge by floating roundoff.
  // Restore just that tiny excess; do not move points outside the accepted K.
  if (weights.some(weight => weight < POINT_INSET)) {
    const clamped = weights.map(weight => Math.max(POINT_INSET, weight));
    const largest = clamped.indexOf(Math.max(...clamped));
    clamped[largest] -= clamped.reduce((sum, weight) => sum + weight, 0) - 1;
    return {r: (1 - clamped[2]) / 2, s: (1 - clamped[0]) / 2};
  }

  return {r, s};
}

self.onmessage = async ({data}) => {
  if (data.type === "solve") {
    const id = data.id;

    try {
      const checked = atlasAttachment.validateCounts(data.counts);
      if (checked.status !== "ok") throw new RangeError(checked.message);

      const {r, s} = compactSourcePoint(data.r, data.s, data.inset);
      const key = data.counts.join(",");
      const previous = pointStateCache.get(key) ?? null;
      const result = atlasAttachment.solveAttachment(r, s, data.counts, previous);

      // Only a complete numerical success may seed the next point. The
      // original fixed geometric anchors remain available for solver retries.
      if (result.status === "ok") {
        pointStateCache.delete(key);
        pointStateCache.set(key, result);

        if (pointStateCache.size > 8) {
          pointStateCache.delete(pointStateCache.keys().next().value);
        }
      }

      self.postMessage({type: "point", id, result});
    } catch (error) {
      self.postMessage({
        type: "pointError", id,
        message: error instanceof Error ? error.message : String(error)
      });
    }

    return;
  }

  if (data.type !== "compute") return;
  const id = data.id;
  try {
    const checked = atlasAttachment.validateCounts(data.counts);
    if (checked.status !== "ok") throw new RangeError(checked.message);
    const mesh = await atlasAttachment.sampleAttachment({
      counts: data.counts,
      resolution: data.resolution ?? 12,
      inset: data.inset ?? 0.04,
      onSample: (_sample, _index, completed, total) => {
        if (completed % 4 === 0 || completed === total) {
          self.postMessage({type: "progress", id, completed, total});
        }
      }
    });
    self.postMessage({type: "done", id, mesh});
  } catch (error) {
    self.postMessage({
      type: "error", id,
      message: error instanceof Error ? error.message : String(error)
    });
  }
};

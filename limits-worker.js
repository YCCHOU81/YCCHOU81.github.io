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

const BRANCHES = Object.freeze([
  { id: 'base', label: 'n = 1', n: 1, dr: 0, ds: 0,
    seed: c(.5904953761444916, .7963943351098978) },
  { id: 'L1', label: 'Attach at L1', n: 2, dr: .5, ds: 0,
    seed: c(.5823785379137855, 1.1341826530853316) },
  { id: 'L2', label: 'Attach at L2', n: 2, dr: 0, ds: .5,
    seed: c(.3544868509081326, .6879997451009817) },
  { id: 'L3', label: 'Attach at L3', n: 2, dr: .5, ds: .5,
    seed: c(.7515095138734595, .6524433021803057) },
]);

/** Safeguarded complex Newton within the fixed Gamma_0(2) domain F0. */
function solveNewton(r, s, n, seed, options = {}) {
  let tau = asComplex(seed), evaluations = 0;
  const maxIterations = options.maxIterations ?? 32;
  const tolerance = options.tolerance ?? ROOT_TOL;
  if (!inF0(tau)) return { status: 'invalid_seed', tau, iterations: 0, residual: Infinity,
    relativeResidual: Infinity, evaluations };
  try {
    let e = evaluate(tau, r, s, n); evaluations++;
    for (let iteration = 0; iteration <= maxIterations; iteration++) {
      if (e.relativeResidual < tolerance) return { status: 'ok', tau, iterations: iteration,
        residual: e.residual, relativeResidual: e.relativeResidual, evaluations };
      if (iteration === maxIterations) break;
      const h = 1e-5 * Math.max(1, abs(tau));
      const fPlus = evaluate(c(tau.re + h, tau.im), r, s, n).value;
      const fMinus = evaluate(c(tau.re - h, tau.im), r, s, n).value;
      evaluations += 2;
      const derivative = scale(sub(fPlus, fMinus), 1 / (2 * h));
      if (!finite(derivative) || abs(derivative) < 1e-15) return { status: 'small_derivative', tau,
        iterations: iteration, residual: e.residual, relativeResidual: e.relativeResidual, evaluations };
      let step = div(e.value, derivative);
      const limit = Math.max(.12, .65 * abs(tau));
      if (abs(step) > limit) step = scale(step, limit / abs(step));
      let accepted = false;
      for (let backtrack = 0; backtrack < 25; backtrack++) {
        const trial = sub(tau, scale(step, 2 ** -backtrack));
        if (!inF0(trial)) continue;
        const candidate = evaluate(trial, r, s, n); evaluations++;
        if (candidate.relativeResidual < e.relativeResidual * (1 - .00005 * 2 ** -backtrack) ||
          candidate.relativeResidual < tolerance) {
          tau = trial; e = candidate; accepted = true; break;
        }
      }
      if (!accepted) return { status: 'line_search_failed', tau, iterations: iteration + 1,
        residual: e.residual, relativeResidual: e.relativeResidual, evaluations };
    }
    return { status: 'iteration_limit', tau, iterations: maxIterations, residual: e.residual,
      relativeResidual: e.relativeResidual, evaluations };
  } catch (error) {
    return { status: 'evaluation_failed', message: error.message, tau, iterations: 0,
      residual: Infinity, relativeResidual: Infinity, evaluations };
  }
}

const FALLBACK_SEEDS = [.15, .35, .55, .75, .9].flatMap(x =>
  [.28, .5, .8, 1.2, 2, 3.5, 6].map(y => c(x, y)).filter(z => inF0(z)));

function continueBranch(branch, r, s, previous) {
  const old = previous?.branches?.find(item => item.id === branch.id && item.status === 'ok');
  const fromR = old ? previous.r : .31, fromS = old ? previous.s : .29;
  let tau = old ? asComplex(old.tau) : asComplex(branch.seed);
  let fraction = 0, iterations = 0, evaluations = 0, continuationSteps = 0;
  const distance = Math.hypot(r - fromR, s - fromS);
  let increment = distance ? Math.min(1, .045 / distance) : 1;
  let last = null;
  for (let attempts = 0; fraction < 1 && attempts < 160; attempts++) {
    const next = Math.min(1, fraction + increment);
    const rr = fromR + next * (r - fromR) + branch.dr;
    const ss = fromS + next * (s - fromS) + branch.ds;
    last = solveNewton(rr, ss, branch.n, tau);
    iterations += last.iterations; evaluations += last.evaluations;
    if (last.status === 'ok') {
      fraction = next; tau = last.tau; continuationSteps++;
      increment = Math.min(1 - fraction, increment * 1.45);
    } else {
      increment *= .5;
      if (increment < 2e-5) break;
    }
  }
  const target = { r: r + branch.dr, s: s + branch.ds };
  if (fraction === 1 && last?.status === 'ok') return { ...last, id: branch.id,
    label: branch.label, n: branch.n, target, iterations, evaluations, continuationSteps,
    method: 'parameter_continuation' };
  for (const seed of [tau, branch.seed, ...FALLBACK_SEEDS]) {
    last = solveNewton(target.r, target.s, branch.n, seed);
    iterations += last.iterations; evaluations += last.evaluations;
    if (last.status === 'ok') return { ...last, id: branch.id, label: branch.label,
      n: branch.n, target, iterations, evaluations, continuationSteps, method: 'reseeded_newton' };
  }
  return { ...(last ?? {}), id: branch.id, label: branch.label, n: branch.n, target,
    status: 'unresolved', iterations, evaluations, continuationSteps,
    message: 'No reliable root was found in the selected chart.' };
}

/** Four independently checked branches. Passing the previous returned state
 * accelerates dragging; every returned root is recomputed at the new target.
 */
function solveAll(r, s, previousState = null) {
  const status = sourceStatus(r, s);
  if (status !== 'ok') return { r, s, status, branches: [], margin: SOURCE_MARGIN,
    message: status === 'near_boundary' ? 'Move slightly inside the triangle; boundary roots are not estimated.' :
      'Choose a point inside the source triangle.' };
  const previous = previousState && sourceStatus(previousState.r, previousState.s) === 'ok' ? previousState : null;
  const branches = BRANCHES.map(branch => continueBranch(branch, r, s, previous));
  return { r, s, status: branches.every(branch => branch.status === 'ok') ? 'ok' : 'partial',
    branches, margin: SOURCE_MARGIN };
}

/** A mesh of the source triangle, inset in barycentric coordinates. */
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

/** Async sampling hook for either a worker or the main thread. The UI decides
 * whether to draw incomplete samples; failures retain explicit status fields.
 */
async function sampleAtlas(options = {}) {
  const mesh = buildSourceMesh(options.resolution ?? 14, options.inset ?? .025);
  const samples = [];
  let previous = null;
  for (let index = 0; index < mesh.points.length; index++) {
    if (options.signal?.aborted) return { ...mesh, samples, status: 'aborted' };
    const point = mesh.points[index];
    const sample = solveAll(point.r, point.s, previous);
    samples.push(sample);
    if (sample.status === 'ok') previous = sample;
    if (options.onSample) await options.onSample(sample, index, mesh.points.length);
    if (index % 8 === 7) await new Promise(resolve => setTimeout(resolve, 0));
  }
  return { ...mesh, samples, status: samples.every(sample => sample.status === 'ok') ? 'ok' : 'partial' };
}

return {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh};
})();
const attachmentMath = (() => {
const {c,add,sub,mul,scale,abs,evaluate,inF0,sourceStatus,buildSourceMesh}=atlasMath;
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

const MAX_N = 12;
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
 const counts=options.counts??[0,0,0],checked=validateCounts(counts);
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
  if(options.signal?.aborted)return{...mesh,counts:[...counts],n:checked.n,samples,status:'aborted',completed};
  const p=mesh.points[index],sample=solveAttachment(p.r,p.s,counts,previous);
  samples[index]=sample;completed++;
  if(sample.status==='ok')previous=sample;
  if(options.onSample)await options.onSample(sample,index,completed,mesh.points.length);
  // Every point yields to queued interactive jobs in a worker.
  await new Promise(resolve=>setTimeout(resolve,0));
 }
 return{...mesh,counts:[...counts],n:checked.n,samples,completed,
  status:samples.every(s=>s?.status==='ok')?'ok':'partial',
  metadata:{sampledInteriorOnly:true,sourceCoordinateSideGap:mesh.inset/2,
   maximumRelativeResidual:Math.max(...samples.flatMap(s=>s.branches.filter(b=>b.status==='ok').map(b=>b.relativeResidual)))}};
}

return {solveAttachment,validateCounts,setAnchorData,MAX_N};
})();
const approachMath = (() => {
const {solveAttachment,validateCounts,setAnchorData,MAX_N}=attachmentMath;
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
 function value(tau,r,s,n){const b=basis(tau,r,s,n);return V(fromBasis(b,n),P(b.factor,n*(n+1)/2));}
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

const clock=()=>globalThis.performance?.now?.()??Date.now();
let storedAnchors=null;
function setApproachAnchorData(data){setAnchorData(data);storedAnchors=data;return{status:'ok',maxN:Math.min(MAX_N,data.maxN)};}

function normalizePath(path={}){
 const kind=path.kind??path.mode??'edge',index=Number(kind==='vertex'?(path.vertex??path.index??0):(path.edge??path.index??0));
 const position=String(kind==='vertex'?(path.direction??path.position??'.5'):(path.position??'.5'));
 if(!['edge','vertex'].includes(kind)||![0,1,2].includes(index))throw new RangeError('Choose an edge or vertex numbered 0, 1, or 2.');
 const t=new Decimal(position);if(!t.gt(0)||!t.lt(1))throw new RangeError('The path position must lie strictly between zero and one.');
 let start;
 if(path.start!==undefined){
  if(!Array.isArray(path.start)||path.start.length!==2)throw new RangeError('Choose an interior starting point with two coordinates.');
  start=path.start.map(value=>{
   const text=String(value);
   if(text.length>60||!/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(text))throw new RangeError('Use finite decimal coordinates for the starting point.');
   const coordinate=new Decimal(text);if(!coordinate.isFinite())throw new RangeError('Use finite decimal coordinates for the starting point.');
   return coordinate.toString();
  });
  const r=new Decimal(start[0]),s=new Decimal(start[1]);
  if(!r.gt(0)||!s.gt(0)||!r.lt('.5')||!s.lt('.5')||!r.plus(s).gt('.5'))throw new RangeError('The starting point must lie strictly inside the base monodromy triangle.');
 }
 return{kind,index,position,...(kind==='edge'?{edge:index}:{vertex:index,direction:position}),...(start?{start}:{})};
}
function coordinates(engine,path,epsilon){
 const {d}=engine,e=d(epsilon),t=d(path.position),b=[];
 if(!e.gt(0)||(path.start?e.gt(1):!e.lt(1)))throw new RangeError(path.start?'Segment epsilon must be positive and at most one.':'Approach epsilon must lie strictly between zero and one.');
 const other=[0,1,2].filter(i=>i!==path.index);
 if(path.start){
  const r=d(path.start[0]),s=d(path.start[1]),start=[d(1).minus(s.times(2)),r.plus(s).times(2).minus(1),d(1).minus(r.times(2))];
  const boundary=[d(0),d(0),d(0)];
  if(path.kind==='edge'){boundary[other[0]]=t;boundary[other[1]]=d(1).minus(t);}else boundary[path.index]=d(1);
  for(let i=0;i<3;i++)b[i]=boundary[i].times(d(1).minus(e)).plus(start[i].times(e));
 }else{
  b[path.index]=path.kind==='edge'?e:d(1).minus(e);
  const remaining=path.kind==='edge'?d(1).minus(e):e;
  b[other[0]]=remaining.times(t);b[other[1]]=remaining.times(d(1).minus(t));
 }
 return{epsilon:e,r:d(1).minus(b[2]).div(2),s:d(1).minus(b[0]).div(2),barycentric:b};
}
function approachCoordinates(path,epsilon,precisionDigits=70){
 const p=createPrecisionEngine(precisionDigits),point=coordinates(p,normalizePath(path),epsilon);
 return{epsilon:point.epsilon.toString(),r:point.r.toNumber(),s:point.s.toNumber(),
  rExact:point.r.toString(),sExact:point.s.toString(),barycentric:point.barycentric.map(v=>v.toString())};
}
function defaultApproachEpsilons(target='1e-6',samplesPerDecade=3){
 const last=new Decimal(target);if(!last.gt(0)||last.gt('.1')||last.lt('1e-10'))throw new RangeError('Choose a target epsilon from 0.1 through 1e-10.');
 const end=-Math.log10(last.toNumber()),out=[];
 for(let exponent=1;exponent<end-1e-12;exponent+=1/samplesPerDecade)
  out.push(Decimal.pow(10,-exponent).toPrecision(18));
 out.push(last.toString());return out;
}

/** Follow one geometrically labeled sheet to actual points approaching a face.
 * Epsilon is barycentric; the ordinary UI's source margin is deliberately not
 * used here. Returned coordinates remain in the fixed original F0 marking. */
async function solveApproachPath(counts,path,epsilons,options={}){
 const started=clock(),data=options.anchorData??storedAnchors;
 if(!data)throw new TypeError('Load the geometric attachment anchors before tracing an approach.');
 setApproachAnchorData(data);const checked=validateCounts(counts,data);
 if(checked.status!=='ok')throw new RangeError(checked.message);
 path=normalizePath(path);epsilons=epsilons??defaultApproachEpsilons(options.target??'1e-6');
 if(!Array.isArray(epsilons)||!epsilons.length||epsilons.length>100)throw new RangeError('Provide between one and 100 approach samples.');
 const precisionDigits=Math.max(60,Math.min(110,options.precisionDigits??70));
 const validationDigits=Math.max(precisionDigits+25,Math.min(160,options.validationDigits??100));
 const p=createPrecisionEngine(precisionDigits),v=createPrecisionEngine(validationDigits);
 const requested=epsilons.map(e=>p.d(String(e)));
 for(let i=0;i<requested.length;i++)if(!requested[i].gte('1e-10')||requested[i].gt(path.start?'1':'.3')||(i&&requested[i].gte(requested[i-1])))
  throw new RangeError(path.start?'Segment epsilons must strictly decrease from at most 1 to at least 1e-10.':'Epsilons must strictly decrease from at most 0.3 to at least 1e-10.');
 const key=counts.join(','),entry=data.anchors[key],baseEntry=data.anchors['0,0,0'];
 if(!entry||!baseEntry)throw new RangeError('This geometric sheet has no validated anchor.');
 const shift=[p.d(counts[0]+counts[2]).div(2),p.d(counts[1]+counts[2]).div(2)];
 let baseTau=p.C(...(baseEntry.highPrecisionTau??baseEntry.tau)),selectedTau=p.C(...(entry.highPrecisionTau??entry.tau));
 let currentR=p.d(String(data.r)),currentS=p.d(String(data.s)),currentEpsilon=null,steps=0,totalIterations=0;
 let stateLogEpsilon=null,previousLogEpsilon=null,previousBaseTau=null,previousSelectedTau=null;
 const nodes=[],cancelled=()=>Boolean(options.signal?.aborted),yieldUI=()=>new Promise(resolve=>setTimeout(resolve,0));
 async function correctAt(point,initial=false){
  let baseSeed=baseTau,selectedSeed=selectedTau;
  const nextLogEpsilon=point.epsilon?.ln()??null;
  if(!initial&&nextLogEpsilon&&stateLogEpsilon&&previousLogEpsilon&&!stateLogEpsilon.eq(previousLogEpsilon)){
   // A local secant predictor in log(epsilon) reduces Newton work near a
   // cusp. The accepted root still obeys the original continuation guard.
   const ratio=nextLogEpsilon.minus(stateLogEpsilon).div(stateLogEpsilon.minus(previousLogEpsilon));
   if(ratio.gt(0)&&ratio.lte(3)){
    const predict=(tau,previous)=>{
     let move=p.K(p.S(tau,previous),ratio),distance=p.abs(move);
     if(distance.gt('.12'))move=p.K(move,p.d('.12').div(distance));
     return p.A(tau,move);
    };
    baseSeed=predict(baseTau,previousBaseTau);selectedSeed=predict(selectedTau,previousSelectedTau);
   }
  }
  const base=p.correct(baseSeed,point.r,point.s,1,{tolerance:'1e-27',trustRadius:initial?.2:.3});
  if(base.status!=='ok')return{status:base.status,branch:'base'};
  const selected=checked.n===1?base:p.correct(selectedSeed,point.r.plus(shift[0]),point.s.plus(shift[1]),checked.n,{tolerance:'1e-27',trustRadius:initial?.2:.3});
  if(selected.status!=='ok')return{status:selected.status,branch:'result'};
  if(p.abs(p.S(base.tau,baseTau)).gt('.16')||p.abs(p.S(selected.tau,selectedTau)).gt('.16'))return{status:'step_too_large'};
  if(nextLogEpsilon){previousLogEpsilon=stateLogEpsilon;stateLogEpsilon=nextLogEpsilon;previousBaseTau=baseTau;previousSelectedTau=selectedTau;}
  baseTau=base.tau;selectedTau=selected.tau;currentR=point.r;currentS=point.s;
  totalIterations+=base.iterations+(checked.n===1?0:selected.iterations);steps++;await yieldUI();return{status:'ok'};
 }
 const first=coordinates(p,path,requested[0]);
 // The existing double engine supplies only an initial root on this exact
 // labeled sheet; multiprecision correction and comparison accept every node.
 const warm=solveAttachment(first.r.toNumber(),first.s.toNumber(),counts,null,data);
 if(warm.status==='ok'){
  baseTau=p.asComplex(warm.base.tau);selectedTau=p.asComplex(warm.result.tau);
  const result=await correctAt(first,true);if(result.status!=='ok')throw new Error('The high-precision initial correction failed: '+result.status);
 }else{
  const originR=currentR,originS=currentS,distance=Math.hypot(first.r.minus(originR).toNumber(),first.s.minus(originS).toNumber());
  let fraction=0,increment=Math.min(1,.012/Math.max(distance,.001));
  for(let attempts=0;fraction<1&&attempts<300;attempts++){
   if(cancelled())return{status:'aborted',counts,path,nodes,elapsedMs:clock()-started};
   const next=Math.min(1,fraction+increment),point={r:originR.plus(first.r.minus(originR).times(next)),s:originS.plus(first.s.minus(originS).times(next))};
   const result=await correctAt(point,true);
   if(result.status==='ok'){fraction=next;increment=Math.min(1-fraction,increment*1.3);}else{increment/=2;if(increment<1e-6)throw new Error('The geometric anchor could not reach the requested approach path.');}
  }
  if(fraction<1)throw new Error('The initial continuation exceeded its safe step budget.');
 }
 async function checkedBranch(tau,n,dr,ds){
  const point=coordinates(v,path,currentEpsilon.toString());
  const initial=v.C(tau[0].toString(),tau[1].toString());
  const corrected=v.correct(initial,point.r.plus(dr),point.s.plus(ds),n,{tolerance:'1e-38',trustRadius:.025,maxIterations:16});
  if(corrected.status!=='ok')return{status:'precision_limited',message:'The independent higher-precision correction did not converge.'};
  const error=v.abs(v.S(corrected.tau,initial));
  if(error.gt('1e-17'))return{status:'precision_limited',crossPrecisionTauError:error.toString(),message:'The two working precisions did not agree to the required accuracy.'};
  // An unchanged higher-precision iterate is not an exact root. Include its
  // remaining Newton correction so zero inter-precision difference cannot
  // create an inflated accuracy label. This remains a numerical estimate.
  const numericalErrorEstimate=v.D.max(error,corrected.correction,v.d(10).pow(-precisionDigits+15));
  const agreementDigits=Math.max(0,Math.min(precisionDigits-15,Math.floor(-numericalErrorEstimate.log(10).toNumber())));
  return{status:'ok',n,tau:{re:corrected.tau[0].toNumber(),im:corrected.tau[1].toNumber()},
   tauExact:corrected.tau.map(z=>z.toPrecision(Math.min(validationDigits,60))),precisionDigits,validationDigits,
   precisionAgreementDigits:agreementDigits,crossPrecisionTauError:error.toString(),numericalErrorEstimate:numericalErrorEstimate.toString(),
   correctionEstimate:corrected.correction.toString(),iterations:corrected.iterations,validation:'two_precisions_not_interval_certified'};
 }
 for(let index=0;index<requested.length;index++){
  if(cancelled())return{status:'aborted',counts:[...counts],path,nodes,elapsedMs:clock()-started};
  const target=requested[index];let failure=null;
  if(currentEpsilon!==null){
   let log=currentEpsilon.ln(),targetLog=target.ln(),increment=p.d('.8');
   for(let attempts=0;log.gt(targetLog)&&attempts<500;attempts++){
    if(cancelled())return{status:'aborted',counts:[...counts],path,nodes,elapsedMs:clock()-started};
    const next=p.D.max(targetLog,log.minus(increment)),point=coordinates(p,path,next.exp());
    const corrected=await correctAt(point);
    if(corrected.status==='ok'){log=next;increment=p.D.min('1.2',increment.times('1.25'));}
    else{increment=increment.div(2);if(increment.lt('1e-5')){failure=corrected;break;}}
   }
   if(log.gt(targetLog)&&!failure)failure={status:'continuation_limit'};
  }
  currentEpsilon=target;const point=coordinates(p,path,target);
  let base,result;
  if(failure){base={status:'unresolved'};result={status:failure.status,message:'Continuation stopped before this node; no other root was substituted.'};}
  else{
   base=await checkedBranch(baseTau,1,0,0);
   result=checked.n===1?{...base}:await checkedBranch(selectedTau,checked.n,shift[0].toString(),shift[1].toString());
  }
  const node={epsilon:target.toString(),r:point.r.toNumber(),s:point.s.toNumber(),rExact:point.r.toString(),sExact:point.s.toString(),
   barycentric:point.barycentric.map(z=>z.toString()),base,result,status:base.status==='ok'&&result.status==='ok'?'ok':'unresolved',
   elapsedMs:clock()-started,continuationSteps:steps};
  nodes.push(node);if(options.onProgress)await options.onProgress(node,index+1,requested.length);await yieldUI();
  if(node.status!=='ok')break;
 }
 return{status:nodes.length===requested.length&&nodes.every(n=>n.status==='ok')?'ok':'partial',counts:[...counts],n:checked.n,path,nodes,
  requestedEpsilons:requested.map(z=>z.toString()),completed:nodes.filter(n=>n.status==='ok').length,total:requested.length,
  stoppedAt:nodes.at(-1)?.status==='ok'?null:nodes.at(-1)?.epsilon,precisionDigits,validationDigits,elapsedMs:clock()-started,
  continuationSteps:steps,iterations:totalIterations,method:'geometric_label_continuation_and_two_decimal_precisions',
  epsilonDefinition:path.start?'Segment fraction: weights=(1-epsilon)*boundaryWeights+epsilon*startWeights; epsilon=1 is the selected start and epsilon=0 is the ideal boundary endpoint.':'Normalized barycentric face or corner distance; weights=(1-2s,2(r+s)-1,1-2r).'};
}

return {setApproachAnchorData,defaultApproachEpsilons,solveApproachPath,normalizePath,approachCoordinates};
})();
const attachmentAnchors = {"r":0.31,"s":0.29,"maxN":12,"anchors":{"0,0,0":{"tau":[0.5904953761444915,0.7963943351098978],"highPrecisionTau":["0.59049537614449157416284192577349226908012778620242941899071673872","0.79639433510989779706922649871403435841001942722694363801613000114"]},"0,0,1":{"tau":[0.7515095138734595,0.6524433021803057],"highPrecisionTau":["0.75150951387345949645384410891502963838365256442799368836445356092","0.65244330218030567273735558173074768316920664224957470615138134323"]},"0,1,0":{"tau":[0.35448685090813264,0.6879997451009816],"highPrecisionTau":["0.35448685090813261961279501711880500203196199222685656403084028255","0.6879997451009816549969185915679019939473208000296978983723724225"]},"1,0,0":{"tau":[0.5823785379137855,1.1341826530853316],"highPrecisionTau":["0.58237853791378546366575728619209822223639032030552707907671234284","1.1341826530853315962403696846543904016643758273288518054099512232"]},"0,0,2":{"tau":[0.8103860604950234,0.5797764171757713],"highPrecisionTau":["0.8103860604950234175731630008717259854202578066353195524099023101","0.57977641717577128277353126631812503383246140420934427464566736476"]},"0,1,1":{"tau":[0.5400832084805506,0.646138229333322],"highPrecisionTau":["0.54008320848055053686221181766791611682008256560754518625148399504","0.64613822933332195609918022468290574950003947134280503981826771478"]},"0,2,0":{"tau":[0.2661530927667447,0.6170110529155463],"highPrecisionTau":["0.26615309276674471073268400652330747298467210193293483517054268986","0.61701105291554622938162772230971657851016530026401569330520689838"]},"1,0,1":{"tau":[0.7554851567438122,0.9044222083961427],"highPrecisionTau":["0.75548515674381221284180199628711477404486436987017071537132345048","0.90442220839614275355237867598541062254703419291474191429440448789"]},"1,1,0":{"tau":[0.33641698617555044,0.9361692214053354],"highPrecisionTau":["0.33641698617555046406636948191940451238830501503669948991252830308","0.93616922140533537502706857312216666531397406482493637466547587488"]},"2,0,0":{"tau":[0.5802309037096132,1.3519173892371166],"highPrecisionTau":["0.58023090370961317901706581078627485590837001710168383595215299674","1.3519173892371164981016024007805815496337166724563258239745718894"]},"0,0,3":{"tau":[0.841728094023958,0.5344069082972436],"highPrecisionTau":["0.84172809402395794437798963984200296198278994039341481977615430925","0.53440690829724356011926677948157995589544993666502620691599366763"]},"0,1,2":{"tau":[0.6299910352020642,0.6065667355063957],"highPrecisionTau":["0.62999103520206424666514873154487720771148631318073206502307391102","0.6065667355063956826424102709933489384973360426935705179468808005"]},"0,2,1":{"tau":[0.4336528965362644,0.6170485958998836],"highPrecisionTau":["0.43365289653626434940007447921084977524877221510063418025648297513","0.61704859589988364076591207303555004368890260561499544131519013268"]},"0,3,0":{"tau":[0.2194990271876052,0.5699249827074794],"highPrecisionTau":["0.21949902718760520093099855513898729088449236729422662883266496726","0.56992498270747942086610151461424252140108037271255773478539646094"]},"1,0,2":{"tau":[0.8190763243004178,0.788631148844644],"highPrecisionTau":["0.81907632430041772892328519180239142067540837434400890714722488429","0.78863114884464407776975156091075380608999434055926951769016416544"]},"1,1,1":{"tau":[0.5345911839088004,0.8413672802455601],"highPrecisionTau":["0.53459118390880048696444534995708777020751403305227293502385883776","0.84136728024556007641866915000768572689992943651461139960498303591"]},"1,2,0":{"tau":[0.24508363041652606,0.820541713936849],"highPrecisionTau":["0.24508363041652606668447818507255165119622906630612429387480143497","0.82054171393684906578724201430520828588042921667795591595711869001"]},"2,0,1":{"tau":[0.7566877827649808,1.0780300107149037],"highPrecisionTau":["0.75668778276498081344247991255363812083684005946317176155617978348","1.0780300107149035645178205761694741260449757878805112991332941373"]},"2,1,0":{"tau":[0.330728071930688,1.1089889089128258],"highPrecisionTau":["0.33072807193068796164601677265754057598583464328590510218069029164","1.1089889089128256834164626285410938482130621636056091194820439256"]},"3,0,0":{"tau":[0.5793596671311855,1.5133680043418907],"highPrecisionTau":["0.57935966713118549646785059370440192075542422181054615870192472471","1.5133680043418905691687554543054482224506433570565979892916896397"]},"0,0,4":{"tau":[0.8615254401911367,0.5026248506350188],"highPrecisionTau":["0.86152544019113679176238261036090362791599906106155313418191988599","0.50262485063501883243903202792683037354349133615644838533614111682"]},"0,1,3":{"tau":[0.6842316944717046,0.575021459683328],"highPrecisionTau":["0.68423169447170460977163152039873638600324844760931760601219117439","0.57502145968332792830429631206393700577703847815327820997756480831"]},"0,2,2":{"tau":[0.5267688263467021,0.6000381495242691],"highPrecisionTau":["0.52676882634670203761321702168536353822548443379150011373399405036","0.60003814952426916341894021800570857749777995150537650109461929878"]},"0,3,1":{"tau":[0.36881306450740986,0.5894886050882322],"highPrecisionTau":["0.36881306450740986600019944766522646292726181838047444235391206642","0.58948860508823214720729626181194754655378803182106548012591917938"]},"0,4,0":{"tau":[0.19030803642972227,0.5361316608784794],"highPrecisionTau":["0.19030803642972226195923039906606751094429657648058034784256067124","0.53613166087847942481143944943590013202057526501703482121631460924"]},"1,0,3":{"tau":[0.8526909974897778,0.7164178524340136],"highPrecisionTau":["0.85269099748977779302243528157438245855071283373839515788272474449","0.71641785243401359520687806770305032026740422441408457852575423185"]},"1,1,2":{"tau":[0.6314478387366457,0.7723654650396221],"highPrecisionTau":["0.63144783873664569897096624501800919017731117667232515026955020756","0.7723654650396220668988983506195457141834741745512073129968633332"]},"1,2,1":{"tau":[0.42284498422152705,0.7812601291251633],"highPrecisionTau":["0.42284498422152703173937940212605588027360784149563578338332795252","0.78126012912516329735344916647031757865170924841893912207017006811"]},"1,3,0":{"tau":[0.1973984260185703,0.7458993031656952],"highPrecisionTau":["0.19739842601857030394926260489959434252067563703764074971139778595","0.74589930316569522719944247157262638098458789659112178159557837239"]},"2,0,2":{"tau":[0.8220118840367037,0.9366554849150815],"highPrecisionTau":["0.82201188403670369127734153917229160521180002468082310926184341498","0.93665548491508158237731357882655183193651114764700543756265658965"]},"2,1,1":{"tau":[0.5326043104373224,0.9845389328712051],"highPrecisionTau":["0.5326043104373223345175149757944280437126504808483973725883114957","0.98453893287120513983018885997343556056975419330537722734883394418"]},"2,2,0":{"tau":[0.23787172685968883,0.9673218781297337],"highPrecisionTau":["0.23787172685968881552678510133220799921056458967044165201797200653","0.96732187812973369267972860672174032570951911869253887776187155073"]},"3,0,1":{"tau":[0.7572248859548885,1.2123248746321915],"highPrecisionTau":["0.75722488595488845319906528414650769917798664631152593431883823719","1.2123248746321916028050942702024680304975404206377622726359490042"]},"3,1,0":{"tau":[0.32814204362632093,1.2429897592333141],"highPrecisionTau":["0.32814204362632094045321879972854872701078605734163594623952594241","1.2429897592333140928882380021799093100987073416475871014131038247"]},"4,0,0":{"tau":[0.5789215292300838,1.64185870386297],"highPrecisionTau":["0.57892152923008382347305331105318138443897892382885769774757269725","1.6418587038629699722056853929650416131181086950589412496293698848"]},"0,0,5":{"tau":[0.87532736448948,0.47873158903313423],"highPrecisionTau":["0.87532736448947998989392931131756188132882542758413995587904279248","0.47873158903313425512431564199566686010994468294081736618363131118"]},"0,1,4":{"tau":[0.7209441947002367,0.5497539960286885],"highPrecisionTau":["0.72094419470023678379441655762213503456221642166883993470250605756","0.5497539960286884913533071351943236046348998525070588908654107854"]},"0,2,3":{"tau":[0.587659408841123,0.5813159221898846],"highPrecisionTau":["0.58765940884112298984260853113034659798112206512704164358904334387","0.58131592218988460152940071150247813259238731326111584413191558218"]},"0,3,2":{"tau":[0.45853298654130986,0.5864119195544156],"highPrecisionTau":["0.45853298654130985413521440762088757942486287528760971471361864912","0.58641191955441559775222197743266137078719362107241107908411946493"]},"0,4,1":{"tau":[0.32481896292121665,0.5658832421781945],"highPrecisionTau":["0.32481896292121666360821188628315927742123723995897777835699154832","0.56588324217819451154632017772738250071382839776731011019832529447"]},"0,5,0":{"tau":[0.17012729074846442,0.5104168701229547],"highPrecisionTau":["0.17012729074846442719830977101281023773357590692605419970522340727","0.51041687012295462879895022007481363011995661392148045537708790572"]},"1,0,4":{"tau":[0.8737146117761868,0.6659979038828044],"highPrecisionTau":["0.87371461177618683686792676045109496013823980693340644678345352783","0.66599790388280449506727724738527214440301055670187071001000178766"]},"1,1,3":{"tau":[0.6897768161485427,0.7213488596924934],"highPrecisionTau":["0.68977681614854273670238952302105090526885629109728594642841292682","0.72134885969249345991885241542859087965666040248495668524284454564"]},"1,2,2":{"tau":[0.5228133583527507,0.7416723275535352],"highPrecisionTau":["0.52281335835275071335372516715059654417468751928851430680670848517","0.74167232755353515308656362678206612559167239852240288279847073706"]},"1,3,1":{"tau":[0.3550011865750481,0.7333908261815947],"highPrecisionTau":["0.35500118657504811105054704307691854997202539677053553183165233342","0.73339082618159470258805910817218290282473827443673705600063241579"]},"1,4,0":{"tau":[0.16791151078300529,0.6931050520647739],"highPrecisionTau":["0.16791151078300528008349169871694909977370254364645526411680998573","0.69310505206477390360081980771087257096288705206867300177110263357"]},"2,0,3":{"tau":[0.856600950061261,0.8473150080377062],"highPrecisionTau":["0.85660095006126107603085759331448131591926957388724166238879057861","0.84731500803770616682614732436409075459340519026671566753033440705"]},"2,1,2":{"tau":[0.632157377106861,0.8969694583022713],"highPrecisionTau":["0.63215737710686109457554655443471075462615074021461477502987149979","0.89696945830227129525532275557527030842890861675509371485537132268"]},"2,2,1":{"tau":[0.41854260474439214,0.905436867652474],"highPrecisionTau":["0.41854260474439212256403006588805169214482855376088308305384923866","0.9054368676524740437391142374833323370833201913735984967688947774"]},"2,3,0":{"tau":[0.18945832036306007,0.875279890371322],"highPrecisionTau":["0.18945832036306005460419402300552539794039105226516813325044573107","0.87527989037132206200438493392664737788090466683572871410143389135"]},"3,0,2":{"tau":[0.8234166474721952,1.0536888952271188],"highPrecisionTau":["0.82341664747219527044209338600618231731810406559927273990913608647","1.0536888952271188137532989429423864912752877788734776048425626928"]},"3,1,1":{"tau":[0.5316145465016989,1.0991956202308],"highPrecisionTau":["0.53161454650169893714610633865980810574296630645367251874081216952","1.0991956202308000481404997115596830952381833549336505136119099664"]},"3,2,0":{"tau":[0.23439333176450308,1.0838565533760522],"highPrecisionTau":["0.23439333176450307736411955181909579165510046821841070533424084074","1.0838565533760521676096512597819470030794662973820065995755979604"]},"4,0,1":{"tau":[0.7575144695144952,1.3223411987468239],"highPrecisionTau":["0.75751446951449524024547708864000577406543475566745071038382009349","1.3223411987468237660502454926151864447380280705288792306802344496"]},"4,1,0":{"tau":[0.3267340104740559,1.3528628092637536],"highPrecisionTau":["0.32673401047405592949129418813456998877711332448204661597705588541","1.3528628092637535868084849863346567031313176593542605399277377056"]},"5,0,0":{"tau":[0.5786704237611745,1.7486344886564826],"highPrecisionTau":["0.57867042376117452490104457756525778178840133017895355965498943366","1.7486344886564825628276727408487155710443017132425650661551689344"]},"0,0,6":{"tau":[0.8855883749676,0.459893337978774],"highPrecisionTau":["0.88558837496760001689701039401265396872705556184812745947278918631","0.45989333797877403982308818181416148616913797206418153515844414766"]},"0,1,5":{"tau":[0.7476547390629024,0.5290741282185445],"highPrecisionTau":["0.74765473906290235728704380927307720678498484041861837770577659765","0.52907412821854450758465974785199194628435627206421227162885795137"]},"0,2,4":{"tau":[0.6310714701953587,0.5638367265401152],"highPrecisionTau":["0.63107147019535872176236397873568064523026928063100004154683640219","0.56383672654011517453895440387700648674074333888088718032256425625"]},"0,3,3":{"tau":[0.5204211943977096,0.5770842731764217],"highPrecisionTau":["0.52042119439770957192456156023878881487418254556011419283712549052","0.57708427317642170967626542290450102697494494132529032617857613792"]},"0,4,2":{"tau":[0.40965039211323434,0.5716922698326926],"highPrecisionTau":["0.40965039211323433337323945301435064843358891821475307217089056163","0.57169226983269256245385693415069817794105394707859884186279865714"]},"0,5,1":{"tau":[0.29281381542373114,0.5458775986357257],"highPrecisionTau":["0.2928138154237311125115774733454388486464267988743328631177230274","0.54587759863572566722113689732651636836904989182044236774066676159"]},"0,6,0":{"tau":[0.1552308493390186,0.49000353533730145],"highPrecisionTau":["0.15523084933901860067231531816440399926991551480064015099819265645","0.49000353533730144852416362697322401358068339166426027403075242072"]},"1,0,5":{"tau":[0.8882196634076015,0.6282589367010324],"highPrecisionTau":["0.88821966340760146444024861666919190735157592365813410524143767743","0.62825893670103244082411810756964209488335829231643078625987129953"]},"1,1,4":{"tau":[0.7290769007922019,0.6820042681833915],"highPrecisionTau":["0.72907690079220191578437785402494521803752743126677680463838969438","0.68200426818339153710338327159037902508138632118593445197064685026"]},"1,2,3":{"tau":[0.5881532493915438,0.7075651684523052],"highPrecisionTau":["0.58815324939154380106576707309851677133732149528119124291864145041","0.70756516845230517135711282291658655520618390872770731883257701078"]},"1,3,2":{"tau":[0.4510576512137605,0.7117896923823168],"highPrecisionTau":["0.45105765121376054069923593773251282600855024608173282334640725286","0.71178969238231680684338984336637680933874415010923337820968652717"]},"1,4,1":{"tau":[0.3092118506224589,0.695216006203829],"highPrecisionTau":["0.30921185062245885994015969729226740965302035887882430105729306513","0.69521600620382894610330745140008141569485121516791446485127257225"]},"1,5,0":{"tau":[0.1477519141514894,0.6533558004801285],"highPrecisionTau":["0.14775191415148940214730895203560951569170286495389917395531618832","0.65335580048012852092837195860983012506322195276232498978091689156"]},"2,0,4":{"tau":[0.8782133522035309,0.7844578175952746],"highPrecisionTau":["0.87821335220353094042726294393664672553156314471002550120013074007","0.78445781759527455599119007985660224061405900219129400475043842081"]},"2,1,3":{"tau":[0.6922991189465529,0.8329147799736228],"highPrecisionTau":["0.69229911894655296676861149736159614062421720434183465611635121365","0.83291477997362284981357378329660109760235747885753783710889870637"]},"2,2,2":{"tau":[0.5212230759002172,0.8511681141462283],"highPrecisionTau":["0.5212230759002172614306987826350007892286443672814369554169548737","0.85116811414622827757140159012376847577678375596443688605848424793"]},"2,3,1":{"tau":[0.34917986947047064,0.8442598735593233],"highPrecisionTau":["0.34917986947047062637633555303720313451862408185632408162946329534","0.84425987355932335021435375726319731171474916828709971677314908112"]},"2,4,0":{"tau":[0.15960719886127595,0.8098782385880396],"highPrecisionTau":["0.1596071988612759526338024768523755539977552356625024048095187703","0.80987823858803967823550265368919752342431644661615530130239384758"]},"3,0,3":{"tau":[0.8585464147633195,0.9521664394926173],"highPrecisionTau":["0.85854641476331944673754273460655262786745359767357990349315470301","0.95216643949261724572267354663318922313081007923077376617509340288"]},"3,1,2":{"tau":[0.6325417195306986,0.9985719839741891],"highPrecisionTau":["0.6325417195306986182747090801291883629494670865767526355407651868","0.99857198397418915262142482918692510875520880015575339466021101146"]},"3,2,1":{"tau":[0.4162831369003377,1.0068608466821334],"highPrecisionTau":["0.41628313690033770660882007957803817225646066870984727250026943391","1.006860846682133340413689058595090881435203678968653282849719714"]},"3,3,0":{"tau":[0.18548711432422435,0.979494178808416],"highPrecisionTau":["0.18548711432422435111335643795778531309896435686984842787575145588","0.9794941788084160989640805981540994458511397355180929070281671434"]},"4,0,2":{"tau":[0.8242124660453741,1.151180299071085],"highPrecisionTau":["0.82421246604537416282967221432402634619356765568250688524567187383","1.1511802990710850586240856728966466896535347463643848840578202466"]},"4,1,1":{"tau":[0.5310399425071837,1.1953234039601657],"highPrecisionTau":["0.53103994250718363680378579454853962668277413214754692536867270272","1.1953234039601656143849404359854521159840082767998228545372767591"]},"4,2,0":{"tau":[0.2324122394543271,1.1810929728116042],"highPrecisionTau":["0.23241223945432709546839307043893110858611149687303660229889116872","1.1810929728116041574948426309424751714771273342904699809303568356"]},"5,0,1":{"tau":[0.7576893917453332,1.415707662218239],"highPrecisionTau":["0.757689391745333247528234198525341799396657791035095473555384135","1.4157076622182388621098222179203795658677349003691051237424452893"]},"5,1,0":{"tau":[0.32587834922456627,1.4461479825119008],"highPrecisionTau":["0.32587834922456628166152581375478867700790983474920457069435268965","1.4461479825119008242495525393020983405971588814043943620480785326"]},"6,0,0":{"tau":[0.578513211933155,1.8400026346529397],"highPrecisionTau":["0.57851321193315493012856392864328060551480893868003335606331070938","1.8400026346529395740690358325599088985257271554023491653570396556"]},"0,0,7":{"tau":[0.8935695170631885,0.4445238550546686],"highPrecisionTau":["0.89356951706318849075508769891489852966763893799898237342261918614","0.44452385505466860120741290521925276103626487376772335258913648391"]},"0,1,6":{"tau":[0.7680816619446168,0.5117859055568753],"highPrecisionTau":["0.76808166194461680712201859559477977622690524208049922407970224875","0.51178590555687528760780020317802243616545967509102477315345214022"]},"0,2,5":{"tau":[0.6638087929376337,0.5481312045758051],"highPrecisionTau":["0.66380879293763374423740883452960200491604101268363617298727852572","0.54813120457580514691695397091290671917269913644764008484439164416"]},"0,3,4":{"tau":[0.5662502583558813,0.5660674216966612],"highPrecisionTau":["0.56625025835588122226192078988337018277252002008060330619115071576","0.56606742169666121806160811835941575812501913284662678530877988483"]},"0,4,3":{"tau":[0.47038431580611356,0.5691191152884885],"highPrecisionTau":["0.47038431580611354786360872914299885857360822710116304628900296598","0.56911911528848842356759777786567364167955659435964057048289114733"]},"0,5,2":{"tau":[0.37270851191302595,0.5575773347273169],"highPrecisionTau":["0.37270851191302597714723449691379776878768339656970990418084405079","0.55757733472731683212656249973699847580085809461235626080737651275"]},"0,6,1":{"tau":[0.2683647278635283,0.52879748719793],"highPrecisionTau":["0.26836472786352825943652763604182126894324183389907044806611236917","0.52879748719793000028373366987019175669553487877778934863617084826"]},"0,7,0":{"tau":[0.14371472207625965,0.4732798229549187],"highPrecisionTau":["0.14371472207625965452142122890917411932368895878540360538890123715","0.47327982295491869472361722979403071729528778884555460327658048054"]},"1,0,6":{"tau":[0.8988946785558809,0.5986474112765273],"highPrecisionTau":["0.89889467855588088476482507030608875844906881159172819300197611728","0.59864741127652726744888854782277050263463337549008763442081842092"]},"1,1,5":{"tau":[0.757515289092806,0.6505645672211018],"highPrecisionTau":["0.75751528909280597565043295356054253809747916470560228477010955912","0.65056456722110181225943745550001991687361517005313930485920498805"]},"1,2,4":{"tau":[0.6345895968123779,0.6786954967391343],"highPrecisionTau":["0.63458959681237787124775009447253622250312407612340239694020994685","0.67869549673913432174000036648635359401118017554777036927023171463"]},"1,3,3":{"tau":[0.5173362701255009,0.6896700070199766],"highPrecisionTau":["0.5173362701255008789636984007982471794580046293471055590681803961","0.68967000701997655653283838374336209148302550892942842229928132056"]},"1,4,2":{"tau":[0.39982624715111464,0.6851433733448178],"highPrecisionTau":["0.39982624715111466252482764020428321786024951551736380761213301904","0.68514337334481773526812707435432375182190702429812505495888157146"]},"1,5,1":{"tau":[0.276090469234877,0.6641428959307979],"highPrecisionTau":["0.27609046923487702449296912278986955822912752301445163465783154524","0.66414289593079784394398336748013410864438610897279802540444625927"]},"1,6,0":{"tau":[0.13302346902452356,0.6220753749174971],"highPrecisionTau":["0.1330234690245235593649361691312010034161408487946490611721650217","0.62207537491749704169103662037455469005599701762019976070917570568"]},"2,0,5":{"tau":[0.8930918423482282,0.7371942141027298],"highPrecisionTau":["0.89309184234822820384087079989037466992686501861297942938348765661","0.73719421410272978828101180246342420590830176787854697483933198902"]},"2,1,4":{"tau":[0.7328497447335313,0.7837854959610272],"highPrecisionTau":["0.73284974473353128358311006140309555124211936797982910951047604633","0.78378549596102716748266199763819421380300956444152877203360952113"]},"2,2,3":{"tau":[0.5885851008708306,0.8064731758072106],"highPrecisionTau":["0.58858510087083060333941607072529884486673386606847125504384227563","0.80647317580721053851314565755508170731956816583397503777328240864"]},"2,3,2":{"tau":[0.447737001043542,0.8104257477661121],"highPrecisionTau":["0.44773700104354204935935945504578860061428127216814639730175223981","0.81042574776611212830034045503909405378152042702940772164159040232"]},"2,4,1":{"tau":[0.3023745099048237,0.7961194599484742],"highPrecisionTau":["0.30237450990482365453672080171501620847786451082741108020321925467","0.79611945994847424258042172955922167626857329996595125285483857958"]},"2,5,0":{"tau":[0.13926952084402813,0.7604944279863987],"highPrecisionTau":["0.1392695208440281220643756277114309463142828721138326913463173153","0.76049442798639864973560588342119045920129937127135056002140654301"]},"3,0,4":{"tau":[0.8805111459535212,0.88014891120271],"highPrecisionTau":["0.88051114595352122095271303893219459268701620479921293821688917677","0.88014891120270998576037563733841175149821439615730837210067285841"]},"3,1,3":{"tau":[0.6936929117459214,0.9249420837481607],"highPrecisionTau":["0.69369291174592142370358703435042566341084249106166355422090702999","0.92494208374816068509230694985157704983128795812717153649895360852"]},"3,2,2":{"tau":[0.5203604410850786,0.9420472610477201],"highPrecisionTau":["0.52036044108507859103368682415975442196995487464893298896143931498","0.94204726104772003961591098857670758612623804713004590225271851619"]},"3,3,1":{"tau":[0.3460120862317194,0.9359838200637809],"highPrecisionTau":["0.34601208623171941889594885436728463319572829310694634678387610005","0.93598382006378090292257396590979919261184353144775735791927728069"]},"3,4,0":{"tau":[0.1553506593979241,0.9048369846987677],"highPrecisionTau":["0.15535065939792410344823686283246145703605557170120958895039749881","0.90483698469876767247081923853466180748106878882758301381509290058"]},"4,0,3":{"tau":[0.8596823845481332,1.040458892557427],"highPrecisionTau":["0.85968238454813314310993681891020962542934931004995104211142327842","1.040458892557426918942268055991296300689459012675318158110228267"]},"4,1,2":{"tau":[0.6327740811462541,1.0849285948795253],"highPrecisionTau":["0.63277408114625406112482043900783098950217704768036632484737604631","1.0849285948795252130086891797011699341213425384076514870570773624"]},"4,2,1":{"tau":[0.41492179738690005,1.0931249845327629],"highPrecisionTau":["0.41492179738690006817033915385071066717175334835301144470402042493","1.0931249845327628327127570471298683493023793817481643128156505858"]},"4,3,0":{"tau":[0.18315877131009872,1.0674498166149522],"highPrecisionTau":["0.1831587713100987090648297757547107449872057702355942301767512777","1.0674498166149522614176492249613412522448428854832789315950362611"]},"5,0,2":{"tau":[0.8247117287312125,1.2350094503071851],"highPrecisionTau":["0.82471172873121243404698884179436968293657667625424122986196646102","1.2350094503071850913706625684913210728359604519277619270930447363"]},"5,1,1":{"tau":[0.5306736884551937,1.2782904040218237],"highPrecisionTau":["0.53067368845519374423756622828343644977136066014540875220950512178","1.2782904040218237009963416085671048823809411191519655999477783855"]},"5,2,0":{"tau":[0.23116478113594885,1.264772313747665],"highPrecisionTau":["0.23116478113594885143927173333528984842476679352740649306009519574","1.2647723137476648760834420769044670699633605902752187034493342103"]},"6,0,1":{"tau":[0.7578034945743788,1.4968886398046868],"highPrecisionTau":["0.75780349457437882579292464108728269126231236265293503927381062953","1.4968886398046867075347062265150151550308777611294659576584592251"]},"6,1,0":{"tau":[0.3253179805035468,1.5272780873967173],"highPrecisionTau":["0.32531798050354681179920532379330847929799973615650668238463649851","1.5272780873967174303502841641138728735512841220412202572335486014"]},"7,0,0":{"tau":[0.5784082889971616,1.9198626778550154],"highPrecisionTau":["0.57840828899716160506181526297230049942770937340638262617192138181","1.9198626778550155068119237124633891221957876946091110865794277259"]},"0,0,8":{"tau":[0.8999886735985638,0.43165749057710767],"highPrecisionTau":["0.8999886735985638164935848730960737450247202396415483724266427115768604283764402816197","0.4316574905771076544906715328770268110639221134798922035891507255165463507430055272952"]},"0,1,7":{"tau":[0.7842843731758058,0.4970654360110489],"highPrecisionTau":["0.7842843731758058397309777888971647632761796975825942988830995551360814827084533385268","0.4970654360110489107820601546828928928563490804352351244219463002985153939018937886154"]},"0,2,6":{"tau":[0.6895018910795871,0.5341374020413163],"highPrecisionTau":["0.6895018910795871753760376912251855668524255078764198502779889149271066593023441955548","0.5341374020413163427617500712603758865632680914192855574974659526749740424487595705793"]},"0,3,5":{"tau":[0.6017919435613742,0.5549293341084078],"highPrecisionTau":["0.6017919435613742434020853736825321015854385747331274463028246695118022883579667339876","0.5549293341084078254887743623922555341929209879360764449086098403692562753231974765941"]},"0,4,4":{"tau":[0.5166510006237209,0.5631702795678296],"highPrecisionTau":["0.516651000623720883116475998713225599500741782219765151749194361921473673608631603524","0.5631702795678295673957416145579126596000363549551665845620286538133578210173621235154"]},"0,5,3":{"tau":[0.43147276876014107,0.5599209902678184],"highPrecisionTau":["0.4314727687601410501564668413641215256018386829718760328072355605671341499577179538788","0.5599209902678183605394205106091168376941611002634460400772134261481596244303357158151"]},"0,6,2":{"tau":[0.34369136170007636,0.5445324732829234],"highPrecisionTau":["0.3436913617000763447889912360811007432833242900553914487043378409227138580589179569249","0.544532473282923393526382737406269702929944984532267787048661614332148386939806384443"]},"0,7,1":{"tau":[0.24900061019025993,0.5140511474271879],"highPrecisionTau":["0.2490006101902599309323841379787908470778136960958963763044606387550853240204880717713","0.5140511474271879565954863614195891279316968821113849215669962668096775906539479356999"]},"0,8,0":{"tau":[0.1345010488513675,0.4592428537479833],"highPrecisionTau":["0.1345010488513674747064786752784890636372475293797493130126073449753058890730886863367","0.4592428537479832671404288337348455128546409697258791351046367885587198212544396588371"]},"1,0,7":{"tau":[0.9071185345899363,0.5746069889852019],"highPrecisionTau":["0.9071185345899363278372097646264172857443786689534955120386637716828501262636686437669","0.5746069889852018900102985434864218834269072382808042934493071311347777359092138986597"]},"1,1,6":{"tau":[0.7791416217273306,0.6247271194210997],"highPrecisionTau":["0.7791416217273306223622639562432323489698525995328403137073075629181197203253756367587","0.62472711942109970928824380957934603279079541275540804549001984800031191345361659353"]},"1,2,5":{"tau":[0.6694686003680538,0.6541175224491951],"highPrecisionTau":["0.6694686003680537390665232793364906056666339918438252431210586334030552931605840825071","0.6541175224491951145640216557642679196714271957653673390161109912225795352751794158835"]},"1,3,4":{"tau":[0.566291131542106,0.6690524526333079],"highPrecisionTau":["0.5662911315421059204394568817858508875384426412826816837973955712903484445083638077717","0.6690524526333079566909013323927705271924132196031198802999164839146771740778318689242"]},"1,4,3":{"tau":[0.464760493774037,0.6715543759066558],"highPrecisionTau":["0.4647604937740370227584414812700978795910449557712499174293962044279250807578501219229","0.6715543759066557788772843701403794451902709308749971146448698837122285558235126842945"]},"1,5,2":{"tau":[0.36126522625152746,0.6618033415718381],"highPrecisionTau":["0.3612652262515274668268158481593664780905506693431912865663737729555101013166657697171","0.6618033415718380212886265194553006401416139769100329039092307002449810382162505803362"]},"1,6,1":{"tau":[0.25093182248719836,0.6383140187936306],"highPrecisionTau":["0.2509318224871983618895031803697801432482922624221840048068551628779876080225521168318","0.6383140187936306350227618950471425062991285846748631712272145209047033796041581332055"]},"1,7,0":{"tau":[0.1217443668562202,0.5966418365607582],"highPrecisionTau":["0.1217443668562201924831472807841784691460039856086937993125849550515639824414951148891","0.5966418365607581937428214995963749352026852257197013559116158636284099306154131443286"]},"2,0,6":{"tau":[0.9040110611158669,0.7000080185294845],"highPrecisionTau":["0.9040110611158669192212899767339517092859912569672296711930418130194009771378797468059","0.7000080185294845207709644371088323823310194951350510459115530919567347543198532930885"]},"2,1,5":{"tau":[0.762178227429426,0.7446696257600844],"highPrecisionTau":["0.762178227429426008334650355692060459044584055421072031145214488255220027521091390815","0.744669625760084408492485931938331135369286863075035741130579751277945534181254931658"]},"2,2,4":{"tau":[0.6365164099219819,0.7694924688840932],"highPrecisionTau":["0.6365164099219819233115329832823904349680124416273358021939568759436547560854796837576","0.7694924688840931250272809559336781813396628888497276728559472137995345283082370502598"]},"2,3,3":{"tau":[0.5160314298136841,0.7793499753523696],"highPrecisionTau":["0.5160314298136840773287430455102130483058568777325736402203904188413316682973284364644","0.7793499753523696120236716195764260907048640834395330090223451491697699167631740106525"]},"2,4,2":{"tau":[0.3952342014435136,0.7754820407072519],"highPrecisionTau":["0.3952342014435136105768065107510190187285360266268170407099615843419028627271402523893","0.7754820407072519053336140151413559192355718141497807549592089822457884095823398926408"]},"2,5,1":{"tau":[0.2685551858462416,0.7572438938373648],"highPrecisionTau":["0.2685551858462416166091166510379829149774321566118218860283119901997337751780057639078","0.7572438938373648147238626332796820571442177928818211650736515791135103186416589479297"]},"2,6,0":{"tau":[0.12446556844497446,0.7215668307205353],"highPrecisionTau":["0.1244655684449744567335898344554386933336740338018613940627148813701859649233727997165","0.7215668307205353291108006898233156799601461323776097474733411934294358772204737073322"]},"3,0,5":{"tau":[0.8956278040791181,0.8256929702307054],"highPrecisionTau":["0.8956278040791181766362410671082806693878148925799627692802619835614552793284561300064","0.8256929702307054407364489756982487154259369643879653123157769864986685012231286177811"]},"3,1,4":{"tau":[0.734983197185718,0.8684124065360601],"highPrecisionTau":["0.734983197185717977415741206254667155314129814427745195094410575275568613134058207791","0.8684124065360600896404854753483447537996260339470407575185347424660983095523652841342"]},"3,2,3":{"tau":[0.5888624815053458,0.8894315983125628],"highPrecisionTau":["0.5888624815053457812234529859582162121033141947974769120714694195576542531956609343246","0.8894315983125628284818622787466645824845646961184169636893083724284428517030951880484"]},"3,3,2":{"tau":[0.4458530141184615,0.8932607935389507],"highPrecisionTau":["0.4458530141184614948909730378671069293379104423876682264704248273395439047754475125425","0.893260793538950667945861595574269455057507511931454750488622554087194283182796446784"]},"3,4,1":{"tau":[0.2985581051723347,0.8803522880710194],"highPrecisionTau":["0.298558105172334699699348692838903794703226393210214394727204514924600933580655644289","0.880352288071019408639009384849866734341486675579347529845509799340532820017056554082"]},"3,5,0":{"tau":[0.13484413406375953,0.8481938901081237],"highPrecisionTau":["0.1348441340637595294487988027615264452300053552253394968450590649594523588427045319184","0.8481938901081236866354956739558237308734144436659394048455927806863835978804137603003"]},"4,0,4":{"tau":[0.8818814555448774,0.9613304545714055],"highPrecisionTau":["0.8818814555448774303777166962619584899489866777458802989156559383799143135179583120558","0.9613304545714055660604808543295819702836673929048036204792860925256891658809841513811"]},"4,1,3":{"tau":[0.694557975951264,1.0038902962802296],"highPrecisionTau":["0.6945579759512639533249365412249316318480445828129750872483750005482214033421956674357","1.003890296280229667582218742599672925347417775283102450071853425599964355055305625045"]},"4,2,2":{"tau":[0.5198248176004105,1.0202785499707214],"highPrecisionTau":["0.519824817600410494676432329601875250765541030744704431707671613941568899641873656566","1.020278549970721335552817636423128537116831120156337869872637584961899666910906207613"]},"4,3,1":{"tau":[0.3440525661399386,1.0147704677679885],"highPrecisionTau":["0.3440525661399386179504371460083575879998307370204140775545541480487990090823671661438","1.014770467767988566836792152646826632321566933477638846470435210795667768780335317458"]},"4,4,0":{"tau":[0.15280393431861738,0.9856218132307234],"highPrecisionTau":["0.1528039343186173880912209608097533022058734750092344218853444296829475871139170608958","0.9856218132307233698969403434951282830923849871186270759267338520071259966171781428124"]},"5,0,3":{"tau":[0.8604126867429136,1.1170602611367253],"highPrecisionTau":["0.8604126867429135954201254710346728968823570240960828252707248194540525927398456111579","1.117060261136725297497179952935359299579514849995811442334319078374838728592265070999"]},"5,1,2":{"tau":[0.632926040787829,1.1602697414617],"highPrecisionTau":["0.6329260407878289625512173544495867123425845453657077828799821023217010887650452831027","1.160269741461699995802906186250831789887795939748422852135522420300440656647329447624"]},"5,2,1":{"tau":[0.41402886114936516,1.1684116115104175],"highPrecisionTau":["0.4140288611493651403721220929424816174758680347367932385744845296641902121203350427462","1.168411611510417498322966742775506754480472624882525594052032778721043914882319927641"]},"5,3,0":{"tau":[0.18165721310513627,1.1438482089549935],"highPrecisionTau":["0.1816572131051362693410927088861038310729204404284640215573081050890424425080507494477","1.143848208954993511769399155950813139119190596055360096171194247232730072972953191538"]},"6,0,2":{"tau":[0.8250473662854825,1.308670545276684],"highPrecisionTau":["0.8250473662854824877592757110301806364515649933252744358139354774599049200077009604375","1.308670545276683954951880192856140977382015735287680324880743042005027259195551576678"]},"6,1,1":{"tau":[0.5304247761233722,1.351368916988763],"highPrecisionTau":["0.5304247761233722929035515273117735734691063121038204917591072161672333160619649915843","1.351368916988762970770294745889279163084223251591872489570127655742303613668144595193"]},"6,2,0":{"tau":[0.23032393881576027,1.338337097595142],"highPrecisionTau":["0.2303239388157602807394385777262604174230204012561773770413000536344732021757366853968","1.338337097595142137458159423329153353790495744771964112051984353575422073599371549918"]},"7,0,1":{"tau":[0.7578821950635308,1.568741781822955],"highPrecisionTau":["0.757882195063530816243557963662140429799946864116175307844386971540554317879462319071","1.568741781822954938435694040696023428476506280080328571933415288907741482138038118804"]},"7,1,0":{"tau":[0.3249304151620226,1.599097138180921],"highPrecisionTau":["0.3249304151620226062317938303992132620307821655213726322458405249021832229235787045152","1.599097138180920899455858452199476350642370889496044188709018979357701815330376975143"]},"8,0,0":{"tau":[0.5783347884054936,1.9907964524404396],"highPrecisionTau":["0.5783347884054935449942169978418890387880981656433728541099263446594450869351617530501","1.990796452440439632300880746142202742461386479495388196960284304927881257137275945998"]},"0,0,9":{"tau":[0.9052863407274091,0.42066837982167543],"highPrecisionTau":["0.9052863407274090933574655407508129107675384959949349928159072893343066456665583657345","0.4206683798216754501742518183810390914028967073829707442709249865630389854379704428759"]},"0,1,8":{"tau":[0.7975003961237308,0.4843360211540446],"highPrecisionTau":["0.7975003961237307925185014346728138893363248495699354904329837911932262889758660018863","0.4843360211540446346988434748030080160755093550880718164236123623353545999493363737415"]},"0,2,7":{"tau":[0.7102810744526495,0.5216545187363899],"highPrecisionTau":["0.7102810744526495547836604462144826759260183481047176257615361329946656200108277121228","0.5216545187363899425828053261378549346801346650131571586022134506060360766383946306552"]},"0,3,6":{"tau":[0.6302881075104035,0.5442438663912446],"highPrecisionTau":["0.6302881075104035903229456678402202537373441938150574962530570038080310451433110931012","0.54424386639124453316982290351997029623721989185919928106236717211215366332431450678"]},"0,4,5":{"tau":[0.5533301666697178,0.5558690621223308],"highPrecisionTau":["0.5533301666697178126144911742925652820883795753852642588814807433782622324798520716104","0.5558690621223307634721158603176359462800825367786263051443129160328631551311961403577"]},"0,5,4":{"tau":[0.47722162225118536,0.5579177435088031],"highPrecisionTau":["0.4772216222511853600858835269279170748359638770973882762194610560980397222435732472427","0.5579177435088031317513698140040679174416197672779382340718145864876419620510972934248"]},"0,6,3":{"tau":[0.4002268745047575,0.5505193291289712],"highPrecisionTau":["0.400226874504757525678208186415780277659221503894824154505797929223420448510859259029","0.5505193291289711242101637705322417783379089512842106339708839947776748937096196802685"]},"0,7,2":{"tau":[0.32022070355933313,0.5326224962894985],"highPrecisionTau":["0.3202207035593331111754600935860983837288505790236501903824677154138703400711322442626","0.5326224962894985510170968657169428341472627844882339418116028707818250276314916765483"]},"0,8,1":{"tau":[0.2332316901885675,0.5011748764210534],"highPrecisionTau":["0.2332316901885674762306188725043922003629316132064274016113584946891834876912515886629","0.5011748764210534281316651678349490925535081635727841955527104148974785157704082629047"]},"0,9,0":{"tau":[0.12693194518398745,0.44723328196147555],"highPrecisionTau":["0.1269319451839874673464908159514107961815052180809058674960209823010786530691501634155","0.4472332819614755266932946235641348228192144149762099758844047448197555924034847402909"]},"1,0,8":{"tau":[0.9136737045620708,0.5545792734711206],"highPrecisionTau":["0.9136737045620708090302697229127568758004394669925635161537679997550649005312924402306","0.5545792734711205896967319504248535500487531465710751478477666849992396305070789808501"]},"1,1,7":{"tau":[0.7962008484973406,0.6030158039280128],"highPrecisionTau":["0.7962008484973406178721941447930601291602255808323300774585965539749463184125226623231","0.603015803928012846324211040985473753242084682657123462525188774150043598689918961324"]},"1,2,6":{"tau":[0.6967284862987426,0.6329664899109504],"highPrecisionTau":["0.6967284862987426641866793036719047699719434886639437814026095067528690999987261328901","0.6329664899109503432732136999176182388449811432040839603155515564190443427464516323853"]},"1,3,5":{"tau":[0.6041306363884229,0.6503436645774036],"highPrecisionTau":["0.6041306363884228557502006920331734344724603529353731706179607306800144795339118769892","0.6503436645774036081017912365743501542354299246502485333885710582198931470298102580479"]},"1,4,4":{"tau":[0.5141192369559593,0.6572773685404725],"highPrecisionTau":["0.5141192369559592792709500114077572264859329963900439733200169500581353076623811810229","0.6572773685404725191778986828762327668829139286897701984468574000011194892605643916367"]},"1,5,3":{"tau":[0.42400586746770635,0.6544141023426473],"highPrecisionTau":["0.4240058674677063535778671070524931226949066990464206820370368509730643787020103804959","0.6544141023426472490498000251083046156077314153838757226308415551784632887905331333022"]},"1,6,2":{"tau":[0.331101689698656,0.6413589740543262],"highPrecisionTau":["0.3311016896986560142877321597604885927890313077809931787475044536816301506431682881646","0.6413589740543261928013104360404355912574671049591463290180222742320911755078117563938"]},"1,7,1":{"tau":[0.23111429379408077,0.6164446751774848],"highPrecisionTau":["0.231114293794080774734709916435274105640496373845170191112186784511365692178417056024","0.6164446751774848557927995217538536054054945620681890807514058669178615449037385721614"]},"1,8,0":{"tau":[0.1127983490862128,0.5754377482669324],"highPrecisionTau":["0.1127983490862128026723207129704993286534318581960549906864635796149972458618263272367","0.5754377482669324021247511003045992902305942659582499942333007252203119200267442514942"]},"2,0,7":{"tau":[0.912397230479871,0.6697708427064374],"highPrecisionTau":["0.9123972304798710233551449392158518420623014388482547438593901095763382398682881544749","0.6697708427064373722550906970451528586060565160155593271002276014466950828663989680567"]},"2,1,6":{"tau":[0.7844552780175824,0.712614200135365],"highPrecisionTau":["0.7844552780175824559523677328156541431641699376798310741052114611927945425980838945711","0.7126142001353650784933042507499629538189504197944506071473225801889392567889974182011"]},"2,2,5":{"tau":[0.6725187345711167,0.7384467602617243],"highPrecisionTau":["0.672518734571116699509648248567620869590052697600849927689976492589033743122088277113","0.7384467602617242638895158750170916804449339309447690274974213314192192288412516703405"]},"2,3,4":{"tau":[0.5665518369109653,0.7517816841966645],"highPrecisionTau":["0.566551836910965376811215515499981989396160740102495819024262764287302550402437209547","0.7517816841966644987920304654669251654429645428328287775084946792219433562280960512687"]},"2,4,3":{"tau":[0.462099501317581,0.7540950926399174],"highPrecisionTau":["0.4620995013175809751365881589646015186172964949128080456782672694646345233863322063767","0.7540950926399174051575433734801414490575100550333929504127231535815505275945543894688"]},"2,5,2":{"tau":[0.35572927723620884,0.7455407913806367],"highPrecisionTau":["0.3557292772362088599462014360184847363117094202723312641418986815804674191840850848558","0.7455407913806366890729194851334020783531036014306584619029509221923032574648917307516"]},"2,6,1":{"tau":[0.2429065728888235,0.7251063477057245],"highPrecisionTau":["0.2429065728888234914520767490887756869121668542065150079358598021691238029292635037643","0.7251063477057244233653477535706407303515192708196040651153239964134703973409902566328"]},"2,7,0":{"tau":[0.11317042601141826,0.6898883177182566],"highPrecisionTau":["0.1131704260114182603773193399040641504065796279324918980588070242665956990051494346158","0.6898883177182566998889399435842244800485389539686358983124274398680990397603885024203"]},"3,0,6":{"tau":[0.9067124091344732,0.7826787101347856],"highPrecisionTau":["0.9067124091344731688970616042450424253422087573943836975189503163502097036248656387906","0.7826787101347856328659255790076509721280010198900544646213013365081847336184223401277"]},"3,1,5":{"tau":[0.7648628516050452,0.8233650732281796],"highPrecisionTau":["0.7648628516050451901649737417484720131456302364043250320530011234460666149009172628365","0.8233650732281796908786256721945353599698313939964480950998495671957493211159273081736"]},"3,2,4":{"tau":[0.6376872271367137,0.8462169074874232],"highPrecisionTau":["0.637687227136713692389101076206900374077254958710176828284816497222860470409771500164","0.8462169074874231820291405633032274111181907960952365772135115521883230929141456037989"]},"3,3,3":{"tau":[0.5152900929628823,0.8554069855901278],"highPrecisionTau":["0.5152900929628822836528581388475037732007921843235275932162599754167869923539633709309","0.8554069855901278287832422052870640178468365049630618050236835470987623251728795540705"]},"3,4,2":{"tau":[0.39255422302925636,0.8519925774317841],"highPrecisionTau":["0.3925542230292563710671258291958081693541369434624626137500580087337046524167407801448","0.8519925774317840853410334584823676540853431354663232449464487151427740029581124192859"]},"3,5,1":{"tau":[0.2642680884070581,0.8354776252924652],"highPrecisionTau":["0.2642680884070581451919624813383072907130963964892052924114005973018074609245077641383","0.8354776252924651517505784017489025830115114615024064362949949022959021922828390982745"]},"3,6,0":{"tau":[0.11994084773381948,0.8033902657430243],"highPrecisionTau":["0.1199408477338194822076144059766371824840521494910512711288580447690026245484636724763","0.8033902657430243138076330885258966277013214562154481744079041095114580665843360007086"]},"4,0,5":{"tau":[0.8971639809366251,0.9011761716805341],"highPrecisionTau":["0.8971639809366250672743042601927497449066134340298483546784207548116755588474590333828","0.9011761716805341006341387697503417771662709037822146382054344951347691544755590878896"]},"4,1,4":{"tau":[0.7363335146569225,0.9414993163681471],"highPrecisionTau":["0.7363335146569225896056095986897527868704162765378335035840805068846183501390171958811","0.9414993163681471045605509404793011004426765444060786528979429861753359784130129819063"]},"4,2,3":{"tau":[0.589046567364497,0.9614434107222732],"highPrecisionTau":["0.5890465673644970253339346805175838287103213932787464968388493550602444130469152502118","0.9614434107222732830256634631536037337458410261825001830071642686004350475668468823742"]},"4,3,2":{"tau":[0.4446503039921923,0.9652054073582298],"highPrecisionTau":["0.4446503039921922773753010806485444721337889831379920696748438854533250544046536426695","0.9652054073582298162383827405868421132685020117824362794574607050477241085062372695125"]},"4,4,1":{"tau":[0.2961510986030259,0.9532245778049533],"highPrecisionTau":["0.2961510986030258968323712103383015971453260966887442406457633320456087390199002912863","0.9532245778049533073255361540309515960503288153590936935387863246660665481194948722813"]},"4,5,0":{"tau":[0.1321565057912849,0.9232359840162704],"highPrecisionTau":["0.1321565057912849159063484985339535710937959438495906816609912167782658319337586904144","0.9232359840162703935026412858810281078165234233917570458407097405231346892818578265263"]},"5,0,4":{"tau":[0.8827780402701756,1.0322177292958412],"highPrecisionTau":["0.8827780402701756350777133174304177959051939955519362855959347478554030785119848518481","1.032217729295841307202897115612417457515307690342044670481650791943150937219642832706"]},"5,1,3":{"tau":[0.6951374123938552,1.073296201030138],"highPrecisionTau":["0.6951374123938551278050364276625081983985782982238312778866440698487107086006684426653","1.07329620103013794432353602281106860970730282782262277430450120688996975906743252629"]},"5,2,2":{"tau":[0.5194643875961141,1.0892017355913224],"highPrecisionTau":["0.519464387596114085212663421340600326707612922598126042320064738646790716376657901353","1.089201735591322339158645736454516754174905175669510787091707552668557524950068709975"]},"5,3,1":{"tau":[0.3427400369916056,1.0840791608266542],"highPrecisionTau":["0.3427400369916055953361965128442098826897607337436190218394253286071881262306928307647","1.084079160826654239184596381051429318681089404693444784370670120726591839633086721442"]},"5,4,0":{"tau":[0.1511330994806695,1.0562653963191864],"highPrecisionTau":["0.1511330994806695059453339788005540858805268067218644916558261871392963397063835393901","1.056265396319186446466899157028464104543469282505696853785902643266826837405147040975"]},"6,0,3":{"tau":[0.8609137234808509,1.1848760155011235],"highPrecisionTau":["0.8609137234808508804462395713111362996433954975936203226683278384115035837740848475941","1.184876015501123553264030995614024506893891229057355377403327855230051540120354729435"]},"6,1,2":{"tau":[0.6330312822893928,1.2272136268889113],"highPrecisionTau":["0.6330312822893927931063097990881707581122506709052680204369328790861868759015100731396","1.227213626888911297497891711506906898884561583949231655321987259356226450167081898486"]},"6,2,1":{"tau":[0.4134077837747312,1.2353204330404501],"highPrecisionTau":["0.4134077837747312055298928357946782724953315515742731611173272239220235684394359034402","1.235320433040450184674122879853663872085784308935157944300098695796774451355206071167"]},"6,3,0":{"tau":[0.18062458321973887,1.2115308098835826],"highPrecisionTau":["0.1806245832197388628016613301944947864290494186961581561406386285046303316419963568514","1.211530809883582598521283222319243504451398101136800193851900794565956289818332369982"]},"7,0,2":{"tau":[0.8252846473115993,1.3744349481415081],"highPrecisionTau":["0.8252846473115993000583257520761830566102880039392885966401551946888580915744941425337","1.374434948141508142007558931006539798467975867422874670639601394769408060263885480287"]},"7,1,1":{"tau":[0.5302474343356233,1.4167199794724403],"highPrecisionTau":["0.530247434335623291370265435945443505852681546543328215391861699974266702319666524666","1.416719979472440375979215389697480104012510454645450885047454196998656734733217679439"]},"7,2,0":{"tau":[0.22972834561982114,1.404035548910081],"highPrecisionTau":["0.2297283456198211473373909361349059383789173447292782244090900027264500810104787817743","1.404035548910081081333490725876124823816072384770694760343358163296941925058929077583"]},"8,0,1":{"tau":[0.7579388342521736,1.6332140920183909],"highPrecisionTau":["0.7579388342521736071410984610155705822045679910995072560583750497771945291323161058004","1.633214092018390830262229809864757523547625202964731672135995441047809910767320523534"]},"8,1,0":{"tau":[0.3246509413617829,1.6635454197877955],"highPrecisionTau":["0.3246509413617829319741417580776544701746355968294362578770034702582368764437007533204","1.663545419787795430809478641798195533619487705548094193136762971460344688364459721837"]},"9,0,0":{"tau":[0.5782813053022888,2.054602837847386],"highPrecisionTau":["0.5782813053022888666224578892797046972378670816647624564493195497899442557309325870068","2.054602837847386278392741945434590330852437132605779492967572371930483245461452510423"]},"0,0,10":{"tau":[0.9097487527704582,0.4111306827603229],"highPrecisionTau":["0.9097487527704582311829954035110662226462386637572822797406447692212361802232315758738","0.4111306827603229112076606665239937209981640300906934810190325052363006172169388998746"]},"0,1,9":{"tau":[0.8085205030604694,0.4731841000501837],"highPrecisionTau":["0.8085205030604694461325300212907446310711489163497119382313925288517186144290727204435","0.4731841000501836812793093361737933746741014432665465208931020638496840365821414556122"]},"0,2,8":{"tau":[0.727484913591236,0.5104695998726397],"highPrecisionTau":["0.7274849135912359935395427827278660431160290699415985677085502471519823022808833463531","0.5104695998726397307961867603762266663607305513313859679442355736331953669971344523805"]},"0,3,7":{"tau":[0.6537225575156654,0.5342070175203936],"highPrecisionTau":["0.6537225575156654015733107991827132789060879075169381315911251374713589078510508778357","0.5342070175203936190463111646672820642209455189936083210614496995384313437659833730767"]},"0,4,6":{"tau":[0.5832557827619747,0.5481280598591305],"highPrecisionTau":["0.5832557827619746525809213563212223583141684552349144256131657007918514102492268615786","0.548128059859130543893974458794355599356452913421144200562674016359558380508406435414"]},"0,5,5":{"tau":[0.5141321200052162,0.5537657245792036],"highPrecisionTau":["0.5141321200052162539130890920804087891171774889777019904260412534715358255103508080571","0.5537657245792036411719277629055726605103724618333836811983801444200981167317507660221"]},"0,6,4":{"tau":[0.4449961846412119,0.551605558803055],"highPrecisionTau":["0.4449961846412118680214941970304425110577902587156706714825489809457540713632857984239","0.5516055588030549979192325176288147562056255097815412790737486782235211234190840691204"]},"0,7,3":{"tau":[0.3745100535182204,0.5413557437284982],"highPrecisionTau":["0.3745100535182204283520743143631119892186922743516572439401310895923405283394766872266","0.5413557437284981652834254377242361187235346496078510938741787571912318084191642101086"]},"0,8,2":{"tau":[0.3007934566085034,0.5217792175341178],"highPrecisionTau":["0.3007934566085033701396994595313023709638607462420587085137884836541699047048467153618","0.5217792175341177899617030766835046022844365339677729222142071619310117329969797166622"]},"0,9,1":{"tau":[0.22010452946094078,0.4898137606886881],"highPrecisionTau":["0.2201045294609407835751005158402713186693271698348910471300566678598983735334896432693","0.4898137606886881105952215977268522842130944006575020736278451175042405721308286921784"]},"0,10,0":{"tau":[0.12058200288812793,0.4367980258823749],"highPrecisionTau":["0.120582002888127924364146423289292650954714253320446667603822165196076530125928182836","0.4367980258823748597960042388667750158583041418653129033955093339879957825050510735667"]},"1,0,9":{"tau":[0.9190384211101494,0.5375538951023546],"highPrecisionTau":["0.9190384211101494734414068157153184941830995069960151844715559069973313292736108285195","0.5375538951023546051355336994158317210659610343969470049771658910507894633127937016922"]},"1,1,8":{"tau":[0.8100412482120526,0.5844412215281697],"highPrecisionTau":["0.8100412482120525647923264364386024363933817660384897794440716628263333836182740328968","0.5844412215281697052887081226217038735876655266253464959166425205695879671939935788096"]},"1,2,7":{"tau":[0.7186835859387825,0.6145586210627302],"highPrecisionTau":["0.7186835859387825008478073229841623513711627998302965995382667425814616059843835336601","0.6145586210627301849283773345826778611031793580435261196197978751268513870922364041437"]},"1,3,6":{"tau":[0.6343628118808657,0.6334855626949449],"highPrecisionTau":["0.6343628118808657171002008285489444703199905391913219967473203546160192391994569092755","0.6334855626949449364723893762855386002688843998754264474544404896866235869189131962729"]},"1,4,5":{"tau":[0.5531333364818894,0.6433577571540114],"highPrecisionTau":["0.5531333364818894526888977167787875419595852753030661198563342677178280107223810770044","0.6433577571540114968711169143304900486208760167469335533524581659919433770517220488689"]},"1,5,4":{"tau":[0.4727572835445658,0.645028601098364],"highPrecisionTau":["0.4727572835445657897711160583205570614008491925779982896904681554112685775000780754485","0.6450286010983640168787039079495773279267674572764228743627301796170928524208020880567"]},"1,6,3":{"tau":[0.3913913054166512,0.6385772047483168],"highPrecisionTau":["0.3913913054166511952285889006601062657019995261046405197283409805910274780007300042907","0.638577204748316862200419338637901220761559077739381215657451965198834989381961572511"]},"1,7,2":{"tau":[0.3068027654273528,0.6233524874653934],"highPrecisionTau":["0.3068027654273527801354760197334907454378698255229963593347090798372196687991845151513","0.623352487465393433118506052669513049986195884241602511886819611303105885739359089508"]},"1,8,1":{"tau":[0.21505973900580513,0.5976350819356401],"highPrecisionTau":["0.2150597390058051379754469266165959426576872160433563899853342716296071882471572163058","0.5976350819356400397629158313504199220273680331012108034666905625106353839353682466911"]},"1,9,0":{"tau":[0.1055076289568248,0.557406650482888],"highPrecisionTau":["0.1055076289568248023481595975531730096860213820736295322227858294291885041810006943752","0.5574066504828880025758448129572641317054210702959268357067886038946615196839498946032"]},"2,0,8":{"tau":[0.9190606076711629,0.644560676162795],"highPrecisionTau":["0.9190606076711629274439796838945805462146440830850808157776831872000570389849940982253","0.6445606761627950106064310819017378035457942738271093642026916446636153861985829935983"]},"2,1,7":{"tau":[0.8020005852137946,0.6857419783728032],"highPrecisionTau":["0.8020005852137945376806602375315539194663056368284454984861325546528423225117853208063","0.6857419783728031833905837069669880707954548569491661238578471078660157550249000674937"]},"2,2,6":{"tau":[0.7006387670582919,0.7119871111940945],"highPrecisionTau":["0.7006387670582918332387796608590107182883483268733914559694400511330150745910301771306","0.7119871111940945587144826743549621709184614929093861202204839700971223586507505702888"]},"2,3,5":{"tau":[0.6056123636276993,0.7274684775197592],"highPrecisionTau":["0.6056123636276993634973917783410737630624409921247238639963484480079692939253212005269","0.7274684775197591397685519673784359392739974716684421876771786701687398925876687246496"]},"2,4,4":{"tau":[0.5130170313391375,0.7337170748264283],"highPrecisionTau":["0.5130170313391374209379406845900595514276440388821862875501239432627502587866916772142","0.7337170748264283347082754095583634994920907111821322568407822311753797436915978352951"]},"2,5,3":{"tau":[0.4202900175632987,0.7312134305812743],"highPrecisionTau":["0.4202900175632987047160961264515075506090430336500358024595057108820127965409537992383","0.7312134305812743003521169461276234399628337575368879709044542445185909295130946274858"]},"2,6,2":{"tau":[0.3248540459501834,0.7196891904709856],"highPrecisionTau":["0.3248540459501834309459843068053545938134694218663553733026307358606862067187822183142","0.7196891904709855471911517553488473634095238095420442251182240258796868392991067839124"]},"2,7,1":{"tau":[0.22273975015311986,0.6980096115053703],"highPrecisionTau":["0.2227397501531198651433313320988675161196668945109576118145108883241770219817143294766","0.6980096115053702959954238328694307944978066457308445010510781240203225356737588234232"]},"2,8,0":{"tau":[0.10424399065456832,0.6634698710594106],"highPrecisionTau":["0.1042439906545683291718902873849519542274907995590766198004319737674953111976726567331","0.6634698710594105606397609789980444569207715013460381069084018044272690332220807364215"]},"3,0,7":{"tau":[0.9152155571535485,0.7476036605229398],"highPrecisionTau":["0.915215557153548468126933152087089212048919177872991104561428040462611983880806250385","0.747603660522939788804215305558641288059563507552877357195140566637251186545585148456"]},"3,1,6":{"tau":[0.7875582513635886,0.7864254368915393],"highPrecisionTau":["0.7875582513635885825426314228104352165205668181688691612984989469762724327962179888028","0.7864254368915393128042702666085320044317383731117380156665638558676087627062384123145"]},"3,2,5":{"tau":[0.674390219373184,0.8101019886409881],"highPrecisionTau":["0.6743902193731839819010930133783194490048500157888232985715857536813346639692518209884","0.8101019886409880915423626376911826816118959019214989170082962195948229277093482708294"]},"3,3,4":{"tau":[0.5667527982850816,0.8224381075576594],"highPrecisionTau":["0.5667527982850816648839931090648323464432090014620580335562483638099708287623314671345","0.8224381075576593991495193455965995652255251101665397296639905228868299429728167606083"]},"3,4,3":{"tau":[0.4605137349559376,0.8246607204324721],"highPrecisionTau":["0.4605137349559376121837480416003691666454035518483808104097466460820842198677965307948","0.8246607204324720256693087191946098115664984552979044745227484783267919009620108539364"]},"3,5,2":{"tau":[0.3524306517251559,0.8169135208837309],"highPrecisionTau":["0.3524306517251559136326060199716910287116877149620732967528125040174791219301267793335","0.8169135208837308901605600787708190100444585703242345740811855986046857886212957903361"]},"3,6,1":{"tau":[0.23827210887070233,0.7984044507900027],"highPrecisionTau":["0.2382721088707023285657601886176418563283061581885434224176641925618732914189152037137","0.7984044507900026490289665757238428598578639553713252111365817273141494737996838400709"]},"3,7,0":{"tau":[0.10858981240248783,0.7668390401184514],"highPrecisionTau":["0.1085898124024878289095360127603484420903187604005999931922523136044288814151233560111","0.7668390401184514229588083289927586075974034227757814445188459228232170669839402104658"]},"4,0,6":{"tau":[0.9083686362080391,0.853471128910256],"highPrecisionTau":["0.9083686362080391457666791940200675722152313361562062249082895557972276924118284631071","0.853471128910255949936534148664034624515691297799938283900417061917484306378069706874"]},"4,1,5":{"tau":[0.7665877408192001,0.8916706712153726],"highPrecisionTau":["0.7665877408192001329905185106113827771956664157203096030944492303042890073970703898582","0.8916706712153725732472175471566519217860650701134100080232073449920759010124546147843"]},"4,2,4":{"tau":[0.6384606435127549,0.9132281048057832],"highPrecisionTau":["0.6384606435127549352801982705834649560645076611794193652814094787908874199039396926855","0.9132281048057831927041283092170043265766826670672219226464690133765713338128212167624"]},"4,3,3":{"tau":[0.5148113446942287,0.9219781522126247],"highPrecisionTau":["0.5148113446942287652480235115864151282683380867359402676207511631921227755297359205891","0.9219781522126247401418852981007161464767758836182221092137687001648563067180422423846"]},"4,4,2":{"tau":[0.3908091918996164,0.918884894011992],"highPrecisionTau":["0.3908091918996163836458659262274802317775619516588918214775814044704967815186201478188","0.9188848940119919800287743156866409467396296417159187139679171533834508250295825794801"]},"4,5,1":{"tau":[0.26152358969543366,0.9035272110690888],"highPrecisionTau":["0.2615235896954336429776629463742982070762532839152471091340575580685031369792294227192","0.9035272110690887969832850014441621721656072796709315364741178719496096650128120569616"]},"4,6,0":{"tau":[0.11716117880608924,0.8737082193732995],"highPrecisionTau":["0.1171611788060892451857549314140909689415255113781764288659893007546743732263561833823","0.873708219373299492559018781740758144569065771819293798520362242239111949563523937017"]},"5,0,5":{"tau":[0.8981825574879845,0.9674058441093298],"highPrecisionTau":["0.8981825574879844968008303508699159729493082402153648136797114476296636247072736240272","0.9674058441093298608671873605620390148584536349497528289507791184403153614642434777689"]},"5,1,4":{"tau":[0.7372529540960431,1.006118644389196],"highPrecisionTau":["0.7372529540960431118379173417386997017693667787705213431904437142440565016271391053633","1.006118644389196154714145333114457346482127526953547296292074479679154083568279297935"]},"5,2,3":{"tau":[0.5891746214659662,1.0253215576534365],"highPrecisionTau":["0.5891746214659662540269123127360659631125538044747944560314951278854532165731587993903","1.02532155765343641257042384978835300428377858234390881499639696812439624915913642686"]},"5,3,2":{"tau":[0.4438246155489609,1.0290425985006992],"highPrecisionTau":["0.4438246155489609235405631714983246615783681388064192360472679579110397061724478507396","1.029042598500699192827885142119258984667613643194400667502949400694164298975229725723"]},"5,4,1":{"tau":[0.29451315034278686,1.0177125624401901],"highPrecisionTau":["0.2945131503427868776148001721857462099366908671280518446110750507678732539068311898454","1.017712562440190132856550601915848999195062467463777237329074818318824177224000077518"]},"5,5,0":{"tau":[0.13037027880383903,0.989190627497009],"highPrecisionTau":["0.1303702788038390392375269262955746047929111977664795859567259612370579644788552600754","0.9891906274970090102436547282805374741391632654015093294324427640920943712698190506554"]},"6,0,4":{"tau":[0.883402440379279,1.0953251469013674],"highPrecisionTau":["0.8834024403792789959274272999066395263440085710054821217378132088959223602246557017956","1.095325146901367339503445097605014389740876032639153409556472191418512794380560134073"]},"6,1,3":{"tau":[0.6955470786137932,1.1353624094097583],"highPrecisionTau":["0.6955470786137932186484695754440777118446226735144491279950848159844484478751774071276","1.13536240940975824371812634597277727453954446849563882286301424200000001944233963053"]},"6,2,2":{"tau":[0.5192082355101879,1.1509252351605825],"highPrecisionTau":["0.5192082355101878586987562359564658629213431705675248341361258774251117101408530897706","1.150925235160582623979969560400804818889613891330726103926416332383358085565399627666"]},"6,3,1":{"tau":[0.34181112754341164,1.1460818676471523],"highPrecisionTau":["0.3418111275434116429031379723650755776158122146642467238990110130692397616301956528584","1.146081867647152239808267507603993618993477696656240909349007403047927359806213269666"]},"6,4,0":{"tau":[0.14996694760821253,1.1192104404698613],"highPrecisionTau":["0.1499669476082125125658653848884694861319420582282715304973718039290242583759344480003","1.119210440469861389612073201726076277602979073268243535774180966814943695342483581426"]},"7,0,3":{"tau":[0.8612740853110933,1.2458070785396733],"highPrecisionTau":["0.8612740853110933211804932394673768494157801333943319613603535929576795724269746790612","1.245807078539673288035104164996755037698840444366572842188697050903851575749902876496"]},"7,1,2":{"tau":[0.6331074021030786,1.2875138975974518],"highPrecisionTau":["0.6331074021030785710283419891071225812777509542159872519948644404427526703392432852283","1.287513897597451828712165164451478867456389610023731712742623509621866226395591879642"]},"7,2,1":{"tau":[0.4129566717215887,1.2955967017886751],"highPrecisionTau":["0.4129566717215886590107519143645261951590381832713242953429149710234336264142508331324","1.295596701788675059412625564017237143508295043317531832958960714583366077844275595403"]},"7,3,0":{"tau":[0.17988051687742626,1.2723691650393238],"highPrecisionTau":["0.1798805168774262682696285910572521078397276674812932367075854438862054030446352100603","1.272369165039323829413312405120135008250489204223627371718507850988920287205181474767"]},"8,0,2":{"tau":[0.8254589686702685,1.4338726835189464],"highPrecisionTau":["0.8254589686702684530890600983682051093715907895346290630049311433994735816854119749681","1.433872683518946311896827847023965149288501880569966830612317177863291293901961316042"]},"8,1,1":{"tau":[0.5301163985559172,1.4758532559781143],"highPrecisionTau":["0.5301163985559171730436915979856548030398756299364931963909465181080504029114021596286","1.47585325597811434419078961162265158574456721548040651872987703420597664938281028863"]},"8,2,0":{"tau":[0.22929014552917062,1.4634259603517479],"highPrecisionTau":["0.2292901455291706287898755210663644004394053033984805609364849069833568307067577726391","1.46342596035174785503170908040767280450708732751666609644746365378068628360584781516"]},"9,0,1":{"tau":[0.7579809867273274,1.6916950098860004],"highPrecisionTau":["0.7579809867273274340394245737183159245765860180673385166831415555930107829854021216083","1.691695009886000421107935827069443226085815593055247133634483862812335516931393463699"]},"9,1,0":{"tau":[0.32444264395471506,1.7220087308402454],"highPrecisionTau":["0.3244426439547150418335005243447246150096517966583535753764461516245644223022983886183","1.722008730840245388050020524575257807829901432857733440353060252871934701972630824612"]},"10,0,0":{"tau":[0.578241177463598,2.112585796401073],"highPrecisionTau":["0.5782411774635980251332596435413468642018407820538830672579350519913895905666842583884","2.112585796401072819440156958687628931747672980445556593406961587102019295346519233618"]},"0,0,11":{"tau":[0.9135705198512508,0.40274319084995835],"highPrecisionTau":["0.9135705198512508562467599637767572126303768773539714967606643564737534529199501803088","0.4027431908499583799739497243638518639527090875514927017746256062012694588057514798861"]},"0,1,10":{"tau":[0.8178748914036359,0.4633057764560679],"highPrecisionTau":["0.8178748914036359330092484774787552186013220338633505644636773410181721796082837790445","0.4633057764560679163945783701572315116570122859502497957648261008274761921686029932926"]},"0,2,9":{"tau":[0.7419995615345809,0.5003921913931301],"highPrecisionTau":["0.7419995615345809607337706927779890021035540272722804221086610593180162286931281295759","0.5003921913931301079879780165454001523106736869464091347533425313685401314329906588707"]},"0,3,8":{"tau":[0.6733854903201331,0.5248610171571945],"highPrecisionTau":["0.6733854903201330830647098645012613107325080345853528708303151699694364994244030372238","0.524861017157194505522759918294313434010238866869863828407647708961446291235610751469"]},"0,4,7":{"tau":[0.6082149762714262,0.5403819371013867],"highPrecisionTau":["0.6082149762714261919572433561154190972248505472803556382347788467162908053644981295088","0.5403819371013866664308420901061839104598229849475971034264565423782795442616787303288"]},"0,5,6":{"tau":[0.5446816420698254,0.5485497483853404],"highPrecisionTau":["0.5446816420698254002648353292563199420447806096737747389266421237310057593862858126666","0.548549748385340430887979951425264355742901157502387647892862556943679236353527525421"]},"0,6,5":{"tau":[0.4816333320402924,0.5500281308875302],"highPrecisionTau":["0.4816333320402923783139300154037572748650955990957191979509067933286525864316144463136","0.5500281308875302720116404501786520424173878835068226151855829525676565878821136979835"]},"0,7,4":{"tau":[0.4180892754139823,0.544884048450245],"highPrecisionTau":["0.4180892754139822880432662587876217148713962526809297802163036529310828775384026233337","0.544884048450244936957198173870475992132206981917353671039463452967242552613265747856"]},"0,8,3":{"tau":[0.3529239050607367,0.5326145495488299],"highPrecisionTau":["0.3529239050607366655454066099317006323253306972639389023586768412831034236155953611566","0.5326145495488299203852863647359907738304041120920926111598882625057846486571514217447"]},"0,9,2":{"tau":[0.2844107691180769,0.5118959676022632],"highPrecisionTau":["0.2844107691180768892572073966021934596492498012139616150078080480046475700445396314516","0.5118959676022632057725436009030211592419540009734332206161802664685760032039710926685"]},"0,10,1":{"tau":[0.2089795013447162,0.4796957634301651],"highPrecisionTau":["0.2089795013447162036851714954156573216162778675977840078032491571880934273857624564103","0.4796957634301650971968214709792949534522231574215224352623510794378931038034110545238"]},"0,11,0":{"tau":[0.11516327836048905,0.42761431079597856],"highPrecisionTau":["0.1151632783604890533432331821087081478349047040535333043184212405714026241269235430836","0.4276143107959785683501567588329211986750411941954409290434918200923864812322377450748"]},"1,0,10":{"tau":[0.923522083580442,0.5228436920813326],"highPrecisionTau":["0.9235220835804419729577901991715419987250549366151670161340742648752549472330305850573","0.5228436920813325850993154657244999545263440604383660804820052713429442048377723484474"]},"1,1,9":{"tau":[0.8215231642780433,0.5683143809096617],"highPrecisionTau":["0.8215231642780433163580964151508197060617665715690054629989497703111854809781069165924","0.5683143809096617399450528390543231460022710504087478773373460556851723344594939528437"]},"1,2,8":{"tau":[0.7367880463090198,0.59836942505751],"highPrecisionTau":["0.7367880463090198132721744920254927279897686533349195982415736543507545514329209996966","0.5983694250575099878466405617122802907415133315223020647848599322052468157456215326461"]},"1,3,7":{"tau":[0.6591386056231198,0.6182916746489376],"highPrecisionTau":["0.6591386056231197271435140457983655938322555943486422364954586633700828542286786548335","0.6182916746489376934191796386591662131170555386956679966406374860962732272471142317643"]},"1,4,6":{"tau":[0.5848640159391224,0.6301805299640165],"highPrecisionTau":["0.5848640159391224500202942814738757445345565111037668203820102124086932036139173125612","0.6301805299640165664700033413926517568730686010027946935039420256944468510877846810768"]},"1,5,5":{"tau":[0.5119820863111095,0.6349826256776173],"highPrecisionTau":["0.511982086311109493507537708981070304293983501675918253805120518826341679194922287459","0.634982625677617368933222150397826033005457464953172667965156308115136494190607712075"]},"1,6,4":{"tau":[0.43905263991071825,0.6330082494050312],"highPrecisionTau":["0.4390526399107182211685785980424288875456818676329292572606016986009538480026142471781","0.6330082494050312081414004218894735169002607728985982453499232710112852169840014532377"]},"1,7,3":{"tau":[0.36463837528057746,0.6240648804722684],"highPrecisionTau":["0.3646383752805774757056101600163852569046239707737847813070260240295458678045085279706","0.6240648804722683448542266968189657173543669460794700737676483979638659290423613522375"]},"1,8,2":{"tau":[0.28676838049023473,0.6073809424113679],"highPrecisionTau":["0.2867683804902347196048770132181734408377724815427078859471807379372273847519093888711","0.607380942411367894129323381781689047775785734647530377864089215402518663188723312206"]},"1,9,1":{"tau":[0.2017602409055539,0.5812406655903369],"highPrecisionTau":["0.2017602409055539024237780505374384754377283992962996719077680066461514799036460578812","0.5812406655903369127449874816711285234110524492044475708270495517964386733075509926213"]},"1,10,0":{"tau":[0.09943620012791543,0.5418265137711081],"highPrecisionTau":["0.09943620012791542854090753386088448607309338262289711771613578994760203294436265508793","0.5418265137711081193891078239328210091146579330393318042662484261089917850757301491147"]},"2,0,9":{"tau":[0.9244965945055537,0.6231238419640506],"highPrecisionTau":["0.9244965945055537033866823112469756256011078600394006562629326277668886820519369298468","0.6231238419640505857391914396910870295710390379396611990125659594416028634580588418946"]},"2,1,8":{"tau":[0.816210312738608,0.6628008359309334],"highPrecisionTau":["0.8162103127386079787121032444016782168771948896389734779048405187939297714940239719734","0.6628008359309334041428724950659488356060212683078733243359898863105974991188074412735"]},"2,2,7":{"tau":[0.7232639602197353,0.6891250629829302],"highPrecisionTau":["0.7232639602197352254594581415315386831025598502118060684430781336504283050946623488299","0.6891250629829301832698850544321089553722294994043497262523444491149964340155231611045"]},"2,3,6":{"tau":[0.6368087234473495,0.7059664572287357],"highPrecisionTau":["0.6368087234473494638359635669485754181674330806923235345901781168427281678922737541681","0.7059664572287357333740441603834797216648977568595418084074630706603682411036702325727"]},"2,4,5":{"tau":[0.5532820844409287,0.7148409867471551],"highPrecisionTau":["0.5532820844409287833357437879776394838444857435115188940921776552241306533686179426903","0.714840986747155085843323734383676438236257483028694181966811019171667380418960605625"]},"2,5,4":{"tau":[0.4705583079566005,0.7163735393510486],"highPrecisionTau":["0.470558307956600474507400414620380593782529567954615469077282222716876704060368922447","0.7163735393510486857681860643708441185343843942072203464329650772129084309593403998324"]},"2,6,3":{"tau":[0.38684796399943494,0.7106298176971176],"highPrecisionTau":["0.3868479639994349430225714719899321473338348149817255007728172931900630467559838116734","0.7106298176971176784786254972163817223954278795314535989028442801686490149947490554444"]},"2,7,2":{"tau":[0.300010590519598,0.6971539737455481],"highPrecisionTau":["0.300010590519597999647105645893867578735816833285816684458976981468653224201451326359","0.6971539737455481075217255949168136349537271229129299777667027749599942349319687208537"]},"2,8,1":{"tau":[0.2064339173398105,0.6747839892977446],"highPrecisionTau":["0.2064339173398105128226331189422485377291326658964226835679961980105183182403894913696","0.6747839892977446495462629976967210131987897901569154818185291436122593194174486047625"]},"2,9,0":{"tau":[0.09699457093004263,0.6410066821773549],"highPrecisionTau":["0.09699457093004263119384513204303825470175925112677859259127493842740811552845596814454","0.6410066821773549290063426871108965278667776990600113582278140323452289289917342030057"]},"3,0,8":{"tau":[0.9219625862618696,0.718300170447209],"highPrecisionTau":["0.9219625862618696061758688595737421905240876791857212140805567502537527809082236774776","0.7183001704472089347864251584625502352725791620628955207762600282659421751650401793775"]},"3,1,7":{"tau":[0.8054263691316562,0.7554462829596336],"highPrecisionTau":["0.8054263691316562442319122204624785224895826152683646303563315767336616436374168645858","0.7554462829596335624871084168398299413593267332432904986218220490608763299456931131483"]},"3,2,6":{"tau":[0.7030651717204842,0.7794188051572261],"highPrecisionTau":["0.7030651717204842244282132142943324863367355149652892781798256180198037492291777877131","0.7794188051572260179160896286115659931530288578172995012257030027097997007607704905724"]},"3,3,5":{"tau":[0.6065802114192296,0.793688038869172],"highPrecisionTau":["0.6065802114192295741928984953593369887253824383555008067435292972915409424512204255707","0.7936880388691720546704724054798711150377555932522159322470428498729507578499856056018"]},"3,4,4":{"tau":[0.5123726481422266,0.7995062578759802],"highPrecisionTau":["0.5123726481422265599740854126566306903558057254710494262822308803553720402195893492211","0.7995062578759801552128594227458599781736672533371420798185933127730759762653207654099"]},"3,5,3":{"tau":[0.4180178046003677,0.7972728654869734],"highPrecisionTau":["0.4180178046003676656672558986699584577526784127736206810543667586321849056077347280767","0.7972728654869733647550574697703155877700121476831940458820577605284989641581452164532"]},"3,6,2":{"tau":[0.3210707964935976,0.7867852240665781],"highPrecisionTau":["0.3210707964935976112350320523427308997437626026791105650375553431048167081070466800616","0.7867852240665781772616183315830587848102322068013026533684236450886326865187558994586"]},"3,7,1":{"tau":[0.21784528541401077,0.7671610731254896],"highPrecisionTau":["0.217845285414010761513443469105781761603094063108797437106675305710648210187395728737","0.7671610731254896729605482785285286417431246238828229964875914080279358042378041904809"]},"3,8,0":{"tau":[0.09963559115192062,0.7363017420446022],"highPrecisionTau":["0.09963559115192061956899548107391425144481966370200223286897439919284118394594493401995","0.7363017420446021480470998456497556090176589653007422197580369205508593401627484394277"]},"4,0,7":{"tau":[0.916960169040881,0.8144531131952817],"highPrecisionTau":["0.9169601690408809186246636828724219343185859652237106020939546142606121149980169840437","0.8144531131952817299749420385298760925556821143594853871666711333635953177386089857408"]},"4,1,6":{"tau":[0.7895758107902655,0.8507391923521915],"highPrecisionTau":["0.7895758107902655512486732280829055876718546644080402701156934234957882013813389147809","0.8507391923521914810130231776890859870621086206062339584841885362154923116651306697272"]},"4,2,5":{"tau":[0.6756412810453973,0.8729810884061157],"highPrecisionTau":["0.6756412810453973904055663821081764168174224819735103757989300866233262020124105056617","0.8729810884061156195313865115240512359438503425280966671173991774064624538167007606631"]},"4,3,4":{"tau":[0.5668984611129122,0.884639177811809],"highPrecisionTau":["0.5668984611129121994589108955069267205636449388407635843461345256383200046194607289302","0.8846391778118089982934913390838681225022673136458098671405745611633748223500631312644"]},"4,4,3":{"tau":[0.45946160702678335,0.8868103076086719],"highPrecisionTau":["0.4594616070267833741943909767192444741196871359463386593822496857776718624092806038259","0.8868103076086718204060716556428265437134748332249803252414822823836682594667293516968"]},"4,5,2":{"tau":[0.35024985048842067,0.8796329618867217],"highPrecisionTau":["0.3502498504884206435063172396714807211849384803266306496700798520331481697482089513944","0.8796329618867216644413114364845278139381497450162649475601183453383681095935270726631"]},"4,6,1":{"tau":[0.23526991624878854,0.8624294725194341],"highPrecisionTau":["0.2352699162487885284366729164465139443476285328176830084648017389138659937554696848992","0.8624294725194341708860914159579891683441220997592890794821697286906786086199163173775"]},"4,7,0":{"tau":[0.10575014763779851,0.8331886558415479],"highPrecisionTau":["0.105750147637798504363794784310411560877099474467030513040625825984552168658528493866","0.8331886558415479054521844678553935594189932959953897031016212910555971498483236217661"]},"5,0,6":{"tau":[0.909478336130215,0.9158150216727714],"highPrecisionTau":["0.9094783361302150557108391902872134878151302511125261183905395470129389570610486454509","0.9158150216727714061180676088322787470162951439676538495345755813981965715875651490709"]},"5,1,5":{"tau":[0.7677772064188196,0.952327367648329],"highPrecisionTau":["0.7677772064188196649284631197180791841321858366345961937191071888776850123902300524408","0.9523273676483290353445104806012995746565741833906869227211291770674063497632723065479"]},"5,2,4":{"tau":[0.6390028604706922,0.9729788345009572],"highPrecisionTau":["0.6390028604706921686131628287473855067271967539930432855073123454006013604195334238048","0.9729788345009571892015887753426938460860487074758568881528461516129652308485569437618"]},"5,3,3":{"tau":[0.5144782751688971,0.9814199027712702],"highPrecisionTau":["0.5144782751688970918713312926290899896154236098105216374807884654669169209476358732632","0.9814199027712701980263237095685257542441618212941891126421480657187929375164556881108"]},"5,4,2":{"tau":[0.38959268539149444,0.9785620709554678],"highPrecisionTau":["0.389592685391494444626287364037887056147695881806248057170551690733193720477209607699","0.9785620709554677761187422520986295121086858261294875040506384487197412080870911810467"]},"5,5,1":{"tau":[0.25963277795822404,0.9640250154428804],"highPrecisionTau":["0.2596327779582240347786605170179396712824517430939205156014848372751262484317789116044","0.9640250154428804291726894752281089905547195297946662794988161739487143233126190082015"]},"5,6,0":{"tau":[0.11529505721141707,0.935752735570142],"highPrecisionTau":["0.1152950572114170703893530136458360072725746142689572847685918920209894903142273661161","0.935752735570142029910630096606450211988774383154740880639740882367062441333513030863"]},"6,0,5":{"tau":[0.898900147925036,1.0266196244194499],"highPrecisionTau":["0.8989001479250360826020908296251193327980945899320215236940797700370075731198803956087","1.026619624419449899129963787214252374200126543433037019167616344311070475509497444695"]},"6,1,4":{"tau":[0.7379121304054201,1.0641875888152226],"highPrecisionTau":["0.7379121304054201072914374133142229745104293681082260911944907498704607030831524366738","1.064187588815222678016022646216693715098806258784928697832668016710333267443837216241"]},"6,2,3":{"tau":[0.5892673990558824,1.0828543465832914],"highPrecisionTau":["0.5892673990558824334081194886208872474692048165176690789166317775794190349145836274304","1.082854346583291306319814091792350454094987147847469199206105812456279483213549876207"]},"6,3,2":{"tau":[0.44322851898129456,1.0865484268975043],"highPrecisionTau":["0.4432285189812945556654389124176920227244013331260448749302817511867544285587546095924","1.086548426897504291271466771835886332481478380282823432225603503940668201792125567091"]},"6,4,1":{"tau":[0.29333838206637347,1.0756945523849286],"highPrecisionTau":["0.2933383820663734876143170777222620996630025901507038580475372079530507079181399650256","1.075694552384928558999578233253509634807386921547107769770771502313579670247796692677"]},"6,5,0":{"tau":[0.12910937368416608,1.0482188756958317],"highPrecisionTau":["0.1291093736841660978415936903634978882903302453792273876973833400043405925700048757768","1.048218875695831706997893658057598133796981413306640672636755721399383301736979324091"]},"7,0,4":{"tau":[0.8838573865685705,1.1523009263231543],"highPrecisionTau":["0.8838573865685704500944383792515242195836393120995967023832957374767500028505003361808","1.152300926323154222214222354062904767781832200894453505682903058107558516692613909228"]},"7,1,3":{"tau":[0.6958486702907646,1.1915745597379706],"highPrecisionTau":["0.6958486702907645589399707376838871480383138477220002282574834491371340761573506991421","1.191574559737970556016219188405528746475727333367392837091317087432202191982400679437"]},"7,2,2":{"tau":[0.5190187547853987,1.206884256897688],"highPrecisionTau":["0.5190187547853986731772602043277033549724554993818830449918504254568060454985503344797","1.206884256897687982650573772914494550978849532151797469442396240524396622418158710142"]},"7,3,1":{"tau":[0.3411264129888459,1.2022500080874907],"highPrecisionTau":["0.3411264129888458942294885644772922262796535119242351116675528934567107442822996262415","1.202250008087490795494728593931030177706114438407810760330002294102075733932591626784"]},"7,4,0":{"tau":[0.14911578112600785,1.1760716818208807],"highPrecisionTau":["0.1491157811260078541039496021823864906414465677516271271250007556208241244717948208489","1.176071681820880760877765400286112572008613113636319176406919723986683905479600152155"]},"8,0,3":{"tau":[0.8615427773268887,1.3011769323276752],"highPrecisionTau":["0.8615427773268886139879252450336681552126440034252649047498804971539866609268446597465","1.301176932327675126633026282892309352708089881706277570105579063232174367412963012843"]},"8,1,2":{"tau":[0.6331643609690645,1.3424113761332421],"highPrecisionTau":["0.6331643609690645133728193714605564506875411490191975167557432316176208189389026460193","1.34241137613324217470805706768958059916785599283596215216046255481264973233131841837"]},"8,2,1":{"tau":[0.4126178490351031,1.350476962940132],"highPrecisionTau":["0.4126178490351031033567431401002192881778207226364897695081527823189070343974142380196","1.350476962940132143661314446148842438090532344296801321282229274506424738406548009729"]},"8,3,0":{"tau":[0.17932493075464095,1.3276715574532627],"highPrecisionTau":["0.1793249307546409464908047704803337719381577931882926663647123159969628568842405782946","1.327671557453262642868664507926703143137033164427689571613915762144484120737301505696"]},"9,0,2":{"tau":[0.8255909930696255,1.4881190103136215],"highPrecisionTau":["0.8255909930696254385697971254395621258546023922022563765842967135839609326471939993589","1.488119010313621415543923385067656097928792988852525057929594667330087782494390661077"]},"9,1,1":{"tau":[0.5300167236508954,1.5298685448852112],"highPrecisionTau":["0.5300167236508953312960493092019995772294675595928494076986599394220193201498810906052","1.529868544885211157035627998543097612071083948731358961551148166444533027239211506622"]},"9,2,0":{"tau":[0.22895789502527408,1.5176370820909544],"highPrecisionTau":["0.2289578950252740884004843355837457432670495654230058020510583388769893642561626099664","1.517637082090954479274640733603402303340614603532317153326775269385697102043976096036"]},"10,0,1":{"tau":[0.7580132224444402,1.745212327495684],"highPrecisionTau":["0.7580132224444402370842221420850855951868034863207474220892669411656637373529934248402","1.745212327495683822505253994103933648187178247744828935507184861821561378726431900062"]},"10,1,0":{"tau":[0.32428317208656704,1.7755127431036775],"highPrecisionTau":["0.3242831720865670354123909088632218890267988839170007772557912146083859335600299058348","1.775512743103677402323272501950754971075152678273899048412310728558364569289836384503"]},"11,0,0":{"tau":[0.578210300862536,2.1657209292992188],"highPrecisionTau":["0.5782103008625360088420900191486876606838751668915027312211958223157391736243088166203","2.165720929299218616175686199131973625932322574539811368854081776190694379553616923172"]}}};
approachMath.setApproachAnchorData(attachmentAnchors);

const TRACE_TIME_BUDGET_MS = 45000;
const ViewDecimal = Decimal.clone({precision: 120, rounding: Decimal.ROUND_HALF_EVEN});

function decorateNode(node) {
  for (const branch of [node.base, node.result]) {
    if (branch?.status !== "ok" || !branch.tauExact) continue;
    const x = new ViewDecimal(branch.tauExact[0]);
    const y = new ViewDecimal(branch.tauExact[1]);
    const denominator0 = x.times(x).plus(y.times(y));
    const u = new ViewDecimal(1).minus(x);
    const denominator1 = u.times(u).plus(y.times(y));
    const error = new ViewDecimal(branch.numericalErrorEstimate);
    branch.views = {
      raw: [...branch.tauExact],
      zero: [new ViewDecimal(1).minus(x.div(denominator0)).toString(), y.div(denominator0).toString()],
      one: [u.div(denominator1).toString(), y.div(denominator1).toString()]
    };
    branch.viewErrors = {
      raw: error.toString(),
      zero: error.div(denominator0).toString(),
      one: error.div(denominator1).toString()
    };
  }
  return node;
}

function validateRequest(data) {
  const counts = data.counts;
  if (!Array.isArray(counts) || counts.length !== 3 ||
      !counts.every(x => Number.isSafeInteger(x) && x >= 0 && x <= 11) ||
      counts.reduce((sum, x) => sum + x, 0) > 11) {
    throw new RangeError("Enter three nonnegative whole numbers with a + b + c <= 11.");
  }
  const path = data.path;
  if (!path || !["edge", "vertex"].includes(path.kind) ||
      !Number.isInteger(path.index) || path.index < 0 || path.index > 2) {
    throw new RangeError("Choose an edge or vertex numbered 0, 1, or 2.");
  }
  const position = String(path.position ?? ".5");
  if (position.length > 60 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(position)) {
    throw new RangeError("Enter a decimal path position strictly between 0 and 1.");
  }
  const p = new Decimal(position);
  if (!p.isFinite() || !p.gt(0) || !p.lt(1)) {
    throw new RangeError("The path position must lie strictly between 0 and 1.");
  }
  const normalizedPath = approachMath.normalizePath({...path,position});
  const target = String(data.target ?? "1e-5");
  if (target.length > 40 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(target)) {
    throw new RangeError("Choose a target epsilon between 0.1 and 1e-10.");
  }
  const targetDecimal = new Decimal(target);
  if (!targetDecimal.isFinite() || targetDecimal.lt("1e-10") || targetDecimal.gt(".1")) {
    throw new RangeError("Choose a target epsilon between 0.1 and 1e-10.");
  }
  const startEpsilon = normalizedPath.start ? String(data.startEpsilon ?? "1") : (data.startEpsilon == null ? null : String(data.startEpsilon));
  if (startEpsilon !== null) {
    if (startEpsilon.length > 40 || !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(startEpsilon)) {
      throw new RangeError(normalizedPath.start ? "A custom segment starts at epsilon 1." : "Choose a starting epsilon above the target and at most 0.3.");
    }
    const startDecimal = new Decimal(startEpsilon);
    if (!startDecimal.isFinite() || !startDecimal.gt(targetDecimal) || (normalizedPath.start ? !startDecimal.eq(1) : startDecimal.gt(".3"))) {
      throw new RangeError(normalizedPath.start ? "A custom segment starts at epsilon 1." : "Choose a starting epsilon above the target and at most 0.3.");
    }
  }
  const samplesPerDecade = data.samplesPerDecade ?? 2;
  if (!Number.isInteger(samplesPerDecade) || samplesPerDecade < 1 || samplesPerDecade > 4) {
    throw new RangeError("Choose between 1 and 4 samples per decade.");
  }
  const tailRefinement = data.tailRefinement ?? 0;
  if (!Number.isInteger(tailRefinement) || tailRefinement < 0 || tailRefinement > 2) {
    throw new RangeError("Choose between 0 and 2 extra tail samples.");
  }
  return {
    counts: [...counts],
    path: normalizedPath,
    target,
    samplesPerDecade,
    tailRefinement,
    startEpsilon
  };
}

let traceRunning = false;
// Only complete, independently checked traces enter this small worker-local
// cache. Each exact request includes its marking, path, range, and density.
const completedTraceCache = new Map();

self.onmessage = async ({data}) => {
  if (!data || data.type !== "trace") return;
  const id = data.id;
  if (traceRunning) {
    self.postMessage({type: "error", id, message: "A trace is already running in this worker."});
    return;
  }
  let request;
  try {
    request = validateRequest(data);
  } catch (error) {
    self.postMessage({type: "error", id, message: error.message});
    return;
  }
  const cacheKey = JSON.stringify([request.counts,request.path.kind,request.path.index,
    new ViewDecimal(request.path.position).toString(),new ViewDecimal(request.target).toString(),
    request.samplesPerDecade,request.startEpsilon===null?null:new ViewDecimal(request.startEpsilon).toString(),
    request.path.start??null,request.tailRefinement]);
  if(completedTraceCache.has(cacheKey)){
    const result=completedTraceCache.get(cacheKey);
    completedTraceCache.delete(cacheKey);completedTraceCache.set(cacheKey,result);
    result.nodes.forEach((node,index)=>self.postMessage({type:'progress',id,node,completed:index+1,total:result.nodes.length}));
    self.postMessage({type:'done',id,result:{...result,fromCache:true}});return;
  }
  traceRunning = true;
  const started = performance.now();
  const nodes = [];
  const epsilons = request.startEpsilon === null
    ? approachMath.defaultApproachEpsilons(request.target, request.samplesPerDecade)
    : (() => {
      const start = new ViewDecimal(request.startEpsilon), target = new ViewDecimal(request.target);
      const out = [start.toString()];
      for (let k = 1; k < 100; k++) {
        const epsilon = start.times(ViewDecimal.pow(10, new ViewDecimal(-k).div(request.samplesPerDecade)));
        if (epsilon.lte(target)) break;
        out.push(epsilon.toPrecision(30));
      }
      out.push(target.toString());
      return out;
    })();
  // Retain every requested sample and both endpoints. Extra samples only
  // resolve the last existing logarithmic interval, where the cusp tail is
  // smallest. They use the same continuation and independent precision checks.
  if (request.tailRefinement && epsilons.length > 1) {
    const target = new ViewDecimal(epsilons.at(-1));
    const previous = new ViewDecimal(epsilons.at(-2));
    const ratio = target.div(previous);
    const extra = Array.from({length: request.tailRefinement}, (_, index) =>
      previous.times(ratio.pow(new ViewDecimal(index + 1).div(request.tailRefinement + 1))).toPrecision(60));
    epsilons.splice(epsilons.length - 1, 0, ...extra);
  }
  const signal = {get aborted() {return performance.now() - started > TRACE_TIME_BUDGET_MS;}};
  try {
    let result = await approachMath.solveApproachPath(request.counts, request.path, epsilons, {
      anchorData: attachmentAnchors,
      precisionDigits: 70,
      validationDigits: 100,
      signal,
      onProgress(node, completed, total) {
        decorateNode(node);
        nodes.push(node);
        self.postMessage({type: "progress", id, node, completed, total});
      }
    });
    if (result.status === "aborted") {
      result = {
        ...result,
        status: "partial",
        n: 1 + request.counts.reduce((sum, x) => sum + x, 0),
        completed: result.nodes.filter(node => node.status === "ok").length,
        total: epsilons.length,
        requestedEpsilons: epsilons,
        precisionDigits: 70,
        validationDigits: 100,
        stopReason: "time_budget",
        message: "The 45-second continuation budget was reached. Validated samples are retained; the target may not have been reached."
      };
    }
    result.nodes.forEach(decorateNode);
    result.tailRefinement = request.tailRefinement;
    if(result.status==='ok'){
      completedTraceCache.set(cacheKey,result);
      if(completedTraceCache.size>12)completedTraceCache.delete(completedTraceCache.keys().next().value);
    }
    self.postMessage({type: "done", id, result});
  } catch (error) {
    if (nodes.length) {
      self.postMessage({type: "done", id, result: {
        status: "partial", counts: request.counts,
        n: 1 + request.counts.reduce((sum, x) => sum + x, 0),
        path: request.path, nodes, completed: nodes.filter(node => node.status === "ok").length,
        total: epsilons.length, requestedEpsilons: epsilons,
        tailRefinement: request.tailRefinement,
        precisionDigits: 70, validationDigits: 100,
        elapsedMs: performance.now() - started,
        stopReason: "numerical_error", message: error.message
      }});
    } else {
      self.postMessage({type: "error", id, message: error.message});
    }
  } finally {
    traceRunning = false;
  }
};

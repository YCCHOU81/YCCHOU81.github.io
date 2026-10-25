/* Pure display coordinates. No elliptic solves and no changes to the input data.
 * Complex arguments accept {re, im} or [re, im]; returned points use {re, im}.
 * equal mode gives unit sample RMS only when the reference and mesh are resolved.
 */
(function (root) {
  "use strict";

  function complex(value) {
    if (value == null) return null;
    const re = Array.isArray(value) ? value[0] : value.re;
    const im = Array.isArray(value) ? value[1] : value.im;
    return Number.isFinite(re) && Number.isFinite(im) ? {re, im} : null;
  }

  function upper(value) {
    const z = complex(value);
    return z && z.im > 0 ? z : null;
  }

  function enumerate(n) {
    if (!Number.isInteger(n) || n < 1 || n > 21) {
      throw new RangeError("Choose an integer n from 1 to 21.");
    }
    const result = [];
    for (let a = 0; a < n; a += 1) {
      for (let b = 0; b < n - a; b += 1) {
        result.push([a, b, n - 1 - a - b]);
      }
    }
    return result;
  }

  function key(counts) {
    const abc = Array.isArray(counts) ? counts : [counts?.a, counts?.b, counts?.c];
    if (abc.length !== 3 || !abc.every(x => Number.isInteger(x) && x >= 0)) {
      throw new TypeError("An attachment is three nonnegative integer counts.");
    }
    return abc.join(",");
  }

  // ζ = (τ − τ*)/(τ − conjugate(τ*)). Scaling before multiplication
  // prevents overflow; differencing before scaling preserves nearby points.
  function disk(tau, reference) {
    const t = upper(tau);
    const p = upper(reference);
    if (!t || !p) return null;
    const delta = t.re - p.re;
    let dx, dy, y, v;
    if (Number.isFinite(delta)) {
      const magnitude = Math.max(Math.abs(delta), t.im, p.im);
      dx = delta / magnitude;
      dy = (t.im - p.im) / magnitude;
      y = t.im / magnitude;
      v = p.im / magnitude;
    } else {
      const magnitude = Math.max(Math.abs(t.re), Math.abs(p.re), t.im, p.im);
      dx = t.re / magnitude - p.re / magnitude;
      dy = (t.im - p.im) / magnitude;
      y = t.im / magnitude;
      v = p.im / magnitude;
    }
    const sum = y + v;
    const denominator = dx * dx + sum * sum;
    const re = (dx * dx + dy * sum) / denominator;
    const im = -2 * v * dx / denominator;
    return Number.isFinite(re) && Number.isFinite(im) ? {re, im} : null;
  }

  // Undo display scaling as well as the disk coordinate. Points at ζ = 1
  // represent infinity and return null rather than an invented finite value.
  function inverseDisk(value, reference, scale = 1) {
    const z = complex(value);
    const p = upper(reference);
    if (!z || !p || !Number.isFinite(scale) || scale <= 0) return null;
    const x = scale * z.re;
    const y = scale * z.im;
    if (!Number.isFinite(x) || !Number.isFinite(y)) return null;
    const c = 1 - x;
    const d = -y;
    const magnitude = Math.max(Math.abs(c), Math.abs(d));
    if (magnitude === 0) return null;
    const cs = c / magnitude;
    const ds = d / magnitude;
    const denominator = cs * cs + ds * ds;
    const a = (1 + x) / magnitude;
    const b = y / magnitude;
    const realRatio = (a * cs + b * ds) / denominator;
    const imaginaryRatio = (b * cs - a * ds) / denominator;
    const re = p.re - p.im * imaginaryRatio;
    const im = p.im * realRatio;
    return Number.isFinite(re) && Number.isFinite(im) ? {re, im} : null;
  }

  function bounds(points) {
    let minRe = Infinity, maxRe = -Infinity, minIm = Infinity, maxIm = -Infinity;
    let maxAbs = 0, count = 0;
    for (const z of points) {
      if (!z) continue;
      minRe = Math.min(minRe, z.re);
      maxRe = Math.max(maxRe, z.re);
      minIm = Math.min(minIm, z.im);
      maxIm = Math.max(maxIm, z.im);
      maxAbs = Math.max(maxAbs, Math.hypot(z.re, z.im));
      count += 1;
    }
    return count ? {minRe, maxRe, minIm, maxIm, maxAbs,
      width: maxRe - minRe, height: maxIm - minIm} : null;
  }

  function normalize(record, options = {}) {
    const mode = options.mode ?? "equal";
    if (mode !== "equal" && mode !== "hyperbolic") {
      throw new RangeError("Normalization mode must be equal or hyperbolic.");
    }
    const tau = Array.isArray(record?.tau) ? record.tau : [];
    const total = tau.length;
    const reference = upper(record?.reference);
    if (!reference) {
      return {points: tau.map(() => null), reference: null, scale: null,
        rms: null, normalizedRms: null, complete: false, resolved: 0, total,
        invalid: total, bounds: null, mode,
        normalizationLabel: "Reference point unresolved",
        message: "The common reference point has not resolved; no substitute center is used."};
    }
    const raw = tau.map(t => disk(t, reference));
    let squareSum = 0, resolved = 0;
    for (const z of raw) {
      if (!z) continue;
      squareSum += z.re * z.re + z.im * z.im;
      resolved += 1;
    }
    const rms = resolved ? Math.sqrt(squareSum / resolved) : null;
    const nonzeroRadius = rms !== null && rms > 0;
    const scale = mode === "equal" && nonzeroRadius ? rms : 1;
    const points = raw.map(z => z ? {re: z.re / scale, im: z.im / scale} : null);
    const complete = total > 0 && resolved === total && record.complete !== false &&
      !record.cancelled && (mode !== "equal" || nonzeroRadius);
    let normalizationLabel;
    let message = "";
    if (resolved === 0) {
      normalizationLabel = "Mesh unresolved";
      message = "No finite upper-half-plane sample has resolved.";
    } else if (mode === "equal" && !nonzeroRadius) {
      normalizationLabel = "Equal-size normalization unavailable";
      message = "The resolved samples have zero disk radius; an equal-size scale cannot be defined.";
    } else if (!complete) {
      normalizationLabel = mode === "equal" ? "Partial sample normalization" : "Partial hyperbolic view";
      message = "Some samples are missing; equal-size comparison is provisional until the common mesh resolves.";
    } else {
      normalizationLabel = mode === "equal" ? "Equal sample RMS" : "Hyperbolic disk";
    }
    return {points, reference, scale, rms,
      normalizedRms: rms === null ? null : rms / scale,
      complete, resolved, total, invalid: total - resolved,
      bounds: bounds(points), normalizationLabel, message, mode};
  }

  root.ComparisonMath = Object.freeze({enumerate, key, disk, inverseDisk, normalize, bounds});
})(typeof globalThis !== "undefined" ? globalThis : this);

/* An independent overlap layer: pixels belonging to just one attachment are
 * transparent. Pass one union Path2D per distinct attachment so that its
 * interior and cusp ribbons cannot count as overlaps with one another.
 * All paths and dimensions use CSS pixels; the destination may use a DPR
 * transform. This module only draws, and does not change numerical data.
 */
(function (root) {
  "use strict";
  const buffers = new WeakMap();
  const DEFAULT_COLOR = [22, 45, 64];
  function boundedStrength(value) {
    const number = Number(value);
    return Number.isFinite(number) ? Math.max(0, Math.min(1, number)) : 0;
  }
  function alphaForCount(count, strength) {
    if (!(count >= 2)) return 0;
    // Strength is the maximum extra opacity. More distinct attachments
    // darken their shared area smoothly, without changing a lone image.
    return boundedStrength(strength) * -Math.expm1(-(count - 1) / 2);
  }
  function makeCanvas(width, height) {
    const canvas = typeof root.OffscreenCanvas === "function"
      ? new root.OffscreenCanvas(width, height)
      : root.document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    return canvas;
  }
  function getBuffers(ctx, width, height, color) {
    let entry = buffers.get(ctx);
    if (!entry || entry.width !== width || entry.height !== height) {
      const countCanvas = makeCanvas(width, height);
      const shadeCanvas = makeCanvas(width, height);
      const countContext = countCanvas.getContext("2d", {
        alpha: false, willReadFrequently: true
      });
      const shadeContext = shadeCanvas.getContext("2d");
      entry = {width, height, countCanvas, shadeCanvas, countContext,
        shadeContext, pixels: shadeContext.createImageData(width, height)};
      buffers.set(ctx, entry);
    }
    const key = color.join(",");
    if (entry.color !== key) {
      const pixels = entry.pixels.data;
      for (let i = 0; i < pixels.length; i += 4) {
        pixels[i] = color[0];
        pixels[i + 1] = color[1];
        pixels[i + 2] = color[2];
      }
      entry.color = key;
    }
    return entry;
  }
  function paint(ctx, areaPaths, options = {}) {
    const width = Math.max(1, Math.ceil(Number(options.width) || 1));
    const height = Math.max(1, Math.ceil(Number(options.height) || 1));
    const strength = boundedStrength(options.strength ?? 0.45);
    const paths = Array.from(areaPaths || []).filter(Boolean);
    if (strength === 0 || paths.length < 2) {
      return {imageCount: paths.length, overlapPixels: 0, maxCoverage: 0};
    }
    const color = Array.isArray(options.color) && options.color.length === 3
      ? options.color.map(value => Math.max(0, Math.min(255, Math.round(Number(value) || 0))))
      : DEFAULT_COLOR;
    const entry = getBuffers(ctx, width, height, color);
    const count = entry.countContext;
    // An opaque black background is essential: using alpha=1/255 instead
    // would un-premultiply the readback and destroy the coverage counts.
    count.globalCompositeOperation = "copy";
    count.fillStyle = "rgb(0,0,0)";
    count.fillRect(0, 0, width, height);
    count.globalCompositeOperation = "lighter";
    count.fillStyle = "rgb(1,0,0)";
    count.globalAlpha = 1;
    for (const path of paths) count.fill(path, "nonzero");
    const counts = count.getImageData(0, 0, width, height).data;
    const pixels = entry.pixels.data;
    const alpha = new Uint8Array(256);
    for (let k = 2; k < alpha.length; k++) {
      alpha[k] = Math.round(255 * alphaForCount(k, strength));
    }
    let overlapPixels = 0, maxCoverage = 0;
    for (let i = 0; i < counts.length; i += 4) {
      const coverage = counts[i];
      pixels[i + 3] = alpha[coverage];
      if (coverage >= 2) overlapPixels++;
      if (coverage > maxCoverage) maxCoverage = coverage;
    }
    entry.shadeContext.putImageData(entry.pixels, 0, 0);
    ctx.save();
    ctx.globalAlpha = 1;
    ctx.globalCompositeOperation = "source-over";
    ctx.drawImage(entry.shadeCanvas, 0, 0, width, height);
    ctx.restore();
    return {imageCount: paths.length, overlapPixels, maxCoverage};
  }
  root.OverlapShading = Object.freeze({paint, alphaForCount});
})(typeof globalThis !== "undefined" ? globalThis : this);

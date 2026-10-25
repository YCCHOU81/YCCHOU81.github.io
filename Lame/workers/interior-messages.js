/* Existing classic-worker message protocol, point cache, and source validation. */
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

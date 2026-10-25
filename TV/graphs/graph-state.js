/* The graph lab expands small base-core pair decorations only. Never silently
   reinterpret a full Explore decoration as a different labelled graph. */
"use strict";
globalThis.GraphPageState = Object.freeze({
  read(params, seed, defaults) {
    const counts = {...defaults};
    const keys = seed === "double" ? ["qt", "qu", "ql", "qm"] :
      seed === "tetra" ? ["qAB", "qAC", "qAD", "qBC", "qBD", "qCD"] :
      ["bq0", "bq1", "bq2", "bq3", "bq4", "bq5"];
    const operation = params.get("decoration") || "repeated";
    const pairPreset = operation === "pair" && ["double", "tetra"].includes(seed);
    if (pairPreset) Object.assign(counts, DecorationModel.preset(seed, "pair"));
    let reason = !["seed", "repeated"].includes(operation) && !pairPreset ?
      "The requested decoration operation is outside this graph lab's base-core pair model." : "";
    if (operation === "repeated" && seed === "double" && ["A", "B"].some(end =>
      Number(params.get("count" + end) || 0) !== 0 || params.get("mode" + end) === "cross"))
      reason = "The requested split or sprout is outside this graph lab's base-core pair model.";
    if (operation === "repeated" && seed === "bouquet" && params.has("bouquetCore") && params.get("bouquetCore") !== "bouquet_B0_T0")
      reason = "The requested three-loop core is outside this graph lab's base-core pair model.";
    for (const key of keys) {
      if (operation === "seed" || pairPreset) continue;
      if (!params.has(key)) continue;
      const raw = params.get(key), value = Number(raw);
      if (!raw.trim() || !Number.isInteger(value) || value < 0 || value > 3)
        reason = "This expanded graph lab allows 0–3 pairs on each chain.";
      else counts[key] = value;
    }
    if (seed === "bouquet" && ["bq3", "bq4", "bq5"].some(key => counts[key] !== 0))
      reason = "The three-loop base core has only three pair controls; the requested extra chains need another core.";
    const hemispheres = {bouquet:1, double:2, tetra:4}[seed] +
      2 * keys.reduce((sum, key) => sum + counts[key], 0);
    if (hemispheres > 16) reason = "This expanded graph lab allows at most 16 downstairs hemispheres.";
    const rKey = params.has("r") ? "r" : "fromR", sKey = params.has("s") ? "s" : "fromS";
    const R = Number(params.get(rKey)), S = Number(params.get(sKey));
    const x = seed === "bouquet" ? .5 : 1, y = seed === "tetra" ? 1 : .5;
    const source = params.has(rKey) && params.has(sKey) && Number.isFinite(R) && Number.isFinite(S) &&
      R > 0 && S > 0 && R < x && S < y && R/x + S/y > 1 ? [R, S] : null;
    return {counts:reason ? {...defaults} : counts, source,
      notice:reason ? reason + " Showing the undecorated base core of the selected seed instead." : ""};
  }
});

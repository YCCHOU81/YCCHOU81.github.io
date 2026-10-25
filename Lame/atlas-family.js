/* Normalize incoming links before any page controller reads them. */
(() => {
  "use strict";
  window.AtlasFamily = Object.freeze({
    id: "lame",
    allowedSeeds: Object.freeze(["lame"]),
    defaultSeed: "lame",
    home: "Lame.html"
  });
  document.documentElement.dataset.atlasFamily = "lame";

  const url = new URL(location.href);
  const query = url.searchParams;
  const requested = query.get("seed");
  query.set("seed", "lame");
  if (query.get("decoration") !== "attachment") query.set("decoration", "attachment");

  // Older links may come from the combined TV atlas. Keep navigation and
  // display choices, but discard controls that do not belong to this folder.
  const tvControls = [
    "qt", "qu", "ql", "qm", "modeA", "modeB", "countA", "countB",
    "qAB", "qAC", "qAD", "qBC", "qBD", "qCD", "bouquetCore",
    "bq0", "bq1", "bq2", "bq3", "bq4", "bq5"
  ];
  for (const key of tvControls) query.delete(key);
  if (requested && requested !== "lame") {
    for (const key of ["r", "s", "fromR", "fromS"]) query.delete(key);
  }
  if (url.href !== location.href) history.replaceState(null, "", url);
})();

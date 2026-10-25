/* Normalize links before controllers read them. This is the single registry of
 * exposed TV families and decoration names; mathematical admissibility belongs
 * to DecorationModel/BouquetModel and is checked separately. */
(() => {
  "use strict";
  const operations = Object.freeze({
    double:Object.freeze(["seed", "pair", "split", "sprout", "repeated"]),
    tetra:Object.freeze(["seed", "pair", "repeated"]),
    bouquet:Object.freeze(["seed", "repeated"])
  });
  const family = Object.freeze({
    id:"tv", allowedSeeds:Object.freeze(["double", "tetra", "bouquet"]),
    defaultSeed:"double", home:"TV.html", operations
  });
  window.AtlasFamily = family;
  document.documentElement.dataset.atlasFamily = family.id;

  const url = new URL(location.href), query = url.searchParams;
  const requested = query.get("seed");
  const seed = family.allowedSeeds.includes(requested) ? requested : family.defaultSeed;
  query.set("seed",seed);
  if (!operations[seed].includes(query.get("decoration"))) query.set("decoration","repeated");
  if (requested && !family.allowedSeeds.includes(requested)) {
    for (const key of ["r", "s", "fromR", "fromS"]) query.delete(key);
  }
  if (url.href !== location.href) history.replaceState(null,"",url);
})();

/* This folder is an independent atlas. Normalize incoming links before controllers read them. */
(() => {
  "use strict";
  const family = {"id": "tv", "allowedSeeds": ["double", "tetra", "bouquet"], "defaultSeed": "double", "home": "TV.html"};
  window.AtlasFamily = Object.freeze({...family,allowedSeeds:Object.freeze(family.allowedSeeds)});
  document.documentElement.dataset.atlasFamily=family.id;
  const url=new URL(location.href),query=url.searchParams;
  const requested=query.get("seed"),seed=family.allowedSeeds.includes(requested)?requested:family.defaultSeed;
  query.set("seed",seed);
  const operations=seed==="lame"?["attachment"]:seed==="double"?["seed","pair","split","sprout","repeated"]:seed==="tetra"?["seed","pair","repeated"]:["seed","repeated"];
  if(!operations.includes(query.get("decoration")))query.set("decoration",seed==="lame"?"attachment":"repeated");
  if(family.id==="lame")for(const key of ["qt","qu","ql","qm","modeA","modeB","countA","countB","qAB","qAC","qAD","qBC","qBD","qCD","bouquetCore","bq0","bq1","bq2","bq3","bq4","bq5"])query.delete(key);
  if(requested&&!family.allowedSeeds.includes(requested))for(const key of ["r","s","fromR","fromS"])query.delete(key);
  if(url.href!==location.href)history.replaceState(null,"",url);
})();

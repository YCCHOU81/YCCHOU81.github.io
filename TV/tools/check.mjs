/* Optional maintenance checks. Node.js built-ins only; the website needs no build.
 * From the TV folder: node tools/check.mjs [--numerics]
 */
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import {fileURLToPath, pathToFileURL} from "node:url";
import {spawnSync} from "node:child_process";

const root = fileURLToPath(new URL("../",import.meta.url));
const list = directory=>fs.readdirSync(directory,{withFileTypes:true}).flatMap(entry=> {
  const filename = path.join(directory,entry.name);
  return entry.isDirectory() ? list(filename) : [filename];
});
const files = list(root);
let references = 0;
function localReference(filename,reference) {
  if (!reference || /^(?:[a-z]+:|\/\/|#)/i.test(reference)) return;
  const relative = decodeURIComponent(reference.split(/[?#]/)[0]);
  if (!relative) return;
  const target = path.resolve(path.dirname(filename),relative);
  assert(target===path.resolve(root)||target.startsWith(root),`${filename}: reference leaves TV: ${reference}`);
  assert(fs.existsSync(target),`${filename}: missing ${reference}`);
  references++;
}
for (const filename of files) {
  const extension = path.extname(filename);
  const source = fs.readFileSync(filename,"utf8");
  if ([".js",".mjs"].includes(extension)) {
    const checked = spawnSync(process.execPath,["--check",filename],{encoding:"utf8"});
    assert.equal(checked.status,0,checked.stderr);
    for (const match of source.matchAll(/(?:\bfrom\s*|\bimport\s*\(\s*|\bimport\s*|\bnew Worker\(\s*|\bnew URL\(\s*)["']([^"']+)["']/g)) {
      localReference(filename,match[1]);
    }
    for (const match of source.matchAll(/importScripts\(([^)]*)\)/g)) {
      for (const argument of match[1].matchAll(/["']([^"']+)["']/g)) localReference(filename,argument[1]);
    }
  }
  if (extension === ".html") {
    for (const match of source.matchAll(/\b(?:src|href)="([^"]+)"/g)) localReference(filename,match[1]);
    const ids=[...source.matchAll(/<[a-z][^>]*\bid="([^"]+)"/gi)].map(match=>match[1]);
    assert.equal(new Set(ids).size,ids.length,`${filename}: duplicate static IDs`);
    assert(!/<details\b[^>]*\sopen(?:[\s>])/i.test(source),`${filename}: help should start collapsed`);
  }
}

// These topology checks are independent of period solving and apply to each
// decoration kind. The three labelled graph faces and the exterior give F=4.
const scope={URLSearchParams};
scope.window=scope;
const context=vm.createContext(scope);
for (const name of ["decoration-model.js","bouquet-data.js","bouquet-model.js"]) {
  vm.runInContext(fs.readFileSync(path.join(root,name),"utf8"),context,{filename:name});
}
const operations={double:["seed","pair","split","sprout","repeated"],
  tetra:["seed","pair","repeated"],bouquet:["seed","repeated"]};
const mixed={
  double:"qt=1&qu=1&ql=1&qm=0&modeA=cross&countA=1&modeB=arm&countB=1",
  tetra:"qAB=2&qAC=1&qAD=0&qBC=1&qBD=0&qCD=0",
  bouquet:"bouquetCore=bouquet_B0_T15_24&bq0=1&bq3=1"
};
function configuration(seed,operation,query="") {
  const params=new URLSearchParams(query);
  if (seed === "bouquet") return scope.BouquetModel.make(operation,scope.BouquetModel.parse(operation,params));
  return scope.DecorationModel.make(seed,operation,scope.DecorationModel.parse(seed,operation,params));
}
let modelChecks=0;
for (const [seed,kinds] of Object.entries(operations)) {
  for (const operation of kinds) {
    const config=configuration(seed,operation,operation==="repeated"?mixed[seed]:"");
    const graph=config.graph;
    assert.equal(config.seed,seed);
    assert(config.weights.every(n=>Number.isInteger(n)&&n>=0));
    assert.equal(graph.vertices.length,config.weights.reduce((sum,n)=>sum+n,0));
    assert.equal(graph.faces.length,4);
    assert.equal(graph.vertices.length-graph.pairs.length+graph.faces.length,2);
    graph.faces.forEach((face,i)=>assert.equal(face.length,2*config.weights[i]+1));
    assert(config.shift.every(x=>Number.isInteger(2*x)));
    modelChecks++;
  }
}

let numericalChecks=0;
if (process.argv.includes("--numerics")) {
  const {handleRequest}=await import(pathToFileURL(path.join(root,"repeated-worker.mjs")));
  await assert.rejects(handleRequest({type:"unknown"}),/Unknown numerical request/);
  numericalChecks++;
  const sources={double:[.58,.31],tetra:[.58,.62],bouquet:[.29,.31]};
  for (const seed of Object.keys(operations)) {
    const config=configuration(seed,"repeated");
    const [R,S]=sources[seed];
    const point=await handleRequest({type:"point",caseId:config.id,caseConfig:config,R,S,budgetMs:30000});
    for (const branch of [point.base,point.result]) {
      assert.equal(branch?.status,"ok",`${seed}: ${point.message??branch?.message??"point unresolved"}`);
      const tau=branch.tauExact??[branch.tau.re,branch.tau.im];
      assert(tau.every(x=>Number.isFinite(Number(x)))&&Number(tau[1])>0);
    }
    numericalChecks++;
  }
}
console.log(JSON.stringify({passed:true,files:files.length,localReferences:references,
  seedFamilies:Object.keys(operations),modelChecks,numericalChecks},null,2));

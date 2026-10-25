/* Optional maintenance checks. Node.js built-ins only; no npm install or build.
 * Run: node tools/check.mjs [--numerics]
 */
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import {fileURLToPath} from "node:url";
import {spawnSync} from "node:child_process";
import {performance} from "node:perf_hooks";

const root = fileURLToPath(new URL("../", import.meta.url));
const files = directory => fs.readdirSync(directory, {withFileTypes:true}).flatMap(entry => {
  const filename = path.join(directory, entry.name);
  return entry.isDirectory() ? files(filename) : [filename];
});
const allFiles = files(root);
let references = 0;

function localReference(filename, reference) {
  if (!reference || /^(?:[a-z]+:|\/\/|#)/i.test(reference)) return;
  const relative = decodeURIComponent(reference.split(/[?#]/)[0]);
  if (!relative) return;
  const target = path.resolve(path.dirname(filename), relative);
  assert(target.startsWith(root), `${filename}: reference leaves the Lamé folder: ${reference}`);
  assert(fs.existsSync(target), `${filename}: missing ${reference}`);
  references++;
}

for (const filename of allFiles) {
  const extension = path.extname(filename);
  if (extension === ".js" || extension === ".mjs") {
    const syntax = spawnSync(process.execPath, ["--check", filename], {encoding:"utf8"});
    assert.equal(syntax.status, 0, syntax.stderr);
    const source = fs.readFileSync(filename,"utf8");
    for (const match of source.matchAll(/importScripts\(([^)]*)\)/g)) {
      for (const argument of match[1].matchAll(/["']([^"']+)["']/g)) localReference(filename,argument[1]);
    }
    for (const match of source.matchAll(/new Worker\(["']([^"']+)["']/g)) localReference(filename,match[1]);
  }
  if (extension === ".html") {
    const source = fs.readFileSync(filename,"utf8");
    for (const match of source.matchAll(/\b(?:src|href)="([^"]+)"/g)) localReference(filename,match[1]);
    const ids = [...source.matchAll(/<[a-z][^>]*\bid="([^"]+)"/gi)].map(match=>match[1]);
    assert.equal(new Set(ids).size,ids.length,`${filename}: duplicate static IDs`);
    assert(!/<details\b[^>]*\sopen(?:[\s>])/i.test(source),`${filename}: help should start collapsed`);
  }
}

// Workers have private global scopes. Imports resolve against their entry URL,
// not against the imported script, just as classic browser workers specify.
function loadWorker(entry) {
  const entryFile = path.join(root,entry);
  let messages = [];
  const scope = {console,performance,setTimeout,clearTimeout,
    postMessage:message=>messages.push(structuredClone(message))};
  scope.self = scope;
  const context = vm.createContext(scope);
  const run = filename=>vm.runInContext(fs.readFileSync(filename,"utf8"),context,{filename});
  scope.importScripts = (...urls)=>urls.forEach(url=> {
    const filename = path.resolve(path.dirname(entryFile),url.split(/[?#]/)[0]);
    assert(filename.startsWith(root),"Worker import leaves the Lamé folder");
    run(filename);
  });
  run(entryFile);
  return async data=> {
    messages = [];
    await scope.onmessage({data});
    assert(messages.length,"Worker returned no messages");
    return messages.at(-1);
  };
}

let numericalChecks = 0;
if (process.argv.includes("--numerics")) {
  const send = loadWorker("numerics-worker.js");
  const point = await send({type:"solve",id:1,counts:[1,0,0],r:.31,s:.29,inset:.04});
  assert.equal(point.type,"point");
  assert.equal(point.result.status,"ok");
  assert.equal(point.result.base.status,"ok");
  assert.equal(point.result.result.status,"ok");
  assert(Math.abs(point.result.base.tau.re-.5904953761444916)<1e-9);
  assert(Math.abs(point.result.base.tau.im-.7963943351098978)<1e-9);
  numericalChecks++;
  const invalid = await send({type:"solve",id:2,counts:[21,0,0],r:.31,s:.29,inset:.04});
  assert.equal(invalid.type,"pointError");
  numericalChecks++;
  const boundary = await send({type:"solve",id:3,counts:[1,0,0],r:.5,s:.29,inset:.04});
  assert.equal(boundary.type,"pointError");
  numericalChecks++;
  // Loading both entry families catches lexical/global collisions in imports.
  const trace = loadWorker("limits-worker.js");
  const rejected = await trace({type:"trace",id:4,counts:[21,0,0],path:{kind:"edge",index:0,position:".5"},target:"1e-3"});
  assert.equal(rejected.type,"error");
  numericalChecks++;
  loadWorker("overlap/normalized-worker.js");
  numericalChecks++;
}
console.log(JSON.stringify({passed:true,files:allFiles.length,localReferences:references,numericalChecks},null,2));

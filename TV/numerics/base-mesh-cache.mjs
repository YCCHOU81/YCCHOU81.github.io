// Reuse only checked canonical base periods. Decorated roots are always
// computed normally. Entries stay in the existing bounded page cache.
export const BASE_MESH_VERSION='tv-base-mesh-checked-1';
const prepared=new WeakMap(),families=new Set(['double','bouquet']);
const key=(r,s)=>Number(r)+'|'+Number(s);
export function checkedBaseBranch(branch,family){
 if(!families.has(family)||branch?.status!=='ok'||!Array.isArray(branch.tauExact)||branch.tauExact.length!==2||
    !branch.tauExact.every(x=>Number.isFinite(Number(x)))||!(Number(branch.tauExact[1])>0)||!Number.isFinite(Number(branch.precisionDigits))||!(Number(branch.precisionDigits)>=50))return false;
 const limit=family==='bouquet'?1e-12:1e-10;
 if(!['raw','zero','one','half'].every(name=>Array.isArray(branch.views?.[name])&&branch.views[name].length===2&&branch.views[name].every(x=>Number.isFinite(Number(x)))&&
     branch.viewErrors?.[name]!=null&&Number.isFinite(Number(branch.viewErrors[name]))&&Number(branch.viewErrors[name])>=0&&Number(branch.viewErrors[name])<=limit))return false;
 if(family==='bouquet'&&!(branch.endpointChecked===true&&Number(branch.lowRequestedDigits)>=18&&Number(branch.highRequestedDigits)>=28))return false;
 return true;
}
function prepare(artifact,family,inset){
 const empty=new Map();
 if(!artifact||typeof artifact!=='object'||artifact.version!==BASE_MESH_VERSION||artifact.family!==family||!families.has(family)||Number(artifact.inset)!==Number(inset))return empty;
 const n=Number(artifact.resolution),delta=Number(artifact.inset);
 if(!Number.isInteger(n)||n<3||n>24||!(delta>=.025&&delta<1/3)||!Array.isArray(artifact.points)||!Array.isArray(artifact.samples)||artifact.points.length!==artifact.samples.length||artifact.points.length>(n+1)*(n+2)/2)return empty;
 const allowed=new Set();
 for(let i=0;i<=n;i++)for(let j=0;j<=n-i;j++){
  const b=[delta+(1-3*delta)*i/n,delta+(1-3*delta)*j/n,delta+(1-3*delta)*(n-i-j)/n];
  allowed.add(key((family==='bouquet'?.5:1)*(b[0]+b[2]),.5*(b[1]+b[2])));
 }
 const entries=new Map();
 for(let i=0;i<artifact.points.length;i++){
  const p=artifact.points[i],source=key(p?.r,p?.s);
  if(!allowed.has(source)||entries.has(source)||!checkedBaseBranch(artifact.samples[i],family))return empty;
  entries.set(source,artifact.samples[i]);
 }
 return entries;
}
export function cachedBaseBranch(artifact,family,source,inset=.04){
 if(!artifact||typeof artifact!=='object')return null;
 let value=prepared.get(artifact);
 if(!value||value.family!==family||value.inset!==Number(inset)){
  value={family,inset:Number(inset),entries:prepare(artifact,family,inset)};prepared.set(artifact,value);
 }
 return value.entries.get(key(source[0],source[1]))??null;
}
const tables=new Map();
export async function verifiedBaseMesh(family,resolution,inset=.04,points=null){
 if(!families.has(family)||Number(inset)!==.04)return null;
 if(!tables.has(family))tables.set(family,(family==='double'
  ?import('./verified-base-double.mjs?v=20261024-tv-audit2')
  :import('./verified-base-bouquet.mjs?v=20261024-tv-audit2')).then(module=>Array.isArray(module.default?.meshes)?module.default.meshes:[]).catch(()=>[]));
 const all=await tables.get(family),n=Number(resolution);
 const candidates=[...all].sort((a,b)=>(Number(b.resolution)===n)-(Number(a.resolution)===n)||b.resolution-a.resolution);
 if(Array.isArray(points))candidates.sort((a,b)=>points.reduce((sum,p)=>sum+!!cachedBaseBranch(b,family,[p.r,p.s],inset)-!!cachedBaseBranch(a,family,[p.r,p.s],inset),0));
 for(const artifact of candidates)if(Array.isArray(artifact?.points)&&prepare(artifact,family,inset).size===artifact.points.length&&artifact.points.length)return artifact;
 return null;
}
export function baseMeshCovers(artifact,family,points,inset=.04){
 return !!artifact&&points.every(p=>cachedBaseBranch(artifact,family,[p.r,p.s],inset));
}

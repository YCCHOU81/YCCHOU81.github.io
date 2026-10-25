// Research-only canonical labelled tetrahedral ribbon graphs. Browser-compatible.
// Generated from the same primitive graph and alternating pair expansion as
// DecorationModel; names and arbitrary case ids are deliberately ignored.
export const CANONICAL_TETRA_VERSION='tetra-canonical-labelled-v1';
export const tetraKeys=['qAB','qAC','qAD','qBC','qBD','qCD'];
const primitive={"vertices":[["AB","AD","AC"],["BC","BD","BA"],["CA","CD","CB"],["DA","DB","DC"]],"face_representatives":["AB","BA","AC","CB"],"P":[1,1,1],"stable_edges":[{"id":"AB","chain":["AB"],"epsilon":1,"faces":[0,1],"weight_per_pair":[1,1,0,0],"shift_twice_per_pair":[0,1],"first_length":[1,0,0],"count":0,"multiplicity":1,"start":"AB","end":"BA"},{"id":"AD","chain":["AD"],"epsilon":1,"faces":[1,2],"weight_per_pair":[0,1,1,0],"shift_twice_per_pair":[1,1],"first_length":[0,0,1],"count":0,"multiplicity":1,"start":"AD","end":"DA"},{"id":"AC","chain":["AC"],"epsilon":1,"faces":[2,0],"weight_per_pair":[1,0,1,0],"shift_twice_per_pair":[1,0],"first_length":[0,1,0],"count":0,"multiplicity":1,"start":"AC","end":"CA"},{"id":"BC","chain":["BC"],"epsilon":1,"faces":[0,3],"weight_per_pair":[1,0,0,1],"shift_twice_per_pair":[1,1],"first_length":[0,0,1],"count":0,"multiplicity":1,"start":"BC","end":"CB"},{"id":"BD","chain":["BD"],"epsilon":1,"faces":[3,1],"weight_per_pair":[0,1,0,1],"shift_twice_per_pair":[1,0],"first_length":[0,1,0],"count":0,"multiplicity":1,"start":"BD","end":"DB"},{"id":"CD","chain":["CD"],"epsilon":1,"faces":[2,3],"weight_per_pair":[0,0,1,1],"shift_twice_per_pair":[0,1],"first_length":[1,0,0],"count":0,"multiplicity":1,"start":"CD","end":"DC"}]};
const json=x=>JSON.stringify(x),same=(a,b)=>json(a)===json(b),cache=new Map(),LIMIT=24;
function readCounts(counts){
 if(!counts||typeof counts!=='object'||Array.isArray(counts))throw Error('A full six-count tetrahedral object is required.');
 const out=Object.fromEntries(tetraKeys.map(k=>[k,counts[k]]));
 if(!Object.values(out).every(v=>Number.isSafeInteger(v)&&v>=0))throw Error('Each labelled tetrahedral count must be a nonnegative safe integer.');
 const total=Object.values(out).reduce((a,b)=>a+b,0);
 if(!Number.isSafeInteger(total)||4+2*total>10000)throw Error('The requested graph exceeds the finite canonical validation allocation limit.');
 return out;
}
function faceCycles(graph){
 const alpha=new Map(graph.pairs.flatMap(([a,b])=>[[a,b],[b,a]]));
 const sigma=new Map(graph.vertices.flatMap(v=>v.map((d,i)=>[d,v[(i+1)%v.length]])));
 const faces=[],seen=new Set();
 for(const pair of graph.pairs)for(const first of pair)if(!seen.has(first)){
  const f=[];let d=first;
  do{if(seen.has(d)||!sigma.has(alpha.get(d)))throw Error('Invalid labelled tetrahedral face.');seen.add(d);f.push(d);d=sigma.get(alpha.get(d));}while(d!==first);
  faces.push(f);
 }
 return graph.face_representatives.map(d=>faces.find(f=>f.includes(d)));
}
/** Fresh canonical expanded graph. Counts are never inferred from weights. */
export function canonicalGraph(input){
 const counts=readCounts(input),graph={vertices:primitive.vertices.map(v=>[...v]),pairs:[],lengths:{},length_coefficients:{},face_representatives:[...primitive.face_representatives],P:[...primitive.P],weights:[1,1,1,1],shift:[0,0],faces:[],stable_edges:[]};
 for(const [ei,original] of primitive.stable_edges.entries()){
  const {count,multiplicity,start:oldStart,end:oldEnd,...e}=original,q=counts['q'+e.id],m=1+2*q,start=e.id,end=e.id.split('').reverse().join('');
  for(let i=0;i<4;i++)graph.weights[i]+=q*e.weight_per_pair[i];
  for(let i=0;i<2;i++)graph.shift[i]=(graph.shift[i]*2+q*e.shift_twice_per_pair[i])%2/2;
  const ds=[start,...Array.from({length:m-1},(_,j)=>`qe${ei}_${j+1}p`)],Ds=[...Array.from({length:m-1},(_,j)=>`qe${ei}_${j}m`),end];
  for(let j=0;j<m;j++){
   const v=j%2?graph.P.map((p,i)=>p-e.first_length[i]):[...e.first_length];
   graph.pairs.push([ds[j],Ds[j]]);graph.length_coefficients[ds[j]]=v;graph.length_coefficients[Ds[j]]=[...v];
   graph.lengths[ds[j]]=v.map((x,i)=>x?`(${x})*${['x','y','z'][i]}`:null).filter(Boolean).join('+')||'0';
   if(j<m-1)graph.vertices.push([Ds[j],ds[j+1]]);
  }
  graph.stable_edges.push({...e,count:q,multiplicity:m,start,end});
 }
 graph.faces=faceCycles(graph);
 if(graph.faces.length!==4||graph.faces.some((f,i)=>!f||(f.length-1)/2!==graph.weights[i])||graph.vertices.length!==graph.weights.reduce((a,b)=>a+b,0))throw Error('Internal canonical tetrahedral graph degree check failed.');
 return graph;
}
function expected(input){
 const counts=readCounts(input),key=json(tetraKeys.map(k=>counts[k]));
 if(cache.has(key)){const old=cache.get(key);cache.delete(key);cache.set(key,old);return old;}
 const graph=canonicalGraph(counts),fields={};
 for(const f of ['vertices','pairs','face_representatives','weights','shift','P','stable_edges','faces','lengths'])fields[f]=json(graph[f]);
 const lengths=Object.entries(graph.length_coefficients).map(([k,v])=>[k,json(v)]),edgeCounts=Object.fromEntries(graph.stable_edges.map(e=>[e.id,e.count]));
 const value={counts,fields,lengths,edgeCounts,weights:[...graph.weights],shift:[...graph.shift],hemispheres:graph.vertices.length};
 cache.set(key,value);while(cache.size>LIMIT)cache.delete(cache.keys().next().value);return value;
}
const invalid=(reason,field)=>({valid:false,canonical:false,reason,field,version:CANONICAL_TETRA_VERSION});
/** Validate labels, cyclic orders, lengths, marking and all derived counts.
 * This is an exact combinatorial certificate, not a numerical root certificate.
 */
export function checkCanonicalTetra(config){
 try{
  if(config?.seed!=='tetra')return invalid('A tetrahedral configuration is required.','seed');
  const e=expected(config.counts),graph=config.graph;
  if(!graph||typeof graph!=='object')return invalid('A full labelled ribbon graph is required.','graph');
  for(const field of ['vertices','pairs','face_representatives','weights','shift','P','stable_edges'])if(json(graph[field])!==e.fields[field])return invalid('The labelled ribbon graph differs from the canonical six-count expansion: '+field+'.','graph.'+field);
  const lc=graph.length_coefficients;
  if(!lc||Object.keys(lc).length!==e.lengths.length||e.lengths.some(([k,v])=>json(lc[k])!==v))return invalid('The labelled spherical edge-length coefficients differ from the canonical expansion.','graph.length_coefficients');
  for(const field of ['faces','lengths'])if(graph[field]!==undefined&&json(graph[field])!==e.fields[field])return invalid('The derived graph '+field+' differ from the canonical expansion.','graph.'+field);
  for(const field of ['weights','shift'])if(!same(config[field],e[field]))return invalid('The ordered configuration '+field+' disagree with the full decorated graph.',field);
  if(config.core!==undefined&&config.core!=='tetra')return invalid('The marked core must be tetra.','core');
  if(config.top!==undefined&&config.top!==1)return invalid('The tetrahedral source triangle has top=1.','top');
  if(config.hemispheres!==undefined&&config.hemispheres!==e.hemispheres)return invalid('The hemisphere count disagrees with the expanded graph.','hemispheres');
  if(config.edgeCounts!==undefined&&(Object.keys(config.edgeCounts).length!==6||Object.entries(e.edgeCounts).some(([k,v])=>config.edgeCounts[k]!==v)))return invalid('The stable edge counts disagree with the six labelled decorations.','edgeCounts');
  const isConnector=tetraKeys.filter(k=>k!=='qAB').every(k=>e.counts[k]===0),connector=isConnector?{family:'tetra',q:e.counts.qAB,edge:'AB',isogenyLift:[0,0]}:undefined;
  if(config.connector!==undefined&&!same(config.connector,connector))return invalid('The connector shortcut metadata do not match this labelled graph.','connector');
  return{valid:true,canonical:true,version:CANONICAL_TETRA_VERSION,counts:{...e.counts},weights:[...e.weights],shift:[...e.shift],hemispheres:e.hemispheres,edgeCounts:{...e.edgeCounts},reason:'The full labelled ribbon graph matches the canonical six-count expansion.'};
 }catch(error){return invalid(String(error?.message??error),'counts');}
}
export const validateCanonicalTetra=checkCanonicalTetra;
export const clearCanonicalCache=()=>cache.clear();
export const canonicalCacheStats=()=>({version:CANONICAL_TETRA_VERSION,entries:cache.size,limit:LIMIT});

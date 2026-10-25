// Cone-centered marked Hodge mesh for decorated three-loop graphs.
// Cotangent weights use cone-centered per-face coordinates, normalized before
// areas. Angular phases are stored in two parts rather than recovered from
// globally rounded XY positions. The labelled topology is unchanged.
const dot=(a,b)=>{let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s};
const rational=x=>typeof x==='string'&&x.includes('/')?x.split('/').map(Number).reduce((a,b)=>a/b):Number(x);
const unitDelta=beta=>[-2*Math.sin(beta/2)**2,Math.sin(beta)];
const nearAngle=(a,b)=>{let x=a-b;while(x>Math.PI)x-=2*Math.PI;while(x< -Math.PI)x+=2*Math.PI;return x;};
function powerTriangle(points,order){let reference=null;return points.map(([x,y])=>{const radius=Math.hypot(x,y);if(!radius)return[0,0];let angle=Math.atan2(y,x);if(reference===null)reference=angle;else angle=reference+nearAngle(angle,reference);const r=Math.exp(Math.log(radius)/order),a=angle/order;return[r*Math.cos(a),r*Math.sin(a)];});}
export function hemisphereMesh(graph,seed,R,S,{sub=6,rings=4,grades=8,regularize=false}={}){
 const alpha=new Map(graph.pairs.flatMap(([a,b])=>[[a,b],[b,a]]));
 const cycles=graph.vertices.map(v=>seed==='double'?[...v].reverse():[...v]);
 const reps=graph.face_representatives.map(d=>seed==='double'?alpha.get(d):d),P=2*Math.PI;
 const values=seed==='bouquet'?[Math.PI*(1-2*R),Math.PI*(1-2*S),Math.PI*(2*R+2*S-1)]:seed==='double'?[P*(1-R),P*(R+2*S-1),Math.PI*(1-2*S)]:[P*(1-R),P*(1-S),P*(R+S-1)];
 const lens=new Map();for(const[a,b]of graph.pairs){const cf=graph.length_coefficients[a];const v=cf.reduce((s,c,j)=>s+rational(c)*values[j],0);if(!(v>0))throw Error('The graph has a nonpositive arc');lens.set(a,v);lens.set(b,v)}
 const xy=[],sheets=[],boundary=[],faces=[],arcs=new Map(),starts=new Map(),phaseParts=[];
 for(let sheet=0;sheet<2;sheet++)for(const cycle of cycles){
  const m=cycle.length*sub,off=xy.length,angles=[];let angle=0;
  for(const d of cycle){starts.set(sheet+':'+d,off+1+(rings-1)*m+angles.length);for(let j=0;j<sub;j++)angles.push(angle+lens.get(d)*j/sub);angle+=lens.get(d)}
  if(Math.abs(angle-P)>1e-8)throw Error('Hemisphere perimeter does not equal 2π');
  xy.push([0,0]);sheets.push(sheet);boundary.push(false);phaseParts.push([0,0]);
  for(let k=1;k<=rings;k++)for(const a of angles){xy.push([k/rings*Math.cos(a),k/rings*Math.sin(a)]);sheets.push(sheet);boundary.push(k===rings);phaseParts.push([a,0]);}
  const node=(k,j)=>off+1+(k-1)*m+(j%m+m)%m;
  for(let j=0;j<m;j++)faces.push([off,node(1,j),node(1,j+1)]);
  for(let k=1;k<rings;k++)for(let j=0;j<m;j++){const a=node(k,j),b=node(k,j+1),c=node(k+1,j),d=node(k+1,j+1);faces.push([a,c,d],[a,d,b])}
  cycle.forEach((d,k)=>arcs.set(sheet+':'+d,Array.from({length:sub+1},(_,j)=>node(rings,k*sub+j))));
 }
 const parent=Array.from({length:xy.length},(_,i)=>i),find=a=>{while(parent[a]!==a){parent[a]=parent[parent[a]];a=parent[a]}return a};
 for(const[a,b]of graph.pairs)for(let sh=0;sh<2;sh++){const aa=arcs.get(sh+':'+a),bb=arcs.get((1-sh)+':'+b);aa.forEach((v,j)=>parent[find(v)]=find(bb[sub-j]))}
 const ids=new Map(),g=parent.map((_,i)=>{const r=find(i);if(!ids.has(r))ids.set(r,ids.size);return ids.get(r)}),cones=reps.map(d=>g[starts.get('0:'+d)]),coneSet=new Set(cones),coneOrders=new Map(cones.map((c,i)=>[c,2*graph.weights[i]+1]));let nv=ids.size,ff=faces;
 // Each cone face has one local frame. Empty cone-data entries are faces
 // already away from all cone points and will never be subdivided again.
 let coords=faces.map(f=>f.map(v=>xy[v])),frames=faces.map(f=>{
  const k=f.findIndex(v=>coneSet.has(g[v]));if(k<0)return null;
  if(f.filter(v=>coneSet.has(g[v])).length!==1)throw Error('Insufficient cone subdivision');
  const c=f[k],theta=phaseParts[c][0],C=xy[c],cs=Math.cos(theta),sn=Math.sin(theta);
  const local=f.map(v=>{if(v===c)return[0,0];if(boundary[v])return unitDelta(nearAngle(phaseParts[v][0],theta));const dx=xy[v][0]-C[0],dy=xy[v][1]-C[1];return[cs*dx+sn*dy,-sn*dx+cs*dy]});
  return{c,theta,local,order:coneOrders.get(g[c]),beta:f.map(v=>boundary[v]?nearAngle(phaseParts[v][0],theta):null)};
 });
 coords=frames.map((frame,i)=>frame?frame.local:coords[i]);let orders=frames.map(frame=>frame?.order??1);
 for(let grade=0;grade<grades;grade++){
  const local=new Map(),global=new Map(),next=[],nextCoords=[],nextFrames=[],nextOrders=[];
  function midpoint(c,v,delta,beta,frame){
   const key=c<v?c+':'+v:v+':'+c;if(local.has(key))return local.get(key);
   const ga=g[c],gb=g[v],edge=ga<gb?ga+':'+gb:gb+':'+ga;if(!global.has(edge))global.set(edge,nv+global.size);
   const bound=boundary[c]&&boundary[v],newBeta=bound?beta/2:null,newDelta=bound?unitDelta(newBeta):delta.map(x=>x/2),cs=Math.cos(frame.theta),sn=Math.sin(frame.theta),C=xy[c];
   const X=[C[0]+cs*newDelta[0]-sn*newDelta[1],C[1]+sn*newDelta[0]+cs*newDelta[1]],i=xy.length;
   xy.push(X);sheets.push(sheets[c]);boundary.push(bound);g.push(global.get(edge));phaseParts.push([frame.theta,Math.atan2(newDelta[1],1+newDelta[0])]);
   const item={i,delta:newDelta,beta:newBeta};local.set(key,item);return item;
  }
  for(let fi=0;fi<ff.length;fi++){
   const f=ff[fi],frame=frames[fi];if(!frame){next.push(f);nextCoords.push(coords[fi]);nextFrames.push(null);nextOrders.push(orders[fi]);continue;}
   const k=f.indexOf(frame.c),c=f[k],a=f[(k+1)%3],b=f[(k+2)%3],da=frame.local[(k+1)%3],db=frame.local[(k+2)%3],ma=midpoint(c,a,da,frame.beta[(k+1)%3],frame),mb=midpoint(c,b,db,frame.beta[(k+2)%3],frame);
   const keep=[[0,0],ma.delta,mb.delta];
   next.push([c,ma.i,mb.i],[ma.i,a,b],[ma.i,b,mb.i]);
   nextCoords.push(keep,[ma.delta,da,db],[ma.delta,db,mb.delta]);
   nextFrames.push({c,theta:frame.theta,local:keep,order:frame.order,beta:[0,ma.beta,mb.beta]},null,null);nextOrders.push(frame.order,frame.order,frame.order);
  }
  nv+=global.size;ff=next;coords=nextCoords;frames=nextFrames;orders=nextOrders;
 }
 if(regularize)coords=coords.map((f,i)=>orders[i]>1?powerTriangle(f,orders[i]):f);
 return{xy,g,faces:ff,faceCoordinates:coords,phaseParts,cones,nv,sheets};
}
export function meshHealth(mesh){
 let bad=0,minArea=Infinity,maxCot=0;
 for(const f of mesh.faceCoordinates){const u=[f[0][0]-f[2][0],f[0][1]-f[2][1]],v=[f[1][0]-f[2][0],f[1][1]-f[2][1]],s=Math.max(...u.map(Math.abs),...v.map(Math.abs)),a=u.map(x=>x/s),b=v.map(x=>x/s),area=Math.abs(a[0]*b[1]-a[1]*b[0]);if(!(area>0)||!Number.isFinite(area)){bad++;continue;}minArea=Math.min(minArea,area);maxCot=Math.max(maxCot,Math.abs((a[0]*b[0]+a[1]*b[1])/area));}
 return{healthy:bad===0,bad,triangles:mesh.faces.length,vertices:mesh.nv,minNormalizedArea:minArea,maxCotangent:maxCot};
}

export function markedHodge(graph,seed,R=.58,S=seed==='double'?.31:.62,options={}){
 const begun=performance.now(),mesh=options.mesh??hemisphereMesh(graph,seed,R,S,options),{xy,g,faces,nv,sheets}=mesh,nf=faces.length;
 const edgeMap=new Map(),edges=[],w=[],edgefaces=[],faceedges=[],facesigns=[],localRep=[];
 for(let fi=0;fi<nf;fi++){
  const f=faces[fi],coords=mesh.faceCoordinates?.[fi]??f.map(v=>xy[v]),es=[],ss=[];
  for(const[aa,bb,cc]of [[0,1,2],[1,2,0],[2,0,1]]){
   const ia=f[aa],ib=f[bb],ic=f[cc],a=g[ia],b=g[ib],key=a<b?a+':'+b:b+':'+a;
   let ei=edgeMap.get(key);if(ei===undefined){ei=edges.length;edgeMap.set(key,ei);edges.push(a<b?[a,b]:[b,a]);w.push(0);edgefaces.push([]);localRep.push([ia,ib,a<b?1:-1])}
   const u0=[coords[aa][0]-coords[cc][0],coords[aa][1]-coords[cc][1]],v0=[coords[bb][0]-coords[cc][0],coords[bb][1]-coords[cc][1]],scale=Math.max(...u0.map(Math.abs),...v0.map(Math.abs)),u=u0.map(x=>x/scale),v=v0.map(x=>x/scale),area=Math.abs(u[0]*v[1]-u[1]*v[0]);
   w[ei]+=.5*(u[0]*v[0]+u[1]*v[1])/area;edgefaces[ei].push(fi);es.push(ei);ss.push(a<b?1:-1);
  }faceedges.push(es);facesigns.push(ss);
 }
 const ne=edges.length;if(nv-ne+nf!==0||edgefaces.some(f=>f.length!==2))throw Error('The doubled graph mesh is not a torus');
 const adj=Array.from({length:nv},()=>[]);edges.forEach(([a,b],e)=>{adj[a].push([b,e]);adj[b].push([a,e])});
 const parent=new Int32Array(nv).fill(-1),parentEdge=new Int32Array(nv),tree=new Set(),queue=[0];parent[0]=0;
 for(let k=0;k<queue.length;k++){const a=queue[k];for(const[b,e]of adj[a])if(parent[b]<0){parent[b]=a;parentEdge[b]=e;queue.push(b);tree.add(e)}}
 const dadj=Array.from({length:nf},()=>[]);edgefaces.forEach(([f,h],e)=>{if(!tree.has(e)){dadj[f].push([h,e]);dadj[h].push([f,e])}});
 const par=new Int32Array(nf).fill(-1),pe=new Int32Array(nf),order=[0],cotree=new Set();par[0]=0;
 for(let k=0;k<order.length;k++){const f=order[k];for(const[h,e]of dadj[f])if(par[h]<0){par[h]=f;pe[h]=e;order.push(h);cotree.add(e)}}
 const generators=edges.map((_,e)=>e).filter(e=>!tree.has(e)&&!cotree.has(e));if(generators.length!==2)throw Error('Could not recover the two period generators');
 const q=[new Float64Array(ne),new Float64Array(ne)];generators.forEach((e,k)=>q[k][e]=1);
 for(let n=order.length-1;n>0;n--){const f=order[n],e=pe[f],es=faceedges[f],ss=facesigns[f],loc=es.indexOf(e);for(let k=0;k<2;k++){let v=0;for(let j=0;j<3;j++)v+=q[k][es[j]]*ss[j];q[k][e]=-v/ss[loc]}}
 let intersection=0;for(let f=0;f<nf;f++){const es=faceedges[f],ss=facesigns[f];intersection+=.5*(q[0][es[0]]*ss[0]*(-q[1][es[2]]*ss[2])-(-q[0][es[2]]*ss[2])*q[1][es[0]]*ss[0])}
 if(Math.abs(Math.abs(intersection)-1)>1e-8)throw Error('The integral period intersection is invalid');if(intersection<0)for(let e=0;e<ne;e++)q[1][e]*=-1;
 const diag=new Float64Array(nv);edges.forEach(([a,b],e)=>{diag[a]+=w[e];diag[b]+=w[e]});diag[0]=1;
 const op=(x,out)=>{out.fill(0);for(let e=0;e<ne;e++){const[a,b]=edges[e],z=w[e]*(x[a]-x[b]);out[a]+=z;out[b]-=z}out[0]=x[0]};
 const U=[],iterations=[];
 for(let k=0;k<2;k++){
  const rhs=new Float64Array(nv);edges.forEach(([a,b],e)=>{const v=w[e]*q[k][e];rhs[a]+=v;rhs[b]-=v});rhs[0]=0;
  const u=new Float64Array(nv),r=rhs.slice(),z=r.map((v,i)=>v/diag[i]),p=z.slice(),Ap=new Float64Array(nv);let rz=dot(r,z),r0=Math.sqrt(dot(r,r)),it=0;
  for(;it<12000;it++){op(p,Ap);const alpha=rz/dot(p,Ap);for(let i=0;i<nv;i++){u[i]+=alpha*p[i];r[i]-=alpha*Ap[i]}if(Math.sqrt(dot(r,r))<2e-10*r0)break;for(let i=0;i<nv;i++)z[i]=r[i]/diag[i];const nr=dot(r,z),beta=nr/rz;for(let i=0;i<nv;i++)p[i]=z[i]+beta*p[i];rz=nr}
  if(it===12000)throw Error('Hodge initializer did not converge');U.push(u);iterations.push(it+1);
 }
 let G00=0,G01=0,G11=0;edges.forEach(([a,b],e)=>{const x=q[0][e]+U[0][b]-U[0][a],y=q[1][e]+U[1][b]-U[1][a];G00+=w[e]*x*x;G01+=w[e]*x*y;G11+=w[e]*y*y});
 function phases(R1,S1){const m=hemisphereMesh(graph,seed,R1,S1,options),parts=m.phaseParts;const phase=localRep.map(([a,b,sign])=>{const v=((parts[b][0]-parts[a][0])+(parts[b][1]-parts[a][1]))/(2*Math.PI)*(1-2*m.sheets[a]);return(v-Math.round(v))*sign}),pot=new Float64Array(nv);for(const v of queue.slice(1)){const pr=parent[v],e=parentEdge[v];pot[v]=pot[pr]+(pr<v?phase[e]:-phase[e])}return generators.map((e,k)=>{const[a,b]=edges[e];return(phase[e]+pot[a]-pot[b])*q[k][e]})}
 const epsilon=1e-5,deriv=[[epsilon,0],[0,epsilon]].map(([dr,ds])=>{const a=phases(R+dr,S+ds),b=phases(R-dr,S-ds);return a.map((v,i)=>{const d=v-b[i];return(d-Math.round(d))/(2*epsilon)})}),mark=[[-deriv[0][0]/2,deriv[1][0]/2],[-deriv[0][1]/2,deriv[1][1]/2]],M=mark.map(row=>row.map(Math.round));
 const error=Math.max(...mark.flatMap((v,i)=>v.map((x,j)=>Math.abs(x-M[i][j]))));if(error>1e-5||Math.abs(M[0][0]*M[1][1]-M[0][1]*M[1][0])!==1)throw Error('Developing phases did not resolve the integral marking');
 const[a,b]=M[0],[c,d]=M[1],H00=a*a*G00+2*a*c*G01+c*c*G11,H01=a*b*G00+(a*d+b*c)*G01+c*d*G11,H11=b*b*G00+2*b*d*G01+d*d*G11,det=H00*H11-H01*H01;
 if(!(det>0))throw Error('The Hodge period matrix is not positive');
 return{tau:[-H01/H11,Math.sqrt(det)/H11],marking:M,markError:error,vertices:nv,iterations,milliseconds:performance.now()-begun};
}

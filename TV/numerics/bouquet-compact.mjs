// Optional binary64 prepolish for the COMPACT source only. Its output is never
// accepted as a displayed root; arbitrary-precision correction follows it.
import * as F from './generic-fast.mjs';
const charts=[[1,0,0,1],[1,-1,1,0],[0,1,-1,1]],perms=[[0,1,2,3],[0,2,3,1],[0,3,1,2]],clock=()=>performance.now();
const mob=(z,[a,b,c,d])=>F.div(F.add(F.mul(a,z),b),F.add(F.mul(c,z),d));
const finite=z=>z.every(Number.isFinite);
function frozen(config,source,initial,edge,deadline,counter){
 const[a,b,c,d]=charts[edge],r=Number(source[1])+config.shift[1],s=Number(source[0])+config.shift[0],red=F.reduceModular(perms[edge].map(i=>config.weights[i]),initial,a*r-b*s,d*s-c*r);
 const get=w=>{if(clock()>deadline)throw Error('Compact prepolish deadline reached');counter.evaluations++;return F.tvMatrix(red.weights,mob(w,red.matrix),red.r,red.s,{reduce:false}).M;};
 const raw=get(initial),N=raw.length,rows=Array.from({length:N},(_,i)=>i),cols=[...rows],rs=raw.map(row=>Math.max(...row.map(F.abs))),M=raw.map((row,i)=>row.map(z=>F.mul(1/rs[i],z))),cs=Array.from({length:N},(_,j)=>Math.max(...M.map(row=>F.abs(row[j]))));
 if([...rs,...cs].some(x=>!Number.isFinite(x)||!x))throw Error('Compact prepolish scaling is unresolved');
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=F.mul(1/cs[j],M[i][j]);
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=F.abs(M[j][j]);for(let i=j;i<N;i++)for(let k=j;k<N;k++){const n=F.abs(M[i][k]);if(n>best){p=i;q=k;best=n;}}
  if(!(best>0))throw Error('Compact prepolish minor is unresolved');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]];}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]];}
  for(let i=j+1;i<N;i++){const t=F.div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=F.sub(M[i][k],F.mul(t,M[j][k]));}
 }
 let first=M[N-1][N-1];
 return w=>{
  if(first&&F.abs(F.sub(w,initial))===0){const v=first;first=null;return v;}
  const raw=get(w),A=rows.map(i=>cols.map(j=>F.mul(1/(rs[i]*cs[j]),raw[i][j])));
  for(let j=0;j<N-1;j++){if(!F.abs(A[j][j]))throw Error('Compact prepolish pivot is unresolved');for(let i=j+1;i<N;i++){const t=F.div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=F.sub(A[i][k],F.mul(t,A[j][k]));}}
  const v=A[N-1][N-1];if(!finite(v))throw Error('Compact prepolish residual is nonfinite');return v;
 };
}
export function compactPrepolish({config,source,guess,edge=0,radius=.08,maxStep=.025,deadline=Infinity,maxIterations=14}={}){
 const begun=clock(),initial=guess.map(Number),counter={evaluations:0};let z=initial,lastStep=null;
 const finish=(converged,message)=>({converged,w:z.map(String),estimatedError:lastStep,milliseconds:clock()-begun,...counter,message,internalOnly:true});
 try{
  if(!finite(z)||z[1]<=.00001)throw Error('Unsafe compact prepolish initializer');
  const f=frozen(config,source,initial,edge,deadline,counter),h=[1e-5*(1+F.abs(z)),2e-5*(1+F.abs(z))];let value=f(z),derivative=F.div(F.sub(f(F.add(z,h)),value),h);
  for(let i=0;i<maxIterations;i++){
   if(clock()>deadline)throw Error('Compact prepolish deadline reached');
   if(!F.abs(value)||!F.abs(derivative))throw Error('Compact prepolish cannot resolve the local residual');
   let step=F.div(value,derivative),length=F.abs(step);lastStep=length;if(!Number.isFinite(length))throw Error('Compact prepolish step is nonfinite');
   if(length<8e-8)return finish(true,'Approximate compact prepolish; arbitrary-precision correction required');
   if(length>maxStep)step=F.mul(maxStep/length,step);let accepted=false,factor=1,next,nextValue;
   for(let k=0;k<8;k++){
    next=F.sub(z,F.mul(factor,step));
    if(finite(next)&&next[1]>.00001&&F.abs(F.sub(next,initial))<=radius){try{nextValue=f(next);if(F.abs(nextValue)<F.abs(value)){accepted=true;break;}}catch{}}
    factor/=2;
   }
   if(!accepted)throw Error('No decreasing compact prepolish step');derivative=F.div(F.sub(nextValue,value),F.sub(next,z));z=next;value=nextValue;
  }
  return finish(false,'Compact prepolish iteration budget reached');
 }catch(error){return finish(false,String(error.message??error));}
}

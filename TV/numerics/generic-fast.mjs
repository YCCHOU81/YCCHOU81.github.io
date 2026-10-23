// Standalone Bloch-jet determinant. No website imports or mutations.
export const C=x=>typeof x==='number'?[x,0]:x;
export const add=(...xs)=>xs.reduce((a,b)=>{b=C(b);return[a[0]+b[0],a[1]+b[1]]},[0,0]);
export const neg=a=>[-a[0],-a[1]];
export const sub=(a,b)=>add(a,neg(C(b)));
export const mul=(a,b)=>{a=C(a);b=C(b);return[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]]};
export const div=(a,b)=>{a=C(a);b=C(b);const d=b[0]**2+b[1]**2;return[(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d]};
export const abs=a=>Math.hypot(...a);
export const exp=a=>{const r=Math.exp(a[0]);return[r*Math.cos(a[1]),r*Math.sin(a[1])]};

export function reduceModular(weights,tau,r,s){
 weights=[...weights];tau=[...tau];let matrix=[1,0,0,1];
 for(let steps=0;steps<60;steps++){
  const m=Math.round(tau[0]);
  if(m){
   tau[0]-=m;r+=m*s;
   if(Math.abs(m)%2)weights=[weights[0],weights[1],weights[3],weights[2]];
   matrix=[matrix[0]-m*matrix[2],matrix[1]-m*matrix[3],matrix[2],matrix[3]];
  }
  if(abs(tau)>=1-1e-13)break;
  tau=div(-1,tau);[r,s]=[s,-r];weights=[weights[0],weights[2],weights[1],weights[3]];
  matrix=[-matrix[2],-matrix[3],matrix[0],matrix[1]];
 }
 r-=Math.round(r);s-=Math.round(s);
 return {weights,tau,r,s,matrix};
}

// Taylor coefficients in v=2*pi*z; theta common q^(1/4) omitted.
function thetaJets(z,tau,N){
 z=C(z);const count=Math.max(4,Math.ceil(Math.sqrt(42/(Math.PI*tau[1]))+Math.abs(z[1])/tau[1])+3);
 const out=Array.from({length:N+1},()=>[0,0]);
 for(let n=0;n<count;n++){
  const k=n+.5,q=exp(mul([0,Math.PI*n*(n+1)],tau));
  const ep=exp(mul([0,2*Math.PI*k],z)),em=div(1,ep),factor=mul([0,-((-1)**n)],q);
  let fp=[1,0],fm=[1,0];
  for(let d=0;d<=N;d++){
   out[d]=add(out[d],mul(factor,sub(mul(fp,ep),mul(fm,em))));
   fp=mul(fp,[0,k/(d+1)]);fm=mul(fm,[0,-k/(d+1)]);
  }
 }
 return out;
}
function polyDiv(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=a[k]??[0,0];
  for(let j=1;j<=Math.min(k,b.length-1);j++)v=sub(v,mul(b[j],out[k-j]));
  out.push(div(v,b[0]));
 }
 return out;
}
function polyMul(a,b,N){
 const out=[];
 for(let k=0;k<=N;k++){
  let v=[0,0];for(let j=Math.max(0,k-b.length+1);j<=Math.min(k,a.length-1);j++)v=add(v,mul(a[j],b[k-j]));
  out.push(v);
 }
 return out;
}
function choose(n,k){let v=1;for(let i=1;i<=k;i++)v=v*(n-k+i)/i;return v}

export function tvMatrix(weights,tau,r,s,{reduce=true}={}){
 if(!(tau[1]>0))throw Error('Upper half-plane required');
 if(reduce)({weights,tau,r,s}=reduceModular(weights,tau,r,s));
 else {r-=Math.round(r);s-=Math.round(s)}
 const sites=[];
 for(let i=0;i<4;i++)for(let h=weights[i]-1;h>=0;h-=2)sites.push([i,h]);
 const degree=2*Math.max(0,...sites.map(x=>x[1])),N=degree+1;
 const alpha=add(r,mul(s,tau)),zero=thetaJets(0,tau,N+1),num=thetaJets(alpha,tau,N),den=zero.slice(1);
 const expjet=[[1,0]];
 for(let d=1;d<=N;d++)expjet.push(mul(expjet[d-1],[0,s/d]));
 const factor=div(zero[1],num[0]),diag=polyMul(expjet,polyDiv(num,den,N),N).map(x=>mul(x,factor));
 const bits=[[0,0],[1,0],[0,1],[1,1]],points=[[0,0],[.5,0],mul(.5,tau),mul(.5,add(1,tau))],jets=new Map();
 const M=sites.map(([i,l])=>sites.map(([j,h])=>{
  const d=l+h;
  if(i===j)return mul(choose(d,l),diag[d+1]);
  const du=bits[i][0]-bits[j][0],dv=bits[i][1]-bits[j][1],hu=(du+2)%2,hv=(dv+2)%2;
  const key=hu+2*hv,pointIndex=key,periodU=(du-hu)/2,periodV=(dv-hv)/2;
  if(!jets.has(key)){
   const z=points[pointIndex],scalar=mul(exp(mul([0,2*Math.PI*s],z)),factor);
   const jet=polyMul(expjet,polyDiv(thetaJets(add(z,alpha),tau,degree),thetaJets(z,tau,degree),degree),degree).map(x=>mul(x,scalar));
   jets.set(key,jet);
  }
  const periodFactor=exp([0,2*Math.PI*(s*periodU-r*periodV)]);
  return mul(choose(d,l),mul(periodFactor,jets.get(key)[d]));
 }));
 return {M,sites,weights,tau,r,s};
}
export function determinant(weights,tau,r,s,options={}){
 // This shortcut selects a square-root character lift. It is NOT the full
 // zero locus of the larger-period equation. Enable only after the graph
 // branch and its lift have been calibrated; full generic is the default.
 if(options.isogeny && weights[0]===weights[1] && weights[2]===weights[3]){
  const [u,v]=options.isogenyLift??[0,0];
  if(weights[0]===weights[2])return determinant([weights[0],0,0,0],tau,(r+u)/2,(s+v)/2,{...options,isogeny:false});
  if(weights[2]===0)return determinant([weights[0],0,0,0],mul(2,tau),r,(s+v)/2,{...options,isogeny:false});
  return determinant([weights[0],weights[2],0,0],div(-1,mul(2,tau)),(s+v)/2,-r,{...options,isogeny:false});
 }
 const {M}=tvMatrix(weights,tau,r,s,options),N=M.length;
 let product=[1,0],minRelativePivot=1;
 // Row normalization keeps intermediate products in the double range.
 for(let i=0;i<N;i++){
  // A 1x1 row is the residual itself: dividing by its absolute value would
  // erase its zero and make secant iteration impossible.
  const scale=N===1?1:Math.max(...M[i].map(abs));
  if(!scale||!Number.isFinite(scale))return {value:[NaN,NaN],minRelativePivot:0,dimension:N};
  for(let j=0;j<N;j++)M[i][j]=mul(1/scale,M[i][j]);
 }
 for(let j=0;j<N;j++){
  let p=j;for(let i=j+1;i<N;i++)if(abs(M[i][j])>abs(M[p][j]))p=i;
  const pivot=abs(M[p][j]);minRelativePivot=Math.min(minRelativePivot,pivot);
  if(pivot<1e-300)return {value:[0,0],minRelativePivot:0,dimension:N};
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];product=neg(product)}
  product=mul(product,M[j][j]);
  for(let i=j+1;i<N;i++){
   const ratio=div(M[i][j],M[j][j]);
   for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]));
  }
 }
 return {value:product,minRelativePivot,dimension:N};
}

export function solve(weights,r,s,guess,{tol=1e-11,maxsteps=32,isogeny=false,isogenyLift=[0,0]}={}){
 let a=guess,b=add(guess,[1e-6,2e-6]),fa=determinant(weights,a,r,s,{isogeny,isogenyLift}).value,fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 for(let i=0;i<maxsteps;i++){
  const denominator=sub(fb,fa);if(abs(denominator)===0)break;
  const step=div(mul(fb,sub(b,a)),denominator);
  if(abs(step)<tol)return {tau:b,iterations:i+2,residual:abs(fb)};
  a=b;fa=fb;b=sub(b,step);
  if(!Number.isFinite(abs(b))||b[1]<=.005)throw Error('Unsafe iterate');
  fb=determinant(weights,b,r,s,{isogeny,isogenyLift}).value;
 }
 throw Error('No convergence');
}

// Fixed row/column equilibration and complete-pivot Schur residual.
// Unlike renormalizing a determinant at every iterate, this preserves a
// holomorphic residual and resolves the formerly unstable count-32 anchor.
export function fixedSchur(E,get,initial){
 const abs=E.abs,div=E.div,mul=E.mul,sub=E.sub;
 const ref=get(initial),N=ref.length;if(!N)throw Error('Empty elliptic matching matrix');
 const greater=(a,b)=>a?.gt?a.gt(b):a>b;
 const maximum=a=>a.reduce((m,x)=>greater(abs(x),m)?abs(x):m,abs(a[0]));
 const scales=ref.map(maximum),M=ref.map((row,i)=>row.map(x=>div(x,scales[i])));
 const cs=Array.from({length:N},(_,j)=>maximum(M.map(row=>row[j])));
 for(let i=0;i<N;i++)for(let j=0;j<N;j++)M[i][j]=div(M[i][j],cs[j]);
 const rows=Array.from({length:N},(_,i)=>i),cols=[...rows];
 for(let j=0;j<N-1;j++){
  let p=j,q=j,best=abs(M[j][j]);
  for(let i=j;i<N;i++)for(let k=j;k<N;k++){const v=abs(M[i][k]);if(greater(v,best)){p=i;q=k;best=v}}
  if(!greater(best,0))throw Error('The Schur complement has a singular fixed pivot');
  if(p!==j){[M[p],M[j]]=[M[j],M[p]];[rows[p],rows[j]]=[rows[j],rows[p]]}
  if(q!==j){for(let i=0;i<N;i++)[M[i][q],M[i][j]]=[M[i][j],M[i][q]];[cols[q],cols[j]]=[cols[j],cols[q]]}
  for(let i=j+1;i<N;i++){const ratio=div(M[i][j],M[j][j]);for(let k=j+1;k<N;k++)M[i][k]=sub(M[i][k],mul(ratio,M[j][k]))}
 }
 return{dimension:N,evaluate:z=>{
  const raw=get(z),A=rows.map(i=>cols.map(j=>div(div(raw[i][j],scales[i]),cs[j])));
  for(let j=0;j<N-1;j++)for(let i=j+1;i<N;i++){const ratio=div(A[i][j],A[j][j]);for(let k=j+1;k<N;k++)A[i][k]=sub(A[i][k],mul(ratio,A[j][k]))}
  return A[N-1][N-1];
 }};
}

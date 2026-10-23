import{createEvaluator,connectorData}from'./decimal-generic.mjs';

export function leadingGuess(E,data){
 const{D,PI}=E,{n,edge,family,delta,x0}=data,tetra=family==='tetra';
 let bin=1n;for(let j=1;j<=n;j++)bin=bin*BigInt(n+j)/BigInt(j);
 const b=new D(bin.toString());
 const K=tetra?[new D(8*n*(n+1)).sqrt(),b.times(2),b.times(2)][edge]:[new D(2*n*(n+1)),b.times(n),b.times(2*n)][edge];
 const alpha=tetra?[new D(1),new D(2).div(n+1),new D(2).div(n+1)][edge]:[new D('.5'),new D(2).div(n),new D(1).div(2*n)][edge];
 return[x0,alpha.div(PI).times(K.div(delta).ln())];
}

export function runCase(spec){
 const {precision=80,guess,digits=precision-25}=spec,E=createEvaluator(precision),data=connectorData(E,spec);
 let evaluations=0;
 const f=w=>{evaluations++;return E.determinant(data.weights,E.inverse(w,data.chart),data.r,data.s,{isogeny:true}).value};
 const started=performance.now();
 try{
  const initial=guess?E.c(guess):leadingGuess(E,data);
  const got=E.solveFunction(f,initial,{digits,maxsteps:35,guard:w=>w[1].gt(0)&&w[0].minus(initial[0]).abs().lt('.4')});
  const correction=got.root[0].minus(data.x0);
  return{...spec,guess:E.strs(initial),w:E.strs(got.root),x0:data.x0.toString(),horizontalCorrection:correction.toString(),delta:data.delta.toString(),characters:[data.r.toString(),data.s.toString()],iterations:got.iterations,evaluations,lastStep:E.strs(got.lastStep),milliseconds:performance.now()-started,success:true};
 }catch(error){return{...spec,error:String(error),evaluations,milliseconds:performance.now()-started,success:false}}
}

if(typeof process!=='undefined'&&process.argv[1]?.endsWith('decimal-case.mjs')){
 const fs=await import('node:fs');
 const spec=JSON.parse(process.argv[2]??'{}'),result=runCase(spec);
 console.log(JSON.stringify(result));
 if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2));
}

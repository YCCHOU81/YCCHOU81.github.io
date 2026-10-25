// The same bounded corrector and stopping tests, with fixed-point arithmetic
// for mixed tetrahedral matrices. Geometry selection remains independent.
import {correct as decimalCorrect} from './double-corrector.mjs?v=20261024-tv-audit2';
import {createEvaluator} from './bouquet-fixed.mjs?v=20261024-tv-audit2';

const evaluators=new Map(),now=()=>performance.now();
// Precision here is the caller's working precision. Thirty further decimal
// places are retained for fixed-point rounding, exactly as in correct().
export function getGuardedEvaluator(precision,{compactFourier=false}={}){
 const p=Math.max(80,Math.min(970,Math.ceil((precision+30)/10)*10));
 const key=p+'|'+Boolean(compactFourier);
 if(!evaluators.has(key)){
  const E=createEvaluator(p,{axisShortcuts:true,cacheMatrixConstants:true,compactFourier,cacheFourierConstants:compactFourier});
  // Only ordering and zero tests use this metric. The raw squared integers
  // have no rounding or underflow; no square root is needed for LU pivots.
  const squared=z=>{const a=E.c(z);return a[0].n*a[0].n+a[1].n*a[1].n;};
  const norm2=z=>{const n=squared(z);return{n,gt(other){return n>other.n;},isZero(){return n===0n;}};};
  const rowScales=M=>M.map(row=>{
   let index=0,best=squared(row[0]);
   for(let j=1;j<row.length;j++){const n=squared(row[j]);if(n>best){best=n;index=j;}}
   return E.abs(row[index]);
  });
  evaluators.set(key,{...E,norm2,rowScales});
 }
 return evaluators.get(key);
}

export function fourierDiagnosticsSince(E,before){
 const after=E.getFourierDiagnostics();
 for(const name of ['recurrenceCalls','directCalls','guardFallbackCalls','coefficientCacheHits','coefficientCacheMisses'])after[name]-=before[name];
 return after;
}

// An under-resolved fixed-point residual cannot pass as an exact root:
// double-corrector rejects zero-at-working-precision and singular derivatives.
// Failures retry the original floating-precision Decimal engine, inside the
// same local neighbourhood and deadline. No accuracy tolerance is relaxed.
export function correct(options={}){
 if(options.config?.seed!=='tetra'||options.config.connector)return decimalCorrect(options);
 const begun=now(),precision=Math.max(Number(options.precision)||80,(Number(options.digits)||25)+30);
 let fixed;
 try{
  const arithmetic=getGuardedEvaluator(precision,{compactFourier:options.compactFourier===true}),before=arithmetic.getFourierDiagnostics();
  fixed=decimalCorrect({...options,arithmetic});
  fixed.diagnostics.arithmetic='guarded fixed-point Fourier jets';
  fixed.diagnostics.fourier=fourierDiagnosticsSince(arithmetic,before);
  if(fixed.converged)return fixed;
 }catch(error){fixed={message:String(error?.message??error),evaluations:0};}
 if(now()>(options.deadline??Infinity)||fixed.failureKind==='time_budget')return fixed;
 const fallback=decimalCorrect(options);
 fallback.elapsedMs=now()-begun;
 fallback.evaluations+=(fixed.evaluations??0);
 fallback.diagnostics.arithmetic='Decimal fallback after fixed-point check';
 fallback.diagnostics.fixedPointFallback={message:fixed.message,precision:fixed.precision,evaluations:fixed.evaluations??0};
 return fallback;
}

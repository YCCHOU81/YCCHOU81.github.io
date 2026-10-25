// Compact-source double-edge transport only. The existing bounded Decimal
// corrector supplies the equation, trust neighbourhood, stopping tests and
// fallback. A bounded same-source binary64 prepolish supplies an initial
// derivative; the unchanged fixed correction target still has to pass.
// Displayed endpoints retain their independent Decimal correction.
import {correct as decimalCorrect} from './double-corrector.mjs?v=20261024-tv-audit2';
import {getGuardedEvaluator,fourierDiagnosticsSince} from './tetra-fixed-corrector.mjs?v=20261024-tv-audit2';

const now=()=>performance.now();

export function correct(options={}){
 if(options.config?.seed!=='double')return decimalCorrect(options);
 const begun=now(),precision=Math.max(Number(options.precision)||80,(Number(options.digits)||25)+30);
 let fixed;
 try{
  const arithmetic=getGuardedEvaluator(precision,{compactFourier:true}),before=arithmetic.getFourierDiagnostics();
  // double-corrector checks the frozen source frame, weights and characters
  // before derivative reuse, then forces a fixed update before convergence.
  // An unsafe binary64 prepolish is caught inside that original corrector and
  // the ordinary fixed secant proceeds at the caller's precision and digits.
  fixed=decimalCorrect({...options,fast:true,reuseFastDerivative:true,arithmetic});
  fixed.diagnostics={...fixed.diagnostics,arithmetic:'guarded fixed-point Fourier jets with bounded binary64 prepolish; compact double-edge transport',fourier:fourierDiagnosticsSince(arithmetic,before),internalOnly:true,requiresDisplayedPrecisionCheck:true};
  if(fixed.converged)return fixed;
 }catch(error){fixed={converged:false,message:String(error?.message??error),evaluations:0};}
 if(now()>(options.deadline??Infinity)||fixed.failureKind==='time_budget')return fixed;
 const fallback=decimalCorrect(options);
 fallback.elapsedMs=now()-begun;
 fallback.evaluations+=(fixed.evaluations??0);
 fallback.diagnostics={...fallback.diagnostics,arithmetic:'Decimal fallback after compact double-edge fixed-point check',internalOnly:true,requiresDisplayedPrecisionCheck:true,
  fixedPointFallback:{message:fixed.message,precision:fixed.precision,evaluations:fixed.evaluations??0}};
 return fallback;
}

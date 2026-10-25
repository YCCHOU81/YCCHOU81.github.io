// Independent compact endpoint Decimal arithmetic, with the same precision
// and formulas as double-corrector's default evaluator. Only repeated exact
// matrix constants are reused; no fixed-point or transport data enter here.
import {createEvaluator} from './decimal-generic.mjs?v=20261024-tv-audit2';
const engines=new Map();
export function getEvaluator(precision=80){
 const p=Math.max(50,Math.min(970,Math.ceil(Number(precision)/10)*10));
 if(!engines.has(p))engines.set(p,createEvaluator(p,{cacheMatrixConstants:true}));
 return engines.get(p);
}

// Tetrahedral checkpoint adapter. An inconclusive geometric computation is a
// retryable audit, never a branch disagreement. Reuse only its checked prefix;
// solver.mjs independently requires accepted start/gate metadata before tails.
import {traceExtension as originalExtension,extendTraceGrid} from './trace-checkpoints.mjs?v=20261024-audit1';
export {extendTraceGrid};
export function traceExtension(request,family){
 const previous=request.extendFrom?.request;
 if(previous&&(!!previous.strictTracking!==!!request.strictTracking||Math.max(.01,Math.min(1,Number(previous.maxDepthStep)||1))!==Math.max(.01,Math.min(1,Number(request.maxDepthStep)||1))))throw Error('The continuation safety policy changed. Compute a separate checked trace.');
 const saved=request.extendFrom?.result;
 if(saved?.stop?.kind!=='geometry_inconclusive')return originalExtension(request,family);
 if(!saved.stop.retryable||saved.stop.branchRejected!==false)throw Error('The saved geometric audit is not a retryable inconclusive computation.');
 const adapted={...request,extendFrom:{...request.extendFrom,result:{...saved,stop:{...saved.stop,kind:'page_change'}}}};
 return originalExtension(adapted,family);
}

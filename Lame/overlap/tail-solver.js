/* Overlap adapter for the shared cusp evaluator. solver.js loads Decimal and
 * atlasAttachment first. This keeps the existing optimized precision policy. */
"use strict";
importScripts("../workers/approach-engine.js?v=20261025-maintenance1");
globalThis.previewApproach=createApproachMath({...atlasAttachment,MAX_N:21},
  {reuseAcceptedEvaluation:true,exposePrecisionEngine:true});

# Lamé numerical workers

The website remains a static site with no build step. Explore, Limits & Zoom,
and Overlap load these shared **classic scripts** using `importScripts`.
Keep the `workers` folder inside the Lamé folder when publishing.

## Entry points

| Entry | Responsibility |
| --- | --- |
| `../numerics-worker.js` | Loads the interior engine and its `solve` / `compute` message protocol. |
| `../limits-worker.js` | Loads the approach engine; owns request validation, checkpoints, and checked depth extension. |
| `../overlap/normalized-worker.js` | Schedules an Overlap case, its compact interior, and all six cusp approaches. |
| `../overlap/solver.js` | Small adapter that selects Overlap's existing interior precision policy. |
| `../overlap/tail-solver.js` | Small adapter that selects Overlap's existing cusp precision policy. |

## Shared implementation

| File | What belongs here |
| --- | --- |
| `atlas-math.js` | Double-precision elliptic evaluation, marked-domain checks, and source-mesh geometry. |
| `precision-engine.js` | Stable Fourier/Decimal evaluation and correction. Each caller gets its own precision cache. |
| `attachment-engine.js` | Interior geometric branch continuation and mesh sampling. |
| `trace-initializer.js` | The distinct initial-root correction policy used by Limits. |
| `approach-engine.js` | Marked log-epsilon continuation and independent two-precision validation. |
| `interior-messages.js` | Compact-source validation, point cache, and the interior worker message handler. |
| `attachment-anchors.js` | Generated geometric anchors for the supported orders. Do not edit individual numbers. |
| `decimal.js` | Vendored Decimal 10.6.0, with its MIT licence. |

The generated anchor table and Decimal library each have **one copy**. Do not
paste them into a worker entry or an Overlap adapter.

`reuseAcceptedEvaluation` preserves an optimization already used by Overlap:
retain a successful Newton trial's value and derivative and skip an unused final
power square. It remains off in Explore and Limits, preserving those pages'
existing arithmetic paths. `trace-initializer.js` also deliberately remains
separate: its correction policy differs from the interactive interior engine.

The refactor does not change equations, branch anchors, continuation thresholds,
two-precision checks, message payloads, or result-cache versions. Each worker
still owns its original private solver state. A sample is usable only after the
existing numerical checks succeed; a missing sample remains missing.

Worker imports resolve relative to the **entry worker URL**. The small entry
files therefore list their root or `../workers/` paths explicitly. Serve the
site over HTTP (including GitHub Pages); opening HTML through `file://` can
prevent workers from loading.

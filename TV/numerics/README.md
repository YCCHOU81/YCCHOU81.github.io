# TV numerical source map

The site is static and has no build step. `../repeated-worker.mjs` is the public
module-worker entry for all three seed families. It accepts `point`, `mesh`, and
`trace` requests and preserves the existing message IDs and result schemas.
Each worker owns its own solver and continuation state.

## Where an algorithm belongs

| Files | Responsibility |
| --- | --- |
| `double-solver.mjs`, `double-geometry.mjs`, `double-corrector.mjs` | Double-edge marked initialization, continuation and checked endpoints. |
| `double-compact-*.mjs`, `double-marked-hodge.mjs` | Compact-image correction and marked geometric checks. |
| `tetra-solver.mjs`, `tetra-geometry.mjs`, `tetra-canonical.mjs` | Tetrahedral branch, geometric initialization and labelled graph validation. |
| `tetra-compact-fast.mjs`, `tetra-fixed-corrector.mjs` | Existing bounded compact-image correction paths. |
| `tetra-lame-reduction.mjs`, `tetra-lame-anchors.mjs` | Exact opposite-pair reduction and its actual marked classical anchors. |
| `bouquet-solver.mjs`, `bouquet-anchor.mjs`, `bouquet-hodge.mjs` | Three-loop entry points and checked geometric initialization. |
| `bouquet-continuation.mjs`, `bouquet-mesh-corrector.mjs`, `bouquet-trace.mjs` | Three-loop point/mesh/trace continuation and endpoint checks. |
| `bouquet-fixed.mjs`, `bouquet-compact.mjs`, `bouquet-verified-anchors.mjs` | Existing stable matching and checked anchor support. |
| `decimal-generic.mjs`, `generic-fast.mjs`, Fourier/BigInt helpers | Shared analytic evaluation at the required working precision. |
| `marked-hodge.mjs` and family Hodge helpers | Independent geometric branch selection. |
| `compact-mesh-pool.mjs` | Nested-worker scheduling, one shared deadline and serial fallback. |
| `compact-mesh-reuse.mjs`, `base-mesh-cache.mjs` | Validation and reuse of checked compact vertices/base artifacts. |
| `trace-checkpoints.mjs`, `tetra-checkpoints.mjs` | Exact graph/path/sampling identity and checked continuation checkpoints. |
| `verified-base-double.mjs`, `verified-base-bouquet.mjs` | Generated, independently checked base grids; loaded only where needed. |

Keep the whole `numerics/` folder and `../decimal.mjs` when publishing. The
tetrahedral Lamé reduction is part of TV and is retained here; it does **not**
depend on the separate Lamé website or its classic workers.

## Changes that require care

- Keep ordered graph labels and decoration counts in cache/checkpoint keys.
  Equal weight vectors do not identify the same decorated branch.
- Retained vertices require their existing graph, grid, source, precision and
  coordinate metadata. A stored root is not a new branch initializer.
- Keep the independent Hodge, bounded correction, repeated-precision and
  endpoint agreement checks. A time stop or failed check leaves partial data.
- The generated anchor/base-grid numbers are data, not hand-edited code.
- `../cusp-formulas.js` and `../cusp-predictions.js` annotate simulation output.
  They are not numerical dependencies and must not initialize solved roots.
- Presentation/packaging changes retain numerical revisions and cache keys.
  Change an appropriate numerical revision only when old results become
  mathematically incompatible.

The old `numerics-worker.js` and `limits-worker.js` were unreachable classical
Lamé page workers left from the combined atlas. They are removed. The TV
entry, its three solver families, the tetrahedral reduction, all numerical
checks, and all generated TV datasets remain in this folder.

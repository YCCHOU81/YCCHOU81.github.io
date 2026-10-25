# Lamé Atlas

Classical one-cone Lamé only. Home: Lame.html. Pages: explore.html, limits.html, mathematics.html, overlap/index.html. TV seed selectors and TV-only datasets are excluded.

This is a self-contained static website. Keep this folder intact. There is no build step and no server-side numerical service. The workers compute in the visitor's browser.

## Open and publish

For GitHub Pages, upload the **contents of this folder** to its own repository (or chosen Pages directory). Enable Pages from that branch and directory. `index.html` opens `Lame.html` automatically. Entry filenames are case-sensitive on hosting.

For local work, start an ordinary static HTTP server in this folder. For example, `python -m http.server 8000 --bind 127.0.0.1`, then visit `http://127.0.0.1:8000/`. Opening the HTML directly as a file can prevent workers and module imports from loading.

## Files to edit

| Change | File |
| --- | --- |
| Home, page text, mathematical background | `Lame.html`, `explore.html`, `limits.html`, `mathematics.html` |
| Appearance and compact controls | `styles.css`, `repeated-controls.css`, `shell.css` |
| Navigation, sticky steps, collapsed progress, copy-link button | `site-ui.js` |
| Incoming URL cleanup and Lamé-only identity | `atlas-family.js` |
| Attachment controls and linked dual-graph highlighting | `seed-cases.js` |
| Explore interactions and computation scheduling | `main.js` |
| Limits source editor, trace scheduling, lenses, tables, export | `limits.js` |
| Common marked-domain and monodromy drawings | `atlas-view.js` |
| Recommended preview settings | `preview-settings.js` |
| Main bounded persistence | `result-cache.js` |
| Explore/Limits result identities, sharing, and depth reuse | `shared-numerics.js` |
| Numerical algorithms and generated anchors | `workers/` — see its `README.md` |
| Overlap interactions and scheduling | `overlap/preview.js` |
| Overlap bounded persistence | `overlap/result-cache.js` |
| Overlap coordinates, symmetry, contours, overlap shading | Named helper files in `overlap/` |

The HTML pages have no generated build output. Edit the corresponding source
file, refresh the browser, and keep the folder structure intact. The Lamé
controllers no longer contain TV family switching or TV decoration controls.
The TV website is independent; changing this folder does not update it.

## Numerical code: one shared copy

The worker entry points are intentionally small. `numerics-worker.js` handles
Explore points and interiors; `limits-worker.js` handles cusp traces and their
checkpoints. Overlap has adapters for its existing precision optimizations.
All load shared classic scripts from `workers/`. The large generated anchor
table and the vendored Decimal library each appear once. Do not paste a shared
engine back into a page or another worker. Do not manually edit anchor numbers.

Each worker still owns its private numerical state. The refactor preserves
the equations, branch choices, precision checks, message formats, and cache
identities. Preserve cache namespaces and numerical revisions when changing
presentation only. Change a numerical revision when a mathematical change
makes previously computed results incompatible.

## Check a maintenance change

No tools are needed for publishing. If Node.js is installed, these optional
checks use only its built-in modules; there is no `npm install`:

```sh
node tools/check.mjs
node tools/check.mjs --numerics
```

The first checks JavaScript syntax, local assets/imports, duplicate static IDs,
and collapsed help. The second also loads the real workers, checks a known
interior value and invalid-input rejection, and loads the Overlap worker.
It is a smoke check, not an exhaustive accuracy test or a browser test.

After an interaction change, open Explore, Limits & Zoom, and Overlap over
HTTP. Check dragging both source points, graph/control highlighting and
click-away clearing, automatic previews, stop/retry, deeper trace reuse,
cross-page saved results, zoom, export, and independent tint/overlap sliders.
Help and cusp progress should remain collapsed until opened. Keep visual
interpolation separate from actual computed points and retain unresolved gaps.

When deploying changed scripts, update their URL version in the HTML and any
changed worker/import entry URLs. This refreshes code without discarding
compatible saved numerical results.

There are no required assets outside this folder, no npm dependencies, and no connection to the old local website. Browser-result caches have a dedicated atlas name, so these copies do not change the original site's saved results. Explore and Limits still share this atlas's checked numerical data. The saved result budget is approximately 16 MB per atlas, with older records evicted; this is not a bound on total browser memory. The Overlap Atlas has its own additional 8 MB result cache.

Computed samples and interpolated contours remain numerical evidence. Missing or unchecked tail portions stay open; display settings do not create new numerical roots.

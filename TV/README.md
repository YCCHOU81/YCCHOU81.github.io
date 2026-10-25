# TV Atlas

A self-contained static website for the **double-edge, tetrahedral, and three-loop seed families**, including their admissible decorations. Home: `TV.html`. Other pages: Explore, Limits & Zoom, Graphs, and Mathematics. Classical Lamé has its own independent website.

There is no build step, npm dependency, numerical service, or required file outside this folder. Numerical work runs in the visitor's browser. Keep the entire folder intact when publishing, including `numerics/` and `graphs/`.

## Open and publish

For GitHub Pages, upload this folder's contents to the repository or chosen Pages directory, then enable Pages for that location. `index.html` redirects to `TV.html`. Hosted filenames are case-sensitive.

For local work, serve this folder over HTTP:

```sh
python -m http.server 8000 --bind 127.0.0.1
```

Open `http://127.0.0.1:8000/`. Opening HTML through `file://` can block module workers and their imports.

## File map

| Change | File |
| --- | --- |
| Home, page text, mathematics | `TV.html`, `explore.html`, `limits.html`, `mathematics.html` |
| Appearance and compact controls | `styles.css`, `repeated-controls.css`, `shell.css` |
| Navigation, sticky steps, collapsed progress, copy-link button | `site-ui.js` |
| Allowed families/operations and incoming URL cleanup | `atlas-family.js` |
| Family/decoration controls, remembered choices, linked dual graphs | `seed-cases.js` |
| Explore mounting, family changes, display persistence | `main.js` |
| Explore point/image/trace scheduling and source interactions | `seed-explore.js` |
| Limits source editor, traces, lenses, playback, evidence, export | `limits.js` |
| Common monodromy and marked-domain drawings | `atlas-view.js` |
| Double-edge and tetrahedral graph/count model | `decoration-model.js` |
| Three-loop registry, model, graph layout | `bouquet-data.js`, `bouquet-model.js`, `bouquet-layout.js` |
| Recommended sampling/precision/time budgets | `preview-settings.js` |
| Bounded browser persistence | `result-cache.js` |
| Explore/Limits sharing, checked partial images, trace extension | `shared-numerics.js` |
| Numerical request routing | `repeated-worker.mjs` |
| Solvers, pools, checkpoints, generated anchors/base roots | `numerics/` — see its `README.md` |
| Conjectural formula annotations | `cusp-formulas.js`, `cusp-predictions.js`, `cusp-formula-checks.json` |
| Hemisphere/dual-graph lab | `graphs/` — see its `README.md` |

Edit the file responsible for the change, then refresh the browser. The TV controllers support all three families; they no longer carry unreachable classical Lamé page controllers or their large unused worker bundles. The tetrahedral solver's actual Lamé reduction is retained in `numerics/`.

## Preserve the important boundaries

- The labelled decorated graph, not just its four weights, identifies a numerical branch. Keep the full case configuration in requests and keys.
- Conjectural cusp formulas annotate solved data. They are not imported by numerical solvers and must not replace roots or create ideal endpoints.
- Numerical changes belong in the corresponding solver module. Preserve repeated-precision, healthy Hodge-mesh, marking, trust-radius, and branch checks when changing scheduling or presentation.
- The three families keep separate numerical state while sharing the bounded store. Keep family switches in place and dispose the old controller's workers before restoring another family.
- Compatible checked partial meshes and trace prefixes retain their actual source/grid/path metadata. A different graph or unchecked branch cannot initialize a new computation from a saved result.
- Generated anchors and verified base meshes are data, not hand-edited code. Keep their precision/provenance metadata and import paths.

Detailed existing solver behavior, cusp references, count limits, slower-device fallbacks and accuracy caveats are retained in [NUMERICAL-NOTES.md](NUMERICAL-NOTES.md).

## Caches and saved settings

The result store is approximately **16 MB / 32 entries**, with older records evicted. IndexedDB retains results when available; a memory fallback keeps the site usable without it. This budget does not bound total browser RAM. The store has a TV-specific name and is independent of the Lamé site.

Explore and Limits share checked results. Family snapshots, display preferences, base meshes, and deeper trace continuation retain their existing identities. Transient drag previews must not overwrite checked released-point results.

Keep numerical revisions and cache names unchanged for presentation-only work. Change a numerical revision when a mathematical change invalidates old roots. Use script/worker URL versions to refresh deployed code without discarding compatible numerical results.

## Check a maintenance change

Publishing requires no development tools. If Node.js is installed, these optional checks need only its built-in modules:

```sh
node tools/check.mjs
node tools/check.mjs --numerics
```

The first checks JavaScript syntax, local imports/assets, duplicate static IDs, collapsed help, and graph topology/face labels for all decoration kinds. The second also loads the real TV worker and computes one simple checked point per seed family. It is a smoke check, not an exhaustive accuracy or browser suite.

After controller changes, check all three families in Explore and Limits: change decorations/counts, drag either linked point, apply the automatic preview, stop/retry, increase depth, return to a computed family, switch pages, zoom, play and export. Verify saved display settings and that help and cusp progress stay collapsed until opened. Test the Graphs page's reciprocal edge/arc highlights and the independent Mathematics page as well.

The expanded Graphs lab intentionally covers only small base-core pair decorations: 0–3 pairs per chain and at most 16 hemispheres. Unsupported incoming configurations are labelled explicitly. Explore retains its full decorations and larger numerical capacity limits.

Missing or unchecked numerical portions remain open in the display.

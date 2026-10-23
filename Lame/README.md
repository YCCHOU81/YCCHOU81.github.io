# Lamé Atlas

Classical one-cone Lamé only. Home: Lame.html. Pages: explore.html, limits.html, mathematics.html, overlap/index.html. TV seed selectors and TV-only datasets are excluded.

This is a self-contained static website. Keep this folder intact. There is no build step and no server-side numerical service. The workers compute in the visitor's browser.

## Open and publish

For GitHub Pages, upload the **contents of this folder** to its own repository (or chosen Pages directory). Enable Pages from that branch and directory. `index.html` opens `Lame.html` automatically. Entry filenames are case-sensitive on hosting.

For local work, start an ordinary static HTTP server in this folder. For example, `python -m http.server 8000 --bind 127.0.0.1`, then visit `http://127.0.0.1:8000/`. Opening the HTML directly as a file can prevent workers and module imports from loading.

## Files to edit

- `Lame.html`: home content.
- `explore.html`, `limits.html`, `mathematics.html`: the corresponding page content.
- `styles.css`, `repeated-controls.css`, `shell.css`: appearance. The shared navigation is the same across this folder's pages.
- `atlas-family.js`: the allowed families, entry page, and incoming URL normalization.
- `seed-cases.js`: graph and decoration controls; this copy restricts initialization, menus, and family switching to this atlas.
- `main.js`, `limits.js`, numerical workers, and supporting modules: existing computation and drawing code. Some shared internal routines remain in both folders; only this atlas's controls are exposed.

There are no required assets outside this folder, no npm dependencies, and no connection to the old local website. Browser-result caches have a dedicated atlas name, so these copies do not change the original site's saved results. Explore and Limits still share this atlas's checked numerical data. The saved result budget is approximately 16 MB per atlas, with older records evicted; this is not a bound on total browser memory. The Overlap Atlas has its own additional 8 MB result cache.

Computed samples and interpolated contours remain numerical evidence. Missing or unchecked tail portions stay open; display settings do not create new numerical roots.

# TV graph lab

Open `index.html` from the TV Atlas navigation. Upload this whole `graphs/`
folder with its parent TV folder; it uses the parent's `atlas-family.js`,
`shell.css`, `decoration-model.js`, `bouquet-data.js`, and `bouquet-model.js`.
There are no external libraries, fonts, services, or localhost dependencies.

The hemisphere graph and dual graph are generated from the same exact ribbon
map. Selecting an edge, a hemisphere/dual vertex, or a marked point/dual face
highlights its matching object. Edge controls and red arc sums are linked in
both directions. Click elsewhere or use Clear to remove a highlight.

The four steps distinguish topology, actual spherical lengths, developed
half-turn axes, and the ordinary-character shift. All displayed equations use
native MathML. The amber insertion diagram is an unapplied operation preview;
the red triangle is the result for the currently applied decoration.

For readable expanded drawings, this graph lab covers the three **base cores**
with **0–3 pairs per chain**, and at most **16 downstairs hemispheres**. It does
not draw the other three-loop cores, sprouts, center splits, or large repeated
decorations supported by Explore. This restriction is shown on the page.

Links from Explore or Limits & Zoom restore the seed, supported pair counts
(including the double-edge and tetrahedral one-pair presets), and source point.
If the requested decoration exceeds this scope, a visible
notice says the undecorated base core is displayed instead, with a link back
to the exact requested full case. The separate “Open this marked case in
Explore” link always carries the graph actually displayed, including its
current counts and source point. It does not launch a period solve automatically.

This page computes graph holonomy and monodromy shifts, not elliptic periods
or conformal images. Numerical graph-product checks are diagnostics rather
than a replacement for a proof of general asymptotic formulas.

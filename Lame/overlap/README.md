# Classical Lamé Overlap Atlas

The Overlap Atlas compares every classical Lamé attachment of a chosen order, including its computed cusp tails. It is part of the Lamé-only website: Home, Explore, Limits & Zoom, Overlap Atlas, and Mathematics share the same navigation. All numerical scripts required by this page are included in this folder.

## Open locally

Serve the parent Lamé folder over HTTP. With Python installed, run from that folder:

```sh
python -m http.server 8775 --bind 127.0.0.1
```

Open **http://127.0.0.1:8775/Lame.html**, then choose **Overlap Atlas**. The page itself is at **http://127.0.0.1:8775/overlap/**. Use HTTP rather than opening `index.html` directly, so the numerical workers can load their scripts. The same folder can be deployed as a static website, including GitHub Pages.

## Start with the overall picture

1. Choose an integer **n = 1,…,21**. The atlas includes all **ordered** nonnegative integer triples with **a + b + c = n − 1**: there are **n(n + 1)/2** labeled images.
2. The atlas computes the interior and **six cusp approaches for each image**: two source paths toward each of the three cusps. Interiors and resolved tail samples appear progressively.
3. The left **Γ(2) overview always shows the original τ-plane**, with equal units on its real and imaginary axes. The imaginary part is not compressed. An infinite cusp necessarily extends beyond every finite window.
4. The right plot defaults to a **common disk coordinate**. It places the three cusps on a finite circle while retaining the relative locations and actual overlaps of all images.

Drag either large plot to pan; use the mouse wheel or **+/−** buttons to zoom. **Fit data** and **Fit images & tails** include the available computed tail points. **Reset** restores the standard Γ(2) window on the left. A common fit may make some narrow tails look very small: zoom into them rather than changing the underlying samples.

## Compare shapes, spreading, and overlaps

Choose **All images overlaid**, **Selected image**, or **Separate images**. Gallery panels use the same coordinate and a shared scale. Select an attachment triple to emphasize its image; **Shift-click** its label to hide or show that case. **Triangle tint** and **Overlap darkness** adjust the single-image fill and extra overlap shading independently.

Interior, cusp-tail, and tail-width layers can be shown independently. Smooth dashed contours pass through the computed points; sample markers are hidden by default and can be enabled with **Sample points**. Four additional transition samples connect each cusp side to its actual compact-boundary point. The intermediate curve is a bounded display interpolation, not additional solved roots or a certified boundary. Missing or unresolved samples remain gaps, and no finite endpoint is substituted for infinity.

Individual images have a light **10% Triangle tint** by default. **Overlap darkness** is a separate slider, initially 55%: it adds shading only where at least two distinct attachment images cover the same pixels, increasing with the number of images. Changing it leaves singly covered regions unchanged. The base tint is filled once across the visible union, so it does not accumulate as cases overlap. At zero overlap darkness, all filled regions share the base tint. Both controls redraw saved points without recomputation.

Each attachment's interior and tail ribbons are counted together once, so a connection seam does not count as an additional overlap. Dashed contours are drawn above the shading. The ribbon's inner edge follows the computed compact boundary rather than a straight join chord. Tail width remains an interpolation between the two source approaches, not a certified enclosure of the full domain.

The pointer readout counts overlaps of the sampled interior triangles and interpolated tail ribbons. It is numerical evidence, not a certified count of overlaps of the complete conformal domains.

## Eleven coordinate choices

The measured τ-values are retained unchanged. Switching coordinates, layouts, visibility, or opacity redraws them without recomputation.

### Common coordinates: inspect actual overlaps

- **Original τ:** no change of variable. Both axes remain linear.
- **Common disk** (default): **ζ = (τ − τ₀)/(τ − conjugate(τ₀))**, with the same fixed **τ₀ = 1/2 + i√3/2** for every attachment. This is one injective holomorphic coordinate on the upper half-plane, so it preserves all actual overlaps. The cusp limits ∞, 0, and 1 lie respectively at **1**, **−1/2 + i√3/2**, and **−1/2 − i√3/2** on its boundary. The drawn tails stop at their computed depths.

### Recenter each branch: compare shapes

For each attachment the reference is **τ* = Fₐᵦ꜀(1/3, 1/3)**.

- **Centered and scaled by Im τ*:** **Z = (τ − τ*)/Im τ***. This is a translation and one positive real scaling applied equally to both axes, including every tail.
- **Centered and equal interior scale:** **Z = (τ − τ*)/R**, where **R** is the sampled root-mean-square distance of the interior points from τ*. Tail points do not determine this scale.
- **Centered disk:** **ζ = (τ − τ*)/(τ − conjugate(τ*))**. Each image has its own disk center.

These three modes use a different coordinate for each branch. They are useful for comparing shapes and scales, but they change the relative positions and overlaps between different images. They do not assert that the full images are equal-sized or congruent.

### Common cusp charts: open one cusp

- **Cusp 0:** **w = −1/τ**.
- **Cusp 1:** **w = −1/(τ − 1)**.
- **Cusp ∞:** **w = τ**.

Each sends the named cusp to infinity, using the same injective coordinate for all images, so actual overlaps are preserved. The real and imaginary axes of each chart remain linear.

### Local cusp coordinates: inspect the end near zero

- **q₀ = exp(πi(−1/τ))**.
- **q₁ = exp(πi/(1 − τ))**.
- **q∞ = exp(πiτ)**.

These send their respective cusp to zero. The exponential identifies points differing by **w → w + 2**. Consequently, overlaps in these views are faithful only within a chosen strip of width two; different translates can otherwise appear superimposed. Extremely small q-values can also underflow at floating-point display precision.

## Sampling, accuracy, and continuation

All attachments use the same interior mesh, whose source barycentric coordinates are at least **0.06**. That compact interior is supplemented by the six actual cusp approaches; the tails are not discarded.

The recommended tail depth is **ε = 10⁻⁵ for n ≤ 8**, and **ε = 10⁻³ for n > 8**, to keep the first comparison responsive as the number and complexity of the cases increase. A user can instead choose **10⁻³**, **10⁻⁵**, **10⁻⁸**, **10⁻¹⁰**, **10⁻¹⁵**, or **10⁻²⁰**, with Quick, Standard, or Detailed tail sampling. The two deepest choices wait for the Compute button and can take substantially longer for many cases. The chosen ε specifies the source approach parameter, not an error estimate for τ.

Tail nodes use marked continuation in log ε and a stable Fourier evaluator. Plotted nodes are refined and compared at **45- and 70-digit working precision**, increasing to **70 and 100 digits** for depths beyond 10; failed comparisons remain gaps. This checks numerical consistency but is **not an interval certificate**, a proof of a complete cusp expansion, or a guarantee of that many correct digits in the displayed τ-value.

Increasing depth reuses previously resolved nodes and continues the same approaches. Existing interior samples are reused too; denser settings add missing samples. **Stop** preserves resolved data. **Compute or extend** resumes the case and retries unresolved portions. Downloaded JSON contains the original τ-values and their approach parameters, rather than just transformed screen coordinates.

## Local computation and bounded storage

Computation runs on the user's device in up to **three background workers**, according to its capabilities. Exact triangle relabelings reuse solves between symmetry-related ordered triples. Larger orders have more independent cases; deeper or denser tail sampling adds work.

The Overlap Atlas uses its **own result cache**, limited to approximately **8 MB** of serialized records and at most 600 records, in the IndexedDB database `lame-atlas-overlap-v1`. IndexedDB retains results when available; otherwise they last for the session. Older entries are evicted as needed. This bound covers the saved numerical-result cache, not the browser's total working memory. Its cache is separate from the single-case Explore and Limits & Zoom cache, and no computation service or upload is required.

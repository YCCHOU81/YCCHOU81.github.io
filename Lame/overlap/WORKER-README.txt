Classical Lame Overlap Atlas: numerical adapter
===============================================

solver.js is byte-for-byte the current atlas numerics-worker.js. Its checked
branch anchors and analytic equations have not been modified. The adapter is
normalized-worker.js, a classic Worker importing solver.js locally. New isolated
tail helpers reuse the atlas's stable Decimal/Fourier evaluator. The original
solver.js and the current website remain unchanged.

Request (one labelled case at a time):
  {type:"run", id, counts:[a,b,c], resolution:4, inset:0.06,
   timeBudgetMs:45000, depth:5, density:1, tailBudgetMs:90000,
   existingRecord:null}

Counts must be nonnegative integers and n=1+a+b+c must be at most 21.
Resolution is rounded and restricted to 2..24. All cases at a fixed resolution
and inset use exactly the same mesh, including source indexing and triangles.
The default compact source has all three barycentric weights at least 0.06.
The common reference point is solved at (r,s)=(1/3,1/3), separately from the
mesh; it is never approximated by its nearest vertex or another sheet.

Messages, in order:
  start     {id,counts,points,triangles,resolution,inset,total}
  reference {id,counts,reference,status,elapsedMs}
  sample    {id,counts,index,source:{r,s},tau,status,completed,total,elapsedMs}
  progress  {id,completed,total,resolved}
  tail-start  {id,counts,tailIndex,edge,t,total,depth,density}
  tail-sample {id,counts,tailIndex,edge,t,node,completed,total,elapsedMs}
  tail-done   {id,counts,tailIndex,edge,t,tail,elapsedMs}
  done      {id,record}
  error     {id,message,elapsedMs}  -- invalid job or unexpected adapter failure

The final record includes counts,n,points,triangles,tau,reference,
referenceSource,resolution,inset,status,resolved,total,completed,stopReason,
elapsedMs,diagnostics,pointDiagnostics. Each tau is {re,im} or null. The client
must not fill any triangle touching null data. A failed reference is null and
does not define a normalized image. No fabricated or replacement roots occur.

Cusp tails are mandatory parts of each run, following the interior image. Each
record includes tails, tailDepth, tailDensity, tailStatus, tailResolved,
tailTotal, tailElapsedMs and transitionVersion:"join-v1". Interior status
remains separate from tailStatus.
There are six trajectories, two edge positions t=0.2 and t=0.8 for each of the
three open source edges. The edge convention uses barycentric coordinates
b=(1-2s, 2(r+s)-1, 1-2r): edge0 approaches infinity, edge1 approaches1,
and edge2 approaches0. Each requested epsilon is strictly positive.

Tail record schema:
  {edge,t,position,path,cusp,nodes,status,resolved,total,requestedEpsilons,
   reused,stopReason,elapsedMs,continuationSteps,method,
   transitionInset,transitionVersion:"join-v1",
   numericalNotIntervalCertified:true}
Tail node schema:
  {epsilon,r,s,source:{r,s},rExact,sExact,barycentric,
   tau:{re,im}|null,status,solverStatus,tauExact,errorEstimate,
   comparisonTauError,precisionDigits:45,validationDigits:70,
   elapsedMs,continuationSteps}
With the default inset0.06, the join-v1 epsilon grid begins
0.06,0.055,0.05,0.045,0.04, then follows logarithmic steps from0.01 through
1e-depth. The first node is an actual corrected root at the compact source
boundary b[edge]=inset, shared with the displayed interior perimeter. The four
additional transition samples remove the previously uncomputed gap between
the inset0.06 interior and the0.04 tail start. If another inset is requested,
the transition grid begins at that inset. Depth is3..20, density is1..3 samples
per decade. Actual nodes only are output; continuation stopping leaves the
tail partial. There is no invented cusp endpoint, extrapolated endpoint, or
schematic substitute for a failed root.

The first tail root is initialized from the same labelled geometric anchor.
Small steps in log(epsilon), with secant prediction and a movement guard, then
carry that root. Every plotted node is corrected using the stable Fourier
evaluator at45 digits and independently checked at70 digits through depth10.
Optional depths15 and20 use70 working digits and100 validation digits;
the default remains depth3 or5. Deep jobs may use a180-second tail budget and
can take substantially longer, especially for all triples at a large order.
Direct large
point jumps are deliberately avoided: at high order they can converge to a
different root even when the residual is small. The older70/100digit full
approach routine is retained as an initialization fallback and accuracy baseline.

When a count permutation fixes the branch, exact geometric relabeling reuses
the corresponding computed tail. These stabilizer images are marked
symmetryDerived/inheritedRepresentativeError. Their decimal coordinates and
Möbius transforms are computed at70 digits. Thus a balanced branch requires
one actual trajectory and its five exact relabelings; a branch with two equal
counts requires three actual trajectories. Every branch still returns all six.

An existingRecord with matching counts, resolution and inset reuses its complete
interior. Compatible validated tail nodes are reused, and a greater depth
continues from the retained exact endpoint. Failure gaps remain visible and
cannot seed continuation. Old cached records beginning at0.04 can be upgraded
by computing only the missing transition prefix: initialization at the compact
boundary still uses the same labelled geometric branch, and every cached
0.04-and-deeper tauExact value is retained unchanged. This cache belongs to the
separate preview.

Each mesh walk starts near the reference, then visits the nearest unvisited
source vertex. Successful complete solver states seed subsequent points.
Failures use only the same labelled branch's fixed anchor through the existing
solver's fallback. A failed point remains a hole. The 45-second default case
budget is checked between point solves; a single solve cannot be interrupted
except by terminating the Worker. Optional timeBudgetMs is limited to5..180s.
Termination is the intended cancellation mechanism; do not overlap jobs in
one Worker.

Accuracy is numerical, not interval-certified. The existing engine checks
its Newton correction and an independently evaluated equation; for n>12,
45- and70-digit arithmetic are compared. Use comparisonTauError and the
solver's explicit status, rather than the high-degree polynomial's raw
relative residual, to assess numerical conditioning. The latter has a
degree-dependent scale and is deliberately not presented by the adapter.

Reproducible checks (outside this deployable directory):
  work/lame-normalized-preview/benchmark-worker.cjs
  work/lame-normalized-preview/worker-benchmark.json
  work/lame-normalized-preview/check-permutations.cjs
  work/lame-normalized-preview/permutation-checks.json

2026-10-23 Node24 check: all120 mesh points and8 references succeeded for
corner/balanced examples at n=4,8,12,20, resolution4 (15 points), inset0.06.
Per-case elapsed times were approximately 0.35–0.43s,0.63–0.70s,
1.18–1.39s,6.45–6.62s respectively. Browser and device timings vary.
Independent cyclic/reflection covariance checks passed at24 comparisons,
maximum absolute tau discrepancy1.1e-10 (double precision tolerance2e-10).
These are arithmetic checks, not by themselves a proof of geometric sheet
correspondence. Exact relabeling is justified separately by the triangle-edge
and period-basis correspondence in symmetry.js.

2026-10-23 join-v1 transition checks:
  check-tail-transition.cjs/tail-transition-verification.json: n6(1,2,2),
  all six initial nodes lie exactly on the inset0.06 source boundary; all42
  depth3 nodes resolved. Adding the missing prefix to an old-style cached
  record retained all18 old nodes exactly and took2.14s; the fresh run took
  5.00s. Extending to depth5 resolved54 nodes, and a fully cached repeat took
  0.090s. Maximum two-precision tau discrepancy was3.22e-22.
  check-tail-worker.cjs: balancedn4 depth3 resolved42nodes in1.99s;
  depth5 extension retained all42 earlier nodes and added12 in0.86s;
  complete54node cache reuse took0.088s.

Earlier tail checks (before the join-v1 transition prefix):
  check-tail-worker.cjs checks events, independent tail status, depth extension,
  exact retention of18 earlier nodes, and full cache reuse.
  tail-worker-verification.json: balancedn4 depth3 first18nodes1.55s,
  depth5 extension12newnodes0.82s, fullrevisit0.067s.
  benchmark-worker-tails.cjs/worker-tail-speed.json: full six-tail worker jobs,
  including six interior samples, took4.76s atn4(depth5),9.20s atn8(depth5),
  9.59s atn20 with repeated counts(depth3),21.26s atn20 with three distinct
  counts(depth3). All requested nodes resolved in those runs. Timings vary by
  machine, parallel load and decoration. Large all-case comparisons take longer.
  tail-accuracy-verification.json:26 nodes acrossn4,n8,n20 agreed with the
  original70/100digit stable approach to the full binary64 display precision.

Optional-depth20 validation:
  n4(0,0,3): all six trajectories resolved120/120 nodes through1e-20.
  n20(0,0,19): both positions on edges0 and1 resolved80/80 nodes through1e-20.
  The infinity-edge n20 trajectories agreed with the original70/100digit
  approach (only its input epsilon guard extended) at all40 stored nodes.
  The n4 infinity endpoints and n20 edge1 endpoints also agreed at display
  precision. Endpoints and full details appear in work/.../deep-tail*.json.
  Extending the balancedn4 worker from depth3 to20 retained18 earlier nodes,
  resolved120 total nodes, and took11.83s; a complete repeat took0.092s.
  These checks establish numerical examples, not a uniform error theorem.

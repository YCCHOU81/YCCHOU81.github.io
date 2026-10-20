"use strict";

(() => {
  const $ = id => document.getElementById(id);

  if (!$("limits-controls")) return;

  const color = "#b55e3c";
  const baseColor = "#285f65";

  const s = {
    job: 0,
    workers: {},
    pending: new Set(),
    nodes: [],
    mesh: null,
    step: -1,
    chart: "raw",
    fit: "image",
    view: null,
    result: null,
    request: null
  };

  const meshCache = new Map();
  const traceCache = new Map();
  const query = new URLSearchParams(location.search);

  ["a", "b", "c"].forEach((name, i) => {
    $("limit-" + name).value = query.get(name) ?? [1, 0, 0][i];
  });

  if (/^(edge|vertex)[012]$/.test(query.get("path"))) {
    $("limit-path").value = query.get("path");
  }

  if (["3", "5", "7", "10"].includes(query.get("depth"))) {
    $("limit-depth").value = query.get("depth");
  }

  if (query.has("t")) {
    $("limit-position").value = query.get("t");
  }

  const good = node =>
    node?.status === "ok" && node.result?.status === "ok";

  const values = (node, branch = "result") =>
    node?.[branch]?.views?.[s.chart];

  const lastGood = () =>
    s.nodes.findLastIndex(good);

  const fmt = (x, digits = 6) =>
    Number.isFinite(Number(x))
      ? Number(Number(x).toPrecision(digits)).toString()
      : "—";

  const escape = x =>
    String(x)
      .replaceAll("&", "&amp;")
      .replaceAll("<", "&lt;")
      .replaceAll('"', "&quot;");

  /*
    Subtract decimal strings BEFORE converting the small difference
    to an ordinary JavaScript number.
  */
  function difference(a, b) {
    function parts(text) {
      const [mantissa, exponent = "0"] =
        String(text).toLowerCase().split("e");

      const [whole, fraction = ""] = mantissa.split(".");

      return {
        n: BigInt(whole + fraction),
        e: Number(exponent) - fraction.length
      };
    }

    const x = parts(a);
    const y = parts(b);
    const e = Math.min(x.e, y.e);

    return Number(
      x.n * 10n ** BigInt(x.e - e) -
      y.n * 10n ** BigInt(y.e - e)
    ) * 10 ** e;
  }

  const offset = (z, origin) => [
    difference(z[0], origin[0]),
    difference(z[1], origin[1])
  ];

  function meshValue(sample) {
    const branch = sample?.result;

    if (branch?.status !== "ok") return null;

    const {re: x, im: y} = branch.tau;

    if (s.chart === "raw") {
      return [String(x), String(y)];
    }

    const u = s.chart === "zero" ? x : 1 - x;
    const den = u * u + y * y;

    return [
      String(s.chart === "zero" ? 1 - x / den : u / den),
      String(y / den)
    ];
  }

  function pathNote() {
    const kind = $("limit-path").value.startsWith("edge")
      ? "edge"
      : "vertex";

    $("limit-path-note").textContent =
      "Source weights are (1−2s, 2(r+s)−1, 1−2r). " +
      (
        kind === "edge"
          ? "The selected weight is ε; the other two are " +
            "(1−ε)t and (1−ε)(1−t)."
          : "The selected weight is 1−ε; the other two are " +
            "εt and ε(1−t)."
      ) +
      " The other weights are taken in increasing index order. " +
      "ε is a source parameter.";
  }

  function status(kind, text, error = false) {
    const el = $("limits-" + kind + "-status");

    el.textContent = text;
    el.classList.toggle("error", error);
  }

  function settle(kind) {
    s.workers[kind]?.terminate();

    delete s.workers[kind];
    s.pending.delete(kind);

    $("limits-run").disabled = s.pending.size > 0;
    $("limits-cancel").disabled = !s.pending.size;
  }

  function cancel() {
    ++s.job;

    for (const kind of [...s.pending]) {
      settle(kind);

      status(
        kind === "mesh" ? "image" : "trace",
        "Cancelled. Computed samples are retained."
      );
    }
  }

  function startWorker(kind, file, data, handle) {
    s.pending.add(kind);

    $("limits-run").disabled = true;
    $("limits-cancel").disabled = false;

    try {
      const worker = new Worker(file);

      s.workers[kind] = worker;

      worker.onmessage = ({data: reply}) => {
        if (reply.id === s.job) handle(reply);
      };

      worker.onerror = () => {
        if (data.id !== s.job) return;

        settle(kind);

        status(
          kind === "mesh" ? "image" : "trace",
          "Could not run " + file +
          ". Check the file is beside limits.html and open " +
          "the site through GitHub Pages or a local server.",
          true
        );
      };

      worker.postMessage(data);
    } catch (error) {
      settle(kind);

      status(
        kind === "mesh" ? "image" : "trace",
        error.message,
        true
      );
    }
  }

  function finishTrace(result) {
    s.result = result;
    s.nodes = result.nodes;
    s.step = lastGood();

    if (
      s.fit !== "image" &&
      Number($("limit-zoom").value) === 0
    ) {
      s.view = null;
    }

    const completed = s.nodes.filter(good).length;

    status(
      "trace",
      (
        result.status === "ok"
          ? "Trace complete: "
          : "Trace stopped early: "
      ) +
      completed + " checked samples. " +
      (result.message ?? "") +
      (
        result.elapsedMs != null
          ? ` ${(result.elapsedMs / 1000).toFixed(1)} seconds.`
          : ""
      )
    );

    render();

    if (Number($("limit-zoom").value) > 0) {
      centerSelected();
    }
  }

  function run(event) {
    event?.preventDefault();

    if (!$("limits-controls").reportValidity()) return;

    const counts = ["a", "b", "c"].map(name =>
      $("limit-" + name).valueAsNumber
    );

    if (
      !counts.every(x =>
        Number.isSafeInteger(x) && x >= 0 && x <= 11
      ) ||
      counts.reduce((a, b) => a + b) > 11
    ) {
      status(
        "trace",
        "Enter nonnegative whole numbers with a + b + c ≤ 11.",
        true
      );
      return;
    }

    const pathValue = $("limit-path").value;

    const path = {
      kind: pathValue.startsWith("edge") ? "edge" : "vertex",
      index: Number(pathValue.at(-1)),
      position: $("limit-position").value
    };

    const target = "1e-" + $("limit-depth").value;
    const density = Number($("limit-density").value);

    cancel();

    const id = ++s.job;

    s.nodes = [];
    s.mesh = null;
    s.result = null;
    s.step = -1;
    s.view = null;

    s.request = {
      counts,
      path,
      target,
      samplesPerDecade: density
    };

    $("limit-zoom").value = 0;

    $("limits-applied").textContent =
      `Applied attachment: (${counts.join(", ")}), ` +
      `n = ${1 + counts.reduce((a, b) => a + b)}. ` +
      `Path: ${pathValue}, t = ${path.position}; ` +
      `requested ε = ${target}.`;

    const url = new URL(location.href);

    ["a", "b", "c"].forEach((name, i) => {
      url.searchParams.set(name, counts[i]);
    });

    url.searchParams.set("path", pathValue);
    url.searchParams.set("t", path.position);
    url.searchParams.set("depth", $("limit-depth").value);

    history.replaceState(null, "", url);

    document.querySelectorAll("a[data-page]").forEach(link => {
      link.href =
        link.dataset.page + "?" +
        new URLSearchParams({
          a: counts[0],
          b: counts[1],
          c: counts[2]
        });
    });

    const meshKey = counts.join(",");
    const traceKey = JSON.stringify([
      counts,
      path,
      target,
      density
    ]);

    status("image", "Computing a 45-point interior image…");

    status(
      "trace",
      "Following the geometric branch and checking each sample " +
      "at two precisions…"
    );

    render();

    if (meshCache.has(meshKey)) {
      s.mesh = meshCache.get(meshKey);

      status(
        "image",
        "Interior image reused from this tab's cache."
      );
    } else {
      startWorker(
        "mesh",
        "numerics-worker.js",
        {
          type: "compute",
          id,
          counts,
          resolution: 8,
          inset: 0.04
        },
        reply => {
          if (reply.type === "progress") {
            status(
              "image",
              `Image samples: ${reply.completed}/${reply.total}.`
            );
          }

          if (reply.type === "error") {
            settle("mesh");
            status("image", reply.message, true);
          }

          if (reply.type === "done") {
            settle("mesh");

            s.mesh = reply.mesh;
            meshCache.set(meshKey, s.mesh);

            if (meshCache.size > 4) {
              meshCache.delete(meshCache.keys().next().value);
            }

            const count = s.mesh.samples.filter(x =>
              x.result?.status === "ok"
            ).length;

            status(
              "image",
              `Interior image: ${count}/${s.mesh.points.length} ` +
              "resolved samples. Unresolved faces are omitted."
            );

            if (
              s.fit === "image" &&
              Number($("limit-zoom").value) === 0
            ) {
              s.view = null;
            }

            render();
          }
        }
      );
    }

    if (traceCache.has(traceKey)) {
      finishTrace(traceCache.get(traceKey));
    } else {
      startWorker(
        "trace",
        "limits-worker.js",
        {
          type: "trace",
          id,
          counts,
          path,
          target,
          samplesPerDecade: density
        },
        reply => {
          if (reply.type === "progress") {
            s.nodes.push(reply.node);
            s.step = lastGood();

            status(
              "trace",
              `Trace samples: ${reply.completed}/${reply.total}; ` +
              `ε = ${reply.node.epsilon}.`
            );

            if (
              (s.fit !== "image" || !s.view) &&
              Number($("limit-zoom").value) === 0
            ) {
              s.view = null;
            }

            render();

            if (Number($("limit-zoom").value) > 0) {
              centerSelected();
            }
          }

          if (reply.type === "error") {
            settle("trace");
            status("trace", reply.message, true);
          }

          if (reply.type === "done") {
            settle("trace");

            if (reply.result.status === "ok") {
              traceCache.set(traceKey, reply.result);

              if (traceCache.size > 4) {
                traceCache.delete(traceCache.keys().next().value);
              }
            }

            finishTrace(reply.result);
          }
        }
      );
    }

    render();
  }

  /* Shared canvas axes */

  function frame(canvas, bounds, xLabel, yLabel) {
    const width = Math.max(
      200,
      canvas.getBoundingClientRect().width
    );

    const height = 340;
    const ratio = Math.min(devicePixelRatio || 1, 2);

    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);

    const c = canvas.getContext("2d");

    c.setTransform(ratio, 0, 0, ratio, 0, 0);

    const box = {
      x: 68,
      y: 18,
      w: width - 88,
      h: height - 68
    };

    const map = p => [
      box.x +
        (p[0] - bounds[0]) /
        (bounds[1] - bounds[0]) * box.w,

      box.y +
        (bounds[3] - p[1]) /
        (bounds[3] - bounds[2]) * box.h
    ];

    c.font = "11px system-ui";
    c.fillStyle = "#596e62";
    c.strokeStyle = "#d6ded2";

    for (let i = 0; i <= 4; i++) {
      const x =
        bounds[0] + (bounds[1] - bounds[0]) * i / 4;

      const y =
        bounds[2] + (bounds[3] - bounds[2]) * i / 4;

      const [px] = map([x, bounds[2]]);
      const [, py] = map([bounds[0], y]);

      c.beginPath();
      c.moveTo(px, box.y);
      c.lineTo(px, box.y + box.h);
      c.moveTo(box.x, py);
      c.lineTo(box.x + box.w, py);
      c.stroke();

      c.textAlign = "center";
      c.fillText(fmt(x), px, box.y + box.h + 18);

      c.textAlign = "right";
      c.fillText(fmt(y), box.x - 7, py + 4);
    }

    c.textAlign = "center";
    c.fillText(
      xLabel,
      box.x + box.w / 2,
      height - 6
    );

    c.save();
    c.translate(12, box.y + box.h / 2);
    c.rotate(-Math.PI / 2);
    c.fillText(yLabel, 0, 0);
    c.restore();

    c.save();
    c.beginPath();
    c.rect(box.x, box.y, box.w, box.h);
    c.clip();

    return {
      c,
      map,
      box,
      finish: () => c.restore()
    };
  }

  function drawSeries(plot, points, tint, selected = -1) {
    const {c, map} = plot;

    c.strokeStyle = tint;
    c.lineWidth = 1.8;
    c.beginPath();

    points.forEach((p, i) => {
      if (!p) return;

      if (i && points[i - 1]) {
        c.lineTo(...map(p));
      } else {
        c.moveTo(...map(p));
      }
    });

    c.stroke();
    c.fillStyle = tint;

    points.forEach((p, i) => {
      if (!p) return;

      c.beginPath();

      c.arc(
        ...map(p),
        i === selected ? 4.5 : 2.5,
        0,
        2 * Math.PI
      );

      c.fill();
    });
  }

  /* Magnified image and trajectory */

  function renderLens() {
    const chosen = s.nodes[s.step];

    const trace = s.nodes.map(node =>
      good(node) ? values(node) : null
    );

    const mesh = s.mesh?.samples.map(meshValue) ?? [];
    const available = [...mesh, ...trace].filter(Boolean);

    const width = Math.max(
      200,
      $("limit-lens").getBoundingClientRect().width
    );

    const aspect = (width - 88) / 272;

    if (!s.view) {
      const origin =
        values(chosen) ??
        trace.find(Boolean) ??
        available[0] ??
        ["0", "0"];

      let fitted = s.fit === "image"
        ? mesh.filter(Boolean)
        : trace.filter(Boolean);

      if (s.fit === "tail") {
        fitted = s.nodes
          .slice(Math.max(0, s.step - 3), s.step + 1)
          .filter(good)
          .map(node => values(node));
      }

      if (!fitted.length) fitted = available;

      const offsets = fitted.map(z => offset(z, origin));
      const xs = offsets.map(p => p[0]);
      const ys = offsets.map(p => p[1]);

      const x0 = xs.length ? Math.min(...xs) : -0.5;
      const x1 = xs.length ? Math.max(...xs) : 0.5;
      const y0 = ys.length ? Math.min(...ys) : -0.5;
      const y1 = ys.length ? Math.max(...ys) : 0.5;

      s.view = {
        origin,
        x: (x0 + x1) / 2,
        y: (y0 + y1) / 2,
        span:
          Math.max(
            y1 - y0,
            (x1 - x0) / aspect
          ) * 1.2 || 0.001
      };
    }

    const v = s.view;

    const span =
      v.span / 10 ** Number($("limit-zoom").value);

    const variable = s.chart === "raw" ? "τ" : "w";

    const bounds = [
      v.x - span * aspect / 2,
      v.x + span * aspect / 2,
      v.y - span / 2,
      v.y + span / 2
    ];

    const plot = frame(
      $("limit-lens"),
      bounds,
      "Δ Re " + variable,
      "Δ Im " + variable
    );

    /*
      Ordinary-precision mesh samples cannot justify arbitrarily
      fine detail. At small scales show the checked trace alone.
    */
    const showMesh = span >= 1e-6;

    if (s.mesh && showMesh) {
      plot.c.fillStyle = color + "28";
      plot.c.strokeStyle = color + "55";
      plot.c.lineWidth = 0.5;

      for (const face of s.mesh.triangles) {
        const points = face.map(i =>
          mesh[i]
            ? plot.map(offset(mesh[i], v.origin))
            : null
        );

        if (points.some(p => !p)) continue;

        plot.c.beginPath();
        plot.c.moveTo(...points[0]);
        plot.c.lineTo(...points[1]);
        plot.c.lineTo(...points[2]);
        plot.c.closePath();
        plot.c.fill();
        plot.c.stroke();
      }
    }

    drawSeries(
      plot,
      trace.map(z => z ? offset(z, v.origin) : null),
      color,
      s.step
    );

    plot.finish();

    $("limit-zoom-value").textContent =
      fmt(10 ** Number($("limit-zoom").value)) + "×";

    const error = chosen?.result?.viewErrors?.[s.chart];

    $("limit-lens-note").textContent =
      `Origin ${variable}* = ` +
      `(${v.origin[0]}, ${v.origin[1]}). ` +
      `Visible height: ${fmt(span)}. ` +
      (
        !showMesh
          ? "The ordinary-precision mesh is hidden at this scale; " +
            "only the checked trajectory is shown. "
          : ""
      ) +
      (
        error && span < 10 * Number(error)
          ? "This scale is below the numerical precision diagnostic; " +
            "further zoom adds no verified detail."
          : "Axes show offsets from this origin."
      );
  }

  /* Behaviour against decreasing epsilon */

  function renderTrend() {
    const chosen = s.nodes[lastGood()];
    const metric = $("limit-metric").value;
    const origin = values(chosen) ?? ["0", "0"];

    function series(branch) {
      return s.nodes.map((node, i) => {
        const z = values(node, branch);

        if (!good(node) || !z) return null;

        let y;

        if (metric === "change") {
          const prior = i
            ? values(s.nodes[i - 1], branch)
            : null;

          if (!prior || !good(s.nodes[i - 1])) {
            return null;
          }

          const distance = Math.hypot(...offset(z, prior));

          if (!(distance > 0)) return null;

          y = Math.log10(distance);
        } else {
          const index = metric === "real" ? 0 : 1;
          y = difference(z[index], origin[index]);
        }

        return [
          -Math.log10(Number(node.epsilon)),
          y
        ];
      });
    }

    const result = series("result");

    const base = $("limit-compare").checked
      ? series("base")
      : [];

    const points = [...result, ...base].filter(Boolean);
    const ys = points.map(p => p[1]);

    let low = ys.length ? Math.min(...ys) : 0;
    let high = ys.length ? Math.max(...ys) : 1;

    const pad =
      (high - low) * 0.12 ||
      Math.max(Math.abs(low) * 0.1, 1e-12);

    low -= pad;
    high += pad;

    const right = Math.max(
      3,
      ...points.map(p => p[0])
    );

    const variable = s.chart === "raw" ? "τ" : "w";

    const label = metric === "change"
      ? `log₁₀ |Δ${variable}|`
      : (metric === "real" ? "Δ Re " : "Δ Im ") + variable;

    const plot = frame(
      $("limit-trend"),
      [1, right, low, high],
      "−log₁₀ ε",
      label
    );

    drawSeries(plot, base, baseColor, s.step);
    drawSeries(plot, result, color, s.step);

    plot.finish();

    $("limit-trend-note").textContent =
      metric === "change"
        ? "Consecutive computed changes, on a logarithmic vertical " +
          "scale. These are not errors from an assumed limit."
        : "Vertical origin is the last checked result: " +
          `${metric === "real" ? "Re" : "Im"} ${variable}* = ` +
          `${origin[metric === "real" ? 0 : 1]}. ` +
          "Orange: result; teal: optional n = 1 comparison.";
  }

  /* Sample selection and table */

  function render() {
    renderLens();
    renderTrend();

    const node = s.nodes[s.step];

    $("limit-step").max = Math.max(0, s.nodes.length - 1);
    $("limit-step").value = Math.max(0, s.step);
    $("limit-step").disabled = !s.nodes.length;

    $("limit-step-label").textContent = node
      ? `Sample ${s.step + 1}/${s.nodes.length}; ε = ${node.epsilon}`
      : "No checked sample selected";

    $("limit-selected").textContent =
      node?.result?.tauExact
        ? `τ = (${node.result.tauExact[0]}, ` +
          `${node.result.tauExact[1]}). ` +
          "Estimated numerical error: " +
          `${node.result.numericalErrorEstimate}.`
        : "";

    $("limit-table").innerHTML = s.nodes.map((n, i) => {
      const resolved = good(n);

      return `
        <tr class="${i === s.step ? "selected" : ""}">
          <td>
            <button type="button"
                    data-step="${i}"
                    ${resolved ? "" : "disabled"}>
              ${escape(n.epsilon)}
            </button>
          </td>
          <td>
            ${resolved ? fmt(n.result.tauExact[0], 14) : "Unresolved"}
          </td>
          <td>
            ${resolved ? fmt(n.result.tauExact[1], 14) : "—"}
          </td>
          <td>
            ${
              resolved
                ? escape(fmt(n.result.numericalErrorEstimate))
                : "—"
            }
          </td>
        </tr>
      `;
    }).join("");

    $("limit-export").disabled = !s.nodes.length;
  }

  /* Zoom, pan, and selection */

  function pan(dx, dy) {
    if (!s.view) return;

    const span =
      s.view.span / 10 ** Number($("limit-zoom").value);

    const aspect =
      (
        Math.max(
          200,
          $("limit-lens").getBoundingClientRect().width
        ) - 88
      ) / 272;

    s.view.x += dx * span * aspect;
    s.view.y += dy * span;

    renderLens();
  }

  function centerSelected() {
    if (!s.view) renderLens();

    const z = values(s.nodes[s.step]);

    if (z && s.view) {
      [s.view.x, s.view.y] = offset(z, s.view.origin);
    }

    renderLens();
  }

  function selectStep(index) {
    s.step = index;

    if (s.fit === "tail") s.view = null;

    render();

    if (Number($("limit-zoom").value) > 0) {
      centerSelected();
    }
  }

  $("limits-controls").addEventListener("submit", run);
  $("limits-cancel").onclick = cancel;
  $("limit-path").onchange = pathNote;

  $("limit-chart").onchange = () => {
    s.chart = $("limit-chart").value;
    s.view = null;
    $("limit-zoom").value = 0;
    render();
  };

  $("limit-zoom").oninput = centerSelected;
  $("limit-center").onclick = centerSelected;
  $("limit-metric").onchange = renderTrend;
  $("limit-compare").onchange = renderTrend;

  document.querySelectorAll("[data-fit]").forEach(button => {
    button.onclick = () => {
      s.fit = button.dataset.fit;
      s.view = null;
      $("limit-zoom").value = 0;
      render();
    };
  });

  document.querySelectorAll("[data-pan]").forEach(button => {
    button.onclick = () => {
      pan(
        ...button.dataset.pan
          .split(",")
          .map(x => Number(x) * 0.15)
      );
    };
  });

  $("limit-step").oninput = () => {
    selectStep(Number($("limit-step").value));
  };

  $("limit-table").onclick = event => {
    const button = event.target.closest("button[data-step]");

    if (button) {
      selectStep(Number(button.dataset.step));
    }
  };

  let drag = null;

  $("limit-lens").onpointerdown = event => {
    $("limit-lens").setPointerCapture(event.pointerId);

    drag = {
      x: event.clientX,
      y: event.clientY
    };
  };

  $("limit-lens").onpointermove = event => {
    if (!drag) return;

    const width = Math.max(
      200,
      $("limit-lens").getBoundingClientRect().width
    );

    pan(
      (drag.x - event.clientX) / (width - 88),
      (event.clientY - drag.y) / 272
    );

    drag = {
      x: event.clientX,
      y: event.clientY
    };
  };

  $("limit-lens").onpointerup =
    $("limit-lens").onpointercancel =
      () => {
        drag = null;
      };

  $("limit-lens").onkeydown = event => {
    const delta = {
      ArrowLeft: [-0.15, 0],
      ArrowRight: [0.15, 0],
      ArrowUp: [0, 0.15],
      ArrowDown: [0, -0.15]
    }[event.key];

    if (delta) {
      event.preventDefault();
      pan(...delta);
    }
  };

  /* Export the stored decimal coordinates and diagnostics */

  $("limit-export").onclick = () => {
    const data = {
      ...s.request,
      ...(s.result ?? {}),
      status:
        s.result?.status ??
        (s.pending.has("trace") ? "in_progress" : "stopped"),
      nodes: s.nodes,
      displayChart: s.chart,
      note:
        "Numerical estimates checked at two working precisions; " +
        "no exact-limit or interval certificate."
    };

    const blob = new Blob(
      [JSON.stringify(data, null, 2)],
      {type: "application/json"}
    );

    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");

    link.href = url;
    link.download = "hemisphere-limit-trace.json";
    link.click();

    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };

  window.addEventListener("resize", () => {
    s.view = null;
    render();
  });

  pathNote();
  render();
  run();
})();

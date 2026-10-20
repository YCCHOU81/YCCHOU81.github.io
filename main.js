"use strict";

/* Shared navigation and attachment selection */

const countNames = ["a", "b", "c"];

const order = counts =>
  1 + counts.reduce((sum, value) => sum + value, 0);

const validCounts = counts =>
  counts.length === 3 &&
  counts.every(value => Number.isSafeInteger(value) && value >= 0) &&
  Number.isSafeInteger(order(counts));

const query = new URLSearchParams(location.search);

let selectedCounts = countNames.map((name, i) =>
  query.has(name) ? Number(query.get(name)) : [1, 0, 0][i]
);

if (!validCounts(selectedCounts)) {
  selectedCounts = [1, 0, 0];
}

function refreshSharedInformation() {
  const params = new URLSearchParams();

  countNames.forEach((name, i) => {
    params.set(name, selectedCounts[i]);
  });

  document.querySelectorAll("a[data-page]").forEach(link => {
    link.href = `${link.dataset.page}?${params}`;
  });

  const context = document.getElementById("selection-context");

  if (context) {
    context.textContent =
      `Selected attachment: (${selectedCounts.join(", ")}). ` +
      `Order n = ${order(selectedCounts)}.`;
  }
}

refreshSharedInformation();

if (document.getElementById("attachment-form")) {
  connectExplorer();
}

/* Explore controls and numerical worker */

function connectExplorer() {
  const get = id => document.getElementById(id);
  const status = get("attachment-status");
  const cache = new Map();

  let worker = null;
  let request = 0;
  let currentMesh = null;

  countNames.forEach((name, i) => {
    get(name).value = selectedCounts[i];
  });

  function message(text, error = false) {
    status.textContent = text;
    status.classList.toggle("error", error);
  }

  function stop() {
    worker?.terminate();
    worker = null;

    get("compute").disabled = false;
    get("cancel").disabled = true;
  }

  function redraw() {
    drawImage(
      get("base-image"),
      currentMesh,
      "base",
      "#285f65"
    );

    drawImage(
      get("result-image"),
      currentMesh,
      "result",
      "#b55e3c"
    );
  }

  function accept(mesh, elapsed, cached = false) {
    currentMesh = mesh;

    const count = name =>
      mesh.samples.filter(sample =>
        sample?.[name]?.status === "ok"
      ).length;

    const total = mesh.points.length;

    message(
      `Resolved samples: base ${count("base")}/${total}; ` +
      `result ${count("result")}/${total}. ` +
      (
        cached
          ? "Loaded from this tab's cache."
          : `${elapsed.toFixed(2)} seconds.`
      ) +
      (
        mesh.status === "ok"
          ? ""
          : " Unresolved samples are omitted."
      )
    );

    redraw();
  }

  function compute(event) {
    event?.preventDefault();

    const form = get("attachment-form");

    if (!form.reportValidity()) {
      message(
        "Enter whole numbers from 0 through 11 in all three fields.",
        true
      );
      return;
    }

    const counts = countNames.map(name =>
      get(name).valueAsNumber
    );

    if (!validCounts(counts) || order(counts) > 12) {
      message(
        "Enter nonnegative whole numbers with a + b + c ≤ 11.",
        true
      );
      return;
    }

    stop();

    const id = ++request;
    const resolution = Number(get("quality").value);
    const key = `${counts.join(",")}/${resolution}`;

    selectedCounts = counts;
    refreshSharedInformation();

    const url = new URL(location.href);

    countNames.forEach((name, i) => {
      url.searchParams.set(name, counts[i]);
    });

    history.replaceState(null, "", url);

    get("result-heading").textContent =
      `Result: (${counts.join(", ")}), n = ${order(counts)}`;

    currentMesh = null;
    redraw();

    if (cache.has(key)) {
      accept(cache.get(key), 0, true);
      return;
    }

    get("compute").disabled = true;
    get("cancel").disabled = false;

    message("Computing the two image regions…");

    const started = performance.now();

    try {
      worker = new Worker("numerics-worker.js");

      worker.onmessage = ({data}) => {
        if (data.id !== request) return;

        if (data.type === "progress") {
          message(
            `Computing samples: ${data.completed}/${data.total}…`
          );
        } else if (data.type === "done") {
          stop();

          cache.set(key, data.mesh);

          if (cache.size > 8) {
            cache.delete(cache.keys().next().value);
          }

          accept(
            data.mesh,
            (performance.now() - started) / 1000
          );
        } else if (data.type === "error") {
          stop();
          message(data.message, true);
        }
      };

      worker.onerror = () => {
        if (id !== request) return;

        stop();

        message(
          "Numerical engine failed to load or run. " +
          "Check that numerics-worker.js is beside explore.html " +
          "and open the site through GitHub Pages or a local web server.",
          true
        );
      };

      worker.postMessage({
        type: "compute",
        id,
        counts,
        resolution,
        inset: 0.04
      });
    } catch (error) {
      stop();

      message(
        `Cannot start the numerical worker: ${error.message}`,
        true
      );
    }
  }

  get("attachment-form").addEventListener("submit", compute);

  get("cancel").addEventListener("click", () => {
    ++request;
    stop();
    message("Computation cancelled.");
  });

  window.addEventListener("resize", redraw);

  redraw();
  compute();
}

/* Draw one conformal image in the marked chart */

function drawImage(canvas, mesh, field, color) {
  const width = canvas.getBoundingClientRect().width;
  const height = 340;

  if (width < 1) return;

  const ratio = Math.min(devicePixelRatio || 1, 2);

  canvas.width = Math.round(width * ratio);
  canvas.height = Math.round(height * ratio);

  const ctx = canvas.getContext("2d");

  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);

  const left = 44;
  const top = 16;
  const w = width - 60;
  const h = height - 62;

  const px = x => left + x * w;
  const py = u => top + (1 - u) * h;

  // Compactify height: y becomes y / (1 + y).
  const point = z => [
    px(z.re),
    py(z.im / (1 + z.im))
  ];

  const line = (a, b) => {
    ctx.moveTo(...a);
    ctx.lineTo(...b);
  };

  // Boundary of the marked 0 < Re τ < 1 chart.
  const chart = new Path2D();

  chart.moveTo(px(0), py(1));

  for (let i = 0; i <= 120; i++) {
    const x = i / 120;
    const y = Math.sqrt(x * (1 - x));

    chart.lineTo(...point({re: x, im: y}));
  }

  chart.lineTo(px(1), py(1));
  chart.closePath();

  ctx.fillStyle = "#eff3e9";
  ctx.fill(chart);

  ctx.save();
  ctx.clip(chart);

  // Coordinate grid.
  ctx.strokeStyle = "#d6ded2";
  ctx.lineWidth = 1;
  ctx.beginPath();

  for (const y of [0.5, 1, 2, 4]) {
    line(
      point({re: 0, im: y}),
      point({re: 1, im: y})
    );
  }

  for (const x of [0.25, 0.5, 0.75]) {
    line(
      [px(x), py(0)],
      [px(x), py(1)]
    );
  }

  ctx.stroke();

  if (mesh) {
    const points = mesh.samples.map(sample => {
      const branch = sample?.[field];

      if (
        branch?.status !== "ok" ||
        !Number.isFinite(branch.tau?.re) ||
        !Number.isFinite(branch.tau?.im)
      ) {
        return null;
      }

      return point(branch.tau);
    });

    ctx.fillStyle = color + "30";
    ctx.strokeStyle = color + "65";
    ctx.lineWidth = 0.5;

    const edges = new Map();

    for (const face of mesh.triangles) {
      // Record source-mesh edges to find its inset perimeter.
      for (let i = 0; i < 3; i++) {
        const a = face[i];
        const b = face[(i + 1) % 3];

        const key = [
          Math.min(a, b),
          Math.max(a, b)
        ].join(",");

        edges.set(key, (edges.get(key) || 0) + 1);
      }

      const p = face.map(i => points[i]);

      // Do not fill across unresolved samples.
      if (p.some(value => !value)) continue;

      ctx.beginPath();
      ctx.moveTo(...p[0]);
      ctx.lineTo(...p[1]);
      ctx.lineTo(...p[2]);
      ctx.closePath();

      ctx.fill();
      ctx.stroke();
    }

    // Dashed perimeter of the sampled interior.
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();

    for (const [key, count] of edges) {
      if (count !== 1) continue;

      const [a, b] = key.split(",").map(Number);

      if (points[a] && points[b]) {
        line(points[a], points[b]);
      }
    }

    ctx.stroke();
  }

  ctx.restore();

  ctx.strokeStyle = "#94a68e";
  ctx.lineWidth = 1;
  ctx.stroke(chart);

  // Axis labels.
  ctx.fillStyle = "#596e62";
  ctx.font = "12px system-ui";
  ctx.textAlign = "center";

  for (const x of [0, 0.5, 1]) {
    ctx.fillText(x, px(x), height - 25);
  }

  ctx.fillText("Re τ", left + w / 2, height - 5);
  ctx.textAlign = "right";

  for (const y of [0.5, 1, 2, 4]) {
    ctx.fillText(
      y,
      left - 8,
      point({re: 0, im: y})[1] + 4
    );
  }

  ctx.fillText("∞", left - 8, top + 5);

  ctx.save();
  ctx.translate(12, top + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText("Im τ", 0, 0);
  ctx.restore();
}

/* Limits retains its settings-only skeleton. */

const limitForm = document.getElementById("limit-form");

if (limitForm) {
  limitForm.addEventListener("submit", event => {
    event.preventDefault();

    const zoom = document.getElementById("zoom").valueAsNumber;
    const status = document.getElementById("limit-status");
    const valid = Number.isSafeInteger(zoom) && zoom >= 1;

    status.classList.toggle("error", !valid);

    status.textContent = valid
      ? "Settings selected. Numerical continuation is not connected yet."
      : "Enter a positive whole-number zoom factor.";

    if (valid) {
      const approach = document.getElementById("approach").value;

      document.getElementById("limit-placeholder").textContent =
        `${approach} approach; requested zoom: ${zoom}×. ` +
        "Visualization pending.";
    }
  });
}

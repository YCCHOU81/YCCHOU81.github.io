"use strict";

/*
  Shared skeleton behavior.

  - Reads attachment counts from the URL.
  - Carries those counts between pages.
  - Handles forms only on pages where they exist.
  - Leaves all numerical calculations and graphs as placeholders.
*/

const parameterNames = ["a", "b", "c"];
const defaultCounts = [1, 0, 0];

function getOrder(counts) {
  return 1 + counts.reduce((sum, value) => sum + value, 0);
}

function validCounts(counts) {
  return (
    counts.length === 3 &&
    counts.every(value => Number.isSafeInteger(value) && value >= 0) &&
    Number.isSafeInteger(getOrder(counts))
  );
}

function readCountsFromURL() {
  const query = new URLSearchParams(window.location.search);

  const counts = parameterNames.map((name, index) => {
    return query.has(name)
      ? Number(query.get(name))
      : defaultCounts[index];
  });

  return validCounts(counts) ? counts : [...defaultCounts];
}

let selectedCounts = readCountsFromURL();

function refreshSharedInformation() {
  const [a, b, c] = selectedCounts;
  const query = new URLSearchParams({
    a: String(a),
    b: String(b),
    c: String(c)
  });

  // Keep the selected attachment when navigating between pages.
  document.querySelectorAll("a[data-page]").forEach(link => {
    link.href = `${link.dataset.page}?${query.toString()}`;
  });

  const context = document.getElementById("selection-context");

  if (context) {
    context.textContent =
      `Selected attachment: (${a}, ${b}, ${c}). ` +
      `Order n = ${getOrder(selectedCounts)}.`;
  }
}

/* Explore page */
const parameterForm = document.getElementById("parameter-form");

if (parameterForm) {
  parameterNames.forEach((name, index) => {
    document.getElementById(name).value = selectedCounts[index];
  });

  parameterForm.addEventListener("submit", event => {
    event.preventDefault();

    const counts = parameterNames.map(name => {
      return document.getElementById(name).valueAsNumber;
    });

    const status = document.getElementById("parameter-status");

    if (!validCounts(counts)) {
      status.textContent =
        "Please enter smaller nonnegative whole numbers.";
      status.classList.add("error");
      return;
    }

    selectedCounts = counts;
    refreshSharedInformation();

    status.classList.remove("error");
    status.textContent = "Attachment selected. Visualization code goes next.";

    updateComparison({
      a: counts[0],
      b: counts[1],
      c: counts[2],
      n: getOrder(counts)
    });
  });
}

function updateComparison(parameters) {
  // TODO: connect the solver and comparison plots here.
  // This placeholder does not compute a period ratio.

  const message = document.getElementById("comparison-message");

  if (message) {
    message.textContent =
      `Ready for the n = ${parameters.n} comparison. ` +
      "Add the mathematical calculation and drawing code here.";
  }
}

/* Limits page */
const limitForm = document.getElementById("limit-form");

if (limitForm) {
  limitForm.addEventListener("submit", event => {
    event.preventDefault();

    const approach = document.getElementById("approach").value;
    const depth = document.getElementById("depth").valueAsNumber;
    const status = document.getElementById("limit-status");

    if (!Number.isInteger(depth) || depth < 3 || depth > 10) {
      status.textContent = "Choose a whole-number depth from 3 to 10.";
      status.classList.add("error");
      return;
    }

    status.classList.remove("error");
    status.textContent =
      `Selected ${approach} approach, depth ${depth}.`;

    prepareLimit({
      counts: [...selectedCounts],
      approach,
      depth
    });
  });
}

function prepareLimit(settings) {
  // TODO: connect numerical continuation, zoom, and trend plots here.
  // This placeholder does not compute a trajectory.

  const message = document.getElementById("limit-message");

  if (message) {
    message.textContent =
      `Ready for an ${settings.approach} approach at depth ` +
      `${settings.depth}. Add the continuation and plotting code here.`;
  }
}

refreshSharedInformation();

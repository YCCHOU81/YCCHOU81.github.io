"use strict";

/*
  Shared by all four pages.

  This skeleton:
  - preserves attachment counts between pages;
  - validates controls;
  - updates placeholder messages.

  It does not compute triangles or conformal structures.
*/

const countNames = ["a", "b", "c"];
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

function readCounts() {
  const query = new URLSearchParams(window.location.search);

  const counts = countNames.map((name, index) => {
    return query.has(name)
      ? Number(query.get(name))
      : defaultCounts[index];
  });

  return validCounts(counts) ? counts : [...defaultCounts];
}

let selectedCounts = readCounts();

function refreshSharedInformation() {
  const [a, b, c] = selectedCounts;

  const query = new URLSearchParams({
    a: String(a),
    b: String(b),
    c: String(c)
  });

  // These remain ordinary links to separate HTML pages.
  document.querySelectorAll("a[data-page]").forEach(link => {
    link.href = `${link.dataset.page}?${query.toString()}`;
  });

  const context = document.getElementById("selection-context");

  if (context) {
    context.textContent =
      `Selected attachment: (${a}, ${b}, ${c}). ` +
      `Order n = ${getOrder(selectedCounts)}.`;
  }

  const result = document.getElementById("result-placeholder");

  if (result) {
    result.textContent =
      `Resulting triangle: n = ${getOrder(selectedCounts)}. ` +
      "Plot pending.";
  }
}

/* Explore page */

const attachmentForm = document.getElementById("attachment-form");

if (attachmentForm) {
  countNames.forEach((name, index) => {
    document.getElementById(name).value = selectedCounts[index];
  });

  attachmentForm.addEventListener("submit", event => {
    event.preventDefault();

    const counts = countNames.map(name => {
      return document.getElementById(name).valueAsNumber;
    });

    const status = document.getElementById("attachment-status");

    if (!validCounts(counts)) {
      status.textContent =
        "Enter nonnegative whole numbers within the supported number range.";
      status.classList.add("error");
      return;
    }

    selectedCounts = counts;
    refreshSharedInformation();

    // Save the selection in this page's address for refresh or sharing.
    const url = new URL(window.location.href);

    countNames.forEach((name, index) => {
      url.searchParams.set(name, String(counts[index]));
    });

    window.history.replaceState(null, "", url.href);

    status.classList.remove("error");
    status.textContent =
      "Attachment selected. The numerical solver is not connected yet.";

    // TODO: call the solver and draw the comparison here.
  });
}

/* Limits page */

const limitForm = document.getElementById("limit-form");

if (limitForm) {
  limitForm.addEventListener("submit", event => {
    event.preventDefault();

    const approach = document.getElementById("approach").value;
    const zoom = document.getElementById("zoom").valueAsNumber;

    const status = document.getElementById("limit-status");
    const placeholder = document.getElementById("limit-placeholder");

    if (!Number.isSafeInteger(zoom) || zoom < 1) {
      status.textContent = "Enter a positive whole-number zoom factor.";
      status.classList.add("error");
      return;
    }

    status.classList.remove("error");
    status.textContent =
      "Settings selected. Numerical continuation is not connected yet.";

    placeholder.textContent =
      `${approach === "edge" ? "Edge" : "Vertex"} approach; ` +
      `requested zoom: ${zoom}×. Visualization pending.`;

    // TODO: compute the approach and draw the magnified view here.
  });
}

refreshSharedInformation();

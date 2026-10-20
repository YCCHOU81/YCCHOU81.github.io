"use strict";

const form = document.getElementById("parameter-form");
const statusElement = document.getElementById("status");
const plotMessage = document.getElementById("plot-message");

form.addEventListener("submit", function (event) {
  event.preventDefault();

  const values = ["a", "b", "c"].map(function (id) {
    return document.getElementById(id).valueAsNumber;
  });

  const valid = values.every(function (value) {
    return Number.isSafeInteger(value) && value >= 0;
  });

  const n = 1 + values.reduce(function (sum, value) {
    return sum + value;
  }, 0);

  if (!valid || !Number.isSafeInteger(n)) {
    statusElement.textContent =
      "Please enter smaller nonnegative whole numbers.";

    statusElement.classList.add("error");
    return;
  }

  const [a, b, c] = values;

  statusElement.classList.remove("error");
  statusElement.textContent =
    `Applied (${a}, ${b}, ${c}). Order n = ${n}.`;

  updateVisualization({ a, b, c, n });
});

function updateVisualization(parameters) {
  // TODO: connect the mathematical solver here.
  // TODO: draw its results inside the #plot container.
  // This placeholder does not compute a period ratio.

  plotMessage.textContent =
    `Parameters received for n = ${parameters.n}. ` +
    "Add your visualization here.";
}

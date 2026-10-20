// EDIT: connect button actions and your mathematical code here.
const form = document.getElementById('parameter-form');
const summary = document.getElementById('parameter-summary');

form.addEventListener('submit', event => {
  event.preventDefault();
  const values = ['a', 'b', 'c'].map(id => Number(document.getElementById(id).value));

  if (!values.every(value => Number.isSafeInteger(value) && value >= 0) ||
      !Number.isSafeInteger(1 + values.reduce((sum, value) => sum + value, 0))) {
    summary.textContent = 'Please enter three nonnegative whole numbers in a smaller range.';
    summary.classList.add('error');
    return;
  }

  const [a, b, c] = values;
  const n = 1 + a + b + c;
  summary.classList.remove('error');
  summary.textContent = `Applied (${a}, ${b}, ${c}). Order n = ${n}.`;
  updateVisualization({ a, b, c, n });
});

function updateVisualization(parameters) {
  // TODO: call your solver and draw its result in the #plot container.
  // This skeleton deliberately performs no period-ratio computation.
  document.getElementById('plot-message').textContent =
    `Parameters received for n = ${parameters.n}. Add your visualization here.`;
}

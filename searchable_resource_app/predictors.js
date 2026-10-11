// Set this to the Zenodo record ID once the predictor FASTAs are uploaded.
// A new Zenodo version gets a new ID, so this is the only line to change.
const ZENODO_RECORD_ID = null;

const EMBEDDER_COLORS = {
  "One-hot":   "#A8A8A8",
  "ESM2 8M":   "#F2B38A",
  "ESM2 650M": "#E07840",
  "ProtT5":    "#1B4F72",
};
const EMBEDDER_ORDER = ["One-hot", "ESM2 8M", "ESM2 650M", "ProtT5"];

const MODEL_NAMES = {
  fiveclass_split1_30pct:   "5-class · split 1 · 30%",
  fiveclass_split1_20pct:   "5-class · split 1 · 20%",
  twelveclass_split2_30pct: "12-class · split 2 · 30%",
  twelveclass_split2_20pct: "12-class · split 2 · 20%",
  binary_clustered:         "Binary · extremophile vs. non-extremophile",
  sixclass_clustered:       "6-class · 5 classes + non-extremophile",
  algae_filtered:           "Eukaryotic algae · test only",
};

const fmt = (v, digits = 3) => (v === null || v === undefined ? "–" : Number(v).toFixed(digits));
const count = (n) => Number(n).toLocaleString("en-US");

function zenodoUrl(file) {
  return ZENODO_RECORD_ID ? `https://zenodo.org/records/${ZENODO_RECORD_ID}/files/${file}?download=1` : null;
}

async function loadJson(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Failed to load ${path}: ${response.status}`);
  return response.json();
}

// ── Top metric cards ─────────────────────────────────────────────────────────
function renderMetrics(summary) {
  const binary = summary.binary_clustered.ProtT5;
  const bacteriaOnly = summary.binary_bacteria_only.ProtT5;
  const algae = summary.algae.ProtT5;
  const homology = summary.homology.binary_clustered;
  const cards = [
    ["Binary MCC (ProtT5)", fmt(binary.mcc.value), `accuracy ${fmt(binary.accuracy.value)} on clustered prokaryote test set; ${fmt(bacteriaOnly.mcc.value)} on bacteria only`],
    ["5-class macro-F1", fmt(summary.fiveclass.fiveclass_split1_30pct.ProtT5.macro_f1.value), "ProtT5, split 1, 30% clustering"],
    ["Algae recall", `${(algae.overall * 100).toFixed(1)}%`, "binary ProtT5 model on 133,689 algae proteins"],
    ["Without any homolog", fmt(homology.prott5_accuracy_without_hit), `ProtT5 accuracy where homology finds nothing (homology overall ${fmt(homology.homology_accuracy)})`],
  ];
  const grid = document.getElementById("predictor-metrics");
  grid.innerHTML = cards.map(([label, value, note]) => `
    <article class="metric-card">
      <p class="metric-label">${label}</p>
      <strong class="metric-value">${value}</strong>
      <p class="metric-note">${note}</p>
    </article>`).join("");
}

// ── Model cards with downloads ───────────────────────────────────────────────
function renderModelCards(manifest) {
  const cards = Object.entries(manifest).map(([name, model]) => {
    const rows = Object.entries(model.splits).map(([split, info]) => {
      const url = zenodoUrl(info.file);
      const button = url
        ? `<a class="button" href="${url}">FASTA</a>`
        : `<span class="button" aria-disabled="true" title="Available after the Zenodo upload">FASTA</span>`;
      const labels = Object.entries(info.labels).map(([l, n]) => `${l.replace("_", " ")} ${count(n)}`).join(" · ");
      return `<tr>
        <td><strong>${split}</strong></td>
        <td>${count(info.n_proteins)}</td>
        <td title="${labels}">${info.size_mb} MB</td>
        <td>${button}</td>
      </tr>`;
    }).join("");
    return `<article class="model-card">
      <h3>${MODEL_NAMES[name] || name}</h3>
      <p>${model.description}</p>
      <table>
        <thead><tr><th>Split</th><th>Proteins</th><th>Size</th><th></th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </article>`;
  });
  document.getElementById("model-cards").innerHTML = cards.join("");
  document.getElementById("zenodo-status").textContent =
    ZENODO_RECORD_ID ? `Zenodo record ${ZENODO_RECORD_ID}` : "Zenodo (upload pending)";
}

// ── Charts ───────────────────────────────────────────────────────────────────
function ciTooltip(ciLookup) {
  return {
    callbacks: {
      label: (ctx) => {
        const ci = ciLookup(ctx);
        const base = `${ctx.dataset.label}: ${fmt(ctx.parsed.y)}`;
        return ci && ci[0] !== null ? `${base}  (95% CI ${fmt(ci[0])}–${fmt(ci[1])})` : base;
      },
    },
  };
}

function barChart(id, labels, datasets, { yLabel, ciLookup, min = 0 } = {}) {
  const canvas = document.getElementById(id);
  if (!canvas || typeof Chart === "undefined") return;
  new Chart(canvas, {
    type: "bar",
    data: { labels, datasets },
    options: {
      responsive: true,
      plugins: {
        legend: { position: "bottom", labels: { boxWidth: 12, font: { size: 11 } } },
        tooltip: ciLookup ? ciTooltip(ciLookup) : {},
      },
      scales: {
        y: { min, max: 1, title: { display: !!yLabel, text: yLabel } },
        x: { grid: { display: false } },
      },
    },
  });
}

function renderCharts(summary) {
  const fiveSplits = Object.keys(summary.fiveclass);
  barChart("chart-fiveclass", fiveSplits.map((s) => MODEL_NAMES[s].replace(" · ", " ")),
    EMBEDDER_ORDER.map((e) => ({
      label: e,
      backgroundColor: EMBEDDER_COLORS[e],
      data: fiveSplits.map((s) => summary.fiveclass[s][e].macro_f1.value),
    })),
    { yLabel: "Macro-F1", ciLookup: (ctx) => {
      const m = summary.fiveclass[fiveSplits[ctx.dataIndex]][ctx.dataset.label].macro_f1;
      return [m.ci_lower, m.ci_upper];
    } });

  const binaryMetrics = ["accuracy", "specificity", "recall", "mcc"];
  barChart("chart-binary", ["Accuracy", "Specificity", "Recall", "MCC"],
    EMBEDDER_ORDER.map((e) => ({
      label: e,
      backgroundColor: EMBEDDER_COLORS[e],
      data: binaryMetrics.map((m) => summary.binary_clustered[e][m].value),
    })),
    { yLabel: "Score", ciLookup: (ctx) => {
      const m = summary.binary_clustered[ctx.dataset.label][binaryMetrics[ctx.dataIndex]];
      return [m.ci_lower, m.ci_upper];
    } });

  const bias = summary.binary_bacteria_only.ProtT5.mcc.value;
  document.getElementById("binary-bias-note").textContent =
    `The non-extremophiles are bacteria only, so archaeal proteins are easy to call extremophile. Without the ` +
    `archaeal extremophiles, ProtT5 reaches an MCC of ${fmt(bias)} instead of ${fmt(summary.binary_clustered.ProtT5.mcc.value)}.`;

  const algaeGroups = ["overall", "acidophile", "halophile", "psychrophile", "thermophile"];
  barChart("chart-algae", algaeGroups.map((g) => g[0].toUpperCase() + g.slice(1)),
    EMBEDDER_ORDER.map((e) => ({
      label: e,
      backgroundColor: EMBEDDER_COLORS[e],
      data: algaeGroups.map((g) => summary.algae[e][g]),
    })),
    { yLabel: "Recall" });

  const detection = summary.sixclass_clustered;
  const detectionModels = Object.keys(detection);
  const detectionColors = ["#C9C9C9", "#9DB9DC", "#3B6EA8"];
  barChart("chart-detection", ["Accuracy", "Specificity", "Recall", "MCC"],
    detectionModels.map((name, i) => ({
      label: name,
      backgroundColor: detectionColors[i],
      data: binaryMetrics.map((m) => detection[name][m].value),
    })),
    { yLabel: "Score", ciLookup: (ctx) => {
      const m = detection[ctx.dataset.label][binaryMetrics[ctx.dataIndex]];
      return [m.ci_lower, m.ci_upper];
    } });

  const hom = summary.homology;
  const tasks = ["fiveclass_split1_30pct", "binary_clustered"];
  barChart("chart-homology", ["5-class test set", "Binary test set"], [
    { label: "Homology-based inference", backgroundColor: "#A8A8A8", data: tasks.map((t) => hom[t].homology_accuracy) },
    { label: "ProtT5 (all proteins)", backgroundColor: "#1B4F72", data: tasks.map((t) => hom[t].prott5_accuracy) },
    { label: "ProtT5 (proteins without homolog)", backgroundColor: "#5C8AE0", data: tasks.map((t) => hom[t].prott5_accuracy_without_hit) },
  ], { yLabel: "Accuracy" });
}

async function init() {
  try {
    const [manifest, summary] = await Promise.all([
      loadJson("/predictors/manifest.json"),
      loadJson("/predictors/results_summary.json"),
    ]);
    renderMetrics(summary);
    renderModelCards(manifest);
    renderCharts(summary);
  } catch (error) {
    document.getElementById("model-cards").textContent = `Could not load predictor data: ${error.message}`;
  }
}

init();

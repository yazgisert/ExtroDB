// Zenodo record of the predictor FASTAs. A new Zenodo version gets a new ID, so this is the only line to change.
const ZENODO_RECORD_ID = "23298950";

// Colors match the thesis figures (scripts/plotting/style.py)
const EMBEDDER_COLORS = {
  "Amino acid composition": "#AEBBCB",
  "ESM2 8M":                "#8095AE",
  "ESM2 650M":              "#55708F",
  "ProtT5":                 "#2E4A6B",
};
const EMBEDDER_ORDER = ["Amino acid composition", "ESM2 8M", "ESM2 650M", "ProtT5"];

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
  // no-cache: the browser revalidates, so updated result files are never mixed with newer page code
  const response = await fetch(path, { cache: "no-cache" });
  if (!response.ok) throw new Error(`Failed to load ${path}: ${response.status}`);
  return response.json();
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
  const doi = `10.5281/zenodo.${ZENODO_RECORD_ID}`;
  document.getElementById("zenodo-status").innerHTML = ZENODO_RECORD_ID
    ? `<a href="https://doi.org/${doi}" target="_blank" rel="noopener">doi:${doi}</a>`
    : "Zenodo (upload pending)";
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
  const detectionColors = ["#AEBBCB", "#8095AE", "#2E4A6B"];
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
    { label: "Homology-based inference", backgroundColor: "#BDBDBD", data: tasks.map((t) => hom[t].homology_accuracy) },
    { label: "ProtT5 (all proteins)", backgroundColor: "#2E4A6B", data: tasks.map((t) => hom[t].prott5_accuracy) },
    { label: "ProtT5 (proteins without homolog)", backgroundColor: "#55708F", data: tasks.map((t) => hom[t].prott5_accuracy_without_hit) },
  ], { yLabel: "Accuracy" });
}

async function init() {
  try {
    const [manifest, summary] = await Promise.all([
      loadJson("/predictors/manifest.json"),
      loadJson("/predictors/results_summary.json"),
    ]);
    renderModelCards(manifest);
    renderCharts(summary);
  } catch (error) {
    document.getElementById("model-cards").textContent = `Could not load predictor data: ${error.message}`;
  }
}

init();

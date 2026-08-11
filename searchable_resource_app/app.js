const state = {
  rows: [],
  summaryRows: [],
  filteredRows: [],
  selectedKey: null,
};

const filterIds = [
  "search",
  "domain-filter",
  "prokaryote-family",
  "prokaryote-class",
  "eukaryote-phenotype",
  "min-protein-count",
  "min-ogt",
  "max-ogt",
  "min-nacl",
  "max-nacl",
  "min-ph",
  "max-ph",
];

const CANONICAL_BRANCHES_BY_FAMILY = {
  acidophile: ["acidophile_expansion", "moderate_acidophile_expansion", "acidotolerant_expansion", "eukaryotic_algae_exploratory"],
  alkaliphile: ["alkaliphile_core"],
  halophile: ["halophile_core", "eukaryotic_algae_exploratory"],
  psychrophile: ["psychrophile_core", "eukaryotic_algae_exploratory"],
  thermophile: ["thermophile_core", "eukaryotic_algae_exploratory"],
};

function isEukaryote(row) {
  return splitMultiValue(row.dataset_branch).some((b) => b.startsWith("eukaryotic_algae"));
}

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Failed to load ${url}: ${response.status}`);
  }
  return response.json();
}

function formatInteger(value) {
  if (value === "" || value == null) {
    return "—";
  }
  return Number(value).toLocaleString();
}

function formatDecimal(value, suffix = "") {
  if (value === "" || value == null) {
    return "—";
  }
  return `${Number(value).toLocaleString(undefined, { maximumFractionDigits: 3 })}${suffix}`;
}

function thermalClass(row) {
  if (row.psychrophile === "yes") {
    return "psychrophile";
  }
  if (row.hyperthermophile === "yes") {
    return "hyperthermophile";
  }
  if (row.thermophile === "yes") {
    return "thermophile";
  }
  return "";
}

function titleCaseLabel(value) {
  if (!value) {
    return "";
  }
  return value
    .split("_")
    .map((part) => part.charAt(0).toUpperCase() + part.slice(1))
    .join(" ");
}

function phenotypeFocusLabel(value) {
  if (!value) {
    return "";
  }
  const labels = {
    acidophile: "Acidophile",
    alkaliphile: "Alkaliphile",
    eukaryotic_algae: "Eukaryotic Algae",
    halophile: "Halophile",
    psychrophile: "Psychrophile",
    thermophile: "Thermophile",
  };
  return labels[value] || titleCaseLabel(value);
}

function splitMultiValue(value) {
  if (!value) {
    return [];
  }
  return value
    .split("|")
    .map((item) => item.trim())
    .filter(Boolean);
}

function splitTextList(value) {
  if (!value) {
    return [];
  }
  return value
    .split(" || ")
    .map((item) => item.trim())
    .filter(Boolean);
}

function containsMultiValue(value, expected) {
  return splitMultiValue(value).includes(expected);
}

function joinedDisplay(value) {
  const values = splitMultiValue(value);
  return values.length ? values.join(", ") : "—";
}

function canonicalBranchesForFamily(family) {
  return CANONICAL_BRANCHES_BY_FAMILY[family] || [];
}

function rowBelongsToCanonicalFamilyBranch(row, family) {
  const canonicalBranches = canonicalBranchesForFamily(family);
  if (!canonicalBranches.length) {
    return false;
  }
  const rowBranches = splitMultiValue(row.dataset_branch);
  return canonicalBranches.some((branch) => rowBranches.includes(branch));
}

function phenotypeLabelValues(row) {
  if (row.phenotype_labels) {
    return splitMultiValue(row.phenotype_labels);
  }

  const acidClasses = splitMultiValue(row.acidophile_class);
  if (acidClasses.includes("acidotolerant")) {
    return ["Acidotolerant"];
  }
  if (acidClasses.includes("moderate_acidophile")) {
    return ["Moderate Acidophile"];
  }
  if (row.acidophile === "yes") {
    const labels = ["Acidophile"];
    if (acidClasses.includes("hyperacidophile")) {
      labels.push("Hyperacidophile");
    }
    return labels;
  }
  if (row.acidophilic_alga === "yes") {
    return ["Acidophilic Alga"];
  }

  const alkaliClasses = splitMultiValue(row.alkaliphile_class);
  if (alkaliClasses.includes("alkalitolerant")) {
    return ["Alkalitolerant"];
  }
  if (row.alkaliphile === "yes") {
    if (alkaliClasses.includes("extreme_alkaliphile")) {
      return ["Extreme Alkaliphile"];
    }
    return ["Alkaliphile"];
  }

  const halophileClasses = splitMultiValue(row.halophile_class);
  if (row.halophile === "yes") {
    if (halophileClasses.length) {
      return halophileClasses.map((value) => `Halophile (${titleCaseLabel(value)})`);
    }
    return ["Halophile"];
  }
  if (row.halophilic_alga === "yes" || row.saline_adapted_alga === "yes") {
    const labels = [];
    if (row.saline_adapted_alga === "yes") {
      labels.push("Saline-Adapted Alga");
    }
    if (row.halophilic_alga === "yes") {
      labels.push("Halophilic Alga");
    }
    return labels;
  }

  const label = thermalClass(row);
  if (label === "psychrophile") {
    return ["Psychrophile"];
  }
  if (row.psychrophilic_alga === "yes") {
    return ["Psychrophilic Alga"];
  }
  if (label === "hyperthermophile") {
    return ["Hyperthermophile"];
  }
  if (label === "thermophile") {
    return ["Thermophile"];
  }
  if (
    row.high_temperature_alga === "yes" ||
    row.thermophilic_alga === "yes" ||
    row.thermoacidophilic_alga === "yes"
  ) {
    const labels = [];
    if (row.high_temperature_alga === "yes") {
      labels.push("High-Temperature Alga");
    }
    if (row.thermophilic_alga === "yes") {
      labels.push("Thermophilic Alga");
    }
    if (row.thermoacidophilic_alga === "yes") {
      labels.push("Thermoacidophilic Alga");
    }
    return labels;
  }
  return row.phenotype_focus ? [row.phenotype_focus] : ["—"];
}

function phenotypeLabel(row) {
  return phenotypeLabelValues(row).join(", ");
}

function phenotypeClassValue(row) {
  return phenotypeLabelValues(row);
}

function rowKey(row) {
  return row.proteome_id || `${row.dataset_branch}::${row.taxon_name}::${row.strain}`;
}

function relativeFileHref(relativePath) {
  if (!relativePath) {
    return "";
  }
  return `/files/${relativePath.split("/").map(encodeURIComponent).join("/")}`;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function controlValue(id) {
  const element = document.getElementById(id);
  return element ? element.value.trim() : "";
}

function populateSelect(id, values, formatter = (value) => value) {
  const select = document.getElementById(id);
  const defaultOption = select.querySelector("option");
  const currentValue = select.value;
  select.innerHTML = "";
  select.append(defaultOption);
  defaultOption.selected = currentValue === "";

  values.forEach((value) => {
    if (!value) {
      return;
    }
    const option = document.createElement("option");
    option.value = value;
    option.textContent = formatter(value);
    if (value === currentValue) {
      option.selected = true;
    }
    select.append(option);
  });
}

function updateFilterVisibility() {
  const domain = controlValue("domain-filter");
  document.querySelectorAll(".prokaryote-field").forEach((el) => {
    el.style.display = domain === "prokaryote" ? "" : "none";
    if (domain !== "prokaryote") el.querySelector("select").value = "";
  });
  document.querySelectorAll(".eukaryote-field").forEach((el) => {
    el.style.display = domain === "eukaryote" ? "" : "none";
    if (domain !== "eukaryote") el.querySelector("select").value = "";
  });
}

function applyFilters() {
  const query = controlValue("search").toLowerCase();
  const prokaryoteFamily = controlValue("prokaryote-family");
  const prokaryoteClass = controlValue("prokaryote-class");
  const eukaryotePhenotype = controlValue("eukaryote-phenotype");
  const minProteinCount = Number(controlValue("min-protein-count") || 0);
  const minOgt = Number(controlValue("min-ogt") || 0);
  const maxOgtRaw = controlValue("max-ogt");
  const maxOgt = maxOgtRaw ? Number(maxOgtRaw) : null;
  const minNacl = Number(controlValue("min-nacl") || 0);
  const maxNaclRaw = controlValue("max-nacl");
  const maxNacl = maxNaclRaw ? Number(maxNaclRaw) : null;
  const minPhRaw = controlValue("min-ph");
  const minPh = minPhRaw ? Number(minPhRaw) : null;
  const maxPhRaw = controlValue("max-ph");
  const maxPh = maxPhRaw ? Number(maxPhRaw) : null;

  const domainFilter = controlValue("domain-filter");

  state.filteredRows = state.rows.filter((row) => {
    if (domainFilter === "prokaryote" && isEukaryote(row)) return false;
    if (domainFilter === "eukaryote" && !isEukaryote(row)) return false;
    if (prokaryoteFamily) {
      if (isEukaryote(row)) return false;
      if (!containsMultiValue(row.phenotype_focus, prokaryoteFamily)) return false;
    }
    if (prokaryoteClass) {
      if (isEukaryote(row)) return false;
      if (!phenotypeClassValue(row).includes(prokaryoteClass)) return false;
    }
    if (eukaryotePhenotype) {
      if (!isEukaryote(row)) return false;
      if (!phenotypeClassValue(row).includes(eukaryotePhenotype)) return false;
    }
    if (minProteinCount && Number(row.protein_count || 0) < minProteinCount) {
      return false;
    }
    if (minOgt && Number(row.ogt_celsius || 0) < minOgt) {
      return false;
    }
    if (maxOgt != null && row.ogt_celsius !== "" && Number(row.ogt_celsius) > maxOgt) {
      return false;
    }
    if (minNacl && Number(row.nacl_optimum_M || 0) < minNacl) {
      return false;
    }
    if (maxNacl != null && row.nacl_optimum_M !== "" && Number(row.nacl_optimum_M) > maxNacl) {
      return false;
    }
    if (minPh != null && row.ph_optimum !== "" && Number(row.ph_optimum) < minPh) {
      return false;
    }
    if (maxPh != null && row.ph_optimum !== "" && Number(row.ph_optimum) > maxPh) {
      return false;
    }
    if (!query) {
      return true;
    }

    const haystack = [
      row.taxon_name,
      row.strain,
      row.dataset_branch,
      row.phenotype_focus,
      row.sequence_source_accession,
      row.retrieval_notes,
      row.notes,
      row.evidence_source,
      row.local_fasta_file,
    ]
      .join(" ")
      .toLowerCase();

    return haystack.includes(query);
  });

  if (!state.filteredRows.some((row) => rowKey(row) === state.selectedKey)) {
    state.selectedKey = state.filteredRows[0] ? rowKey(state.filteredRows[0]) : null;
  }
}

function renderMetrics() {
  const container = document.getElementById("metrics-grid");
  const template = document.getElementById("metric-card-template");
  const phenotypeFamily = controlValue("prokaryote-family");
  const downloaded = state.filteredRows.filter((row) => row.proteome_status === "downloaded");
  const proteinSum = downloaded.reduce((sum, row) => sum + Number(row.protein_count || 0), 0);
  const branchCount = new Set(state.filteredRows.flatMap((row) => splitMultiValue(row.dataset_branch))).size;
  const phenotypeCount = new Set(state.filteredRows.flatMap((row) => splitMultiValue(row.phenotype_focus))).size;
  const phenotypeClassCount = new Set(state.filteredRows.flatMap((row) => phenotypeClassValue(row))).size;

  const metrics = [
    {
      label: "Visible Taxa",
      value: formatInteger(state.filteredRows.length),
      note: "Matching your filters — all with downloaded proteomes",
    },
    {
      label: "Proteins in View",
      value: formatInteger(proteinSum),
      note: "Summed protein count across downloaded rows",
    },
    {
      label: phenotypeFamily
        ? `${phenotypeFocusLabel(phenotypeFamily)} Family Branch Proteins`
        : "Family Branch Proteins",
      value: phenotypeFamily
        ? formatInteger(
            downloaded
              .filter((row) => rowBelongsToCanonicalFamilyBranch(row, phenotypeFamily))
              .reduce((sum, row) => sum + Number(row.protein_count || 0), 0)
          )
        : "—",
      note: phenotypeFamily
        ? "Within current filters, restricted to curation branches that belong to the selected phenotype family."
        : "Select a phenotype family to compare overlap against its linked curation branches.",
    },
    {
      label: "Curation Branch Coverage",
      value: formatInteger(branchCount),
      note: `${phenotypeCount} phenotype famil${phenotypeCount === 1 ? "y" : "ies"}, ${phenotypeClassCount} class label(s) represented`,
    },
  ];

  container.innerHTML = "";
  metrics.forEach((metric) => {
    const node = template.content.firstElementChild.cloneNode(true);
    node.querySelector(".metric-label").textContent = metric.label;
    node.querySelector(".metric-value").textContent = metric.value;
    node.querySelector(".metric-note").textContent = metric.note;
    container.append(node);
  });
}

function renderSummary() {
  const container = document.getElementById("summary-table");
  const classContainer = document.getElementById("class-summary-table");
  const rows = state.summaryRows
    .filter((row) => row.proteome_status === "downloaded")
    .map((row) => {
      return `
        <tr>
          <td>${escapeHtml(row.dataset_branch)}</td>
          <td>${escapeHtml(row.row_count)}</td>
          <td>${formatInteger(row.protein_count_sum)}</td>
        </tr>
      `;
    })
    .join("");

  container.innerHTML = `
    <p class="table-section-label">Curation branches</p>
    <table class="summary-table">
      <thead>
        <tr>
          <th>Curation Branch</th>
          <th>Rows</th>
          <th>Protein Count Sum</th>
        </tr>
      </thead>
      <tbody>${rows}</tbody>
    </table>
  `;

  const grouped = new Map();
  state.filteredRows.forEach((row) => {
    phenotypeClassValue(row).forEach((label) => {
      if (!grouped.has(label)) {
        grouped.set(label, { rows: 0, downloaded: 0, proteins: 0 });
      }
      const bucket = grouped.get(label);
      bucket.rows += 1;
      if (row.proteome_status === "downloaded") {
        bucket.downloaded += 1;
        bucket.proteins += Number(row.protein_count || 0);
      }
    });
  });

  const classRows = [...grouped.entries()]
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([label, values]) => {
      return `
        <tr>
          <td>${escapeHtml(label)}</td>
          <td>${formatInteger(values.rows)}</td>
          <td>${formatInteger(values.downloaded)}</td>
          <td>${formatInteger(values.proteins)}</td>
        </tr>
      `;
    })
    .join("");

  classContainer.innerHTML = `
    <p class="table-section-label">Phenotype classes in current view</p>
    <table class="summary-table">
      <thead>
        <tr>
          <th>Phenotype Class</th>
          <th>Rows</th>
          <th>Downloaded Proteomes</th>
          <th>Protein Count Sum</th>
        </tr>
      </thead>
      <tbody>${classRows}</tbody>
    </table>
  `;
}

function statusClass(status) {
  if (status === "downloaded") {
    return "status-download";
  }
  return "status-missing";
}

function renderResults() {
  const body = document.getElementById("results-body");
  const resultsTitle = document.getElementById("results-title");
  const downloadAllBtn = document.getElementById("download-all");
  const n = state.filteredRows.length;
  resultsTitle.textContent = `${n} matching entr${n === 1 ? "y" : "ies"}`;
  if (downloadAllBtn) {
    downloadAllBtn.textContent = n ? `Download All FASTAs (${n})` : "Download All FASTAs";
    downloadAllBtn.disabled = n === 0;
  }

  if (!state.filteredRows.length) {
    body.innerHTML = `
      <tr>
        <td colspan="7">
          <div class="empty-state">No rows match the current filters.</div>
        </td>
      </tr>
    `;
    return;
  }

  body.innerHTML = state.filteredRows
    .map((row) => {
      const selectedClass = rowKey(row) === state.selectedKey ? "is-selected" : "";
      return `
        <tr class="${selectedClass}" data-row-key="${escapeHtml(rowKey(row))}">
          <td>
            <strong>${escapeHtml(row.taxon_name)}</strong><br />
            <span>${escapeHtml(row.strain || "—")}</span>
          </td>
          <td>${escapeHtml(joinedDisplay(row.dataset_branch))}</td>
          <td>${escapeHtml(phenotypeLabel(row))}</td>
          <td>${formatInteger(row.protein_count)}</td>
          <td style="white-space:nowrap">${formatDecimal(row.ogt_celsius, " C")}</td>
          <td style="white-space:nowrap">${formatDecimal(row.nacl_optimum_M, " M")}</td>
          <td>${formatDecimal(row.ph_optimum)}</td>
        </tr>
      `;
    })
    .join("");

  body.querySelectorAll("tr[data-row-key]").forEach((row) => {
    row.addEventListener("click", () => {
      state.selectedKey = row.dataset.rowKey;
      renderResults();
      renderDetail();
    });
  });
}

function detailItem(label, value) {
  return `
    <div class="detail-item">
      <span class="detail-item-label">${escapeHtml(label)}</span>
      <div>${value || "—"}</div>
    </div>
  `;
}

function renderDetail() {
  const container = document.getElementById("detail-card");
  const selectedRow = state.filteredRows.find((row) => rowKey(row) === state.selectedKey);

  if (!selectedRow) {
    container.className = "detail-empty";
    container.textContent = "Select a row to inspect evidence, sequence source, notes, and local file links.";
    return;
  }

  const links = [];
  splitTextList(selectedRow.evidence_source).forEach((value, index) => {
    links.push(
      `<a class="detail-link" href="${escapeHtml(value)}" target="_blank" rel="noreferrer">Open evidence source${index ? ` ${index + 1}` : ""}</a>`
    );
  });
  if (selectedRow.proteome_id) {
    links.push(
      `<a class="button" style="text-decoration:none" href="/api/download/${encodeURIComponent(selectedRow.proteome_id)}">Download FASTA</a>`
    );
  }

  container.className = "detail-card";
  container.innerHTML = `
    <h3 class="detail-title">${escapeHtml(selectedRow.taxon_name)}</h3>
    <p class="detail-subtitle">${escapeHtml(selectedRow.strain || "No strain recorded")} • Curation branch: ${escapeHtml(joinedDisplay(selectedRow.dataset_branch))}</p>
    <div class="detail-grid">
      ${detailItem("Phenotype", escapeHtml(phenotypeLabel(selectedRow)))}
      ${detailItem("Proteome ID", escapeHtml(selectedRow.proteome_id || "—"))}
      ${detailItem("Domain Scope", escapeHtml(joinedDisplay(selectedRow.domain_scope)))}
      ${detailItem("Curation Branches", escapeHtml(joinedDisplay(selectedRow.dataset_branch)))}
      ${detailItem("Identity Level", escapeHtml(joinedDisplay(selectedRow.identity_level)))}
      ${detailItem("Collection Match", escapeHtml(joinedDisplay(selectedRow.collection_match)))}
      ${detailItem("Collection Source", escapeHtml(joinedDisplay(selectedRow.collection_source)))}
      ${detailItem("Strain Aliases", escapeHtml(joinedDisplay(selectedRow.strain_aliases)))}
      ${detailItem("Protein Count", escapeHtml(formatInteger(selectedRow.protein_count)))}
      ${detailItem("Sequence Source", escapeHtml(`${selectedRow.sequence_source_db || "—"} ${selectedRow.sequence_source_accession || ""}`.trim()))}
      ${detailItem("Sequence Level", escapeHtml(selectedRow.sequence_level || "—"))}
      ${detailItem("OGT", escapeHtml(formatDecimal(selectedRow.ogt_celsius, " C")))}
      ${detailItem("NaCl Optimum", escapeHtml(formatDecimal(selectedRow.nacl_optimum_M, " M")))}
      ${detailItem("pH Optimum", escapeHtml(formatDecimal(selectedRow.ph_optimum)))}
      ${detailItem("Retrieval Notes", escapeHtml(selectedRow.retrieval_notes || "—"))}
      ${detailItem("All Local FASTA Paths", escapeHtml(joinedDisplay(selectedRow.local_fasta_paths)))}
      ${detailItem("Notes", escapeHtml(selectedRow.notes || "—"))}
    </div>
    <div class="detail-links">${links.join("") || '<span class="empty-state">No direct links available for this row.</span>'}</div>
  `;
}

function render() {
  applyFilters();
  renderMetrics();
  renderResults();
  renderDetail();
}

function resetFilters() {
  filterIds.forEach((id) => {
    const element = document.getElementById(id);
    if (!element) {
      return;
    }
    element.value = "";
  });
  render();
}

async function downloadAllFastas() {
  const ids = state.filteredRows.map((row) => row.proteome_id).filter(Boolean);
  if (!ids.length) return;

  const btn = document.getElementById("download-all");
  const originalText = btn.textContent;
  btn.disabled = true;
  btn.textContent = "Building zip…";

  try {
    const response = await fetch("/api/download-bulk", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ids }),
    });

    if (!response.ok) {
      const detail = await response.json().then((d) => d.detail).catch(() => response.statusText);
      alert(`Download failed: ${detail}`);
      return;
    }

    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "extremophile_fastas.zip";
    a.click();
    URL.revokeObjectURL(url);
  } catch (err) {
    alert(`Download error: ${err.message}`);
  } finally {
    btn.disabled = false;
    btn.textContent = originalText;
  }
}

function exportResults() {
  const columns = state.rows.length ? Object.keys(state.rows[0]) : [];
  const lines = [columns.join("\t")];
  state.filteredRows.forEach((row) => {
    lines.push(columns.map((column) => row[column] || "").join("\t"));
  });

  const blob = new Blob([`${lines.join("\n")}\n`], { type: "text/tab-separated-values;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "extremophile_filtered_results.tsv";
  link.click();
  URL.revokeObjectURL(url);
}

function attachEvents() {
  filterIds.forEach((id) => {
    document.getElementById(id).addEventListener("input", render);
    document.getElementById(id).addEventListener("change", render);
  });
  document.getElementById("domain-filter").addEventListener("change", updateFilterVisibility);
  document.getElementById("download-all").addEventListener("click", downloadAllFastas);
  document.getElementById("reset-filters").addEventListener("click", () => {
    resetFilters();
    updateFilterVisibility();
  });
  document.getElementById("export-results").addEventListener("click", exportResults);
}

async function init() {
  try {
    const indexPayload = await fetchJson("/api/index");

    // Only show organisms that are actually accessible (have a downloaded FASTA).
    // Missing and not_started entries are excluded from the table entirely.
    state.rows = indexPayload.rows.filter((row) => row.proteome_status === "downloaded");
    state.filteredRows = state.rows.slice();
    state.selectedKey = state.rows[0] ? rowKey(state.rows[0]) : null;

    const updatedElement = document.getElementById("updated-at");
    if (updatedElement) {
      updatedElement.textContent = new Date(indexPayload.meta.updated_at).toLocaleString();
    }

    const prokaryoteRows = state.rows.filter((row) => !isEukaryote(row));
    const eukaryoteRows  = state.rows.filter((row) => isEukaryote(row));

    populateSelect(
      "prokaryote-family",
      [...new Set(prokaryoteRows.flatMap((row) => splitMultiValue(row.phenotype_focus)))].sort(),
      phenotypeFocusLabel
    );
    populateSelect(
      "prokaryote-class",
      [...new Set(prokaryoteRows.flatMap((row) => phenotypeClassValue(row)))].sort()
    );
    populateSelect(
      "eukaryote-phenotype",
      [...new Set(eukaryoteRows.flatMap((row) => phenotypeClassValue(row)))].sort()
    );

    attachEvents();
    updateFilterVisibility();
    render();
  } catch (error) {
    document.getElementById("results-body").innerHTML = `
      <tr>
        <td colspan="8">
          <div class="empty-state">${escapeHtml(error.message)}</div>
        </td>
      </tr>
    `;
  }
}

init();

# Normalised confusion matrix of the 5-class model on the leakage-filtered algae.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT  = Path("results/training_results/eukaryote_fiveclass_filtered_result.json")
OUT     = Path("results/plots/eukaryote_5class_confusion_filtered.png")
OUT.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":  "sans-serif",
    "font.size":    12,
    "figure.dpi":   150,
})

with open(RESULT) as f:
    data = json.load(f)

classes = data["classes"]           # alphabetical
cm = np.array(data["eukaryote_eval"]["confusion_matrix_normalised"])

# Remove only the alkaliphile row (no algae samples); keep its column, since the model still predicts it
present = [i for i, c in enumerate(classes) if c != "alkaliphile"]
cm_filt = cm[present, :]
row_labels = [classes[i].capitalize() for i in present]
col_labels = [c.capitalize() for c in classes]

fig, ax = plt.subplots(figsize=(8, 5))

im = ax.imshow(cm_filt, cmap="Blues", vmin=0, vmax=1)
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Fraction of true class")

for i, row in enumerate(present):
    for j in range(len(col_labels)):
        val = cm_filt[i, j]
        color = "white" if val > 0.55 else "black"
        weight = "bold" if row == j else "normal"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                fontsize=10, color=color, fontweight=weight)

ax.set_xticks(range(len(col_labels)))
ax.set_yticks(range(len(row_labels)))
ax.set_xticklabels(col_labels, rotation=30, ha="right")
ax.set_yticklabels(row_labels)
ax.set_xlabel("Predicted class (trained on prokaryotes)")
ax.set_ylabel("True class (eukaryote algae label)")

ax.set_title(
    "5-Class Model on Eukaryote Algae — Normalised Confusion Matrix\n"
    f"Leakage-filtered algae (n = {data['eukaryote_eval']['n']:,}) · no alkaliphile algae, so no alkaliphile row",
    fontsize=11
)

plt.tight_layout()
plt.savefig(OUT, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}")

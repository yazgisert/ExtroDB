# generate_table_and_curves.py
#
# Generates two additional result figures:
#   plots/results_table.png      — summary table (one row per split)
#   plots/training_curves.png    — train vs val loss curves for all 4 splits
#
# Run:
#   conda activate biocentral
#   cd thesis_data_collection
#   python generate_table_and_curves.py

import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

os.makedirs("plots", exist_ok=True)

# ── Load results from JSON files ──────────────────────────────────────────────

splits = ['split1_30pct', 'split2_30pct', 'split1_20pct', 'split2_20pct']
results = {}
for s in splits:
    with open(f'training_results/{s}_result.json') as f:
        results[s] = json.load(f)

# ── Palette ───────────────────────────────────────────────────────────────────

NAVY   = '#0D1B2A'
BLUE   = '#1B4F72'
TEAL   = '#1A9485'
WHITE  = '#FFFFFF'
OFFWHT = '#F5F7FA'
LGRAY  = '#DDE3EC'
DGRAY  = '#5A6A7A'
MGRAY  = '#9AAABB'

# ══════════════════════════════════════════════════════════════════════════════
# Plot 1: Summary results table
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(11, 3.2))
fig.patch.set_facecolor(OFFWHT)
ax.set_facecolor(OFFWHT)
ax.axis('off')

fig.suptitle('Classification Results — Summary Table',
             fontsize=13, fontweight='bold', color=NAVY, x=0.02, ha='left', y=0.97)
fig.text(0.02, 0.88, 'ProtT5 embeddings + FNN · test set metrics',
         fontsize=8, color=DGRAY)

# Table data
col_labels = ['Split', 'Label set', 'Threshold', 'n_train', 'n_val', 'n_test',
              'Classes', 'Accuracy', 'MCC', 'Early stop\n(epoch)']

row_data = [
    ['Split 1', 'Coarse (5)',     '30%', '99,964', '14,289', '28,576', '5',  '0.728', '0.662', '27'],
    ['Split 1', 'Coarse (5)',     '20%', '68,557',  '9,791', '19,583', '5',  '0.716', '0.646', '27'],
    ['Split 2', 'Fine-grained (12)', '30%', '43,294', '6,182', '12,374', '12', '0.550', '0.511', '39'],
    ['Split 2', 'Fine-grained (12)', '20%', '28,586', '4,081',  '8,172', '12', '0.537', '0.497', '35'],
]

table = ax.table(
    cellText=row_data,
    colLabels=col_labels,
    loc='center',
    cellLoc='center',
)

table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2.0)

# Style header row
for j in range(len(col_labels)):
    cell = table[0, j]
    cell.set_facecolor(NAVY)
    cell.set_text_props(color=WHITE, fontweight='bold', fontsize=8.5)
    cell.set_edgecolor(NAVY)

# Style data rows with alternating background
row_colors = [WHITE, LGRAY, WHITE, LGRAY]
metric_cols = {7, 8}  # Accuracy and MCC columns — highlight in blue
for i, row in enumerate(row_data):
    for j in range(len(col_labels)):
        cell = table[i + 1, j]
        cell.set_facecolor(row_colors[i])
        cell.set_edgecolor(LGRAY)
        if j in metric_cols:
            cell.set_text_props(color=BLUE, fontweight='bold')
        else:
            cell.set_text_props(color=DGRAY)

fig.subplots_adjust(top=0.82, bottom=0.02, left=0.01, right=0.99)
fig.savefig('plots/results_table.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/results_table.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2: Training curves (2×2 grid, one per split)
# ══════════════════════════════════════════════════════════════════════════════

split_labels = {
    'split1_30pct': 'Split 1 — 30% threshold (5 classes)',
    'split2_30pct': 'Split 2 — 30% threshold (12 classes)',
    'split1_20pct': 'Split 1 — 20% threshold (5 classes)',
    'split2_20pct': 'Split 2 — 20% threshold (12 classes)',
}

fig, axes = plt.subplots(2, 2, figsize=(11, 7))
fig.patch.set_facecolor(OFFWHT)

fig.suptitle('Training & Validation Loss Curves',
             fontsize=13, fontweight='bold', color=NAVY, x=0.02, ha='left', y=0.99)
fig.text(0.02, 0.963, 'Adam optimizer · CrossEntropyLoss · early stopping (patience=10)',
         fontsize=8, color=DGRAY)

for ax, split in zip(axes.flat, splits):
    r = results[split]
    train_loss = r['training_loss']
    val_loss   = r['validation_loss']
    epochs     = list(range(1, len(train_loss) + 1))
    best_epoch = int(np.argmin(val_loss)) + 1

    ax.set_facecolor(OFFWHT)
    ax.plot(epochs, train_loss, color=BLUE, linewidth=1.8, label='Train loss')
    ax.plot(epochs, val_loss,   color=TEAL, linewidth=1.8, label='Val loss',   linestyle='--')
    ax.axvline(best_epoch, color=MGRAY, linewidth=1.0, linestyle=':', zorder=0)
    ax.text(best_epoch + 0.3, max(train_loss) * 0.97,
            f'best\nepoch {best_epoch}', fontsize=7, color=MGRAY, va='top')

    ax.set_title(split_labels[split], fontsize=9.5, fontweight='bold',
                 color=NAVY, pad=6, loc='left')
    ax.set_xlabel('Epoch', fontsize=8, color=DGRAY)
    ax.set_ylabel('Loss', fontsize=8, color=DGRAY)
    ax.tick_params(labelsize=8, colors=DGRAY, length=0)
    ax.yaxis.grid(True, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color(LGRAY)
    ax.legend(fontsize=8, frameon=False, labelcolor=DGRAY)

fig.subplots_adjust(top=0.92, bottom=0.08, left=0.08, right=0.97,
                    hspace=0.42, wspace=0.28)
fig.savefig('plots/training_curves.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/training_curves.png")

print("\nDone. Both files saved to thesis_data_collection/plots/")

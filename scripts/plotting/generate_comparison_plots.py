# Comparison plots for ProtT5, ESM-2 650M, ESM-2 8M and one-hot on the 5-class and 12-class splits.

import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from style import EMBEDDER_COLORS

os.makedirs("results/plots", exist_ok=True)


SPLIT_KEYS = ['split1_30pct', 'split1_20pct', 'split2_30pct', 'split2_20pct']
def load(suffix):
    results = {}
    for s in SPLIT_KEYS:
        with open(f'results/training_results/{s}{suffix}_result.json') as f:
            results[s] = json.load(f)
    return results

prot    = load('')
esm2    = load('_esm2')
esm2_8m = load('_esm2_8m')
onehot  = load('_onehot')


NAVY   = '#0D1B2A'
OFFWHT = '#FFFFFF'
LGRAY  = '#DDE3EC'
BGRAY  = '#EEF1F6'
DGRAY  = '#5A6A7A'

def style_ax_v(ax):
    ax.set_facecolor(OFFWHT)
    ax.yaxis.grid(True, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color(LGRAY)
    ax.tick_params(length=0)

def style_ax_h(ax):
    ax.set_facecolor(OFFWHT)
    ax.xaxis.grid(True, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    ax.spines['left'].set_color(LGRAY)
    ax.tick_params(length=0)

def ci_err(r, metric, kind='overall'):
    if kind == 'overall':
        val   = r['test_metrics'][metric]
        lower = r['test_metrics_bootstrap'][metric]['lower']
        upper = r['test_metrics_bootstrap'][metric]['upper']
    else:
        cls, met = kind
        val   = r['per_class_metrics'][cls][met]
        lower = r['per_class_metrics_bootstrap'][cls][met]['lower']
        upper = r['per_class_metrics_bootstrap'][cls][met]['upper']
    return val - lower, upper - val

err_kw  = dict(fmt='none', elinewidth=1.5, capsize=4, capthick=1.5, ecolor=DGRAY, zorder=5)
herr_kw = dict(fmt='none', elinewidth=1.4, capsize=4, capthick=1.4, ecolor=DGRAY, zorder=5)

MODELS = [('ProtT5', prot, EMBEDDER_COLORS['prott5']), ('ESM-2 650M', esm2, EMBEDDER_COLORS['esm2_650m']),
          ('ESM-2 8M', esm2_8m, EMBEDDER_COLORS['esm2_8m']), ('Composition', onehot, EMBEDDER_COLORS['onehot'])]
W = 0.19
OFFSETS = [(i - 1.5) * (W + 0.02) for i in range(len(MODELS))]


def grouped_bars(ax, x, values, errors, decimals, fontsize, horizontal=False, rotate=False):
    for (label, _, color), off, vals, err in zip(MODELS, OFFSETS, values, errors):
        if horizontal:
            pos = x - off
            ax.barh(pos, vals, W, label=label, color=color, zorder=3, linewidth=0)
            ax.errorbar(vals, pos, xerr=err, **herr_kw)
            for p, v, e in zip(pos, vals, err[1]):
                ax.text(v + e + 0.012, p, f'{v:.{decimals}f}', va='center', fontsize=fontsize, color=DGRAY, fontweight='bold')
        else:
            pos = x + off
            ax.bar(pos, vals, W, label=label, color=color, zorder=3, linewidth=0)
            ax.errorbar(pos, vals, yerr=err, **err_kw)
            for p, v, e in zip(pos, vals, err[1]):
                ax.text(p, v + e + 0.014, f'{v:.{decimals}f}', ha='center', va='bottom', fontsize=fontsize, color=DGRAY, fontweight='bold',
                        rotation=90 if rotate else 0)


# Plot 1: overview, accuracy and MCC

x_labels = ['Split 1\n30%', 'Split 1\n20%', 'Split 2\n30%', 'Split 2\n20%']
x   = np.arange(len(SPLIT_KEYS))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), sharey=False)
fig.patch.set_facecolor(OFFWHT)

for ax, metric, ylabel in [(ax1, 'accuracy', 'Accuracy'),
                            (ax2, 'matthews-corr-coeff', 'MCC')]:
    style_ax_v(ax)

    values = [[r[s]['test_metrics'][metric] for s in SPLIT_KEYS] for _, r, _ in MODELS]
    errors = [np.array([ci_err(r[s], metric) for s in SPLIT_KEYS]).T for _, r, _ in MODELS]
    grouped_bars(ax, x, values, errors, decimals=2, fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10, color=DGRAY)
    ax.set_ylim(0, 0.95)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=9)
    ax.set_ylabel(ylabel, fontsize=11, color=DGRAY, labelpad=10)

    ax.axvspan(-0.5, 1.5, color=BGRAY, zorder=0, linewidth=0)
    ax.text(0.5, 0.895, '5 classes', ha='center', va='top',
            transform=ax.get_xaxis_transform(), fontsize=8, color=DGRAY, style='italic')
    ax.text(2.5, 0.895, '12 classes', ha='center', va='top',
            transform=ax.get_xaxis_transform(), fontsize=8, color=DGRAY, style='italic')

handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=10, frameon=False,
           labelcolor=DGRAY, handlelength=1.4, bbox_to_anchor=(0.5, -0.04))
fig.suptitle('ProtT5 vs ESM-2 650M vs ESM-2 8M vs Composition — Overview Comparison',
             fontsize=15, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.88, bottom=0.13, left=0.06, right=0.98, wspace=0.28)
fig.savefig('results/plots/overview_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: results/plots/overview_comparison.png")


# Plot 2: per-class F1, 5 classes (split 1, 30%)

rp = prot['split1_30pct']
classes5 = sorted(rp['per_class_metrics'],
                  key=lambda c: rp['per_class_metrics'][c]['f1_score'])
labels5  = [c.replace('_', ' ').capitalize() for c in classes5]

fig, ax = plt.subplots(figsize=(11, 5.5))
fig.patch.set_facecolor(OFFWHT)
style_ax_v(ax)

x = np.arange(len(classes5))
values = [[r['split1_30pct']['per_class_metrics'][c]['f1_score'] for c in classes5] for _, r, _ in MODELS]
errors = [np.array([ci_err(r['split1_30pct'], None, (c, 'f1_score')) for c in classes5]).T for _, r, _ in MODELS]
grouped_bars(ax, x, values, errors, decimals=2, fontsize=7.5)

ax.set_xticks(x); ax.set_xticklabels(labels5, fontsize=11, color=DGRAY)
ax.set_ylim(0, 1.08); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('F1 Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper left', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('5-Class Per-class F1 — ProtT5 vs ESM-2 650M vs ESM-2 8M vs Composition (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.10, left=0.09, right=0.97)
fig.savefig('results/plots/split1_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: results/plots/split1_comparison.png")


# Plot 3: per-class F1, 12 classes, horizontal (split 2, 30%)

rp = prot['split2_30pct']
classes12 = sorted(rp['per_class_metrics'],
                   key=lambda c: rp['per_class_metrics'][c]['f1_score'])
labels12  = [c.replace('_', ' ').capitalize() for c in classes12]
values12 = [[r['split2_30pct']['per_class_metrics'][c]['f1_score'] for c in classes12] for _, r, _ in MODELS]
errors12 = [np.array([ci_err(r['split2_30pct'], None, (c, 'f1_score')) for c in classes12]).T for _, r, _ in MODELS]

fig, ax = plt.subplots(figsize=(9, 12))
fig.patch.set_facecolor(OFFWHT)
style_ax_h(ax)

y = np.arange(len(classes12))
grouped_bars(ax, y, values12, errors12, decimals=2, fontsize=7.5, horizontal=True)

ax.set_yticks(y); ax.set_yticklabels(labels12, fontsize=11, color=DGRAY)
ax.set_xlim(0, 1.18); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.xaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_xlabel('F1 Score', fontsize=11, color=DGRAY, labelpad=10)
ax.legend(fontsize=10, frameon=False, loc='lower right', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('12-Class Per-class F1 — ProtT5 vs ESM-2 650M vs ESM-2 8M vs Composition (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.995)
fig.subplots_adjust(top=0.96, bottom=0.07, left=0.20, right=0.97)
fig.savefig('results/plots/split2_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: results/plots/split2_comparison.png")


# Plot 4: per-class F1, 12 classes, vertical (split 2, 30%)

labels12_2line = [c.replace('_', ' ').capitalize().replace(' ', '\n', 1) for c in classes12]

fig, ax = plt.subplots(figsize=(16, 6))
fig.patch.set_facecolor(OFFWHT)
style_ax_v(ax)

x = np.arange(len(classes12))
grouped_bars(ax, x, values12, errors12, decimals=2, fontsize=7, rotate=True)

ax.set_xticks(x); ax.set_xticklabels(labels12_2line, fontsize=9, color=DGRAY)
ax.set_ylim(0, 1.10); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('F1 Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper left', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('12-Class Per-class F1 — ProtT5 vs ESM-2 650M vs ESM-2 8M vs Composition (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.15, left=0.06, right=0.98)
fig.savefig('results/plots/split2_comparison_vertical.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: results/plots/split2_comparison_vertical.png")


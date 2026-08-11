# generate_comparison_plots.py
#
# 3-way comparison plots: ProtT5 vs ESM-2 650M vs one-hot encoding.
# Old plots are preserved — these are saved with _comparison suffix.
#
# Run from inside the results/ folder:
#   cd thesis_data_collection/results
#   python generate_comparison_plots.py
#
# Output:
#   plots/overview_comparison.png        — Accuracy & MCC, all three models
#   plots/split1_comparison.png          — per-class F1, 5 classes (vertical)
#   plots/split2_comparison.png          — per-class F1, 12 classes (horizontal)
#   plots/split2_comparison_vertical.png — per-class F1, 12 classes (vertical)

import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("plots", exist_ok=True)

# ── Load results ──────────────────────────────────────────────────────────────

SPLIT_KEYS = ['split1_30pct', 'split1_20pct', 'split2_30pct', 'split2_20pct']
prot   = {}
esm2   = {}
onehot = {}
for s in SPLIT_KEYS:
    with open(f'training_results/{s}_result.json') as f:
        prot[s] = json.load(f)
    with open(f'training_results/{s}_esm2_result.json') as f:
        esm2[s] = json.load(f)
    with open(f'training_results/{s}_onehot_result.json') as f:
        onehot[s] = json.load(f)

# ── Palette ───────────────────────────────────────────────────────────────────

NAVY   = '#0D1B2A'
BLUE   = '#1B4F72'   # ProtT5
GREEN  = '#1A7A4A'   # ESM-2 650M
AMBER  = '#C07A15'   # One-hot
OFFWHT = '#F5F7FA'
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


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1: Overview — Accuracy & MCC, two side-by-side subplots
# ══════════════════════════════════════════════════════════════════════════════

x_labels = ['Split 1\n30%', 'Split 1\n20%', 'Split 2\n30%', 'Split 2\n20%']
x   = np.arange(len(SPLIT_KEYS))
w   = 0.22
off = w + 0.03

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), sharey=False)
fig.patch.set_facecolor(OFFWHT)

for ax, metric, ylabel in [(ax1, 'accuracy', 'Accuracy'),
                            (ax2, 'matthews-corr-coeff', 'MCC')]:
    style_ax_v(ax)

    p_vals = [prot[s]['test_metrics'][metric]   for s in SPLIT_KEYS]
    e_vals = [esm2[s]['test_metrics'][metric]   for s in SPLIT_KEYS]
    o_vals = [onehot[s]['test_metrics'][metric] for s in SPLIT_KEYS]
    p_err  = np.array([ci_err(prot[s],   metric) for s in SPLIT_KEYS]).T
    e_err  = np.array([ci_err(esm2[s],   metric) for s in SPLIT_KEYS]).T
    o_err  = np.array([ci_err(onehot[s], metric) for s in SPLIT_KEYS]).T

    ax.bar(x - off, p_vals, w, label='ProtT5',     color=BLUE,  zorder=3, linewidth=0)
    ax.errorbar(x - off, p_vals, yerr=p_err, **err_kw)
    ax.bar(x,       e_vals, w, label='ESM-2 650M', color=GREEN, zorder=3, linewidth=0)
    ax.errorbar(x,       e_vals, yerr=e_err, **err_kw)
    ax.bar(x + off, o_vals, w, label='One-hot',    color=AMBER, zorder=3, linewidth=0)
    ax.errorbar(x + off, o_vals, yerr=o_err, **err_kw)

    for i, (pv, ev, ov, pe, ee, oe) in enumerate(
            zip(p_vals, e_vals, o_vals, p_err[1], e_err[1], o_err[1])):
        ax.text(x[i]-off, pv+pe+0.018, f'{pv:.3f}',
                ha='center', va='bottom', fontsize=7.5, color=BLUE,  fontweight='bold')
        ax.text(x[i],     ev+ee+0.018, f'{ev:.3f}',
                ha='center', va='bottom', fontsize=7.5, color=GREEN, fontweight='bold')
        ax.text(x[i]+off, ov+oe+0.018, f'{ov:.3f}',
                ha='center', va='bottom', fontsize=7.5, color=AMBER, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10, color=DGRAY)
    ax.set_ylim(0, 0.95)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=9)
    ax.set_ylabel(ylabel, fontsize=11, color=DGRAY, labelpad=10)
    ax.legend(fontsize=9, frameon=False, loc='upper right', labelcolor=DGRAY, handlelength=1.4)

    ax.axvspan(-0.5, 1.5, color=BGRAY, zorder=0, linewidth=0)
    ax.text(0.5, 0.895, '5 classes', ha='center', va='top',
            transform=ax.get_xaxis_transform(), fontsize=8, color=DGRAY, style='italic')
    ax.text(2.5, 0.895, '12 classes', ha='center', va='top',
            transform=ax.get_xaxis_transform(), fontsize=8, color=DGRAY, style='italic')

fig.suptitle('ProtT5 vs ESM-2 650M vs One-hot — Overview Comparison',
             fontsize=15, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.88, bottom=0.13, left=0.06, right=0.98, wspace=0.28)
fig.savefig('plots/overview_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/overview_comparison.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2: Split 1 per-class F1 (5 classes, vertical, 30% threshold)
# ══════════════════════════════════════════════════════════════════════════════

rp = prot['split1_30pct']; re = esm2['split1_30pct']; ro = onehot['split1_30pct']
classes5 = sorted(rp['per_class_metrics'],
                  key=lambda c: rp['per_class_metrics'][c]['f1_score'])
labels5  = [c.replace('_', ' ').capitalize() for c in classes5]

p_f5   = [rp['per_class_metrics'][c]['f1_score'] for c in classes5]
e_f5   = [re['per_class_metrics'][c]['f1_score'] for c in classes5]
o_f5   = [ro['per_class_metrics'][c]['f1_score'] for c in classes5]
p_f5_e = np.array([ci_err(rp, None, (c, 'f1_score')) for c in classes5]).T
e_f5_e = np.array([ci_err(re, None, (c, 'f1_score')) for c in classes5]).T
o_f5_e = np.array([ci_err(ro, None, (c, 'f1_score')) for c in classes5]).T

fig, ax = plt.subplots(figsize=(10, 5.5))
fig.patch.set_facecolor(OFFWHT)
style_ax_v(ax)

x = np.arange(len(classes5)); w = 0.22; off = w + 0.03
ax.bar(x-off, p_f5, w, label='ProtT5',     color=BLUE,  zorder=3, linewidth=0)
ax.errorbar(x-off, p_f5, yerr=p_f5_e, **err_kw)
ax.bar(x,    e_f5, w, label='ESM-2 650M',  color=GREEN, zorder=3, linewidth=0)
ax.errorbar(x,    e_f5, yerr=e_f5_e, **err_kw)
ax.bar(x+off, o_f5, w, label='One-hot',   color=AMBER, zorder=3, linewidth=0)
ax.errorbar(x+off, o_f5, yerr=o_f5_e, **err_kw)

for i, (pv, ev, ov, pe, ee, oe) in enumerate(
        zip(p_f5, e_f5, o_f5, p_f5_e[1], e_f5_e[1], o_f5_e[1])):
    ax.text(x[i]-off, pv+pe+0.016, f'{pv:.2f}',
            ha='center', va='bottom', fontsize=8, color=BLUE,  fontweight='bold')
    ax.text(x[i],     ev+ee+0.016, f'{ev:.2f}',
            ha='center', va='bottom', fontsize=8, color=GREEN, fontweight='bold')
    ax.text(x[i]+off, ov+oe+0.016, f'{ov:.2f}',
            ha='center', va='bottom', fontsize=8, color=AMBER, fontweight='bold')

ax.set_xticks(x); ax.set_xticklabels(labels5, fontsize=11, color=DGRAY)
ax.set_ylim(0, 1.08); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('F1 Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper left', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('5-Class Per-class F1 — ProtT5 vs ESM-2 vs One-hot (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.10, left=0.09, right=0.97)
fig.savefig('plots/split1_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/split1_comparison.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3: Split 2 per-class F1 (12 classes, horizontal, 30% threshold)
# ══════════════════════════════════════════════════════════════════════════════

rp = prot['split2_30pct']; re = esm2['split2_30pct']; ro = onehot['split2_30pct']
classes12 = sorted(rp['per_class_metrics'],
                   key=lambda c: rp['per_class_metrics'][c]['f1_score'])
labels12  = [c.replace('_', ' ').capitalize() for c in classes12]

p_f12   = [rp['per_class_metrics'][c]['f1_score'] for c in classes12]
e_f12   = [re['per_class_metrics'][c]['f1_score'] for c in classes12]
o_f12   = [ro['per_class_metrics'][c]['f1_score'] for c in classes12]
p_f12_e = np.array([ci_err(rp, None, (c, 'f1_score')) for c in classes12]).T
e_f12_e = np.array([ci_err(re, None, (c, 'f1_score')) for c in classes12]).T
o_f12_e = np.array([ci_err(ro, None, (c, 'f1_score')) for c in classes12]).T

fig, ax = plt.subplots(figsize=(9, 10))
fig.patch.set_facecolor(OFFWHT)
style_ax_h(ax)

y = np.arange(len(classes12)); h = 0.22; off = h + 0.03
ax.barh(y+off, p_f12, h, label='ProtT5',     color=BLUE,  zorder=3, linewidth=0)
ax.errorbar(p_f12, y+off, xerr=p_f12_e, **herr_kw)
ax.barh(y,    e_f12, h, label='ESM-2 650M',  color=GREEN, zorder=3, linewidth=0)
ax.errorbar(e_f12, y,    xerr=e_f12_e, **herr_kw)
ax.barh(y-off, o_f12, h, label='One-hot',   color=AMBER, zorder=3, linewidth=0)
ax.errorbar(o_f12, y-off, xerr=o_f12_e, **herr_kw)

for i, (pv, ev, ov, pe, ee, oe) in enumerate(
        zip(p_f12, e_f12, o_f12, p_f12_e[1], e_f12_e[1], o_f12_e[1])):
    ax.text(pv+pe+0.016, i+off, f'{pv:.2f}',
            va='center', fontsize=8.5, color=BLUE,  fontweight='bold')
    ax.text(ev+ee+0.016, i,     f'{ev:.2f}',
            va='center', fontsize=8.5, color=GREEN, fontweight='bold')
    ax.text(ov+oe+0.016, i-off, f'{ov:.2f}',
            va='center', fontsize=8.5, color=AMBER, fontweight='bold')

ax.set_yticks(y); ax.set_yticklabels(labels12, fontsize=11, color=DGRAY)
ax.set_xlim(0, 1.18); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.xaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_xlabel('F1 Score', fontsize=11, color=DGRAY, labelpad=10)
ax.legend(fontsize=10, frameon=False, loc='lower right', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('12-Class Per-class F1 — ProtT5 vs ESM-2 vs One-hot (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.995)
fig.subplots_adjust(top=0.96, bottom=0.07, left=0.20, right=0.97)
fig.savefig('plots/split2_comparison.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/split2_comparison.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4: Split 2 per-class F1 (12 classes, vertical, 30% threshold)
# ══════════════════════════════════════════════════════════════════════════════

labels12_2line = [c.replace('_', ' ').capitalize().replace(' ', '\n', 1) for c in classes12]

fig, ax = plt.subplots(figsize=(16, 6))
fig.patch.set_facecolor(OFFWHT)
style_ax_v(ax)

x = np.arange(len(classes12)); w = 0.22; off = w + 0.03
ax.bar(x-off, p_f12, w, label='ProtT5',     color=BLUE,  zorder=3, linewidth=0)
ax.errorbar(x-off, p_f12, yerr=p_f12_e, **err_kw)
ax.bar(x,    e_f12, w, label='ESM-2 650M',  color=GREEN, zorder=3, linewidth=0)
ax.errorbar(x,    e_f12, yerr=e_f12_e, **err_kw)
ax.bar(x+off, o_f12, w, label='One-hot',   color=AMBER, zorder=3, linewidth=0)
ax.errorbar(x+off, o_f12, yerr=o_f12_e, **err_kw)

for i, (pv, ev, ov, pe, ee, oe) in enumerate(
        zip(p_f12, e_f12, o_f12, p_f12_e[1], e_f12_e[1], o_f12_e[1])):
    ax.text(x[i]-off, pv+pe+0.016, f'{pv:.2f}',
            ha='center', va='bottom', fontsize=7.5, color=BLUE,  fontweight='bold')
    ax.text(x[i],     ev+ee+0.016, f'{ev:.2f}',
            ha='center', va='bottom', fontsize=7.5, color=GREEN, fontweight='bold')
    ax.text(x[i]+off, ov+oe+0.016, f'{ov:.2f}',
            ha='center', va='bottom', fontsize=7.5, color=AMBER, fontweight='bold')

ax.set_xticks(x); ax.set_xticklabels(labels12_2line, fontsize=9, color=DGRAY)
ax.set_ylim(0, 1.10); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('F1 Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper left', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('12-Class Per-class F1 — ProtT5 vs ESM-2 vs One-hot (30% threshold)',
             fontsize=14, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.15, left=0.06, right=0.98)
fig.savefig('plots/split2_comparison_vertical.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/split2_comparison_vertical.png")

print("\nAll comparison plots saved to results/plots/")

# generate_esm2_plots.py
#
# Standalone plots for ESM-2 650M results — same format as generate_plots.py.
# Run from inside the results/ folder:
#   cd thesis_data_collection/results
#   python generate_esm2_plots.py
#
# Output:
#   plots/esm2_overview_metrics.png
#   plots/esm2_split1_per_class_f1.png
#   plots/esm2_split2_per_class_f1.png

import json, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("plots", exist_ok=True)

SPLIT_KEYS = ['split1_30pct', 'split1_20pct', 'split2_30pct', 'split2_20pct']
results = {}
for s in SPLIT_KEYS:
    with open(f'training_results/{s}_esm2_result.json') as f:
        results[s] = json.load(f)

NAVY='#0D1B2A'; BLUE='#1B4F72'; TEAL='#1A9485'
OFFWHT='#F5F7FA'; LGRAY='#DDE3EC'; BGRAY='#EEF1F6'; DGRAY='#5A6A7A'

def style_ax_v(ax):
    ax.set_facecolor(OFFWHT)
    ax.yaxis.grid(True, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top','right','left']].set_visible(False)
    ax.spines['bottom'].set_color(LGRAY)
    ax.tick_params(length=0)

def style_ax_h(ax):
    ax.set_facecolor(OFFWHT)
    ax.xaxis.grid(True, color=LGRAY, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top','right','bottom']].set_visible(False)
    ax.spines['left'].set_color(LGRAY)
    ax.tick_params(length=0)

def overall_ci_err(r, metric):
    val   = r['test_metrics'][metric]
    lower = r['test_metrics_bootstrap'][metric]['lower']
    upper = r['test_metrics_bootstrap'][metric]['upper']
    return val - lower, upper - val

def per_class_ci_err(r, c, metric='f1_score'):
    val   = r['per_class_metrics'][c][metric]
    lower = r['per_class_metrics_bootstrap'][c][metric]['lower']
    upper = r['per_class_metrics_bootstrap'][c][metric]['upper']
    return val - lower, upper - val

err_kw  = dict(fmt='none', elinewidth=1.6, capsize=6,  capthick=1.6, ecolor=DGRAY, zorder=5)
herr_kw = dict(fmt='none', elinewidth=1.4, capsize=5,  capthick=1.4, ecolor=DGRAY, zorder=5)


# ── Plot 1: Overview ──────────────────────────────────────────────────────────

x_labels = ['Split 1\n30% threshold', 'Split 1\n20% threshold',
            'Split 2\n30% threshold', 'Split 2\n20% threshold']
acc = [results[s]['test_metrics']['accuracy']            for s in SPLIT_KEYS]
mcc = [results[s]['test_metrics']['matthews-corr-coeff'] for s in SPLIT_KEYS]
acc_err = np.array([overall_ci_err(results[s], 'accuracy')            for s in SPLIT_KEYS]).T
mcc_err = np.array([overall_ci_err(results[s], 'matthews-corr-coeff') for s in SPLIT_KEYS]).T

fig, ax = plt.subplots(figsize=(11, 5.8))
fig.patch.set_facecolor(OFFWHT)
style_ax_v(ax)
x = np.arange(len(x_labels)); width = 0.30

ax.axvspan(-0.5, 1.5, color=BGRAY, zorder=0, linewidth=0)
ax.text(0.5, 0.865, 'Split 1 · 5 classes', ha='center', va='top',
        transform=ax.get_xaxis_transform(), fontsize=8.5, color=DGRAY, style='italic')
ax.text(2.5, 0.865, 'Split 2 · 12 classes', ha='center', va='top',
        transform=ax.get_xaxis_transform(), fontsize=8.5, color=DGRAY, style='italic')

ax.bar(x-width/2, acc, width, label='Accuracy', color=BLUE, zorder=3, linewidth=0)
ax.errorbar(x-width/2, acc, yerr=acc_err, **err_kw)
ax.bar(x+width/2, mcc, width, label='MCC',      color=TEAL, zorder=3, linewidth=0)
ax.errorbar(x+width/2, mcc, yerr=mcc_err, **err_kw)

for i, (a, m) in enumerate(zip(acc, mcc)):
    ax.text(x[i]-width/2, a+acc_err[1][i]+0.022, f'{a:.3f}',
            ha='center', va='bottom', fontsize=8.5, color=BLUE, fontweight='bold')
    ax.text(x[i]+width/2, m+mcc_err[1][i]+0.022, f'{m:.3f}',
            ha='center', va='bottom', fontsize=8.5, color=TEAL, fontweight='bold')

ax.set_xticks(x); ax.set_xticklabels(x_labels, fontsize=11, color=DGRAY)
ax.set_ylim(0, 0.92); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper right', labelcolor=DGRAY,
          handlelength=1.4, handleheight=1.1)
fig.suptitle('Classification Results — ESM-2 650M — Overview',
             fontsize=16, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.13, left=0.08, right=0.97)
fig.savefig('plots/esm2_overview_metrics.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/esm2_overview_metrics.png")


# ── Plot 2: Split 1 per-class F1 (5 classes, vertical) ───────────────────────

r30 = results['split1_30pct']; r20 = results['split1_20pct']
classes5  = sorted(r30['per_class_metrics'], key=lambda c: r30['per_class_metrics'][c]['f1_score'])
f5_30     = [r30['per_class_metrics'][c]['f1_score'] for c in classes5]
f5_20     = [r20['per_class_metrics'][c]['f1_score'] for c in classes5]
labels5   = [c.replace('_', ' ').capitalize() for c in classes5]
f5_30_err = np.array([per_class_ci_err(r30, c) for c in classes5]).T
f5_20_err = np.array([per_class_ci_err(r20, c) for c in classes5]).T

fig, ax = plt.subplots(figsize=(9, 5.5))
fig.patch.set_facecolor(OFFWHT); style_ax_v(ax)
x = np.arange(len(classes5)); width = 0.30

ax.bar(x-width/2, f5_30, width, label='30% threshold', color=BLUE, zorder=3, linewidth=0)
ax.errorbar(x-width/2, f5_30, yerr=f5_30_err, **err_kw)
ax.bar(x+width/2, f5_20, width, label='20% threshold', color=TEAL, zorder=3, linewidth=0)
ax.errorbar(x+width/2, f5_20, yerr=f5_20_err, **err_kw)

for i, (v30, v20, e30, e20) in enumerate(zip(f5_30, f5_20, f5_30_err[1], f5_20_err[1])):
    ax.text(x[i]-width/2, v30+e30+0.018, f'{v30:.2f}',
            ha='center', va='bottom', fontsize=8.5, color=BLUE, fontweight='bold')
    ax.text(x[i]+width/2, v20+e20+0.018, f'{v20:.2f}',
            ha='center', va='bottom', fontsize=8.5, color=TEAL, fontweight='bold')

ax.set_xticks(x); ax.set_xticklabels(labels5, fontsize=11, color=DGRAY)
ax.set_ylim(0, 1.05); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.yaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_ylabel('F1 Score', fontsize=11, color=DGRAY, labelpad=12)
ax.legend(fontsize=10, frameon=False, loc='upper left', labelcolor=DGRAY,
          handlelength=1.4, handleheight=1.1)
fig.suptitle('5-Class Classification — ESM-2 650M — Per-class F1 Score',
             fontsize=15, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.97)
fig.subplots_adjust(top=0.90, bottom=0.10, left=0.09, right=0.97)
fig.savefig('plots/esm2_split1_per_class_f1.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/esm2_split1_per_class_f1.png")


# ── Plot 3: Split 2 per-class F1 (12 classes, horizontal) ────────────────────

r30 = results['split2_30pct']; r20 = results['split2_20pct']
classes12  = sorted(r30['per_class_metrics'], key=lambda c: r30['per_class_metrics'][c]['f1_score'])
f12_30     = [r30['per_class_metrics'][c]['f1_score'] for c in classes12]
f12_20     = [r20['per_class_metrics'][c]['f1_score'] for c in classes12]
labels12   = [c.replace('_', ' ').capitalize() for c in classes12]
f12_30_err = np.array([per_class_ci_err(r30, c) for c in classes12]).T
f12_20_err = np.array([per_class_ci_err(r20, c) for c in classes12]).T

fig, ax = plt.subplots(figsize=(9, 8.5))
fig.patch.set_facecolor(OFFWHT); style_ax_h(ax)
y = np.arange(len(classes12)); h = 0.30

ax.barh(y+h/2, f12_30, h, label='30% threshold', color=BLUE, zorder=3, linewidth=0)
ax.errorbar(f12_30, y+h/2, xerr=f12_30_err, **herr_kw)
ax.barh(y-h/2, f12_20, h, label='20% threshold', color=TEAL, zorder=3, linewidth=0)
ax.errorbar(f12_20, y-h/2, xerr=f12_20_err, **herr_kw)

for i, (v30, v20, e30, e20) in enumerate(zip(f12_30, f12_20, f12_30_err[1], f12_20_err[1])):
    ax.text(v30+e30+0.018, i+h/2, f'{v30:.2f}',
            va='center', fontsize=9, color=BLUE, fontweight='bold')
    ax.text(v20+e20+0.018, i-h/2, f'{v20:.2f}',
            va='center', fontsize=9, color=TEAL, fontweight='bold')

ax.set_yticks(y); ax.set_yticklabels(labels12, fontsize=11, color=DGRAY)
ax.set_xlim(0, 1.15); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.xaxis.set_tick_params(labelcolor=DGRAY, labelsize=10)
ax.set_xlabel('F1 Score', fontsize=11, color=DGRAY, labelpad=10)
ax.legend(fontsize=10, frameon=False, loc='lower right', labelcolor=DGRAY, handlelength=1.4)
fig.suptitle('12-Class Classification — ESM-2 650M — Per-class F1 Score',
             fontsize=15, fontweight='bold', color=NAVY, x=0.025, ha='left', y=0.99)
fig.subplots_adjust(top=0.95, bottom=0.07, left=0.20, right=0.97)
fig.savefig('plots/esm2_split2_per_class_f1.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: plots/esm2_split2_per_class_f1.png")

print("\nAll ESM-2 plots saved to results/plots/")

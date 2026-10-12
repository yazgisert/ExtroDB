# Shared plot colours: one fixed colour per extremophile class, and a blue-grey ramp for comparing models or groups
# (light = baseline / simpler, dark = main model), so class colours and model colours never clash.

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

CLASS_COLORS = {
    "thermophile":  "#E05C5C",
    "halophile":    "#E0A040",
    "psychrophile": "#5C8AE0",
    "acidophile":   "#6ABF69",
    "alkaliphile":  "#A85CE0",
}
# Neutral greys: non-extremophiles, and external baselines or aggregates (homology, THERMO, accuracy bars)
NON_EXTREMOPHILE = "#8C96A3"
BASELINE = "#BDBDBD"

# Blue-grey ramp, light to dark
RAMP = ["#AEBBCB", "#8095AE", "#55708F", "#2E4A6B"]
MAIN, SECONDARY = RAMP[3], RAMP[1]

EMBEDDER_COLORS = {
    "onehot":    RAMP[0],
    "esm2_8m":   RAMP[1],
    "esm2_650m": RAMP[2],
    "prott5":    RAMP[3],
}


def light(color, amount=0.5):
    # Tint of a colour (mixed with white), for the second bar of a pair
    r, g, b = mcolors.to_rgb(color)
    return mcolors.to_hex((r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount))


def apply():
    plt.rcParams.update({
        "font.family":       "sans-serif",
        "font.size":         12,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.facecolor":    "white",
        "figure.facecolor":  "white",
        "figure.dpi":        150,
    })

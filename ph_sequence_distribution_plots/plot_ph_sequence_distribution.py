from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "thesis_data_collection"
OUTPUT_DIR = Path(__file__).resolve().parent

# These are the current manifest files we want to visualize.
# The acid-side plot combines strict acidophiles, moderate acidophiles,
# and acidotolerant entries. The alkaline-side plot uses the full current
# alkaliphile manifest, which already includes alkalitolerant entries.
PLOT_CONFIGS = {
    "acid_side": {
        "title": "Acid-Side pH Distribution by Sequence Count",
        "subtitle": "Combined acidophile, moderate-acidophile, and acidotolerant manifests",
        "manifests": [
            DATA_DIR
            / "acidophile_manifests"
            / "sequences_manifest_combined_acidophiles_36.tsv",
            DATA_DIR
            / "acidotolerant_manifests"
            / "sequences_manifest_acidotolerant_10.tsv",
        ],
        "x_min": 0.0,
        "x_max": 7.5,
        "thresholds": [3.0, 5.0, 7.5],
        "intervals": [
            (0.0, 3.0, "<= 3.0"),
            (3.0, 5.0, "3.0-5.0"),
            (5.0, 7.5, "5.0-7.5"),
        ],
        "bar_color": "#c75b39",
        "shade_colors": ["#fde8e2", "#fff4df", "#eef6e9"],
        "svg_name": "acid_side_sequence_distribution.svg",
        "exact_csv_name": "acid_side_sequence_counts_by_ph.csv",
        "interval_csv_name": "acid_side_sequence_counts_by_interval.csv",
    },
    "alkaline_side": {
        "title": "Alkaline-Side pH Distribution by Sequence Count",
        "subtitle": "Current alkaliphile manifest, including entries below pH 9.0",
        "manifests": [
            DATA_DIR
            / "alkaliphile_manifests"
            / "sequences_manifest_alkaliphile_16.tsv"
        ],
        "x_min": 7.5,
        "x_max": 10.75,
        "thresholds": [9.0, 10.0],
        "intervals": [
            (7.5, 9.0, "7.5-<9.0"),
            (9.0, 10.0, "9.0-<10.0"),
            (10.0, 10.75, ">= 10.0"),
        ],
        "bar_color": "#2d6d8a",
        "shade_colors": ["#e8f3f8", "#edf4ea", "#f7ecef"],
        "svg_name": "alkaline_side_sequence_distribution.svg",
        "exact_csv_name": "alkaline_side_sequence_counts_by_ph.csv",
        "interval_csv_name": "alkaline_side_sequence_counts_by_interval.csv",
    },
}


def load_manifest_rows(paths: Iterable[Path]) -> list[dict[str, str]]:
    rows_by_accession: dict[str, dict[str, str]] = {}

    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                accession = row["assembly_accession"].strip()
                if accession in rows_by_accession:
                    continue
                rows_by_accession[accession] = row

    rows = list(rows_by_accession.values())
    rows.sort(key=lambda row: (float(row["ph_optimum"]), row["assembly_accession"]))
    return rows


def aggregate_exact_counts(rows: Iterable[dict[str, str]]) -> list[tuple[float, int, int]]:
    counts: dict[float, int] = {}
    proteome_counts: dict[float, int] = {}

    for row in rows:
        ph_value = float(row["ph_optimum"])
        protein_count = int(row["protein_count"])
        counts[ph_value] = counts.get(ph_value, 0) + protein_count
        proteome_counts[ph_value] = proteome_counts.get(ph_value, 0) + 1

    return [
        (ph_value, counts[ph_value], proteome_counts[ph_value])
        for ph_value in sorted(counts)
    ]


def aggregate_interval_counts(
    rows: Iterable[dict[str, str]],
    intervals: list[tuple[float, float, str]],
) -> list[tuple[str, int, int]]:
    interval_totals = {label: 0 for _, _, label in intervals}
    interval_proteomes = {label: 0 for _, _, label in intervals}

    for row in rows:
        ph_value = float(row["ph_optimum"])
        protein_count = int(row["protein_count"])
        for start, end, label in intervals:
            # The final interval is inclusive on the right edge.
            in_interval = start <= ph_value < end or (
                label == intervals[-1][2] and start <= ph_value <= end
            )
            if in_interval:
                interval_totals[label] += protein_count
                interval_proteomes[label] += 1
                break

    return [
        (label, interval_totals[label], interval_proteomes[label])
        for _, _, label in intervals
    ]


def write_exact_csv(path: Path, rows: list[tuple[float, int, int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ph_optimum", "sequence_count", "proteome_count"])
        for ph_value, sequence_count, proteome_count in rows:
            writer.writerow([f"{ph_value:.3f}", sequence_count, proteome_count])


def write_interval_csv(path: Path, rows: list[tuple[str, int, int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["interval", "sequence_count", "proteome_count"])
        for label, sequence_count, proteome_count in rows:
            writer.writerow([label, sequence_count, proteome_count])


def nice_y_max(value: int) -> int:
    if value <= 10000:
        step = 1000
    elif value <= 50000:
        step = 5000
    else:
        step = 10000
    return ((value + step - 1) // step) * step


def clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(value, maximum))


def svg_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def build_svg(
    title: str,
    subtitle: str,
    exact_rows: list[tuple[float, int, int]],
    interval_rows: list[tuple[str, int, int]],
    intervals: list[tuple[float, float, str]],
    thresholds: list[float],
    x_min: float,
    x_max: float,
    bar_color: str,
    shade_colors: list[str],
) -> str:
    width = 1280
    height = 980
    margin_left = 100
    margin_right = 60
    margin_top = 118
    margin_bottom = 250
    plot_width = width - margin_left - margin_right
    plot_height = 470
    plot_left = margin_left
    plot_top = margin_top
    plot_bottom = plot_top + plot_height

    max_count = max(sequence_count for _, sequence_count, _ in exact_rows)
    y_max = nice_y_max(max_count)
    y_ticks = 6

    def x_to_px(x_value: float) -> float:
        return plot_left + ((x_value - x_min) / (x_max - x_min)) * plot_width

    def y_to_px(y_value: float) -> float:
        return plot_bottom - (y_value / y_max) * plot_height

    unique_ph_values = [ph for ph, _, _ in exact_rows]
    if len(unique_ph_values) > 1:
        min_gap = min(
            b - a for a, b in zip(unique_ph_values, unique_ph_values[1:]) if b > a
        )
    else:
        min_gap = 0.25
    half_bar_width = min_gap * 0.32
    font_family = "Arial, Helvetica, sans-serif"
    axis_color = "#203040"
    grid_color = "#d9dee5"
    text_muted = "#5d6b7a"
    summary_x = plot_left
    summary_y = plot_bottom + 150
    summary_width = plot_width
    summary_height = 162

    svg_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fcfcfb"/>',
        f'<rect x="{plot_left - 24}" y="24" width="{plot_width + 48}" height="{summary_y + summary_height - 4}" rx="24" ry="24" fill="#ffffff" stroke="#e6e9ee"/>',
        f'<text x="{plot_left}" y="54" font-family="{font_family}" font-size="29" font-weight="700" fill="#16212b">{svg_escape(title)}</text>',
        f'<text x="{plot_left}" y="84" font-family="{font_family}" font-size="16" fill="{text_muted}">{svg_escape(subtitle)}</text>',
    ]

    # Background shading for the chosen pH intervals.
    for index, (start, end, label) in enumerate(intervals):
        x_start = x_to_px(start)
        x_end = x_to_px(end)
        svg_parts.append(
            f'<rect x="{x_start:.2f}" y="{plot_top}" width="{(x_end - x_start):.2f}" height="{plot_height}" fill="{shade_colors[index % len(shade_colors)]}"/>'
        )
        pill_width = max(88, len(label) * 8 + 20)
        pill_x = clamp(
            (x_start + x_end - pill_width) / 2,
            plot_left + 6,
            plot_left + plot_width - pill_width - 6,
        )
        pill_y = plot_top - 34
        svg_parts.append(
            f'<rect x="{pill_x:.2f}" y="{pill_y:.2f}" width="{pill_width:.2f}" height="24" rx="12" ry="12" fill="#ffffff" stroke="#d8dee6"/>'
        )
        svg_parts.append(
            f'<text x="{(pill_x + pill_width / 2):.2f}" y="{pill_y + 16:.2f}" text-anchor="middle" font-family="{font_family}" font-size="13" font-weight="700" fill="#415466">{svg_escape(label)}</text>'
        )

    # Horizontal gridlines and y tick labels.
    for tick_index in range(y_ticks + 1):
        y_value = (y_max / y_ticks) * tick_index
        y_px = y_to_px(y_value)
        svg_parts.append(
            f'<line x1="{plot_left}" y1="{y_px:.2f}" x2="{plot_left + plot_width}" y2="{y_px:.2f}" stroke="{grid_color}" stroke-width="1"/>'
        )
        svg_parts.append(
            f'<text x="{plot_left - 12}" y="{y_px + 5:.2f}" text-anchor="end" font-family="{font_family}" font-size="13" fill="{text_muted}">{int(y_value):,}</text>'
        )

    # Threshold lines.
    for threshold in thresholds:
        if x_min <= threshold <= x_max:
            x_px = x_to_px(threshold)
            label = f"pH {threshold:g}"
            label_width = max(58, len(label) * 8 + 18)
            label_x = clamp(
                x_px - label_width / 2,
                plot_left + 6,
                plot_left + plot_width - label_width - 6,
            )
            svg_parts.append(
                f'<line x1="{x_px:.2f}" y1="{plot_top}" x2="{x_px:.2f}" y2="{plot_bottom}" stroke="#49606f" stroke-width="2" stroke-dasharray="7,6"/>'
            )
            svg_parts.append(
                f'<rect x="{label_x:.2f}" y="{plot_top + 10:.2f}" width="{label_width:.2f}" height="22" rx="11" ry="11" fill="#ffffff" stroke="#cad2dc"/>'
            )
            svg_parts.append(
                f'<text x="{(label_x + label_width / 2):.2f}" y="{plot_top + 25:.2f}" text-anchor="middle" font-family="{font_family}" font-size="12" font-weight="700" fill="#49606f">{label}</text>'
            )

    # Bars for exact pH optima.
    for index, (ph_value, sequence_count, _) in enumerate(exact_rows):
        x_left = x_to_px(ph_value - half_bar_width)
        x_right = x_to_px(ph_value + half_bar_width)
        bar_top = y_to_px(sequence_count)
        bar_height = plot_bottom - bar_top
        svg_parts.append(
            f'<rect x="{x_left:.2f}" y="{bar_top:.2f}" width="{(x_right - x_left):.2f}" height="{bar_height:.2f}" rx="3" ry="3" fill="{bar_color}" stroke="#1e2b36" stroke-opacity="0.25" opacity="0.94"/>'
        )
        svg_parts.append(
            f'<line x1="{((x_left + x_right) / 2):.2f}" y1="{plot_bottom:.2f}" x2="{((x_left + x_right) / 2):.2f}" y2="{plot_bottom + 8:.2f}" stroke="#7e8d9a" stroke-width="1"/>'
        )
        label_y = plot_bottom + 34 + (index % 2) * 16
        label_x = (x_left + x_right) / 2
        svg_parts.append(
            f'<text x="{label_x:.2f}" y="{label_y:.2f}" transform="rotate(-45 {label_x:.2f} {label_y:.2f})" text-anchor="end" font-family="{font_family}" font-size="12" fill="#334455">{ph_value:.2f}</text>'
        )

    # Axes.
    svg_parts.append(
        f'<line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_left + plot_width}" y2="{plot_bottom}" stroke="{axis_color}" stroke-width="2.2"/>'
    )
    svg_parts.append(
        f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" y2="{plot_bottom}" stroke="{axis_color}" stroke-width="2.2"/>'
    )
    svg_parts.append(
        f'<text x="{plot_left + plot_width / 2:.2f}" y="{plot_bottom + 118:.2f}" text-anchor="middle" font-family="{font_family}" font-size="16" font-weight="700" fill="{axis_color}">Stored pH optimum in the manifest</text>'
    )
    svg_parts.append(
        f'<text x="28" y="{plot_top + plot_height / 2:.2f}" transform="rotate(-90 28 {plot_top + plot_height / 2:.2f})" text-anchor="middle" font-family="{font_family}" font-size="16" font-weight="700" fill="{axis_color}">Number of protein sequences</text>'
    )

    total_sequences = sum(sequence_count for _, sequence_count, _ in exact_rows)
    total_proteomes = sum(proteome_count for _, _, proteome_count in exact_rows)
    svg_parts.append(
        f'<rect x="{summary_x}" y="{summary_y}" width="{summary_width}" height="{summary_height}" rx="16" ry="16" fill="#f7f9fb" stroke="#d8dee7"/>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 20}" y="{summary_y + 30}" font-family="{font_family}" font-size="16" font-weight="700" fill="#16212b">Summary</text>'
    )
    svg_parts.append(
        f'<rect x="{summary_x + 20}" y="{summary_y + 44}" width="180" height="42" rx="10" ry="10" fill="#ffffff" stroke="#d5dce5"/>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 34}" y="{summary_y + 62}" font-family="{font_family}" font-size="12" fill="{text_muted}">Total sequences</text>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 34}" y="{summary_y + 80}" font-family="{font_family}" font-size="19" font-weight="700" fill="#16212b">{total_sequences:,}</text>'
    )
    svg_parts.append(
        f'<rect x="{summary_x + 216}" y="{summary_y + 44}" width="180" height="42" rx="10" ry="10" fill="#ffffff" stroke="#d5dce5"/>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 230}" y="{summary_y + 62}" font-family="{font_family}" font-size="12" fill="{text_muted}">Total proteomes</text>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 230}" y="{summary_y + 80}" font-family="{font_family}" font-size="19" font-weight="700" fill="#16212b">{total_proteomes}</text>'
    )
    svg_parts.append(
        f'<rect x="{summary_x + 412}" y="{summary_y + 44}" width="{summary_width - 432}" height="42" rx="10" ry="10" fill="#ffffff" stroke="#d5dce5"/>'
    )
    svg_parts.append(
        f'<text x="{summary_x + 428}" y="{summary_y + 69}" font-family="{font_family}" font-size="13" fill="{text_muted}">Sequence counts come from the manifest protein_count values grouped by stored pH optimum.</text>'
    )

    card_y = summary_y + 102
    card_gap = 16
    card_width = (summary_width - 40 - card_gap * (len(interval_rows) - 1)) / len(interval_rows)
    for index, (label, sequence_count, proteome_count) in enumerate(interval_rows):
        card_x = summary_x + 20 + index * (card_width + card_gap)
        svg_parts.append(
            f'<rect x="{card_x:.2f}" y="{card_y:.2f}" width="{card_width:.2f}" height="44" rx="10" ry="10" fill="#ffffff" stroke="#d5dce5"/>'
        )
        svg_parts.append(
            f'<text x="{card_x + 14:.2f}" y="{card_y + 18:.2f}" font-family="{font_family}" font-size="12" font-weight="700" fill="#415466">{svg_escape(label)}</text>'
        )
        svg_parts.append(
            f'<text x="{card_x + 14:.2f}" y="{card_y + 34:.2f}" font-family="{font_family}" font-size="13" fill="#253443">{sequence_count:,} sequences across {proteome_count} proteomes</text>'
        )

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for plot_name, config in PLOT_CONFIGS.items():
        rows = load_manifest_rows(config["manifests"])
        exact_rows = aggregate_exact_counts(rows)
        interval_rows = aggregate_interval_counts(rows, config["intervals"])

        write_exact_csv(OUTPUT_DIR / config["exact_csv_name"], exact_rows)
        write_interval_csv(OUTPUT_DIR / config["interval_csv_name"], interval_rows)

        svg_text = build_svg(
            title=config["title"],
            subtitle=config["subtitle"],
            exact_rows=exact_rows,
            interval_rows=interval_rows,
            intervals=config["intervals"],
            thresholds=config["thresholds"],
            x_min=config["x_min"],
            x_max=config["x_max"],
            bar_color=config["bar_color"],
            shade_colors=config["shade_colors"],
        )
        (OUTPUT_DIR / config["svg_name"]).write_text(svg_text, encoding="utf-8")
        print(f"Wrote {plot_name} outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

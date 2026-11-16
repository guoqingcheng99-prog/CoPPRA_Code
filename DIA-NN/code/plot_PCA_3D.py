"""Render saved CoPPRA PCA scores in the bioRxiv Figure 2A style.

The Figure 2A view emphasizes arithmetic group means, while showing every
biological replicate in gray. The companion view labels all replicates.
Coordinates, PC signs, sample membership, and PCA fitting are never changed.
Both views use the same orthographic camera and equal physical units on all axes.
Only the title, axes, legend, and (in the replicate view) sample labels are drawn.
Group-mean definitions and camera provenance are retained in the saved metadata.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

GROUPS = ["Ctrl", "R", "U", "RU"]
PCS = ["PC1", "PC2", "PC3"]
COLORS = dict(Ctrl="#1F77B4", R="#D62728", U="#2CA02C", RU="#FF7F0E")
MARKERS = dict(Ctrl="o", R="^", U="s", RU="D")
FIGURE_SIZE_INCHES = (4.5, 4.5)
FONT_SIZES = dict(title=18, axis_title=14, axis_tick=13,
                  legend_title=14, legend_text=13)
LEGEND_ANCHOR = (0.73, 0.65)


def read_inputs(scores_path, variance_path):
    with open(scores_path) as handle:
        rows = list(csv.DictReader(handle))
    with open(variance_path) as handle:
        variance = {row["PC"]: float(row["Variance_percent"])
                    for row in csv.DictReader(handle)}
    xyz = np.array([[float(row[pc]) for pc in PCS] for row in rows])
    groups = np.array([row["PCA_group"] for row in rows])
    if len(rows) != 12 or any(np.sum(groups == g) != 3 for g in GROUPS):
        raise ValueError("Expected 3 biological replicates per group (12 total).")
    if len({r["label"] for r in rows}) != len(rows) or not np.isfinite(xyz).all():
        raise ValueError("Sample labels must be unique and PCA scores finite.")
    if not np.isfinite([variance[p] for p in PCS]).all():
        raise ValueError("Missing or nonfinite explained variance.")
    return rows, xyz, groups, variance


def draw(rows, xyz, groups, variance, elev=30, azim=165,
         show_means=True):
    plt.rcParams.update({"font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 12,
        "axes.labelsize": FONT_SIZES["axis_title"], "savefig.facecolor": "white"})
    fig = plt.figure(figsize=FIGURE_SIZE_INCHES, facecolor="white")
    ax = fig.add_axes([-0.01, 0.09, 0.80, 0.78], projection="3d")
    ax.set_proj_type("ortho")
    ax.view_init(elev=elev, azim=azim)
    spans = np.ptp(xyz, axis=0)
    if np.any(spans == 0):
        raise ValueError("Each displayed PC must have a nonzero score range.")
    limits = np.column_stack((xyz.min(0) - 0.12 * spans,
                              xyz.max(0) + 0.12 * spans))
    ax.set_xlim(*limits[0]); ax.set_ylim(*limits[1]); ax.set_zlim(*limits[2])
    ax.set_box_aspect(np.ptp(limits, axis=1), zoom=0.95)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
        axis.set_pane_color((1, 1, 1, 1))
        axis.pane.set_edgecolor("#777777")
        axis.pane.set_linewidth(0.7)
        axis.line.set_color("#777777")
        axis.line.set_linewidth(0.7)
        axis._axinfo["grid"].update(color="#E5E5E5", linewidth=0.8, linestyle="-")
        axis._axinfo["tick"].update(inward_factor=0, outward_factor=0.15)
    ax.tick_params(labelsize=FONT_SIZES["axis_tick"], pad=1.5, colors="#555555")
    ax.set_xlabel(f"PC1 ({variance['PC1']:.2f}% variance)", labelpad=20)
    ax.set_ylabel(f"PC2 ({variance['PC2']:.2f}% variance)", labelpad=10)
    ax.set_zlabel(f"PC3 ({variance['PC3']:.2f}% variance)", labelpad=22)
    means = np.array([xyz[groups == g].mean(0) for g in GROUPS])
    for i, group in enumerate(GROUPS):
        values = xyz[groups == group]
        if show_means:
            ax.scatter(*values.T, marker=MARKERS[group], c="#B8B8B8", s=32,
                       alpha=0.65, edgecolors="none", depthshade=False)
            ax.scatter(*means[i], marker=MARKERS[group], c=COLORS[group], s=155,
                       edgecolors="white", linewidths=0.5, depthshade=False)
        else:
            ax.add_collection3d(Poly3DCollection([values],
                facecolors=COLORS[group], alpha=0.09, edgecolors="none"))
            closed = np.vstack([values, values[0]])
            ax.plot(*closed.T, color=COLORS[group], linewidth=1.6,
                    linestyle=(0, (5, 4)), alpha=0.8)
            ax.scatter(*values.T, marker=MARKERS[group], c=COLORS[group], s=100,
                       edgecolors="white", linewidths=0.4, depthshade=False)
    labels = []
    if not show_means:
        # Screen-space placement changes only text, never sample coordinates.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        projected = np.array(proj3d.proj_transform(*xyz.T, ax.get_proj())[:2]).T
        display_points = ax.transData.transform(projected)
        occupied = [Bbox.from_bounds(x-7, y-7, 14, 14) for x, y in display_points]
        candidates = [(9, 10), (9, -11), (-9, 10), (-9, -11),
                      (15, 0), (-15, 0), (0, 18), (0, -18), (20, 16), (-20, -16)]
        for row, point in zip(rows, projected):
            label = ax.annotate(row["label"], point, xytext=(0, 0),
                textcoords="offset points", fontsize=11, color=COLORS[row["PCA_group"]],
                va="center", annotation_clip=False, zorder=100,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.1))
            choices = []
            for dx, dy in candidates:
                label.set_position((dx, dy))
                label.set_ha("left" if dx > 0 else "right" if dx < 0 else "center")
                bounds = label.get_window_extent(renderer).expanded(1.04, 1.15)
                choices.append((sum(bounds.overlaps(b) for b in occupied), dx, dy, bounds))
                if choices[-1][0] == 0:
                    break
            _, dx, dy, bounds = min(choices, key=lambda item: item[0])
            label.set_position((dx, dy))
            label.set_ha("left" if dx > 0 else "right" if dx < 0 else "center")
            occupied.append(bounds)
            labels.append(label)
    fig.suptitle("Peptide-level PCA\nby Experimental Group", y=0.98,
                 fontsize=FONT_SIZES["title"], fontweight="bold")
    handles = [Line2D([], [], marker=MARKERS[g], color=COLORS[g],
                     linestyle="none", markersize=8, label=g) for g in GROUPS]
    fig.legend(handles=handles, title="Experiment", frameon=False,
               loc="center left", bbox_to_anchor=LEGEND_ANCHOR,
               fontsize=FONT_SIZES["legend_text"], title_fontsize=FONT_SIZES["legend_title"],
               ncol=1, handlelength=0.8, handletextpad=0.35,
               borderaxespad=0, borderpad=0, labelspacing=0.25)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    collisions = [(a.get_text(), b.get_text()) for i, a in enumerate(labels)
                  for b in labels[i+1:] if a.get_window_extent(renderer).overlaps(
                      b.get_window_extent(renderer))]
    return fig, means, limits, collisions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--variance", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tables-dir", required=True, type=Path)
    parser.add_argument("--label", default="CoPPRA")
    parser.add_argument("--elev", type=float, default=30)
    parser.add_argument("--azim", type=float, default=165)
    args = parser.parse_args()
    rows, xyz, groups, variance = read_inputs(args.scores, args.variance)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.tables_dir.mkdir(parents=True, exist_ok=True)
    for show_means, name in [(True, "PCA_peptide_3D_Figure2A_style.png"),
                             (False, "PCA_peptide_3D_grouped.png")]:
        fig, means, limits, collisions = draw(rows, xyz, groups, variance,
            args.elev, args.azim, show_means)
        if collisions:
            raise RuntimeError(f"Overlapping sample labels: {collisions}")
        fig.savefig(args.output_dir / name, dpi=600)
        plt.close(fig)
    with (args.tables_dir / "PCA_3D_group_centroids.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["PCA_group", "n", *PCS])
        writer.writerows([[g, int(np.sum(groups == g)), *m] for g, m in zip(GROUPS, means)])
    settings = {"analysis_label": args.label,
        "views": {
            "Figure2A_style": {"elevation_degrees": args.elev,
                "azimuth_degrees": args.azim},
            "grouped": {"elevation_degrees": args.elev, "azimuth_degrees": args.azim}},
        "view_selection_note": "Camera selected to foreshorten control spread in projection; actual variability is unchanged.",
        "visible_text": "Title, axis labels and ticks, legend; sample labels in the grouped view",
        "figure_size_inches": list(FIGURE_SIZE_INCHES),
        "font_sizes_pt": FONT_SIZES,
        "legend_anchor_figure_fraction": list(LEGEND_ANCHOR),
        "legend_position": "right side, compact single column",
        "projection": "orthographic",
        "equal_data_units": True, "score_coordinates_modified": False,
        "group_means": "Arithmetic means of saved PC scores, n=3 per group",
        "limits": dict(zip(PCS, limits.tolist())),
        "variance_percent": {p: variance[p] for p in PCS},
        "style_reference": "https://doi.org/10.64898/2026.09.24.752349 (v1, Figure 2A)",
        "inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [args.scores, args.variance]},
        "matplotlib": matplotlib.__version__, "numpy": np.__version__,
        "sample_label_overlaps": collisions}
    (args.tables_dir / "PCA_3D_plot_settings.json").write_text(json.dumps(settings, indent=2))
    print(f"Saved {args.label} Figure 2A-style and replicate PNGs; 12 samples; no label overlaps.")


if __name__ == "__main__":
    main()

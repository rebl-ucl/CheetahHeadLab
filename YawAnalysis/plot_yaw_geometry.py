#!/usr/bin/env python3
"""Create exactly two figures from the matched all-18-run yaw/geometry tables."""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'outputs/yaw/geometry'
OUT = SOURCE = HERE
CONTEXTS = ["pitch tracking", "steady lateral", "close pass-by"]
GEOMETRY_TITLES = ["Frontal pursuit", "Lateral pursuit", "Close pass-by pursuit"]
GEOMETRY_PREFIXES = ["F", "L", "C"]
COLORS = ["#0072B2", "#D55E00", "#882255"]
MARKERS = ["o", "s", "^"]
PREFIXES = ["P", "L", "C"]
GREEN = "#009E73"
GRID = np.arange(-45, 200.01, 0.25)
BANDWIDTH = 5.0

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.direction": "out", "ytick.direction": "out",
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def reference_band(ax, horizontal=False):
    if horizontal:
        ax.axhspan(-10, 10, color=GREEN, alpha=0.09, zorder=0)
        ax.axhline(0, color="0.4", lw=0.8, ls=":", zorder=1)
    else:
        ax.axvspan(-10, 10, color=GREEN, alpha=0.09, zorder=0)
        ax.axvline(0, color="0.4", lw=0.8, ls=":", zorder=1)


def density(values):
    z = (GRID[:, None] - np.asarray(values)[None, :]) / BANDWIDTH
    return np.mean(np.exp(-0.5 * z**2), axis=1) / (BANDWIDTH * np.sqrt(2*np.pi))


def save_figure(fig, stem):
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(HERE / f"{stem}.{extension}", dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def context_figure(frames, runs):
    # Match the geometry export's actual canvas, including its tight-crop margins.
    png = (HERE / "02_yaw_lead_vs_geometry.png").read_bytes()
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    target_pixels = struct.unpack(">I", png[16:20])[0]
    svg = ET.parse(HERE / "02_yaw_lead_vs_geometry.svg").getroot()
    vector_width = float(svg.attrib["width"].removesuffix("pt")) / 72
    runs = runs.copy()
    for context, prefix in zip(CONTEXTS, GEOMETRY_PREFIXES):
        indexes = runs.index[runs.context == context]
        runs.loc[indexes, "runId"] = [f"{prefix}{i+1}" for i in range(len(indexes))]
    with plt.rc_context({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8}):
        _context_figure(frames, runs, target_pixels, vector_width)


def _context_figure(frames, runs, target_pixels, vector_width):
    height = 3.2
    fig = plt.figure(figsize=(target_pixels / 300, height))
    gs = fig.add_gridspec(1, 3, wspace=0.28,
                         left=0.105, right=0.975, bottom=0.20, top=0.73)
    xlimits = [(-35, 35), (-25, 110), (-15, 195)]
    xticks = [[-30, -15, 0, 15, 30], [0, 30, 60, 90], [0, 60, 120, 180]]
    density_rows = []
    all_densities = {name: density(f.signedPathGammaDeg) for name, f in frames.groupby("dataset", sort=False)}
    ymax = max(v.max() for v in all_densities.values()) * 1.12
    for k, (context, title, color) in enumerate(zip(CONTEXTS, GEOMETRY_TITLES, COLORS)):
        group = runs[runs.context == context]
        ax = fig.add_subplot(gs[0, k])
        reference_band(ax)
        values = np.stack([all_densities[name] for name in group.dataset], axis=1)
        mean = values.mean(axis=1)
        se = values.std(axis=1, ddof=1) / np.sqrt(len(group))
        for j, name in enumerate(group.dataset):
            ax.plot(GRID, values[:, j], color=color, alpha=0.30, lw=0.65)
        ax.fill_between(GRID, np.maximum(0, mean-se), mean+se, color=color, alpha=0.14, lw=0)
        ax.plot(GRID, mean, color=color, lw=1.8)
        ax.set(xlim=xlimits[k], ylim=(0, ymax))
        ax.set_title(f"{title}\n{int(group.nYawFrames.sum()):,} frames", pad=7)
        ax.grid(axis="y", alpha=0.15, lw=0.6)
        ax.set_xticks(xticks[k])
        ax.yaxis.set_major_locator(MultipleLocator(0.02))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
        if k == 0:
            ax.set_ylabel("Probability density (deg⁻¹)")
        else:
            ax.tick_params(axis="y", labelleft=False)
        observed = frames.loc[frames.context == context, "signedPathGammaDeg"]
        assert observed.between(*xlimits[k]).all(), "Context zoom clips observed yaw values."
        assert (group.yawReconLoDeg >= xlimits[k][0]).all()
        assert (group.yawReconHiDeg <= xlimits[k][1]).all()
        for i, x in enumerate(GRID):
            density_rows.append({"context": context, "signedYawDeg": x,
                                 "runEqualDensityPerDeg": mean[i], "runSE": se[i]})
    fig.supxlabel("Travel-directed head–lure yaw offset (°)", x=0.54, y=0.04, fontsize=9)
    handles = [Line2D([], [], color="0.25", lw=1.8, label="Mean across runs"),
               Line2D([], [], color="0.55", lw=0.65, label="Individual run"),
               Patch(facecolor="0.4", alpha=0.14, label="±1 SE across runs"),
               Patch(facecolor=GREEN, alpha=0.09, edgecolor=GREEN, label="±10° alignment")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, 0.985),
               ncol=4, frameon=False, fontsize=8, handlelength=2.1, columnspacing=1.4)
    pd.DataFrame(density_rows).to_csv(OUT / "context_densities.csv", index=False)
    for extension in ["png", "pdf", "svg"]:
        width = (target_pixels + 1e-6) / 300 if extension == "png" else vector_width
        fig.savefig(HERE / f"01_yaw_lead_by_context.{extension}", dpi=300,
                    bbox_inches=Bbox.from_bounds(0, 0, width, height), pad_inches=0)
    plt.close(fig)


def draw_geometry_points(ax, runs, norm, labels=True, detail=False):
    cmap = plt.get_cmap("viridis")
    # Keep true point coordinates; move only labels for the crowded low-angle runs.
    label_positions = {
        "F1": (12, 6.4), "F2": (8, -2.2),
        "F4": (23, 2.0), "F5": (17.5, -3.8), "F6": (1.0, 1.5),
        "F7": (5, -5.9), "F8": (27, -1.8), "F9": (1.0, 6.8),
        "F10": (23, 6.2), "F11": (19, 9.8),
    }
    updated_run = runs.loc[runs.runId == "F3"].iloc[0]
    label_positions["F3"] = (updated_run.medianThetaDeg + 2.5,
                             updated_run.medianSignedYawDeg - 1.4)
    reference_band(ax, horizontal=True)
    for k, context in enumerate(CONTEXTS):
        group = runs[runs.context == context]
        for row in group.itertuples():
            color = cmap(norm(row.medianRangeM))
            ax.hlines(row.medianSignedYawDeg, row.thetaReconLoDeg, row.thetaReconHiDeg,
                      color=color, alpha=0.5, lw=1, zorder=2)
            ax.vlines(row.medianThetaDeg, row.yawReconLoDeg, row.yawReconHiDeg,
                      color=color, alpha=0.5, lw=1, zorder=2)
            # Fill every marker with distance colour; bars show reconstruction uncertainty.
            ax.scatter(row.medianThetaDeg, row.medianSignedYawDeg, marker=MARKERS[k],
                       facecolors=[color], s=34, edgecolors="0.18",
                       linewidth=0.6, zorder=3)
            if labels and ((detail and ax.get_xlim()[0] <= row.medianThetaDeg <= ax.get_xlim()[1]
                            and ax.get_ylim()[0] <= row.medianSignedYawDeg <= ax.get_ylim()[1])
                           or (not detail and row.context == "close pass-by")):
                label = row.runId
                if detail and row.runId in label_positions:
                    ax.annotate(label, (row.medianThetaDeg, row.medianSignedYawDeg),
                                xytext=label_positions[row.runId], textcoords="data",
                                fontsize=6.5, color="0.18", annotation_clip=True, zorder=4,
                                arrowprops={"arrowstyle": "-", "color": "0.5", "lw": 0.5,
                                            "shrinkA": 2, "shrinkB": 4})
                else:
                    ax.annotate(label, (row.medianThetaDeg, row.medianSignedYawDeg),
                                xytext=(5, 5), textcoords="offset points", fontsize=6.5, color="0.18",
                                annotation_clip=True, zorder=4)
    ax.grid(alpha=0.14, lw=0.6)


def geometry_figure(runs):
    runs = runs.copy()
    for context, prefix, title in zip(CONTEXTS, GEOMETRY_PREFIXES, GEOMETRY_TITLES):
        indexes = runs.index[runs.context == context]
        runs.loc[indexes, "runId"] = [f"{prefix}{i+1}" for i in range(len(indexes))]
        runs.loc[indexes, "displayContext"] = title
    runs[["runId", "dataset", "animal", "context", "displayContext", "nYawFrames",
          "medianThetaDeg", "medianRangeM", "medianSignedYawDeg"]].to_csv(
              HERE / "geometry_run_key.csv", index=False)
    with plt.rc_context({"font.size": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7}):
        _geometry_figure(runs)


def _geometry_figure(runs):
    norm = LogNorm(vmin=runs.medianRangeM.min(), vmax=runs.medianRangeM.max())
    fig = plt.figure(figsize=(7.4, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 0.045], wspace=0.29,
                         left=0.10, right=0.91, bottom=0.18, top=0.86)
    main = fig.add_subplot(gs[0, 0])
    zoom = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])
    low = runs[runs.context != "close pass-by"]
    xmin = max(0, np.floor((low.medianThetaDeg.min()-5)/5)*5)
    xmax = min(180, np.ceil((low.medianThetaDeg.max()+5)/5)*5)
    ymin = np.floor((low.medianSignedYawDeg.min()-4)/5)*5
    ymax = np.ceil((low.medianSignedYawDeg.max()+4)/5)*5
    main.set(xlim=(0, 180), ylim=(-30, 155))
    zoom.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    draw_geometry_points(main, runs, norm)
    draw_geometry_points(zoom, runs, norm, detail=True)
    main.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,
                            fill=False, edgecolor="0.45", lw=0.8, ls="--", zorder=1))
    main.set_xticks([0, 45, 90, 135, 180])
    for ax in [main, zoom]:
        ax.set_xlabel("Lure-travel angle (°)", labelpad=6)
        ax.axhline(10, color=GREEN, alpha=0.35, lw=0.5)
        ax.axhline(-10, color=GREEN, alpha=0.35, lw=0.5)
        ax.text(0.96, 5, "aligned", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8, color="#28785e")
    main.set_ylabel("Travel-directed head–lure yaw offset (°)", labelpad=6)
    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Median horizontal head–lure distance (m)", labelpad=7)
    ticks = [v for v in [0.25, 0.5, 1, 2, 5, 10, 20, 50, 100] if norm.vmin <= v <= norm.vmax]
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{v:g}" for v in ticks])
    cb.ax.minorticks_off()
    handles = [Line2D([], [], marker=m, linestyle="none", markersize=5, color="0.35",
                      label=f"{title} ({int(runs.loc[runs.context == context, 'nYawFrames'].sum()):,} frames)")
               for m, title, context in zip(MARKERS, GEOMETRY_TITLES, CONTEXTS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.51, 1.005),
               ncol=3, frameon=False, fontsize=7.5, columnspacing=1.4)
    save_figure(fig, "02_yaw_lead_vs_geometry")


def main(geometry_only=False, context_only=False):
    frames = pd.read_csv(SOURCE / "signed_yaw_by_frame_all18.csv")
    runs = pd.read_csv(OUT / "geometry_by_run.csv")
    geometry = pd.read_csv(OUT / "geometry_by_frame.csv")
    assert len(runs) == 18 and len(frames) == len(geometry) == 1708
    assert not frames.duplicated(["dataset", "absFrame"]).any()
    joined = frames.merge(geometry, on=["dataset", "absFrame"], validate="one_to_one", suffixes=("", "_geometry"))
    assert len(joined) == len(frames)
    assert np.max(np.abs(joined.signedPathGammaDeg-joined.signedPathGammaDeg_geometry)) < 1e-9
    assert np.all(runs.nYawFrames == runs.nGeometryFrames)
    assert np.all(runs.nDraws == 500)
    assert ((geometry.thetaMedianDeg >= 0) & (geometry.thetaMedianDeg <= 180)).all()
    assert (geometry.rangeMedianM > 0).all()
    for context, prefix in zip(CONTEXTS, PREFIXES):
        indexes = runs.index[runs.context == context]
        runs.loc[indexes, "runId"] = [f"{prefix}{i+1}" for i in range(len(indexes))]
    key = runs[["runId", "dataset", "animal", "context", "nYawFrames", "medianThetaDeg", "medianRangeM", "medianSignedYawDeg"]]
    if not geometry_only and not context_only:
        key.to_csv(HERE / "run_key.csv", index=False)
        runs.to_csv(OUT / "geometry_by_run_with_ids.csv", index=False)
    if not context_only:
        geometry_figure(runs)
    if not geometry_only:
        context_figure(frames, runs)
    key.to_csv(HERE / 'run_key.csv', index=False)
    print('Geometry plots: 18 pursuits, 1708 frames, 500 joint draws.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--geometry-only", action="store_true",
                        help="Regenerate only the geometry figure; preserve the context figure and data tables.")
    selection.add_argument("--context-only", action="store_true",
                           help="Regenerate only the context density panels at the geometry export width.")
    args = parser.parse_args()
    main(args.geometry_only, args.context_only)

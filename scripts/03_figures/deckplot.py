"""
Plot helpers for the slide deck.

Written after three figures came back with clipped titles and bars running
off the top of their axes. Both failures had the same cause: text and data
were placed at coordinates chosen by eye, with nothing reserving room for
them. The helpers here make that room structurally.

  panel()    returns axes sitting below a reserved title band, so a title
             can never overlap the plot or fall outside the canvas
  headroom() sets limits with padding in both directions, so a positive
             bar in an otherwise negative series has somewhere to go
  card()     fixed-width blocks; width carries no meaning, which is what
             lets a 2.6% category keep a readable label
  note()     statistics inside the axes, because panel() sets titles in
             caps and would turn r into R
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch

EMBER, GOLD, FLAME = "#E8402A", "#F2A623", "#F27127"
BONE, ASH, DUSK = "#F4EFEA", "#A29591", "#7A6C68"
GREEN, RUST, SLATE = "#6F9E44", "#8E3020", "#5A4744"
GRID, INK = "#3A2F2D", "#241A18"

_have = {f.name for f in fm.fontManager.ttflist}
HEAD = next((f for f in ["TeX Gyre Heros Cn", "Liberation Sans Narrow",
                         "DejaVu Sans Condensed", "Liberation Sans",
                         "DejaVu Sans"] if f in _have), "DejaVu Sans")
BODY = next((f for f in ["Carlito", "Open Sans", "Liberation Sans",
                         "DejaVu Sans"] if f in _have), "DejaVu Sans")

DPI = 300
TITLE_BAND = 0.13


def figure(w, h):
    f = plt.figure(figsize=(w, h), dpi=DPI)
    f.patch.set_alpha(0)
    return f


def panel(fig, left, width, bottom=0.16, height=None, title=None):
    if height is None:
        height = 1.0 - bottom - TITLE_BAND - 0.04
    ax = fig.add_axes([left, bottom, width, height])
    ax.set_facecolor((0, 0, 0, 0))
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(colors=DUSK, labelsize=11, length=0)
    if title:
        fig.text(left, bottom + height + 0.055, title.upper(), family=HEAD,
                 fontsize=14, color=BONE, fontweight="bold", va="bottom")
    return ax


def finish(ax, xlabel=None, ylabel=None, grid="y"):
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontfamily(BODY)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel, color=ASH, fontsize=12, family=BODY, labelpad=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=ASH, fontsize=12, family=BODY, labelpad=9)


def headroom(ax, values, pad=0.16, axis="y", zero=True):
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    lo, hi = v.min(), v.max()
    if zero:
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    span = hi - lo or 1.0
    lim = (lo - span * pad, hi + span * pad)
    (ax.set_ylim if axis == "y" else ax.set_xlim)(*lim)
    return lim


def bar_labels(ax, xs, values, fmt="{:.3f}", size=11.5, color=BONE):
    lo, hi = ax.get_ylim()
    off = (hi - lo) * 0.028
    for x, v in zip(xs, values):
        ax.text(x, v + (off if v >= 0 else -off), fmt.format(v),
                family=BODY, fontsize=size, color=color, ha="center",
                va="bottom" if v >= 0 else "top", zorder=6)


def card(ax, x, y, w, h, title, value, accent=EMBER, fill=INK,
         tsize=11.5, vsize=13):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=0.018",
                                facecolor=fill, edgecolor=accent, lw=1.6,
                                zorder=3))
    ax.add_patch(FancyBboxPatch((x, y), min(0.012, w * 0.05), h,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                facecolor=accent, edgecolor="none", zorder=4))
    ax.text(x + w * 0.09, y + h * 0.68, title.upper(), family=HEAD,
            fontsize=tsize, color=ASH, va="center", fontweight="bold",
            zorder=5)
    ax.text(x + w * 0.09, y + h * 0.28, value, family=BODY, fontsize=vsize,
            color=BONE, va="center", zorder=5)


def note(ax, text, loc="upper left", size=12.5, color=GOLD, pad=0.045):
    xa = pad if "left" in loc else 1 - pad
    ya = 1 - pad if "upper" in loc else pad
    ax.text(xa, ya, text, transform=ax.transAxes, family=BODY,
            fontsize=size, color=color, fontweight="bold",
            ha="left" if "left" in loc else "right",
            va="top" if "upper" in loc else "bottom", zorder=8)


def blank(fig, rect=(0.03, 0.04, 0.94, 0.92)):
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax


def save(fig, path):
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=DPI, transparent=True, bbox_inches="tight",
                pad_inches=0.16)
    plt.close(fig)
    return path

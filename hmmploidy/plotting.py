"""Plotting utilities: translation of the original ``hmmPlotting`` R
function into matplotlib."""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np


def plot_ploidy(
    window_start: np.ndarray,
    window_end: np.ndarray,
    states: np.ndarray,
    posterior: np.ndarray,
    depth: np.ndarray,
    mu: np.ndarray,
    out_path: str,
    true_ploidy: Optional[Sequence[int]] = None,
    title: str = "Inferred ploidies",
) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [2, 1]}
    )

    ax = axes[0]
    x = np.arange(len(states))
    ax.step(x, states, where="mid", color="navyblue" if False else "#1f3d7a", lw=1.5, label="Inferred ploidy")
    if true_ploidy is not None:
        ax.plot(x, true_ploidy, "s", color="coral", ms=3, label="True ploidy")

    prop_states = posterior.mean(axis=0)
    present_states = np.unique(states)
    for s in present_states:
        y0 = s
        ax.fill_between(
            x, y0, y0 + posterior[:, s - 1] * 0.5, step="mid", color="deepskyblue", alpha=0.5
        )

    ax.set_ylabel("Ploidy")
    ax.set_title(title)
    ax.legend(loc="upper left")

    ax2 = axes[1]
    ax2.plot(x, depth, ".", color="deepskyblue", ms=3, label="Window mean depth")
    for s in present_states:
        ax2.axhline(mu[s - 1], color="coral", lw=1, alpha=0.7)
    ax2.set_ylabel("Depth")
    ax2.set_xlabel("Window index")
    ax2.set_title("Window-mean depth")
    ax2.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)

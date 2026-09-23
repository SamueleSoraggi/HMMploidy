"""Simulation benchmark of the HMM models under different kinds of
overdispersion (used by ``notebooks/demo_HMMploidy.ipynb``).

Each scenario simulates a population with known ploidies, computes genotype
likelihoods and fits every requested model to every individual, scoring the
fraction of windows whose ploidy is called correctly.
"""
from __future__ import annotations

import argparse
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import data as dm
from . import genotype_likelihoods as glmod
from . import hmm
from . import simulate as sim


@dataclass
class Scenario:
    description: str
    sim_kwargs: Dict = field(default_factory=dict)


SCENARIOS: Dict[str, Scenario] = {
    "poisson": Scenario("Poisson depth, Binomial allele counts (no overdispersion)",
                        dict(depth_model="poisson")),
    "nb1": Scenario("Strongly overdispersed per-site depth (NB size 1)",
                    dict(depth_model="negbinom", dispersion=1.0)),
    "cv30": Scenario("Regional coverage bias, CV 30% over 50-locus blocks",
                     dict(depth_model="negbinom", dispersion=5.0, depth_cv=0.3)),
    "ab05": Scenario("Beta-Binomial allele counts, rho = 0.05",
                     dict(depth_model="negbinom", dispersion=5.0, allele_overdispersion=0.05)),
    "ab15": Scenario("Beta-Binomial allele counts, rho = 0.15",
                     dict(depth_model="negbinom", dispersion=5.0, allele_overdispersion=0.15)),
    "refbias": Scenario("Reference bias: non-reference reads under-sampled (odds 0.7)",
                        dict(depth_model="negbinom", dispersion=5.0, ref_bias=0.7)),
    "hard": Scenario("All combined: NB size 2 + regional CV 30% + rho 0.1",
                     dict(depth_model="negbinom", dispersion=2.0, depth_cv=0.3,
                          allele_overdispersion=0.1)),
    "hardrb": Scenario("'hard' + reference bias (odds 0.8)",
                       dict(depth_model="negbinom", dispersion=2.0, depth_cv=0.3,
                            allele_overdispersion=0.1, ref_bias=0.8)),
}


def ploidy_layout(layout: str, n_samples: int) -> List[List[int]]:
    """Ploidy of each sample in each of 3 genome segments.

    ``changing``: every sample changes ploidy along the genome (2/4/3 for the
    first half of the samples, 3/2/4 for the second).
    ``uniform``: every sample keeps one ploidy (2, 3 or 4) genome-wide, so
    depth alone cannot tell ploidies apart.
    """
    half = n_samples // 2
    if layout == "uniform":
        per_sample = [2 + i % 3 for i in range(n_samples)]
        return [per_sample] * 3
    if layout == "changing":
        return [
            [2] * half + [3] * (n_samples - half),
            [4] * half + [2] * (n_samples - half),
            [3] * half + [4] * (n_samples - half),
        ]
    raise ValueError(f"Unknown layout {layout!r}")


def simulate_scenario(
    name: str,
    outdir: str,
    depth: float = 5.0,
    layout: str = "changing",
    loci_per_segment: int = 15000,
    n_samples: int = 6,
    seed: int = 1,
    threads: int = 1,
) -> str:
    """Simulate one scenario and compute its genotype likelihoods; returns
    the ``.genolikes.gz`` path (cached: not recomputed if present)."""
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.join(outdir, f"{name}.{layout}")
    segments = ploidy_layout(layout, n_samples)
    gl_path = f"{prefix}.DP{float(depth)}.genolikes.gz"
    if os.path.exists(gl_path):
        return gl_path
    ploidy_file = f"{prefix}.ploidy"
    with open(ploidy_file, "w") as fh:
        fh.write("\n".join(",".join(map(str, s)) for s in segments) + "\n")
    sim.simulate_dataset(
        ploidy_file, [loci_per_segment] * 3, [depth], prefix, seed=seed,
        **SCENARIOS[name].sim_kwargs,
    )
    args = argparse.Namespace(
        min_quality_score=1, min_non_major_freq=0.2, max_minor2_freq=0.1,
        max_minor3_freq=0.1, min_global_depth=0, min_ind_depth=-1,
        downsampling=1, random_seed=seed, threads=threads, chunk_size=20000,
    )
    glmod.run_on_file(f"{prefix}.DP{float(depth)}.mpileup.gz", gl_path, [1, 2, 3, 4, 5, 6], args)
    return gl_path


def evaluate(
    gl_path: str,
    segments: Sequence[Sequence[int]],
    loci_per_segment: int,
    models: Sequence[str] = ("overdispersed", "negbinom", "em"),
    wind: int = 100,
    max_ploidy: int = 6,
    model_kwargs: Optional[Dict[str, dict]] = None,
) -> pd.DataFrame:
    """Fit every model to every individual; one row per (individual, model)."""
    model_kwargs = model_kwargs or {}
    df = dm.read_genolikes(gl_path)
    pop = dm.population_counts(df)
    rows = []
    for ind in range(1, int(df["ind"].max()) + 1):
        iw = dm.prepare_individual(df, ind, tuple(range(1, 7)), max_ploidy, wind, pop=pop)
        centre = (iw.window_start + iw.window_end) / 2
        seg = np.clip((centre - 1) // loci_per_segment, 0, len(segments) - 1).astype(int)
        truth = np.array([segments[k][ind - 1] for k in seg])
        for model in models:
            t0 = time.perf_counter()
            res = hmm.fit(model, iw.depth, iw.log_emission, snps=iw.snps, **model_kwargs.get(model, {}))
            rows.append(dict(
                individual=ind, model=model,
                accuracy=float(np.mean(res["states"] == truth)),
                fit_seconds=time.perf_counter() - t0, windows=len(iw.depth),
                rho=res.get("rho", np.nan),
            ))
    return pd.DataFrame(rows)


def run_benchmark(
    scenarios: Sequence[str],
    outdir: str,
    depth: float = 5.0,
    layout: str = "changing",
    loci_per_segment: int = 15000,
    n_samples: int = 6,
    models: Sequence[str] = ("overdispersed", "negbinom", "em"),
    threads: int = 1,
    seed: int = 1,
) -> pd.DataFrame:
    """Simulate + evaluate several scenarios; one row per (scenario,
    individual, model)."""
    frames = []
    for name in scenarios:
        gl_path = simulate_scenario(
            name, outdir, depth, layout, loci_per_segment, n_samples, seed, threads
        )
        res = evaluate(gl_path, ploidy_layout(layout, n_samples), loci_per_segment, models)
        res.insert(0, "scenario", name)
        res.insert(1, "depth", depth)
        res.insert(2, "layout", layout)
        frames.append(res)
    return pd.concat(frames, ignore_index=True)

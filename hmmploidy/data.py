"""Loading and windowing of ``.genolikes`` files.

Pure-Python/numpy translation of the data-preparation logic that used to
live inline in ``HMMploidy.R`` (population allele-frequency estimation,
per-individual SNP/­depth filtering, window construction and per-window
aggregation of depth and genotype likelihoods). All the per-window
aggregation ``sapply`` loops in the R script have been replaced with
vectorised ``numpy.searchsorted`` + cumulative-sum lookups.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import dataclasses
import os

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import binom


def resolve_genolikes_path(path: str) -> str:
    """Return ``path``, or ``path + '.gz'`` if only the compressed file exists."""
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        return path + ".gz"
    return path


def read_genolikes(path: str, ploidy_list: Sequence[int] = (1, 2, 3, 4, 5, 6)) -> pd.DataFrame:
    """Read a ``.genolikes`` or ``.genolikes.gz`` file into a DataFrame with
    named columns: chrom, pos, ind, ref, depth, major, minor, major_count,
    minor_count, and one ``gl_<ploidy>_<dosage>`` column per GL entry.

    Genotype likelihoods are stored as float32 and counts as int32 to halve
    memory use on genome-scale files."""
    gl_names = []
    for p in ploidy_list:
        gl_names.extend([f"gl_{p}_{k}" for k in range(p + 1)])
    names = [
        "chrom", "pos", "ind", "ref", "depth", "major", "minor",
        "major_count", "minor_count",
    ] + gl_names
    dtypes = {
        "chrom": "category", "pos": np.int64, "ind": np.int32, "ref": "category",
        "depth": np.int32, "major": "category", "minor": "category",
        "major_count": np.int32, "minor_count": np.int32,
        **{c: np.float32 for c in gl_names},
    }
    path = resolve_genolikes_path(path)
    df = pd.read_csv(
        path, sep="\t", header=None, names=names, usecols=range(len(names)), dtype=dtypes
    )
    check_site_major_layout(df)
    return df


def check_site_major_layout(df: pd.DataFrame) -> None:
    """The file must list every individual (1..n, in order) at every site,
    as written by ``hmmploidy genolikes``; the per-site reshapes below rely
    on it."""
    ind = df["ind"].to_numpy()
    n_ind = int(ind.max()) if len(ind) else 0
    if len(ind) % max(n_ind, 1) or not np.array_equal(
        ind.reshape(-1, n_ind), np.broadcast_to(np.arange(1, n_ind + 1), (len(ind) // n_ind, n_ind))
    ):
        raise ValueError(
            "genolikes rows must list individuals 1..n in order at every site "
            "(as written by `hmmploidy genolikes`)"
        )


@dataclass
class PopulationCounts:
    """Per-site read counts summed over all individuals of a file."""

    minor: np.ndarray      # (n_sites,) sum of minor-allele reads
    total: np.ndarray      # (n_sites,) sum of major + minor reads
    n_covered: np.ndarray  # (n_sites,) nr of individuals with >= 1 read
    n_ind: int


def population_counts(df: pd.DataFrame) -> PopulationCounts:
    """Compute once per file and pass to :func:`prepare_individual`."""
    n_ind = int(df["ind"].max())
    major = df["major_count"].to_numpy().reshape(-1, n_ind).astype(np.int64)
    minor = df["minor_count"].to_numpy().reshape(-1, n_ind).astype(np.int64)
    nt = major + minor
    return PopulationCounts(minor.sum(axis=1), nt.sum(axis=1), (nt > 0).sum(axis=1), n_ind)


def _other_counts(pop: PopulationCounts, minor_i: np.ndarray, total_i: np.ndarray):
    if pop.n_ind < 2:
        return pop.minor, pop.total, pop.n_covered
    return pop.minor - minor_i, pop.total - total_i, pop.n_covered - (total_i > 0)


def frequency_beta_prior(
    pop: PopulationCounts,
    minor_i: np.ndarray,
    total_i: np.ndarray,
    freq: np.ndarray,
    chromosomes_per_ind: float = 2.0,
):
    """Beta(a, b) parameters expressing the uncertainty of the leave-one-out
    frequency ``freq``. Reads of one individual are not independent allele
    draws, so the effective sample size is capped at
    ``chromosomes_per_ind`` per covered individual (plus a Jeffreys 1/2).
    Used by the ``overdispersed`` model to put a Beta-Binomial rather than
    Binomial prior on the allele dosage, which keeps noisy frequencies from
    few individuals from biasing the ploidy calls."""
    _, total_o, covered_o = _other_counts(pop, minor_i, total_i)
    ess = np.minimum(total_o, chromosomes_per_ind * covered_o).astype(float)
    f = np.clip(freq, 1e-3, 1 - 1e-3)
    return f * ess + 0.5, (1 - f) * ess + 0.5


def leave_one_out_frequencies(
    pop: PopulationCounts,
    minor_i: np.ndarray,
    total_i: np.ndarray,
    min_ind: int = 1,
    eps: float = 0.0005,
) -> np.ndarray:
    """Population minor-allele frequency per site *excluding* one
    individual (same read-weighted estimator as
    :func:`estimate_allele_frequencies`).

    The frequency is the prior on the individual's allele dosages, so it
    must not be estimated from the individual's own reads: in small
    samples a high-depth (e.g. high-ploidy) individual dominates the
    estimate, the prior then tracks its own allele ratios, and the finer
    dosage grid of the next ploidy up fits them better - a systematic
    upward bias. Falls back to including the individual when the file has
    a single individual. Sites with fewer than ``min_ind`` informative
    other individuals get -1.
    """
    minor_o, total_o, covered_o = _other_counts(pop, minor_i, total_i)
    with np.errstate(divide="ignore", invalid="ignore"):
        freq = (minor_o - eps * total_o) / (total_o * (1 - 2 * eps))
    return np.where((covered_o >= min_ind) & (total_o > 0), freq, -1.0)


def estimate_allele_frequencies(
    major_counts: np.ndarray, minor_counts: np.ndarray, min_ind: int, eps: float = 0.0005
) -> np.ndarray:
    """Population allele-frequency estimate per site, pooling over
    individuals with a read-depth weighted average. ``major_counts`` and
    ``minor_counts`` have shape (n_sites, n_ind). Sites with fewer than
    ``min_ind`` informative individuals get frequency -1 (filtered out
    downstream), matching the original C++ ``alleleFrequencies``.
    :func:`leave_one_out_frequencies` computes the same estimator from
    per-file sums (see :func:`population_counts`)."""
    nt = major_counts + minor_counts
    tot = nt.sum(axis=1, keepdims=True)
    tot_safe = np.where(tot == 0, 1, tot)
    with np.errstate(divide="ignore", invalid="ignore"):
        pis = np.where(nt > 0, (minor_counts - eps * nt) / (nt * (1 - 2 * eps)), 0.0)
    wis = np.where(nt > 0, nt / tot_safe, 0.0)
    ind_with_data = (nt > 0).sum(axis=1)
    freq = (wis * pis).sum(axis=1)
    freq = np.where(ind_with_data >= min_ind, freq, -1.0)
    return freq


def logrescale(x: np.ndarray) -> np.ndarray:
    """Log-normalise log-likelihoods along the last axis (log-sum-exp),
    matching the original ``logRescale``."""
    out = x - logsumexp(x, axis=-1, keepdims=True)
    return np.clip(out, None, -1e-5)


def windows_builder(min_size: int, loci: np.ndarray, loci_snp: np.ndarray) -> np.ndarray:
    """Build windows of (approximately) ``min_size`` loci along a
    chromosome, merging windows that end up without any SNP so that every
    window contains at least one SNP. Returns an (n_windows, 2) array of
    (start, end) boundaries (in the same coordinate as `loci`)."""
    loci = np.asarray(loci)
    loci_snp = np.asarray(loci_snp)

    if len(loci) > 2:
        edges = np.arange(loci[0], loci[-2] + 1, min_size, dtype=float)
        s = np.concatenate([edges, [float(loci[-1])]])
    elif len(loci) == 2:
        s = np.array([float(loci[0]), float(loci[1])])
    else:
        s = np.array([float(loci[0]), float(loci[0])])

    if loci.max() - loci.min() < min_size:
        s = np.array([float(loci.min()), float(loci.max())])

    # Merge windows without SNPs into the following window (the last ones
    # into the preceding window). Equivalent to deleting boundaries one at a
    # time, but in closed form: an interior boundary survives iff the
    # original window ending at it, (s[j-1], s[j]], contains a SNP.
    n_le = np.searchsorted(np.sort(loci_snp), s, side="right")  # nr SNPs <= boundary
    keep = np.ones(len(s), dtype=bool)
    keep[1:-1] = n_le[1:-1] > n_le[:-2]
    kept = np.flatnonzero(keep)
    if len(kept) > 2 and n_le[-1] == n_le[kept[-2]]:
        keep[kept[-2]] = False  # trailing empty window merges backwards
    s = s[keep]
    return np.stack([s[:-1], s[1:]], axis=1)


def window_sum(values: np.ndarray, windows: np.ndarray, loci: np.ndarray) -> np.ndarray:
    """Sum ``values`` (aligned with sorted ``loci``) within each window,
    inclusive of both boundaries. Vectorised with a cumulative sum + binary
    search instead of a per-window scan."""
    order = np.argsort(loci)
    loci_sorted = loci[order]
    cum = np.concatenate([[0.0], np.nancumsum(values[order])])
    left = np.searchsorted(loci_sorted, windows[:, 0], side="left")
    right = np.searchsorted(loci_sorted, windows[:, 1], side="right")
    return cum[right] - cum[left]


def window_mean(values: np.ndarray, windows: np.ndarray, loci: np.ndarray) -> np.ndarray:
    order = np.argsort(loci)
    loci_sorted = loci[order]
    vals_sorted = values[order]
    cum = np.concatenate([[0.0], np.nancumsum(vals_sorted)])
    left = np.searchsorted(loci_sorted, windows[:, 0], side="left")
    right = np.searchsorted(loci_sorted, windows[:, 1], side="right")
    counts = np.maximum(right - left, 1)
    return (cum[right] - cum[left]) / counts


def geno_loglik_for_ploidy(freq: np.ndarray, gl: np.ndarray) -> np.ndarray:
    """log P(reads | ploidy) at each SNP, marginalising the (unknown)
    allele dosage under a Binomial(ploidy, freq) prior, i.e. the
    translation of the original ``pGenoDataAll``.

    ``gl`` has shape (n_snp, ploidy + 1): natural-log genotype
    likelihoods for dosage 0..ploidy. ``freq`` has shape (n_snp,).
    """
    ploidy = gl.shape[1] - 1
    k = np.arange(ploidy + 1)
    log_prior = binom.logpmf(k[None, :], ploidy, np.clip(freq, 1e-12, 1 - 1e-12)[:, None])
    log_prior = np.where(np.isfinite(log_prior), log_prior, -1000.0)
    return logsumexp(gl + log_prior, axis=1)


@dataclass
class SnpData:
    """Per-SNP allele counts of one individual, for emission models that
    compute their own genotype likelihoods (see
    :mod:`hmmploidy.hmm.overdispersed`)."""

    window: np.ndarray         # (n_snp,) index of the window containing the SNP
    minor: np.ndarray          # (n_snp,) minor-allele read count
    total: np.ndarray          # (n_snp,) major + minor read count
    freq: np.ndarray           # (n_snp,) population minor-allele frequency
    error_rate: float          # per-read probability of reading a given wrong allele
    # Beta(prior_a, prior_b) uncertainty of `freq`, from the other individuals
    prior_a: Optional[np.ndarray] = None
    prior_b: Optional[np.ndarray] = None


@dataclass
class IndividualWindows:
    chrom: np.ndarray          # window -> chromosome name
    window_start: np.ndarray
    window_end: np.ndarray
    depth: np.ndarray          # (n_windows,) mean depth per window
    log_emission: np.ndarray   # (n_windows, max_ploidy) rescaled log P(data|ploidy)
    sites_snp: np.ndarray      # positions of SNPs actually used
    snps: Optional[SnpData] = None


def prepare_individual(
    df: pd.DataFrame,
    which_ind: int,
    ploidy_list: Sequence[int],
    max_ploidy: int,
    wind: int,
    min_ind: int = 1,
    eps: float = 0.0005,
    quantile_trim: Tuple[float, float] = (0.02, 0.98),
    snp_trim: Tuple[float, float] = (0.1, 0.9),
    pop: Optional[PopulationCounts] = None,
) -> IndividualWindows:
    """Build per-window depth + genotype-likelihood matrices for one
    individual, following the per-file / per-individual block of
    ``HMMploidy.R`` (but with leave-one-out population allele frequencies,
    see :func:`leave_one_out_frequencies`). ``pop`` (from
    :func:`population_counts`) can be passed to avoid recomputing it for
    every individual."""
    pop = population_counts(df) if pop is None else pop

    n_ind = int(df["ind"].max())
    ind_df = df.iloc[which_ind - 1 :: n_ind].reset_index(drop=True)
    minor_all = ind_df["minor_count"].to_numpy()
    total_all = ind_df["major_count"].to_numpy() + minor_all
    # SNP selection and the precomputed genotype likelihoods of the legacy
    # models use the all-individual frequency, as in HMMploidy.R; the
    # overdispersed model gets leave-one-out frequencies with uncertainty
    freqs_all = leave_one_out_frequencies(
        dataclasses.replace(pop, n_ind=1), minor_all, total_all, min_ind, eps
    )
    loo_all = leave_one_out_frequencies(pop, minor_all, total_all, min_ind, eps)
    prior_a_all, prior_b_all = frequency_beta_prior(pop, minor_all, total_all, loo_all)
    depth = ind_df["depth"].to_numpy(dtype=float)
    pos = ind_df["pos"].to_numpy(dtype=float)
    chrom = ind_df["chrom"].to_numpy().astype(str)

    lo, hi = np.quantile(depth, quantile_trim)
    keep = (depth >= lo) & (depth <= hi)
    depth, pos, chrom = depth[keep], pos[keep], chrom[keep]
    freqs_ind = freqs_all[keep]
    loo_ind = loo_all[keep]
    prior_a_snp = prior_a_all[keep]
    prior_b_snp = prior_b_all[keep]
    major_i = ind_df["major_count"].to_numpy()[keep]
    minor_i = ind_df["minor_count"].to_numpy()[keep]

    # SNPs are selected on the *population* allele frequency, as in the
    # original HMMploidy.R. Selecting on the individual's own allele ratio
    # (i.e. keeping only sites that look heterozygous) biases the genotype
    # likelihoods towards high ploidies, whose heterozygosity is higher.
    is_snp = (freqs_ind > snp_trim[0]) & (freqs_ind < snp_trim[1])

    pos_snp = pos[is_snp]
    freq_snp = freqs_ind[is_snp]
    loo_snp = loo_ind[is_snp]

    gl_cols = {p: [c for c in df.columns if c.startswith(f"gl_{p}_")] for p in ploidy_list}
    for pp in range(1, max_ploidy + 1):
        if pp not in gl_cols:
            raise ValueError(f"Ploidy {pp} (<= maxPloidy={max_ploidy}) not present in genolikes columns")
    snp_rows = np.flatnonzero(keep)[is_snp]
    geno_snp = np.stack(
        [
            geno_loglik_for_ploidy(freq_snp, ind_df[gl_cols[pp]].to_numpy(dtype=float)[snp_rows])
            for pp in range(1, max_ploidy + 1)
        ]
    )  # (max_ploidy, n_snp)

    # contigs are contiguous blocks of rows (the input is a sorted pileup)
    block_starts = np.flatnonzero(np.r_[True, chrom[1:] != chrom[:-1]])
    block_ends = np.r_[block_starts[1:], len(chrom)]
    snp_index = np.flatnonzero(is_snp)
    snp_starts = np.searchsorted(snp_index, block_starts)
    snp_ends = np.searchsorted(snp_index, block_ends)
    windows_chrom, w_start, w_end, w_depth = [], [], [], []
    log_emission_rows: List[np.ndarray] = []
    used_sites: List[float] = []
    snp_window: List[np.ndarray] = []
    snp_minor: List[np.ndarray] = []
    snp_total: List[np.ndarray] = []
    snp_freq: List[np.ndarray] = []
    snp_a: List[np.ndarray] = []
    snp_b: List[np.ndarray] = []
    prior_a_snp, prior_b_snp = prior_a_snp[is_snp], prior_b_snp[is_snp]
    major_snp = major_i[is_snp]
    minor_snp = minor_i[is_snp]
    n_windows_done = 0

    for b0, b1, sb0, sb1 in zip(block_starts, block_ends, snp_starts, snp_ends):
        c = chrom[b0]
        mask_indiv = slice(b0, b1)
        mask_snp = slice(sb0, sb1)
        loci_indiv = pos[mask_indiv]
        loci_snp = pos_snp[mask_snp]
        if len(loci_snp) == 0 or len(loci_indiv) < 2:
            continue
        wt = windows_builder(wind, loci_indiv, loci_snp)
        geno_per_ploidy = geno_snp[:, mask_snp]

        geno_win = np.stack(
            [window_sum(geno_per_ploidy[pp], wt, loci_snp) for pp in range(max_ploidy)], axis=1
        )
        depth_win = window_mean(depth[mask_indiv], wt, loci_indiv)

        keep_win = np.all(np.isfinite(geno_win), axis=1)
        geno_win = geno_win[keep_win]
        depth_win = depth_win[keep_win]

        # assign each SNP to exactly one window: (end_{i-1}, end_i]
        win_of_snp = np.searchsorted(wt[:, 1], loci_snp, side="left")
        win_of_snp = np.minimum(win_of_snp, len(wt) - 1)
        new_index = np.cumsum(keep_win) - 1
        snp_keep = keep_win[win_of_snp]
        snp_window.append(new_index[win_of_snp[snp_keep]] + n_windows_done)
        snp_minor.append(minor_snp[mask_snp][snp_keep])
        snp_total.append(minor_snp[mask_snp][snp_keep] + major_snp[mask_snp][snp_keep])
        snp_freq.append(loo_snp[mask_snp][snp_keep])
        snp_a.append(prior_a_snp[mask_snp][snp_keep])
        snp_b.append(prior_b_snp[mask_snp][snp_keep])
        wt = wt[keep_win]
        n_windows_done += len(wt)

        geno_win_rescaled = logrescale(geno_win)

        windows_chrom.extend([c] * len(wt))
        w_start.extend(wt[:, 0].tolist())
        w_end.extend(wt[:, 1].tolist())
        w_depth.extend(depth_win.tolist())
        log_emission_rows.append(geno_win_rescaled)
        used_sites.extend(loci_snp.tolist())

    log_emission = (
        np.concatenate(log_emission_rows, axis=0) if log_emission_rows else np.zeros((0, max_ploidy))
    )

    def _cat(parts, dtype):
        return np.concatenate(parts).astype(dtype) if parts else np.zeros(0, dtype=dtype)

    # reads that are neither major nor minor are sequencing errors into one
    # of the two other bases; the rate into one specific base is half that
    depth_i = ind_df["depth"].to_numpy(dtype=float)[keep]
    other = np.clip(depth_i - major_i - minor_i, 0, None).sum()
    error_rate = float(np.clip(0.5 * other / max(depth_i.sum(), 1.0), 1e-4, 0.05))
    snps = SnpData(
        window=_cat(snp_window, np.int64),
        minor=_cat(snp_minor, float),
        total=_cat(snp_total, float),
        freq=_cat(snp_freq, float),
        error_rate=error_rate,
        prior_a=_cat(snp_a, float),
        prior_b=_cat(snp_b, float),
    )
    return IndividualWindows(
        chrom=np.array(windows_chrom),
        window_start=np.array(w_start),
        window_end=np.array(w_end),
        depth=np.array(w_depth),
        log_emission=log_emission,
        sites_snp=np.array(used_sites),
        snps=snps,
    )

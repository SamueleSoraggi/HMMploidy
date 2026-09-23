"""Simulate polyploid sequencing data in mpileup format.

Pure-Python, vectorised re-implementation of ``simulMpileup.R`` +
``simulationScript.sh``. For every site the read-base sampling loop (the
performance bottleneck in the original R script) is vectorised with numpy
across all reads of a sample at once.

In addition to the original Poisson read-depth model, a Negative-Binomial
depth model is available (``--depth-model negbinom``) to generate
over-dispersed depths -- useful to test how the three HMM emission models
(``negbinom``, ``em``, ``bayesian``, see :mod:`hmmploidy.hmm`) behave under
matched or mismatched depth-generating processes. ``--depth-cv`` adds a
regional coverage bias and ``--allele-overdispersion`` Beta-Binomial allele
counts, the two kinds of overdispersion that do not average out in windows.
"""
from __future__ import annotations

import argparse
import gzip
import os
from dataclasses import dataclass
from typing import List, Sequence, TextIO

import numpy as np

REF_ALLELE = "A"
NONREF_ALLELE = "C"
OTHER_ALLELES = ["G", "T"]


@dataclass
class SegmentConfig:
    """One segment of the simulated genome: ploidy per sample + nr loci."""

    ncopy: np.ndarray  # ploidy per sample, shape (n_samples,)
    n_loci: int


def parse_copy_spec(spec: str) -> np.ndarray:
    """Parse a ploidy spec like '2x10,4x5' -> array([2]*10 + [4]*5)."""
    values: List[int] = []
    for token in spec.split(","):
        parts = token.split("x")
        if len(parts) == 2:
            values.extend([int(parts[0])] * int(parts[1]))
        else:
            values.append(int(parts[0]))
    return np.asarray(values, dtype=int)


def parse_ploidy_file(path: str) -> List[np.ndarray]:
    with open(path) as fh:
        return [parse_copy_spec(line.strip()) for line in fh if line.strip()]


def sample_depths(
    n_loci: int,
    haploid_depth: float,
    ncopy: np.ndarray,
    err_depth: float,
    depth_model: str,
    dispersion: float,
    rng: np.random.Generator,
    depth_cv: float = 0.0,
    depth_block: int = 50,
) -> np.ndarray:
    """Depths of shape (n_samples, n_loci) for one segment."""
    n_samples = len(ncopy)
    lam = haploid_depth * ncopy[:, None]
    lo = lam * (1 - err_depth)
    hi = lam * (1 + err_depth)
    lam_sampled = rng.uniform(lo, hi, size=(n_samples, n_loci))
    lam_sampled = np.clip(lam_sampled, 0.01, None)
    if depth_cv > 0:
        # Regional coverage bias (GC content, mappability, ...): a log-normal
        # multiplicative factor shared by blocks of `depth_block` consecutive
        # loci. Unlike per-site noise, this does not average out in windows.
        n_blocks = -(-n_loci // depth_block)
        sigma = np.sqrt(np.log1p(depth_cv**2))
        factor = rng.lognormal(-sigma**2 / 2, sigma, size=(n_samples, n_blocks))
        lam_sampled = lam_sampled * np.repeat(factor, depth_block, axis=1)[:, :n_loci]
    if depth_model == "poisson":
        return rng.poisson(lam_sampled)
    if depth_model == "negbinom":
        r = dispersion
        p = r / (r + lam_sampled)
        return rng.negative_binomial(r, p)
    raise ValueError(f"Unknown depth model: {depth_model}")


def sample_allele_frequencies(
    n_loci: int, ne: int, ksfs: float, pvar: float, panc: float, rng: np.random.Generator
) -> np.ndarray:
    """Derived allele frequency per site, following the same neutral SFS
    parameterisation as the original R script."""
    ee = 1.0 / (np.arange(1, ne)) ** (1.0 / ksfs)
    ee = ee / ee.sum() * pvar
    pder = np.concatenate(([(1 - pvar) * 0.9], ee, [(1 - pvar) * 0.1]))
    pder = pder / pder.sum()
    states = np.arange(0, ne + 1)
    qq = rng.choice(states, size=n_loci, p=pder) / ne
    flip_prob = 1 - panc
    flip = rng.random(n_loci) < flip_prob
    qq[flip] = 1 - qq[flip]
    return qq


def _sample_bases_for_sample(
    dosage: int,
    ploidy: int,
    depth: int,
    mean_qual: int,
    rng: np.random.Generator,
    allele_overdispersion: float = 0.0,
    ref_bias: float = 1.0,
) -> tuple[str, str]:
    """Vectorised read sampling for one sample at one site.

    With ``allele_overdispersion`` (rho) > 0, the expected non-reference read
    fraction at heterozygous sites is drawn from a Beta distribution with
    mean ``dosage/ploidy`` and intra-class correlation rho, so that allele
    counts are Beta-Binomial instead of Binomial. With ``ref_bias`` < 1,
    non-reference reads at heterozygous sites are under-sampled by that
    odds ratio (reference mapping bias: a shift of the mean allele ratio).
    """
    if depth <= 0:
        return ",", "!"

    q = np.round(rng.normal(loc=mean_qual, scale=2.0, size=depth))
    q = np.clip(q, 0, None)
    p_err = 10.0 ** (-(q / 10.0))

    frac_nonref = dosage / ploidy
    if ref_bias != 1.0 and 0 < dosage < ploidy:
        frac_nonref = frac_nonref * ref_bias / (frac_nonref * ref_bias + 1 - frac_nonref)
    if allele_overdispersion > 0 and 0 < dosage < ploidy:
        conc = (1 - allele_overdispersion) / allele_overdispersion
        frac_nonref = rng.beta(frac_nonref * conc, (1 - frac_nonref) * conc)
    frac_ref = 1.0 - frac_nonref
    prob_nonref = frac_nonref * (1 - p_err) + frac_ref * (p_err / 3.0)
    prob_ref = frac_ref * (1 - p_err) + frac_nonref * (p_err / 3.0)
    prob_other = p_err / 3.0

    # cumulative probs over [ref, nonref, other1, other2] per read
    cum = np.stack(
        [prob_ref, prob_ref + prob_nonref, prob_ref + prob_nonref + prob_other],
        axis=1,
    )
    u = rng.random(depth)
    choice = np.sum(u[:, None] >= cum, axis=1)  # 0=ref,1=nonref,2=other1,3=other2
    symbol_map = np.array([".", NONREF_ALLELE, OTHER_ALLELES[0], OTHER_ALLELES[1]])
    bases = "".join(symbol_map[choice])
    quals = "".join(chr(33 + int(x)) for x in q)
    return bases, quals


def simulate_segment(
    ncopy: np.ndarray,
    n_loci: int,
    haploid_depth: float,
    mean_qual: int,
    ksfs: float,
    pvar: float,
    panc: float,
    ne: int,
    err_depth: float,
    depth_model: str,
    dispersion: float,
    offset: int,
    rng: np.random.Generator,
    depth_cv: float = 0.0,
    depth_block: int = 50,
    allele_overdispersion: float = 0.0,
    ref_bias: float = 1.0,
) -> List[str]:
    """Simulate one segment; returns a list of mpileup text lines."""
    n_samples = len(ncopy)
    depths = sample_depths(
        n_loci, haploid_depth, ncopy, err_depth, depth_model, dispersion, rng,
        depth_cv=depth_cv, depth_block=depth_block,
    )
    qq = sample_allele_frequencies(n_loci, ne, ksfs, pvar, panc, rng)

    lines: List[str] = []
    for i in range(n_loci):
        fields = [f"copy_{'_'.join(map(str, ncopy))}", str(offset + i), REF_ALLELE]
        for n in range(n_samples):
            dosage = int(rng.binomial(int(ncopy[n]), qq[i]))
            depth = int(depths[n, i])
            bases, quals = _sample_bases_for_sample(
                dosage, int(ncopy[n]), depth, mean_qual, rng, allele_overdispersion, ref_bias
            )
            fields.extend([str(max(depth, 1)), bases, quals])
        lines.append("\t".join(fields))
    return lines


def simulate_dataset(
    ploidy_file: str,
    loci: Sequence[int],
    depths: Sequence[float],
    out_prefix: str,
    qual: int = 20,
    ksfs: float = 1.0,
    pvar: float = 1.0,
    panc: float = 1.0,
    ne: int = 10000,
    err_depth: float = 0.05,
    depth_model: str = "poisson",
    dispersion: float = 5.0,
    seed: int = 180218,
    depth_cv: float = 0.0,
    depth_block: int = 50,
    allele_overdispersion: float = 0.0,
    ref_bias: float = 1.0,
) -> List[str]:
    """Simulate a full dataset (all segments x all depths), mirroring
    ``simulationScript.sh``. Returns the list of generated file basenames
    (one per depth), each with a corresponding ``.mpileup.gz`` and a
    ``<out_prefix>.filelist``.
    """
    segments = parse_ploidy_file(ploidy_file)
    if len(segments) != len(loci):
        raise ValueError(
            f"ploidy file has {len(segments)} segments but {len(loci)} loci counts were given"
        )

    rng = np.random.default_rng(seed)
    out_names = []
    for depth in depths:
        out_name = f"{out_prefix}.DP{depth}"
        out_path = f"{out_name}.mpileup.gz"
        offset = 1
        with gzip.open(out_path, "wt") as fh:
            for ncopy, n_loci in zip(segments, loci):
                lines = simulate_segment(
                    ncopy,
                    n_loci,
                    depth,
                    qual,
                    ksfs,
                    pvar,
                    panc,
                    ne,
                    err_depth,
                    depth_model,
                    dispersion,
                    offset,
                    rng,
                    depth_cv=depth_cv,
                    depth_block=depth_block,
                    allele_overdispersion=allele_overdispersion,
                    ref_bias=ref_bias,
                )
                fh.write("\n".join(lines) + "\n")
                offset += n_loci
        out_names.append(out_name)

    filelist_path = f"{out_prefix}.filelist"
    with open(filelist_path, "w") as fh:
        for name in out_names:
            fh.write(name + ".mpileup.gz\n")
    return out_names


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hmmploidy simulate",
        description="Simulate polyploid mpileup data for testing HMMploidy",
    )
    parser.add_argument("-p", "--ploidyFile", required=True, help="Segment ploidy file")
    parser.add_argument(
        "-l", "--loci", required=True, help="Comma-separated nr of loci per segment"
    )
    parser.add_argument(
        "-d", "--depth", required=True, help="Comma-separated haploid depths to simulate"
    )
    parser.add_argument("-o", "--out", default="out", help="Output file prefix")
    parser.add_argument("-q", "--qual", type=int, default=20, help="Mean base quality")
    parser.add_argument("-k", "--ksfs", type=float, default=1.0)
    parser.add_argument("-r", "--pvar", type=float, default=1.0)
    parser.add_argument("-a", "--panc", type=float, default=1.0)
    parser.add_argument("-n", "--ne", type=int, default=10000)
    parser.add_argument("-e", "--errdepth", type=float, default=0.05)
    parser.add_argument(
        "--depth-model",
        choices=["poisson", "negbinom"],
        default="poisson",
        help="Distribution used to generate per-site sequencing depth",
    )
    parser.add_argument(
        "--dispersion",
        type=float,
        default=5.0,
        help="Negative-Binomial dispersion (size) parameter, used only with --depth-model negbinom",
    )
    parser.add_argument(
        "--depth-cv",
        type=float,
        default=0.0,
        help="Coefficient of variation of a regional (block-wise) log-normal coverage bias",
    )
    parser.add_argument(
        "--depth-block", type=int, default=50, help="Nr of loci sharing one coverage-bias factor"
    )
    parser.add_argument(
        "--allele-overdispersion",
        type=float,
        default=0.0,
        help="Beta-Binomial intra-class correlation (rho) of allele read counts at heterozygous sites",
    )
    parser.add_argument(
        "--ref-bias",
        type=float,
        default=1.0,
        help="Odds ratio (<1) by which non-reference reads are under-sampled at heterozygous sites",
    )
    parser.add_argument("-u", "--seed", type=int, default=180218)
    return parser


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    loci = [int(x) for x in args.loci.split(",")]
    depths = [float(x) for x in args.depth.split(",")]
    out_names = simulate_dataset(
        args.ploidyFile,
        loci,
        depths,
        args.out,
        qual=args.qual,
        ksfs=args.ksfs,
        pvar=args.pvar,
        panc=args.panc,
        ne=args.ne,
        err_depth=args.errdepth,
        depth_model=args.depth_model,
        dispersion=args.dispersion,
        seed=args.seed,
        depth_cv=args.depth_cv,
        depth_block=args.depth_block,
        allele_overdispersion=args.allele_overdispersion,
        ref_bias=args.ref_bias,
    )
    print(f"Generated: {', '.join(n + '.mpileup.gz' for n in out_names)}")
    print(f"File list: {args.out}.filelist")


if __name__ == "__main__":
    main()

"""Main ploidy-inference pipeline: replaces the top-level driver logic of
``HMMploidy.R``, dispatching to one of three HMM fitting techniques.
"""
from __future__ import annotations

import argparse
import multiprocessing
import os
from typing import List, Optional, Sequence

import numpy as np

from . import data as dm
from . import hmm as hmm_module
from .plotting import plot_ploidy


def analyze_file(
    genolikes_path: str,
    out_prefix: str,
    model: str = "overdispersed",
    wind: int = 100,
    max_ploidy: int = 6,
    ploidy_list: Sequence[int] = (1, 2, 3, 4, 5, 6),
    chosen_ind: Optional[Sequence[int]] = None,
    min_ind: int = 1,
    eps: float = 0.0005,
    quantile_trim=(0.02, 0.98),
    snp_trim=(0.1, 0.9),
    names: Optional[Sequence[str]] = None,
    true_ploidy: Optional[int] = None,
    model_kwargs: Optional[dict] = None,
    threads: int = 1,
    plot: bool = True,
) -> List[dict]:
    """Run ploidy inference on one ``.genolikes`` (or ``.genolikes.gz``)
    file, for all (or the chosen) individuals. Writes
    ``<out_prefix>.<name>.HMMploidy`` (text) and ``<out_prefix>.<name>.pdf``
    (plot) per individual, and returns the list of result dictionaries (one
    per individual) for programmatic use (e.g. from the demo notebook).

    With ``threads > 1`` individuals are fitted in parallel worker
    processes; the genolikes table itself stays in the main process, which
    prepares one individual at a time and sends only its windows."""
    model_kwargs = model_kwargs or {}
    df = dm.read_genolikes(genolikes_path, ploidy_list=ploidy_list)
    pop = dm.population_counts(df)
    n_ind = int(df["ind"].max())
    individuals = list(chosen_ind) if chosen_ind else list(range(1, n_ind + 1))

    def tasks():
        for ind in individuals:
            name = names[ind - 1] if names else f"ind{ind}"
            iw = dm.prepare_individual(
                df, ind, ploidy_list, max_ploidy, wind, min_ind=min_ind, eps=eps,
                quantile_trim=quantile_trim, snp_trim=snp_trim, pop=pop,
            )
            if len(iw.depth) == 0:
                print(f"[{name}] no usable windows, skipping")
                continue
            yield (name, iw, model, model_kwargs, out_prefix, max_ploidy, true_ploidy, plot)

    if threads <= 1 or model == "bayesian":
        return [_fit_individual(task) for task in tasks()]
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(threads) as pool:
        return list(pool.imap(_fit_individual, tasks(), chunksize=1))


def _fit_individual(task) -> dict:
    """Fit, write and plot one individual (runs in a worker with threads > 1)."""
    name, iw, model, model_kwargs, out_prefix, max_ploidy, true_ploidy, plot = task
    if multiprocessing.current_process().name != "MainProcess":
        import matplotlib

        matplotlib.use("Agg")  # no GUI / notebook backend in worker processes
    res = hmm_module.fit(model, iw.depth, iw.log_emission, snps=iw.snps, **model_kwargs)
    res["name"] = name
    res["window_start"] = iw.window_start
    res["window_end"] = iw.window_end
    res["chrom"] = iw.chrom
    res["depth"] = iw.depth

    _write_text_output(f"{out_prefix}.{name}.HMMploidy", name, res, iw, true_ploidy)
    if plot:
        true_series = np.full(len(res["states"]), true_ploidy) if true_ploidy else None
        plot_ploidy(
            iw.window_start,
            iw.window_end,
            res["states"],
            res.get("posterior", np.zeros((len(res["states"]), max_ploidy)) + 1 / max_ploidy),
            iw.depth,
            res["mu"],
            f"{out_prefix}.{name}.pdf",
            true_ploidy=true_series,
            title=f"{name} ({model} model) - Inferred ploidies",
        )
    print(f"[{name}] done: {len(iw.depth)} windows")
    return res


def _write_text_output(path, name, res, iw, true_ploidy):
    with open(path, "w") as fh:
        fh.write(f"# sample: {name}\n")
        fh.write(f"# model: {res['model']}\n")
        fh.write("delta\t" + "\t".join(f"{x:.6g}" for x in res["delta"]) + "\n")
        fh.write("trans\t" + "\t".join(f"{x:.6g}" for x in np.asarray(res["trans"]).ravel()) + "\n")
        fh.write("mu\t" + "\t".join(f"{x:.6g}" for x in res["mu"]) + "\n")
        if "rho" in res:
            fh.write(f"haploid_depth\t{res['haploid_depth']:.6g}\n")
            fh.write(f"depth_shape\t{res['depth_shape']:.6g}\n")
            fh.write(f"allele_overdispersion_rho\t{res['rho']:.6g}\n")
            fh.write(f"error_rate\t{res['error_rate']:.6g}\n")
        if "alpha" in res:
            fh.write("alpha\t" + "\t".join(f"{x:.6g}" for x in res["alpha"]) + "\n")
            fh.write("beta\t" + "\t".join(f"{x:.6g}" for x in res["beta"]) + "\n")
        if "loglik" in res and len(np.atleast_1d(res["loglik"])):
            fh.write(f"final_loglik\t{np.atleast_1d(res['loglik'])[-1]:.6g}\n")
        fh.write("ploidy\t" + "\t".join(str(int(x)) for x in res["states"]) + "\n")
        if true_ploidy is not None:
            fh.write(f"true_ploidy\t{true_ploidy}\n")
        prop = np.bincount(res["states"], minlength=len(res["mu"]) + 1)[1:] / len(res["states"])
        fh.write("proportions\t" + "\t".join(f"{x:.4g}" for x in prop) + "\n")
        fh.write("window_start\t" + "\t".join(str(int(x)) for x in iw.window_start) + "\n")
        fh.write("window_end\t" + "\t".join(str(int(x)) for x in iw.window_end) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hmmploidy run", description="Infer ploidy levels from .genolikes files"
    )
    parser.add_argument("--file", help="Single basename of a .genolikes(.gz) file")
    parser.add_argument("--fileList", help="File with a list of basenames")
    parser.add_argument(
        "--model", choices=list(hmm_module.MODEL_CHOICES), default="overdispersed",
        help="HMM fitting technique to use",
    )
    parser.add_argument("--wind", type=int, required=True, help="Window size (bp)")
    parser.add_argument("--maxPloidy", type=int, default=6)
    parser.add_argument("--ploidyList", default=None, help="Comma separated ploidies present in the genolikes columns (default 1..maxPloidy)")
    parser.add_argument("--minInd", type=int, default=1)
    parser.add_argument("--eps", type=float, default=0.0005)
    parser.add_argument("--quantileTrim", default="0.02,0.98")
    parser.add_argument("--SNPtrim", default="0.1,0.9")
    parser.add_argument("--chosenInd", default=None, help="Comma separated 1-indexed individuals")
    parser.add_argument("--nameFile", default=None)
    parser.add_argument("--truePl", type=int, default=None)
    parser.add_argument("--outSuffix", default=None)
    parser.add_argument(
        "-t", "--threads", type=int, default=min(8, os.cpu_count() or 1),
        help="Individuals fitted in parallel (default: min(8, nr of CPUs))",
    )
    parser.add_argument("--noPlot", action="store_true", help="Skip the per-individual PDF plots")
    return parser


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    if args.file is None and args.fileList is None:
        raise SystemExit("Provide either --file or --fileList")

    basenames = [args.file] if args.file else [
        line.strip() for line in open(args.fileList) if line.strip()
    ]

    ploidy_list = (
        tuple(int(x) for x in args.ploidyList.split(","))
        if args.ploidyList
        else tuple(range(1, args.maxPloidy + 1))
    )
    chosen_ind = (
        [int(x) for x in args.chosenInd.split(",")] if args.chosenInd else None
    )
    names = (
        [line.strip() for line in open(args.nameFile) if line.strip()]
        if args.nameFile
        else None
    )
    quantile_trim = tuple(float(x) for x in args.quantileTrim.split(","))
    snp_trim = tuple(float(x) for x in args.SNPtrim.split(","))

    for base in basenames:
        # accept genolikes basenames as well as the mpileup list given to `genolikes`
        for suffix in (".genolikes.gz", ".genolikes", ".mpileup.gz", ".mpileup"):
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break
        genolikes_path = dm.resolve_genolikes_path(f"{base}.genolikes")
        out_prefix = base if not args.outSuffix else f"{base}.{args.outSuffix}"
        print(f"==> Analyze {base}")
        analyze_file(
            genolikes_path,
            out_prefix,
            model=args.model,
            wind=args.wind,
            max_ploidy=args.maxPloidy,
            ploidy_list=ploidy_list,
            chosen_ind=chosen_ind,
            min_ind=args.minInd,
            eps=args.eps,
            quantile_trim=quantile_trim,
            snp_trim=snp_trim,
            names=names,
            true_ploidy=args.truePl,
            threads=args.threads,
            plot=not args.noPlot,
        )


if __name__ == "__main__":
    main()

"""Genotype-likelihood calculation from mpileup/BAM input.

Python re-implementation of the original ``Genotype_Likelihoods.py``. The
output file format is unchanged: one ``<name>.genolikes.gz`` file per input
mpileup (or one for a whole list of BAM files, see below), tab separated,
with columns::

    chrom  pos  ind  ref  depth  major  minor  major_count  minor_count \
        GL_ploidy1(2 cols) GL_ploidy2(3 cols) ... GL_ploidyN(N+1 cols)

where the set of ploidies matches ``--ploidyFile`` (default ``1..6``).

Input can be:

* multi-sample mpileup / mpileup.gz files (one output per file), or
* a list of BAM files (``--bamList``), one per sample. These are piled up
  *jointly* with ``samtools mpileup`` (streamed, no temporary file), so the
  population allele frequencies used downstream are estimated from all
  samples. When the BAMs are indexed, the genome is split into regions that
  are piled up and processed in parallel.

Each site is processed by the compiled kernel in :mod:`hmmploidy._gl_kernel`;
:func:`process_site` below is the equivalent, readable pure-Python
reference implementation. Sites are processed in parallel chunks
(``--threads``), each written as an independent gzip member.

Note: the ``--Inbreeding`` option is accepted for backward CLI
compatibility with the original script but (as in the original) is not
currently folded into the likelihood model.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import math
import multiprocessing as mp_ctx
import os
import random
import shlex
import shutil
import subprocess
import sys
from typing import Iterable, Iterator, List, Optional, Sequence, TextIO, Tuple

import numpy as np

from . import mpileup as mp

ALLELES = mp.ALLELES


def process_site(
    line: str,
    nsams: int,
    ploidies: Sequence[int],
    min_quality_score: int,
    min_non_major_freq: float,
    max_minor2_freq: float,
    max_minor3_freq: float,
    min_global_depth: float,
    min_ind_depth: float,
    downsampling: float,
    rng: random.Random,
) -> List[str]:
    """Return the list of output lines (one per sample) for one site, or
    an empty list if the site is filtered out."""
    site = mp.parse_pileup_line(line, nsams)

    decoded = [
        mp.decode_pileup_field(b, q, site.reference)
        for b, q in zip(site.sample_bases, site.sample_quals)
    ]
    if min_quality_score > 0:
        decoded = [
            mp.filter_by_quality(b, q, min_quality_score) for b, q in decoded
        ]

    individual_depth = np.array([len(b) for b, _ in decoded], dtype=float)
    pooled_bases = "".join(b for b, _ in decoded)
    pooled_quals = "".join(q for _, q in decoded)
    global_depth = len(pooled_bases)

    if not (global_depth > min_global_depth and individual_depth.min() > min_ind_depth):
        return []

    non_major_count = global_depth - mp.base_counts(pooled_bases).max()
    non_major_prop = non_major_count / global_depth
    if non_major_prop <= min_non_major_freq:
        return []

    haploid = mp.haploid_log_likelihoods(pooled_bases, pooled_quals)[:4]
    # rank on likelihoods rounded to 1e-8 so that exact ties (e.g. two
    # alleles with the same reads) do not depend on floating-point
    # summation order; ties go to the lower allele index (A < C < G < T)
    order = sorted(range(4), key=lambda i: -math.floor(haploid[i] * 1e8 + 0.5))
    major, minor, minor2, minor3 = order

    minor2_prop = mp.base_counts(pooled_bases)[minor2] / global_depth
    minor3_prop = mp.base_counts(pooled_bases)[minor3] / global_depth
    if not (minor2_prop < max_minor2_freq and minor3_prop < max_minor3_freq):
        return []

    out_lines = []
    for n in range(nsams):
        bases, quals = decoded[n]
        if downsampling < 1:
            k = math.ceil(len(bases) * downsampling)
            if k < len(bases):
                idx = rng.sample(range(len(bases)), k)
                bases = "".join(bases[i] for i in idx)
                quals = "".join(quals[i] for i in idx)

        counts = mp.base_counts(bases)
        major_count = int(counts[major])
        minor_count = int(counts[minor])
        sample_depth = len(bases)

        fields = [
            site.chrom,
            str(site.position),
            str(n + 1),
            site.reference,
            str(sample_depth),
            ALLELES[major],
            ALLELES[minor],
            str(major_count),
            str(minor_count),
        ]
        for ploidy in ploidies:
            gl = mp.genotype_log_likelihoods(ploidy, bases, quals, major, minor)
            fields.extend(f"{x}" for x in gl)
        out_lines.append("\t".join(fields) + "\n")
    return out_lines


def _kernel_params(args: argparse.Namespace) -> tuple:
    return (
        int(args.min_quality_score),
        float(args.min_non_major_freq),
        float(args.max_minor2_freq),
        float(args.max_minor3_freq),
        float(args.min_global_depth),
        float(args.min_ind_depth),
        float(args.downsampling),
    )


def _process_lines(task) -> bytes:
    """Worker: turn a chunk of raw mpileup lines into a gzip member of
    ``.genolikes`` lines."""
    lines, nsams, ploidies, params, seed = task
    from ._gl_kernel import process_line, seed_kernel_rng

    if seed is not None:
        seed_kernel_rng(seed)
    ploidies = np.asarray(ploidies, dtype=np.int64)
    out: List[str] = []
    for raw in lines:
        raw = raw.rstrip(b"\r\n")
        if not raw.strip():
            continue
        keep, major, minor, depth, major_c, minor_c, gls = process_line(
            np.frombuffer(raw, dtype=np.uint8), nsams, ploidies, *params
        )
        if not keep:
            continue
        chrom, pos, ref = raw.split(b"\t", 3)[:3]
        prefix = f"{chrom.decode()}\t{pos.decode()}\t"
        ref = ref.decode()
        a_major, a_minor = ALLELES[major], ALLELES[minor]
        for n in range(nsams):
            out.append(
                f"{prefix}{n + 1}\t{ref}\t{depth[n]}\t{a_major}\t{a_minor}\t"
                f"{major_c[n]}\t{minor_c[n]}\t" + "\t".join(map(str, gls[n].tolist()))
            )
    if not out:
        return b""
    return gzip.compress(("\n".join(out) + "\n").encode(), compresslevel=4)


def _chunks(lines: Iterable[bytes], size: int) -> Iterator[List[bytes]]:
    chunk: List[bytes] = []
    for line in lines:
        chunk.append(line)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _chunk_seed(args: argparse.Namespace, index: int) -> Optional[int]:
    # one seed per chunk, so results do not depend on the number of threads
    return None if args.random_seed is None else (int(args.random_seed) * 1_000_003 + index) % 2**32


def _write_members(tasks: Iterable[tuple], output_path: str, threads: int) -> None:
    with open(output_path, "wb") as out:
        if threads <= 1:
            for task in tasks:
                out.write(_process_lines(task))
        else:
            ctx = mp_ctx.get_context("spawn")
            with ctx.Pool(threads) as pool:
                for member in pool.imap(_process_lines, tasks, chunksize=1):
                    out.write(member)


def _count_samples_bytes(line: bytes) -> int:
    return (len(line.rstrip(b"\r\n").split(b"\t")) - 3) // 3


def process_stream(
    lines: Iterable[bytes],
    output_path: str,
    ploidies: Sequence[int],
    args: argparse.Namespace,
    nsams: Optional[int] = None,
) -> None:
    """Compute genotype likelihoods for a stream of raw mpileup lines."""
    lines = iter(lines)
    first = next(lines, None)
    if first is None:
        open(output_path, "wb").close()
        return
    if nsams is None:
        nsams = _count_samples_bytes(first)
    params = _kernel_params(args)
    threads = int(getattr(args, "threads", 1) or 1)
    chunk_size = int(getattr(args, "chunk_size", 20000) or 20000)

    def tasks():
        for i, chunk in enumerate(_chunks(itertools.chain([first], lines), chunk_size)):
            yield (chunk, nsams, tuple(ploidies), params, _chunk_seed(args, i))

    _write_members(tasks(), output_path, threads)


def run_on_file(
    input_path: str,
    output_path: str,
    ploidies: Sequence[int],
    args: argparse.Namespace,
) -> None:
    """Compute genotype likelihoods for one (multi-sample) mpileup file."""
    opener = gzip.open if input_path.endswith(".gz") else open
    with opener(input_path, "rb") as fh:
        process_stream(fh, output_path, ploidies, args)


# --------------------------------------------------------------------------
# BAM input
# --------------------------------------------------------------------------


def _require_samtools() -> str:
    # also look next to the Python interpreter, i.e. in a non-activated conda env
    exe = shutil.which("samtools") or shutil.which(
        "samtools", path=os.path.dirname(sys.executable)
    )
    if exe is None:
        raise SystemExit(
            "'samtools' not found; it is needed for BAM input. Install it with "
            "`conda install -c bioconda samtools`."
        )
    return exe


def mpileup_command(
    bams: Sequence[str], ref_fasta: str, args: argparse.Namespace, region: Optional[str] = None
) -> List[str]:
    cmd = [_require_samtools(), "mpileup", "-f", ref_fasta]
    if getattr(args, "min_mapq", 0):
        cmd += ["-q", str(args.min_mapq)]
    if getattr(args, "mpileup_args", None):
        cmd += shlex.split(args.mpileup_args)
    if region is not None:
        cmd += ["-r", region]
    return cmd + list(bams)


def _is_indexed(bam: str) -> bool:
    return any(os.path.exists(bam + ext) for ext in (".bai", ".csi")) or os.path.exists(
        bam[: -len(".bam")] + ".bai"
    )


def bam_contigs(bam: str) -> List[Tuple[str, int]]:
    """(name, length) of the reference sequences in a BAM header."""
    header = subprocess.run(
        [_require_samtools(), "view", "-H", bam], check=True, capture_output=True, text=True
    ).stdout
    contigs = []
    for line in header.splitlines():
        if line.startswith("@SQ"):
            tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
            contigs.append((tags["SN"], int(tags["LN"])))
    return contigs


def split_regions(contigs: Sequence[Tuple[str, int]], region_size: int) -> List[str]:
    regions = []
    for name, length in contigs:
        label = f"{{{name}}}" if ":" in name else name
        for start in range(1, length + 1, region_size):
            regions.append(f"{label}:{start}-{min(start + region_size - 1, length)}")
    return regions


def _process_region(task) -> bytes:
    """Worker: pile up one region of all BAMs and process it."""
    bams, ref_fasta, args, region, nsams, ploidies, index = task
    cmd = mpileup_command(bams, ref_fasta, args, region)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    params = _kernel_params(args)
    members = []
    for i, chunk in enumerate(_chunks(proc.stdout, int(args.chunk_size))):
        seed = _chunk_seed(args, index * 100_000 + i)
        members.append(_process_lines((chunk, nsams, tuple(ploidies), params, seed)))
    proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"samtools mpileup failed on region {region}: {' '.join(cmd)}")
    return b"".join(members)


def run_on_bams(
    bams: Sequence[str],
    ref_fasta: str,
    output_path: str,
    ploidies: Sequence[int],
    args: argparse.Namespace,
) -> None:
    """Pile up all BAMs jointly (one sample per BAM) and compute genotype
    likelihoods. Uses region-level parallelism when every BAM is indexed."""
    threads = int(getattr(args, "threads", 1) or 1)
    nsams = len(bams)
    if threads > 1 and all(_is_indexed(b) for b in bams):
        regions = split_regions(bam_contigs(bams[0]), int(args.region_size))
        print(f"Piling up {nsams} BAM files in {len(regions)} regions with {threads} threads")
        tasks = [
            (list(bams), ref_fasta, args, r, nsams, tuple(ploidies), i)
            for i, r in enumerate(regions)
        ]
        ctx = mp_ctx.get_context("spawn")
        with open(output_path, "wb") as out, ctx.Pool(threads) as pool:
            for member in pool.imap(_process_region, tasks, chunksize=1):
                out.write(member)
        return

    if threads > 1:
        print("BAM files are not all indexed (samtools index): piling up as a single stream")
    cmd = mpileup_command(bams, ref_fasta, args)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    try:
        process_stream(proc.stdout, output_path, ploidies, args, nsams=nsams)
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise SystemExit(f"samtools mpileup failed: {' '.join(cmd)}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hmmploidy genolikes",
        description="Calculate genotype likelihoods from mpileup or BAM files",
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="File with the list of mpileup / mpileup.gz file names to use (one output per file). "
        "A list made only of .bam files is treated like --bamList.",
    )
    parser.add_argument(
        "-b", "--bamList",
        help="File with the list of BAM files, one per sample, piled up jointly into ONE "
        "genolikes file (requires samtools and --ref_fasta)",
    )
    parser.add_argument("-r", "--ref_fasta", help="Reference fasta, required for BAM input")
    parser.add_argument(
        "-O", "--outName", default=None,
        help="Output basename for --bamList input (default: the list file name without extension)",
    )
    parser.add_argument("-p", "--ploidyFile", help="File listing ploidy levels to test")
    parser.add_argument("-o", "--outFolder", default=None, help="Output folder")
    parser.add_argument(
        "-i", "--Inbreeding", help="Per-sample inbreeding coefficients (reserved, unused)"
    )
    parser.add_argument("-d", "--downsampling", default=1, help="Fraction of reads to keep")
    parser.add_argument(
        "-m", "--min_non_major_freq", type=float, default=0.2,
        help="Minimum frequency of non-major alleles for a site to be kept",
    )
    parser.add_argument(
        "-q", "--min_quality_score", type=int, default=1,
        help="Minimum base quality (phred) for a read to be kept",
    )
    parser.add_argument("-dp", "--min_global_depth", type=float, default=0)
    parser.add_argument(
        "-dpInd", "--min_ind_depth", type=float, default=-1,
        help="Sites are kept only if every individual has depth > this value. The default (-1) "
        "keeps sites where some individuals have no reads, which is essential for many "
        "low-coverage samples; their genotype likelihoods are then uninformative (all 0)",
    )
    parser.add_argument("-M2", "--max_minor2_freq", type=float, default=0.1)
    parser.add_argument("-M3", "--max_minor3_freq", type=float, default=0.1)
    parser.add_argument("-s", "--random_seed", type=int, default=None)
    parser.add_argument(
        "-t", "--threads", type=int, default=min(8, os.cpu_count() or 1),
        help="Number of worker processes (default: min(8, nr of CPUs))",
    )
    parser.add_argument(
        "--chunk_size", type=int, default=20000, help="Nr of mpileup lines per work unit"
    )
    parser.add_argument(
        "--region_size", type=int, default=10_000_000,
        help="Size (bp) of the genomic regions piled up in parallel for indexed BAMs",
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help="Minimum mapping quality (samtools mpileup -q)"
    )
    parser.add_argument(
        "--mpileup_args", default=None,
        help="Extra arguments passed verbatim to samtools mpileup, e.g. \"-B -Q 0\"",
    )
    return parser


def _read_ploidies(args) -> List[int]:
    if args.ploidyFile:
        with open(args.ploidyFile) as fh:
            ploidies = [int(line.strip()) for line in fh if line.strip()]
        print(f"Ploidy levels to be tested in analysis are: {ploidies}")
    else:
        ploidies = [1, 2, 3, 4, 5, 6]
        print(f"Default ploidy levels to be tested in analysis are: {ploidies}")
    return ploidies


def _read_list(path: str) -> List[str]:
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _run_bam_list(bams: List[str], list_path: str, ploidies, args) -> None:
    if not args.ref_fasta:
        sys.exit("--ref_fasta is required for BAM input")
    base = args.outName or os.path.splitext(os.path.basename(list_path))[0]
    folder = args.outFolder or os.path.dirname(os.path.abspath(list_path))
    output_path = os.path.join(folder, base + ".genolikes.gz")
    names_path = os.path.join(folder, base + ".names")
    with open(names_path, "w") as fh:
        for bam in bams:
            fh.write(os.path.basename(bam)[: -len(".bam")] + "\n")
    print(f"{len(bams)} BAM file(s), piled up jointly (one sample each)")
    print(f"Output file is: {output_path} (sample names in {names_path}, usable as --nameFile)")
    run_on_bams(bams, args.ref_fasta, output_path, ploidies, args)


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    if args.input is None and args.bamList is None:
        sys.exit("Provide a list of mpileup files, or --bamList")

    if args.random_seed is not None:
        print(f"Seed for calculations using random function is: {args.random_seed}")
    else:
        print("Seed is not set.")
    ploidies = _read_ploidies(args)

    if args.bamList:
        _run_bam_list(_read_list(args.bamList), args.bamList, ploidies, args)
        return

    names = _read_list(args.input)
    print(f"{len(names)} file(s) found")
    print(names)
    if names and all(n.endswith(".bam") for n in names):
        _run_bam_list(names, args.input, ploidies, args)
        return

    for name in names:
        if name.endswith(".mpileup.gz"):
            base = name[: -len(".mpileup.gz")]
        elif name.endswith(".mpileup"):
            base = name[: -len(".mpileup")]
        elif name.endswith(".bam"):
            sys.exit(
                f"{name}: BAM files cannot be mixed with mpileup files; list them with --bamList"
            )
        else:
            sys.exit(f"Unsupported file type: {name}")

        out_name = os.path.basename(base) + ".genolikes.gz"
        output_path = os.path.join(args.outFolder, out_name) if args.outFolder else base + ".genolikes.gz"
        print(f"Output file is: {output_path}")
        run_on_file(name, output_path, ploidies, args)


if __name__ == "__main__":
    main()

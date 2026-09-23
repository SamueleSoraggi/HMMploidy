"""Shared fixtures: small simulated datasets, built once per test session."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from hmmploidy import benchmark

LOCI_PER_SEGMENT = 5000
BAM_PLOIDY = [2, 2, 4, 3]


@pytest.fixture(scope="session")
def sim_dir(tmp_path_factory):
    return str(tmp_path_factory.mktemp("sim"))


@pytest.fixture(scope="session")
def hard_dataset(sim_dir):
    """'hard' scenario (all kinds of overdispersion), 6 individuals whose
    ploidy changes along the genome. Returns (genolikes path, segments)."""
    path = benchmark.simulate_scenario(
        "hard", sim_dir, depth=5.0, loci_per_segment=LOCI_PER_SEGMENT, seed=11
    )
    return path, benchmark.ploidy_layout("changing", 6)


@pytest.fixture(scope="session")
def mpileup_path(hard_dataset):
    return hard_dataset[0].replace(".genolikes.gz", ".mpileup.gz")


def find_samtools():
    return shutil.which("samtools") or shutil.which(
        "samtools", path=os.path.dirname(sys.executable)
    )


@pytest.fixture(scope="session")
def bam_dataset(tmp_path_factory):
    """Four indexed BAMs (ploidy 2, 2, 4, 3) on a reference whose second
    contig is partly soft-masked. Returns (bam list file, reference)."""
    samtools = find_samtools()
    if samtools is None:
        pytest.skip("samtools not available")
    d = tmp_path_factory.mktemp("bam")
    rng = np.random.default_rng(7)
    contigs = {"chr1": 60_000, "chr2": 40_000}
    B = np.array(list("ACGT"))
    ref = {c: rng.integers(0, 4, L) for c, L in contigs.items()}
    with open(d / "ref.fa", "w") as fh:
        for c, seq in ref.items():
            s = "".join(B[seq])
            if c == "chr2":
                s = s[:10_000] + s[10_000:30_000].lower() + s[30_000:]
            fh.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    subprocess.run([samtools, "faidx", str(d / "ref.fa")], check=True)
    snp_pos = {c: np.sort(rng.choice(L, L // 50, replace=False)) for c, L in contigs.items()}
    freq = {c: rng.uniform(0.1, 0.9, len(snp_pos[c])) for c in contigs}
    alt = {c: (ref[c][snp_pos[c]] + rng.integers(1, 4, len(snp_pos[c]))) % 4 for c in contigs}
    bams = []
    for si, k in enumerate(BAM_PLOIDY):
        sam = ["@HD\tVN:1.6"] + [f"@SQ\tSN:{c}\tLN:{L}" for c, L in contigs.items()]
        for c, L in contigs.items():
            haps = np.tile(ref[c], (k, 1))
            carriers = rng.random((k, len(snp_pos[c]))) < freq[c]
            for h in range(k):
                haps[h, snp_pos[c][carriers[h]]] = alt[c][carriers[h]]
            n_reads = L * k * 4 // 100
            for r, (st, h) in enumerate(zip(rng.integers(0, L - 100, n_reads), rng.integers(0, k, n_reads))):
                seq = "".join(B[haps[h, st:st + 100]])
                sam.append(f"{c}_{r}\t0\t{c}\t{st + 1}\t60\t100M\t*\t0\t0\t{seq}\t{'I' * 100}")
        sam_path = d / f"sample{si + 1}.sam"
        sam_path.write_text("\n".join(sam) + "\n")
        bam = str(d / f"sample{si + 1}.bam")
        subprocess.run([samtools, "sort", "-o", bam, str(sam_path)], check=True, capture_output=True)
        subprocess.run([samtools, "index", bam], check=True)
        bams.append(bam)
    (d / "bams.txt").write_text("\n".join(bams) + "\n")
    return str(d / "bams.txt"), str(d / "ref.fa")

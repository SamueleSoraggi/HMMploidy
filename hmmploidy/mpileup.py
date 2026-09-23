"""Fast mpileup parsing and genotype-likelihood computation.

This module is a from-scratch, vectorised rewrite of the logic that used to
live in ``Genotype_Likelihoods.py`` / ``generics.py``. The output format
(``.genolikes`` columns) is unchanged so existing downstream tooling keeps
working, but the implementation avoids the O(n^2) character-by-character
string building of the original script:

* :func:`decode_pileup_field` decodes one mpileup "read bases" field in a
  single left-to-right pass using list buffers instead of repeated string
  concatenation / ``str.find`` scans.
* :func:`genotype_log_likelihoods` computes, for every requested ploidy
  level, the log-likelihoods of all ``N+1`` possible dosages using numpy
  array operations across all reads at once instead of a Python loop per
  read.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

ALLELES = ["A", "C", "G", "T"]
PHRED_OFFSET = 33
_DIGITS = set("0123456789")


def decode_pileup_field(bases: str, quals: str, ref: str) -> Tuple[str, str]:
    """Decode one mpileup 'read bases' field into plain ACGT bases.

    Equivalent to the original ``convertSyms`` + 'X' stripping + quality
    filter pipeline, but runs in O(len(bases)) time.
    """
    out_bases: List[str] = []
    out_quals: List[str] = []
    i = 0
    n = len(bases)
    qi = 0  # index into the quality string (one char per read)
    while i < n:
        c = bases[i]
        if c == "." or c == ",":
            out_bases.append(ref)
            out_quals.append(quals[qi])
            i += 1
            qi += 1
        elif c in "ACGTacgt":
            out_bases.append(c.upper())
            out_quals.append(quals[qi])
            i += 1
            qi += 1
        elif c == "^":
            i += 2  # skip caret + mapping-quality char, no read consumed
        elif c == "$":
            i += 1  # end of read marker, no read consumed
        elif c in "*Nn><":
            # deleted base / undefined / ref skip: consumes one quality
            # slot but is dropped from the output (matches original 'X'
            # removal behaviour)
            i += 1
            qi += 1
        elif c == "+" or c == "-":
            j = i + 1
            while j < n and bases[j] in _DIGITS:
                j += 1
            length = int(bases[i + 1 : j])
            i = j + length
        else:
            i += 1
    return "".join(out_bases), "".join(out_quals)


def filter_by_quality(bases: str, quals: str, min_quality: int) -> Tuple[str, str]:
    if min_quality <= 0:
        return bases, quals
    keep_b: List[str] = []
    keep_q: List[str] = []
    for b, q in zip(bases, quals):
        if (ord(q) - PHRED_OFFSET) > min_quality:
            keep_b.append(b)
            keep_q.append(q)
    return "".join(keep_b), "".join(keep_q)


def base_counts(bases: str) -> np.ndarray:
    """Return counts of A,C,G,T in `bases` as a length-4 numpy array."""
    counts = np.zeros(4, dtype=np.int64)
    for i, a in enumerate(ALLELES):
        counts[i] = bases.count(a)
    return counts


def error_probs(quals: str) -> np.ndarray:
    """Per-read error probability implied by phred quality characters."""
    q = np.frombuffer(quals.encode("ascii"), dtype=np.uint8).astype(np.float64)
    return 10.0 ** (-(q - PHRED_OFFSET) / 10.0)


def haploid_log_likelihoods(bases: str, quals: str) -> List[float]:
    """Per-allele haploid genotype log-likelihoods (+ a 'triallelic ref'
    slot), matching the original ``calcGenoLogLike1``."""
    p = error_probs(quals)
    base_arr = np.frombuffer(bases.encode("ascii"), dtype=np.uint8)
    log_likes = [0.0, 0.0, 0.0, 0.0, 0.0]
    for j, allele in enumerate(ALLELES):
        match = base_arr == ord(allele)
        sublike = np.where(match, 1.0 - p, p / 3.0)
        log_likes[j] = float(np.sum(np.log(sublike)))
        if j == 0:
            log_likes[4] = float(np.sum(np.log(p / 3.0)))
    return log_likes


def genotype_log_likelihoods(
    ploidy: int, bases: str, quals: str, major: int, minor: int
) -> List[float]:
    """Log-likelihoods of the ``ploidy + 1`` possible dosages of the minor
    allele (0..ploidy copies), vectorised over all reads at a site.

    Equivalent to the original ``calcGenoLogLikeN_MajorMinor`` but without
    the per-read Python loop: the dosage loop only runs ``ploidy + 1``
    times (<= 9 for the ploidies supported here) while all reads are
    processed together with numpy.
    """
    if len(bases) == 0:
        return [0.0] * (ploidy + 1)

    p = error_probs(quals)
    base_arr = np.frombuffer(bases.encode("ascii"), dtype=np.uint8)
    is_major = base_arr == ord(ALLELES[major])
    is_minor = base_arr == ord(ALLELES[minor])

    match_prob = 1.0 - p
    mismatch_prob = p / 3.0

    log_likes = [0.0] * (ploidy + 1)
    # k = number of major-allele copies among the `ploidy` gene copies;
    # order matches itertools.combinations_with_replacement([major, minor], N)
    for k in range(ploidy, -1, -1):
        frac_major = k / ploidy
        frac_minor = (ploidy - k) / ploidy
        sublike = np.where(
            is_major,
            frac_major * match_prob + frac_minor * mismatch_prob,
            np.where(
                is_minor,
                frac_major * mismatch_prob + frac_minor * match_prob,
                frac_major * mismatch_prob + frac_minor * mismatch_prob,
            ),
        )
        log_likes[ploidy - k] = float(np.sum(np.log(sublike)))
    return log_likes


@dataclass
class PileupSite:
    chrom: str
    position: int
    reference: str
    sample_bases: List[str]
    sample_quals: List[str]


def parse_pileup_line(line: str, nsams: int) -> PileupSite:
    fields = line.rstrip("\n").split("\t")
    chrom, pos, ref = fields[0], int(fields[1]), fields[2]
    sample_bases = []
    sample_quals = []
    for n in range(nsams):
        base_field = fields[3 + n * 3 + 1]
        qual_field = fields[3 + n * 3 + 2]
        sample_bases.append(base_field)
        sample_quals.append(qual_field)
    return PileupSite(chrom, pos, ref, sample_bases, sample_quals)


def parse_inbreeding(spec: str, nsams: int) -> List[float]:
    values: List[float] = []
    for token in spec.split(","):
        parts = token.split("x")
        if len(parts) == 2:
            values.extend([float(parts[0])] * int(parts[1]))
        else:
            values.append(float(parts[0]))
    if len(values) != nsams:
        raise ValueError(
            f"Inbreeding coefficients ({len(values)}) do not match number of samples ({nsams})"
        )
    return values

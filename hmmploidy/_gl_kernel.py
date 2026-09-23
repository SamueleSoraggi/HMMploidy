"""numba kernel computing one ``.genolikes`` site from one raw mpileup line.

This is the compiled counterpart of :func:`hmmploidy.genotype_likelihoods.process_site`
(which is kept as the readable reference implementation): the line is
decoded, quality-filtered, site-filtered and turned into genotype
likelihoods for every requested ploidy in a single pass over its bytes, so
no per-read or per-sample Python objects are created.
"""
from __future__ import annotations

import math

import numpy as np

from .hmm.common import njit

# byte codes
_TAB = 9
_DOT = 46
_COMMA = 44
_CARET = 94
_DOLLAR = 36
_PLUS = 43
_MINUS = 45
_STAR = 42
_GT = 62
_LT = 60


@njit(cache=True)
def _allele_code(c):
    # A C G T (either case) -> 0..3, anything else -> 4
    if c == 65 or c == 97:
        return 0
    if c == 67 or c == 99:
        return 1
    if c == 71 or c == 103:
        return 2
    if c == 84 or c == 116:
        return 3
    return 4


@njit(cache=True)
def seed_kernel_rng(seed):
    np.random.seed(seed)


@njit(cache=True)
def process_line(
    line,
    nsams,
    ploidies,
    min_quality_score,
    min_non_major_freq,
    max_minor2_freq,
    max_minor3_freq,
    min_global_depth,
    min_ind_depth,
    downsampling,
):
    """Process one mpileup line (uint8 array without trailing newline).

    Returns ``(keep, major, minor, depth, major_count, minor_count, gls)``
    where the per-sample arrays have length ``nsams`` and ``gls`` has shape
    ``(nsams, sum(ploidies + 1))``.
    """
    n_gl = 0
    for p in ploidies:
        n_gl += p + 1
    depth = np.zeros(nsams, dtype=np.int64)
    major_count = np.zeros(nsams, dtype=np.int64)
    minor_count = np.zeros(nsams, dtype=np.int64)
    gls = np.zeros((nsams, n_gl))

    # locate tab-separated fields
    L = line.shape[0]
    n_fields = 3 + 3 * nsams
    starts = np.zeros(n_fields, dtype=np.int64)
    ends = np.zeros(n_fields, dtype=np.int64)
    f = 0
    starts[0] = 0
    for i in range(L):
        if line[i] == _TAB:
            ends[f] = i
            f += 1
            if f >= n_fields:
                break
            starts[f] = i + 1
    if f < n_fields:
        ends[f] = L
    if f < n_fields - 1:
        return False, 0, 0, depth, major_count, minor_count, gls

    ref = _allele_code(line[starts[2]])

    # decode every sample's reads into allele codes + error probabilities
    max_reads = 0
    for n in range(nsams):
        q_len = ends[3 + 3 * n + 2] - starts[3 + 3 * n + 2]
        if q_len > max_reads:
            max_reads = q_len
    codes = np.full((nsams, max(max_reads, 1)), 4, dtype=np.int64)
    perr = np.zeros((nsams, max(max_reads, 1)))

    for n in range(nsams):
        b0 = starts[3 + 3 * n + 1]
        b1 = ends[3 + 3 * n + 1]
        q0 = starts[3 + 3 * n + 2]
        q1 = ends[3 + 3 * n + 2]
        i = b0
        qi = q0
        k = 0
        while i < b1:
            c = line[i]
            code = -1
            if c == _DOT or c == _COMMA:
                code = ref
                i += 1
            elif _allele_code(c) < 4:
                code = _allele_code(c)
                i += 1
            elif c == _CARET:
                i += 2
                continue
            elif c == _DOLLAR:
                i += 1
                continue
            elif c == _STAR or c == 78 or c == 110 or c == _GT or c == _LT:
                # deletion / N / reference skip: consumes a quality, dropped
                i += 1
                qi += 1
                continue
            elif c == _PLUS or c == _MINUS:
                j = i + 1
                length = 0
                while j < b1 and line[j] >= 48 and line[j] <= 57:
                    length = length * 10 + (line[j] - 48)
                    j += 1
                i = j + length
                continue
            else:
                i += 1
                continue
            if qi >= q1:
                break
            q = line[qi] - 33
            qi += 1
            if min_quality_score > 0 and q <= min_quality_score:
                continue
            codes[n, k] = code
            perr[n, k] = 10.0 ** (-q / 10.0)
            k += 1
        depth[n] = k

    global_depth = 0
    min_depth = depth[0] if nsams > 0 else 0
    for n in range(nsams):
        global_depth += depth[n]
        if depth[n] < min_depth:
            min_depth = depth[n]
    if not (global_depth > min_global_depth and min_depth > min_ind_depth):
        return False, 0, 0, depth, major_count, minor_count, gls

    pooled = np.zeros(5, dtype=np.int64)
    haploid = np.zeros(4)
    for n in range(nsams):
        for r in range(depth[n]):
            c = codes[n, r]
            pooled[c] += 1
            p = perr[n, r]
            for a in range(4):
                if c == a:
                    haploid[a] += math.log(1.0 - p)
                else:
                    haploid[a] += math.log(p / 3.0)

    max_count = max(pooled[0], pooled[1], pooled[2], pooled[3])
    if (global_depth - max_count) / global_depth <= min_non_major_freq:
        return False, 0, 0, depth, major_count, minor_count, gls

    # descending order of the haploid likelihoods rounded to 1e-8, ties to the
    # lower allele index (same rule as the reference process_site)
    for a in range(4):
        haploid[a] = math.floor(haploid[a] * 1e8 + 0.5)
    order = np.zeros(4, dtype=np.int64)
    used = np.zeros(4, dtype=np.bool_)
    for r in range(4):
        best = -1
        for a in range(4):
            if not used[a] and (best < 0 or haploid[a] > haploid[best]):
                best = a
        order[r] = best
        used[best] = True
    major, minor, minor2, minor3 = order[0], order[1], order[2], order[3]
    if not (
        pooled[minor2] / global_depth < max_minor2_freq
        and pooled[minor3] / global_depth < max_minor3_freq
    ):
        return False, 0, 0, depth, major_count, minor_count, gls

    for n in range(nsams):
        d = depth[n]
        if downsampling < 1.0:
            keep_n = int(math.ceil(d * downsampling))
            if keep_n < d:
                # partial Fisher-Yates: move a random subset to the front
                for r in range(keep_n):
                    j = r + int(np.random.random() * (d - r))
                    tc = codes[n, r]
                    codes[n, r] = codes[n, j]
                    codes[n, j] = tc
                    tp = perr[n, r]
                    perr[n, r] = perr[n, j]
                    perr[n, j] = tp
                d = keep_n
                depth[n] = d
        for r in range(d):
            if codes[n, r] == major:
                major_count[n] += 1
            elif codes[n, r] == minor:
                minor_count[n] += 1

        col = 0
        for p in ploidies:
            if d > 0:
                # dosage of the minor allele = ploidy - (#major copies)
                for kmaj in range(p, -1, -1):
                    fmaj = kmaj / p
                    fmin = (p - kmaj) / p
                    s = 0.0
                    for r in range(d):
                        c = codes[n, r]
                        match = 1.0 - perr[n, r]
                        mis = perr[n, r] / 3.0
                        if c == major:
                            s += math.log(fmaj * match + fmin * mis)
                        elif c == minor:
                            s += math.log(fmaj * mis + fmin * match)
                        else:
                            s += math.log(fmaj * mis + fmin * mis)
                    gls[n, col + p - kmaj] = s
            col += p + 1

    return True, major, minor, depth, major_count, minor_count, gls

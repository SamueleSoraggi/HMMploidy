"""Small numerical helpers shared across the package.

This module replaces the original ``generics.py``. The read/base parsing
helpers used by :mod:`hmmploidy.genotype_likelihoods` were rewritten to
avoid quadratic-time string operations (the original implementation
repeatedly called ``str.find``/slicing in a loop for every filtered
character); see :mod:`hmmploidy.mpileup` for the new implementation.
"""
from __future__ import annotations

import math
from itertools import combinations_with_replacement
from typing import Iterable, List, Sequence

ALLELES = ["A", "C", "G", "T"]


def exp_or_zero(x: float) -> float:
    return 0.0 if x == 0 else math.exp(x)


def log_or_zero(x: float) -> float:
    return 0.0 if x == 0 else math.log(x)


def dist(ploidies: Sequence[int], max_ploidy: int = 8) -> List[float]:
    """Empirical distribution of inferred ploidies over 1..max_ploidy."""
    n = len(ploidies)
    return [ploidies.count(i) / n for i in range(1, max_ploidy + 1)]


def combinations_with_rep(iterable: Iterable, r: int):
    """Kept for backward compatibility; use itertools directly in new code."""
    return combinations_with_replacement(iterable, r)

"""Efficient Poisson-emission HMM, fit by Baum-Welch EM with closed-form
M-step updates.

This is the new "``em``" model choice requested in addition to the
original Negative-Binomial method: by using a Poisson emission for depth
instead of a Negative Binomial, the M-step for the emission parameters has
a closed-form solution (the weighted mean), so no Newton-Raphson
inner-loop is needed. This is both simpler and noticeably faster to
converge than :mod:`hmmploidy.hmm.negbinom`, at the cost of not modelling
depth over-dispersion explicitly (over-dispersion is instead partly
absorbed by the window-averaging of depth, which already reduces variance
relative to per-site read counts).
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.special import gammaln

from .common import forward_backward, transition_xi_sum, viterbi


def _log_poisson_density(count: np.ndarray, lam: np.ndarray) -> np.ndarray:
    lam = np.clip(lam, 1e-6, None)
    count = count[:, None]
    return count * np.log(lam)[None, :] - lam[None, :] - gammaln(count + 1)


def fit_poisson_hmm(
    depth: np.ndarray,
    log_emission: np.ndarray,
    max_iter: int = 500,
    tol: float = 1e-6,
) -> Dict[str, np.ndarray]:
    count = np.asarray(depth, dtype=float)
    T, N = log_emission.shape

    positive = count[count > 0]
    if positive.size == 0:
        positive = np.array([1.0])
    qs = np.linspace(0, 1, N + 1)[1:]
    lam = np.clip(np.quantile(positive, qs), 1e-3, None)
    delta = np.full(N, 1.0 / N)
    trans = np.full((N, N), 1.0 / N)

    loglik_hist = []
    for it in range(max_iter):
        log_dens = _log_poisson_density(count, lam) + log_emission
        dens = np.exp(np.clip(log_dens, -700, 700))

        loglik, gamma, forward, backward, scale = forward_backward(dens, delta, trans)
        loglik_hist.append(loglik)

        new_delta = gamma[0].copy()
        new_trans = transition_xi_sum(forward, backward, dens, trans, scale)

        eq_count = gamma.sum(axis=0)
        eq_count_safe = np.where(eq_count > 0, eq_count, 1e-300)
        new_lam = np.clip((count @ gamma) / eq_count_safe, 1e-6, None)

        lam, delta, trans = new_lam, new_delta, new_trans
        if it > 0 and abs((loglik_hist[-1] - loglik_hist[-2]) / (loglik_hist[-2] or 1e-300)) < tol:
            break

    log_dens = _log_poisson_density(count, lam) + log_emission
    path, _nu = viterbi(log_dens, delta, trans)
    _, gamma, *_ = forward_backward(np.exp(np.clip(log_dens, -700, 700)), delta, trans)

    return {
        "model": "em",
        "mu": lam,
        "delta": delta,
        "trans": trans,
        "loglik": np.array(loglik_hist),
        "states": path + 1,
        "posterior": gamma,
    }

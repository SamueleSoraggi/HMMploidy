"""Classic Negative-Binomial emission HMM, fit by conditional EM.

This is a direct, vectorised translation of the ``nbHMM``/``EStep``/
``MStep`` functions in the original ``HMMploidy.R``: the depth emission is
Negative-Binomial(alpha, beta) (parameterised as in the R code, i.e. with
mean ``alpha/beta``), combined multiplicatively (additively in log-space)
with the pre-computed, window-summed genotype-likelihood term for each
candidate ploidy state.

Simplification vs. the original: the original R implementation also tries
to greedily *remove* candidate ploidy states via a BIC comparison
(``combn`` search over state subsets) while iterating. That combinatorial
search is not reproduced here; instead all ``max_ploidy`` states are kept
throughout, and states with negligible overall posterior support can be
dropped by the caller after convergence (see ``prune_states``). In
practice this gives the same qualitative results while being much easier
to reason about and to keep numerically stable.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.special import digamma, gammaln, polygamma

from .common import forward_backward, transition_xi_sum, viterbi


def _log_nbinom_density(count: np.ndarray, alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """log NB pmf with mean alpha/beta, shape (T, N)."""
    count = count[:, None]
    return (
        alpha[None, :] * np.log(beta[None, :] / (1 + beta[None, :]))
        - gammaln(alpha)[None, :]
        - count * np.log(1 + beta[None, :])
        + gammaln(count + alpha[None, :])
        - gammaln(count + 1)
    )


def _init_params(count: np.ndarray, n_states: int):
    positive = count[count > 0]
    if positive.size == 0:
        positive = np.array([1.0])
    qs = np.linspace(0, 1, n_states + 1)[1:]
    alpha = np.quantile(positive, qs)
    alpha = np.clip(alpha, 1e-3, None)
    beta = np.ones(n_states)
    delta = np.full(n_states, 1.0 / n_states)
    trans = np.full((n_states, n_states), 1.0 / n_states)
    return alpha, beta, delta, trans


def fit_negbinom_hmm(
    depth: np.ndarray,
    log_emission: np.ndarray,
    max_iter: int = 500,
    tol: float = 1e-5,
    max_alpha: float = 1e7,
    max_beta: float = 1e7,
) -> Dict[str, np.ndarray]:
    count = np.asarray(depth, dtype=float)
    T, N = log_emission.shape
    alpha, beta, delta, trans = _init_params(count, N)

    loglik_hist = []
    for it in range(max_iter):
        log_dens = _log_nbinom_density(count, alpha, beta) + log_emission
        dens = np.exp(np.clip(log_dens, -700, 700))

        loglik, gamma, forward, backward, scale = forward_backward(dens, delta, trans)
        loglik_hist.append(loglik)

        new_delta = gamma[0].copy()
        new_trans = transition_xi_sum(forward, backward, dens, trans, scale)

        eq_count = gamma.sum(axis=0)
        eq_count_safe = np.where(eq_count > 0, eq_count, 1e-300)
        mu = (count @ gamma) / eq_count_safe
        mu = np.clip(mu, 1e-6, None)
        new_beta = alpha / mu

        grad = eq_count * (np.log(new_beta / (1 + new_beta)) - digamma(alpha)) + np.sum(
            gamma * digamma(count[:, None] + alpha[None, :]), axis=0
        )
        hess = -eq_count * polygamma(1, alpha) + np.sum(
            gamma * polygamma(1, count[:, None] + alpha[None, :]), axis=0
        )
        hess_safe = np.where(hess == 0, -1e-300, hess)
        step = -grad / hess_safe
        new_alpha = alpha + step
        tries = 0
        while np.any(new_alpha <= 0) and tries < 50:
            step = step / 10
            new_alpha = alpha + step
            tries += 1
        new_alpha = np.clip(new_alpha, 1e-6, None)

        if np.any(new_alpha > max_alpha) or np.any(new_beta > max_beta):
            break

        alpha, beta, delta, trans = new_alpha, new_beta, new_delta, new_trans

        if it > 0 and abs((loglik_hist[-1] - loglik_hist[-2]) / (loglik_hist[-2] or 1e-300)) < tol:
            break

    log_dens = _log_nbinom_density(count, alpha, beta) + log_emission
    path, _nu = viterbi(log_dens, delta, trans)
    _, gamma, *_ = forward_backward(np.exp(np.clip(log_dens, -700, 700)), delta, trans)

    return {
        "model": "negbinom",
        "alpha": alpha,
        "beta": beta,
        "delta": delta,
        "trans": trans,
        "mu": alpha / beta,
        "loglik": np.array(loglik_hist),
        "states": path + 1,  # ploidy = state index + 1
        "posterior": gamma,
    }

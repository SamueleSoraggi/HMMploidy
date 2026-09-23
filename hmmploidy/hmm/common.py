"""Shared HMM utilities: scaled forward-backward, Viterbi decoding.

These are used by every emission model in :mod:`hmmploidy.hmm` so that only
the emission-density computation differs between models.

The recursions are compiled with numba when it is installed (a large
speed-up, since they are sequential loops over windows that numpy cannot
vectorise); otherwise the same code runs as plain Python/numpy.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np

try:
    from numba import njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _forward_backward_kernel(dens, delta, trans):
    T, N = dens.shape
    forward = np.zeros((T, N))
    backward = np.zeros((T, N))
    scale = np.zeros(T)

    s = 0.0
    for j in range(N):
        forward[0, j] = delta[j] * dens[0, j]
        s += forward[0, j]
    scale[0] = s if s > 0 else 1e-300
    for j in range(N):
        forward[0, j] /= scale[0]

    for t in range(1, T):
        s = 0.0
        for j in range(N):
            acc = 0.0
            for i in range(N):
                acc += forward[t - 1, i] * trans[i, j]
            forward[t, j] = acc * dens[t, j]
            s += forward[t, j]
        scale[t] = s if s > 0 else 1e-300
        for j in range(N):
            forward[t, j] /= scale[t]

    for j in range(N):
        backward[T - 1, j] = 1.0
    for t in range(T - 2, -1, -1):
        for i in range(N):
            acc = 0.0
            for j in range(N):
                acc += trans[i, j] * dens[t + 1, j] * backward[t + 1, j]
            backward[t, i] = acc / scale[t + 1]
    return forward, backward, scale


def forward_backward(
    dens: np.ndarray, delta: np.ndarray, trans: np.ndarray
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Scaled forward-backward algorithm.

    Parameters
    ----------
    dens: (T, N) array of (non-log) emission densities per state.
    delta: (N,) initial state distribution.
    trans: (N, N) transition matrix.

    Returns
    -------
    loglik, gamma (T,N posteriors), forward (T,N), backward (T,N), scale (T,)
    """
    forward, backward, scale = _forward_backward_kernel(
        np.ascontiguousarray(dens, dtype=np.float64),
        np.ascontiguousarray(delta, dtype=np.float64),
        np.ascontiguousarray(trans, dtype=np.float64),
    )
    loglik = float(np.sum(np.log(scale)))
    gamma = forward * backward
    gamma_sum = gamma.sum(axis=1, keepdims=True)
    gamma_sum[gamma_sum == 0] = 1e-300
    gamma /= gamma_sum
    return loglik, gamma, forward, backward, scale


def forward_backward_log(
    log_dens: np.ndarray, delta: np.ndarray, trans: np.ndarray
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Forward-backward on log emission densities of arbitrary magnitude.

    Each row is shifted by its maximum before exponentiating (and the shift
    added back to the log-likelihood), so emissions that are sums of many
    per-site log-likelihoods do not underflow.

    Returns loglik, gamma (T,N posteriors) and the (N,N) matrix of expected
    transition counts sum_t P(s_t=i, s_{t+1}=j | data).
    """
    row_max = log_dens.max(axis=1, keepdims=True)
    dens = np.exp(log_dens - row_max)
    loglik, gamma, forward, backward, scale = forward_backward(dens, delta, trans)
    loglik += float(row_max.sum())
    xi = trans * (forward[:-1].T @ (dens[1:] * backward[1:] / scale[1:, None]))
    return loglik, gamma, xi


@njit(cache=True)
def _viterbi_kernel(log_dens, log_delta, log_trans):
    T, N = log_dens.shape
    nu = np.zeros((T, N))
    psi = np.zeros((T, N), dtype=np.int64)
    for j in range(N):
        nu[0, j] = log_delta[j] + log_dens[0, j]
    for t in range(1, T):
        for j in range(N):
            best = -np.inf
            arg = 0
            for i in range(N):
                v = nu[t - 1, i] + log_trans[i, j]
                if v > best:
                    best = v
                    arg = i
            psi[t, j] = arg
            nu[t, j] = best + log_dens[t, j]
    path = np.zeros(T, dtype=np.int64)
    path[T - 1] = np.argmax(nu[T - 1])
    for t in range(T - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]
    return path, nu


def viterbi(log_dens: np.ndarray, delta: np.ndarray, trans: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Standard log-space Viterbi decoding.

    Returns the most likely state path (0-indexed) and the DP table.
    """
    log_delta = np.log(np.clip(delta, 1e-300, None))
    log_trans = np.log(np.clip(trans, 1e-300, None))
    return _viterbi_kernel(
        np.ascontiguousarray(log_dens, dtype=np.float64),
        np.ascontiguousarray(log_delta, dtype=np.float64),
        np.ascontiguousarray(log_trans, dtype=np.float64),
    )


def transition_xi_sum(
    forward: np.ndarray,
    backward: np.ndarray,
    dens: np.ndarray,
    trans: np.ndarray,
    scale: np.ndarray,
) -> np.ndarray:
    """Row-normalised sum over t of xi(t) = P(state_t=i, state_{t+1}=j |
    data), i.e. the transition-matrix M-step update."""
    xi_sum = trans * (forward[:-1].T @ (dens[1:] * backward[1:] / scale[1:, None]))
    row_sum = xi_sum.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1e-300
    return xi_sum / row_sum

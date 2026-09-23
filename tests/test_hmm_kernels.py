"""The compiled forward-backward / Viterbi recursions against brute-force
enumeration of every state path of a small HMM."""
import itertools

import numpy as np
import pytest
from scipy.special import logsumexp

from hmmploidy.hmm.common import forward_backward, forward_backward_log, transition_xi_sum, viterbi


def brute_force(log_dens, delta, trans):
    T, N = log_dens.shape
    paths = list(itertools.product(range(N), repeat=T))
    lp = np.array([
        np.log(delta[p[0]]) + log_dens[0, p[0]]
        + sum(np.log(trans[p[t - 1], p[t]]) + log_dens[t, p[t]] for t in range(1, T))
        for p in paths
    ])
    loglik = logsumexp(lp)
    w = np.exp(lp - loglik)
    gamma = np.zeros((T, N))
    xi = np.zeros((N, N))
    for p, wi in zip(paths, w):
        for t in range(T):
            gamma[t, p[t]] += wi
        for t in range(T - 1):
            xi[p[t], p[t + 1]] += wi
    return loglik, gamma, xi, np.array(paths[int(np.argmax(lp))])


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_forward_backward_and_viterbi_match_brute_force(seed):
    rng = np.random.default_rng(seed)
    T, N = 6, 3
    # emissions far below exp() range, as with window-summed likelihoods
    log_dens = rng.normal(size=(T, N)) * 20 - 800
    delta = rng.dirichlet(np.ones(N))
    trans = rng.dirichlet(np.ones(N), size=N)

    ll_ref, gamma_ref, xi_ref, path_ref = brute_force(log_dens, delta, trans)
    ll, gamma, xi = forward_backward_log(log_dens, delta, trans)

    assert ll == pytest.approx(ll_ref)
    np.testing.assert_allclose(gamma, gamma_ref, atol=1e-10)
    np.testing.assert_allclose(xi, xi_ref, atol=1e-10)
    np.testing.assert_array_equal(viterbi(log_dens, delta, trans)[0], path_ref)


def test_transition_update_is_row_normalised_expected_counts():
    rng = np.random.default_rng(3)
    T, N = 6, 3
    log_dens = rng.normal(size=(T, N))
    delta = rng.dirichlet(np.ones(N))
    trans = rng.dirichlet(np.ones(N), size=N)
    _, _, xi_ref, _ = brute_force(log_dens, delta, trans)

    _, _, forward, backward, scale = forward_backward(np.exp(log_dens), delta, trans)
    new_trans = transition_xi_sum(forward, backward, np.exp(log_dens), trans, scale)
    np.testing.assert_allclose(new_trans, xi_ref / xi_ref.sum(axis=1, keepdims=True), atol=1e-10)

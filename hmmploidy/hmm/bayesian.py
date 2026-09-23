"""Hierarchical Bayesian HMM for ploidy inference, fit with PyMC.

Depth emissions are modelled as Negative-Binomial with a mean that is
constrained to be an integer multiple of a single shared "haploid depth"
parameter (``mu_k = haploid_depth * k`` for candidate ploidy ``k``), which
encodes the same biological assumption as the other two models (depth
scales with copy number) while sharing statistical strength across states
and giving a principled, checkable prior. The transition matrix rows and
the initial-state distribution get Dirichlet priors. The genotype
likelihoods computed from the data (``log_emission``) enter the model as a
fixed additive offset to the depth log-density, exactly as in the other
two models.

Because PyMC has no built-in discrete-HMM likelihood, the forward
algorithm is implemented as a ``pytensor.scan`` and exposed to the sampler
through a ``pm.Potential`` (a common pattern for marginalising discrete
HMM states in a fully continuous, NUTS-friendly model). After sampling,
the posterior-mean parameters are plugged into the same Viterbi decoder
used by the other two models to obtain a single best ploidy path.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from .common import viterbi


def _forward_loglik(log_dens, log_trans, log_delta):
    import pytensor
    import pytensor.tensor as pt

    def step(log_e_t, log_alpha_prev, log_trans_ns):
        m = log_alpha_prev[:, None] + log_trans_ns  # (from, to)
        return pt.logsumexp(m, axis=0) + log_e_t

    init = log_delta + log_dens[0]

    result, _ = pytensor.scan(
        fn=step,
        sequences=log_dens[1:],
        outputs_info=init,
        non_sequences=[log_trans],
    )
    log_alpha_T = result[-1]
    return pt.logsumexp(log_alpha_T)


def fit_bayesian_hmm(
    depth: np.ndarray,
    log_emission: np.ndarray,
    draws: int = 500,
    tune: int = 500,
    chains: int = 2,
    target_accept: float = 0.9,
    seed: int = 1,
) -> Dict[str, np.ndarray]:
    try:
        import pymc as pm
        import pytensor.tensor as pt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "The 'bayesian' model requires pymc and pytensor. Install with "
            "`conda env create -f environment.yml` or `pip install pymc`."
        ) from exc

    count = np.asarray(depth, dtype=float)
    T, N = log_emission.shape
    log_emission_const = pt.as_tensor_variable(log_emission)
    mean_depth_guess = max(float(count.mean()), 1e-3)

    with pm.Model() as model:
        haploid_depth = pm.Gamma("haploid_depth", alpha=2.0, beta=2.0 / mean_depth_guess)
        dispersion = pm.Gamma("dispersion", alpha=2.0, beta=1.0)

        ploidies = pt.arange(1, N + 1, dtype="float64")
        mu = haploid_depth * ploidies  # (N,)
        p_nb = dispersion / (dispersion + mu)

        count_col = pt.as_tensor_variable(count)[:, None]
        log_nb = (
            pt.gammaln(count_col + dispersion)
            - pt.gammaln(dispersion)
            - pt.gammaln(count_col + 1)
            + dispersion * pt.log(p_nb)[None, :]
            + count_col * pt.log(1 - p_nb)[None, :]
        )
        log_dens = log_nb + log_emission_const  # (T, N)

        trans_conc = np.ones(N) + np.eye(N) * (5.0 * N)
        trans = pm.Dirichlet("trans", a=trans_conc, shape=(N, N))
        delta = pm.Dirichlet("delta", a=np.ones(N))

        loglik = _forward_loglik(log_dens, pt.log(trans + 1e-300), pt.log(delta + 1e-300))
        pm.Potential("hmm_loglik", loglik)

        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            random_seed=seed,
            target_accept=target_accept,
            progressbar=True,
        )

    post = idata.posterior
    haploid_depth_hat = float(post["haploid_depth"].mean())
    dispersion_hat = float(post["dispersion"].mean())
    trans_hat = post["trans"].mean(dim=("chain", "draw")).values
    delta_hat = post["delta"].mean(dim=("chain", "draw")).values

    mu_hat = haploid_depth_hat * np.arange(1, N + 1)
    p_hat = dispersion_hat / (dispersion_hat + mu_hat)
    from scipy.special import gammaln

    log_nb_hat = (
        gammaln(count[:, None] + dispersion_hat)
        - gammaln(dispersion_hat)
        - gammaln(count[:, None] + 1)
        + dispersion_hat * np.log(p_hat)[None, :]
        + count[:, None] * np.log(1 - p_hat)[None, :]
    )
    log_dens_hat = log_nb_hat + log_emission
    path, _nu = viterbi(log_dens_hat, delta_hat, trans_hat)

    return {
        "model": "bayesian",
        "haploid_depth": haploid_depth_hat,
        "dispersion": dispersion_hat,
        "mu": mu_hat,
        "delta": delta_hat,
        "trans": trans_hat,
        "states": path + 1,
        "idata": idata,
    }

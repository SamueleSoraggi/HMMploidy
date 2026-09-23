"""Overdispersion-aware HMM: Beta-Binomial allele counts + Gamma depth.

The other models take the genotype-likelihood term as fixed, pre-computed
input (a Binomial model of read sampling, summed over all SNPs of a window).
Real data are overdispersed relative to that model -- reference bias,
mapping artefacts and PCR duplicates spread the allele ratios of
heterozygous sites well beyond Binomial noise -- and because the per-SNP
terms of a window are summed, even mild misspecification turns into very
confident, wrong ploidy calls. This model instead computes its own
genotype likelihoods from the per-SNP allele counts:

* **Allele counts**: minor-allele reads at a SNP with dosage ``g`` out of
  ``k`` copies are Beta-Binomial with mean ``g/k`` (corrected for sequencing
  error) and intra-class correlation ``rho``; the dosage is marginalised
  under a Beta-Binomial prior: Binomial(k, population frequency) with the
  frequency, estimated leave-one-out from the *other* individuals,
  integrated over its Beta uncertainty (see
  :func:`hmmploidy.data.frequency_beta_prior`). With few individuals a
  point-estimate frequency biases the calls. ``rho`` is estimated by
  maximising the HMM marginal likelihood over a grid; ``rho = 0`` recovers
  the Binomial model.
* **Window depth**: Gamma with mean ``k * haploid_depth`` and a shared shape,
  i.e. depth is proportional to ploidy, with a variance that grows with the
  square of the mean (as regional coverage biases such as GC content do). A
  small outlier component absorbs windows with aberrant depth (repeats,
  collapsed duplications). The M-step is closed-form for the haploid depth
  and a 1-D Newton solve for the shape.
* **Transitions**: a Dirichlet ("sticky") prior adds ``sticky`` pseudo-counts
  to self-transitions. Without it, the model can use rapid switching
  between neighbouring ploidies (e.g. 4 <-> 5) to explain regional coverage
  noise when the allele data only weakly separate them.

Because depth only identifies ploidy up to the haploid-depth scale (2 copies
at depth h look like 4 copies at h/2), the fit at the selected ``rho`` is
restarted from every "whole genome at ploidy k" configuration and the
restart with the highest likelihood is kept.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
from scipy.special import betaln, digamma, gammaln, polygamma

from .common import forward_backward_log, viterbi

DEFAULT_RHO_GRID = (0.0, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3)


def snp_log_likelihoods(snps, n_states: int, rho: float) -> np.ndarray:
    """log P(allele counts | ploidy k), shape (n_snp, n_states), up to a
    constant (the binomial coefficient) that does not depend on k or rho."""
    m = snps.minor[:, None]
    n = snps.total[:, None]
    f = np.clip(snps.freq, 1e-6, 1 - 1e-6)[:, None]
    uncertain_f = snps.prior_a is not None
    if uncertain_f:
        pa, pb = snps.prior_a[:, None], snps.prior_b[:, None]
    err = snps.error_rate
    conc = (1.0 - rho) / rho if rho > 0 else None

    out = np.empty((len(snps.minor), n_states))
    for k in range(1, n_states + 1):
        g = np.arange(k + 1)[None, :]
        frac = g / k
        # P(read shows minor | read shows major or minor)
        p_minor = frac * (1 - 3 * err) + (1 - frac) * err
        p_major = (1 - frac) * (1 - 3 * err) + frac * err
        pi = p_minor / (p_minor + p_major)

        lik = m * np.log(pi) + (n - m) * np.log1p(-pi)
        if conc is not None and k > 1:
            het = slice(1, k)
            a, b = pi[:, het] * conc, (1 - pi[:, het]) * conc
            lik[:, het] = betaln(m + a, n - m + b) - betaln(a, b)

        log_prior = gammaln(k + 1) - gammaln(g + 1) - gammaln(k - g + 1)
        if uncertain_f:
            log_prior = log_prior + betaln(g + pa, k - g + pb) - betaln(pa, pb)
        else:
            log_prior = log_prior + g * np.log(f) + (k - g) * np.log1p(-f)
        x = lik + log_prior
        x_max = x.max(axis=1, keepdims=True)
        out[:, k - 1] = (x_max + np.log(np.exp(x - x_max).sum(axis=1, keepdims=True)))[:, 0]
    return out


def window_log_likelihoods(snps, n_windows: int, n_states: int, rho: float) -> np.ndarray:
    """Sum of :func:`snp_log_likelihoods` over the SNPs of each window."""
    per_snp = snp_log_likelihoods(snps, n_states, rho)
    out = np.zeros((n_windows, n_states))
    for k in range(n_states):
        out[:, k] = np.bincount(snps.window, weights=per_snp[:, k], minlength=n_windows)
    return out


def _gamma_logpdf(x: np.ndarray, mean: np.ndarray, shape: float) -> np.ndarray:
    x = x[:, None]
    mean = mean[None, :]
    return (
        shape * np.log(shape / mean)
        + (shape - 1) * np.log(x)
        - gammaln(shape)
        - shape * x / mean
    )


def _solve_gamma_shape(s: float, shape: float) -> float:
    """Solve log(a) - digamma(a) = s for the Gamma shape a (Newton on log a)."""
    if not np.isfinite(s) or s <= 0:
        return 1e4
    a = (3 - s + np.sqrt((s - 3) ** 2 + 24 * s)) / (12 * s)
    for _ in range(50):
        f = np.log(a) - digamma(a) - s
        df = 1 / a - polygamma(1, a)
        step = f / (df * a)
        a = a * np.exp(-np.clip(step, -2, 2))
        if abs(step) < 1e-10:
            break
    return float(np.clip(a, 1e-2, 1e4))


class _DepthModel:
    def __init__(self, depth: np.ndarray, n_states: int, outlier: float):
        self.x = np.clip(np.asarray(depth, dtype=float), 1e-3, None)
        self.k = np.arange(1, n_states + 1, dtype=float)
        self.outlier = outlier
        # outlier component: uniform over [0, max depth]
        self.log_outlier = np.log(outlier) - np.log(self.x.max() * 1.01) if outlier > 0 else None

    def init(self, states: np.ndarray):
        ratio = self.x / self.k[states]
        self.h = float(np.median(ratio))
        cv2 = np.var(ratio / self.h)
        self.shape = float(np.clip(1.0 / max(cv2, 1e-3), 1.0, 1e4))

    def logpdf(self):
        """Return (T, N) log density and the (T, N) inlier responsibility."""
        inlier = _gamma_logpdf(self.x, self.k * self.h, self.shape)
        if self.log_outlier is None:
            return inlier, None
        inlier = inlier + np.log1p(-self.outlier)
        total = np.logaddexp(inlier, self.log_outlier)
        return total, np.exp(inlier - total)

    def m_step(self, post: np.ndarray, resp: Optional[np.ndarray]):
        w = post if resp is None else post * resp
        wsum = w.sum()
        if wsum <= 0:
            return
        # closed-form haploid depth for a Gamma with mean k*h and fixed shape
        self.h = float((w * (self.x[:, None] / self.k[None, :])).sum() / wsum)
        z = self.x[:, None] / (self.k[None, :] * self.h)
        s = float((w * (z - np.log(z) - 1)).sum() / wsum)
        self.shape = _solve_gamma_shape(s, self.shape)


def _em(depth_model, gl_win, delta, trans, sticky, max_iter, tol):
    N = gl_win.shape[1]
    loglik_hist = []
    for _ in range(max_iter):
        log_depth, resp = depth_model.logpdf()
        ll, post, xi = forward_backward_log(log_depth + gl_win, delta, trans)
        loglik_hist.append(ll)

        delta = np.clip(post[0], 1e-6, None)
        delta /= delta.sum()
        counts = xi + sticky * np.eye(N) + 1e-6
        trans = counts / counts.sum(axis=1, keepdims=True)
        depth_model.m_step(post, resp)

        if len(loglik_hist) > 1 and abs(loglik_hist[-1] - loglik_hist[-2]) < tol * abs(loglik_hist[-2]):
            break
    log_depth, _ = depth_model.logpdf()
    ll, post, _ = forward_backward_log(log_depth + gl_win, delta, trans)
    loglik_hist.append(ll)
    # EM maximises the posterior under the sticky prior, so fits are
    # compared on that objective rather than on the raw likelihood
    objective = ll + sticky * float(np.log(np.diag(trans)).sum())
    return dict(
        ll=ll, objective=objective, delta=delta, trans=trans, post=post,
        log_dens=log_depth + gl_win, h=depth_model.h, shape=depth_model.shape,
        hist=loglik_hist,
    )


def _initial_transitions(N):
    trans = np.full((N, N), 0.1 / max(N - 1, 1))
    np.fill_diagonal(trans, 0.9)
    return np.full(N, 1.0 / N), trans


def fit_overdispersed_hmm(
    depth: np.ndarray,
    log_emission: np.ndarray,
    snps=None,
    rho: Optional[float] = None,
    rho_grid: Sequence[float] = DEFAULT_RHO_GRID,
    sticky: Optional[float] = None,
    outlier: float = 0.01,
    max_iter: int = 200,
    tol: float = 1e-7,
) -> Dict[str, np.ndarray]:
    """Fit the overdispersion-aware HMM.

    ``snps`` (a :class:`hmmploidy.data.SnpData`) is required; ``log_emission``
    is only used for its shape (number of windows x ploidy states). Pass a
    fixed ``rho`` to skip its estimation. ``sticky`` defaults to 20% of the
    number of windows.
    """
    if snps is None:
        raise ValueError("The 'overdispersed' model needs per-SNP allele counts (snps=...)")
    T, N = log_emission.shape
    if sticky is None:
        sticky = 0.2 * T
    grid = [rho] if rho is not None else list(rho_grid)

    depth_model = _DepthModel(depth, N, outlier)
    best = None
    profile = []
    for r in grid:
        gl_win = window_log_likelihoods(snps, T, N, r)
        if best is None:
            # initialise from the per-window genotype-likelihood argmax
            depth_model.init(np.argmax(gl_win, axis=1))
            delta, trans = _initial_transitions(N)
        else:
            # warm start from the previous grid point
            delta, trans = fit["delta"], fit["trans"]
        fit = _em(depth_model, gl_win, delta, trans, sticky, max_iter, tol)
        fit["rho"] = r
        fit["gl_win"] = gl_win
        profile.append(fit["objective"])
        if best is None or fit["objective"] > best["objective"]:
            best = fit

    # multi-start at the selected rho, against haploid-depth scale ambiguity
    gl_win = best["gl_win"]
    for k0 in range(N):
        depth_model.init(np.full(T, k0))
        fit = _em(depth_model, gl_win, *_initial_transitions(N), sticky, max_iter, tol)
        if fit["objective"] > best["objective"] + 1e-6:
            fit["rho"] = best["rho"]
            best = fit
    profile[grid.index(best["rho"])] = best["objective"]

    path, _nu = viterbi(best["log_dens"], best["delta"], best["trans"])
    k = np.arange(1, N + 1)
    return {
        "model": "overdispersed",
        "haploid_depth": best["h"],
        "depth_shape": best["shape"],
        "rho": best["rho"],
        "rho_grid": np.asarray(grid, dtype=float),
        "rho_profile_objective": np.asarray(profile),
        "error_rate": snps.error_rate,
        "mu": best["h"] * k,
        "delta": best["delta"],
        "trans": best["trans"],
        "loglik": np.asarray(best["hist"]),
        "states": path + 1,
        "posterior": best["post"],
    }

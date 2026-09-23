from .negbinom import fit_negbinom_hmm
from .overdispersed import fit_overdispersed_hmm
from .poisson_em import fit_poisson_hmm
from .common import viterbi, forward_backward

__all__ = [
    "fit_negbinom_hmm",
    "fit_overdispersed_hmm",
    "fit_poisson_hmm",
    "viterbi",
    "forward_backward",
]

MODEL_CHOICES = ("overdispersed", "negbinom", "em", "bayesian")


def fit(model: str, depth, log_emission, snps=None, **kwargs):
    """Dispatch to one of the supported HMM fitting techniques.

    * ``"overdispersed"`` - Beta-Binomial allele counts + Gamma depth with a
      shared haploid depth, overdispersion estimated from the data (needs
      the per-SNP counts ``snps``; see :mod:`hmmploidy.hmm.overdispersed`).

    * ``"negbinom"`` - classic conditional Newton-Raphson EM with a
      Negative-Binomial emission (translation of the original R method).
    * ``"em"``       - Poisson-emission Baum-Welch EM with closed-form
      M-step updates (faster, no Newton-Raphson).
    * ``"bayesian"`` - hierarchical Bayesian HMM fit with PyMC/NUTS.
    """
    if model == "overdispersed":
        return fit_overdispersed_hmm(depth, log_emission, snps=snps, **kwargs)
    if model == "negbinom":
        return fit_negbinom_hmm(depth, log_emission, **kwargs)
    if model == "em":
        return fit_poisson_hmm(depth, log_emission, **kwargs)
    if model == "bayesian":
        from .bayesian import fit_bayesian_hmm

        return fit_bayesian_hmm(depth, log_emission, **kwargs)
    raise ValueError(f"Unknown model '{model}', choose one of {MODEL_CHOICES}")

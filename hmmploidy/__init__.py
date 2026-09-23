"""HMMploidy: genotype-likelihood- and depth-based inference of ploidy levels.

Pure-Python reimplementation of the original R/Bash/Python toolchain, with a
choice of emission/inference models for the underlying Hidden Markov
Model (see ``hmmploidy.hmm``):

* ``overdispersed`` (default) - Beta-Binomial allele counts with estimated
  overdispersion and a Gamma depth emission with a shared haploid depth;
  robust to allelic and coverage overdispersion.

* ``negbinom``  - classic conditional Newton-Raphson EM on a Negative
  Binomial emission (direct translation of the original ``HMMploidy.R``).
* ``em``        - a Poisson-emission Baum-Welch EM with closed-form M-step
  updates (no Newton-Raphson needed), which is faster and numerically more
  stable for the typical low/medium depth regime used in ploidy calling.
* ``bayesian``  - a hierarchical Bayesian HMM fit with PyMC, placing priors
  on the depth-emission parameters and the transition matrix and returning
  full posterior uncertainty in addition to a MAP Viterbi path.
"""

__version__ = "2.0.0"

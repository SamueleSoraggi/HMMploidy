"""Model behaviour on simulated data with known ploidy."""
import numpy as np
import pandas as pd
import pytest

from hmmploidy import benchmark, data as dm, hmm
from hmmploidy.data import SnpData
from hmmploidy.hmm.overdispersed import snp_log_likelihoods

BENCH_LOCI = 15000  # the genome size of the notebook benchmark (450 windows)


@pytest.fixture(scope="module")
def fitted(sim_dir):
    """The notebook's 'hard' benchmark scenario. On genomes a third of this
    size (150 windows) accuracy drops to ~0.94, with single individuals ~0.8."""
    path = benchmark.simulate_scenario("hard", sim_dir + "/bench", depth=5.0,
                                       loci_per_segment=BENCH_LOCI, seed=1)
    return benchmark.evaluate(path, benchmark.ploidy_layout("changing", 6), BENCH_LOCI)


def test_overdispersed_model_calls_ploidy_under_overdispersion(fitted):
    acc = fitted[fitted.model == "overdispersed"].accuracy
    assert acc.mean() > 0.95
    assert acc.min() > 0.85


def test_overdispersed_model_beats_binomial_models(fitted):
    mean = fitted.groupby("model").accuracy.mean()
    assert mean["overdispersed"] > mean["negbinom"] + 0.1
    assert mean["overdispersed"] > mean["em"] + 0.1


def test_overdispersion_is_estimated(fitted):
    # the 'hard' scenario simulates rho = 0.1
    rho = fitted[fitted.model == "overdispersed"].rho
    assert 0.05 <= rho.median() <= 0.15


def test_haploid_depth_is_recovered(hard_dataset):
    df = dm.read_genolikes(hard_dataset[0])
    iw = dm.prepare_individual(df, 2, tuple(range(1, 7)), 6, 100)
    res = hmm.fit("overdispersed", iw.depth, iw.log_emission, snps=iw.snps)
    assert res["haploid_depth"] == pytest.approx(5.0, rel=0.1)
    np.testing.assert_allclose(res["posterior"].sum(axis=1), 1.0)
    assert set(np.unique(res["states"])) <= set(range(1, 7))


def test_overdispersed_model_requires_snp_data():
    with pytest.raises(ValueError, match="per-SNP allele counts"):
        hmm.fit("overdispersed", np.ones(5), np.zeros((5, 3)))


def idealised_counts(ploidies, seed=0, n_sites=20000, haploid_depth=4):
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.1, 0.9, n_sites)
    minor, total = [], []
    for k in ploidies:
        n = rng.poisson(k * haploid_depth, n_sites)
        minor.append(rng.binomial(n, rng.binomial(k, f) / k))
        total.append(n)
    return np.array(minor).T, np.array(total).T


@pytest.mark.parametrize("ploidies", [[2, 2, 4, 3], [2, 4], [2, 3, 4, 5, 6, 2, 3, 4, 5, 6]])
def test_genotype_likelihoods_unbiased_in_small_cohorts(ploidies):
    """Leave-one-out frequencies with Beta uncertainty: the allele evidence
    alone picks the true ploidy of every individual, even with 2 individuals
    of very different depth (a self-inclusive point estimate calls the
    tetraploids of [2, 2, 4, 3] pentaploid)."""
    minor, total = idealised_counts(ploidies)
    df = pd.DataFrame({
        "ind": np.tile(np.arange(1, len(ploidies) + 1), len(minor)),
        "major_count": (total - minor).ravel(), "minor_count": minor.ravel(),
    })
    pop = dm.population_counts(df)
    for i, k in enumerate(ploidies):
        freq = dm.leave_one_out_frequencies(pop, minor[:, i], total[:, i])
        a, b = dm.frequency_beta_prior(pop, minor[:, i], total[:, i], freq)
        sel = (freq > 0.1) & (freq < 0.9)
        snps = SnpData(np.zeros(sel.sum(), dtype=int), minor[sel, i].astype(float),
                       total[sel, i].astype(float), freq[sel], 1e-3, a[sel], b[sel])
        assert np.argmax(snp_log_likelihoods(snps, 6, 0.0).sum(axis=0)) + 1 == k

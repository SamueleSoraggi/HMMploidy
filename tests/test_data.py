"""Windowing, allele frequencies and per-individual data preparation."""
import dataclasses

import numpy as np
import pandas as pd
import pytest

from hmmploidy import data as dm


def windows_builder_loop(min_size, loci, loci_snp):
    """The original boundary-deletion loop (translation of HMMploidy.R's
    windowsBuilder), kept here as the reference for the vectorised version."""
    loci = np.asarray(loci)
    if len(loci) > 2:
        s = np.concatenate([np.arange(loci[0], loci[-2] + 1, min_size, dtype=float), [float(loci[-1])]])
    elif len(loci) == 2:
        s = np.array([float(loci[0]), float(loci[1])])
    else:
        s = np.array([float(loci[0]), float(loci[0])])
    if loci.max() - loci.min() < min_size:
        s = np.array([float(loci.min()), float(loci.max())])
    s = list(s)
    i = 0
    while i < len(s) - 1:
        if np.count_nonzero((loci_snp > s[i]) & (loci_snp <= s[i + 1])) == 0:
            if i < len(s) - 2:
                del s[i + 1]
            else:
                del s[i]
            continue
        i += 1
    return np.stack([np.asarray(s[:-1]), np.asarray(s[1:])], axis=1)


def test_vectorised_windows_match_original_loop():
    rng = np.random.default_rng(1)
    checked = 0
    for _ in range(1500):
        n_loci = int(rng.integers(2, 300))
        loci = np.sort(rng.choice(5000, n_loci, replace=False)).astype(float)
        n_snp = int(rng.integers(1, max(2, n_loci // int(rng.integers(1, 30)))))
        snp = np.sort(rng.choice(loci, min(n_snp, n_loci), replace=False))
        if np.all(snp == loci[0]):
            continue  # documented edge case: the only SNP is on the first locus
        wind = int(rng.choice([1, 5, 20, 100, 1000]))
        expected = windows_builder_loop(wind, loci, snp)
        got = dm.windows_builder(wind, loci, snp)
        assert got.shape == expected.shape
        np.testing.assert_allclose(got, expected)
        checked += 1
    assert checked > 1000


def test_every_window_contains_a_snp():
    rng = np.random.default_rng(2)
    loci = np.sort(rng.choice(100_000, 5000, replace=False)).astype(float)
    snp = np.sort(rng.choice(loci[1:], 60, replace=False))
    w = dm.windows_builder(500, loci, snp)
    counts = np.searchsorted(snp, w[:, 1], "right") - np.searchsorted(snp, w[:, 0], "right")
    assert np.all(counts >= 1)


def random_counts(seed=0, n_sites=500, n_ind=5):
    rng = np.random.default_rng(seed)
    total = rng.poisson(8, (n_sites, n_ind))
    total[rng.random((n_sites, n_ind)) < 0.1] = 0
    minor = rng.binomial(total, 0.3)
    return total - minor, minor


def counts_frame(major, minor):
    n_sites, n_ind = major.shape
    return pd.DataFrame({
        "ind": np.tile(np.arange(1, n_ind + 1), n_sites),
        "major_count": major.ravel(), "minor_count": minor.ravel(),
    })


def test_all_individual_frequencies_match_original_estimator():
    major, minor = random_counts()
    pop = dm.population_counts(counts_frame(major, minor))
    got = dm.leave_one_out_frequencies(dataclasses.replace(pop, n_ind=1), minor[:, 0],
                                       (major + minor)[:, 0], min_ind=2)
    np.testing.assert_allclose(got, dm.estimate_allele_frequencies(major, minor, min_ind=2))


def test_leave_one_out_frequencies_exclude_the_individual():
    major, minor = random_counts()
    pop = dm.population_counts(counts_frame(major, minor))
    for i in range(major.shape[1]):
        others = [j for j in range(major.shape[1]) if j != i]
        expected = dm.estimate_allele_frequencies(major[:, others], minor[:, others], min_ind=1)
        got = dm.leave_one_out_frequencies(pop, minor[:, i], (major + minor)[:, i], min_ind=1)
        np.testing.assert_allclose(got, expected)


def test_frequency_prior_caps_evidence_per_individual():
    major, minor = random_counts()
    pop = dm.population_counts(counts_frame(major, minor))
    tot0 = (major + minor)[:, 0]
    freq = dm.leave_one_out_frequencies(pop, minor[:, 0], tot0)
    a, b = dm.frequency_beta_prior(pop, minor[:, 0], tot0, freq)
    covered_others = ((major + minor)[:, 1:] > 0).sum(axis=1)
    assert np.all(a + b <= 2 * covered_others + 1 + 1e-9)
    ok = freq > 0
    np.testing.assert_allclose(((a - 0.5) / (a + b - 1))[ok & (covered_others > 0)],
                               np.clip(freq, 1e-3, 1 - 1e-3)[ok & (covered_others > 0)])


def test_read_genolikes_rejects_wrong_row_layout(tmp_path, hard_dataset):
    df = dm.read_genolikes(hard_dataset[0])
    bad = df.iloc[[1, 0] + list(range(2, len(df)))]
    path = tmp_path / "bad.genolikes"
    bad.to_csv(path, sep="\t", header=False, index=False)
    with pytest.raises(ValueError, match="individuals 1..n in order"):
        dm.read_genolikes(str(path))


def test_prepare_individual_snp_data_is_consistent(hard_dataset):
    df = dm.read_genolikes(hard_dataset[0])
    pop = dm.population_counts(df)
    iw = dm.prepare_individual(df, 1, tuple(range(1, 7)), 6, 100, pop=pop)
    snps = iw.snps
    n_windows = len(iw.depth)
    assert iw.log_emission.shape == (n_windows, 6)
    assert np.all(np.diff(snps.window) >= 0)
    assert snps.window.min() == 0 and snps.window.max() == n_windows - 1
    assert np.all(np.bincount(snps.window, minlength=n_windows) >= 1)
    assert np.all((snps.minor >= 0) & (snps.minor <= snps.total))
    assert 1e-4 <= snps.error_rate <= 0.05
    # passing precomputed population counts gives the same result
    iw2 = dm.prepare_individual(df, 1, tuple(range(1, 7)), 6, 100)
    np.testing.assert_allclose(iw.log_emission, iw2.log_emission)

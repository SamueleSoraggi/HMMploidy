"""The compiled genotype-likelihood kernel against the pure-Python
reference implementation (`process_site`)."""
import gzip
import random

import numpy as np
import pytest

from hmmploidy import genotype_likelihoods as gl

PLOIDIES = [1, 2, 3, 4, 5, 6]


def options(**overrides):
    args = gl.build_parser().parse_args(["-s", "1", "-t", "1", "--chunk_size", "500"])
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def reference_lines(lines, nsams, args):
    out = []
    for line in lines:
        out += gl.process_site(
            line, nsams, PLOIDIES, args.min_quality_score, args.min_non_major_freq,
            args.max_minor2_freq, args.max_minor3_freq, args.min_global_depth,
            args.min_ind_depth, float(args.downsampling), random.Random(1),
        )
    return out


def assert_same_genolikes(expected, got):
    assert len(expected) == len(got)
    exp = [line.rstrip("\n").split("\t") for line in expected]
    new = [line.rstrip("\n").split("\t") for line in got]
    # chrom, pos, ind, ref, depth, major, minor, major/minor counts: identical
    assert [e[:9] for e in exp] == [n[:9] for n in new]
    # likelihoods: equal up to floating-point summation order
    np.testing.assert_allclose(
        np.array([[float(v) for v in e[9:]] for e in exp]),
        np.array([[float(v) for v in n[9:]] for n in new]),
        rtol=0, atol=1e-9,
    )


@pytest.fixture(scope="module")
def pileup_lines(mpileup_path, tmp_path_factory):
    lines = gzip.open(mpileup_path, "rt").read().splitlines()[:3000]
    path = tmp_path_factory.mktemp("gl") / "sub.mpileup"
    path.write_text("\n".join(lines) + "\n")
    return lines, str(path)


@pytest.mark.parametrize(
    "overrides",
    [{}, {"min_quality_score": 20}, {"min_ind_depth": 0}, {"min_non_major_freq": 0.05}],
    ids=["defaults", "minQ20", "min_ind_depth0", "min_non_major0.05"],
)
def test_kernel_matches_reference(pileup_lines, tmp_path, overrides):
    lines, path = pileup_lines
    args = options(**overrides)
    nsams = (len(lines[0].split("\t")) - 3) // 3
    out = tmp_path / "out.genolikes.gz"
    gl.run_on_file(path, str(out), PLOIDIES, args)
    got = gzip.open(out, "rt").read().splitlines(True)
    assert len(got) > 0
    assert_same_genolikes(reference_lines(lines, nsams, args), got)


def test_parallel_output_identical_to_serial(pileup_lines, tmp_path):
    _, path = pileup_lines
    outputs = []
    for threads in (1, 3):
        out = tmp_path / f"t{threads}.genolikes.gz"
        gl.run_on_file(path, str(out), PLOIDIES, options(threads=threads))
        outputs.append(gzip.open(out).read())
    assert outputs[0] == outputs[1]


def test_downsampling_reproducible_across_thread_counts(pileup_lines, tmp_path):
    _, path = pileup_lines
    outputs = []
    for threads in (1, 3):
        out = tmp_path / f"ds{threads}.genolikes.gz"
        gl.run_on_file(path, str(out), PLOIDIES, options(threads=threads, downsampling=0.5))
        outputs.append(gzip.open(out).read())
    assert outputs[0] == outputs[1]
    full = tmp_path / "full.genolikes.gz"
    gl.run_on_file(path, str(full), PLOIDIES, options())
    assert outputs[0] != gzip.open(full).read()


def test_lowercase_reference_counts_reference_reads(tmp_path):
    # '.'/',' reads match a soft-masked (lowercase) reference base
    line = "chr1\t10\ta\t6\t...CCC\tIIIIII\t6\t,,,ccc\tIIIIII"
    path = tmp_path / "lower.mpileup"
    path.write_text(line + "\n")
    out = tmp_path / "lower.genolikes.gz"
    gl.run_on_file(str(path), str(out), [2], options())
    rows = [r.split("\t") for r in gzip.open(out, "rt").read().splitlines()]
    assert len(rows) == 2
    for row in rows:
        major_count, minor_count = int(row[7]), int(row[8])
        assert major_count + minor_count == 6
        assert {row[5], row[6]} == {"A", "C"}


def test_sites_with_uncovered_individuals_are_kept_by_default(tmp_path):
    # the second individual has no reads at this site
    line = "chr1\t10\tA\t6\t...CCC\tIIIIII\t0\t*\t*"
    path = tmp_path / "gap.mpileup"
    path.write_text(line + "\n")
    out = tmp_path / "gap.genolikes.gz"
    gl.run_on_file(str(path), str(out), [2], options())
    rows = [r.split("\t") for r in gzip.open(out, "rt").read().splitlines()]
    assert [r[4] for r in rows] == ["6", "0"]
    assert [float(x) for x in rows[1][9:]] == [0.0, 0.0, 0.0]
    # the old default (min_ind_depth 0) drops the site
    gl.run_on_file(str(path), str(out), [2], options(min_ind_depth=0))
    assert gzip.open(out, "rt").read() == ""

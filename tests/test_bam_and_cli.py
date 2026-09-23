"""End-to-end runs through the command line, from BAM and mpileup input."""
import gzip

import numpy as np
import pytest

from hmmploidy import cli
from hmmploidy import run as hmm_run

from conftest import BAM_PLOIDY


def read_output(path):
    fields = {}
    for line in open(path):
        if not line.startswith("#"):
            key, *values = line.rstrip("\n").split("\t")
            fields[key] = values
    return fields


@pytest.fixture(scope="module")
def bam_genolikes(bam_dataset, tmp_path_factory):
    bam_list, ref = bam_dataset
    out = tmp_path_factory.mktemp("bam_out")
    for name, threads in (("stream", "1"), ("regions", "4")):
        cli.main(["genolikes", "--bamList", bam_list, "-r", ref, "-O", name, "-o", str(out),
                  "-t", threads, "--region_size", "15000"])
    return out


def test_bam_region_parallel_output_identical_to_stream(bam_genolikes):
    stream = gzip.open(bam_genolikes / "stream.genolikes.gz").read()
    regions = gzip.open(bam_genolikes / "regions.genolikes.gz").read()
    assert len(stream) > 0
    assert stream == regions


def test_bam_list_writes_sample_names(bam_genolikes):
    names = (bam_genolikes / "stream.names").read_text().split()
    assert names == [f"sample{i + 1}" for i in range(len(BAM_PLOIDY))]


def test_soft_masked_region_is_used(bam_genolikes):
    rows = [line.split("\t") for line in gzip.open(bam_genolikes / "stream.genolikes.gz", "rt")]
    masked = [r for r in rows if r[0] == "chr2" and 10_000 < int(r[1]) <= 30_000]
    chr2 = [r for r in rows if r[0] == "chr2"]
    # the masked stretch is half of chr2 and must hold a similar share of sites
    assert 0.35 < len(masked) / len(chr2) < 0.65
    assert all(int(r[7]) + int(r[8]) > 0 for r in masked if int(r[4]) > 0)


def test_bam_workflow_calls_true_ploidy(bam_genolikes):
    prefix = str(bam_genolikes / "stream")
    cli.main(["run", "--file", prefix, "--wind", "100", "--nameFile", prefix + ".names",
              "-t", "2", "--noPlot"])
    for i, k in enumerate(BAM_PLOIDY):
        out = read_output(f"{prefix}.sample{i + 1}.HMMploidy")
        calls = np.array(out["ploidy"], dtype=int)
        assert np.mean(calls == k) > 0.95, f"sample{i + 1}: expected ploidy {k}"


def test_mpileup_cli_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "ploidy.file").write_text("2x3\n4x3\n")
    cli.main(["simulate", "-p", "ploidy.file", "-l", "3000,3000", "-d", "5", "-o", "sim"])
    cli.main(["genolikes", "sim.filelist", "-t", "2"])
    # the mpileup list is reused for `run`, which reads the .genolikes.gz directly
    cli.main(["run", "--fileList", "sim.filelist", "--wind", "100", "-t", "1"])
    for ind in (1, 2, 3):
        out = read_output(f"sim.DP5.0.ind{ind}.HMMploidy")
        calls = np.array(out["ploidy"], dtype=int)
        starts = np.array(out["window_start"], dtype=float)
        truth = np.where(starts < 3000, 2, 4)
        assert np.mean(calls == truth) > 0.9
        assert (tmp_path / f"sim.DP5.0.ind{ind}.pdf").exists()


def test_run_accepts_genolikes_basenames(hard_dataset, tmp_path):
    path, _ = hard_dataset
    res = hmm_run.analyze_file(path, str(tmp_path / "x"), wind=100, chosen_ind=[1], plot=False)
    assert len(res) == 1 and res[0]["model"] == "overdispersed"
    assert (tmp_path / "x.ind1.HMMploidy").exists()

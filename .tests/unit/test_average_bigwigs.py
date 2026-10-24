"""Unit tests of average_bigwigs.py on tiny bigWigs."""

import shutil
import subprocess

import numpy as np
import polars as pl
import pytest
from figwig import BigWigReader, BigWigWriter

from average_bigwigs import finish, main

SIZES = [("chrA", 250), ("chrB", 100)]


def write_bigwig(path, values):
    """Write dense values per contig, NaN marking uncovered bases."""
    sizes = [(c, len(v)) for c, v in values.items()]
    with BigWigWriter(path, sizes, zooms=0) as w:
        for chrom, v in values.items():
            w.write(chrom, [0], np.asarray(v, np.float32)[None, :], missing=np.nan)


def read_dense(path, chrom, size):
    return BigWigReader(path).read(chrom, [0], size, missing=np.nan)[0]


@pytest.fixture
def group(tmp_path):
    rng = np.random.default_rng(0)
    a1 = rng.uniform(0, 5, 250).round(2)
    a1[10:40] = np.nan  # uncovered
    b1 = rng.uniform(0, 5, 100).round(2)
    a2 = rng.uniform(0, 5, 250).round(2)
    a3 = rng.uniform(0, 5, 250).round(2)
    b3 = np.zeros(100)
    write_bigwig(tmp_path / "1.bw", {"chrA": a1, "chrB": b1})
    write_bigwig(tmp_path / "2.bw", {"chrA": a2})  # lacks chrB
    write_bigwig(tmp_path / "3.bw", {"chrA": a3, "chrB": b3})
    pl.DataFrame(
        {
            "group": ["skinDNase"] * 3,
            "path": [str(tmp_path / f"{i}.bw") for i in (1, 2, 3)],
            "in_tp": [True, False, True],
        }
    ).write_csv(tmp_path / "inputs.tsv", separator="\t")
    pl.DataFrame(
        {
            "group": ["skinDNase", "skinDNase"],
            "set": ["all", "tp"],
            "file": ["skinDNase.bw", "tpSkinDNase.bw"],
        }
    ).write_csv(tmp_path / "tracks.tsv", separator="\t")
    (tmp_path / "chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in SIZES))
    zero = np.nan_to_num
    expected = {
        "all": {
            "chrA": (zero(a1) + a2 + a3) / 3,
            "chrB": (b1 + b3) / 3,
        },
        "tp": {"chrA": (zero(a1) + a3) / 2, "chrB": (b1 + b3) / 2},
    }
    return tmp_path, expected


def run(tmp_path, tile=64):
    main(
        [
            "--inputs", str(tmp_path / "inputs.tsv"),
            "--tracks", str(tmp_path / "tracks.tsv"),
            "--group", "skinDNase",
            "--chrom-sizes", str(tmp_path / "chrom.sizes"),
            "--outdir", str(tmp_path / "out"),
            "--threads", "2",
            "--tile", str(tile),
        ]
    )  # fmt: skip


@pytest.mark.parametrize("tile", [64, 1000])
def test_means_match_float64_reference(group, tile):
    tmp_path, expected = group
    run(tmp_path, tile)
    for name, s in [("skinDNase.bw", "all"), ("tpSkinDNase.bw", "tp")]:
        for chrom, size in SIZES:
            got = read_dense(tmp_path / "out" / name, chrom, size)
            want = np.round(expected[s][chrom], 3).astype(np.float32)
            np.testing.assert_array_equal(got, want)


def test_every_base_is_written(group):
    tmp_path, _ = group
    run(tmp_path)
    for chrom, size in SIZES:
        assert not np.isnan(
            read_dense(tmp_path / "out" / "skinDNase.bw", chrom, size)
        ).any()


@pytest.mark.skipif(shutil.which("bigWigInfo") is None, reason="needs bigWigInfo")
def test_kent_reads_full_coverage(group):
    tmp_path, _ = group
    run(tmp_path)
    info = subprocess.run(
        ["bigWigInfo", str(tmp_path / "out" / "skinDNase.bw")],
        capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip
    assert "basesCovered: 350" in info
    assert "chromCount: 2" in info


def test_negative_zero_is_written_as_zero(tmp_path):
    write_bigwig(tmp_path / "1.bw", {"chrA": np.full(250, -0.0001)})
    pl.DataFrame(
        {"group": ["g"], "path": [str(tmp_path / "1.bw")], "in_tp": [True]}
    ).write_csv(tmp_path / "inputs.tsv", separator="\t")
    pl.DataFrame({"group": ["g"], "set": ["all"], "file": ["g.bw"]}).write_csv(
        tmp_path / "tracks.tsv", separator="\t"
    )
    (tmp_path / "chrom.sizes").write_text("chrA\t250\n")
    main(
        [
            "--inputs", str(tmp_path / "inputs.tsv"),
            "--tracks", str(tmp_path / "tracks.tsv"),
            "--group", "g",
            "--chrom-sizes", str(tmp_path / "chrom.sizes"),
            "--outdir", str(tmp_path / "out"),
            "--threads", "1",
        ]
    )  # fmt: skip
    got = read_dense(tmp_path / "out" / "g.bw", "chrA", 250)
    assert not np.signbit(got).any()
    assert (got == 0).all()


def test_contig_length_mismatch_raises(group):
    tmp_path, _ = group
    (tmp_path / "chrom.sizes").write_text("chrA\t251\nchrB\t100\n")
    with pytest.raises(ValueError, match="chrA has 250 bp"):
        run(tmp_path)


def test_finish_matches_numpy_round():
    rng = np.random.default_rng(2)
    total = np.concatenate(
        [
            rng.uniform(-50, 5000, 100_000),
            np.array([0.0005, 0.0015, 0.0025, 1.2345, -0.0004, 2.0004999, 1e-9]) * 7,
        ]
    )
    for n in (1, 3, 7):
        out = np.empty(total.size + 5, np.float32)
        finish(np.concatenate([total, np.zeros(5)]), n, 1000.0, total.size, out)
        want = (np.round(total / n, 3) + 0.0).astype(np.float32)
        np.testing.assert_array_equal(out[: total.size], want)
        assert np.isnan(out[total.size :]).all()

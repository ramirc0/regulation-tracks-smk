"""The rule-based selection against Gao's 2022 Registry lists."""

import polars as pl
import pytest

# Minimum share of Gao's files that the rules pick again. Misses are files
# archived or revoked since 2022 and DNase replicates that ENCODE no longer
# prefers.
SAME_FILE = {"DNase": 0.93, "ATAC": 0.99, "H3K4me3": 0.98, "H3K27ac": 0.99, "CTCF": 1.0}


@pytest.mark.parametrize("assay", list(SAME_FILE))
def test_rules_pick_gaos_files(gao, ours, assay):
    theirs = gao.filter(pl.col("assay") == assay)
    picked = theirs.join(ours, on="file", how="semi").height
    assert picked / theirs.height >= SAME_FILE[assay], f"{picked} of {theirs.height}"


@pytest.mark.parametrize("assay", list(SAME_FILE))
def test_rules_keep_gaos_experiments(gao, ours, assay):
    theirs = gao.filter(pl.col("assay") == assay)
    kept = theirs.join(ours, on="experiment", how="semi").height
    assert kept / theirs.height >= 0.99, f"{kept} of {theirs.height}"


def test_organs_agree_with_registry(gao, ours):
    both = gao.join(ours, on="experiment", suffix="_ours").filter(
        pl.col("organ") != "UNASSIGNED"
    )
    agree = both.filter(pl.col("organ") == pl.col("organ_ours")).height
    assert agree / both.height >= 0.99, both.filter(
        pl.col("organ") != pl.col("organ_ours")
    )

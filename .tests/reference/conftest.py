"""Shared paths and tables for the reference tests.

These tests read other users' trees and the full run's outputs, so they run
only on the lab cluster. Set RUN_ID to test a run other than the default.
"""

import os
import sys
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

RUN = ROOT / "results" / os.environ.get("RUN_ID", "2026-10-24")
GAO = Path("/home/gaomingshi/Project-slurm/encode-track")
HUB = Path("/zata/public_html/users/gaomingshi/ENCODE_Reg")
CHROM_SIZES = Path(
    "/data/zusers/ramirezc/static/hg38/genome/"
    "GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes"
)
ASSAYS = ["DNase", "ATAC", "H3K4me3", "H3K27ac", "CTCF"]


@pytest.fixture(scope="session")
def gao():
    """Gao's input lists: assay, experiment, file, classification and organ."""
    return pl.concat(
        pl.read_csv(
            GAO / f"{a}-List.txt",
            separator="\t",
            has_header=False,
            infer_schema_length=0,
            truncate_ragged_lines=True,
        ).select(
            pl.lit(a).alias("assay"),
            pl.col("column_0").alias("experiment"),
            pl.col("column_1").alias("file"),
            pl.col("column_5").alias("type"),
            pl.col("column_7").str.replace_all("_", " ").alias("organ"),
        )
        for a in ASSAYS
    )


@pytest.fixture(scope="session")
def ours():
    """The frozen selection of the run under test."""
    path = RUN / "selection" / "inputs.tsv"
    if not path.exists():
        pytest.skip(f"{path} is missing; run `snakemake freeze`")
    return pl.read_csv(path, separator="\t")

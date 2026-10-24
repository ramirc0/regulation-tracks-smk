"""Rebuild the toy inputs of the integration test.

Cuts the metadata snapshot down to a few experiments and writes small random
bigWigs for their files on a toy genome, in the mirror's layout.
Run from the repository root:
`pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/make_data.py`.
"""

from pathlib import Path

import numpy as np
import polars as pl
from figwig import BigWigWriter

SNAPSHOT = Path(
    "/home/ramirezc/Projects/encode-metadata-smk/results/2026-10-09/parquets"
)
OUT = Path(__file__).parent / "resources"
# Skin DNase: 2 tissue or primary cell, 1 cell line. Thyroid ATAC: 2 tissue.
# Urinary bladder total RNA-seq: 2 experiments, both strands.
EXPERIMENTS = ["ENCSR000EPP", "ENCSR092OFK", "ENCSR000EIN", "ENCSR201FIW",
               "ENCSR474XFV", "ENCSR000AAC", "ENCSR000AEV"]  # fmt: skip
SIZES = [("chr1", 2000), ("chr2", 1000), ("chrM", 100)]
TYPES = [
    "read-depth normalized signal",
    "fold change over control",
    "plus strand signal of unique reads",
    "minus strand signal of unique reads",
]


def main():
    """Write the toy metadata, mirror and chrom sizes."""
    meta = OUT / "metadata"
    meta.mkdir(parents=True, exist_ok=True)
    exp = pl.read_parquet(SNAPSHOT / "experiment.parquet").filter(
        pl.col("accession").is_in(EXPERIMENTS)
    )
    assert exp.height == len(EXPERIMENTS), set(EXPERIMENTS) - set(exp["accession"])
    files = pl.read_parquet(SNAPSHOT / "file.parquet").filter(
        pl.col("dataset").is_in(exp["id"].implode())
        & (pl.col("file_format") == "bigWig")
    )
    types = pl.read_parquet(SNAPSHOT / "biosample_type.parquet").filter(
        pl.col("id").is_in(exp["biosample_ontology"].implode())
    )
    exp.write_parquet(meta / "experiment.parquet")
    files.write_parquet(meta / "file.parquet")
    types.write_parquet(meta / "biosample_type.parquet")
    (OUT / "chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in SIZES))

    rng = np.random.default_rng(1)
    toy = files.filter(
        (pl.col("genome_assembly") == "GRCh38") & pl.col("output_type").is_in(TYPES)
    )
    for row in toy.iter_rows(named=True):
        path = (
            OUT
            / "encode"
            / row["dataset"].split("/")[2]
            / row["download_url"].split("/")[-1]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        with BigWigWriter(path, SIZES[:2], zooms=0) as w:  # no chrM, as in ATAC files
            for chrom, n in SIZES[:2]:
                v = rng.uniform(0, 3, n).round(2).astype(np.float32)
                v[rng.uniform(size=n) < 0.2] = np.nan
                w.write(chrom, [0], v[None, :], missing=np.nan)


if __name__ == "__main__":
    main()

"""Download the selected bigWigs that the local ENCODE mirror lacks.

Runs outside Snakemake, after `snakemake freeze` and before the main run.
Files come from ENCODE's public S3 bucket and must match their MD5.
"""

import argparse
import hashlib
import os
import subprocess
from pathlib import Path

import polars as pl


def build_parser():
    """Build the fetch_missing.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the fetch_missing.py flags.
    """
    p = argparse.ArgumentParser(description="Download selected files the mirror lacks.")
    p.add_argument(
        "--inputs", required=True, help="Inputs TSV from `snakemake freeze`."
    )
    p.add_argument("--metadata", required=True, help="encode-metadata-smk parquet dir.")
    p.add_argument("--dry-run", action="store_true", help="Only list the files.")
    return p


def md5(path):
    """Return the MD5 hex digest of a file.

    Parameters
    ----------
    path : pathlib.Path
        File to hash.

    Returns
    -------
    str
        Hex digest.
    """
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    """Download every selected file whose path does not exist.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    inputs = pl.read_csv(args.inputs, separator="\t")
    missing = inputs.filter(
        ~pl.col("path").map_elements(os.path.exists, return_dtype=bool)
    )
    files = pl.read_parquet(
        Path(args.metadata) / "file.parquet", columns=["accession", "s3_uri", "md5sum"]
    )
    todo = missing.join(files, left_on="file", right_on="accession", how="left")
    print(f"{todo.height} of {inputs.height} files are missing.")
    for r in todo.iter_rows(named=True):
        print(r["s3_uri"], "->", r["path"])
        if args.dry_run:
            continue
        dest = Path(r["path"])
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_name(dest.name + ".part")
        subprocess.run(
            ["s5cmd", "--no-sign-request", "cp", r["s3_uri"], str(part)], check=True
        )
        if md5(part) != r["md5sum"]:
            part.unlink()
            raise ValueError(f"MD5 mismatch for {r['file']}")
        part.rename(dest)


if __name__ == "__main__":
    main()

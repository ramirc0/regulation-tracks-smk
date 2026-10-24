"""Select the ENCODE bigWigs that feed each organ-average track.

Reads an encode-metadata-smk snapshot and writes two tables. The inputs table
has one row per selected file. The tracks table has one row per output track.
"""

import argparse
import re
import sys
from pathlib import Path

import polars as pl

# Assay, track suffix, experiment filter and the output type of the chosen file.
# "fc" takes the fold change over control file that matches the preferred
# signal p-value file's analysis and replicates.
ASSAYS = [
    ("DNase", "DNase", pl.col("assay_title") == "DNase-seq", "preferred"),
    ("ATAC", "ATAC", pl.col("assay_title") == "ATAC-seq", "fc"),
    (
        "H3K4me3",
        "H3K4me3",
        pl.col("assay_title").is_in(["Histone ChIP-seq", "Mint-ChIP-seq"])
        & (pl.col("target_of_assay") == "H3K4me3"),
        "fc",
    ),
    (
        "H3K27ac",
        "H3K27ac",
        pl.col("assay_title").is_in(["Histone ChIP-seq", "Mint-ChIP-seq"])
        & (pl.col("target_of_assay") == "H3K27ac"),
        "fc",
    ),
    (
        "CTCF",
        "CTCF",
        (pl.col("assay_title") == "TF ChIP-seq")
        & (pl.col("target_of_assay") == "CTCF"),
        "fc",
    ),
    (
        "RNA",
        "",
        (pl.col("assay_title") == "total RNA-seq")
        & (pl.col("cellular_component").fill_null("") == ""),
        "preferred",
    ),
]

PREFERRED_TYPES = {
    "read-depth normalized signal": "",
    "plus strand signal of unique reads": "Plus",
    "minus strand signal of unique reads": "Minus",
}
TP_CLASSES = ["tissue", "primary cell"]


def build_parser():
    """Build the select_inputs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the select_inputs.py flags.
    """
    p = argparse.ArgumentParser(description="Select inputs of organ-average tracks.")
    p.add_argument("--metadata", required=True, help="encode-metadata-smk parquet dir.")
    p.add_argument("--organs", required=True, help="Organ table (config/organs.tsv).")
    p.add_argument("--encode-root", required=True, help="Local ENCODE mirror.")
    p.add_argument("--download-dir", required=True, help="Where missing files go.")
    p.add_argument("--inputs", required=True, help="Inputs TSV to write.")
    p.add_argument("--tracks", required=True, help="Tracks TSV to write.")
    return p


def read_organs(path):
    """Read the organ table.

    Parameters
    ----------
    path : str
        Organ TSV with columns organ, slims, fallback, pattern and color.

    Returns
    -------
    polars.DataFrame
        One row per organ in priority order, with `slims` and `fallback` as
        lists of str.
    """
    organs = pl.read_csv(
        path, separator="\t", comment_prefix="#", infer_schema_length=0
    )
    return organs.with_columns(
        pl.col("slims", "fallback")
        .fill_null("")
        .str.split(",")
        .list.eval(pl.element().filter(pl.element() != "")),
        pl.col("pattern").fill_null(""),
    )


def assign_organ(slims, fallback, name, organs):
    """Assign one organ to a biosample type.

    Parameters
    ----------
    slims : str or None
        Comma-joined ENCODE organ_slims.
    fallback : str or None
        Comma-joined ENCODE system_slims and cell slims.
    name : str
        Biosample term name.
    organs : polars.DataFrame
        Organ table from `read_organs`.

    Returns
    -------
    str or None
        The first organ whose slims match, else the first whose fallback
        terms match, else the first whose pattern matches the name. None
        when nothing matches.
    """
    rows = list(organs.iter_rows(named=True))
    for terms, column in [(slims, "slims"), (fallback, "fallback")]:
        found = set((terms or "").split(","))
        for row in rows:
            if found.intersection(row[column]):
                return row["organ"]
    for row in rows:
        if row["pattern"] and re.search(row["pattern"], name):
            return row["organ"]
    return None


def stem(organ):
    """Return the camel-case file stem of an organ.

    Parameters
    ----------
    organ : str
        Organ name, such as "large intestine".

    Returns
    -------
    str
        Stem, such as "largeIntestine".
    """
    first, *rest = organ.split()
    return first + "".join(w.capitalize() for w in rest)


def encode_path(encode_root, download_dir, experiment, name):
    """Return where an ENCODE file is read from.

    Parameters
    ----------
    encode_root : str
        Local mirror, laid out as `<experiment>/<name>`.
    download_dir : str
        Same layout, for files the mirror lacks.
    experiment : str
        Experiment accession.
    name : str
        File name, such as `ENCFF000AAA.bigWig`.

    Returns
    -------
    str
        The mirror path if it exists, else the download path.
    """
    mirror = Path(encode_root, experiment, name)
    return str(mirror) if mirror.exists() else f"{download_dir}/{experiment}/{name}"


def select_files(experiments, files):
    """Choose the bigWig of each experiment and strand.

    Parameters
    ----------
    experiments : polars.DataFrame
        Released human experiments with an `assay` and a `rule` column.
    files : polars.DataFrame
        Released GRCh38 bigWigs of experiments.

    Returns
    -------
    polars.DataFrame
        One row per experiment and strand: assay, suffix, experiment, file and
        download_url.

    Raises
    ------
    ValueError
        If an experiment has more than one candidate file for a strand.
    """
    f = files.join(experiments, left_on="dataset", right_on="id").filter(
        pl.col("analyses").str.split(",").list.contains(pl.col("default_analysis"))
    )
    preferred = f.filter(
        (pl.col("rule") == "preferred")
        & pl.col("preferred_default").fill_null(False)
        & pl.col("output_type").is_in(list(PREFERRED_TYPES))
    ).with_columns(
        pl.col("suffix") + pl.col("output_type").replace_strict(PREFERRED_TYPES)
    )
    keys = ["dataset", "analyses", "biological_replicates"]
    pvalue = f.filter(
        (pl.col("rule") == "fc")
        & pl.col("preferred_default").fill_null(False)
        & (pl.col("output_type") == "signal p-value")
    ).select(keys)
    fc = f.filter(
        (pl.col("rule") == "fc") & (pl.col("output_type") == "fold change over control")
    ).join(pvalue, on=keys)
    cols = ["assay", "suffix", "accession_right", "accession", "download_url"]
    chosen = (
        pl.concat([preferred.select(cols), fc.select(cols)])
        .rename({"accession_right": "experiment", "accession": "file"})
        .sort("assay", "suffix", "experiment")
    )
    dup = chosen.filter(pl.struct("experiment", "suffix").is_duplicated())
    if dup.height:
        raise ValueError(f"Several files for one experiment and strand:\n{dup}")
    return chosen


def main(argv=None):
    """Select inputs and write the inputs and tracks tables.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    meta = Path(args.metadata)
    organs = read_organs(args.organs)

    experiments = (
        pl.read_parquet(meta / "experiment.parquet")
        .filter(
            (pl.col("organism") == "Homo sapiens") & (pl.col("status") == "released")
        )
        .with_columns(
            pl.coalesce(
                pl.when(cond).then(
                    pl.struct(assay=pl.lit(a), suffix=pl.lit(s), rule=pl.lit(r))
                )
                for a, s, cond, r in ASSAYS
            ).alias("spec")
        )
        .filter(pl.col("spec").is_not_null())
        .unnest("spec")
        .select(
            "id",
            "accession",
            "assay",
            "suffix",
            "rule",
            "default_analysis",
            "biosample_ontology",
            "biosample_term_name",
            "biosample_classification",
        )
    )
    files = pl.read_parquet(
        meta / "file.parquet",
        columns=[
            "accession",
            "dataset",
            "file_format",
            "genome_assembly",
            "status",
            "output_type",
            "preferred_default",
            "analyses",
            "biological_replicates",
            "restricted_file",
            "no_file_available",
            "download_url",
        ],
    ).filter(
        (pl.col("file_format") == "bigWig")
        & (pl.col("genome_assembly") == "GRCh38")
        & (pl.col("status") == "released")
        & ~pl.col("restricted_file").fill_null(False)
        & ~pl.col("no_file_available").fill_null(False)
    )
    chosen = select_files(experiments, files)

    types = pl.read_parquet(
        meta / "biosample_type.parquet",
        columns=["id", "organ", "system_slims", "cell", "biosample_name"],
    )
    organ_of = {
        r["id"]: assign_organ(
            r["organ"],
            ",".join(filter(None, [r["system_slims"], r["cell"]])),
            r["biosample_name"],
            organs,
        )
        for r in types.iter_rows(named=True)
    }

    inputs = chosen.join(
        experiments.select(
            pl.col("accession").alias("experiment"),
            "biosample_ontology",
            pl.col("biosample_term_name").alias("biosample"),
            pl.col("biosample_classification").alias("classification"),
        ),
        on="experiment",
    ).with_columns(
        pl.col("biosample_ontology")
        .replace_strict(organ_of, default=None)
        .alias("organ"),
        pl.col("classification").is_in(TP_CLASSES).alias("in_tp"),
        pl.col("download_url").str.extract(r"([^/]+)$").alias("name"),
    )
    unassigned = inputs.filter(pl.col("organ").is_null())
    print(f"Dropped {unassigned.height} experiments with no organ:", file=sys.stderr)
    print(
        unassigned.group_by("assay", "biosample").len().sort("len", descending=True),
        file=sys.stderr,
    )
    inputs = inputs.filter(pl.col("organ").is_not_null())

    inputs = (
        inputs.with_columns(
            pl.struct("experiment", "name")
            .map_elements(
                lambda r: encode_path(
                    args.encode_root, args.download_dir, r["experiment"], r["name"]
                ),
                return_dtype=pl.String,
            )
            .alias("path"),
            (
                pl.col("organ").map_elements(stem, return_dtype=pl.String)
                + pl.col("suffix")
            ).alias("group"),
        )
        .select(
            "group",
            "assay",
            "experiment",
            "file",
            "biosample",
            "classification",
            "organ",
            "in_tp",
            "path",
        )
        .sort("group", "experiment", "file")
    )

    tracks = (
        inputs.group_by("group", "assay", "organ")
        .agg(n_all=pl.len(), n_tp=pl.col("in_tp").sum())
        .with_columns(pl.col("n_all", "n_tp").cast(pl.Int64))
    )
    all_sets = tracks.select(
        "group",
        "assay",
        "organ",
        pl.lit("all").alias("set"),
        pl.col("n_all").alias("n"),
        (pl.col("group") + ".bw").alias("file"),
    )
    tp_sets = tracks.filter(
        (pl.col("n_tp") > 0) & (pl.col("n_tp") < pl.col("n_all"))
    ).select(
        "group",
        "assay",
        "organ",
        pl.lit("tp").alias("set"),
        pl.col("n_tp").alias("n"),
        (
            "tp"
            + pl.col("group").str.slice(0, 1).str.to_uppercase()
            + pl.col("group").str.slice(1)
            + ".bw"
        ).alias("file"),
    )
    tracks = pl.concat([all_sets, tp_sets]).sort("group", "set")

    inputs.write_csv(args.inputs, separator="\t")
    tracks.write_csv(args.tracks, separator="\t")
    print(
        inputs.group_by("assay")
        .agg(
            experiments=pl.col("experiment").n_unique(),
            files=pl.len(),
            groups=pl.col("group").n_unique(),
            to_download=pl.col("path").str.starts_with(args.download_dir).sum(),
        )
        .sort("assay"),
        file=sys.stderr,
    )
    print(tracks.group_by("assay", "set").len().sort("assay", "set"), file=sys.stderr)


if __name__ == "__main__":
    main()

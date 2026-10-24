"""Write a one-file UCSC track hub of the organ-average tracks.

Containers, labels, colors and default visibility follow UCSC's ENCODE
Regulation v4 track (`wgEncodeReg4`).
"""

import argparse

import polars as pl

SUPER = "wgEncodeReg4"

# Assay: container suffix, measure in labels, container settings.
CONTAINERS = {
    "DNase": (
        "Dnase",
        "DNase",
        {
            "shortLabel": "DNase (Layered)",
            "longLabel": "Chromatin accessibility from DNase-seq signal, "
            "averaged by organ/tissue",
            "visibility": "hide",
            "viewLimits": "0:50",
            "priority": "1.1",
        },
    ),
    "ATAC": (
        "Atac",
        "ATAC",
        {
            "shortLabel": "ATAC (Layered)",
            "longLabel": "Chromatin accessibility from ATAC-seq signal, "
            "averaged by organ/tissue",
            "visibility": "hide",
            "viewLimits": "0:50",
            "priority": "1.2",
        },
    ),
    "H3K4me3": (
        "MarkH3k4me3",
        "H3K4me3",
        {
            "shortLabel": "H3K4me3 (Layered)",
            "longLabel": "H3K4me3 signal marking active and poised promoters, "
            "averaged by organ/tissue",
            "visibility": "hide",
            "viewLimits": "0:100",
            "priority": "1.3",
        },
    ),
    "H3K27ac": (
        "MarkH3k27ac",
        "H3K27ac",
        {
            "shortLabel": "H3K27ac (Layered)",
            "longLabel": "H3K27ac signal marking active enhancers and promoters, "
            "averaged by organ/tissue",
            "visibility": "full",
            "viewLimits": "0:100",
            "priority": "1.4",
        },
    ),
    "CTCF": (
        "MarkCtcf",
        "CTCF",
        {
            "shortLabel": "CTCF (Layered)",
            "longLabel": "CTCF binding signal from ChIP-seq, averaged by organ/tissue",
            "visibility": "hide",
            "viewLimits": "0:100",
            "priority": "1.5",
        },
    ),
    "RNA": (
        "Txn",
        "total RNA-seq",
        {
            "shortLabel": "Transcription (Layered)",
            "longLabel": "Strand-specific transcription signal measured by total "
            "RNA-seq, averaged by organ/tissue",
            "visibility": "hide",
            "priority": "1.6",
        },
    ),
}
SHARED = {
    "type": "bigWig",
    "container": "multiWig",
    "aggregate": "transparentOverlay",
    "showSubtrackColorOnUi": "on",
    "maxHeightPixels": "100:50:11",
    "autoScale": "on",
    "allButtonPair": "on",
    "dragAndDrop": "subtracks",
    "noInherit": "on",
}
SHOWN = {"blood", "brain", "kidney", "liver", "muscle"}
SET_LABEL = {"all": "all biosamples", "tp": "tissues and primary cells only"}
STRANDS = {"Plus": "+", "Minus": "-"}


def build_parser():
    """Build the make_hub.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the make_hub.py flags.
    """
    p = argparse.ArgumentParser(description="Write the track hub.")
    p.add_argument("--tracks", required=True, help="Tracks TSV from select_inputs.py.")
    p.add_argument("--organs", required=True, help="Organ table with colors.")
    p.add_argument("--email", required=True, help="Hub contact email.")
    p.add_argument("--output", required=True, help="hub.txt to write.")
    return p


def long_label(measure, n, organ, track_set):
    """Return a leaf track's long label.

    Parameters
    ----------
    measure : str
        What is averaged, such as "DNase" or "+ strand total RNA-seq".
    n : int
        Number of experiments.
    organ : str
        Organ name.
    track_set : str
        "all" or "tp".

    Returns
    -------
    str
        For example "Avg. DNase level of 46 skin experiments (all biosamples)".
    """
    avg, plural = ("Avg. ", "s") if n > 1 else ("", "")
    return f"{avg}{measure} level of {n} {organ} experiment{plural} ({SET_LABEL[track_set]})"


def leaves(tracks, colors):
    """Yield the stanza of each leaf track, in display order.

    Parameters
    ----------
    tracks : polars.DataFrame
        Tracks table: group, assay, organ, set, n and file.
    colors : dict of str to str
        RGB color of each organ.

    Yields
    ------
    tuple of (str, dict of str to str)
        Assay and the settings of one leaf.
    """
    tracks = tracks.with_columns(
        pl.col("group").str.extract(r"(Plus|Minus)$").fill_null("").alias("strand"),
        (pl.col("set") == "all").alias("is_all"),
    ).with_columns(
        # A tp track leads its organ. Otherwise the all track does.
        (pl.col("is_all") & (pl.len().over("group") > 1)).alias("secondary")
    )
    tracks = tracks.sort("assay", "secondary", "organ", pl.col("strand") == "Minus")
    priority = {}
    for r in tracks.iter_rows(named=True):
        container, measure, _ = CONTAINERS[r["assay"]]
        organ = r["organ"]
        camel = "".join(w.capitalize() for w in organ.split())
        sign = STRANDS.get(r["strand"], "")
        if sign:
            measure = f"{sign} strand {measure}"
        short = organ.capitalize() + (f" {sign}" if sign else "")
        name = f"{SUPER}{container}{'All' if r['is_all'] else ''}{camel}{r['strand']}"
        if r["is_all"]:
            short += " (all biosamples)"
        shown = organ in SHOWN and not r["secondary"]
        priority[container] = priority.get(container, 0) + 1
        leaf = {
            "track": name,
            "parent": f"{SUPER}{container}" + ("" if shown else " off"),
            "shortLabel": short,
            "longLabel": long_label(measure, r["n"], organ, r["set"]),
            "type": "bigWig",
            "bigDataUrl": f"organAve/{r['group']}/{r['file']}",
            "color": colors[organ],
            "priority": str(priority[container]),
        }
        if r["strand"] == "Minus":
            leaf["negateValues"] = "on"
        yield r["assay"], leaf


def stanza(settings, indent=""):
    """Format one trackDb stanza.

    Parameters
    ----------
    settings : dict of str to str
        Setting names and values, in order.
    indent : str, optional
        Prefix of every line.

    Returns
    -------
    str
        The stanza, ending in a blank line.
    """
    return "".join(f"{indent}{k} {v}\n" for k, v in settings.items()) + "\n"


def main(argv=None):
    """Write hub.txt.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    tracks = pl.read_csv(args.tracks, separator="\t")
    organs = pl.read_csv(args.organs, separator="\t", comment_prefix="#")
    colors = dict(zip(organs["organ"], organs["color"]))

    out = [
        stanza(
            {
                "hub": "ENCODE4RegulationAverages",
                "shortLabel": "ENCODE4 organ averages",
                "longLabel": "ENCODE 4 regulation signal averaged by organ/tissue",
                "useOneFile": "on",
                "email": args.email,
            }
        ),
        stanza({"genome": "hg38"}),
        stanza(
            {
                "track": SUPER,
                "superTrack": "on show",
                "shortLabel": "ENC4 Regulation",
                "longLabel": "ENCODE 4 regulation signal averaged by organ/tissue",
            }
        ),
    ]
    by_assay = {}
    for assay, leaf in leaves(tracks, colors):
        by_assay.setdefault(assay, []).append(leaf)
    for assay, (container, _, settings) in CONTAINERS.items():
        if assay not in by_assay:
            continue
        out.append(
            stanza(
                {"track": f"{SUPER}{container}", "parent": SUPER} | settings | SHARED
            )
        )
        out += [stanza(leaf, "    ") for leaf in by_assay[assay]]
    with open(args.output, "w") as fh:
        fh.write("".join(out))


if __name__ == "__main__":
    main()

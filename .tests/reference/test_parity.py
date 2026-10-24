"""Averages against Gao's hub files, which UCSC serves as wgEncodeReg4.

`test_rebuild_matches_hub` reruns the averaging on Gao's own inputs. The
other test compares the full run's tracks where our inputs equal Gao's.
"""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from conftest import CHROM_SIZES, HUB, RUN
from figwig import BigWigReader

from average_bigwigs import main

MIRROR = Path("/data/projects/encode/data")
CHROMS = ["chr21", "chrX"]
TP_TYPES = ["tissue", "primary_cell"]

# Hub files with known defects (see .claude/CLAUDE.md).
DEFECTS = {
    "ave.skin.CTCF", "ave.mouth.DNase", "ave.brain.ATAC", "ave.heart.ATAC",
    "ave.brain.CTCF", "ave.brain.H3K27ac", "ave.lung.H3K27ac", "ave.skin.H3K27ac",
    "ave.brain.H3K4me3", "ave.embryo.H3K4me3", "ave.lung.H3K4me3", "ave.skin.H3K4me3",
}  # fmt: skip


def read(path, chrom):
    r = BigWigReader(path)
    return r.read(chrom, [0], r.chrom_sizes[chrom], missing=0.0)[0]


def rebuild(tmp_path, paths):
    """Average `paths` on CHROMS and return the written track."""
    for p in paths:
        if not Path(p).exists():
            pytest.skip(f"{p} is not on the mirror")
    sizes = dict(line.split()[:2] for line in CHROM_SIZES.read_text().splitlines())
    (tmp_path / "chrom.sizes").write_text("".join(f"{c}\t{sizes[c]}\n" for c in CHROMS))
    pl.DataFrame({"group": "g", "path": paths, "in_tp": True}).write_csv(
        tmp_path / "inputs.tsv", separator="\t"
    )
    pl.DataFrame({"group": ["g"], "set": ["all"], "file": ["g.bw"]}).write_csv(
        tmp_path / "tracks.tsv", separator="\t"
    )
    main(
        [
            "--inputs", str(tmp_path / "inputs.tsv"),
            "--tracks", str(tmp_path / "tracks.tsv"),
            "--group", "g",
            "--chrom-sizes", str(tmp_path / "chrom.sizes"),
            "--outdir", str(tmp_path),
        ]
    )  # fmt: skip
    return tmp_path / "g.bw"


@pytest.mark.parametrize(
    "assay, organ, hub_file",
    [
        ("DNase", "uterus", "Epigenetics/average_bigWig/tp.uterus.DNase.bigWig"),
        (
            "H3K4me3",
            "prostate",
            "Epigenetics/average_bigWig/tp.prostate.H3K4me3.bigWig",
        ),
        ("RNA", None, "Transcription/average_bigWig/urinary_bladder.plus.bigWig"),
    ],
)
def test_rebuild_matches_hub(gao, tmp_path, assay, organ, hub_file):
    if assay == "RNA":
        files = [("ENCSR000AAC", "ENCFF442VHH"), ("ENCSR000AEV", "ENCFF775XPO")]
    else:
        rows = gao.filter(
            (pl.col("assay") == assay)
            & (pl.col("organ") == organ)
            & pl.col("type").is_in(TP_TYPES)
        )
        files = rows.select("experiment", "file").rows()
    track = rebuild(tmp_path, [str(MIRROR / e / f"{a}.bigWig") for e, a in files])
    for chrom in CHROMS:
        np.testing.assert_array_equal(read(track, chrom), read(HUB / hub_file, chrom))


def hub_name(set_, organ, assay):
    return f"{'tp' if set_ == 'tp' else 'ave'}.{organ.replace(' ', '_')}.{assay}"


def unchanged_tracks(gao, ours):
    """List tracks whose input files equal Gao's."""
    tracks = pl.read_csv(RUN / "selection" / "tracks.tsv", separator="\t")
    out = []
    for r in tracks.filter(pl.col("assay") != "RNA").iter_rows(named=True):
        mine = ours.filter(pl.col("group") == r["group"])
        theirs = gao.filter(
            (pl.col("assay") == r["assay"]) & (pl.col("organ") == r["organ"])
        )
        if r["set"] == "tp":
            mine = mine.filter(pl.col("in_tp"))
            theirs = theirs.filter(pl.col("type").is_in(TP_TYPES))
        name = hub_name(r["set"], r["organ"], r["assay"])
        if set(mine["file"]) == set(theirs["file"]) and name not in DEFECTS:
            out.append((r["group"], r["file"], name))
    return out


def test_run_matches_hub_where_inputs_are_unchanged(gao, ours):
    hub_dir = RUN / "hub" / "organAve"
    if not hub_dir.exists():
        pytest.skip(f"{hub_dir} is missing; run the workflow")
    checked = 0
    for group, file, name in unchanged_tracks(gao, ours):
        theirs = HUB / "Epigenetics" / "average_bigWig" / f"{name}.bigWig"
        if not theirs.exists():
            continue
        for chrom in CHROMS:
            got, want = read(hub_dir / group / file, chrom), read(theirs, chrom)
            if name.startswith("tp."):
                np.testing.assert_array_equal(got, want, err_msg=name)
            else:
                # ave files are bigWigMerge sums over n, unrounded, about 6 digits.
                np.testing.assert_allclose(
                    got, want, rtol=1e-5, atol=5e-4, err_msg=name
                )
        checked += 1
    assert checked > 0

"""Unit tests of make_hub.py."""

import polars as pl

from make_hub import leaves, long_label, main

TRACKS = pl.DataFrame(
    {
        "group": ["bloodDNase", "bloodDNase", "boneDNase", "skinPlus", "skinMinus"],
        "assay": ["DNase", "DNase", "DNase", "RNA", "RNA"],
        "organ": ["blood", "blood", "bone", "skin", "skin"],
        "set": ["all", "tp", "all", "all", "all"],
        "n": [401, 359, 1, 18, 18],
        "file": [
            "bloodDNase.bw", "tpBloodDNase.bw", "boneDNase.bw", "skinPlus.bw", "skinMinus.bw",
        ],
    }
)  # fmt: skip
COLORS = {"blood": "254,75,173", "bone": "121,147,150", "skin": "127,133,209"}


def test_long_label_matches_ucsc_wording():
    assert (
        long_label("DNase", 46, "skin", "all")
        == "Avg. DNase level of 46 skin experiments (all biosamples)"
    )
    assert (
        long_label("DNase", 1, "bone", "tp")
        == "DNase level of 1 bone experiment (tissues and primary cells only)"
    )


def test_leaves_follow_ucsc_names_order_and_visibility():
    got = {leaf["track"]: (assay, leaf) for assay, leaf in leaves(TRACKS, COLORS)}
    dnase = [t for t, (a, _) in got.items() if a == "DNase"]
    # tp tracks and sole all tracks lead; secondary all tracks follow.
    assert dnase == [
        "wgEncodeReg4DnaseBlood",
        "wgEncodeReg4DnaseAllBone",
        "wgEncodeReg4DnaseAllBlood",
    ]
    blood_tp = got["wgEncodeReg4DnaseBlood"][1]
    assert blood_tp["parent"] == "wgEncodeReg4Dnase"
    assert blood_tp["shortLabel"] == "Blood"
    assert blood_tp["bigDataUrl"] == "organAve/bloodDNase/tpBloodDNase.bw"
    assert got["wgEncodeReg4DnaseAllBlood"][1]["parent"] == "wgEncodeReg4Dnase off"
    minus = got["wgEncodeReg4TxnAllSkinMinus"][1]
    assert minus["negateValues"] == "on"
    assert minus["shortLabel"] == "Skin - (all biosamples)"
    assert minus["longLabel"].startswith("Avg. - strand total RNA-seq level of 18 skin")


def test_main_writes_one_file_hub(tmp_path):
    TRACKS.write_csv(tmp_path / "tracks.tsv", separator="\t")
    pl.DataFrame({"organ": list(COLORS), "color": list(COLORS.values())}).write_csv(
        tmp_path / "organs.tsv", separator="\t"
    )
    main(
        [
            "--tracks", str(tmp_path / "tracks.tsv"),
            "--organs", str(tmp_path / "organs.tsv"),
            "--email", "a@b.org",
            "--output", str(tmp_path / "hub.txt"),
        ]
    )  # fmt: skip
    hub = (tmp_path / "hub.txt").read_text()
    assert hub.startswith("hub ENCODE4RegulationAverages\n")
    assert "useOneFile on" in hub
    assert hub.count("\n    track ") == 5
    assert "track wgEncodeReg4Txn\nparent wgEncodeReg4\n" in hub
    assert "track wgEncodeReg4Atac" not in hub

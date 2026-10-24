"""Unit tests of select_inputs.py."""

import polars as pl
import pytest

from select_inputs import assign_organ, read_organs, select_files, stem

ORGANS = """organ\tslims\tfallback\tpattern\tcolor
skin\tskin of body\t\t\t1,1,1
brain\tbrain\tcentral nervous system\t\t2,2,2
blood\tblood\timmune system\t^GM\\d+$\t3,3,3
epithelium\tepithelium\t\t\t4,4,4
"""


@pytest.fixture
def organs(tmp_path):
    path = tmp_path / "organs.tsv"
    path.write_text(ORGANS)
    return read_organs(path)


def experiment(acc, assay, suffix, rule, analysis="/analyses/A/"):
    return {
        "id": f"/experiments/{acc}/",
        "accession": acc,
        "assay": assay,
        "suffix": suffix,
        "rule": rule,
        "default_analysis": analysis,
    }


def bigwig(acc, exp, output_type, reps="1,2", preferred=None, analyses="/analyses/A/"):
    return {
        "accession": acc,
        "dataset": f"/experiments/{exp}/",
        "output_type": output_type,
        "preferred_default": preferred,
        "analyses": analyses,
        "biological_replicates": reps,
        "download_url": f"/files/{acc}/@@download/{acc}.bigWig",
    }


def test_fc_matches_preferred_pvalue_replicates():
    exps = pl.DataFrame([experiment("ENCSR1", "ATAC", "ATAC", "fc")])
    files = pl.DataFrame(
        [
            bigwig("P12", "ENCSR1", "signal p-value", "1,2", True),
            bigwig("P1", "ENCSR1", "signal p-value", "1"),
            bigwig("F1", "ENCSR1", "fold change over control", "1"),
            bigwig("F12", "ENCSR1", "fold change over control", "1,2"),
            bigwig(
                "OLD", "ENCSR1", "fold change over control", "1,2", None, "/analyses/B/"
            ),
        ]
    )
    chosen = select_files(exps, files)
    assert chosen["file"].to_list() == ["F12"]
    assert chosen["suffix"].to_list() == ["ATAC"]


def test_preferred_rule_keeps_default_analysis_only():
    exps = pl.DataFrame([experiment("ENCSR2", "DNase", "DNase", "preferred")])
    files = pl.DataFrame(
        [
            bigwig("D1", "ENCSR2", "read-depth normalized signal", "1", True),
            bigwig(
                "D2",
                "ENCSR2",
                "read-depth normalized signal",
                "1",
                True,
                "/analyses/B/",
            ),
            bigwig("D3", "ENCSR2", "read-depth normalized signal", "2", False),
        ]
    )
    assert select_files(exps, files)["file"].to_list() == ["D1"]


def test_rna_gets_one_file_per_strand():
    exps = pl.DataFrame([experiment("ENCSR3", "RNA", "", "preferred")])
    files = pl.DataFrame(
        [
            bigwig("PLUS", "ENCSR3", "plus strand signal of unique reads", "1", True),
            bigwig("MINUS", "ENCSR3", "minus strand signal of unique reads", "1", True),
            bigwig("ALL", "ENCSR3", "plus strand signal of all reads", "1"),
        ]
    )
    chosen = select_files(exps, files).sort("suffix")
    assert chosen.select("file", "suffix").rows() == [
        ("MINUS", "Minus"),
        ("PLUS", "Plus"),
    ]


def test_two_candidates_raise():
    exps = pl.DataFrame([experiment("ENCSR4", "DNase", "DNase", "preferred")])
    files = pl.DataFrame(
        [
            bigwig("D1", "ENCSR4", "read-depth normalized signal", "1", True),
            bigwig("D2", "ENCSR4", "read-depth normalized signal", "2", True),
        ]
    )
    with pytest.raises(ValueError, match="Several files"):
        select_files(exps, files)


@pytest.mark.parametrize(
    "slims, fallback, name, expected",
    [
        ("epithelium,skin of body", "", "keratinocyte", "skin"),
        ("brain,blood", "", "x", "brain"),
        (None, "immune system,T cell", "T cell", "blood"),
        (None, "", "GM12878", "blood"),
        (None, "reproductive system", "fallopian tube", None),
        ("epithelium", "central nervous system", "x", "epithelium"),
    ],
)
def test_assign_organ(organs, slims, fallback, name, expected):
    assert assign_organ(slims, fallback, name, organs) == expected


def test_stem():
    assert stem("large intestine") == "largeIntestine"
    assert stem("skin") == "skin"

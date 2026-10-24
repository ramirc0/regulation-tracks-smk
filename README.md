# regulation-tracks-smk

Rebuilds the organ-average tracks of UCSC's ENCODE Regulation v4 track
(`wgEncodeReg4`) from a current ENCODE metadata snapshot. Each track is the
per-base mean of ENCODE signal bigWigs from one organ:

| Assay | Input file per experiment |
| --- | --- |
| DNase-seq | preferred read-depth normalized signal |
| ATAC-seq | fold change over control matching the preferred signal p-value |
| H3K4me3, H3K27ac (Histone and Mint-ChIP) | same as ATAC-seq |
| CTCF (TF ChIP-seq) | same as ATAC-seq |
| total RNA-seq, whole cell | preferred + and − strand signal of unique reads |

Inputs come from the default analysis of released GRCh38 experiments. Each
organ gets a track over all biosamples. It also gets a track over tissues and
primary cells when that set differs and is not empty. Missing data counts as
0, means are rounded to 3 decimals, and every base of the 195 no-alt GRCh38
contigs is written. The output is a one-file UCSC track hub.

## Layout

```
config/
  config.yaml.template   paths and run id; copy to config.yaml
  organs.tsv             organ vocabulary: priority, ENCODE slim terms, colors
workflow/
  Snakefile
  rules/                 common, select (freeze), average, hub
  scripts/               select_inputs.py, average_bigwigs.py, make_hub.py
  envs/pixi.toml         rule env
scripts/fetch_missing.py downloads selected files the mirror lacks
profiles/                local and SLURM
.tests/                  unit, integration and reference tests
```

Outputs go under `results/<run_id>/`:

```
selection/inputs.tsv     one row per input file
selection/tracks.tsv     one row per track
hub/hub.txt              the track hub
hub/organAve/<group>/    the tracks of one assay, strand and organ
```

## Organs

Each biosample type gets one organ: the first row of `config/organs.tsv` whose
`slims` match its ENCODE `organ_slims`. Without a match, the `fallback` terms
are checked against its system and cell slims, then the `pattern` against its
term name. Biosamples that match nothing are dropped and listed in the
`select` log.

## Running

```bash
cp config/config.yaml.template config/config.yaml   # then edit
pixi run snakemake --profile profiles/slurm freeze  # writes the selection
pixi run fetch --inputs results/<run_id>/selection/inputs.tsv \
    --metadata <snapshot>                            # only if files are missing
pixi run snakemake --profile profiles/slurm         # averages and the hub
```

`freeze` must run first: the selection decides which jobs exist. Rerun it
after changing the snapshot or `organs.tsv`.

Each `average` job reads every input of its group once with
[figwig](https://github.com/jmschrei/figwig), on 8 threads. It writes the
tracks to local disk, then copies them with `O_DSYNC`.

## Tests

```bash
pixi run unit        # scripts on tiny inputs
pixi run test        # whole workflow on toy data in .tests/integration
pixi run reference   # against the 2022 lists and hub files (lab cluster only)
pixi run lint
```

`.tests/integration/make_data.py` rebuilds the toy data from a snapshot.

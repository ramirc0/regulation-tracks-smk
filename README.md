# regulation-tracks-smk

A Snakemake workflow that rebuilds the organ-average tracks of UCSC's
[ENCODE Regulation v4](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=wgEncodeReg4)
track from current ENCODE metadata. The output is a one-file UCSC track hub.

Each track is the per-base mean of one assay's signal bigWigs from one organ:

| Assay | Input per experiment |
| --- | --- |
| DNase-seq | preferred read-depth normalized signal |
| ATAC-seq, H3K4me3, H3K27ac, CTCF | fold change over control, matching the preferred signal p-value |
| total RNA-seq (whole cell) | preferred + and − strand signal of unique reads |

Every organ gets a track over all biosamples. It also gets one over tissues
and primary cells when that set differs. Missing data counts as 0, means are
rounded to 3 decimals, and every base of the 195 no-alt GRCh38 contigs is
written. `config/organs.tsv` maps ENCODE organ slims to one organ per
biosample.

## Running

Needs [pixi](https://pixi.sh), an
[encode-metadata-smk](https://github.com/ramirc0/encode-metadata-smk) snapshot
and a local ENCODE mirror.

```bash
cp config/config.yaml.template config/config.yaml   # then edit
pixi run snakemake --profile profiles/slurm freeze  # select inputs
pixi run fetch --inputs results/<run_id>/selection/inputs.tsv \
    --metadata <snapshot>                            # only for files the mirror lacks
pixi run snakemake --profile profiles/slurm         # averages and the hub
```

Outputs land in `results/<run_id>/`: the selection tables in `selection/`,
and `hub/hub.txt` with its tracks in `hub/organAve/`.

## Tests

```bash
pixi run unit        # scripts on tiny inputs
pixi run test        # whole workflow on toy data
pixi run lint
```

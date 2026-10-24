# Organ averages: one job per assay, strand and organ.


rule average:
    input:
        inputs=f"{SELECTION}/inputs.tsv",
        tracks=f"{SELECTION}/tracks.tsv",
        chrom_sizes=config["chrom_sizes"],
        bigwigs=group_bigwigs,
    output:
        directory(f"{HUB}/organAve/{{group}}"),
    log:
        f"{LOGDIR}/average/{{group}}.txt",
    benchmark:
        f"{BENCHDIR}/average/{{group}}.tsv",
    threads: 8
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # numba caches are CPU-specific, so each job keeps its own. Tracks are
        # written to local disk, then copied with O_DSYNC: big buffered writes
        # to CephFS stall in memory-limited jobs.
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        export NUMBA_CACHE_DIR="$tmp/numba"

        python {workflow.basedir:q}/scripts/average_bigwigs.py \
            --inputs {input.inputs:q} \
            --tracks {input.tracks:q} \
            --group {wildcards.group:q} \
            --chrom-sizes {input.chrom_sizes:q} \
            --outdir "$tmp/out" \
            --threads {threads}

        mkdir -p {output:q}
        for f in "$tmp"/out/*.bw; do
            dd if="$f" of={output:q}/"${{f##*/}}" bs=64M oflag=dsync status=none
        done

        sync {output:q}
        """

# Input selection from the ENCODE metadata snapshot.


rule select:
    input:
        metadata=ancient(config["encode_metadata"]),
        organs=config["organs"],
    output:
        inputs=f"{SELECTION}/inputs.tsv",
        tracks=f"{SELECTION}/tracks.tsv",
    params:
        encode_root=config["encode_root"],
        download_dir=config["download_dir"],
    log:
        f"{LOGDIR}/select.txt",
    benchmark:
        f"{BENCHDIR}/select.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/select_inputs.py \
            --metadata {input.metadata:q} \
            --organs {input.organs:q} \
            --encode-root {params.encode_root:q} \
            --download-dir {params.download_dir:q} \
            --inputs {output.inputs:q} \
            --tracks {output.tracks:q}

        sync {output:q}
        """


rule freeze:
    input:
        f"{SELECTION}/inputs.tsv",
        f"{SELECTION}/tracks.tsv",

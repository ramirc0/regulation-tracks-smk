# One-file track hub over the averages.


rule hub:
    input:
        tracks=f"{SELECTION}/tracks.tsv",
        organs=config["organs"],
        groups=lambda w: expand(f"{HUB}/organAve/{{group}}", group=groups()),
    output:
        f"{HUB}/hub.txt",
    params:
        email=config["hub_email"],
    log:
        f"{LOGDIR}/hub.txt",
    benchmark:
        f"{BENCHDIR}/hub.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/make_hub.py \
            --tracks {input.tracks:q} \
            --organs {input.organs:q} \
            --email {params.email:q} \
            --output {output:q}

        sync {output:q}
        """

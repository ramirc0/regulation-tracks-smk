# Config, the frozen selection, and shared helpers.

from functools import cache
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config["run_id"]
OUTDIR = f"{RESULTS}/{RUN_ID}"
SELECTION = f"{OUTDIR}/selection"
HUB = f"{OUTDIR}/hub"
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")


wildcard_constraints:
    group=r"[A-Za-z0-9]+",


@cache
def frozen(name):
    """Read a table that `snakemake freeze` writes.

    These tables set which jobs exist, so they must exist before the DAG is
    built.

    Parameters
    ----------
    name : str
        Table name, "inputs" or "tracks".

    Returns
    -------
    polars.DataFrame
        The table.
    """
    path = Path(f"{SELECTION}/{name}.tsv")
    if not path.exists():
        raise WorkflowError(f"{path} is missing. Run `snakemake freeze` first.")
    return pl.read_csv(path, separator="\t")


def groups():
    """List the groups, one averaging job each.

    Returns
    -------
    list of str
        Group names, such as "skinDNase".
    """
    return frozen("tracks")["group"].unique().sort().to_list()


def group_bigwigs(wildcards):
    """List the input bigWigs of a group.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Wildcards with `group`.

    Returns
    -------
    list of str
        Input paths, flagged ancient so mirror timestamps never trigger reruns.
    """
    inputs = frozen("inputs").filter(pl.col("group") == wildcards.group)
    return [ancient(p) for p in inputs["path"]]

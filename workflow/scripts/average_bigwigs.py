"""Average the bigWigs of one group base by base and write its tracks.

A group is one assay, strand and organ. Each input is read once. It adds to
the sum over all biosamples and, for tissue and primary cell inputs, to the
tissue and primary cell sum. Each base of a track is the mean of its inputs,
rounded to `--decimals`, with missing data as 0. Every base of every contig
is written, zeros included.
"""

import argparse
from pathlib import Path

import numba
import numpy as np
import polars as pl
from figwig import BigWigReader, BigWigWriter


@numba.njit(parallel=True, nogil=True, cache=True)
def add(total, values):
    """Add `values` into `total` in place.

    Parameters
    ----------
    total : numpy.ndarray of float64
        Running sum.
    values : numpy.ndarray of float32
        Values of one input, as long as `total`.
    """
    for i in numba.prange(total.size):
        total[i] += values[i]


def build_parser():
    """Build the average_bigwigs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the average_bigwigs.py flags.
    """
    p = argparse.ArgumentParser(description="Average the bigWigs of one group.")
    p.add_argument("--inputs", required=True, help="Inputs TSV from select_inputs.py.")
    p.add_argument("--tracks", required=True, help="Tracks TSV from select_inputs.py.")
    p.add_argument("--group", required=True, help="Group to average.")
    p.add_argument("--chrom-sizes", required=True, help="Contigs to write, in order.")
    p.add_argument("--outdir", required=True, help="Directory for the group's tracks.")
    p.add_argument("--threads", type=int, default=8, help="Threads.")
    p.add_argument("--tile", type=int, default=100_000, help="Bases per read window.")
    p.add_argument("--decimals", type=int, default=3, help="Decimals to round to.")
    return p


def read_chrom_sizes(path):
    """Read a chrom sizes file.

    Parameters
    ----------
    path : str
        Two-column TSV of contig and length.

    Returns
    -------
    list of tuple of (str, int)
        Contigs and lengths, in file order.
    """
    with open(path) as fh:
        return [
            (c, int(n)) for c, n in (line.split()[:2] for line in fh if line.strip())
        ]


def average(paths, in_tp, sizes, writers, tile, threads, decimals):
    """Sum the inputs per contig and write each set's mean.

    Parameters
    ----------
    paths : list of str
        Input bigWigs, summed in this order.
    in_tp : list of bool
        Whether each input belongs to the tissue and primary cell set.
    sizes : list of tuple of (str, int)
        Contigs to write, in order.
    writers : dict of str to figwig.BigWigWriter
        Open writer of each set, "all" and optionally "tp".
    tile : int
        Bases per read window.
    threads : int
        Threads for reading.
    decimals : int
        Decimals to round the mean to.

    Raises
    ------
    ValueError
        If an input gives a contig a length other than `sizes`.
    """
    readers = [BigWigReader(p) for p in paths]
    length = dict(sizes)
    for path, reader in zip(paths, readers):
        for chrom, n in reader.chrom_sizes.items():
            if chrom in length and length[chrom] != n:
                raise ValueError(
                    f"{path}: {chrom} has {n} bp, expected {length[chrom]}"
                )
    n = {"all": len(paths), "tp": sum(in_tp)}
    longest = -(-max(length.values()) // tile)
    buf = np.empty((longest, tile), np.float32)
    for chrom, size in sizes:
        k = -(-size // tile)
        starts = np.arange(k, dtype=np.int64) * tile
        sums = {s: np.zeros(k * tile) for s in writers}
        for reader, tp in zip(readers, in_tp):
            if chrom not in reader.chrom_sizes:
                continue
            reader.read(chrom, starts, tile, out=buf[:k], n_jobs=threads)
            values = buf[:k].reshape(-1)
            add(sums["all"], values)
            if tp and "tp" in sums:
                add(sums["tp"], values)
        for s, total in sums.items():
            total /= n[s]
            np.round(total, decimals, out=total)
            total += 0.0  # turns -0.0 into 0.0
            mean = total.astype(np.float32)
            mean[size:] = np.nan
            writers[s].write(chrom, starts, mean.reshape(k, tile), missing=np.nan)


def main(argv=None):
    """Average one group and write its tracks.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    numba.set_num_threads(args.threads)
    inputs = pl.read_csv(args.inputs, separator="\t").filter(
        pl.col("group") == args.group
    )
    tracks = pl.read_csv(args.tracks, separator="\t").filter(
        pl.col("group") == args.group
    )
    if inputs.is_empty():
        raise ValueError(f"No inputs for group {args.group}")
    sizes = read_chrom_sizes(args.chrom_sizes)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    writers = {
        r["set"]: BigWigWriter(Path(args.outdir, r["file"]), sizes, n_jobs=args.threads)
        for r in tracks.iter_rows(named=True)
    }
    try:
        average(
            inputs["path"].to_list(),
            inputs["in_tp"].to_list(),
            sizes,
            writers,
            args.tile,
            args.threads,
            args.decimals,
        )
    finally:
        for w in writers.values():
            w.close()


if __name__ == "__main__":
    main()

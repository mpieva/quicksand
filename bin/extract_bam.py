#!/usr/bin/env python3

from argparse import ArgumentParser, RawTextHelpFormatter, FileType
from pysam import AlignmentFile


def cli():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        '-f', '--family', type=str,
        help="Taxonimical family for which to extract data"
    )
    parser.add_argument(
        '-k', '--krakenfile', type=FileType('r'),
        help="Output from kraken-translate"
    )
    parser.add_argument(
        '-c', '--compression-level', type=int, metavar='N',
        choices=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), default=6,
        help="compression level for output BAM files [0..9]\n"
             "(default: %(default)s)"
    )
    parser.add_argument(
        '-o', '--outfile', type=str,
        help="Output file"
    )
    parser.add_argument(
        'bamfile', type=FileType('rb'),
        help="BAM file from which to extract"
    )
    return parser.parse_args()


if __name__ == '__main__':

    cfg = cli()

    with cfg.krakenfile as kf:
        seq_ids = [l.split('\t', 1)[0] for l in kf
                   if "c__Mammalia" in l and f"f__{cfg.family}" in l]

    with AlignmentFile(cfg.bamfile, 'rb') as bf:
        bf.add_hts_options([f"level={cfg.compression_level}".encode('UTF-8')])
        with AlignmentFile(cfg.outfile, 'wb', template=bf) as of:
            for read in bf.fetch(until_eof=True):
                if read.query_name in seq_ids:
                    of.write(read)

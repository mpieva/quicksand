#!/usr/bin/env python3

"""Count duplicate sequences in a BAM file.

Output a BAM file with unique sequences. The number of times each sequence
appeared in the input file is encoded in the 'XP' field. Statistics may
optionally be written to a separate text output file.

"""

import sys
from argparse import ArgumentParser, FileType, RawTextHelpFormatter
from collections import Counter

from pysam import AlignmentFile, AlignedSegment


def cli():
    """Parse command-line arguments."""
    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        '-o', '--outputfile', metavar='FILE',
        type=FileType('wb'), nargs='?', default=None, const=sys.stdout,
        help="Write output in bam format (default: STDOUT)"
    )
    parser.add_argument(
        '-s', '--statsfile', metavar='FILE', type=FileType('w'),
        help="Write statistics to file"
    )
    parser.add_argument(
        '-c', '--cutoff', metavar='N', type=int, default=35,
        help='Cutoff length for calculating statistics (default: %(default)s)'
    )
    parser.add_argument(
        '-l', '--compression-level', type=int, metavar='N',
        choices=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), default=6,
        help="compression level for output BAM files [0..9]\n"
             "(default: %(default)s)"
    )
    parser.add_argument(
        'inputfile', nargs='?', type=FileType('rb'), default=sys.stdin,
        help="Input file (default: STDIN)"
    )
    return parser.parse_args()


if __name__ == '__main__':

    CFG = cli()

    # Build Counter object to contain sequence occurrence counts
    #
    with AlignmentFile(CFG.inputfile, 'rb') as inbam:
        hdr = inbam.header          # save BAM header
        ctr = Counter(
            a.get_forward_sequence() if a.is_reverse else a.query_sequence
            for a in inbam.fetch(until_eof=True)
        )

    # Build list of sequences above length cutoff
    #
    long_seqs = [s for s in ctr if len(s) >= CFG.cutoff]

    # Create BAM file with newly-generated entries for sequences.
    # Record occurrence count in 'XP' extra field.
    # (This is the "simplest thing that works". Alternatively one could
    # write one of the existing alignment records to the output BAM;
    # that would require a decision as to *which one* to write.)
    #
    if CFG.outputfile:
        with AlignmentFile(CFG.outputfile, 'wb', header=hdr) as outbam:
            outbam.add_hts_options([f"level={CFG.compression_level}".encode('UTF-8')])
            for n, seq in enumerate(long_seqs):
                s = AlignedSegment()
                s.query_name = f'seq_{n}'
                s.query_sequence = seq
                s.tags = (('XP', ctr[seq]),)
                outbam.write(s)

    # Write stats to stats file (if required)
    #
    if CFG.statsfile:

        # Calculate average occurrence count
        try:
            avg = sum(ctr[seq] for seq in long_seqs) / len(long_seqs)
        except ZeroDivisionError:
            avg = 0

        print(len(ctr), len(long_seqs), f"{avg:.1f}",
              sep='\t', file=CFG.statsfile)

#! /usr/bin/env python3

import pysam
import sys

def main(bamfile, doublestranded):

    #open the files
    infile = pysam.AlignmentFile(bamfile, 'rb')
    outfile = pysam.AlignmentFile(f"masked_{bamfile}", 'wb', template=infile)

    #main loop
    for read in infile:
        #alter quality scores
        qual = read.query_qualities
        seq = read.query_sequence
        for n in [0,1,2,-3,-2,-1]:
            #check if that position is deaminated
            if read.is_reverse == False:
                if n >= 0:
                    if seq[n]=='T':
                        qual[n] = 0
                else:
                    if (doublestranded and seq[n]=='A') or (not doublestranded and seq[n]=='T'):
                        qual[n] = 0
            else:
                if n >= 0:
                    if seq[n]=='A':
                        qual[n] = 0
                else:
                    if (doublestranded and seq[n]=='T') or (not doublestranded and seq[n]=='A'):
                        qual[n] = 0

        read.query_qualities = qual
        outfile.write(read)

    infile.close()
    outfile.close()

if __name__ == "__main__":
    bamfile = sys.argv[1]
    doublestranded = 'doublestranded' in sys.argv
    main(bamfile, doublestranded)



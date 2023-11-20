#! /usr/bin/env python3

import pysam
import sys

def main(bamfile):

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
            if (read.is_reverse == False and seq[n]=='T'):
                #set to 0
                qual[n] = 0
            elif(read.is_reverse and seq[n]=='A'):
                qual[n] = 0
            else:
                pass
        read.query_qualities = qual
        outfile.write(read)

    infile.close()
    outfile.close()

if __name__ == "__main__":
    bamfile = sys.argv[1]
    main(bamfile)



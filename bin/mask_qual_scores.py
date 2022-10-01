#! /usr/bin/env python3

import pysam
import sys

def main(bamfile):

    #open the files
    infile = pysam.AlignmentFile(bamfile, 'rb')
    outfile = pysam.AlignmentFile(f"{str(bamfile).rsplit('.')[0]}.masked.bam", 'wb', template=infile)
    
    #main loop
    for read in infile:
        #alter quality scores
        qual = read.query_qualities
        for n in [0,1,2,-3,-2,-1]:
            #set to 15 as mpileup by default ignores anything below 13
            qual[n] = 15
        read.query_qualities = qual
        outfile.write(read)

    infile.close()
    outfile.close()

if __name__ == "__main__":
    bamfile = sys.argv[1]
    main(bamfile)



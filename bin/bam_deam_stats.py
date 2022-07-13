#! /usr/bin/env python3

import pysam
import pandas as pd
from collections import defaultdict
import re
from pathlib import Path
from scipy import stats
import sys


def binomial_ci(x, n, alpha=0.05):
    #Clopper-Pearson interval for binomial distribution
    #x is number of successes, n is number of trials, alpha the percent
    if n==0:
        return "N/A","N/A"
    lower = stats.beta.interval(1-alpha, x,n-x+1)[0] if x != 0 else 0
    upper = stats.beta.interval(1-alpha, x+1,n-x)[1] if x != n else 1  
    return round(lower*100,2), round(upper*100,2)


def convert_reverse(df):
    df["SeqBP_corr"] = df.apply(lambda x: x[0].translate(str.maketrans('ACGT', 'TGCA') if x[2] else x[0]), axis=1)
    df["RefBP_corr"] = df.apply(lambda x: x[1].translate(str.maketrans('ACGT', 'TGCA') if x[2] else x[1]), axis=1)
    return df


def get_data(read):
    df = pd.DataFrame()
    #x,y,z --> x=index of read, y=index on reference, z=base (Upper = match, lower=substitution)
    #if x,y (and thus z) = none --> there is a deletion or insertion on read or reference respecively
    #remove data for gaps in read, since we need substitutions 
    tmp = [(x,y,z) for x,y,z in read.get_aligned_pairs(with_seq=True) if x != None]
    #reconstruct alignment (if ref = gap, count as match(no substitution))
    reference = "".join([z.upper() if None not in (x,y,z) else "-" for x,y,z in tmp])
    df["SeqBP"] = [read.seq[0]+read.seq[-1]]
    df["RefBP"] = [reference[0]+reference[-1]]
    df["Reverse"] = [read.is_reverse]
    df["Name"] = [read.query_name]
    return df


def import_data(bamfile):
    reads = pysam.AlignmentFile(open(bamfile))
    df = pd.DataFrame(columns=["SeqBP","RefBP","Reverse","Name"])
    for read in reads:
        df = df.append(get_data(read), ignore_index=True)
    reads.close()
    return df


def print_header():
    print("\t".join([
        "Ancient",
        "ReadsDeaminated",
        "Deam5",
        "Deam3",
        "Deam5Cond",
        "Deam3Cond",
    ]), file=sys.stdout)


def extract_reads(bamfile, df):
    infile = pysam.AlignmentFile(bamfile, 'rb')
    reads = set(df[(df['5deam1']==1) | (df['3deam1']==1)]['Name'])
    #make empty bam-file in case of 0 deaminated sequences 
    with open('output.deaminated.bam', 'wb') as header:
        header.write(pysam.view(bamfile,"-H","-b"))

    with pysam.AlignmentFile('output.deaminated.bam', 'wb', template=infile) as outfile:
        for read in infile:
            if read.query_name in reads:
                outfile.write(read)
    infile.close()
    return len(reads)

#
#
# MAIN
#
#

def main(bamfile):

    df = import_data(bamfile)
    if(len(df)==0):
        print_header()
        exit()
    df = convert_reverse(df)

    # Create all the stats 
    # Is there an easier way than hardcoding??
    # deam gets a 1 if there is a C-T substitution at the 5' or 3' end
    df["5deam1"] = df.apply(lambda x: 1 if x[4][0]=="T" and x[5][0]=="C" else 0,axis=1)
    df["3deam1"] = df.apply(lambda x: 1 if x[4][-1]=="T" and x[5][-1]=="C" else 0,axis=1)

    # cond gets a 1 if both ends are deaminated
    df["cond1"] = df.apply(lambda x: 1 if x[6]==1 and x[7]==1 else 0,axis=1)

    #subsets of possible deaminations
    deam_subset51 = df[df["RefBP_corr"].str.contains("^C")]
    deam_subset31 = df[df["RefBP_corr"].str.contains("C$")]

    #subsets of possible conditional substitions
    cond_subset51 = df[(df["RefBP_corr"].str.contains("^C"))&(df["3deam1"]==1)]
    cond_subset31 = df[(df["RefBP_corr"].str.contains("C$"))&(df["5deam1"]==1)]

    #Now calculate the percentages and other stats
    #percentage of deaminated C
    deam51 = round(sum(deam_subset51["5deam1"])/len(deam_subset51)*100,2) if len(deam_subset51)>0 else 'N/A'
    deam31 = round(sum(deam_subset31["3deam1"])/len(deam_subset31)*100,2) if len(deam_subset31)>0 else 'N/A'

    #and the confidence interval for that
    deam51_95ci = binomial_ci(sum(deam_subset51["5deam1"]),len(deam_subset51))
    deam31_95ci = binomial_ci(sum(deam_subset31["3deam1"]),len(deam_subset31))

    # percentage of deaminated conditional reads
    deam51cond = round(sum(cond_subset51["cond1"])/len(cond_subset51)*100,2) if len(cond_subset51)>0 else 'N/A'
    deam31cond = round(sum(cond_subset31["cond1"])/len(cond_subset31)*100,2) if len(cond_subset31)>0 else 'N/A'

    # the 95% confidence interval
    deam51cond_95ci = binomial_ci(sum(cond_subset51["cond1"]),len(cond_subset51))
    deam31cond_95ci = binomial_ci(sum(cond_subset31["cond1"]),len(cond_subset31))


    if (type(deam51_95ci[0]) != str and type(deam31_95ci[0]) != str):
        if (deam51_95ci[0] > 9.5 and deam31_95ci[0] > 9.5):
            ancient_string = "++"
        elif(deam51_95ci[0] > 9.5 or deam31_95ci[0] > 9.5):
            ancient_string = "+"
        else:
            ancient_string = "-"
    else:
        ancient_string = "-"

    #write the reads with deamination on terminal 5' or 3' end to file
    n_deam = extract_reads(bamfile,df)

    #header
    print_header()

    #and row
    print("\t".join([
    f"{ancient_string}",
    f"{n_deam}",
    f"{deam51}{deam51_95ci}".replace(" ",""),
    f"{deam31}{deam31_95ci}".replace(" ",""),
    f"{deam51cond}{deam51cond_95ci}".replace(" ",""),
    f"{deam31cond}{deam31cond_95ci}".replace(" ",""),
    ]
    ), file=sys.stdout)


if __name__ == "__main__":
    bamfile = sys.argv[1]
    main(bamfile)



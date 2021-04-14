#!/usr/bin/env python3

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
        return "NA","NA"
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
    #reconstruct alignment (if ref = gap, count as match)
    reference = "".join([z.upper() if None not in (x,y,z) else "-" for x,y,z in tmp])
    
    df["SeqBP"] = [read.seq[0:3]+read.seq[-3:]]
    df["RefBP"] = [reference[0:3]+reference[-3:]]
    df["Reverse"] = [read.is_reverse]
    df["MD"] = [read.get_tag("MD")]
    return df

def import_data(bamfile):
    reads = pysam.AlignmentFile(open(bamfile))
    df = pd.DataFrame(columns=["SeqBP","RefBP","Reverse"])
    for n,read in enumerate(reads):
        df = df.append(get_data(read), ignore_index=True)
    return df

#
#
# MAIN
#
#

bamfile = sys.argv[1]
rg = sys.argv[2]
family = sys.argv[3]
species = sys.argv[4]

df = import_data(bamfile)
df = convert_reverse(df)

#
# Create all the stats 
# Is there an easier way than hardcoding??
# deam gets a 1 if there is a C-T substitution at the 5' or 3' end
df["5deam1"] = df.apply(lambda x: 1 if x[4][0]=="T" and x[5][0]=="C" else 0,axis=1)
df["5deam2"] = df.apply(lambda x: 1 if x[4][1]=="T" and x[5][1]=="C" else 0,axis=1)
df["5deam3"] = df.apply(lambda x: 1 if x[4][2]=="T" and x[5][2]=="C" else 0,axis=1)
df["3deam3"] = df.apply(lambda x: 1 if x[4][3]=="T" and x[5][3]=="C" else 0,axis=1)
df["3deam2"] = df.apply(lambda x: 1 if x[4][4]=="T" and x[5][4]=="C" else 0,axis=1)
df["3deam1"] = df.apply(lambda x: 1 if x[4][5]=="T" and x[5][5]=="C" else 0,axis=1)

# cond gets a 1 if both ends are deaminated
df["cond1"] = df.apply(lambda x: 1 if x[6]==1 and x[11]==1 else 0,axis=1)
df["cond2"] = df.apply(lambda x: 1 if x[7]==1 and x[10]==1 else 0,axis=1)
df["cond3"] = df.apply(lambda x: 1 if x[8]==1 and x[9]==1 else 0,axis=1)

#subsets of possible deaminations
deam_subset51 = df[df["RefBP_corr"].str.contains("^C")]
deam_subset52 = df[df["RefBP_corr"].str.contains("^.C")]
deam_subset53 = df[df["RefBP_corr"].str.contains("^..C")]
deam_subset33 = df[df["RefBP_corr"].str.contains("C..$")]
deam_subset32 = df[df["RefBP_corr"].str.contains("C.$")]
deam_subset31 = df[df["RefBP_corr"].str.contains("C$")]

#subsets of possible conditional substitions
cond_subset51 = df[(df["RefBP_corr"].str.contains("^C"))&(df["3deam1"]==1)]
cond_subset52 = df[(df["RefBP_corr"].str.contains("^.C"))&(df["3deam2"]==1)]
cond_subset53 = df[(df["RefBP_corr"].str.contains("^..C"))&(df["3deam3"]==1)]
cond_subset33 = df[(df["RefBP_corr"].str.contains("C..$"))&(df["5deam3"]==1)]
cond_subset32 = df[(df["RefBP_corr"].str.contains("C.$"))&(df["5deam2"]==1)]
cond_subset31 = df[(df["RefBP_corr"].str.contains("C$"))&(df["5deam1"]==1)]

#Now calculate the percentages and other stats
#percentage of deaminated C
deam51 = round(sum(deam_subset51["5deam1"])/len(deam_subset51)*100,2) if len(deam_subset51)>0 else 0
deam52 = round(sum(deam_subset52["5deam2"])/len(deam_subset52)*100,2) if len(deam_subset52)>0 else 0
deam53 = round(sum(deam_subset53["5deam3"])/len(deam_subset53)*100,2) if len(deam_subset53)>0 else 0
deam33 = round(sum(deam_subset33["3deam3"])/len(deam_subset33)*100,2) if len(deam_subset33)>0 else 0
deam32 = round(sum(deam_subset32["3deam2"])/len(deam_subset32)*100,2) if len(deam_subset32)>0 else 0
deam31 = round(sum(deam_subset31["3deam1"])/len(deam_subset31)*100,2) if len(deam_subset31)>0 else 0
#and the confidence interval for that
deam51_95ci = binomial_ci(sum(deam_subset51["5deam1"]),len(deam_subset51))
deam52_95ci = binomial_ci(sum(deam_subset52["5deam2"]),len(deam_subset52))
deam53_95ci = binomial_ci(sum(deam_subset53["5deam3"]),len(deam_subset53))
deam33_95ci = binomial_ci(sum(deam_subset33["3deam3"]),len(deam_subset33))
deam32_95ci = binomial_ci(sum(deam_subset32["3deam2"]),len(deam_subset32))
deam31_95ci = binomial_ci(sum(deam_subset31["3deam1"]),len(deam_subset31))
# the number of Cs that could have been deaminated
nRef51 = len(deam_subset51)
nRef52 = len(deam_subset52)
nRef53 = len(deam_subset53)
nRef33 = len(deam_subset33)
nRef32 = len(deam_subset32)
nRef31 = len(deam_subset31)
# percentage of deaminated conditional reads
deam51cond = round(sum(cond_subset51["cond1"])/len(cond_subset51)*100,2) if len(cond_subset51)>0 else 0
deam52cond = round(sum(cond_subset52["cond2"])/len(cond_subset52)*100,2) if len(cond_subset52)>0 else 0
deam53cond = round(sum(cond_subset53["cond3"])/len(cond_subset53)*100,2) if len(cond_subset53)>0 else 0
deam33cond = round(sum(cond_subset33["cond3"])/len(cond_subset33)*100,2) if len(cond_subset33)>0 else 0
deam32cond = round(sum(cond_subset32["cond2"])/len(cond_subset32)*100,2) if len(cond_subset32)>0 else 0
deam31cond = round(sum(cond_subset31["cond1"])/len(cond_subset31)*100,2) if len(cond_subset31)>0 else 0
# the 95% confidence interval
deam51cond_95ci = binomial_ci(sum(cond_subset51["cond1"]),len(cond_subset51))
deam52cond_95ci = binomial_ci(sum(cond_subset52["cond2"]),len(cond_subset52))
deam53cond_95ci = binomial_ci(sum(cond_subset53["cond3"]),len(cond_subset53))
deam33cond_95ci = binomial_ci(sum(cond_subset33["cond3"]),len(cond_subset33))
deam32cond_95ci = binomial_ci(sum(cond_subset32["cond2"]),len(cond_subset32))
deam31cond_95ci = binomial_ci(sum(cond_subset31["cond1"]),len(cond_subset31))
# the number of conditional Cs that could have been deaminated
nRef51cond = len(cond_subset51)
nRef52cond = len(cond_subset52)
nRef53cond = len(cond_subset53)
nRef33cond = len(cond_subset33)
nRef32cond = len(cond_subset32)
nRef31cond = len(cond_subset31)

if(deam51_95ci[0] > 9.5 and deam31_95ci[0] > 9.5):
    ancient_string = "++"
elif(deam51_95ci[0] > 9.5 or deam31_95ci[0] > 9.5):
    ancient_string = "+"
else:
    ancient_string = "-"

#header
print("\t".join([
    "RG",
    "Ancient",
    "Family",
    "Species",
    "Freq51(95CI)",
    "Freq52(95CI)",
    "Freq53(95CI)",
    "Freq33(95CI)",
    "Freq32(95CI)",
    "Freq31(95CI)",
    "nRef51",
    "nRef52",
    "nRef53",
    "nRef33",
    "nRef32",
    "nRef31",
    "Freq51cond(95CI)",
    "Freq52cond(95CI)",
    "Freq53cond(95CI)",
    "Freq33cond(95CI)",
    "Freq32cond(95CI)",
    "Freq31cond(95CI)",
    "nRef51cond",
    "nRef52cond",
    "nRef53cond",
    "nRef33cond",
    "nRef32cond",
    "nRef31cond",
]
),file=sys.stdout)

#and row
print("\t".join([
    f"{rg}",
    f"{ancient_string}",
    f"{family}",
    f"{species}",
    f"{deam51} {deam51_95ci}",
    f"{deam52} {deam52_95ci}",
    f"{deam53} {deam53_95ci}",
    f"{deam33} {deam33_95ci}",
    f"{deam32} {deam32_95ci}",
    f"{deam31} {deam31_95ci}",
    f"{nRef51}",
    f"{nRef52}",
    f"{nRef53}",
    f"{nRef33}",
    f"{nRef32}",
    f"{nRef31}",
    f"{deam51cond} {deam51cond_95ci}",
    f"{deam52cond} {deam52cond_95ci}",
    f"{deam53cond} {deam53cond_95ci}",
    f"{deam33cond} {deam33cond_95ci}",
    f"{deam32cond} {deam32cond_95ci}",
    f"{deam31cond} {deam31cond_95ci}",
    f"{nRef51cond}",
    f"{nRef52cond}",
    f"{nRef53cond}",
    f"{nRef33cond}",
    f"{nRef32cond}",
    f"{nRef31cond}",
    ]
), file=sys.stdout)

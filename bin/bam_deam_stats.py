#! /usr/bin/env python3

import pysam
from scipy import stats
import sys

def binomial_ci(x, n, alpha=0.05):
    #Clopper-Pearson interval for binomial distribution
    #x is number of successes, n is number of trials, alpha the percent
    if n==0:
        return None,None
    lower = stats.beta.interval(1-alpha, x,n-x+1)[0] if x != 0 else 0
    upper = stats.beta.interval(1-alpha, x+1,n-x)[1] if x != n else 1  
    return round(lower*100,2), round(upper*100,2)


def main(bamfile, stats_only=False):
    #store all reference bases
    all_first = []
    all_last = []

    #store all bases from the other end of the observed deamination
    cond5_last = []
    cond3_first = []
    
    n_deam_1 = 0
    n_deam_3 = 0
    #count the individual events
    #first or last base deaminated
    n_deam51 = 0
    n_deam31 = 0
    #first or last 3 bases deaminated
    n_deam53 = 0
    n_deam33 = 0
    #first and last base deaminated
    n_cond = 0

    #open the files
    infile = pysam.AlignmentFile(bamfile, 'rb')
    
    if not only_stats:
        out1term = pysam.AlignmentFile('output.deaminated1.bam', 'wb', template=infile)
        out3term = pysam.AlignmentFile('output.deaminated3.bam', 'wb', template=infile)
    
    #main loop
    for read in infile:
        # get read and ref sequences, reverse complement if necessary
        seq = read.query_sequence
        ref = ''.join([str(z).replace('None','-') for x,y,z in read.get_aligned_pairs(with_seq=True) if x != None])

        if read.is_reverse:
            ref = "".join([x.translate(str.maketrans('NACGTnacgt-', 'NTGCAntgca-')) for x in ref])[::-1]
            seq = read.get_forward_sequence()

        rlen = len(ref)-1 #0based required
 
        all_first.append(ref[0].upper())
        all_last.append(ref[-1].upper())

        # check if ct substitutions in the first or last 3 postions and save indices
        mism = [
            n for n,x in enumerate(ref)                             
            if ((x, seq[n]) == ('c','T')) and (n<=2 or n>=rlen-2)   
        ]

        # skip non-damaged reads
        if len(mism) == 0:
            continue
        
        # increase the deamination stats counts
         
        deam51 = 0 in mism
        deam31 = rlen in mism
        deam53 = any(x<=2 for x in mism) 
        deam33 = any(x>=rlen-2 for x in mism)
        cond = deam51 and deam31

        n_deam51 += int(deam51)
        n_deam31 += int(deam31)
        n_deam53 += int(deam53)
        n_deam33 += int(deam33)
        n_cond   += int(cond) 

        if deam51:
            cond5_last.append(ref[-1].upper())
        if deam31:
            cond3_first.append(ref[0].upper())

        #write read to the file(s)
        if deam51 or deam31:
            n_deam_1 += 1
            if not only_stats:
                out1term.write(read)
        if deam53 or deam33:
            n_deam_3 += 1
            if not only_stats:
                out3term.write(read)

    infile.close()
    if not only_stats:
        out1term.close()
        out3term.close()

    ## Calculate the stats

    all_5c = all_first.count('C')
    all_3c = all_last.count('C')
    all_5c_cond3 = cond5_last.count('C')
    all_3c_cond5 = cond3_first.count('C')

    #calculate deamination percentage
    p_deam51 = round(n_deam51/all_5c * 100,2) if all_5c>0 else None
    p_deam31 = round(n_deam31/all_3c * 100,2) if all_3c>0 else None
    
    #and the confidence interval for that
    p_deam51_95ci = binomial_ci(n_deam51,all_5c)
    p_deam31_95ci = binomial_ci(n_deam31,all_3c)

    # percentage of deaminated conditional reads
    p_deam51_cond = round(n_cond/all_3c_cond5 *100,2) if all_3c_cond5>0 else None
    p_deam31_cond = round(n_cond/all_5c_cond3 *100,2) if all_5c_cond3>0 else None
    
    # the 95% confidence interval
    p_deam51_cond_95ci = binomial_ci(n_cond,all_3c_cond5)
    p_deam31_cond_95ci = binomial_ci(n_cond,all_5c_cond3)
    
    ancientness = '-'
    test51,test31 = p_deam51_95ci[0],p_deam31_95ci[0]
    if test51 and test31:
        if test51 > 9.5 or test31 > 9.5:
            ancientness = '+'
        if test51 > 9.5 and test31 > 9.5:
            ancientness = '++'

    #And the report
    #print header
    print(
        "Ancientness", "ReadsDeam(1term)","ReadsDeam(3term)", 
        "Deam5(95ci)", "Deam3(95ci)", "Deam5Cond(95ci)",
        "Deam3Cond(95ci)", sep='\t',file=sys.stdout
        )

    #and row
    print(
        ancientness, n_deam_1, n_deam_3, 
        f"{p_deam51} {p_deam51_95ci}".replace(", ",",").replace("None","N/A"),
        f"{p_deam31} {p_deam31_95ci}".replace(", ",",").replace("None","N/A"), 
        f"{p_deam51_cond} {p_deam51_cond_95ci}".replace(", ",",").replace("None","N/A"),
        f"{p_deam31_cond} {p_deam31_cond_95ci}".replace(", ",",").replace("None","N/A"), 
        sep='\t', file=sys.stdout
        )


if __name__ == "__main__":
    bamfile = sys.argv[1]
    only_stats = False
    try:
        if sys.argv[2]=='only_stats':
            only_stats = True
    except IndexError:
        pass
    main(bamfile, only_stats)



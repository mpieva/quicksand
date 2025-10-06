#! /usr/bin/env python3

import pysam
from scipy import stats
import numpy as np
import sys
import statistics

def binomial_ci(x, n, alpha=0.05):
    #Clopper-Pearson interval for binomial distribution
    #x is number of successes, n is number of trials, alpha the percent
    if n==0:
        return "N/A,N/A"
    lower = stats.beta.interval(1-alpha, x,n-x+1)[0] if x != 0 else 0
    upper = stats.beta.interval(1-alpha, x+1,n-x)[1] if x != n else 1
    return f"{round(lower*100,1):.1f},{round(upper*100,1):.1f}"

def substitution_to_list(reference, sequence, from_char, to_char):
    """
    Compares two strings character by character and returns a list where:
    - None  : the reference character is not ref_char
    - 0     : the reference character is ref_char and matches the test character
    - 1     : the reference character is ref_char and the test character matches test_char
    """
    ref_arr = np.array(list(reference.upper()))
    test_arr = np.array(list(sequence.upper()))

    result = np.full_like(ref_arr, None, dtype=object)  # default None
    mask_ref = ref_arr == from_char
    result[mask_ref] = np.where(test_arr[mask_ref] == to_char, 1,
                                np.where(test_arr[mask_ref] == from_char, 0, None))
    
    # I need to reduce the array to 20bp (10 on each site)
    # to make 2D-calculations later!
    first_last = []
    first_last.extend(result[:10]) 
    first_last.extend(result[-10:])                           
    return first_last

def summarize_arrays_np(arrays):
    """
    Summarize a list of arrays with values [None, 0, 1] by computing,
    for the first 10 and last 10 positions, the fraction of 1's among non-None values.
    Returns a 20-element summary array.
    """
    # Convert to NumPy array, using object dtype to handle None
    # Ensure all arrays are lists and form a 2D object array
    arr = np.array(arrays, dtype=object)

    # Function to compute ratio of 1's / non-None for a slice
    mask = arr != None
    ones = (arr == 1) & mask
    counts_ones = ones.sum(axis=0)
    counts_non_none = mask.sum(axis=0)

    # Avoid division by zero
    percentages = np.where(counts_non_none > 0, np.round((counts_ones / counts_non_none)*100, 3), None)
    
    # Confidence intervals (position-wise using list comprehension)
    ci = [
        binomial_ci(x, n) if n > 0 else None
        for x, n in zip(counts_ones, counts_non_none)
    ]

    return percentages, ci

def main(bamfile):
    #store all subsitutions of all sequences in separate arrays
    #these are the main arrays
    C_T = []
    C_G = []
    C_A = []
    G_C = []
    G_T = []
    G_A = []

    #open the files
    infile = pysam.AlignmentFile(bamfile, 'rb')

    #
    # main loop
    #

    # now iterate the bamfile
    for read in infile:
        # get read and ref sequences, reverse complement if necessary
        seq = read.query_sequence
        ref = ''.join([str(z).replace('None','-') for x,y,z in read.get_aligned_pairs(with_seq=True) if x != None])

        if read.is_reverse:
            ref = "".join([x.translate(str.maketrans('NACGTnacgt-', 'NTGCAntgca-')) for x in ref])[::-1]
            seq = read.get_forward_sequence()


        C_T.append(substitution_to_list(ref,seq,'C','T'))
        C_G.append(substitution_to_list(ref,seq,'C','G'))
        C_A.append(substitution_to_list(ref,seq,'C','A'))
        G_T.append(substitution_to_list(ref,seq,'G','T'))
        G_C.append(substitution_to_list(ref,seq,'G','C'))
        G_A.append(substitution_to_list(ref,seq,'G','A'))

    ## Write the stats to file

    with open('substitutions.tsv', 'w') as outfile:
        print('Sub', '\t'.join([str(x) for x in range(1,11)]), '\t'.join([str(-1*x) for x in range(10,0,-1)]), sep='\t', file=outfile)
        print('C->T', '\t'.join([str(x) for x in summarize_arrays_np(C_T)[0]]), sep='\t', file=outfile)
        print('C->G', '\t'.join([str(x) for x in summarize_arrays_np(C_G)[0]]), sep='\t', file=outfile)
        print('C->A', '\t'.join([str(x) for x in summarize_arrays_np(C_A)[0]]), sep='\t', file=outfile)
        print('G->T', '\t'.join([str(x) for x in summarize_arrays_np(G_T)[0]]), sep='\t', file=outfile)
        print('G->C', '\t'.join([str(x) for x in summarize_arrays_np(G_C)[0]]), sep='\t', file=outfile)
        print('G->A', '\t'.join([str(x) for x in summarize_arrays_np(G_A)[0]]), sep='\t', file=outfile)

    with open('confidence.tsv', 'w') as outfile:
        print('Sub', '\t'.join([str(x) for x in range(1,11)]), '\t'.join([str(-1*x) for x in range(10,0,-1)]), sep='\t', file=outfile)
        print('C->T', '\t'.join([str(x) for x in summarize_arrays_np(C_T)[1]]), sep='\t', file=outfile)
        print('C->G', '\t'.join([str(x) for x in summarize_arrays_np(C_G)[1]]), sep='\t', file=outfile)
        print('C->A', '\t'.join([str(x) for x in summarize_arrays_np(C_A)[1]]), sep='\t', file=outfile)
        print('G->T', '\t'.join([str(x) for x in summarize_arrays_np(G_T)[1]]), sep='\t', file=outfile)
        print('G->C', '\t'.join([str(x) for x in summarize_arrays_np(G_C)[1]]), sep='\t', file=outfile)
        print('G->A', '\t'.join([str(x) for x in summarize_arrays_np(G_A)[1]]), sep='\t', file=outfile)

if __name__ == "__main__":
    bamfile = sys.argv[1]   
    #TODO: 
    # - output also the conditional subsitutions??
    #   - thats not easy, because in double-stranded context the condition is on G->A for C->T
    #   - and I dont know HOW exactly. Condition on the terminal base? better to leave it as it is right now...

    main(bamfile)
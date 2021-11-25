#! /usr/bin/env python3 
import sys
from collections import defaultdict

def find_best(tree):
    return max(tree.keys(), key=lambda x: tree[x]['kmers'])

results = {}

sp = 'placeholder'


min_kmers = int(sys.argv[2])
min_reads = int(sys.argv[3])

#parse the report
for row in open(sys.argv[1],'r'):
    try:
        perc, reads, taxReads, kmers, dup, cov, taxid, level, name = row.split('\t', 8)
    except ValueError: #headerline in krakenUniq
        continue
    if level == 'order':
        order = name.strip()
    elif level == 'family':
        fam = name.strip()
        results[(fam,order)] = {
                'id':taxid,'kmers':kmers, 'cov':cov,
                'dup':dup,'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'genus':
        gen = name.strip()
        results[(fam,order)]['tree'][gen] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'species':
        sp = name.strip()
        results[(fam,order)]['tree'][gen]['tree'][sp] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'subspecies':
        ssp = name.strip()
        results[(fam,order)]['tree'][gen]['tree'][sp]['tree'][ssp] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers)
        }
    else:
        continue

real_results = {}

for (fam,order) in results:
    if(results[(fam,order)]["counts"]<min_reads or results[(fam,order)]["kmers"]<min_kmers):
        continue
    try:
        gen = find_best(results[(fam,order)]['tree'])
        try:
            sp = find_best(results[(fam,order)]['tree'][gen]['tree'])
            try:
                ssp = find_best(results[(fam,order)]['tree'][gen]['tree'][sp]['tree'])
                real_results[(fam,order)] = results[(fam,order)]['tree'][gen]['tree'][sp]['tree'][ssp]['id']
            except:
                real_results[(fam,order)] = results[(fam,order)]['tree'][gen]['tree'][sp]['id']
        except:
            real_results[(fam,order)] = results[(fam,order)]['tree'][gen]['id']
    except:
        real_results[(fam,order)] = results[(fam,order)]['id']

with open('parsed_record.tsv', 'w') as outfile:
    print('Family','Order','BestTaxID','FamReads','FamKmers','FamKmerCov','FamKmerDup', sep='\t', file=outfile)
    for (fam,order),taxid  in real_results.items():
        print(fam, order, taxid,
              results[(fam,order)]['counts'],
              results[(fam,order)]['kmers'], 
              results[(fam,order)]['cov'], 
              results[(fam,order)]['dup'], 
              sep='\t', file=outfile)


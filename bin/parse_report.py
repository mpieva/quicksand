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
        _, reads, _, kmers, dup, cov, taxid, level, name = row.split('\t', 8)
    except ValueError: #headerline in krakenUniq
        continue
    if level == 'family':
        fam = name.strip()
        results[fam] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'genus':
        gen = name.strip()
        results[fam]['tree'][gen] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'species':
        sp = name.strip()
        results[fam]['tree'][gen]['tree'][sp] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers), 'tree':{}
        }
    elif level == 'subspecies':
        ssp = name.strip()
        results[fam]['tree'][gen]['tree'][sp]['tree'][ssp] = {
            'id':taxid, 'counts':int(reads), 'kmers':float(kmers)
        }
    else:
        continue

real_results = {}

for fam in results:
    if(results[fam]["counts"]<min_reads or results[fam]["kmers"]<min_kmers):
        continue
    try:
        gen = find_best(results[fam]['tree'])
        try:
            sp = find_best(results[fam]['tree'][gen]['tree'])
            try:
                ssp = find_best(results[fam]['tree'][gen]['tree'][sp]['tree'])
                real_results[fam] = results[fam]['tree'][gen]['tree'][sp]['tree'][ssp]['id']
            except:
                real_results[fam] = results[fam]['tree'][gen]['tree'][sp]['id']
        except:
            real_results[fam] = results[fam]['tree'][gen]['id']
    except:
        real_results[fam] = results[fam]['id']

with open('parsed_record.tsv', 'w') as outfile:
    for fam,taxid  in real_results.items():
        print(fam, taxid, sep='\t', file=outfile)


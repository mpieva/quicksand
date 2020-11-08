#! /usr/bin/env python3 
import sys
from collections import defaultdict

def find_best(tree):
    return max(tree.keys(), key=lambda x: tree[x]['counts'])

results = {}
sp = 'placeholder'

for row in open(sys.argv[1],'r'):
    _, reads, _, level, taxid, name = row.split('\t', 5)
    if level == 'F':
        fam = name.strip()
        results[fam] = {'id':taxid, 'counts':int(reads), 'tree':{}}
    elif level == 'G':
        gen = name.strip()
        results[fam]['tree'][gen] = {
            'id':taxid, 'counts':int(reads), 'tree':{}
        }
    elif level == 'S':
        sp = name.strip()
        results[fam]['tree'][gen]['tree'][sp] = {
            'id':taxid, 'counts':int(reads), 'tree':{}
            }
    elif level == '-' and name.strip().startswith(sp):
        ssp = name.strip()
        results[fam]['tree'][gen]['tree'][sp]['tree'][ssp] = {
            'id':taxid, 'counts':int(reads)
            }
    else:
        continue

real_results = {}

for fam in results:
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


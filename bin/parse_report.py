#! /usr/bin/env python3 
import sys

report = sys.argv[1]
min_kmers = int(sys.argv[2])
min_reads = int(sys.argv[3])

# parse the report
# Assume the following structure:
#
# family(400kmers)
# |  |  
# |  +- genus1 (300kmers)
# +--genus2 (5 kmers)
#    |
#    +-species1 (2 Kmers)
#
# The "easy" approach would be to take the lowest assigned node (species1)
# However, its more likely, that the read belongs to (a species from) genus1
# so we need to find in each node the path with the highest remaining kmer-count
# For the bwa-step later, we use only order, family, genus, species, subspecies
# So ignore the intermediate clades!

class Node:
    def __init__ (self,parent,name,level,taxid,coverage,kmers,dup,reads):
        self.parent = parent
        self.name = name
        self.level = level
        self.taxid = taxid
        self.coverage = coverage
        self.kmers = int(kmers)
        self.dup = float(dup)
        self.reads = int(reads)
        self.children = []

    @property
    def order(self):
        if self.level == 'order':
            return self.name
        else:
            return self.parent.order if self.parent else None

    @property
    def best_child(self):
        return  max(self.children, key=lambda x: x.kmers).best_child if len(self.children) > 0  else self


hierarchie = ['order','family','genus','species','subspecies']

parents = {
    'order':None,    
    'family':None,
    'genus':None,
    'species':None,
    'subspecies':None
    }

family_nodes = []

for row in open(report,'r'):
    try:
        perc, reads, taxReads, kmers, dup, cov, taxid, level, name = (x.strip() for x in row.split('\t', 8))
    except ValueError: #headerline in krakenUniq
        continue

    if level not in hierarchie:
        continue

    parent = parents[level]
    node = Node(parent,name,level,taxid,cov,kmers,dup,reads)
    if parent!=None:
        parent.children.append(node)

    for lower in hierarchie[hierarchie.index(level)+1:]:
        parents[lower] = node

    if level=='family':
        family_nodes.append(node)


with open('parsed_record.tsv', 'w') as outfile:
    print('Family','Order','BestTaxID','FamReads','FamilyKmers','KmerCoverage','KmerDupRate', sep='\t', file=outfile)
    for node in family_nodes:
        best = node.best_child
        order = node.order

        ## Now apply the filters ##
        if node.reads < min_reads or node.kmers < min_kmers:
            continue

        print(node.name,order,best.taxid,node.reads,node.kmers,node.coverage,node.dup, sep='\t', file=outfile)

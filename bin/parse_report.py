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
    def best_child(self):
        """Find the direct best children"""
        return max(self.children, key=lambda x: x.kmers).best_child if len(self.children) > 0  else self


order_nodes = []

for row in open(report,'r'):
    try:
        perc, reads, taxReads, kmers, dup, cov, taxid, level, name = (x.strip() for x in row.split('\t', 8))
    except ValueError: #headerline in krakenUniq
        continue

    if level == 'order':
        current_node = Node(None,name,level,taxid,cov,kmers,dup,reads)
        order_nodes.append(current_node)

    elif level in ['family', 'genus','species','subspecies']:
        node = Node(current_node, name,level,taxid,cov,kmers,dup,reads)
        current_node.children.append(node)
        current_node = node

    else:
        continue


with open('parsed_record.tsv', 'w') as outfile:
    print('Family','Order','BestTaxID','FamReads','FamKmers','FamKmerCov','FamKmerDup', sep='\t', file=outfile)
    for node in order_nodes:
        for fam_node in node.children:
            if fam_node.level != 'family':
                family = 'Undefined'
            else:
                family = fam_node.name
            best = fam_node.best_child

            ## Now apply the filters ##
            if fam_node.reads < min_reads or fam_node.kmers < min_kmers:
                continue

            print(family, node.name, best.taxid, fam_node.reads, fam_node.kmers, fam_node.coverage, fam_node.dup, sep='\t', file=outfile)

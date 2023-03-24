#! /usr/bin/env pypy3

import sys
import json

def get_names_dict(names_file):
    names_dict = {}
    
    with open(names_file) as f:
        for line in f:
            if 'scientific name' in line:
                taxid,name,_ = line.strip().split("\t|\t", 2)
                names_dict[int(taxid)]=name
    return names_dict

def get_nodes_dict(taxonomy_file):
    taxonomy_dict = {}
    
    # Read the taxonomy file and create a dictionary of taxonomic IDs and their parent and ranks
    with open(taxonomy_file) as f:
        for line in f:
            fields = line.strip().split('\t|\t')
            taxid = int(fields[0])
            parent = int(fields[1].strip())
            rank = fields[2].strip()

            # create or update entry
            if taxid in taxonomy_dict:
                taxonomy_dict[taxid].update({'rank':rank, 'parent':parent})
            else:
                taxonomy_dict[taxid] = {'rank': rank, 'parent': parent, 'children': []}
            
            # add the children link
            if parent in taxonomy_dict:
                taxonomy_dict[parent]['children'].append(taxid)
            else:
                taxonomy_dict[parent] = {'children':[taxid]}
    return taxonomy_dict

def parse_content(seqid2taxid):
    genomes = []
    with open(seqid2taxid) as f:
        for line in f:
            seqid, taxid = line.split('\t',1)
            genomes.append(int(taxid))
    return set(genomes)

def traverse_tree(taxid, taxdict, genomes):
    results = []
    record = taxdict[taxid]
    if taxid in genomes:
        results.append(taxid)
    for child in record['children']:
        if child != taxid:
            results.extend( traverse_tree(child, taxdict, genomes) )
    return results
            

if __name__ == '__main__':
    nodefile = sys.argv[1]
    namesfile = sys.argv[2]
    seqid2taxid = sys.argv[3]

    names = get_names_dict(namesfile)
    nodes = get_nodes_dict(nodefile)
    genomes = parse_content(seqid2taxid)
    
    final_json = {}
    for taxid in nodes:
        tmp = [names[x].replace(' ','_') for x in traverse_tree(taxid, nodes, genomes)]
        if len(tmp)>0:
            final_json[taxid] = tmp
    json.dump(final_json, sys.stdout)



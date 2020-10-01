#!/usr/bin/env python3
import sys
import re

taxa = {
    tuple(
    	[y for y in x.replace("\n","").split("\t")[0].split("|") if y[0] in "fgs"]
    ): int(x.replace("\n",'').split("\t")[1]) 
    for x in open(sys.argv[1])
}
families = {entry[0]:{"fam_count":taxa[entry], "gen_count":0, "sp_count":0,
	"genus":None, "species":None} for entry in taxa.keys() if len(entry) == 1}

for fam,gen in [x for x in taxa if len(x)==2] : 
     if (int(families[fam]["gen_count"]) < taxa[(fam, gen)]) or (families[fam]['genus']==None): 
         families[fam]["gen_count"]=taxa[(fam, gen)] 
         families[fam]["genus"]=gen.split("__")[1] 

for fam,gen,sp in [x for x in taxa if len(x)==3] : 
     if (int(families[fam]["sp_count"]) < taxa[(fam, gen,sp)]) or (families[fam]['species']==None): 
         families[fam]["sp_count"]=taxa[(fam, gen, sp)] 
         families[fam]["species"]=re.search(r"\w*?_..", sp.split("__")[1]).group() 

# This is a minimum 3 reads per family filter --> Hardcoded :/
families = {k:v for k,v in families.items() if v["fam_count"] > 3}

for record in families:
	if families[record]['sp_count']:
		with open(f"s_{families[record]['species']}__{record.split('__')[1]}.txt", 'w') as outfile:
			outfile.write("touch")
	elif families[record]['gen_count']:
		with open(f"g_{families[record]['genus']}__{record.split('__')[1]}.txt", 'w') as outfile:
			outfile.write("touch")
	else:
		with open(f"f_{record.split('__')[1]}__{record.split('__')[1]}.txt", 'w') as outfile:
			outfile.write("touch")

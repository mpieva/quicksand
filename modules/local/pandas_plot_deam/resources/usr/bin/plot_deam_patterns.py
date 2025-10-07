#! /usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

df = pd.read_csv(sys.argv[1], sep='\t')
title = sys.argv[2]

df = df.melt(id_vars=['Sub'])

bases = list(set(df['Sub']))
palette = {x:'lightgrey' for x in bases}
palette.update({'C->T':'red', "G->A":'black'})

#now plot the rate
# fix the positions first

df['Position in Sequence'] = df.apply(lambda x: int(x.variable) if int(x.variable) > 0 else int(x.variable)+21, axis=1)
df['Frequency'] = df['value'].apply(lambda x: float(x))
df['Substitution'] = df['Sub']

g = sns.lineplot(
    data = df, 
    x='Position in Sequence',
    y='Frequency', 
    hue='Substitution',
    palette=palette
)

ticks = list(range(0,11))
ticks.extend(list(range(-10,0)))

# Update the x-ticks
g.set_ylabel("Frequency [%]")
g.set_xticks(list(range(0,21)))
g.set_xticklabels(ticks)
g.axvline(x=10.5, ymax=0.5, ls="--", color='grey')
g.set_title(title)

plt.tight_layout()
plt.savefig(f"{title}_plot.jpg", dpi=300)
import numpy
import matplotlib.pyplot as plt
import pandas as pd

colors = {'blue': tuple(x/256 for x in (46,37,133)),
          'red': tuple(x/256 for x in (194,106,119)),
          'lgreen': tuple(x/256 for x in (93,168,153)),
          'gold': tuple(x/256 for x in (220,205,125)),
          'green': tuple(x/256 for x in (51, 117,56)),
          'lblue': tuple(x/256 for x in (148,203,236)),
          'magenta': tuple(x/256 for x in (159,74,150)),
          'wine': tuple(x/256 for x in (126, 41, 84)), 
          'black': tuple(x/256 for x in (0,0,0)), 
}

x = pd.read_csv("table_s3.tsv", sep="\t").rename(columns={'Unnamed: 0': 'type'})
e = x[x.type == 'edges']
r = x[x.type == 'runtime']

s =  {'orig': {'c': 'black', 'l': 'dash'},
      'I': {'c': 'lgreen', 'l': 'dot'},
      'IE': {'c': 'wine', 'l': None},
      'IS': {'c': 'lblue', 'l': 'dash'},
      'ISE': {'c': 'magenta', 'l': None},
}

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.set_xlabel("number of edges")
ax.set_ylabel("runtime (s)")

for f in x.file.unique():
    ax.plot(x[x.file == f].iloc[0, 2:], x[x.file == f].iloc[1, 2:], c='grey', zorder=0)

for name in ['orig', 'I', 'IE', 'IS', 'ISE']:
    ax.scatter(e[name], r[name], color=colors[s[name]['c']], label=name)

ax.legend()

plt.tight_layout()
plt.savefig("fig_s5.pdf")

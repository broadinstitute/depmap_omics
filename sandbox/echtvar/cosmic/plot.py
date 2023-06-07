import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
test=pd.read_table("cosmic_somatic_gnomad.stat", header=None)

fig, ax = plt.subplots()
fig.set_size_inches(9.2, 5.6)
sns.violinplot(x=2, y=3, data=test, ax=ax)
ax.set_ylim(0, 0.2)
ax.set_ylabel("Gnomad Allele frequency")
fig.autofmt_xdate(rotation=45)
fig.subplots_adjust(bottom=0.5) # or whatever
plt.savefig('cosmic_stat.pdf')

fig, ax = plt.subplots()
fig.set_size_inches(9.2, 5.6)
sns.violinplot(x=2, y=3, data=test, ax=ax)
##ax.set_ylim(0, 0.2)
ax.set_ylabel("Gnomad Allele frequency")
fig.autofmt_xdate(rotation=45)
fig.subplots_adjust(bottom=0.5) # or whatever
fig.savefig('cosmic_stat_all-scales.pdf')

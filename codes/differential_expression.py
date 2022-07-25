import pandas as pd
import numpy as np
import time as tm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

#Volcano and GO term enrichment plot (Fig. 6a and b
er = pd.read_csv("/Users/qingbowang/Desktop/covid_rna_de/n465_edger_lrt_12_vs_34_results.tsv", sep='\t') #from edgeR
plt.figure(figsize=(6,4))
plt.scatter(er.logFC, -1*np.log10(er.PValue), c=er.logCPM, cmap='viridis')
cbar = plt.colorbar(shrink=0.75)
cbar.set_label('log10(CPM+1)')
plt.xlabel("log2(Severe/Non-severe)")
plt.ylabel("-log10(raw p)")
plt.xlim([-2,6])
plt.savefig("/Users/qingbowang/Desktop/tmp/volcano_edger_real.png", bbox_inches='tight', dpi=400)
plt.clf()

#GO term
go = pd.read_csv("/Users/qingbowang/Desktop/covid_rna_de/n465_upreg_goenr.csv", sep=',')
go.source.value_counts()
go = go[go.source=="REAC"].reset_index() #reactome is the easiest
go["frac"] = go.intersection_size/go.query_size #% of the term that is covered
fig, ax = plt.subplots(1, 1, figsize=(9, 3))
scat = ax.scatter(go.negative_log10_of_adjusted_p_value, go.index[::-1], c=go.intersection_size,
                  cmap=cm.autumn, s=150, edgecolor="black")
ax.set_yticks(go.index)
ax.set_yticklabels(np.array(go.term_name[::-1].str.replace("Events associated with p","P"))) #making it shorter
ax.set_xlabel("-log10(adjusted p)")
cbar = fig.colorbar(scat, ax=ax)
cbar.set_label('Number of genes')
plt.savefig("/Users/qingbowang/Desktop/tmp/goterm_enr.png", bbox_inches='tight', dpi=400)
plt.clf()

#Diff. Intron Usage (Fig. S17):
lrt = pd.read_csv('/Users/qingbowang/Desktop/covid_rna_de/n465_edger_lrt_12_vs_34_introns.tsv', sep='\t')
plt.figure(figsize=(6,4))
plt.scatter(lrt.logFC, -1*np.log10(lrt.PValue), c=lrt.logCPM, cmap='viridis')
cbar = plt.colorbar(shrink=0.75)
cbar.set_label('log10(CPM+1)')
plt.xlabel("log2(Severe/Non-severe)")
plt.ylabel("-log10(raw p)")
plt.savefig("/Users/qingbowang/Desktop/tmp/volcano_edger_introns.png", bbox_inches='tight', dpi=400)
plt.clf()


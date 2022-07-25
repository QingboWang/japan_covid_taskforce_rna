import pandas as pd
import numpy as np
import time as tm

#reading BBJ and UKB PIP
M = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_fm/bbj_hemat_maxpip.tsv.gz", sep='\t')
Mb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_fm/ukb_hemat_maxpip.tsv.gz", sep='\t')
#annotate those
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    v = df.variant_id.str.split(":").str[4:].apply(lambda x: ":".join(x))
    #1. no flip
    df.index = v
    df.index.names = ["variant"]
    M.index = M.variant
    M.rename(columns = {"0":"bbj_pip"}, inplace=True)
    df = df.join(M.bbj_pip, how="left")
    Mb.index = Mb.variant.str.replace("chr","")
    Mb.index.names = ["variant"]
    Mb.rename(columns = {"pip":"ukb_pip"}, inplace=True)
    df = df.join(Mb.ukb_pip, how="left")
    print (sum(~df.bbj_pip.isna()), df.shape[0])
    print (sum(~df.ukb_pip.isna()), df.shape[0])
    cols = ["variant_id","gene_id","pip_fm","pip_susie","bbj_pip","ukb_pip",'pp_susie_gtex','pp_fm_gtex']
    df[cols].to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/chunk{0}_bbj_ukb.txt.gz".format(chk), sep='\t')
    df[(df.bbj_pip>0.1) | (df.ukb_pip>0.1)].to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/chunk{0}_bbj_ukb_01.txt.gz".format(chk), sep='\t')
    print ("done {0}, {1}".format(chk, tm.ctime()))


#and then get the enrichment and plot
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/chunk{0}_bbj_ukb.txt.gz".format(chk), sep='\t')

    #remove the ones gtex var is missing, for fair comparison
    df = df[(~df.pp_susie_gtex.isna()) & (~df.pp_fm_gtex.isna())]
    #probably taking the max per variant is the way to go
    df["eqtl_pip"] = np.minimum(df.pip_fm.fillna(-1), df.pip_susie.fillna(-1))
    df["min_pip_gtex"] = np.minimum(df.pp_susie_gtex.replace(0, -1),
                                    df.pp_fm_gtex.replace(0, -1))  # -1: dummy for missing ones
    df["comp_pip01"] = "Neither"
    df.loc[df.ukb_pip>0.1, "comp_pip01"] = "UKB_only"
    df.loc[df.bbj_pip>0.1, "comp_pip01"] = "BBJ_only"
    df.loc[((df.bbj_pip>0.1)&(df.ukb_pip>0.1)), "comp_pip01"] = "Both"
    df["eqtl_pip01"] = "Neither"
    df.loc[df.min_pip_gtex > 0.1, "eqtl_pip01"] = "gtex_only"
    df.loc[df.eqtl_pip > 0.1, "eqtl_pip01"] = "jtex_only"
    df.loc[((df.min_pip_gtex > 0.1) & (df.eqtl_pip > 0.1)), "eqtl_pip01"] = "Both"
    tb = df.groupby(["eqtl_pip01","comp_pip01"]).size().unstack(level=-1)

    #also do the per variant
    df.fillna(-1, inplace=True)
    M = df.groupby("variant_id")[["min_pip_gtex", "eqtl_pip","ukb_pip","bbj_pip"]].max()
    M["comp_pip01"] = "Neither"
    M.loc[M.ukb_pip>0.1, "comp_pip01"] = "UKB_only"
    M.loc[M.bbj_pip>0.1, "comp_pip01"] = "BBJ_only"
    M.loc[((M.bbj_pip>0.1)&(M.ukb_pip>0.1)), "comp_pip01"] = "Both"
    M["eqtl_pip01"] = "Neither"
    M.loc[M.min_pip_gtex > 0.1, "eqtl_pip01"] = "gtex_only"
    M.loc[M.eqtl_pip > 0.1, "eqtl_pip01"] = "jtex_only"
    M.loc[((M.min_pip_gtex > 0.1) & (M.eqtl_pip > 0.1)), "eqtl_pip01"] = "Both"
    tb2 = M.groupby(["eqtl_pip01", "comp_pip01"]).size().unstack(level=-1)
    if chk==0:
        tbsum = tb
        tb2sum = tb2
    else:
        tbsum = tbsum.add(tb, fill_value=0)
        tb2sum = tb2sum.add(tb2, fill_value=0)
    print ("done {0}, {1}".format(chk, tm.ctime()))
tbsum = tbsum.astype(int)
tb2sum = tb2sum.astype(int)
tbsum.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/coloc_per_pip_table_jg_intersect.tsv", sep="\t")
tb2sum.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/coloc_per_pip_table_jg_intersect_perv.tsv", sep="\t")
#perv:
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/coloc_per_pip_table_jg_intersect_perv.tsv",
                 sep="\t", index_col=0)
tb.columns = tb.columns.str.replace("_only", "")
y = (tb.T/tb.sum(axis=1)).T
yerr = (np.sqrt((y*(1-y)).T/tb.sum(axis=1))).T
#remove the 2 rows
y = y.iloc[[2,3],[0,3]].T
yerr = yerr.iloc[[2,3],[0,3]].T
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
x = np.arange(2)
c4 = {}
c4["UKB"] = "tab:blue"
c4["BBJ"] = "tab:red"
fig, ax = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(3.5,3))
cnt = 0
for i in range(2):#order manual
    label = y.index[i]
    ax[cnt].errorbar(x, y.iloc[i,:]*100, yerr.iloc[i,:]*100, color=c4[label], fmt="o")
    ax[cnt].set_ylabel(label)
    cnt += 1
ax[1].set_xlabel("eQTL PIP>0.1 only in:")
ax[1].set_xticks([0,1])
ax[1].set_xticklabels(["GTEx\n(n=9,715)","JCTF\n(n=10,776)"])
fig.text(-0.2, 0.5, '%Hematopoietic trait\n     PIP>0.1 only in:', va='center', rotation='vertical')
plt.subplots_adjust(hspace = .2)
plt.savefig('/Users/qingbowang/Desktop/plots/coloc_ukb_bbj_table.png', bbox_inches='tight', dpi=500)
plt.clf()

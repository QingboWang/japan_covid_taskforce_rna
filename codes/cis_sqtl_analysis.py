import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
c6 = ["tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]


#collect the fine-mapping results:
sus = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_n465_imputed_introns_0001.tsv", sep='\t', index_col=0)
fm = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_n465_imputed_introns_0001.txt", sep=' ', index_col=0)
#join the two to get the min -> keep only the ones with PIP>0.001 for both
print (sus.shape)
sus = sus.join(fm, how='inner')
print (sus.shape)
sus.columns = ["pip_susie", "pip_fm"]
sus.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465_imputed_introns_0001.tsv", sep='\t')


#1. get the value count; (later add the background value)
sus["min_pip"] = np.minimum(sus.pip_susie, sus.pip_fm)
df = sus
def bin_pip(pip_vector):
    out = pip_vector.copy(deep=True)
    out[pip_vector >= 0] = "(0,0.001]"
    out[pip_vector >= 0.001] = "(0.001,0.01]"
    out[pip_vector > 0.01] = "(0.01,0.1]"
    out[pip_vector > 0.1] = "(0.1,0.5]"
    out[pip_vector > 0.5] = "(0.5,0.9]"
    out[pip_vector > 0.9] = "(0.9,1]"
    return (out)
df["min_pip_bin"] = bin_pip(df.min_pip)
pervi = df.min_pip_bin.value_counts()
#perv and peri;
df["v"] = df.index.str.split("_").str[0]
df["i"] = df.index.str.split("_").str[1] + "_" + df.index.str.split("_").str[2]
perv_pip = df.groupby("v").min_pip.agg(max)
peri_pip = df.groupby("i").min_pip.agg(max)
perv = bin_pip(perv_pip).value_counts()
peri = bin_pip(peri_pip).value_counts()
tb = pd.concat([pervi, perv, peri], axis=1).fillna(0).astype(int)
tb.columns = ["vi", "v", "i"]
n_intron_all = 19499 #from wc -l
n_vi_all = 106020550 #from wc -l fm_n465_imputed_introns.txt
n_v_all  = 5872211 #from awk -> uniq -> wc -l -> -1 for the header
tb.iloc[-1,:] = np.array([n_vi_all, n_v_all, n_intron_all]) - tb.sum(axis=0) #subtracting those already counted
tb = tb.iloc[[5,0,1,2,3,4]]
tb.to_csv('/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_table.tsv',sep='\t')
#awk -F':' '{print $1"_"$2"_"$3"_"$4}' fm_n465_imputed_introns.txt > fm_n465_imputed_introns_variants.txt
#sort fm_n465_imputed_introns_variants.txt | uniq > fm_n465_imputed_introns_variants_unique.txt


#TSS distance:
df["g"] = df.i.str.split(":").str[-1]
dfg = df.sort_values(by="min_pip", ascending=False).drop_duplicates(subset=["g", "v"], keep="first")
dfg.set_index(["g", "v"], inplace=True, drop=True)
eq0 = []
for chk in range(26):
    eq = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465.imputed.pairs.chunk{0}.txt".format(chk), sep='\t', index_col=[0,1])
    eq = eq.loc[np.intersect1d(eq.index, dfg.index),"tss_distance"]
    eq0.append(eq)
    print ("done chk{0}, {1}".format(chk, tm.ctime()))
eq0 = pd.concat(eq0)
eq0.index.names = ["g","v"]
dfg = dfg.join(eq0, how='left')
dfg.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_vg_tssdist.tsv", sep="\t")


#CLPP with eQTL PIP (Fig. 3d)
eq0 = []
for chk in range(26):
    eq = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),
                     sep='\t', index_col=[1,0])
    eq = eq.loc[np.intersect1d(eq.index, dfg.index),["pip_susie", "pip_fm"]]
    eq0.append(eq)
    print ("done chk{0}, {1} rows, {2}".format(chk, len(eq), tm.ctime()))
eq0 = pd.concat(eq0)
eq0.index.names = ["g","v"]
eq0.columns = ["pip_susie_e", "pip_fm_e"]
dfg = dfg.join(eq0, how='left')
dfg.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_vg_eqtlpip.tsv", sep="\t")
#the distribution;
#remove those that are in the lowest bin, to begin with
dfg=dfg[dfg.min_pip_bin!="(0,0.001]"]
dfg["pip_min_e"] = np.minimum(dfg.pip_susie_e, dfg.pip_fm_e)
dfg["min_pip_bin_e"] = bin_pip(dfg.pip_min_e)
dfg["clpp"] = dfg.min_pip * dfg.pip_min_e
vg = bin_pip(dfg.clpp).value_counts()
perv_clpp = dfg.groupby("v").clpp.agg(max)
perg_clpp = dfg.groupby("g").clpp.agg(max)
perv = bin_pip(perv_clpp).value_counts()
perg = bin_pip(perg_clpp).value_counts()
tb = pd.concat([vg, perv, perg], axis=1).fillna(0).astype(int)
tb.columns = ["vg", "v", "g"]
tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_vg_eqtlpip_stats.tsv", sep="\t")
#plot this
tbnorm = tb/tb.sum(axis=0)
tbnorm = tbnorm.iloc[[0,1,2,3,5,4],] #manualに順番整える.
c6 = ["tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]
y = [0]
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(12, 6))
for i in range(tbnorm.shape[0]):
    ax[0].barh(y, tbnorm.iloc[i,0], left=tbnorm.iloc[:i,0].sum(), color=c6[i])
for i in range(tbnorm.shape[0]):
    ax[1].barh(y, tbnorm.iloc[i,1], left=tbnorm.iloc[:i,1].sum(), color=c6[i], label=tbnorm.index[i])
for i in range(tbnorm.shape[0]):
    ax[2].barh(y, tbnorm.iloc[i,2], left=tbnorm.iloc[:i,2].sum(), color=c6[i], label=tbnorm.index[i])
ax[2].set_xlim([0,1])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', fontsize=18, title="(max) CLPP bin")
ax[2].set_xlabel("Fraction", fontsize=18)
ax[0].set_title("n= 558,597 variant-genes", fontsize=18) #manual count sum(vg)
ax[1].set_title("n= 469,452 variants", fontsize=18) #manual count sum(v)
ax[2].set_title("n= 3,450 genes", fontsize=18) #manual count sum(g)
ax[0].set_yticks([])
ax[1].set_yticks([])
ax[2].set_yticks([])
plt.subplots_adjust(hspace = .5)
plt.savefig("/Users/qingbowang/Desktop/plots/jtex_sqtl_clpp_nums.png", bbox_inches='tight', dpi=500)
plt.clf()


#Splice donor/accepter, Splice AI score enrichment (Fig. 3b and c)
v = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_pip0001_vep.txt", sep='\t')
v.index = v.Location.str.split("-").str[0] + ":" + v.Allele + ":" + v.Gene
sp = v[v.Consequence.str.contains("splice")] #~10k
sp.Consequence.value_counts()
spdn_v = sp[sp.Consequence.str.contains("splice_donor")].index.unique()
spac_v = sp[sp.Consequence.str.contains("splice_acceptor")].index.unique()
sprg_v = sp[sp.Consequence.str.contains("splice_region")].index.unique()

dfg = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sqtl_vg_eqtlpip.tsv", sep="\t")
dfg.index = dfg.v.str.split(":").str[0].str.replace("chr","") + ":" + \
            dfg.v.str.split(":").str[1] + ":" + dfg.v.str.split(":").str[3] + ":" + dfg.g.str.split("\\.").str[0] #strange key of chr pos alt gene
dfg["acceptor"] = False
dfg["donor"] = False
dfg["region"] = False
dfg.loc[np.intersect1d(dfg.index, spac_v),"acceptor"] = True
dfg.loc[np.intersect1d(dfg.index, spdn_v),"donor"] = True
dfg.loc[np.intersect1d(dfg.index, sprg_v),"region"] = True
act = dfg.groupby("min_pip_bin").acceptor.agg(["count", "sum", "mean", "sem"]).fillna(0)
dnt = dfg.groupby("min_pip_bin").donor.agg(["count", "sum", "mean", "sem"]).fillna(0)
rgt = dfg.groupby("min_pip_bin").region.agg(["count", "sum", "mean", "sem"]).fillna(0)

#spliceAI score as well
sp = v[(v.SpliceAI_pred_DS_AG!="-") | (v.SpliceAI_pred_DS_AL!="-") | \
       (v.SpliceAI_pred_DS_DG!="-") | (v.SpliceAI_pred_DS_DL!="-")].iloc[:,-5:-1]
sp = sp.replace("-","0").astype(float)
sp = sp[sp.sum(axis=1)>0]
print (sp.shape[0], len(sp.index.unique())) #would like to reduce to the max score
sp = sp[~sp.index.duplicated(keep='first')]
#and merge to original
dfg = dfg.join(sp, how="left")
dfg.fillna(0, inplace=True)
splice_cols = dfg.columns[-4:]
spt = dfg.groupby("min_pip_bin")[splice_cols].agg(["mean", "sem"]).fillna(0)
dfg["spliceai_nonzero"] = dfg.loc[:,splice_cols].sum(axis=1)>0
dfg.groupby("min_pip_bin").spliceai_nonzero.agg(["mean", "sem"]).fillna(0)
dfg["spliceai_high"] = dfg.loc[:,splice_cols].max(axis=1)>0.5
dfg.groupby("min_pip_bin").spliceai_high.agg(["mean", "sem"]).fillna(0)
dfg["spliceai_high"] = dfg.loc[:,splice_cols].max(axis=1)>0.9
dfg.groupby("min_pip_bin").spliceai_high.agg(["mean", "sem"]).fillna(0)
def bin_pip_loose(pip_vector):
    out = pip_vector.copy(deep=True)
    out[pip_vector >= 0.001] = "(0.001,0.01]"
    out[pip_vector > 0.01] = "(0.01,0.1]"
    out[pip_vector > 0.1] = "(0.1,0.9]"
    out[pip_vector > 0.9] = "(0.9,1]"
    return (out)
dfg["min_pip_bin_loose"] = bin_pip_loose(dfg.min_pip)
act = dfg.groupby("min_pip_bin_loose").acceptor.agg(["count", "sum", "mean", "sem"]).fillna(0)
dnt = dfg.groupby("min_pip_bin_loose").donor.agg(["count", "sum", "mean", "sem"]).fillna(0)
act = act.iloc[1:,:]
dnt = dnt.iloc[1:,:]
#plot
plt.rcParams.update({'font.size': 8})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
c4 = ["tab:blue", "tab:green", "tab:orange", "tab:red"]
plt.figure(figsize=(3.6,1.5))
plt.errorbar(-0.1,dnt.iloc[0,:]["mean"], dnt.iloc[0,:]["sem"], fmt="o", color="black", markersize=3.5, label="Canonical splice donor") #for legend
plt.errorbar(0.1,act.iloc[0,:]["mean"], act.iloc[0,:]["sem"], fmt="D", color="black", markersize=3.5, label="Canonical splice acceptor")
for i in range(len(c4)):
    plt.errorbar(i-0.1, dnt.iloc[i,:]["mean"], dnt.iloc[i,:]["sem"], fmt="o", color=c4[i], markersize=4)
    plt.errorbar(i+0.1, act.iloc[i,:]["mean"], act.iloc[i,:]["sem"], fmt="D", color=c4[i],markersize=4)
plt.xlabel("sQTL PIP bin")
plt.ylabel("Fraction")
plt.legend(loc="upper left")
plt.xticks(np.arange(len(c4)),act.index, rotation=15)
plt.savefig("/Users/qingbowang/Desktop/plots/canonical_transcirpt_enr.png", bbox_inches='tight', dpi=500)
plt.clf()

#splice AI score
dfg["spliceai_max"] = dfg.loc[:,splice_cols].max(axis=1)
tb = []
for i in np.arange(101) / 100:
    j = (dfg[dfg.min_pip >= i].spliceai_max.replace(0, -1)//0.1).value_counts()
    #j = (np.floor(dfg[dfg.min_pip >= i].spliceai_max.replace(0,-1)*10)/10).value_counts() #exact one does not really exist
    tb.append(j)
tb = pd.DataFrame(tb)
tb.index = range(101)
tb = tb.fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/spliceai_scores.tsv", sep='\t')

tbnorm = tb.apply(lambda x: x / sum(x), axis=1)
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
xlb = ["0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]",
       "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]"]
from matplotlib import cm
vir = cm.viridis
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]}) #the total num is almost obvious
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-\ngene)", fontsize=13)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label="=0", color="tab:gray", edgecolor="#888888ff", linewidth=0.2)
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1),
                      label=xlb[i], color=vir(int(i/(tbnorm.shape[1]-1)*256)),edgecolor="#888888ff", linewidth=0.2)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="SpliceAI score bin")
ax[1].set_xlabel("sQTL significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/spliceai_scores.png', bbox_inches='tight', dpi=500)
plt.clf()

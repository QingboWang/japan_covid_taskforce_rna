import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from sklearn.linear_model import LinearRegression
from scipy import stats
from matplotlib import cm
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
c4 = ["tab:gray", "tab:blue", "tab:orange", "tab:red"]
c6 = ["tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]
c8 = ["tab:gray", "tab:brown", "tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]



#marginal beta agreement (Fig. 2f):
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465.imputed.pairs.gtex_beta.txt", sep='\t')
plt.scatter(df.slope, df.beta_gtex, alpha=0.1, color="tab:gray")
plt.scatter(df[df.pval_nominal<5*10**-8].slope, df[df.pval_nominal<5*10**-8].beta_gtex, alpha=0.1, color="tab:red")
plt.xlabel("beta (JCT)")
plt.ylabel("beta (GTEx)")
plt.savefig("/Users/qingbowang/Desktop/plots/scatter_beta_p.png", dpi=500)
plt.clf()
df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),
                           df.pip_fm.fillna(-1))  # -1 = those that are not fine-mapped since it is not egene
df["min_pip_bin"] = "Undefined"
df["pip_bin"] = "(0,0.001]"  # to begin with
df.loc[df.min_pip >= 0, "pip_bin"] = "(0,0.001]"
df.loc[df.min_pip > 0.001, "pip_bin"] = "(0.001,0.01]"
df.loc[df.min_pip > 0.01, "pip_bin"] = "(0.01,0.1]"
df.loc[df.min_pip > 0.1, "pip_bin"] = "(0.1,0.5]"
df.loc[df.min_pip > 0.5, "pip_bin"] = "(0.5,0.9]"
df.loc[df.min_pip > 0.9, "pip_bin"] = "(0.9,1]"
bins = ["(0,0.001]","(0.001,0.01]","(0.01,0.1]","(0.1,0.5]","(0.5,0.9]","(0.9,1]"]
plt.figure(figsize=(5,4.75))
for i in range(len(bins)):
    plt.scatter(df[df.pip_bin==bins[i]].slope, df[df.pip_bin==bins[i]].beta_gtex, alpha=1, color=c6[i], label=bins[i])
plt.axvline(x=0, linewidth=0.5, color="black", linestyle="--")
plt.axhline(y=0, linewidth=0.5, color="black", linestyle="--")
plt.plot([-3,3], [-3,3], linewidth=0.5, color="black", linestyle="--")
plt.legend(title="PIP bin", loc="lower right", fontsize=10)
plt.xlabel(r"Effect size $\beta$ (JCT)")
plt.ylabel(r"Effect size $\beta$ (GTEx)")
plt.xlim([-2.95,2.95])
plt.ylim([-2.95,2.95])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/scatter_beta_pip.png", dpi=500)

#pearson R:
v = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_beta_sign_tb.tsv",sep='\t', index_col=0)
err = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_beta_sign_err_tb.tsv",sep='\t', index_col=0)
v = v.loc[bins, bins]
err = err.loc[bins, bins]
pear = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_beta_sign_pearsonr.tsv",sep='\t', index_col=0)
pear_u = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_beta_sign_pearsonr_upper.tsv",sep='\t', index_col=0)
pear_l = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_beta_sign_pearsonr_lower.tsv",sep='\t', index_col=0)
pear = pear.loc[bins, bins]
pear_u = pear_u.loc[bins, bins]
pear_l = pear_l.loc[bins, bins]

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7,4))
ax[0].set_ylim([0.78,1.01])
ax[1].set_ylim([0.78,1.01])
for i in range(len(bins)):
    ax[0].errorbar(np.arange(len(bins))-0.25+i*0.1, v.loc[bins[i], :], err.loc[bins[i],:], color=c6[i], fmt="o", markersize=6)
    ax[1].errorbar(np.arange(len(bins))-0.25+i*0.1, pear.loc[bins[i], :],
                   [list(pear.loc[bins[i], :]-pear_l.loc[bins[i],:]),list(pear_u.loc[bins[i],:]-pear.loc[bins[i], :])], color=c6[i], fmt="o", markersize=6)
i = 5
ax[1].errorbar([-0.25+i*0.1], pear.loc[bins[i], bins[0]],
                   [[pear.loc[bins[i], bins[0]]-0.79], [pear_u.loc[bins[i],bins[0]]-pear.loc[bins[i], bins[0]]]],
               color=c6[i], fmt="o", markersize=6, uplims=True)
ax[0].axhline(y=1, linewidth=0.5, linestyle='--', zorder=2, color='black')
ax[1].axhline(y=1, linewidth=0.5, linestyle='--', zorder=2, color='black')
ax[1].set_xlabel("PIP bin (JCTF)")
ax[0].set_ylabel("Prob. same\neffect direction")
ax[1].set_ylabel("Effect sizes corr.\n(pearson r)")
ax[1].set_xticks(np.arange(len(bins)))
ax[1].set_xticklabels(bins, rotation=30)
for xtick, color in zip(ax[1].get_xticklabels(), c6):
    xtick.set_color(color)
plt.tight_layout()
plt.subplots_adjust(hspace = .05)
plt.savefig("/Users/qingbowang/Desktop/plots/eff_size_agreement.png", dpi=500)



#overview of JCTF vs GTEx (Fig. 2 b,d etc)
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_vs_gtex_pval_combined.tsv",sep='\t', index_col=0)
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c4[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c4[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="P-value in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (-log10p)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_p_vs_gtex_p.png', bbox_inches='tight', dpi=500)
plt.clf()
#In terms of PIP as well:
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_gtex_pip_combined.tsv", sep='\t', index_col=0)
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c8[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c8[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (-log10p)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_p_vs_gtex_pip_full.png', bbox_inches='tight', dpi=500)
plt.clf()

#removing undefined and missing
tbrm = tb.iloc[:,2:]
tbrmnorm = tbrm.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tbrm.index, tbrm.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tbrm.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbrmnorm.index, tbrmnorm.iloc[:,0], label=tbrmnorm.columns[0], color=c6[0])
for i in range(1,tbrmnorm.shape[1]):
    ax[1].bar(tbrmnorm.index, tbrmnorm.iloc[:,i], bottom=tbrmnorm.iloc[:, :i].sum(axis=1), label=tbrmnorm.columns[i], color=c6[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (-log10p)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_p_vs_gtex_pip_main.png', bbox_inches='tight', dpi=500)
plt.clf()

#JCTF PIP vs gtex P
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pval_combined.tsv",sep='\t', index_col=0)
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c4[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c4[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="P-value in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_vs_gtex_p.png', bbox_inches='tight', dpi=500)
plt.clf()

#JCTF PIP vs GTEx PIP, full
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_combined.tsv",sep='\t', index_col=0)
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c8[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c8[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_vs_gtex_pip_full.png', bbox_inches='tight', dpi=500)
plt.clf()

#removing undefined and missing ones
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_combined.tsv",sep='\t', index_col=0)
tb = tb.iloc[:,2:]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c6[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c6[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in GTEx:", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_vs_gtex_pip_main.png', bbox_inches='tight', dpi=500)
plt.clf()

#The expression modifier score (EMS) (Fig. 2g)
tbp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_ems_combined.tsv",sep='\t', index_col=0)
tbpip = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_ems_combined.tsv", sep='\t', index_col=0)
tbp.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
tbpip.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
from matplotlib import cm
vir = cm.viridis
tb = tbpip
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=vir(0))
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=vir(int(i/(tbnorm.shape[1]-1)*256)))
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Expression Modifier Score (EMS):", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_vs_ems_main.png', bbox_inches='tight', dpi=500)
plt.clf()

#FM vs SuSiE (Fig. S5)

tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/susie_vs_pip.tsv", sep='\t', index_col=0)
log_norm = LogNorm()
plt.figure(figsize=(8,8))
sns.heatmap(tb.iloc[::-1,:]+1, annot=tb.iloc[::-1,:], fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=20, fontsize=15)
plt.xticks(rotation=20, fontsize=15)
plt.xlabel("FINEMAP PIP", fontsize=16)
plt.ylabel("SuSiE PIP", fontsize=16)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_fm_vs_susie.png', bbox_inches='tight', dpi=500)
plt.clf()
#FM vs SuSiE in terms of EMS distribution
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/susie_vs_pip_ems.tsv", sep='\t', index_col=[0,1]) #unstack again when plotting
tb = tb.astype(int)
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
from matplotlib.patches import Rectangle
cols = ["(0,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"]
vir = cm.viridis
fig, ax = plt.subplots(len(cols), len(cols), sharex=True, sharey=True, figsize=(8, 8))
for i in range(len(cols)): #jtex
    for j in range(len(cols)): #gtex
        v = tb.loc[(tb.index.get_level_values(level=0)==cols[::-1][i]) & (tb.index.get_level_values(level=1)==cols[j]),:]
        v = np.array(v.iloc[0, :]) #to make it an array
        n = sum(v) #to keep it later
        v = np.cumsum(v/sum(v))
        cnt = len(v)-1
        for c in range(len(v))[::-1]:
            if v[c]==0: l = 0
            else: l = np.sqrt(v[c]) #so that the n is the fraction
            pch = Rectangle((-l/2,-l/2), l, l, facecolor=vir(int(cnt/(len(v)-1)*256)))
            ax[i,j].add_patch(pch)
            cnt -= 1
        ax[i,j].set_xlim([-0.5,0.5])
        ax[i,j].set_ylim([-0.5, 0.5])
        if i==(len(cols)-1):
            ax[i,j].set_xlabel(cols[j])
        if j==0:
            ax[i,j].set_ylabel(cols[::-1][i])
        #remove ticks
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        #add title
        ax[i,j].set_title("n={0}".format(n), fontsize=13, pad=0.3)
fig.text(0.5, 0.04, 'FINEMAP PIP bin', ha='center')
fig.text(0.04, 0.5, 'SuSiE PIP bin', va='center', rotation='vertical')
plt.savefig("/Users/qingbowang/Desktop/plots/fm_v_susie_ems_heatmap.png", bbox_inches='tight', dpi=500)
plt.clf()



#ROC for predicting GTEx putative causal variants (Fig. 2e)
from sklearn import metrics
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

st = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/predicting_gtex_putative_causal.tsv", sep='\t')
y = st.gtex_causal
x1 = st.pval
x2 = st.pip
fpr1, tpr1, threshold1 = metrics.roc_curve(y, x1)
roc_auc1 = metrics.auc(fpr1, tpr1)
fpr2, tpr2, threshold2 = metrics.roc_curve(y, x2)
roc_auc2 = metrics.auc(fpr2, tpr2)
plt.figure(figsize=(5.5,4.5))
plt.plot(fpr1, tpr1, c='tab:green', label = 'P-value, AUROC = %0.3f' % roc_auc1)
plt.plot(fpr2, tpr2, c='tab:purple', label = 'PIP, AUROC = %0.3f' % roc_auc2)
plt.title("Task: Prioritization of PIP>0.9 in GTEx")
plt.ylabel("True positive rate (TPR)", fontsize=14)
plt.xlabel("False positive rate (FPR)", fontsize=14)
plt.legend()
plt.savefig('/Users/qingbowang/Desktop/plots/gtexcausal_prediction_roc.png', bbox_inches='tight', dpi=500)
plt.clf()

precision1, recall1, thresholds1 = metrics.precision_recall_curve(y, x1)
precision2, recall2, thresholds2 = metrics.precision_recall_curve(y, x2)
auc1 = metrics.average_precision_score(y, x1)
auc2 = metrics.average_precision_score(y, x2)
plt.figure(figsize=(5.5,4.5))
plt.plot(recall1, precision1, c='tab:green', label = 'P-value, AUPRC = %0.3f' % auc1)
plt.plot(recall2, precision2, c='tab:purple', label = 'PIP, AUPRC = %0.3f' % auc2)
plt.title("Task: Prioritization of PIP>0.9 in GTEx", fontsize=14)
plt.xlabel("Recall", fontsize=14)
plt.ylabel("Precision", fontsize=14)
plt.legend()
plt.savefig('/Users/qingbowang/Desktop/plots/gtexcausal_prediction_prc.png', bbox_inches='tight', dpi=500)
plt.clf()


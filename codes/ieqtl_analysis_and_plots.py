import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans


#ieQTL calls (Fig. 7):
phenotype_bed_file = '/home/qwang/n465_imputed_eqtl_call_files/n465.expression.bed.gz'
covariates_file = "/home/qwang/n465_imputed_eqtl_call_files/n465.combined_covariates_idadded.txt"
iterm_file = "/home/qwang/n465_imputed_eqtl_call_files/interaction_f.txt" #header no need
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates
iterm_df = pd.read_csv(iterm_file, sep='\t', index_col=0, header=None, squeeze=True) #this needs to be series in advance
#perform ieQTL call:
import time as tm
for chr in range(1,23):
    # read in the genotype matrix, just in pandas
    print ("starting chr{0}, {1}".format(chr, tm.ctime()))
    prefix = "ie_chr{0}".format(chr)
    vcfpath = "/home/qwang/covid_imputed_genotypes/ct_imputed_hg38_sorted_chr{0}.vcf.gz".format(chr)
    df = pd.read_csv(vcfpath, sep="\t", skiprows=3384)
    print("read the df, chr{0}, {1}".format(chr, tm.ctime()))
    df.index = df.ID
    genotype_df = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[-1]))
    variant_df = df[["#CHROM","POS"]]
    variant_df.columns = ["chrom", "pos"]
    #also phenotype df needs to be subsetted...
    phenotype_df_sub = phenotype_df.loc[phenotype_pos_df.chr=="chr{0}".format(chr),:]
    phenotype_pos_df_sub = phenotype_pos_df.loc[phenotype_pos_df.chr=="chr{0}".format(chr),:]
    ie_df = cis.map_nominal(genotype_df, variant_df, phenotype_df_sub, phenotype_pos_df_sub, prefix,
                    covariates_df=covariates_df, maf_threshold=0.01,
                    interaction_s=iterm_df, maf_threshold_interaction=0.05,
                    run_eigenmt=True, output_dir='.', write_top=True, write_stats=True)
    print("done chr{0}, {1}".format(chr, tm.ctime()))
#collect the results into single file:
df0 = []
for chr in range(1,23):
    df = pd.read_csv("/home/qwang/ie_chr{0}.cis_qtl_top_assoc.txt.gz".format(chr),sep='\t')
    df0.append(df)
df0 = pd.concat(df0)
df0.to_csv("/home/qwang/ie_cis_qtl_top_assoc.txt",sep='\t', index=False)





#cibersort cell type fraction eQTLs (Fig. 8):
cs = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_n500_CIBERSORT.txt", sep='\t', index_col=0)
pval = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/cibersort_ie_table.tsv", sep='\t', index_col=0)
dir_concord = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/cibersort_ie_dir_table.tsv", sep='\t', index_col=0)
pval.columns = pval.columns.str.replace("pval_gi_","")
dir_concord.columns = dir_concord.columns.str.replace("b_gi_","")
#conversion to gene name;
gns = pd.read_csv("~/Desktop/resources/gene_names.txt", sep='\t', index_col=0)
gene_names = gns.loc[pval.index.str.split("\\.").str[0],"Gene name"]
#bonf = 617
bonf = (13*22) #n(genes)*n(celltypes) are tested essentially
pval.index = gene_names
dir_concord.index = gene_names
colors = pval.copy(deep=True)
colors[(dir_concord)] = 0.5
colors[(~dir_concord)] = 0.5#baseline: white is fine
colors[(pval<0.05)&(dir_concord)] = 0.55
colors[(pval<0.05/bonf)&(dir_concord)] = 0.75
colors[(pval<0.05)&(~dir_concord)] = 0.45
colors[(pval<0.05/bonf)&(~dir_concord)] = 0.25
#put interesting things to left top
roworder = (((pval<0.05/bonf)&(dir_concord)).sum(axis=1) - ((pval<0.05/bonf)&(~dir_concord)).sum(axis=1) + ((pval<0.05/bonf)&(dir_concord)).sum(axis=1)*0.1).sort_values(ascending=False).index
colorder = (((pval<0.05/bonf)&(dir_concord)).sum() - ((pval<0.05/bonf)&(~dir_concord)).sum() + ((pval<0.05/bonf)&(~dir_concord)).sum()*0.1).sort_values(ascending=False).index
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib.colors import LogNorm
import seaborn as sns
cmap=plt.cm.seismic
plt.figure(figsize=(8,8))
sns.heatmap(colors.loc[roworder, colorder], annot=-np.log10(pval.loc[roworder, colorder]), fmt=".1f", square=True, linewidths=.5, linecolor="black",
            cmap="seismic", cbar=False, annot_kws={"size": 10}, vmin=0, vmax=1) #legend; later
plt.yticks(rotation=0, fontsize=14, style='italic')
plt.xticks(rotation=90, fontsize=14)
plt.xlabel("Cell type", fontsize=14)
plt.ylabel("Gene name", fontsize=14)
plt.title("-log$_{10}$(p) for cell type composition-ieQTL effect")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_revision_fig/cibersort_heatmap.png', bbox_inches='tight', dpi=500)
plt.clf()


#comparison with ImmunextUT (Fig. 7b etc)
fs = glob.glob("/Users/qingbowang/Desktop/resources/immunexut/*.txt")
ns = pd.Series(fs).str.split("immunexut/").str[-1].str.replace("_gwsig.txt","")
import numpy as np
#per cell type
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/immunexut_signif_n_per_celltype_per_p.tsv", sep='\t', index_col=0)
tb0 = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/n_per_jtex_p_base.tsv", sep='\t', index_col=0)
tbnorm = (tb.T/tb0.iloc[:,0]).T
fig, ax = plt.subplots(7, 4, sharex=True, sharey=True, figsize=(8, 11.5)) #the total num is almost obvious
for i in range(7):
    for j in range(4):
        if j==0:
            ax[i][j].set_ylim([-0.01, 1.01])
        if i==0:
            ax[i][j].set_xlim([-1,101])
        ax[i][j].bar(tbnorm.index, tbnorm.iloc[:, i*4 + j], color="tab:red")
        ax[i][j].bar(tbnorm.index, 1-tbnorm.iloc[:, i*4 + j], bottom=tbnorm.iloc[:,i*4 + j], color="tab:gray")
        ax[i][j].set_title(ns[i*4 + j], pad=5, fontsize=12)
#for legend inner
ax[0][0].bar(tbnorm.index, tbnorm.iloc[:, i*4 + j], color="tab:red", label="True")
ax[0][0].bar(tbnorm.index, 1-tbnorm.iloc[:, i*4 + j], bottom=tbnorm.iloc[:,i*4 + j], color="tab:gray", label="False")
ax[0][0].legend(loc='upper left', title="p<5e-8\nin ImmuNextUT", fontsize=6, title_fontsize=6)
ax[6][1].set_xlabel("                        Significance threshold (-log10p)", fontsize=15)
ax[3][0].set_ylabel("Fraction", fontsize=15)
plt.subplots_adjust(hspace = .35, wspace=.15)
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pval_vs_immunexut_per_tissue.png', bbox_inches='tight', dpi=500)
plt.clf()
#same thing with PIP
tb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/immunexut_signif_n_per_celltype_per_pip.tsv", sep='\t', index_col=0)
tb0 = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/n_per_jtex_pip_base.tsv", sep='\t', index_col=0)
tbnorm = (tb.T/tb0.iloc[:,0]).T
fig, ax = plt.subplots(7, 4, sharex=True, sharey=True, figsize=(8, 11.5)) #the total num is almost obvious
for i in range(7):
    for j in range(4):
        if j==0:
            ax[i][j].set_ylim([-0.01, 1.01])
        if i==0:
            ax[i][j].set_xlim([-1,101])
            ax[i][j].set_xticklabels(list(np.array(ax[i][j].get_xticks().tolist()) / 100))
        ax[i][j].bar(tbnorm.index, tbnorm.iloc[:, i*4 + j], color="tab:red")
        ax[i][j].bar(tbnorm.index, 1-tbnorm.iloc[:, i*4 + j], bottom=tbnorm.iloc[:,i*4 + j], color="tab:gray")
        ax[i][j].set_title(ns[i*4 + j], pad=5, fontsize=12)
#for legend inner
ax[0][0].bar(tbnorm.index, tbnorm.iloc[:, i*4 + j], color="tab:red", label="True")
ax[0][0].bar(tbnorm.index, 1-tbnorm.iloc[:, i*4 + j], bottom=tbnorm.iloc[:,i*4 + j], color="tab:gray", label="False")
ax[0][0].legend(loc='upper left', title="p<5e-8\nin ImmuNextUT", fontsize=6, title_fontsize=6)
ax[6][1].set_xlabel("                      Significance threshold (PIP)", fontsize=15)
ax[3][0].set_ylabel("Fraction", fontsize=15)
plt.subplots_adjust(hspace = .35, wspace=.15)
plt.savefig('/Users/qingbowang/Desktop/plots/jtex_pip_vs_immunexut_per_tissue.png', bbox_inches='tight', dpi=500)
plt.clf()


import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


#Collecting the result of trans-eQTL mapping:
trans_df = []
for chr in range(1,23):
    df = pd.read_csv("/home/qwang/n465_imputed_eqtl_call_files/chr{0}_trans_eqtls.tsv".format(chr), sep="\t")
    trans_df.append(df)
trans_df = pd.concat(trans_df).iloc[:,1:]
trans_df.to_csv("/home/qwang/n465_imputed_eqtl_call_files/trans_eqtls_all.tsv", sep="\t")

#co-loc with cis-eQTLs:
#get the list of putative cis-eQTLs at PIP>0.9
dfs = []
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),sep='\t')
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),df.pip_fm.fillna(-1))
    df = df[df.min_pip>0.1] #loose filtering
    dfs.append(df)
    print ("done chk {0}, {1}".format(chk, tm.ctime()))
pd.concat(dfs, axis=0).to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/cis_pip01_pairs.tsv", sep="\t", index=False)

dfc = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/cis_pip01_pairs.tsv", sep="\t")
dfc["hg19"] = dfc.variant_id.str.split(":").str[-4] + ":" + \
              dfc.variant_id.str.split(":").str[-3] + ":" + \
              dfc.variant_id.str.split(":").str[-2] + ":" + \
              dfc.variant_id.str.split(":").str[-1]
co = np.intersect1d(trans_df.variant_id, dfc[dfc.min_pip>0.1].hg19) #491
trans_df.index = trans_df.variant_id
trans_df.index.names = ["hg19"]
cotran = trans_df.loc[co,:]
dfc.index = dfc.hg19
cocis = dfc.loc[co,:]
genes = cotran.groupby("hg19").phenotype_id.apply(list)
ps = cotran.groupby("hg19").pval.apply(list)
betas = cotran.groupby("hg19").b.apply(list)
ses = cotran.groupby("hg19").b_se.apply(list)

cocis = cocis.join(genes, how='outer') #This how shouldn't matter
cocis = cocis.join(ps, how='outer')
cocis = cocis.join(betas, how='outer')
cocis = cocis.join(ses, how='outer')
cosis = cocis.rename(columns={"phenotype_id":"trans_genes", "pval":"trans_pval",
                      "b":"trans_beta", "b_se":"trans_beta_se"})
#add gene names
gn = pd.read_csv("~/Desktop/resources/gene_names.txt", sep="\t")
#first filter to those that exists in our cosis
g = np.union1d(pd.Series(cosis.gene_id.unique()).apply(lambda x: x.split(".")[0]),
               pd.Series(cotran.phenotype_id.unique()).apply(lambda x: x.split(".")[0])) #genes of interest
gn.index = gn['Gene stable ID']
gnuse = gn.loc[np.intersect1d(gn.index, g),:]
def find_gn(x, df=gnuse):
    try: return (df.loc[x,"Gene name"])
    except: return ("")
def find_gn_in_list(l, df=gnuse):
    out = []
    for g in l:
        g = g.split(".")[0]
        try: out.append(df.loc[g,"Gene name"])
        except: out.append("")
    return (out)
cosis["gene_name"] = cosis.gene_id.str.split("\\.").str[0].apply(lambda x: find_gn(x))
cosis["trans_gene_names"] = cosis.trans_genes.apply(lambda x: find_gn_in_list(x))
cosis.to_csv('~/Desktop/covid_eqtl_imputed_fm/trans_colocalizing_putative_causal_vars_imputed_pip01_p08.tsv', sep='\t') #loose threshold


#next, check per PIP trans-eQTL enrichment (Fig. 5b, c)
up_base = pd.read_csv('/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/transe_enr_pip_background.tsv', sep="\t", index_col=0)
up = pd.read_csv('/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/transe_enr_bonf_pip.tsv', sep="\t", index_col=0)
#and the baseline (random) is;
#for trans eQTL in eQTLgen;
rd_trans = u.frac[0] #frac for random_J_var
up["frac_enr"] = up.frac / rd_trans
#for random var
rd = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/trans_eff_hg19unique_vars.tsv", sep="\t", index_col=0)
bk = pd.read_csv("~/Downloads/eqtlgen_testedsnps.csv", sep=',')
bk["hg19"] = bk.Chr.astype(str) + ":" +  \
             bk.Position.astype(str) + ":" +  \
             bk.Alleles.str.split("/").str[0] + ":" +  \
             bk.Alleles.str.split("/").str[1]
bk = bk.drop_duplicates(subset=["hg19"])
bk.index = bk.hg19
rd = rd.join(bk.Trait, how="left")
rd_intestset = (~rd.Trait.isna()).mean()
(~rd.Trait.isna()).value_counts()
up_base["frac_enr"] = up_base.frac / rd_intestset
#add the error bar and plot - this one is right.
up["n"] = up.iloc[:,:2].sum(axis=1)
up_base["n"] = up_base.iloc[:,:2].sum(axis=1)
up["err"] = np.sqrt(up.frac*(1-up.frac)/up.n)
up_base["err"] = np.sqrt(up_base.frac*(1-up_base.frac)/up_base.n)
up["err_enr"] = up.err / rd_trans
up_base["err_enr"] = up_base.err / rd_intestset
up.reset_index(inplace=True)
up_base.reset_index(inplace=True)
plt.figure(figsize=(4,3))
c4 = ["mediumseagreen", "tab:olive", "tab:orange", "tab:red"]
i = 0
plt.errorbar([i-0.1], up.frac_enr[i], up.err_enr[i], fmt="o", color="black", label="trans-eQTL\n(FDR<0.05)",mec='white', mew=0.25)
plt.errorbar([i+0.1], up_base.frac_enr[i], up_base.err_enr[i], fmt="D", color="black", label="Tested SNP\n(background)",mec='white', mew=0.25)
plt.legend(title="In eQTLgen:", loc="upper left")
for i in range(4):
    plt.errorbar([i-0.1], up.frac_enr[i], up.err_enr[i], fmt="o", color=c4[i],mec='black', mew=0.25)
    plt.errorbar([i+0.1], up_base.frac_enr[i], up_base.err_enr[i], fmt="D", color=c4[i],mec='black', mew=0.25)
plt.xlabel("cis-eQTL PIP threshold")
plt.ylabel("Enrichment")
plt.xticks([0,1,2,3],["0.01<PIP\n(n={0})".format(up.n[0]),
                      "0.1<PIP\n(n={0})".format(up.n[1]),
                      "0.5<PIP\n(n={0})".format(up.n[2]),
                      "0.9<PIP\n(n={0})".format(up.n[3])], fontsize=10, rotation=30)
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/trans_enr_pip.png", dpi=500)
plt.clf()

#sign match with eQTLgen (Fig. 5a)
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
gen = pd.read_csv("~/Downloads/eqtlgen_trans_eqtls.csv", sep=',') #Not only unique, for now
gen["hg19"] = gen["SNP-Chr"].astype(str) + ":" + gen["SNP-Pos"].astype(str) + ":" + \
              gen["Alleles"].str.split("/").str[0] + ":" + \
              gen["Alleles"].str.split("/").str[1]
trans_df = pd.read_csv("~/Downloads/trans_eqtls_all.tsv", sep="\t", index_col=0)
trans_df.index = [trans_df.variant_id, trans_df.phenotype_id.str.split("\\.").str[0]]
gen.index = [gen.hg19,gen.Gene]
gen.index.names = trans_df.index.names
intersection = trans_df.join(gen, how="inner", rsuffix="eqtlgen")
intersection["z"] = intersection.Z
intersection.loc[intersection.AlleleAssessed==intersection.Alleles.str.split("/").str[0],"z"] = intersection.loc[intersection.AlleleAssessed==intersection.Alleles.str.split("/").str[0],"Z"]*-1
#flipping the minor major difference
intersection["sign_match"] = np.nan
intersection.loc[(intersection.b * intersection.z)>0,"sign_match"] = True
intersection.loc[(intersection.b * intersection.z)<0,"sign_match"] = False
intersection.sign_match.value_counts()
#Perfect 35/35 match in terms of the exact variant-gene
#maybe take the corr then?
intersection["jct_z"] = intersection.b / intersection.b_se
from scipy import stats
import matplotlib.cm as cm
linreg = stats.linregress(intersection.jct_z, intersection.z)
plt.figure(figsize=(4,3))
sc = plt.scatter(intersection.jct_z, intersection.z, edgecolors="black", c=np.log10(intersection.pval)*-1, cmap=cm.viridis)
plt.axvline(x=0, linewidth=0.5, linestyle="--", color="black")
plt.axhline(y=0, linewidth=0.5, linestyle="--", color="black")
plt.text(1, 60, "r={0}".format(linreg.rvalue.round(3)))
plt.xlabel("Trans-eQTL Z score (JCTF)")
plt.ylabel("Trans-eQTL Z score (eQTLgen)")
plt.xlim([-15,15])
plt.ylim([-80,80])
plt.colorbar(sc, label='-log$_{10}$(p) in JCTF')
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/trans_enr_beta.png", dpi=500)
plt.clf()



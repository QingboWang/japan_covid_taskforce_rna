import pandas as pd
import numpy as np
import time as tm




#merge the GTEx p-value information
gpip = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip.tsv.gz", sep="\t")
gtex = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz", sep='\t')
#combine these two first. pval is enough
gpip.set_index(["variant_id","gene_id"], inplace=True)
gtex.set_index(["variant_id","gene_id"], inplace=True)
gpip = gpip.join(gtex.pval_nominal, how="left") #takes some time..
gpip.pval_nominal.fillna(1, inplace=True) #1 = non-significant
gpip.to_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip_and_signif_p.tsv.gz",
            sep="\t", compression="gzip")

gpip = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip_and_signif_p.tsv.gz",sep="\t")
print ("Done reading gpip, {0}".format(tm.ctime()))
#gpip["vid"] = gpip.variant_id.str.replace("_b38","")
gpip["vid"] = gpip.variant_id.str[:-4]#faster than above
print ("Done vid, {0}".format(tm.ctime()))
gpip["gid"] = gpip.gene_id.str.split("\\.").str[0] #to make sure about the match
print ("Done gid, {0}".format(tm.ctime()))
gpip.set_index(["vid", "gid"], inplace=True)
gpip.index.names = ["v","g"]
print ("Done set index, {0}".format(tm.ctime()))
#keep the chr names for later filtering
gpip_chrs = pd.Series(gpip.variant_id.str.split("_").str[0])
del gpip["variant_id"]
del gpip["gene_id"]
gpip.columns = gpip.columns + "_gtex"
print ("Done parsing gpip, {0}".format(tm.ctime()))

for n in range(25+1):
    print ("start read chunk {0}, {1}".format(n, tm.ctime()))
    jtex = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465.imputed.pairs.chunk{0}.txt".format(n),sep='\t')
    #filter to maf>0.01 as in GTEx protocol
    jtex = jtex[jtex.maf>0.01]
    jtex["vid"] = jtex.variant_id.str.split(":").apply(lambda x: "_".join(x[:4]))
    jtex["gid"] = jtex.gene_id.str.split("\\.").str[0]
    jtex.set_index(["vid", "gid"], inplace=True)
    jtex.index.names = ["v","g"]

    #filter gpip to small subset based on at least the chromosome
    chrs = jtex.variant_id.str.split(":").str[0].unique()
    filt_by_chr = gpip_chrs.apply(lambda x: x in chrs)
    gpip_subset = gpip[filt_by_chr]
    #and join
    jtex = jtex.join(gpip_subset, how="left")
    print ("done join, {0}".format(tm.ctime()))
    #also flip join (although not so many..)
    jtex["vid"] = jtex.variant_id.str.split(":").str[0] + "_" + \
                    jtex.variant_id.str.split(":").str[1] + "_" + \
                    jtex.variant_id.str.split(":").str[3] + "_" + \
                    jtex.variant_id.str.split(":").str[2]
    jtex["gid"] = jtex.gene_id.str.split("\\.").str[0]
    jtex.set_index(["vid", "gid"], inplace=True)
    jtex.index.names = ["v","g"]
    jtex = jtex.join(gpip_subset, how="left", rsuffix="_flip")
    jtex["flipped"] = ~jtex.pp_susie_gtex_flip.isna()
    print("done flip join, {0}".format(tm.ctime()))
    #now merge the original vs flip
    for col in jtex.columns:
        if (col + "_flip") in jtex.columns:
            jtex.loc[jtex[col].isna(),col] = jtex.loc[jtex[col].isna(), col+"_flip"]
            del jtex[col+"_flip"]
    print ("stats; ")
    print ("n(GTEx missing); {0}".format(sum(jtex.pp_susie_gtex.isna())))
    print ("n(total); {0}".format(jtex.shape[0]))
    print ("%GTEx missing; {0}".format(sum(jtex.pp_susie_gtex.isna()) / jtex.shape[0]))
    #and write
    jtex.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot.txt.gz".format(n),sep='\t', index=False, compression="gzip")

#check the concordance with GTEx:
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot.txt.gz".format(chk), sep='\t')
    df["pval_bin_gtex"] = "FDR threshold<p" #since we only have that data
    missing_in_gtex = (df.pp_susie_gtex.isna() & df.pp_fm_gtex.isna())
    df.loc[missing_in_gtex,"pval_bin_gtex"] = "missing"
    df.loc[(~missing_in_gtex) & (df.pval_nominal_gtex<5*10**-8),"pval_bin_gtex"] = "p<5e-8"
    df.loc[(~missing_in_gtex) & (df.pval_nominal_gtex>5*10**-8),"pval_bin_gtex"] = "5e-8<p<FDR threshold"
    tb = []
    for i in range(101):
        j = df[df.pval_nominal<10**-i].pval_bin_gtex.value_counts()
        tb.append(j)
        print ("done {0}".format(i))
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.iloc[:,[1,0,2,3]]#ordering manually
    tb.rename(columns = {"5e-8<p<FDR_thres":"5e-8<p<FDR threshold"}, inplace=True)
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_vs_gtex_pval_chunk{0}.tsv".format(chk), sep='\t')
    # vs gtex PIP
    df["min_pip_gtex"] = np.minimum(df.pp_susie_gtex.replace(0, 100),
                                    df.pp_fm_gtex.replace(0, 100))  # 100: just as a dummy for the undefined ones
    df["min_pip_gtex_bin"] = "missing"
    df.loc[df.min_pip_gtex >= 0, "min_pip_gtex_bin"] = "(0,0.001]"
    df.loc[df.min_pip_gtex > 0.001, "min_pip_gtex_bin"] = "(0.001,0.01]"
    df.loc[df.min_pip_gtex > 0.01, "min_pip_gtex_bin"] = "(0.01,0.1]"
    df.loc[df.min_pip_gtex > 0.1, "min_pip_gtex_bin"] = "(0.1,0.5]"
    df.loc[df.min_pip_gtex > 0.5, "min_pip_gtex_bin"] = "(0.5,0.9]"
    df.loc[df.min_pip_gtex > 0.9, "min_pip_gtex_bin"] = "(0.9,1]"
    df.loc[df.min_pip_gtex == 100, "min_pip_gtex_bin"] = "undefined"
    tb = []
    for i in range(101):
        j = df[df.pval_nominal < 10 ** -i].min_pip_gtex_bin.value_counts()
        tb.append(j)
        print("done {0}".format(i))
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.iloc[:, [0, 1, 3, 2, 4, 5, 6, 7]]  # ordering manually
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_gtex_pip_chunk{0}.tsv".format(chk), sep='\t')
#Summary:
chk = 0
tbp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_vs_gtex_pval_chunk{0}.tsv".format(chk),sep='\t', index_col=0)
tbpip = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_gtex_pip_chunk{0}.tsv".format(chk), sep='\t', index_col=0)
colp = tbp.columns
colpip = tbpip.columns
for chk in range(1,26):
    tbpsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_vs_gtex_pval_chunk{0}.tsv".format(chk), sep='\t')
    tbpipsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_gtex_pip_chunk{0}.tsv".format(chk),sep='\t')
    tbp = tbp + tbpsub.loc[tbp.index, colp]
    tbpip = tbpip + tbpipsub.loc[tbpip.index, colpip]
    print ("done chunk {0}, {1}".format(chk, tm.ctime()))
tbp.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_vs_gtex_pval_combined.tsv",sep='\t')
tbpip.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_gtex_pip_combined.tsv", sep='\t')
#sanity check;
print (tbp.apply(lambda x: x/sum(x), axis=1))
print (tbpip.apply(lambda x: x/sum(x), axis=1))


#merge the GTEx PIP comparison:
fm = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_n465_imputed_all_v.txt", sep=' ', index_col=0)
print ("read fm, {0} lines, {1}".format(fm.shape[0], tm.ctime()))
sus = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_n465_imputed_all_v.tsv", sep='\t', index_col=0)
print ("read sus, {0} lines, {1}".format(sus.shape[0], tm.ctime()))
pips = fm.join(sus, how="left")
print ("concatenated the two {0}".format(tm.ctime()))
pips.index = pips.index.str.split("_", expand=True)
pips.index.names = ["variant_id","gene_id"]
pips.columns = ["pip_fm", "pip_susie"]
print ("done index split. pip file is ready. {0}".format(tm.ctime()))
pips.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_all_v.tsv", sep='\t')
print ("and wrote the file. {0}".format(tm.ctime())) #do this later probably?
#checking NAs in SuSiE;
pips[pips.pip_susie.isna()]
pips.loc[pips.pip_susie.isna(), "pip_susie"] = pips.loc[pips.pip_susie.isna(), "pip_fm"] #if SiSiE was missing, impute with FM result
#and go annotate:
for chk in range(26):
    print("reading the chunk {0}, {1}".format(chk, tm.ctime()))
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot.txt.gz".format(chk), sep='\t')
    print ("read the chunk {0}, {1}".format(chk, tm.ctime()))
    df.set_index(["variant_id", "gene_id"], inplace=True)
    df = df.join(pips, how="left")
    print("done join, {0}".format(tm.ctime()))
    if chk==0: #sanity check
        print (df.head(10))
        print (df.tail(10))
    df.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')


#Checking the concordance in the level of PIP as well:
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["pval_bin_gtex"] = "FDR threshold<p" #since we only have that data
    missing_in_gtex = (df.pp_susie_gtex.isna() & df.pp_fm_gtex.isna())
    df.loc[missing_in_gtex,"pval_bin_gtex"] = "missing"
    df.loc[(~missing_in_gtex) & (df.pval_nominal_gtex<5*10**-8),"pval_bin_gtex"] = "p<5e-8"
    df.loc[(~missing_in_gtex) & (df.pval_nominal_gtex>5*10**-8),"pval_bin_gtex"] = "5e-8<p<FDR threshold"
    #j PIP vs gtex pval
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),df.pip_fm.fillna(-1))  # -1 = those that are not fine-mapped since it is not egene
    tb = []
    for i in np.arange(101) / 100:
        j = df[df.min_pip >= i].pval_bin_gtex.value_counts()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.loc[:,["missing","FDR threshold<p","5e-8<p<FDR threshold","p<5e-8"]]#ordering manually
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pval_chunk{0}.tsv".format(chk), sep='\t')

    # vs gtex PIP
    df["min_pip_gtex"] = np.minimum(df.pp_susie_gtex.replace(0, 100),
                                    df.pp_fm_gtex.replace(0, 100))  # 100: Dummy
    df["min_pip_gtex_bin"] = "missing"
    df.loc[df.min_pip_gtex >= 0, "min_pip_gtex_bin"] = "(0,0.001]"
    df.loc[df.min_pip_gtex > 0.001, "min_pip_gtex_bin"] = "(0.001,0.01]"
    df.loc[df.min_pip_gtex > 0.01, "min_pip_gtex_bin"] = "(0.01,0.1]"
    df.loc[df.min_pip_gtex > 0.1, "min_pip_gtex_bin"] = "(0.1,0.5]"
    df.loc[df.min_pip_gtex > 0.5, "min_pip_gtex_bin"] = "(0.5,0.9]"
    df.loc[df.min_pip_gtex > 0.9, "min_pip_gtex_bin"] = "(0.9,1]"
    df.loc[df.min_pip_gtex == 100, "min_pip_gtex_bin"] = "undefined"
    tb = []
    for i in np.arange(101) / 100:
        j = df[df.min_pip >= i].min_pip_gtex_bin.value_counts()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.loc[:, ["missing","undefined","(0,0.001]","(0.001,0.01]","(0.01,0.1]","(0.1,0.5]","(0.5,0.9]","(0.9,1]"]]  # ordering manually
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_chunk{0}.tsv".format(chk), sep='\t')

    #get the min and bin the jtex p-val, for j vs j comparison
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),df.pip_fm.fillna(-1))  # -1 = those that are not into fine-mapping pipeline
    df["min_pip_bin"] = "undefined"
    df.loc[df.min_pip >= 0, "min_pip_bin"] = "(0,0.001]"
    df.loc[df.min_pip > 0.001, "min_pip_bin"] = "(0.001,0.01]"
    df.loc[df.min_pip > 0.01, "min_pip_bin"] = "(0.01,0.1]"
    df.loc[df.min_pip > 0.1, "min_pip_bin"] = "(0.1,0.5]"
    df.loc[df.min_pip > 0.5, "min_pip_bin"] = "(0.5,0.9]"
    df.loc[df.min_pip > 0.9, "min_pip_bin"] = "(0.9,1]"
    #j pval vs j pip
    tb = []
    for i in range(101):
        j = df[df.pval_nominal < 10 ** -i].min_pip_bin.value_counts()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.loc[:, ["undefined","(0,0.001]","(0.001,0.01]","(0.01,0.1]","(0.1,0.5]","(0.5,0.9]","(0.9,1]"]]  # ordering manually
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_jtex_pip_chunk{0}.tsv".format(chk), sep='\t')
    #j pip vs j pval
    df["pval_bin"] = "[0.05,1]"
    df.loc[df.pval_nominal < 5*10**-2,"pval_bin"] = "[5e-8,0.05)"
    df.loc[df.pval_nominal_gtex<5*10**-8,"pval_bin"] = "< 5e-8"
    tb = []
    for i in np.arange(101) / 100:
        j = df[df.min_pip >= i].pval_bin.value_counts()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.loc[:, ["[0.05,1]","[5e-8,0.05)","< 5e-8"]]  # ordering manually -- check and change things
    tb = tb.fillna(0).astype(int)
    print (chk)
    print (tb)
    print (tb.apply(lambda x: x/sum(x), axis=1))
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_jtex_p_chunk{0}.tsv".format(chk), sep='\t')


#merging stats for jtex pip vs gtex comparison results
chk = 0
tbp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pval_chunk{0}.tsv".format(chk),sep='\t', index_col=0)
tbpip = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_chunk{0}.tsv".format(chk), sep='\t', index_col=0)
colp = tbp.columns
colpip = tbpip.columns
for chk in range(1,26):
    tbpsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pval_chunk{0}.tsv".format(chk), sep='\t')
    tbpipsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_chunk{0}.tsv".format(chk),sep='\t')
    tbp = tbp + tbpsub.loc[tbp.index, colp]
    tbpip = tbpip + tbpipsub.loc[tbpip.index, colpip]
    print ("done chunk {0}, {1}".format(chk, tm.ctime()))
tbp.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pval_combined.tsv",sep='\t')
tbpip.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_gtex_pip_combined.tsv", sep='\t')
#sanity check;
print (tbp.apply(lambda x: x/sum(x), axis=1))
print (tbpip.apply(lambda x: x/sum(x), axis=1))

#merging stats for jtex vs jtex itself
chk = 0
tbp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_jtex_p_chunk{0}.tsv".format(chk),sep='\t', index_col=0)
tbpip = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_jtex_pip_chunk{0}.tsv".format(chk), sep='\t', index_col=0)
colp = tbp.columns
colpip = tbpip.columns
for chk in range(1,26):
    tbpsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_jtex_p_chunk{0}.tsv".format(chk), sep='\t')
    tbpipsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_jtex_pip_chunk{0}.tsv".format(chk),sep='\t')
    tbp = tbp + tbpsub.loc[tbp.index, colp]
    tbpip = tbpip + tbpipsub.loc[tbpip.index, colpip]
    print ("done chunk {0}, {1}".format(chk, tm.ctime()))
tbp.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_jtex_p_combined.tsv",sep='\t')
tbpip.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_jtex_pip_combined.tsv", sep='\t')


#Expression Modifier Score comparison (JTEx p vs PIP)
import numpy as np
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    #j PIP vs EMS
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),df.pip_fm.fillna(-1))  # -1 = those that are not fine-mapped since it is not egene
    tb = []
    for i in np.arange(101) / 100:
        j = df[df.min_pip >= i].ems_bin_gtex.value_counts().sort_index()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.fillna(0).astype(int)
    print (chk)
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_ems_chunk{0}.tsv".format(chk), sep='\t')
    #j pval vs EMS
    tb = []
    for i in range(101):
        j = df[df.pval_nominal < 10 ** -i].ems_bin_gtex.value_counts().sort_index()
        tb.append(j)
    tb = pd.DataFrame(tb)
    tb.index = range(101)
    tb = tb.fillna(0).astype(int)
    print (chk)
    tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_ems_chunk{0}.tsv".format(chk), sep='\t')
#merge the result
chk = 0
tbp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_ems_chunk{0}.tsv".format(chk),sep='\t', index_col=0)
tbpip = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_ems_chunk{0}.tsv".format(chk), sep='\t', index_col=0)
colp = tbp.columns
colpip = tbpip.columns
for chk in range(1,26):
    tbpsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_ems_chunk{0}.tsv".format(chk), sep='\t',index_col=0)
    tbpipsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_ems_chunk{0}.tsv".format(chk),sep='\t', index_col=0)
    tbp = tbp + tbpsub.loc[tbp.index, colp]
    tbpip = tbpip + tbpipsub.loc[tbpip.index, colpip]
    print ("done chunk {0}, {1}".format(chk, tm.ctime()))
tbp.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_pip_vs_ems_combined.tsv",sep='\t')
tbpip.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_ems_combined.tsv", sep='\t')


#agreement between SuSiE and FM
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["pip_fm_bin"] = 0
    df.loc[df.pip_fm>=0,"pip_fm_bin"] = "(0,0.001]"
    df.loc[df.pip_fm>0.001,"pip_fm_bin"] = "(0.001,0.01]"
    df.loc[df.pip_fm>0.01,"pip_fm_bin"] = "(0.01,0.1]"
    df.loc[df.pip_fm>0.1,"pip_fm_bin"] = "(0.1,0.5]"
    df.loc[df.pip_fm>0.5,"pip_fm_bin"] = "(0.5,0.9]"
    df.loc[df.pip_fm>0.9,"pip_fm_bin"] = "(0.9,1]"
    df["pip_susie_bin"] = 0
    df.loc[df.pip_susie>=0,"pip_susie_bin"] = "(0,0.001]"
    df.loc[df.pip_susie>0.001,"pip_susie_bin"] = "(0.001,0.01]"
    df.loc[df.pip_susie>0.01,"pip_susie_bin"] = "(0.01,0.1]"
    df.loc[df.pip_susie>0.1,"pip_susie_bin"] = "(0.1,0.5]"
    df.loc[df.pip_susie>0.5,"pip_susie_bin"] = "(0.5,0.9]"
    df.loc[df.pip_susie>0.9,"pip_susie_bin"] = "(0.9,1]"
    tb = df.groupby(["pip_susie_bin","pip_fm_bin"]).size().unstack(level=-1).fillna(0) #0 = failed, or not fine mapped in that platform (somehow..)
    tb = tb.iloc[1:,1:].astype(int)
    if chk==0:
        tbsum = tb
    else:
        tbsum = tbsum + tb
    print("done {0}, {1}".format(chk, tm.ctime()))
tbsum.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/susie_vs_pip.tsv", sep='\t')

#FM vs SuSiE in terms of EMS as well
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["pip_fm_bin"] = 0
    df.loc[df.pip_fm>=0,"pip_fm_bin"] = "(0,0.01]"
    df.loc[df.pip_fm>0.01,"pip_fm_bin"] = "(0.01,0.1]"
    df.loc[df.pip_fm>0.1,"pip_fm_bin"] = "(0.1,0.9]"
    df.loc[df.pip_fm>0.9,"pip_fm_bin"] = "(0.9,1]"
    df["pip_susie_bin"] = 0
    df.loc[df.pip_susie>=0,"pip_susie_bin"] = "(0,0.01]"
    df.loc[df.pip_susie>0.01,"pip_susie_bin"] = "(0.01,0.1]"
    df.loc[df.pip_susie>0.1,"pip_susie_bin"] = "(0.1,0.9]"
    df.loc[df.pip_susie>0.9,"pip_susie_bin"] = "(0.9,1]"
    tb = df.groupby(["pip_susie_bin","pip_fm_bin","ems_bin_gtex"]).size().unstack(level=-1).fillna(0) #0 = failed, or not fine mapped in that platform (somehow..)
    tb = tb.iloc[1:,:].astype(int)
    if chk==0:
        tbsum = tb
    else:
        tbsum = tbsum.add(tb, fill_value=0) #to account for missingness
    print("done {0}, {1}".format(chk, tm.ctime()))
tbsum.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/susie_vs_pip_ems.tsv", sep='\t') #unstack again when plotting


#Effect size concordance:
gtex = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz", sep='\t')
gtex["gene_id_unv"] = gtex.gene_id.str.split("\\.").str[0]
gtex.index = gtex[["variant_id", "gene_id_unv"]]
#for chk in range(26):
for chk in range(23,26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["gene_id_unv"] = df.gene_id.str.split("\\.").str[0]
    df["variant_b38"] = df.variant_id.str.split(":").str[0] + "_" + \
                        df.variant_id.str.split(":").str[1] + "_" + \
                        df.variant_id.str.split(":").str[2] + "_" + \
                        df.variant_id.str.split(":").str[3] + "_b38"
    df["variant_b38_flip"] = df.variant_id.str.split(":").str[0] + "_" + \
                        df.variant_id.str.split(":").str[1] + "_" + \
                        df.variant_id.str.split(":").str[3] + "_" + \
                        df.variant_id.str.split(":").str[2] + "_b38"
    df.index = df[["variant_b38","gene_id_unv"]]
    df.index.names = gtex.index.names
    df = df.join(gtex.slope, how="left", rsuffix="_gtex")
    df.index = df[["variant_b38_flip","gene_id_unv"]]
    df.index.names = gtex.index.names
    df = df.join(gtex.slope, how="left", rsuffix="_gtex_flip")
    df = df[["variant_id","gene_id","slope_gtex", "slope_gtex_flip"]]
    df.rename(columns = {"slope_gtex":"beta_gtex", "slope_gtex_flip":"beta_gtex_flip"}, inplace=True)
    df.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtex_beta.txt".format(chk), sep='\t', index=False)
    print ("done chk{0}, {1}".format(chk, tm.ctime()))
#merge and get the eff, pval_nominal for both, and pip for both
dfbeta = []
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    dfb = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtex_beta.txt".format(chk), sep='\t')
    dfb.loc[dfb.beta_gtex.isna(), "beta_gtex"] = dfb.loc[dfb.beta_gtex.isna(), "beta_gtex_flip"]*-1
    df["beta_gtex"] = np.array(dfb.beta_gtex)
    dfbeta.append(df[~df.beta_gtex.isna()])
    print ("done merging beta chk{0}, {1}".format(chk, tm.ctime()))
pd.concat(dfbeta).to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465.imputed.pairs.gtex_beta.txt", sep='\t', index=False)
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465.imputed.pairs.gtex_beta.txt", sep='\t')
#actually, let's also restrict to our favorite p<5e-8 for GTEx and do the comparison
df["sign_same"] = (df.beta_gtex*df.slope>0)
tb = []
for i in range(101):
    j = df[df.pval_nominal < 10 ** -i].sign_same.value_counts()
    tb.append(j)
tb = pd.DataFrame(tb)
tb.index = range(101)
tb = tb.fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_beta_fdr005.tsv".format(chk),sep='\t')
tb = []
df08 = df[df.pval_nominal_gtex<5*10**-8]
for i in range(101):
    j = df08[df08.pval_nominal < 10 ** -i].sign_same.value_counts()
    tb.append(j)
tb = pd.DataFrame(tb)
tb.index = range(101)
tb = tb.fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
tb.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/jtex_p_vs_beta_5e8.tsv".format(chk),sep='\t')


#SuRE enrichment:
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt
gtex_sure = {}
jtex_sure_miss = {}
jtex_sure_exists = {}
import time as tm
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["hg19"] = df.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))
    df["min_pip_gtex"] = np.minimum(df.pp_susie_gtex.replace(0, 100),
                                    df.pp_fm_gtex.replace(0, 100))  # 100はGTExに存在するけどPIP対象になっていないgenes.
    df["min_pip_gtex_bin"] = "missing"
    df.loc[df.min_pip_gtex >= 0, "min_pip_gtex_bin"] = "(0,0.01]"
    df.loc[df.min_pip_gtex > 0.01, "min_pip_gtex_bin"] = "(0.01,0.1]"
    df.loc[df.min_pip_gtex > 0.1, "min_pip_gtex_bin"] = "(0.1,0.9]"
    df.loc[df.min_pip_gtex > 0.9, "min_pip_gtex_bin"] = "(0.9,1]"
    df.loc[df.min_pip_gtex == 100, "min_pip_gtex_bin"] = "Undefined"
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),
                               df.pip_fm.fillna(-1))  # -1 = those that are not fine-mapped since it is not egene
    df["min_pip_bin"] = "Undefined"
    df.loc[df.min_pip >= 0, "min_pip_bin"] = "(0,0.01]"
    df.loc[df.min_pip > 0.01, "min_pip_bin"] = "(0.01,0.1]"
    df.loc[df.min_pip > 0.1, "min_pip_bin"] = "(0.1,0.9]"
    df.loc[df.min_pip > 0.9, "min_pip_bin"] = "(0.9,1]"
    dfmiss = df[df.min_pip_gtex_bin=="missing"][["min_pip_bin","min_pip_gtex_bin","hg19"]]
    dfex = df[~(df.min_pip_gtex_bin=="missing")][["min_pip_bin","min_pip_gtex_bin","hg19"]]
    for b in df.min_pip_gtex_bin.unique():
        v = df[df.min_pip_gtex_bin==b].hg19
        n = len(np.intersect1d(v, raqtls))
        if chk==0:
            gtex_sure[b] = n
        else:
            gtex_sure[b] = gtex_sure[b]+n
        #jtex, missing in gtex
        v = dfmiss[dfmiss.min_pip_bin==b].hg19
        n = len(np.intersect1d(v, raqtls))
        if chk==0:
            jtex_sure_miss[b] = n
        else:
            jtex_sure_miss[b] = jtex_sure_miss[b]+n
        #jtex, exist in gtex
        v = dfex[dfex.min_pip_bin==b].hg19
        n = len(np.intersect1d(v, raqtls))
        if chk==0:
            jtex_sure_exists[b] = n
        else:
            jtex_sure_exists[b] = jtex_sure_exists[b]+n
    print ("done {0}, {1}".format(chk, tm.ctime()))
df = pd.concat([pd.DataFrame(gtex_sure, index=["g"]),
    pd.DataFrame(jtex_sure_exists, index=["j_g_exist"]),
    pd.DataFrame(jtex_sure_miss, index=["j_g_missing"])], axis=0)
df.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/sure_enr.tsv", sep="\t")


#ROC and PRC for predicting GTEx putative causal ones
st = []
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    df["min_pip_gtex"] = np.minimum(df.pp_susie_gtex.replace(0, 100),
                                    df.pp_fm_gtex.replace(0, 100))  # 100はGTExに存在するけどPIP対象になっていないgenes.
    df["min_pip"] = np.minimum(df.pip_susie.fillna(-1),df.pip_fm.fillna(-1))  # -1 = those that are not fine-mapped since it is not egene
    df_pippos = df[(df.min_pip_gtex>0.9)&(df.min_pip_gtex!=100)]
    df_pipneg = df[(df.min_pip_gtex<0.001)&(df.min_pip_gtex!=100)]
    y = [1]*df_pippos.shape[0] + [0]*df_pipneg.shape[0]
    x1 = np.array(list(df_pippos.pval_nominal) + list(df_pipneg.pval_nominal))#pval
    x2 = np.array(list(df_pippos.min_pip) + list(df_pipneg.min_pip))#pip
    y = np.array(y)
    x1 = x1[~np.isnan(x2)]
    y = y[~np.isnan(x2)]
    x2 = x2[~np.isnan(x2)] #removing those where PIP is not defined at all in jtex
    x1 = -1*np.log10(x1) #-log10, the higher the better for prediction
    st.append(pd.DataFrame({"gtex_causal":y, "pip":x2, "pval":x1}))
    print("done {0}, {1}".format(chk, tm.ctime()))
st = pd.concat(st)
st.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/basic_stats/predicting_gtex_putative_causal.tsv", sep='\t', index=False)

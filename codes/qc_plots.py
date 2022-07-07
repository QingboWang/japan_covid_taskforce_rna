
qc = pd.read_csv("~/Downloads/rnaseqc_metrics_n500.tsv", sep='\t', index_col=0) #Replace this with QC metrics matrix generated with RNAseQC for your own data
r = pd.read_csv("~/Downloads/cnt_n500_corrmat.tsv", sep='\t', index_col=0) #Replace this with the expression correlation between counts for your own data

from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

plt.figure(figsize=(4.5,3))
p = qc.Mapped.hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x() < 10000000:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=10000000, linestyle="--", linewidth=1, color="black")
plt.xlabel("Number of mapped reads")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_mappedreads.png", dpi=500)
plt.clf()

plt.figure(figsize=(4.5,3))
p = qc["Mapping Rate"].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x() < 0.8:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.8, linestyle="--", linewidth=1, color="black")
plt.xlabel("Mapping rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_maprate.png", dpi=500)
plt.clf()

plt.figure(figsize=(4.5,3))
p = qc['Intergenic Rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x() > 0.1:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.1, linestyle="--", linewidth=1, color="black")
plt.xlabel("Intergenic rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_intgenic.png", dpi=500)
plt.clf()

plt.figure(figsize=(4.5,3))
p = qc['Base Mismatch Rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x() > 0.01:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.01, linestyle="--", linewidth=1, color="black")
plt.xlabel("Base mismatch rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_basemismatch.png", dpi=500)
plt.clf()

plt.figure(figsize=(4.5,3))
p = qc['rRNA rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x() > 0.3:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.3, linestyle="--", linewidth=1, color="black")
plt.xlabel("Ribosomal RNA rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_rrna.png", dpi=500)
plt.clf()

#and the r prime
r_i = (r.sum(axis=0)-1)/(r.shape[1]-1)
r_barbar = (r.sum().sum()-r.shape[0])/(r.shape[0]**2-r.shape[0])
D_i = (r_i - r_barbar)/(np.median(r_i - r_barbar))

plt.figure(figsize=(4.5,3))
p = D_i.hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if rectangle.get_x()<-15:
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=-15, linestyle="--", linewidth=1, color="black")
plt.xlabel("Intersample correlation deviation ($D_i$)")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_di.png", dpi=500)
plt.clf()

#Below are the filtering criteria applied:
"""
filtered_by_n_mappedread = dfs.Mapped<10000000
filtered_by_mapping_rate = dfs['Mapping Rate']<0.2
filtered_by_intergenic_mapping_rate = dfs['Intergenic Rate']>0.3
filtered_by_base_mismatch_rate = dfs['Base Mismatch Rate']>0.01
filtered_by_rrna_rate = dfs['rRNA rate']>0.3
filtered_by_r_i= D_i < -15
"""

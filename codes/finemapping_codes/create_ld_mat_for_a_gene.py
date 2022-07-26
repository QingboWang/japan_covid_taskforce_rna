#input = the file, and also the gene z file as a reference
import sys
import pandas as pd
import time as tm
gene = sys.argv[1]
zsc_dir = sys.argv[2]
df_dir = sys.argv[3]

z = pd.read_csv(zsc_dir, sep=' ')

df = pd.read_csv(df_dir, sep='\t', header=None)
df.index = df.iloc[:,0] + ":" + df.iloc[:,1].astype(str) + ":" + \
           df.iloc[:,3] + ":" + df.iloc[:,4] + ":" + df.iloc[:,2]
#remove maf0 and set the columns identical as the .z file
df = df.loc[z.rsid.str.split("_").str[0],:]
#and turn it to dosage;
df = df.iloc[:,9:]
df_dosage = df.applymap(lambda x: float(x.split(":")[-1]))
print ("done getting dosage mat for gene{0}: {1}".format(gene, tm.ctime()))

#then calculate the LD
M = pd.read_csv("~/Dropbox/my_shared_resources/n465_imputed_M_for_covadj.tsv", sep="\t", index_col=0)

X = df_dosage.T #the matrix X we want, pre-normalization etc.
X = X.apply(lambda l: (l - l.mean())/l.std(), axis=0) #standardization
X.index = M.index
X_adj = M @ X
#X_adj.T.to_csv("/Users/qingbowang/Desktop/tmp/{0}_X_adj_T.tsv.gz".format(gene), sep="\t", compression="gzip") #no need. Takes time.

x = X_adj.T
n = x.shape[1]#number of samples = 465
R = x @ x.T / n
print ("done calculating LD mat for gene{0}: {1}".format(gene, tm.ctime()))
print (R.shape)
print (R.iloc[:3,:3])
R.to_csv("~/Desktop/imputed_ld_mat/{0}.ld".format(gene), sep=" ", index=False, header=None) #結局律速はここ、read/write

print ("done writing adj-LD mat for gene{0}: {1}".format(gene, tm.ctime()))

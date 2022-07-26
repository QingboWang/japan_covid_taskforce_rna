import pandas as pd
import numpy as np
cov = pd.read_csv("~/Downloads/n465.imputed.combined_covariates_idadded.txt", sep="\t", index_col=0)
cov = cov.T #row = samples, col = covariates. For Xs too, row = samples
I = pd.DataFrame(np.identity(cov.shape[0]), cov.index, cov.index)
inside_invterm = cov.T@cov
invterm = pd.DataFrame(np.linalg.inv(inside_invterm.values), inside_invterm.columns, inside_invterm.index)
#sanity check
print (invterm @ inside_invterm) #and this is also right
M = I - cov @ invterm @ cov.T
M.to_csv("~/Dropbox/my_shared_resources/n465_imputed_M_for_covadj.tsv", sep="\t")

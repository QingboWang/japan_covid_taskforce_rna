library(susieR)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
gn = args[1]

ld = paste0("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_inputs_n465_imputed/",gn,"/",gn,".ld")
R = as.matrix(fread(ld))
print ("done reading ld in R")
print (Sys.time())
z = paste0("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_inputs_n465_imputed/",gn,"/",gn,"_fminput.z")
st = read.table(z, header=T)
st$z = st$beta / st$se

#filter to maf>0.01
R = R[st$maf>0.01,st$maf>0.01]
st = st[st$maf>0.01,]

z_scores = st$z
print ("Starting SuSiE")
print (Sys.time())
fitted_rss = susie_rss(z_scores, R, L=10, refine=F)
print ("Done SuSiE, starting writing results")
print (Sys.time())
#write the results
cs_out = paste0("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_inputs_n465_imputed/",gn,"/",gn,"_susie_cs.txt")
write.table(summary(fitted_rss)$cs , cs_out, quote=F, sep='\t')
pip_out = paste0("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_inputs_n465_imputed/",gn,"/",gn,"_susie_pip.tsv")
st$pip = fitted_rss$pip
write.table(st[c("rsid","pip")], pip_out, quote=F, sep='\t')
print ("Done writing")
print (Sys.time())


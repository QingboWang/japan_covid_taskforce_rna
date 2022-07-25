library(DESeq2)

#read in the data
df = read.table("/Users/qingbowang/Desktop/covid_rna_de/cnt_n500_int.tsv.gz", sep='\t', header=T)#created integer version for convenience in python.
pheno = read.table("/Users/qingbowang/Desktop/covid_rna_de/Pheno_data_SAIGE.txt", sep='\t', header=T)
samples = read.table("/Users/qingbowang/Desktop/covid_rna_de/n465_samples.txt", sep="\t", header=F)
samples = samples$V1
#also the list of genes to test
gtex_exp = read.table("/Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz", sep='\t')
genes_to_use = gtex_exp$V4
rownames(df) = df$gene_id
df = df[genes_to_use,samples]
rownames(pheno) = pheno$IID
pheno = pheno[samples,]
sample_df = pheno[,c("IID", "Pca_4", "Sex","Age")]
colnames(sample_df) = c("samplename","case_status", "Sex", "Age")

#run edgeR
library(edgeR)
y = DGEList(counts=df, genes=rownames(df))
y <- calcNormFactors(y, method="TMM")

#first, binary comparison of 12 vs 34
sex = pheno$Sex
age = pheno$Age
case = pheno$Pca_1
design = model.matrix(~sex+age+case) #the interesting variable needs to be the last
rownames(design) = colnames(y)
y = estimateDisp(y, design, robust=TRUE)
fit = glmFit(y, design)
#checking the result
lrt = glmLRT(fit)
write.table(lrt$table, "/Users/qingbowang/Desktop/covid_rna_de/n465_edger_lrt_12_vs_34_results.tsv", quote=F, sep='\t')
#then, comparison of 1 vs 4
y = DGEList(counts=df[,to_use], genes=rownames(df))
y <- calcNormFactors(y, method="TMM")
sex = pheno[to_use,]$Sex
age = pheno[to_use,]$Age
case = pheno[to_use,]$Pca_1
design = model.matrix(~sex+age+case) #the interesting variable needs to be the last
rownames(design) = colnames(y)
y = estimateDisp(y, design, robust=TRUE)
fit = glmFit(y, design)
#checking the result
lrt = glmLRT(fit)
write.table(lrt$table, "/Users/qingbowang/Desktop/covid_rna_de/n465_edger_lrt_1_vs_4_results.tsv", quote=F, sep='\t')

#continuous
y = DGEList(counts=df, genes=rownames(df))
y <- calcNormFactors(y, method="TMM")
sex = pheno$Sex
age = pheno$Age
case = pheno$Pca_4
design = model.matrix(~sex+age+case) #the interesting variable needs to be the last
rownames(design) = colnames(y)
y = estimateDisp(y, design, robust=TRUE)
fit = glmFit(y, design)
#checking the result
lrt = glmLRT(fit)
write.table(lrt$table, "/Users/qingbowang/Desktop/covid_rna_de/n465_edger_lrt_cont_results.tsv", quote=F, sep='\t')

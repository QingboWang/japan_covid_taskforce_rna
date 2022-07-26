l=$(gunzip -c /Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz | head -n $n | tail -n 1 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Dropbox/ct_filtered_vcf/ct_imputed_hg38_sorted_"$chr".vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/vcf_"$gn".tsv #output
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
gene_fn=/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/teststats_per_gene_n465/"$gn"_fminput.z #z score file as a reference to match the index
python /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/create_ld_mat_for_a_gene.py $gn $gene_fn $fn_tmp
rm $fn_tmp

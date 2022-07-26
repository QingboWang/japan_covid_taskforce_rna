
cd /Users/qingbowang/Desktop/covid_eqtl_imputed_fm
gene_names=gwsig_genes_ids.tsv


#writing the header:
head -n 1 fm_inputs_n465_imputed/ENSG00000000460.17/ENSG00000000460.17.snp | cut -d " " -f2,11 > fm_n465_imputed_all_v.txt
head -n 1 fm_inputs_n465_imputed/ENSG00000000460.17/ENSG00000000460.17_susie_pip.tsv > susie_n465_imputed_all_v.tsv
#write everything into single file
i=0
while IFS= read -r line
do
    tail -n +2 fm_inputs_n465_imputed/"$line"/"$line".snp | cut -d " " -f2,11 >> fm_n465_imputed_all_v.txt
    tail -n +2 fm_inputs_n465_imputed/"$line"/"$line"_susie_pip.tsv | cut -f2,3 >> susie_n465_imputed_all_v.tsv
    i=$(($i+1))
    if [ $(($i%100)) -eq 0 ]
    then
      now=$(date +"%T")
      echo "done $i $now"
    fi
done < "$gene_names"

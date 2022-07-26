#get the $1-th gene name etc in the expression file

#print the num
echo "Starting $1-th line"

cd /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/
l=$(gunzip -c /Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz | head -n $1 | tail -n 1 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
gene_fn_beforefilt=/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/teststats_per_gene_n465/"$gn"_fminput.z #z score file as a reference to match the index
#if the z score file is missing, that means I didn't make the z-score file since the min(P)>5e-8 -> skip
#also for re-doing, only proceed when the fm and susie output files are missing
if [ ! -f ./fm_inputs_n465_imputed/"$gn"/"$gn".snp ] || [ ! -f ./fm_inputs_n465_imputed/"$gn"/"$gn"_susie_pip.tsv ]
then
  if test -f "$gene_fn_beforefilt"
  then
    now=$(date +"%T")
    echo "Starting $gn : $now" >> /Users/qingbowang/Desktop/imputed_finemap_startlog.txt
    echo "Starting to filter the .z file based on maf>0.01 : $now"
    awk -F'[ ]' ' $6 >= 0.01 ' "$gene_fn_beforefilt" > /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/zfile_maffilt/"$gn"_fminput.z
    gene_fn=/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/zfile_maffilt/"$gn"_fminput.z #this is the filtered .z file.
    now=$(date +"%T")
    echo "Starting to create LD : $now"
    fn=/Users/qingbowang/Dropbox/ct_filtered_vcf/ct_imputed_hg38_sorted_"$chr".vcf.gz #the full vcf file to create LD mat from
    fn_tmp=/Users/qingbowang/Desktop/tmp/vcf_"$gn".tsv #to cut the vcf,
    /Users/qingbowang/samtools-1.13/htslib-1.13/tabix $fn "$chr":$(($tss-1000001))-$(($tss+1000001)) > $fn_tmp #for +-1Mb of TSS
    #wait if the memory usage is high. the next step, creating LD mat, requires tons of memory
    while :
    do
      mem_avail=$(memory_pressure | tail -n 1 | cut -d " " -f5 | cut -d "%" -f1)
      if [ "$mem_avail" -gt 25 ]; then #wait when less than 25% System-wide memory free percentage
        break
      fi
      now=$(date +"%T")
      echo "Waiting $1 th line. Only $mem_avail percent is currently available.. $(date +"%T")"
      sleep 30 #30 second sleep.
    done
    python /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/create_ld_mat_for_a_gene.py $gn $gene_fn $fn_tmp #and create LD matrix
    rm $fn_tmp #remove since it is so heavy
    #put things to the right file and run FINEMAP
    mkdir -p ./fm_inputs_n465_imputed/"$gn"
    mv zfile_maffilt/"$gn"_fminput.z ./fm_inputs_n465_imputed/"$gn" #the one filtered to maf>001
    mv /Users/qingbowang/Desktop/imputed_ld_mat/"$gn".ld ./fm_inputs_n465_imputed/"$gn"
    cd ./fm_inputs_n465_imputed/"$gn"
    echo "z;ld;snp;config;cred;log;n_samples" > ./master
    echo "${gn}_fminput.z;${gn}.ld;${gn}.snp;${gn}.config;${gn}.cred;${gn}.log;465" >> ./master
    ~/Downloads/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX n --sss --in-files ./master
    now=$(date +"%T")
    echo "Done FINEMAP $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/imputed_finemap_endlog_fm.txt
    cd ../../
    #run SuSiE
    #but wait if the memory usage is high.
    while :
    do
      mem_avail=$(memory_pressure | tail -n 1 | cut -d " " -f5 | cut -d "%" -f1)
      if [ "$mem_avail" -gt 25 ]; then #wait when less than 25% System-wide memory free percentage
        break
      fi
      now=$(date +"%T")
      echo "Waiting $1 th line. Only $mem_avail percent is currently available.. $(date +"%T")"
      sleep 30 #30 second sleep.
    done
    echo "Starting to run SuSiE : $now"
    Rscript ./run_susie_for_a_gene.R $gn #put this in a proper path = /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/
    now=$(date +"%T")
    echo "Done running SuSiE (l=$1) : $now"
    echo "Done SuSiE $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/imputed_finemap_endlog_susie.txt
    #remove the LD mat since it is so heavy
    rm ./fm_inputs_n465_imputed/"$gn"/"$gn".ld
    #print things to the log
  else
    now=$(date +"%T")
    echo "$gn is not a eGene, skipping (l=$1)"
    echo "$gn is not a eGene, skipping (l=$1): $now" >> /Users/qingbowang/Desktop/imputed_finemap_skiplog.txt
  fi
fi

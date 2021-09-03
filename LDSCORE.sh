#!/bin/bash

#https://github.com/bulik/ldsc
#https://data.broadinstitute.org/alkesgroup/LDSCORE/
#https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
#https://github.com/bulik/ldsc/wiki/FAQ

#Path to software
ldsc_soft=$1

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2

#Prepare the files
input_file=$2
zcat $input_file |  awk 'OFS="\t"{print $2,$3,$4,$5,$6,$10,$14,$15,$16,$19}' >  forldsc_AllOA.txt
#Repeat for all phenos

#Path to anacoda
anacoda=$3
PATH=$anacoda:$PATH
$ldsc_soft/munge_sumstats.py --sumstats forldsc_AllOA.txt --out AllOA.ldsc --merge-alleles w_hm3.snplist
#Repeat for all phenos

#Running the analyses
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
$ldsc_soft/ldsc.py --rg AllOA.sumstats.gz,SpineOA.ldsc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out AllOA-SpineOA_correlation
#Repeat for all pairs of OA phenos

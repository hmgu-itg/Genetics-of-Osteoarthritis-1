#!/bin/bash

#https://github.com/bulik/ldsc
#https://data.broadinstitute.org/alkesgroup/LDSCORE/
#https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
#https://github.com/bulik/ldsc/wiki/FAQ

#Path to software
/nfs/team144/software/ldsc/

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2

#Prepare the files
zcat HD_hip.allchr.snptest.results.withkey.passedSNPQC_overallMAF.HWE.INFO.tables.allchr_headers |  awk 'OFS="\t"{print $2,$3,$4,$5,$6,$10,$14,$15,$16,$19}' >  forldsc_UKBB_HD_hip.txt
PATH=/nfs/team144/software/anaconda/bin:$PATH
/nfs/team144/software/ldsc/munge_sumstats.py --sumstats forldsc_UKBB_HD_hip.txt --out UKBB_HD_hip.ldsc --merge-alleles w_hm3.snplist

#Running the analyses

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2

/nfs/team144/software/ldsc/ldsc.py --rg DDH.ldsc.sumstats.gz,UKBB_HD_hip.ldsc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out DDH-UKBB-HD-hip_correlation

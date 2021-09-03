#!/bin/bash
#Path to PRSice
prs=$1
#Path to base file
base=$2
#Path to files used for ld calculation
ld=$3
#Path to target phenotype files
target=$4

for PHENO in HipOA AllOA FingerOA HandOA SpineOA KneeOA KneeHipOA  THR ThumbOA TJR TKR ; do
echo "$prs \
--base $base/GO.FILTER.GW.${PHENO}.noUKBB.13082019.txt.gz \
--snp CPTID \
--chr CHR  \
--A1 EA  \
--A2 NEA  \
--or  \
--stat OR  \
--pvalue P  \
--target $ld/chr#-noduplicates-missnp,$target/${PHENO}.fam \
--binary-target T \
--type bed \
--cov $target/Covar.txt \
--cov-col @PC[1-20] \
--clump-kb 1M \
--clump-r2 0.1 \
--clump-p 1e-04 \
--ld $ld/chr#-noduplicates-missnp,$target/${PHENO}.fam \
--ld-keep $ld/LD.fam \
--lower 1.3e-08  \
--upper 1e-4  \
--interval 1e-07 \
--logit-perm \
--perm 10000 \
--all-score \
--print-snp \
--keep-ambig \
--thread 20 \
--out ${PHENO}_thread20_Nopinterval_logit-perm10k-keep-ambig" > ${PHENO}_thread20_Nopinterval_logit-perm10k-keep-ambig.cmds
echo "bsub -G t144_nargwas -q basement -n20 -R 'span[hosts=1]' -M20000 -R 'select[mem>20000] rusage[mem=20000]' -o ${PHENO}_thread20_Nopinterval_logit-perm10k-keep-ambig.o -e ${PHENO}_thread20_Nopinterval_logit-perm10k-keep-ambig.e ./${PHENO}_thread20_Nopinterval_logit-perm10k-keep-ambig.cmds" >> bsub.thread20_Nopinterval_logit-perm10k-keep-ambig.cmds ; done
chmod 770 *.cmds
./bsub.thread20_Nopinterval_logit-perm10k-keep-ambig.cmds

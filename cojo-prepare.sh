#!/bin/bash

#Arguments
trait=$1 
phenotype=$2 # Names to be used: AllOA, HipOA, KneeOA, KneeHipOA, THR, TKR, TJR, SpineOA, HandOA, FingerOA, ThumbOA 
chr=$3
data=$4 # path to data
plink_files=$5

mkdir $phenotype
cd $phenotype

#DATA PATHS
IN=$data/$phenotype

echo "Clupming..."

awk '{print $1,$2,$NF}' ${trait}.sample| tail -n+3 > ${trait}_pheno.txt

plink \
	--bed $plink_files/chr${chr}-noduplicates-missnp.bed \
	--bim $plink_files/chr${chr}-noduplicates-missnp.bim \
	--fam $plink_files/chr${chr}-noduplicates-missnp.fam \
	--pheno ${trait}_pheno.txt \
	--1 \
	--clump-field P \
	--clump-snp-field CPTID \
	--clump ${IN}/GO.RAW.final.meta.results.$phenotype.txt \
	--clump-p1 1.3e-8 --clump-p2 5e-2 --clump-r2 0.1 --clump-kb 1000 \
	--clump-verbose \
	--out clumping.$phenotype.chr${chr} 
	


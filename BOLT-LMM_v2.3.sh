#!/bin/bash

# PROGRAM PATHS
bolt=$1
# DATA PATHS
Geno_bed=$2
Geno_bim=$3
Geno_fam=$4
pheno_cov=$5
sample=$6

$bolt \
        --bed=$Geno_bed \
        --bim=$Geno_bim \
        --fam=$Geno_fam \
        --remove=bolt.in_plink_but_not_imputed.FID_IID.968.txt \
        --exclude=autosome_maf_lt_0.001.txt \
	--exclude=autosome_missing_gt_0.1.txt \
        --phenoFile=$pheno_cov \
        --phenoCol=OA_phenotype \
        --covarFile=$pheno_cov \
        --qCovarCol=PC{1:20} \
        --LDscoresFile=LDSCORE.1000G_EUR.tab.gz \
        --LDscoresMatchBp \
        --geneticMapFile=genetic_map_hg19_withX.txt.gz \
        --lmmForceNonInf \
        --numThreads=8 \
        --statsFile=chr${LSB_JOBINDEX}_bolt_HD_hipknee-OA.stats.gz \
        --bgenFile=ukb_imp_chr${LSB_JOBINDEX}_v3.bgen \
        --bgenMinINFO=0.3 \
        --sampleFile=$sample \
        --statsFileBgenSnps=chr${LSB_JOBINDEX}_bolt_HD_hipknee-OA.bgen.stats.gz \
        --verboseStats

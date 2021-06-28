#!/bin/bash

/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/BOLT-LMM_v2.3/bolt \
        --bed=/lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_cal_chr{1:22}_v2.bed \
        --bim=/lustre/scratch115/projects/ukbiobank/FullRelease/Genotypes/001/ukb_snp_chr{1:22}_v2.bim \
        --fam=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/HD_hipknee/ukb9979_cal_chr1_v2_s488363_HD_hipknee_FixCol6.fam \
        --remove=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/SR/bolt.in_plink_but_not_imputed.FID_IID.968.txt \
        --exclude=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/SR/autosome_maf_lt_0.001.txt \
		--exclude=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/SR/autosome_missing_gt_0.1.txt \
        --phenoFile=/lustre/scratch115/projects/ukbiobank_oa/OA_sample-file_v3_12042018/HD_hipknee_ukb9979_imp_chr1_v3_s487395_pheno-covar.sample \
        --phenoCol=OA_phenotype \
        --covarFile=/lustre/scratch115/projects/ukbiobank_oa/OA_sample-file_v3_12042018/HD_hipknee_ukb9979_imp_chr1_v3_s487395_pheno-covar.sample \
        --qCovarCol=PC{1:20} \
        --LDscoresFile=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/BOLT-LMM_v2.3/tables/LDSCORE.1000G_EUR.tab.gz \
        --LDscoresMatchBp \
        --geneticMapFile=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/BOLT-LMM_v2.3/tables/genetic_map_hg19_withX.txt.gz \
        --lmmForceNonInf \
        --numThreads=8 \
        --statsFile=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/HD_hipknee/chr${LSB_JOBINDEX}_bolt_HD_hipknee-OA.stats.gz \
        --bgenFile=/lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/EGAD00010001474/ukb_imp_chr${LSB_JOBINDEX}_v3.bgen \
        --bgenMinINFO=0.3 \
        --sampleFile=/lustre/scratch115/projects/ukbiobank_oa/OA_sample-file_v3_12042018/HD_hipknee_ukb9979_imp_chr1_v3_s487395.sample \
        --statsFileBgenSnps=/lustre/scratch115/projects/ukbiobank/Kostas/FULL_RELEASE/BOLT-LMM_v2.3/VERSION_3/HD_hipknee/chr${LSB_JOBINDEX}_bolt_HD_hipknee-OA.bgen.stats.gz \
        --verboseStats

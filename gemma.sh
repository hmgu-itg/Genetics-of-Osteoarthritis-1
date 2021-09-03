#!/bin/bash
#PROGRAM PATHS
#plink=" "
tabix=" "
gemma=" "
# DATA PATHS
input_bimbamfiles=" "
matrix_bimbamfile=" "

#Usage
#./gemma.bash <phenotype> <number of column in the phenotype file> <location of the phenotype file> <location of the annotation file> <name of the cohort> <location of the covariates file>

#Arguments
phenotype=$1
n=$2
phenotype_file=$3
annotation_file=$4
cohort=$5
cov=$6

# Exports
export CF9_R_LIBS=""
export PATH=" ":$PATH
export R_LIBS=" "

mkdir $phenotype
cd $phenotype

# Associate
# Here just generate the GEMMA jobs!
echo "Associating..."
for i in `ls $input_bimbamfiles | sed 's/.*\///'`
do
    echo $gemma -g ${input_bimbamfiles}/$i -p $phenotype_file -n $n -a $annotation_file -notsnp -maf 0 -miss 1 -km 1 -k $matrix_bimbamfile -lmm 4 -c $cov -o $i
done | ./array 30g asc | sed 's/red/green/'> assoc.command
chmod +x assoc.command
jobid=$(./assoc.command | sed 's/Job <//;s/> is.*//')

echo "Watching for association job array $jobid to finish..."
sleep 5
njobs=$(bjobs | grep -w $jobid | wc -l)

while [ $njobs -gt 0  ]
do
    njobs=$(bjobs | grep -w $jobid | wc -l)
    sleep 5
done


xitstatus=$(grep xited asc*.o )
xited=$(grep xited asc*.o | wc -l)
if [ $xited -gt 0 ]
then
    echo "Some association jobs have failed:"
    echo $xitstatus
    exit
fi
rm asc*[eo]

echo "Concatenating..."

head -n1 $(ls output/*.assoc.txt | head -n1) > $cohort.$phenotype.assoc.txt
for i in {1..22}
do
    for j in `ls output/*chr${i}_*.assoc.txt`
    do
	tail -n+2 $j
    done | sort -k3,3n
done >> $cohort.$phenotype.assoc.txt


echo "Building graphs.."

bgzip $cohort.$phenotype.assoc.txt
$tabix -s 1 -b 3 -e 3 -S 1 $cohort.$phenotype.assoc.txt.gz
gsub 10g -I -q yesterday ./man_qq_annotate --chr-col chr --pos-col ps --auto-label --pval-col p_score --title "${cohort}-$phenotype" --sig-thresh 1e-08 --sig-thresh-line 1e-08 $cohort.$phenotype.assoc.txt.gz $cohort.$phenotype.assoc

maf=0.05
zcat $cohort.$phenotype.assoc.txt.gz | awk '$7>'$maf  | bgzip > $cohort.$phenotype.maf$maf.assoc.txt.gz
$tabix -s 1 -b 3 -e 3 -S 1 $cohort.$phenotype.maf$maf.assoc.txt.gz
gsub 10g -I -q yesterday ./man_qq_annotate --chr-col chr --pos-col ps --auto-label --pval-col p_score --title "${cohort}-$phenotype-MAF$maf" --sig-thresh 5e-08 --sig-thresh-line 5e-08 $cohort.$phenotype.maf$maf.assoc.txt.gz $cohort.$phenotype.maf$maf.assoc

maf=0.01
zcat $cohort.$phenotype.assoc.txt.gz | awk '$7>'$maf  | bgzip > $cohort.$phenotype.maf$maf.assoc.txt.gz
$tabix -s 1 -b 3 -e 3 -S 1 $cohort.$phenotype.maf$maf.assoc.txt.gz
gsub 10g -I -q yesterday ./man_qq_annotate --chr-col chr --pos-col ps --auto-label --pval-col p_score --title "${cohort}-$phenotype-MAF$maf" --sig-thresh 5e-08 --sig-thresh-line 5e-08 $cohort.$phenotype.maf$maf.assoc.txt.gz $cohort.$phenotype.maf$maf.assoc


cd ..

echo "Finished"

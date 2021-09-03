# GWAMA.sh
# 30th March 2020

# Path to working directory
wd=$1

# Data paths
input=$2

# Install latest version of GWAMA:
mkdir $wd/GWAMA
mkdir $wd/GWAMA/gwama
cd $wd/GWAMA/gwama
wget http://www.geenivaramu.ee/tools/GWAMA_v2.2.2.zip
unzip GWAMA_v2.2.2.zip
make
# Excecutable is called: $wd/GWAMA/gwama/GWAMA 

# Create a list of the files for each phenotype to input into GWAMA
for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do
echo "$input/Female_${PHENO}/GO.FILTER.GW.final.meta.results.Female_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz F" > ${PHENO}.gwamafile.in
echo "$input/Male_${PHENO}/GO.FILTER.GW.final.meta.results.Male_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz M" >> ${PHENO}.gwamafile.in
done

# Notes for GWAMA
# https://genomics.ut.ee/en/tools/gwama-tutorial

# Input requirements are:
# 1) MARKERNAME – snp name
# 2) EA – effect allele
# 3) NEA – non effect allele
# 4) OR - odds ratio
# 5) OR_95L - lower confidence interval of OR
# 6) OR_95U - upper confidence interval of OR
# In case of quantitative trait:
# 4) BETA – beta
# 5) SE – std. error
# Study files might also contain columns:
# 7) N - number of samples
# 8) EAF – effect allele frequency
# 9) STRAND – marker strand (if the column is missing then program expects all markers being on positive strand)
# 10) IMPUTED – if marker is imputed or not (if the column is missing then all markers are counted as directly genotyped ones)

# Print the commands into a single file to run the analysis:
for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do
echo "$wd/GWAMA/gwama/GWAMA --filelist ${PHENO}.gwamafile.in --sex --name_marker CPTID --name_or OR --name_or_95l OR_L95 --name_or_95u OR_U95 --output ${PHENO}.gwama.out" >>gwama.cmds
done

chmod 770 gwama.cmds
bsub -G t144_nargwas -M20000 -R 'select[mem>20000] rusage[mem=20000]'  -o gwama.o -e gwama.e ./gwama.cmds


############################################################################
##   Tabix the files to produce regional plots    						  ##
############################################################################

#!/bin/sh
# Path to hgi profile
hgi=$3
./$hgi
# Path to team144 software
team=$4
export PATH=$team/bin:$PATH
module add $(module avail 2>&1 | grep '/tabix/' | grep latest | sed 's/.latest.//')


#### To prepare the bgzipped files
#Path to tabix directory
tabix=$5
#Females
for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do cat <(echo \#"$(zcat $input/Female_${PHENO}/GO.FILTER.GW.final.meta.results.Female_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz | head -1)") <(zcat $input/Female_${PHENO}/GO.FILTER.GW.final.meta.results.Female_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz | tail -n+2 | sort -k20,20n -k21,21n) | bgzip > $tabix/GO-meta1.Female_${PHENO}.txt.bgz & #; done
done
#Males
for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do cat <(echo \#"$(zcat $input/Male_${PHENO}/GO.FILTER.GW.final.meta.results.Male_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz | head -1)") <(zcat $input/Male_${PHENO}/GO.FILTER.GW.final.meta.results.Male_${PHENO}.p1e-06.n2.MAF0.0001.txt.gz | tail -n+2 | sort -k20,20n -k21,21n) | bgzip > $tabix/GO-meta1.Male_${PHENO}.txt.bgz & #; done
done
####### Then to tabix:
for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do tabix -c "#" -s 20 -b 21 -e 21 $tabix/GO-meta1.Female_${PHENO}.txt.bgz &
done

for PHENO in AllOA HipOA KneeOA KneeHipOA HandOA ThumbOA FingerOA SpineOA THR TKR TJR ; do tabix -c "#" -s 20 -b 21 -e 21 $tabix/GO-meta1.Male_${PHENO}.txt.bgz &
done

#Produce some regional plots for ALDH1A2
#Path to regional plots output directory
rp=$6
mkdir $rp
cd $rp

# Prepare a list of the regions to plot:

# prepare the input files for locus zoom:
# (1) GWAS summary stats
mkdir $rg/input
cat ALHD1A2.list.txt  | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo " tabix -h $tabix/GO-meta1.Female_${PHENO}.txt.bgz ${Window1Mb} | cut -f1,8-10,20,21 | sed '/^#/ d' | awk 'BEGIN{OFS=\"\\t\"; print \"CPTID\",\"BETA\",\"SE\",\"P\",\"chromosome\",\"position\",\"rsid\"}{if(\$1==\$1){newid=\"chr\"\$5\":\"\$6} print \$0,newid;}' > input/${PHENO}_Female_${CPTID}_1Mbwindow.input"  >> Prepare.input.cmds ; done
cat ALHD1A2.list.txt  | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo " tabix -h $tabix/GO-meta1.Female_${PHENO}.txt.bgz ${Window2Mb} | cut -f1,8-10,20,21 | sed '/^#/ d' | awk 'BEGIN{OFS=\"\\t\"; print \"CPTID\",\"BETA\",\"SE\",\"P\",\"chromosome\",\"position\",\"rsid\"}{if(\$1==\$1){newid=\"chr\"\$5\":\"\$6} print \$0,newid;}' > input/${PHENO}_Female_${CPTID}_2Mbwindow.input"  >> Prepare.input.cmds ; done
cat ALHD1A2.list.txt  | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo " tabix -h $tabix/GO-meta1.Male_${PHENO}.txt.bgz ${Window1Mb} | cut -f1,8-10,20,21 | sed '/^#/ d' | awk 'BEGIN{OFS=\"\\t\"; print \"CPTID\",\"BETA\",\"SE\",\"P\",\"chromosome\",\"position\",\"rsid\"}{if(\$1==\$1){newid=\"chr\"\$5\":\"\$6} print \$0,newid;}' > input/${PHENO}_Male_${CPTID}_1Mbwindow.input"  >> Prepare.input.cmds ; done
cat ALHD1A2.list.txt  | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo " tabix -h $tabix/GO-meta1.Male_${PHENO}.txt.bgz ${Window2Mb} | cut -f1,8-10,20,21 | sed '/^#/ d' | awk 'BEGIN{OFS=\"\\t\"; print \"CPTID\",\"BETA\",\"SE\",\"P\",\"chromosome\",\"position\",\"rsid\"}{if(\$1==\$1){newid=\"chr\"\$5\":\"\$6} print \$0,newid;}' > input/${PHENO}_Male_${CPTID}_2Mbwindow.input"  >> Prepare.input.cmds ; done 

./$hgi
export PATH=$team/bin:$PATH
module add $(module avail 2>&1 | grep '/tabix/' | grep latest | sed 's/.latest.//')
chmod 770 Prepare.input.cmds
bsub -G t144_nargwas -o Prepare.input.cmds.o -e Prepare.input.cmds.e ./Prepare.input.cmds

# (2) The input SNPlist for the LD files
cd $rp/input
for FILE in `ls *1Mbwindow.input` ; do cat ${FILE} | awk 'BEGIN{OFS="\t"}{if($1==$1) print $7,$5,$6;}' | sed 's/rsid/snp/g'| sed 's/chromosome/chr/g' | sed 's/position/pos/g' > ${FILE}.snp_pos.txt ; done
for FILE in `ls *2Mbwindow.input` ; do cat ${FILE} | awk 'BEGIN{OFS="\t"}{if($1==$1) print $7,$5,$6;}' | sed 's/rsid/snp/g'| sed 's/chromosome/chr/g' | sed 's/position/pos/g' > ${FILE}.snp_pos.txt ; done


# Prepare running commands for the Locus Zoom script:
# Files to use for LD calculation
ld=$7
cd $rp/input
cat ../ALHD1A2.list.txt | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo "./locusZoomForGO_1MbWindow-Loz.sh ${PHENO}_Female_${CPTID}_1Mbwindow.input ${rsid} $ld/chr${CHR}-noduplicates " >> GO_ALDH1A2.1MbWindow.LZ.cmds ; done 
cat ../ALHD1A2.list.txt | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo "./locusZoomForGO_2MbWindow-Loz.sh ${PHENO}_Female_${CPTID}_2Mbwindow.input ${rsid} $ld/chr${CHR}-noduplicates " >> GO_ALDH1A2.2MbWindow.LZ.cmds ; done 
cat ../ALHD1A2.list.txt | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo "./locusZoomForGO_1MbWindow-Loz.sh ${PHENO}_Male_${CPTID}_1Mbwindow.input ${rsid} $ld/chr${CHR}-noduplicates " >> GO_ALDH1A2.1MbWindow.LZ.cmds ; done 
cat ../ALHD1A2.list.txt | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo "./locusZoomForGO_2MbWindow-Loz.sh ${PHENO}_Male_${CPTID}_2Mbwindow.input ${rsid} $ld/chr${CHR}-noduplicates " >> GO_ALDH1A2.2MbWindow.LZ.cmds ; done 


chmod 770 *cmds
./GO_ALDH1A2.1MbWindow.LZ.cmds
./GO_ALDH1A2.2MbWindow.LZ.cmds

#########################################################################################
# Extract the information for each of the 100 GO signals for the appropriate phenotype
#########################################################################################
# Working dierctory
GO_signals=$8
cd $GO_signals

cat Signals.txt | tail -n+2 | while read PHENO CPTID EA NEA SNP CHR POS rsid Window1Mb Window2Mb ; do echo "zgrep  ${CPTID} ${PHENO}.gwama.out.out " >> GWAMAsexResultsfor100GOindepSignals.cmds ; done 
chmod 770 *cmds
./GWAMAsexResultsfor100GOindepSignals.cmds >> GWAMAsexResultsfor100GOindepSignals.results

# Plot a qq for the Pdiff, Phet and a Manhat for th meta-analysis of the sex specific results only
cd $GO_signals
# Column header is:
head -n 1 TKR.gwama.out.out | sed 's/\t/\n/g' | awk '{print NR, $0}'

# get a list of the files:
ls *out.out | sed 's/.out.out//g' > GWAMA.filelist.txt 
more GWAMA.filelist.txt 

# Add in chromosome and position from the first colum (3:52236363_C_T)
for FILE in `cat GWAMA.filelist.txt` ; do cat ${FILE}.out.out | awk '{if(NR >1){split($1,a,"_"); chrpos=a[1]; print chrpos,$0;}}' | awk 'BEGIN{OFS="\t"; print "chrpos","rs_number","reference_allele","other_allele","eaf","OR","OR_se","OR_95L","OR_95U","z","p-value","_-log10_p-value","q_statistic","q_p-value","i2","n_studies","n_samples","effects","male_eaf","male_OR","male_OR_se","male_OR_95L","male_OR_95U","male_z","male_p-value","male_n_studies","male_n_samples","female_eaf","female_OR","female_OR_se","female_OR_95L","female_OR_95U","female_z","female_p-value","female_n_studies","female_n_samples","gender_differentiated_p-value","gender_heterogeneity_p-value","CHR","POS"}{if(NR >=1){split($1,a,":"); chr=a[1]; pos=a[2]; print chrpos,$0,chr,pos;}}' | gzip > ${FILE}.plotting.gz ; done
for file in `ls *.plotting.gz ` ; do zcat ${file} | wc -l ; done

ls *.plotting.gz

for file in `cat GWAMA.filelist.txt` ; do wc -l ${file}.out.out ; done

grep NA AllOA.gwama.out.out #<-no nothing printed to screen!
grep "\-9" AllOA.gwama.out.out > temp
head temp

ls *.plotting.gz > GWAMA.plotfilelist.txt
for FILE in `cat GWAMA.plotfilelist.txt` ; do zgrep -v "\-9" ${FILE} | gzip > ${FILE}.minus9sout.gz ; done
ls *.minus9sout.gz > GWAMA.plotfilelist.txt
for FILE in `cat GWAMA.plotfilelist.txt` ; do zcat ${FILE} | wc -l  ; done

# Running on gen1 interactive
#Path to manqq
run_manqq=$9
#5e-8
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col gender_differentiated_p-value  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.gender_differentiated.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col gender_heterogeneity_p-value  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.gender_heterogeneity.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col p-value  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.metanalysis.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col male_p-value  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.Male-metanalysis.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col female_p-value  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.Female-metanalysis.plot ; done
#1.3e-8
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col gender_differentiated_p-value  --sig 1.3e-8 --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.gender_differentiated_pval-1-3e8.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col gender_heterogeneity_p-value --sig 1.3e-8  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.gender_heterogeneity_pval-1-3e8.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col p-value  --sig 1.3e-8  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.metanalysis_pval-1-3e8.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col male_p-value --sig 1.3e-8  --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.Male-metanalysis_pval-1-3e8.plot ; done
cat GWAMA.plotfilelist.txt | while read FILE ; do ./Rscript $run_manqq/run_manqq.R --chr-col CHR --pos-col POS --pval-col female_p-value  --sig 1.3e-8 --image pdf --a2 reference_allele --a1 other_allele --build 37 --af-col eaf $FILE $FILE.Female-metanalysis_pval-1-3e8.plot ; done

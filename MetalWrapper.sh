#########################################
##       METALWrapper.README         ###
#########################################
### Creator: Cindy G. Boer
### Last Updated: 27/07/2018
### For questions: c.boer@erasmumc.nl
###########################################
# Step 3 in GO-meta-analysis: METAL meta-analysis
#
###########################################
### Dependancies:
# METAL
# CLEAND.easyQC output files
#
### RUN script with:
# ./MetalWrapper.sh [Prefix.name] [PATH/DIR/Cleaned.EASYQC.files]
# 
# Run METAL: [Will do automatically]
# bsub -G t144_oagwas_meta -J METAL.[prefix] -o [prefix].METAL.LOG -e [Prefix].METAL.ERR -R "select[mem>2000] rusage[mem=2000]" -M2000 < [Prefix].meta.METAL.cmd

# script generates metal.pa and metal.cmd file for meta-analysis0
# Meta analysis will sart with all the files in the PATH/DIR/ given
#
#!/usr/local/bin/bash
#########################################################
#### Genetics of Osteoarthritis MetalWrapper        #####
#########################################################
### GO consortium
### By Cindy G. Boer Erasmus MC
### c.boer@erasmusmc.nl
### Version 1.0
### Last-updated: 26/7/2018 
##########################################
## INPUT
# $1 = Prefix Name Meta-analysis
# $2 = PATH/DIR files for meta-analysis
# $3 = OUTPUT DIR
########################################## 
cd $3
touch $1.meta.METAL.par

## Make METAL.par file
echo " 
#########################################
###  Genetics of Osteoarthritis METAL ###
#########################################
CLEAR
# Input columns:
MARKER CPTID
ALLELE EA NEA
EFFECT BETA
STDERR SE
FREQLABEL EAF
PVALUE P
WEIGHT N

# Custom Variables
CUSTOMVARIABLE NCASES
LABEL NCASES AS NCASES
CUSTOMVARIABLE NCONTROLS
LABEL NCONTROLS AS NCONTROLS
CUSTOMVARIABLE N
LABEL N AS N
  
#Input file seperators
SEPARATOR TAB
 
# Metal Options:
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
GENOMICCONTROL ON"  >> $1.meta.METAL.par

for i in $(ls -d -1 $2/**)
do
echo PROCESS $i >> $1.meta.METAL.par
done

echo "
OUTFILE MetaAnalysis.GO-meta1.$1 .TBL

ANALYZE HETEROGENEITY
QUIT" >> $1.meta.METAL.par

###########################################
### Make Metal .CMD
#Path to metal files
metal_files=$4
touch $1.meta.METAL.cmd
echo -e "#!/usr/local/bin/bash
$metal_files/metal $3$1.meta.METAL.par > $1.meta.METAL.log" > $1.meta.METAL.cmd

chmod +x $1.meta.METAL.cmd
chmod +x $1.meta.METAL.par

###########################################
### Run METAL

./gsub 19G -R"span[hosts=1]" -q yesterday -G t144_oagwas_meta  -o $1.METAL.o  -e $1.METAL.e $3$1.meta.METAL.cmd

###########################################
exit


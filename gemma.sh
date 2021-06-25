#!/usr/bin/bash
  
panel=$1

# protein

prot=$2

# chromosome 
     
chr=$3

# position of the SNP 

pos=$4

# file with the "UniprotID rsID" pairs from the external study 

extfile=$5

# output directory for auxillary files; GEMMA saves its output in the separate "output" directory

outdir=$6

# logfile 

logfile=$7

if [ ! -d $outdir ];then
    mkdir $outdir
fi

# prefix for auxillary files 
prefix="$outdir/$panel"."$prot"."$chr"."$pos"

# prefix for GEMMA output files
outfile="$panel"."$prot"."$chr"."$pos".out

echo "INPUT : PANEL: $panel" >> ${logfile}
echo "INPUT : PROT: $prot" >> ${logfile}
echo "INPUT : CHR: $chr" >> ${logfile}
echo "INPUT : POS: $pos" >> ${logfile}
echo "INPUT : EXTFILE: $extfile" >> ${logfile}
echo "INPUT : OUTDIR: $outdir" >> ${logfile}
echo "PREFIX: $prefix" >> ${logfile}

# file with genotypes
genofile="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/single_point/input/whole_genome/autosomal.correctsnames"

# GEMMA executable
GEMMA="/storage/hmgu/software/gemma-0.94"

# relatedness matrix
RELMAT="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/relmat/final/single-point/merge.cXX.txt"

# script to get most recent (chr pos) for a given rsID
#script="/home/andrei/rsid.to.pos.b38.sh"
script="/home/andrei/scripts/getSNPb38Position.py"

# uniprot file which gives an UniprotID to a protein
ufile="/storage/sanger/team144/OLINK/ALL_ensembl_uniprot"

# uniprot ID for the input
#uID=$(fgrep -w $prot $ufile | cut -f 5)
uID=$(awk -v p=$prot 'BEGIN{FS="\t";}$2==p{print $5;exit 0;}' $ufile)

echo "current UniprotID : " $uID >> ${logfile}

# results of single-point association are in here
assocfile="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/OLINK/single_point/"$panel"/"$panel"."$prot"/MANOLIS."$panel"."$prot".assoc.txt.gz"

# contains only the current SNP
totestIDs="$prefix".IDs.totest
echo "chr"${chr}:${pos}"[b38]" > $totestIDs

# SNPs to condition upon: they correspond to the SNPs (from the external study) that regulate the same protein as the input SNP
tocondIDs="$prefix".IDs.tocond
if [ -f $tocondIDs ];then
        rm $tocondIDs
fi
touch $tocondIDs

# loop over all variants from the external study that correspond to the given uniprot ID, convert their rsID to b38 coordinates
# and if they are on the same chromosome and within 1Mbp window around the input SNP, add their single-point association results to the $tocondIDs

wk -v id=$uID '$1==id{print $2;}' $extfile | while read rs;
do
    read -r cc pp <<<$($script $rs 2>>$logfile)
    echo "SNP from the external study: $rs ($cc $pp) Input SNP: ($chr $pos)" >> ${logfile}
    #read -r cc pp <<<$($script $rs|sed 's/:/ /'| sed 's/^chr//')
    if [[ ${cc} == ${chr} ]];then # if they are on the same chromosome
        let "d=$pp-$pos"
        #echo $cc $chr $pp $pos $d
        if [[ $d -lt 1000000 && $d -gt -1000000 ]];then # and if they are within 1Mbp window from each other
            echo "chr"${cc}:${pp}"[b38]" >> $tocondIDs
            echo "$cc==$chr; |$pp-$pos|<1000000; output $rs into $tocondIDs" >> ${logfile}
        fi
    fi
done

# file with phenotypes
phenofile="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/phenotypes/OLINK/MANOLIS"."$panel"."$prot"."txt"

totestfile="$prefix".totest
tocondfile="$prefix".tocond
covarfile="$prefix".covar

# is the input SNP one the conditioning SNPs ? If yes, we don't have to test anything
if [[ $(fgrep "${chr}:${pos}[b38]" $tocondIDs | wc -l) -gt 0 ]]; then
    echo "MORE THAN 1 occurrence of ${chr}:${pos} : same variant" >> ${logfile}
    echo >> ${logfile}
    touch "$prefix".samevariant
else
    if [[ $(cat $tocondIDs | wc -l) -lt 1 ]]; then
        echo "No variants of condition study present in our dataset for $panel $prot" >> ${logfile}
        echo >> ${logfile}
    else
        # preparing input files for GEMMA
        echo >> ${logfile}
        echo "PLINK: selecting the current SNP" >> ${logfile}
        plink --bfile "$genofile" --extract "$totestIDs" --make-bed --out "$totestfile" --pheno <(awk '{print $1,$1, $3}' $phenofile) --allow-no-sex 2>&1 >>$logfile
        echo >> ${logfile}
        echo "PLINK: selecting the conditional SNP(s)" >> ${logfile}
        plink --bfile "$genofile" --extract "$tocondIDs" --recode A --out "$tocondfile" 2>&1 >>$logfile
        echo >> ${logfile}

        if [[ -f "$tocondfile".raw ]];then # sometimes the "conditioning" variant is not in our dataset, in which case the *.raw file doesn't exist
            tail -n +2 "$tocondfile".raw | cut -d ' ' -f 7- |awk '{print(1, $0)}' > "$covarfile"

            # calling GEMMA
            echo "GEMMA:" >> ${logfile}
            "$GEMMA" -bfile "$totestfile" -n 1 -notsnp  -maf 0  -miss 1  -km 1 -k "$RELMAT" -lmm 4 -c "$covarfile" -o "$outfile" 2>&1 >>$logfile
        else
            echo "No variants of condition study present in our dataset for $panel $prot" >> ${logfile}
            echo >> ${logfile}
        fi
    fi
fi

echo "=============================================================================" >> ${logfile}


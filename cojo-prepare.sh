#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -i <input table>"
    echo "          -o <output table>"
    echo "          -k <known signals>"
    echo "          -w <bp window; default: 1000000>"
    exit 0
}

# default
window=1000000

OPTIND=1
while getopts "i:o:k:w:" optname; do
    case "$optname" in
        "i" ) input="${OPTARG}";;
        "o" ) output="${OPTARG}";;
        "k" ) known="${OPTARG}";;
        "w" ) window="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

winstr=$(echo $window|perl -lne 'print sprintf("%0.2fMb",$_/1000000);')
logfile="$PWD"/"cojo-prepare_$winstr"."log"

# path to meta-analysis results
ma_path="/storage/hmgu/projects/helic/OLINK/meta_analysis"

# PLINK files
ha_plink="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/single_point/input/whole_genome/autosomal.correctsnames"
hp_plink="/storage/sanger/projects/helic/t144_helic_15x/analysis/HP/missingness/gw/merged.shortnames"

# pheno files prefixes
ha_pheno="/storage/sanger/projects/helic/t144_helic_15x/analysis/HA/phenotypes/OLINK/MANOLIS"
hp_pheno="/storage/hmgu/projects/helic/OLINK/HP/phenotypes"

echo "window        : $window" > "$logfile"
echo "known signals : $known" >> "$logfile"
echo "input table   : $input" >> "$logfile"
echo "output table  : $output" >> "$logfile"
echo >> "$logfile"

> "$output"

tmpfile=$(mktemp "$PWD"/tmp_cojo.XXXX)
tmpfile2=$(mktemp "$PWD"/tmp_cojo2.XXXX)
echo "Temp file   : $tmpfile" >> "$logfile"
echo "Temp file 2 : $tmpfile2" >> "$logfile"
echo >> "$logfile"

# reading the input table
cut -f 1,2,5,6,8,10-12,16-18,24 "$input"| tail -n +2|tr '\t' ' '|while read panel prot uniprot chr pos a1 a2 f1 b se p nMiss; do
prefix="$panel"."$prot"
id="$chr:$pos"
suffix="$chr"_"$pos"
varid="$panel"_"$prot"_"$suffix"

echo "=========================================" >> "$logfile"
echo "Variant: $prefix $id" >> "$logfile"
echo >> "$logfile"

#------------------- total number of samples ------------------------------
phenofile="$ha_pheno"."$prefix"."txt"
if [[ ! -f "$phenofile" ]];then
    echo "ERROR: $phenofile doesn't exist" >> "$logfile"
    continue
fi
ha_samples=$(cat <(cut -f 1 "$phenofile") <(cut -f 2 -d ' ' "$ha_plink"."fam") | sort|uniq -d| tr '\n' ',')
ha_N=$(cat <(cut -f 1 "$phenofile") <(cut -f 2 -d ' ' "$ha_plink"."fam") | sort|uniq -d|wc -l)

phenofile="$hp_pheno"/"$panel"/"POMAK"."$panel"."$prot"."txt"
if [[ ! -f "$phenofile" ]];then
    echo "ERROR: $phenofile doesn't exist" >> "$logfile"
    continue
fi
hp_samples=$(cat <(cut -f 1 "$phenofile") <(cut -f 2 -d ' ' "$hp_plink"."fam") | sort|uniq -d| tr '\n' ',')
hp_N=$(cat <(cut -f 1 "$phenofile") <(cut -f 2 -d ' ' "$hp_plink"."fam") | sort|uniq -d|wc -l)
all_samples="$ha_samples""$hp_samples"
N=$((ha_N+hp_N-nMiss))
#--------------------------------------------------------------------------

echo "Total samples: $N" >> "$logfile"
echo >> "$logfile"

# sort if there are more than one ID in the uniprot string
uniprot=$(echo "$uniprot"| perl -lne '@a=split(/,/);print join(",",sort @a);')

start=$((pos-window))
if [[ "$start" -lt 0 ]];then
    start=0
fi
end=$((pos+window))
> "$tmpfile"
echo "intersectBed using chr=$chr start=$start end=$end" >> "$logfile"
intersectBed -wb -a <(echo "$chr $start $end"| tr ' ' '\t') -b "$known" | awk -v x="$uniprot" 'BEGIN{FS="\t";OFS="\t";}$8==x{print $0;}' | cut -f 1,3| tr '\t' ':' | sort | uniq > "$tmpfile" 2>>"$logfile"
echo >> "$logfile"

# if the tested variant is known
c=$(grep -c "$id" "$tmpfile")
if [[ "$c" -ne 0 ]];then
    echo "$id : is known" >> "$logfile"
    continue
fi

# if there are no known signals in the bp window
nKnown=$(cat "$tmpfile"| wc -l)
if [[ "$nKnown" -eq 0 ]];then
    echo "$id : no known signals" >> "$logfile"
    continue
fi

echo "$nKnown known signals found" >> "$logfile"
echo >> "$logfile"

# extract m/a results for known signals and output them
echo -n "Extracting m/a stats for known signals ... " >> "$logfile"
> "$tmpfile2"
cat "$tmpfile"| tr ':' ' '|while read cr ps;do
    tabix "$ma_path/$panel/METAL/$panel.$prot.metal.bgz" $cr:$ps-$ps| cut -f 1-5,9-11| awk -v id=$varid -v s=$all_samples -v c=$cr -v p=$ps -v n=$N 'BEGIN{FS="\t";OFS="\t";}$1==c && $2==p{print id,c,p,toupper($3),toupper($4),$5,$6,$7,$8,n,s;}'  >> "$tmpfile2"
done
echo "Done" >> "$logfile"
echo >> "$logfile"

# check if we have anything in the output
c=$(cat "$tmpfile2"| wc -l)
if [[ "$c" -eq 0 ]];then
    echo "$id : no m/a results for known signals" >> "$logfile"
    continue
fi

echo "$id: output m/a results for $c known signals" >> "$logfile"
cat "$tmpfile2" >> "$output"

# output details of the current variant
echo "$varid $chr $pos $a1 $a2 $f1 $b $se $p $N $all_samples" | tr ' ' '\t' >> "$output"
echo >> "$logfile"

done

rm "$tmpfile"
rm "$tmpfile2"

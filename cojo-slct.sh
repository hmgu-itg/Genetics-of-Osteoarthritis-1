#!/bin/bash

# a simple cojo-slct wrapper

function usage {
    echo ""
    echo "Usage: $0 -i <signal ID>"
    echo "          -n <number of samples in meta analysis>"
    echo "          -f <bfile prefix>"
    echo "          -t <optional: number of threads; default: 1>"
    echo "          -c <optional: collinearity threshold; default: 0.9>"
    echo "          -p <optional: p-value threshold; default: 5e-8>"
    exit 0
}

# default
threads=1
ct=0.9
pt=5e-8

OPTIND=1
while getopts "i:n:f:t:c:p:" optname; do
    case "$optname" in
        "i" ) id="${OPTARG}";;
        "n" ) N="${OPTARG}";;
        "f" ) bfile="${OPTARG}";;
        "t" ) threads="${OPTARG}";;
        "c" ) ct="${OPTARG}";;
        "p" ) pt="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

# m/a results
ma_results="/storage/hmgu/projects/helic/OLINK/meta_analysis"

logfile="cojo-slct.log"

echo "" >> "$logfile"
date >> "$logfile"
echo "" >> "$logfile"

echo "ID                 : $id" >> "$logfile"
echo "N                  : $N" >> "$logfile"
echo "bfile              : $bfile" >> "$logfile"
echo "collinearity       : $ct" >> "$logfile"
echo "p-value threshold  : $pt" >> "$logfile"
echo "threads            : $threads" >> "$logfile"

read -r panel prot chr pos <<<$(echo $id|tr '_' ' ')

gcta64 --bfile "$bfile" --cojo-slct --out "$id".slct --cojo-file <(zcat "$ma_results"/"$panel"/METAL/"$panel"."$prot".metal.bgz| cut -f 1-5,9-11| awk -v n=$N 'BEGIN{FS="\t";OFS="\t";}{if (NR==1){print "SNP","A1","A2","freq","b","se","p","N";}else{ print $1":"$2,$3,$4,$5,$6,$7,$8,n;}}') --extract "$id".cond --threads "$threads" --cojo-collinear "$ct" --cojo-p "$pt"

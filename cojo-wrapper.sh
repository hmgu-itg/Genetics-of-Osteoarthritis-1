#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -f <PLINK file prefix>"
    echo "          -o <output dir>"
    echo "          -i <input table>"
    echo "          -c <optional: collinearity threshold, default: 0.9>"
    exit 0
}

# default
ct=0.9

OPTIND=1
while getopts "f:o:i:c:" optname; do
    case "$optname" in
        "f" ) bfile="${OPTARG}";;
        "o" ) out="${OPTARG}";;
        "i" ) input="${OPTARG}";;
        "c" ) ct="${OPTARG}";;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

# required fields in the input table (tab separated)
# variant/protein association ID: panel_prot_chr_pos; chr and pos correspond to the variant being tested
#
# chr, pos, a1, a2, freq1, beta, SE, p, N: fields correspond to either the variant being tested, 
# or to conditioning variants, they have the same meaning as fields in a cojo file

out=${out/%/}

mkdir -p "$out"

logfile="$out"/"cojo-wrapper.log"
errfile="$out"/"cojo-wrapper.err"

date > "$logfile"
echo "working dir     : $PWD" >> "$logfile"
echo "bfile           : $bfile" >> "$logfile"
echo "input table     : $input" >> "$logfile"
echo "output dir      : $out" >> "$logfile"
echo "collinearity    : $ct" >> "$logfile"

# reading the input table
cut -f 1,11 "$input"|sort|uniq| while read varid samples; do

echo "================== $varid =========================" >> "$logfile"

id=$(echo $varid | cut -f 3- -d '_'| tr '_' ':')

# extract variants
plinkout="$out"/"$varid"
echo -n "Extracting $id and known signals ... " >> "$logfile"
plink --make-bed --bfile "$bfile" --out "$plinkout" --extract <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') --keep <(echo -n "$samples" | tr ',' '\n'| awk '{print $1,$1;}') --allow-no-sex
echo "Done " >> "$logfile"
echo >> "$logfile"

# check if there are enough variants
c=$(cat "$plinkout".bim | wc -l)
if [[ "$c" -lt 2 ]];then
    echo "$varid : not enough variants" >> "$logfile"
    echo >> "$logfile"
    continue
fi

# our variant should be in the PLINK output
c=$(cat "$plinkout".bim | grep "$id" | wc -l)
if [[ "$c" -eq 0 ]];then
    echo "$id : ERROR : not in bfile" >> "$logfile"
    echo >> "$logfile"
    continue
fi

# known signals which are not in the bfile
echo "$id: known signals not in the bfile: " >> "$logfile"
cut -f 2 "$plinkout"."bim" | cat - <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') | sort|uniq -u  >> "$logfile"
echo >> "$logfile"

# GCTA input files
cojofile="$plinkout"."ma"
condfile="$plinkout"."cond"
echo "SNP A1 A2 freq b se p N" | tr ' ' '\t' > "$cojofile"
fgrep -w "$varid" "$input"|cut -f 2-10|sed 's/\t/:/'| grep -v -f <(cut -f 2 "$plinkout"."bim" | cat - <(fgrep -w "$varid" "$input"|cut -f 2,3| tr '\t' ':') | sort|uniq -u) >> "$cojofile"
fgrep -v -w "$id" "$plinkout"."bim"| cut -f 2 > "$condfile"

# calling GCTA
echo -n "Calling GCTA ... " >> "$logfile"
gcta64 --bfile "$plinkout" --cojo-file "$cojofile" --cojo-cond "$condfile" --cojo-collinear "$ct" --out "$plinkout"."out" 1>> "$logfile" 2>> "$logfile"
echo "Done " >> "$logfile"
echo >> "$logfile"

done
date >> "$logfile"

# report signals for which GCTA failed
cat "$logfile" | perl -lne 'BEGIN{$cur=undef;$err=undef;%H=();}{if (/^=+\s+(\w.*\w)\s+=+$/){if (!defined($cur)){$cur=$1;}else{if (defined($err)){$H{$cur}=$err;} $cur=$1;$err=undef;}} if (/^Error:\s+(\S.*)$/){$err=$1;} }END{if (defined($cur) && defined($err)){$H{$cur}=$err;}foreach $x (keys %H){print $x."\t".$H{$x};}}' > "$errfile"


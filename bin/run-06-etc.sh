#!/bin/bash


cd "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross" ||
    {
        echo "Exiting: cd failed."
        exit 1
    }

# shellcheck disable=1091
. ./bin/auxiliary/functions-preprocessing.sh
. ./bin/auxiliary/functions-in-progress.sh

parallelize=4
indir="./data/files_bam"
outdir="./data/files_bam"
chunk="100000"

typeset -A bam=(
    ["${indir}/Disteche_sample_1.dedup.CAST.corrected.bam"]="CAST"
    ["${indir}/Disteche_sample_1.dedup.mm10.corrected.bam"]="mm10"
)
# for key value in "${(@kv)bam}"; do
#     echo "Key: ${key}"
#     echo "Value: ${value}"
# done

# parallel --header : -k -j "${parallelize}" echo \
parallel --header : -k -j "${parallelize}" \
"Rscript bin/workflow/06-get-AS-per-qname.R \
-b {bam} \
-i {bai} \
-o {outdir} \
-s {strain} \
-c {chunk}" \
::: bam "${(k)bam[@]}" \
:::+ bai "${(k)bam[@]}.bai" \
:::+ strain "${(v)bam[@]}" \
::: outdir "${outdir}" \
::: chunk "${chunk}"
# -b, --bam     bam infile, including path <chr>
# -i, --bai     bam index, including path <chr>
# -o, --outdir  directory for saving rds outfile, including path <chr>
# -s, --strain  strain name to be appended to rds outfile columns <chr>
# -c, --chunk   number of records to read into memory at a single time
#               <even int> [default: 100000]

parallel --header : -k -j "${parallelize}" \
""

typeset -A txt_in=(
    ["${indir}/Disteche_sample_1.dedup.CAST.corrected.bam"]="CAST"
    ["${indir}/Disteche_sample_1.dedup.mm10.corrected.bam"]="mm10"
)

# parallel --header : -k -j "${parallelize}" echo \
parallel --header : -k -j "${parallelize}" \
"zcat -df {txt_in} \
| cut -f 1 \
| sort \
| uniq -c \
| sort -nr \
| cut -c 7- \
| awk '$1 > 1 {print $0}' \
| tr ' ' \\t \
| gzip \
> {txt_out}" \
::: txt_in "${txt_in[@]}" \
:::+ txt_out "${txt_out[@]}" \


bash ./bin/workflow/03-remove-duplicate-qnames.sh \
-u FALSE \
-i "./data/files_bam/Disteche_sample_1.dedup.mm10.corrected.mm10.AS.txt.gz" \
-o "./data" \
-c TRUE \
-t TRUE \
-r FALSE

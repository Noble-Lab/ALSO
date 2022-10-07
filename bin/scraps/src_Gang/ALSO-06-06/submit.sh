#!/bin/bash

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }

module add java/1.8.0 
export PATH="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/jdk-18/bin/:$PATH"

module load picard/2.26.4
export PATH="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/subread-2.0.3-Linux-x86_64/bin/utilities/:$PATH"
module load samtools/1.14

module add pcre2/10.39
module add R/4.1.2

module add parallel/20200922 
module add gzip/1.12 

sample_id=$1
echo "${sample_id}"
#strain=$2 # "CAST-EiJ"
#echo "${strain}"
path_1="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO/input/mm10-preprocessed/"
path_2="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO/input/CAST-EiJ-preprocessed/"
bam_input="Disteche_sample_${sample_id}.dedup.corrected.corrected.bam" 
#input/mm10-preprocessed/Disteche_sample_1.dedup.corrected.corrected.bam
outpath="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO/results/"
#infile="Disteche_sample_${sample_id}.dedup.corrected.bam"

chmod 751 *.sh
chmod 751 bin/*.sh
chmod 751 bin/auxiliary/*.sh

#  Call from 2021_kga0_4dn-mouse-cross or a directory containing
#+ functions-in-progress.sh, functions-preprocessing-HPC.sh,
#+ get-AS-per-qname.R, find-set-intersection-set-complement.sh,
#+ generate-assignment-lists.R, and filter-qnames-by-assignment.sh
# bash ./bin/workflow/driver_allelic-segregation.sh \
# -u FALSE \
# -l FALSE \
# -d TRUE \
# -m "512m" \
# -x "4096m" \
# -r "mm10" \
# -s "CAST" \
# -1 "path/to/mm10/Disteche_sample_N.dedup.corrected.bam" \
# -2 "path/to/CAST/Disteche_sample_N.dedup.corrected.bam" \
# -p "Disteche_sample_1_downsampled_test-again" \
# -o "./path/to/directory/for/results" \
# -b 100000 \
# -c 1000000 \
# -t 0 \
# -a TRUE \
# -n 4

# Arguments:
# -h  print this help message and exit
# -u  use safe mode: TRUE or FALSE [logical; default: FALSE]
# -l  run on GS HPC: TRUE or FALSE [logical; default: FALSE]
# -d  run pipeline in ${TMPDIR}: TRUE or FALSE [logical; default:
#     TRUE]
# -m  initial memory allocation pool for JVM [chr; default: "512m"]
# -x  maximum memory allocation pool for JVM [chr; default: "4096m"]
# -r  string for "sample #1" [chr]
# -s  string for "sample #2" [chr]
# -1  bam infile #1, including path [chr]
# -2  bam infile #2, including path [chr]
# -p  prefix for outfiles [chr]
# -o  results directory for outfiles [chr]; path will be made if it
#     does not exist
# -b  number of records to read into memory at one time when running
#     the script for Part #1, get-AS-per-qname.R [int > 0; default:
#     100000]
# -c  number of records to read into memory at one time when running
#     the script for Part #3, generate-assignment-lists.R [int > 0;
#     default: 1000000]
# -t  alignment score threshold [int >= 0; default: 0]; the absolute
#     value of the difference in alignment scores between "sample
#     #1" and "sample #2" must be greater than this value in order
#     for a sample-specific assignment to be made; if not greater than
#     this value, then the assignment will be "ambiguous"
# -a  count lines: TRUE or FALSE [logical; default: TRUE]
# -n  step in pipeline to run up to [int 1-4; default: 4]

## 05.23
#  Call from 2021_kga0_4dn-mouse-cross or a directory containing
#+ functions-in-progress.sh and functions-preprocessing-HPC.sh
# bash ./remove-duplicate-qnames.sh \
# -u FALSE \
# -c TRUE \
# -m "2g" \
# -x "10g" \
# -i "${outpath}${infile}" \
# -o "${outpath}" \
# -n TRUE \
# -t FALSE \
# -e TRUE \
# -r TRUE \
# -p 4 \
# > "${outpath}/rm-dup-qnames_${strain}_${sample_id}.o.txt" \
# 2> "${outpath}/rm-dup-qnames_${strain}_${sample_id}.e.txt"

## 06.06
bash ./bin/workflow/driver_allelic-segregation.sh \
-u FALSE \
-l TRUE \
-d FALSE \
-m "512m" \
-x "82g" \
-r "mm10" \
-s "CAST" \
-1 "${path_1}${bam_input}" \
-2 "${path_2}${bam_input}" \
-p "Disteche_sample_${sample_id}_ALSO" \
-o "${outpath}" \
-b 100000 \
-c 1000000 \
-t 0 \
-a TRUE \
-n 4 \
> "${outpath}/ALSO_${sample_id}.o.txt" \
2> "${outpath}/ALSO_${sample_id}.e.txt"

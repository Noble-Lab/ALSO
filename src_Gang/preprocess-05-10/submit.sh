#!/bin/bash

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-05-03-preprocess ||
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

sample_id=$1
echo "${sample_id}"
strain=$2 # "CAST-EiJ"
echo "${strain}"
bam_input="Disteche_sample_${sample_id}.dedup.bam"
outpath="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-05-03-preprocess/results/${strain}-preprocessed/"

chmod 751 *.sh

#file="Disteche_sample_1.dedup.bam"
bash ./filter-qnames.sh \
-u FALSE \
-c FALSE \
-l TRUE \
-i ./"${strain}"/get_unique_fragments/"${bam_input}" \
-o "${outpath}" \
-f TRUE \
-r FALSE \
-p 4 > "${outpath}"preprocess_${strain}_${sample_id}.o 2>"${outpath}"preprocess_${strain}_${sample_id}.e




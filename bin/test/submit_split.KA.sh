#!/bin/bash

# split mm10
# GL 2022.03.24 / KA 2022-0326

qlogin -l mfree=2G

dir_storage="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results"
dir_test="2022-0326_mm10-split_sample-22_test"
cd "${dir_storage}" ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }

mkdir -p "${dir_test}"

# cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-24-split-mm10/ ||
#     {
#         echo "Exiting: Directory not found."
#         exit 1
#     }
cd "${TMPDIR}" ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }

module add bedtools/2.29.2 
module add samtools/1.14
module add parallel/20200922
export PATH="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/subread-2.0.3-Linux-x86_64/bin/utilities/:$PATH"

sample_id="22"
echo "${sample_id}"
safe_mode="FALSE"
infile="./Disteche_sample_${sample_id}.repair.bam"
outpath="."
prefix="Disteche_sample_${sample_id}"
chromosome="all"
mode="M"
parallelize=4

cp "/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-24-split-mm10/mm10-preprocessed/Disteche_sample_${sample_id}.repair.bam" .
cp "/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-24-split-mm10/split-index-repair-bam.sh" .

bash split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}" > "${outpath}"/split_sample_"${sample_id}".o 2>"${outpath}"/split_sample_"${sample_id}".e

outpath="${dir_storage}/${dir_test}"
mv -f "split_sample_${sample_id}."* "${prefix}."* "${outpath}"

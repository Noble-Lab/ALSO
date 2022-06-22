#!/bin/bash

#  test_generate-downsampled-bam.sh
#  KA


#  mm10 test ------------------------------------------------------------------
start="$(date +%s)"

sample="mm10"
infile="./data/files_bam_test/Disteche_sample_1.${sample}.dedup.bam"
outpath="./data/files_bam_test"
downsample=300000
prefix="test.${sample}.${downsample}.bam"
seed=24
parallelize=4

bash bin/generate-downsampled-bam.sh \
-u "FALSE" \
-i "${infile}" \
-o "${outpath}" \
-d "${downsample}" \
-x "${prefix}" \
-s "${seed}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""


#  CAST test ------------------------------------------------------------------
start="$(date +%s)"

sample="CAST-EiJ"
infile="./data/files_bam_test/Disteche_sample_1.${sample}.dedup.bam"
outpath="./data/files_bam_test"
downsample=300000
prefix="test.${sample}.${downsample}.bam"
seed=24
parallelize=4

bash bin/generate-downsampled-bam.sh \
-u "FALSE" \
-i "${infile}" \
-o "${outpath}" \
-d "${downsample}" \
-x "${prefix}" \
-s "${seed}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

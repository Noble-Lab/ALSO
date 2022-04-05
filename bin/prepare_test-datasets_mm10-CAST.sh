#!/bin/bash

#  prepare_test-datasets_mm10-CAST.sh
#  KA


#  Functions ------------------------------------------------------------------
# :param $1: a single file to grab, including path
sGrab() { scp -p nexus:"${1}" .; }


#  Check Disteche_sample_6.chr11.* --------------------------------------------
bam_HPC_mm10="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/get_unique_fragments/Disteche_sample_6.dedup.bam"
bam_HPC_CAST="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-23/get_unique_fragments/Disteche_sample_6.dedup.bam"
bam_mm10="Disteche_sample_6.dedup.mm10.bam"
bam_CAST="Disteche_sample_6.dedup.CAST.bam"
# bam_mm10_intermed="${bam_mm10/.bam/.intermediate.bam}"
# bam_CAST_intermed="${bam_CAST/.bam/.intermediate.bam}"
bam_mm10_prepro="${bam_mm10/.bam/.prepro.bam}"
bam_CAST_prepro="${bam_CAST/.bam/.prepro.bam}"
parallelize="${1:-"4"}"

4dn-data
mkdir -p ./data/2022-0404_prepare_test-datasets_mm10-CAST
cd ./data/2022-0404_prepare_test-datasets_mm10-CAST || ! echo "Warning: cd failed."

#  Get mm10 test bam
sGrab "${bam_HPC_mm10}" && \
mv -f "Disteche_sample_6.dedup.bam" "${bam_mm10}"

#  Get CAST test bam
sGrab ${bam_HPC_CAST} && \
mv -f "Disteche_sample_6.dedup.bam" "${bam_CAST}"

#  Check flag statistics for each bam file
samtools flagstat -@ "${parallelize}" "${bam_mm10}" > "${bam_mm10/.bam/.flagstat.txt}"
samtools flagstat -@ "${parallelize}" "${bam_CAST}" > "${bam_CAST/.bam/.flagstat.txt}"


#  Perform preprocessing ------------------------------------------------------
#  mm10
start="$(date +%s)"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_mm10}" \
| samtools sort -n -@ "${parallelize}" - \
| samtools fixmate -@ "${parallelize}" - - \
| samtools sort -@ "${parallelize}" - \
> "${bam_mm10_prepro}"
end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "mm10 run time: ${run_time} seconds."

samtools flagstat -@ "${parallelize}" "${bam_mm10_prepro}" > "${bam_mm10_prepro/.bam/.flagstat.txt}"

#  CAST
start="$(date +%s)"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_CAST}" \
| samtools sort -n -@ "${parallelize}" - \
| samtools fixmate -@ "${parallelize}" - - \
| samtools sort -@ "${parallelize}" - \
> "${bam_CAST_prepro}"
end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "CAST run time: ${run_time} seconds."

samtools flagstat -@ "${parallelize}" "${bam_CAST_prepro}" > "${bam_CAST_prepro/.bam/.flagstat.txt}"


#  Create split bams for tests ------------------------------------------------
samtools index -@ "${parallelize}" "${bam_mm10_prepro}"
samtools index -@ "${parallelize}" "${bam_CAST_prepro}"

samtools view -@ "${parallelize}" -b "${bam_mm10_prepro}" "chr19" > "${bam_mm10_prepro/.bam/.chr19.bam}"
samtools view -@ "${parallelize}" -b "${bam_CAST_prepro}" "chr19" > "${bam_CAST_prepro/.bam/.chr19.bam}"

samtools index -@ "${parallelize}" "${bam_mm10_prepro/.bam/.chr19.bam}"
samtools index -@ "${parallelize}" "${bam_CAST_prepro/.bam/.chr19.bam}"

repair -d -T "${parallelize}" -c \
-i "${bam_mm10_prepro/.bam/.chr19.bam}" \
-o "${bam_mm10_prepro/.bam/.chr19.repair.bam}"

repair -d -T "${parallelize}" -c \
-i "${bam_CAST_prepro/.bam/.chr19.bam}" \
-o "${bam_CAST_prepro/.bam/.chr19.repair.bam}"

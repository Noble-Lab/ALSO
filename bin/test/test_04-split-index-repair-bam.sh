#!/bin/bash

#  test_04-split-index-repair-bam.sh
#  KA


#  mm10 test ------------------------------------------------------------------
#  chr19
start="$(date +%s)"

file="./data/files_bam_test/test.mm10.300000.bam"
infile="./data/files_bam_test/test.mm10.300000.repair.bam"
safe_mode="FALSE"
outpath="./data/2022-0326_test_04_mm10_non-sorted_chr19"
prefix="test.mm10.300000"
chromosome="chr19"
mode="M"
parallelize=4

repair -T 4 -d -i "${file}" -o "${infile}"

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

rm "${infile}" "${infile}.bai"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

#  all
start="$(date +%s)"

file="./data/files_bam_test/test.mm10.300000.bam"
infile="./data/files_bam_test/test.mm10.300000.repair.bam"
safe_mode="FALSE"
outpath="./data/2022-0326_test_04_mm10_non-sorted_all"
prefix="test.mm10.300000"
chromosome="all"
mode="M"
parallelize=4

repair -T 4 -d -i "${file}" -o "${infile}"

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

rm "${infile}" "${infile}.bai"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 5 seconds
echo ""


#  CAST test ------------------------------------------------------------------
#  chr19
start="$(date +%s)"

file="./data/files_bam_test/test.CAST-EiJ.300000.bam"
infile="./data/files_bam_test/test.CAST-EiJ.300000.repair.bam"
safe_mode="FALSE"
outpath="./data/2022-0326_test_04_CAST-EiJ_non-sorted_chr19"
prefix="test.CAST-EiJ.300000"
chromosome="chr19"
mode="S"
parallelize=4

repair -T 4 -d -i "${file}" -o "${infile}"

#  Run in "default mode", which assumes fastqs were aligned to a non-mm10
#+ reference
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

rm "${infile}" "${infile}.bai"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

#  all
start="$(date +%s)"

file="./data/files_bam_test/test.CAST-EiJ.300000.bam"
infile="./data/files_bam_test/test.CAST-EiJ.300000.repair.bam"
safe_mode="FALSE"
outpath="./data/2022-0326_test_04_CAST-EiJ_non-sorted_all"
prefix="test.CAST-EiJ.300000"
chromosome="all"
mode="S"
parallelize=4

repair -T 4 -d -i "${file}" -o "${infile}"

#  Run in "default mode", which assumes fastqs were aligned to a non-mm10
#+ reference
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 9 seconds
echo ""


#  CAST-EiJ test for unsorted file ------------------------------------------------
#  chr19
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/Disteche_sample_1.CAST-EiJ.downsampled.shuffle.bam"
outpath="./data/2022-0327_test_04_CAST-EiJ_non-sorted_chr19"
prefix="Disteche_sample_1.CAST-EiJ.downsampled.shuffle"
chromosome="chr19"
mode="S"
parallelize=4

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""


#  CAST-EiJ test for sorted file, test removal of singletons ------------------
#  Just see whether or not the script runs to completion (even though there are
#+ no signletons in the sample data)

#  chr17
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/Disteche_sample_1.CAST-EiJ.downsampled.sort.bam"
outpath="./data/2022-0328_test_04_CAST-EiJ_rm-singletons_chr17"
prefix="Disteche_sample_1.CAST-EiJ.downsampled"
chromosome="chr17"
mode="S"
parallelize=4

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

#  all
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/Disteche_sample_1.CAST-EiJ.downsampled.sort.bam"
outpath="./data/2022-0328_test_04_CAST-EiJ_rm-singletons_all"
prefix="Disteche_sample_1.CAST-EiJ.downsampled"
chromosome="all"
mode="S"
parallelize=4

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mode}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

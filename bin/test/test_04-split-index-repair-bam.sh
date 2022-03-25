#!/bin/bash

#  test_04-split-index-repair-bam.sh
#  KA


#  mm10 test ------------------------------------------------------------------
#  chr19
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/test.mm10.300000.bam"
outpath="./data/2022-0324_test_04_chr19"
prefix="test.mm10.300000"
chromosome="chr19"
mm10="TRUE"
parallelize=4

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mm10}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

#  all
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/test.mm10.300000.bam"
outpath="./data/2022-0324_test_04_all"
prefix="test.mm10.300000"
chromosome="all"
mm10="TRUE"
parallelize=4

#  Run in "mm10 mode"
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-m "${mm10}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 5 seconds
echo ""


#  CAST test ------------------------------------------------------------------
#  chr19
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/test.CAST-EiJ.300000.bam"
outpath="./data/2022-0324_test_04_chr19"
prefix="test.CAST-EiJ.300000"
chromosome="chr19"
parallelize=4

#  Run in "default mode", which assumes fastqs were aligned to a non-mm10
#+ reference
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."
echo ""

#  all
start="$(date +%s)"

safe_mode="FALSE"
infile="./data/files_bam_test/test.CAST-EiJ.300000.bam"
outpath="./data/2022-0324_test_04_all"
prefix="test.CAST-EiJ.300000"
chromosome="all"
parallelize=4

#  Run in "default mode", which assumes fastqs were aligned to a non-mm10
#+ reference
bash bin/workflow/04-split-index-repair-bam.sh \
-u "${safe_mode}" \
-i "${infile}" \
-o "${outpath}" \
-x "${prefix}" \
-c "${chromosome}" \
-p "${parallelize}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo "${0} run time: ${run_time} seconds."  # run time: 9 seconds
echo ""

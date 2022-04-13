#!/bin/bash

#  test-preprocessing-gangliuw_2022-0411.sh
#  KA


#  Set up functions -----------------------------------------------------------
displaySpinningIcon() {
    # Display "spinning icon" while a background process runs
    # 
    # :param 1: PID of the last program the shell ran in the background (int)
    # :param 2: message to be displayed next to the spinning icon (chr)
    spin="/|\\–"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


#  Set up variables -----------------------------------------------------------
parallelize=6

bam_sciatac="Disteche_sample_6.dedup.CAST.bam"
bam_sort_n_1="${bam_sciatac/.bam/.sort-n.bam}"
bam_multiple="${bam_sort_n_1/.bam/.multiple-QNAME.bam}"
bam_sort_c_1="${bam_sort_n_1/.bam/.sort-c.bam}"
bam_filter="${bam_sort_c_1/.bam/.filter-QNAME.bam}"
bam_sort_n_2="${bam_filter/.bam/.sort-n.bam}"
bam_fixmate_1="${bam_sort_n_2/.bam/.fixmate.bam}"
bam_sort_c_2="${bam_fixmate_1/.bam/.sort-c.bam}"
bam_process="${bam_sort_c_2/.bam/.-f3-F12-q30.bam}"
bam_sort_n_3="${bam_process/.bam/.sort-n.bam}"
bam_fixmate_2="${bam_sort_n_3/.bam/.fixmate.bam}"
bam_sort_c_3="${bam_fixmate_2/.bam/.sort-c.bam}"

txt_multiple="${bam_sort_n_1/.bam/.multiple.txt}"
txt_multiple_QNAME="${txt_multiple/.txt/.QNAME.txt}"

flagstat_0="flagstat.0.before-processing.txt"
flagstat_1="flagstat.1.sort-n.txt"
flagstat_2="flagstat.2.only-multiple-QNAME.txt"
flagstat_3="flagstat.3.sort-n.sort-c.txt"
flagstat_4="flagstat.4.filter.txt"
flagstat_5="flagstat.5.filter.sort-n.txt"
flagstat_6="flagstat.6.filter.sort-n.fixmate.txt"
flagstat_7="flagstat.7.filter.sort-n.fixmate.sort-c.txt"
flagstat_8="flagstat.8.process.txt"
flagstat_9="flagstat.9.process.sort-n.fixmate.txt"
flagstat_10="flagstat.10.process.sort-n.fixmate.sort-c.txt"


#  Run flagstat on bam file output by the Shendure-Lab pipeline ---------------
samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}" &
displaySpinningIcon $! "Running flagstat #0 (bam output by Shendure-Lab pipeline)... "

#  1 - QNAME-sort bam file output by the Shendure-Lab pipeline ----------------
#+ ...then run flagstat again
start="$(date +%s)"

samtools sort -n -@ "${parallelize}" "${bam_sciatac}" > "${bam_sort_n_1}" &
displaySpinningIcon $! "Running samtools sort -n #1... "
samtools flagstat -@ "${parallelize}" "${bam_sort_n_1}" > "${flagstat_1}" &
displaySpinningIcon $! "Running flagstat #1 (after doing first sort -n)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 2281 seconds, 38.02 minutes
echo ""


#  2 - Identify QNAMEs with >2 entries ----------------------------------------
start="$(date +%s)"

samtools view "${bam_sort_n_1}" \
| cut -f 1 \
| sort \
| uniq -c \
| sort -nr \
| awk '$1 > 2 {print $0}' \
> "${txt_multiple}" &
displaySpinningIcon $! "Running samtools view to identify abnormal numbers of QNAME entries... "

#  Create a single-column version of the txt file containing QNAMEs with >2
#+ entries
cut -c 6- "${txt_multiple}" > "${txt_multiple_QNAME}"

#  gzip the txt file of QNAMEs with >2 entries
gzip "${txt_multiple}"

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 2281 seconds, 38.02 minutes
echo ""

#  3 - Filter bam file to include only QNAMEs with >2 entries -----------------
start="$(date +%s)"

samtools view -hN "${txt_multiple_QNAME}" "${bam_sort_n_1}" > "${bam_multiple}" &
displaySpinningIcon $! "Running samtools view -N (filter bam file to include only QNAMEs with >2 entries)... "

samtools flagstat -@ "${parallelize}" "${bam_multiple}" > "${flagstat_2}" &
displaySpinningIcon $! "Running flagstat #2 (for bam file that includes only QNAMEs with >2 entries)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 60 seconds
echo ""


#  4 - Filter bam file to exclude QNAMEs with >2 entries --------------------------
#  Coordinate-sort the QNAME-sorted bam file
start="$(date +%s)"

samtools sort -@ "${parallelize}" "${bam_sort_n_1}" > "${bam_sort_c_1}" &
displaySpinningIcon $! "Sorting bam by coordinate in preparation for use with picard FilterSamReads... "

samtools flagstat -@ "${parallelize}" "${bam_sort_c_1}" > "${flagstat_3}"
displaySpinningIcon $! "Running flagstat #3 (sort by QNAME-sorted bam by coordinates again)... "

#  Filter the coordinate-sorted bam file
picard FilterSamReads \
I="${bam_sort_c_1}" \
O="${bam_filter}" \
READ_LIST_FILE="${txt_multiple_QNAME}" \
FILTER="excludeReadList" &
displaySpinningIcon $! "Using picard to filter bam file for abnormal QNAME entries... "

samtools flagstat -@ "${parallelize}" "${bam_filter}" > "${flagstat_4}" &
displaySpinningIcon $! "Running flagstat #4 (duplicate-QNAME-filtered bam file)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 3584 seconds, 59.7 minutes
echo ""


#  QNAME-sort the filtered bam file -------------------------------------------
start="$(date +%s)"

samtools sort -n -@ "${parallelize}" "${bam_filter}" > "${bam_sort_n_2}" &
displaySpinningIcon $! "Running samtools sort -n #2... "

samtools flagstat -@ "${parallelize}" "${bam_sort_n_2}" > "${flagstat_5}" &
displaySpinningIcon $! "Running flagstat #5 (QNAME-sorted duplicate-QNAME-filtered bam file)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 755 seconds, 12.58 minutes
echo ""

#  Free up some space
echo "Removing ${bam_sort_n_1} and ${bam_sort_c_1} to free up space... "
rm "${bam_sort_n_1}" "${bam_sort_c_1}"


#  Perform fixmate on the filtered, QNAME-sorted bam file ---------------------
start="$(date +%s)"

samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_1}" &
displaySpinningIcon $! "Running samtools fixmate on the filtered, QNAME-sorted bam file... "

samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_6}" &
displaySpinningIcon $! "Running flagstat #6 (fixmate has been performed on filtered, QNAME-sorted bam file)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 940 seconds, 15.67 minutes
echo ""


#  Sort by coordinate again prior to filtering out reads... ------------------- 
#+ ...based on status and MAPQ
start="$(date +%s)"

samtools sort -@ "${parallelize}" "${bam_fixmate_1}" > "${bam_sort_c_2}" &
displaySpinningIcon $! "Sorting bam by coordinate in preparation for filtering reads for paired status and MAPQ... "

samtools flagstat -@ "${parallelize}" "${bam_sort_c_2}" > "${flagstat_7}" &
displaySpinningIcon $! "Running flagstat #7 (coordinate-sort the fixmate-status bam)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # Not recorded
echo ""

#  Free up some space
echo "Removing ${bam_sort_n_2} to free up space... "
rm "${bam_sort_n_2}"


#  Filter out reads based on paired status and MAPQ ---------------------------
start="$(date +%s)"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sort_c_2}" -o "${bam_process}" &
displaySpinningIcon $! "Running samtools view (-f 3 -F 12 -q 30)... "

samtools flagstat -@ "${parallelize}" "${bam_process}" > "${flagstat_8}" &
displaySpinningIcon $! "Running flagstat #8 (bam has been processed for properly paired, MAPQ >= 30 reads)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 78 seconds
echo ""

#  Free up some space
echo "Removing ${bam_sort_c_2} to free up space... "
rm "${bam_sort_c_2}"


#  QNAME-sort bam file and run fixmate one last time --------------------------
start="$(date +%s)"
samtools sort -n -@ "${parallelize}" "${bam_process}" > "${bam_sort_n_3}" &
displaySpinningIcon $! "Running samtools sort -n #2... "

samtools fixmate -@ "${parallelize}" "${bam_sort_n_3}" "${bam_fixmate_2}" &
displaySpinningIcon $! "Running samtools fixmate #2... "

samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_9}" &
displaySpinningIcon $! "Running flagstat #9 (processed bam has been QNAME-sorted and subjected to fixmate)... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 702 seconds
echo ""

#  Free up some space
echo "Removing ${bam_sort_n_3} to free up space... "
rm "${bam_sort_n_3}"


#  Coordinate-sort the processed fixmate-status bam ---------------------------
start="$(date +%s)"

samtools sort -@ "${parallelize}" "${bam_fixmate_2}" > "${bam_sort_c_3}" &
displaySpinningIcon $! "Running samtools sort (by coordinate) #2... "

samtools flagstat -@ "${parallelize}" "${bam_sort_c_3}" > "${flagstat_10}" &
displaySpinningIcon $! "Running flagstat #5... "

end="$(date +%s)"
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds."  # 274 seconds
echo ""

#  Free up some space
echo "Removing ${bam_fixmate_1} and ${bam_fixmate_2} to free up space... "
rm "${bam_fixmate_1}" "${bam_fixmate_2}"


# #  Example input, output ------------------------------------------------------
# parallelize=6
# bam_sciatac="${1:-"Disteche_sample_6.dedup.bam"}"
# bam_preprocessed="${2:-"Disteche_sample_6.dedup.preprocessed.bam"}"
#
# #  Pipeline
# samtools sort -n -@ "${parallelize}" "${bam_sciatac}" \
# | samtools fixmate -@ "${parallelize}" - \
# | samtools sort -@ "${parallelize}" - \
# | samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 - \
# | samtools sort -n -@ "${parallelize}" - \
# | samtools fixmate -@ "${parallelize}" - \
# | samtools sort -@ "${parallelize}" "${bam_preprocessed}"


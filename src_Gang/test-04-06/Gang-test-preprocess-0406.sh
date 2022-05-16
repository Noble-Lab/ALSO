#!/bin/bash

#  test-preprocessing-gangliuw_2022-0404.sh
#  KA
# GL : try f1

#  Start recording time
start="$(date +%s)"

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools ||
    {
        echo "Exiting: Directory not found."
        exit 1
    }
module add samtools/1.14

#  Function
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

# #  Example input, output
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

#  Example input, output
parallelize=4

bam_sciatac="mm10/get_unique_fragments/Disteche_sample_6.dedup.bam"
bam_preprocessed="mm10-output/test/Disteche_sample_6.dedup.preprocessed.bam"
output_path="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools/"

# bam_sciatac="Disteche_sample_6.dedup.bam"
bam_sort_n_1="${output_path}mm10-output/test-0406/sort-n.bam"
bam_fixmate_1="${output_path}mm10-output/test-0406/fixmate.bam"
bam_sort_c_1="${output_path}mm10-output/test-0406/sort-c.bam"
bam_process="${output_path}mm10-output/test-0406/f1-F12-q30.bam"
bam_sort_n_2="${output_path}mm10-output/test-0406/sort-n_2.bam"
bam_fixmate_2="${output_path}mm10-output/test-0406/fixmate_2.bam"
bam_sort_c_2="${output_path}mm10-output/test-0406/sort-c_2.bam"

flagstat_0="mm10-output/test-0406/flagstat.0.txt"
flagstat_1="mm10-output/test-0406/flagstat.1.txt"
flagstat_2="mm10-output/test-0406/flagstat.2.txt"
flagstat_3="mm10-output/test-0406/flagstat.3.txt"
flagstat_4="mm10-output/test-0406/flagstat.4.txt"
flagstat_5="mm10-output/test-0406/flagstat.5.txt"

#  Pipeline
# samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}" &
# displaySpinningIcon $! "Running flagstat #0... "

samtools sort -n -@ "${parallelize}" "${bam_sciatac}" -o ${bam_sort_n_1} &
displaySpinningIcon $! "Running samtools sort -n #1... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_1}" "${bam_fixmate_1}" &
displaySpinningIcon $! "Running samtools fixmate #1... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_1}" &
displaySpinningIcon $! "Running flagstat #1... "

# rm "${bam_sort_n_1}"

samtools sort -@ "${parallelize}" "${bam_fixmate_1}" -o "${bam_sort_c_1}" &
displaySpinningIcon $! "Running samtools sort (by coordinate) #1... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_1}" > "${flagstat_2}" &
displaySpinningIcon $! "Running flagstat #2... "

samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sort_c_1}" -o "${bam_process}" &
displaySpinningIcon $! "Running samtools view (-f 3 -F 12 -q 30)... "
samtools flagstat -@ "${parallelize}" "${bam_process}" > "${flagstat_3}" &
displaySpinningIcon $! "Running flagstat #3... "

# rm "${bam_sort_c_1}"

samtools sort -n -@ "${parallelize}" "${bam_process}" -o "${bam_sort_n_2}" &
displaySpinningIcon $! "Running samtools sort -n #2... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_2}" &
displaySpinningIcon $! "Running samtools fixmate #2... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_4}" &
displaySpinningIcon $! "Running flagstat #4... "

# rm "${bam_sort_n_2}"

samtools sort -@ "${parallelize}" "${bam_fixmate_2}" -o "${bam_sort_c_2}" &
displaySpinningIcon $! "Running samtools sort (by coordinate) #2... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_2}" > "${flagstat_5}" &
displaySpinningIcon $! "Running flagstat #5... "

#  End recording time
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""
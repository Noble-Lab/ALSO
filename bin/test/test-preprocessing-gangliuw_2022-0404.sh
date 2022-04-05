#!/bin/bash

#  test-preprocessing-gangliuw_2022-0404.sh
#  KA

#  Start recording time
start="$(date +%s)"

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
parallelize=6

bam_sciatac="Disteche_sample_6.dedup.bam"
bam_sort_n_1="${bam_sciatac/.bam/.sort-n.bam}"
bam_fixmate_1="${bam_sort_n_1/.bam/.fixmate.bam}"
bam_sort_c_1="${bam_fixmate_1/.bam/.sort-c.bam}"
bam_process="${bam_sort_c_1/.bam/.-f3-F12-q30.bam}"
bam_sort_n_2="${bam_process/.bam/.sort-n.bam}"
bam_fixmate_2="${bam_sort_n_2/.bam/.fixmate.bam}"
bam_sort_c_2="${bam_fixmate_2/.bam/.sort-c.bam}"

flagstat_0="flagstat.0.txt"
flagstat_1="flagstat.1.txt"
flagstat_2="flagstat.2.txt"
flagstat_3="flagstat.3.txt"
flagstat_4="flagstat.4.txt"
flagstat_5="flagstat.5.txt"

#  Pipeline
samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}" &
displaySpinningIcon $! "Running flagstat #0... "

samtools sort -n -@ "${parallelize}" "${bam_sciatac}" > "${bam_sort_n_1}" &
displaySpinningIcon $! "Running samtools sort -n #1... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_1}" "${bam_fixmate_1}" &
displaySpinningIcon $! "Running samtools fixmate #1... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_1}" &
displaySpinningIcon $! "Running flagstat #1... "

rm "${bam_sort_n_1}"

samtools sort -@ "${parallelize}" "${bam_fixmate_1}" > "${bam_sort_c_1}" &
displaySpinningIcon $! "Running samtools sort (by coordinate) #1... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_1}" > "${flagstat_2}" &
displaySpinningIcon $! "Running flagstat #2... "

samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sort_c_1}" -o "${bam_process}" &
displaySpinningIcon $! "Running samtools view (-f 3 -F 12 -q 30)... "
samtools flagstat -@ "${parallelize}" "${bam_process}" > "${flagstat_3}" &
displaySpinningIcon $! "Running flagstat #3... "

rm "${bam_sort_c_1}"

samtools sort -n -@ "${parallelize}" "${bam_process}" > "${bam_sort_n_2}" &
displaySpinningIcon $! "Running samtools sort -n #2... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_2}" &
displaySpinningIcon $! "Running samtools fixmate #2... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_4}" &
displaySpinningIcon $! "Running flagstat #4... "

rm "${bam_sort_n_2}"

samtools sort -@ "${parallelize}" "${bam_fixmate_2}" > "${bam_sort_c_2}" &
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

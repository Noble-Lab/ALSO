#!/bin/bash

#  test-preprocessing-gangliuw_2022-0411.sh
#  KA


#  Start recording time (total) -----------------------------------------------
total_start="$(date +%s)"


#  Set up functions -----------------------------------------------------------
calculateRunTime() {
    # Calculate run time for processes
    # 
    # :param 1: start time in $(date +%s) format
    # :param 2: end time in $(date +%s) format
    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo "Run time: ${run_time} seconds."
    echo ""
}


checkDependency() {
    # Check if program is available in "${PATH}"; exit if not
    # 
    # :param 1: program to be checked (chr)
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            exit 1
        }
}


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


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Run bam file through the steps of the test 2021_kga0_4dn-mouse-cross"
    echo "preprocessing module as of 2022-0411-0413."
    echo ""
    echo "Preprocessing is made up of the following steps:"
    echo "    1. QNAME-sort bam file output by the Shendure-Lab pipeline (this"
    echo "       is not necessary but speeds up the following step"
    echo "       considerably)"
    echo "    2. Identify QNAMEs with >2 entries in the bam file (this step is"
    echo "       slow; there's likely room for optimization)"
    echo "    3. Create a bam file comprised of only duplicated QNAME entries"
    echo "       (this step is optional; it's not strictly necessary)"
    echo "    4. Filter bam file to exclude QNAMEs with >2 entries by"
    echo "       (a) coordinate-sorting the QNAME-sorted bam file and"
    echo "       (b) filtering the coordinate-sorted bam file with picard"
    echo "           FilterSamReads (picard FilterSamReads takes only"
    echo "           coordinate-sorted bam files as input; also, picard"
    echo "           FilterSamReads is exponentially faster than filtering"
    echo "           with grep)"
    echo "    5. QNAME-sort the filtered bam file; then, to update flag"
    echo "       information in the bam file, perform samtools fixmate on the"
    echo "       bam file"
    echo "    6. Prior to filtering out reads based on pairing status and MAPQ"
    echo "       values, sort by coordinate again"
    echo "    7. Filter out reads based on pairing status and MAPQ: Run" 
    echo "       samtools view with flags -f 3 -F 12 -q 30"
    echo "    8. Sort bam by QNAME and perform samtools fixmate to update flag"
    echo "       information again"
    echo "    9. Coordinate-sort and index the processed bam file, which"
    echo "       should now be ready for the subsequent module"
    echo ""
    echo "Dependencies:"
    echo "    - picard >= 2.26.4"
    echo "    - samtools >= 1.13"
    echo ""
    echo "Arguments:"
    echo "    -h <print this help message and exit>"
    echo "    -u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
    echo "    -i <bam infile, including path (chr)>"
    echo "    -o <directory for outfile, including path (chr)>"
    echo "    -p <number of cores for parallelization (int >= 1); default: 1>"
    echo ""
    exit
}


while getopts "h:u:i:o:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${parallelize}" ]] && parallelize=6


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Check for necessary dependencies; exit if not found
checkDependency picard
checkDependency samtools

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-u: Safe mode is on."
        set -Eeuxo pipefail
        ;;
    false | f) \
        echo -e "-u: Safe mode is off."
        :
        ;;
    *) \
        echo -e "Exiting: -u safe-mode argument must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

#  Check that "${infile}" exists
[[ -f "${infile}" ]] ||
    {
        echo -e "Exiting: -i ${infile} does not exist.\n"
        exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo "-o: Directory ${outpath} does not exist; making the directory."
        mkdir -p "${outpath}"
    }

#  Check "${parallelize}"
[[ ! "${parallelize}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -p parallelize argument must be an integer.\n"
        exit 1
    }

[[ ! $((parallelize)) -ge 1 ]] &&
    {
        echo -e "Exiting: -p parallelize argument must be an integer >= 1.\n"
        exit 1
    }


#  Set up variables -----------------------------------------------------------
bam_sciatac="$(basename "${infile}")"
bam_sort_n_1="${outpath}/${bam_sciatac/.bam/.sort-n.bam}"
bam_multiple="${outpath}/${bam_sort_n_1/.bam/.multiple-QNAME.bam}"
bam_sort_c_1="${outpath}/${bam_sort_n_1/.bam/.sort-c.bam}"
bam_filter="${outpath}/${bam_sort_c_1/.bam/.filter-QNAME.bam}"
bam_sort_n_2="${outpath}/${bam_filter/.bam/.sort-n.bam}"
bam_fixmate_1="${outpath}/${bam_sort_n_2/.bam/.fixmate.bam}"
bam_sort_c_2="${outpath}/${bam_fixmate_1/.bam/.sort-c.bam}"
bam_process="${outpath}/${bam_sort_c_2/.bam/.-f3-F12-q30.bam}"
bam_sort_n_3="${outpath}/${bam_process/.bam/.sort-n.bam}"
bam_fixmate_2="${outpath}/${bam_sort_n_3/.bam/.fixmate.bam}"
bam_sort_c_3="${outpath}/${bam_fixmate_2/.bam/.sort-c.bam}"

txt_multiple="${outpath}/${bam_sort_n_1/.bam/.multiple.txt}"
txt_multiple_QNAME="${outpath}/${txt_multiple/.txt/.QNAME.txt}"

flagstat_0="${outpath}/flagstat.0.before-processing.txt"
flagstat_1="${outpath}/flagstat.1.sort-n.txt"
flagstat_2="${outpath}/flagstat.2.only-multiple-QNAME.txt"
flagstat_3="${outpath}/flagstat.3.sort-n.sort-c.txt"
flagstat_4="${outpath}/flagstat.4.filter.txt"
flagstat_5="${outpath}/flagstat.5.filter.sort-n.txt"
flagstat_6="${outpath}/flagstat.6.filter.sort-n.fixmate.txt"
flagstat_7="${outpath}/flagstat.7.filter.sort-n.fixmate.sort-c.txt"
flagstat_8="${outpath}/flagstat.8.process.txt"
flagstat_9="${outpath}/flagstat.9.process.sort-n.fixmate.txt"
flagstat_10="${outpath}/flagstat.10.process.sort-n.fixmate.sort-c.txt"


#  Step 0 - Run flagstat on bam infile ----------------------------------------
start="$(date +%s)"

samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}" &
displaySpinningIcon $! "Running flagstat #0 (bam output by Shendure-Lab pipeline)... "

end="$(date +%s)"
echo ""
echo "Step 0 - Run flagstat on bam infile"
calculateRunTime "${start}" "${end}"  # 43 seconds


#  Step 1 - QNAME-sort bam infile ---------------------------------------------
#+ ...then run flagstat again
start="$(date +%s)"

samtools sort -n -@ "${parallelize}" "${bam_sciatac}" > "${bam_sort_n_1}" &
displaySpinningIcon $! "Running samtools sort -n #1... "
samtools flagstat -@ "${parallelize}" "${bam_sort_n_1}" > "${flagstat_1}" &
displaySpinningIcon $! "Running flagstat #1 (after doing first sort -n)... "

end="$(date +%s)"
echo ""
echo "Step 1 - QNAME-sort bam infile"
calculateRunTime "${start}" "${end}"  # 665 seconds, 11.08 minutes (HPC)


#  Step 2 - Identify QNAMEs with >2 entries -----------------------------------
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
echo ""
echo "Step 2 - Identify QNAMEs with >2 entries"
calculateRunTime "${start}" "${end}"  # 2281 seconds, 38.02 minutes (2020 MacBook Pro M1); 4070 seconds, 67.83 minutes (HPC)


#  Step 3 - Filter bam file to include only QNAMEs with >2 entries ------------
start="$(date +%s)"

samtools view -hN "${txt_multiple_QNAME}" "${bam_sort_n_1}" > "${bam_multiple}" &
displaySpinningIcon $! "Running samtools view -N (filter bam file to include only QNAMEs with >2 entries)... "

samtools flagstat -@ "${parallelize}" "${bam_multiple}" > "${flagstat_2}" &
displaySpinningIcon $! "Running flagstat #2 (for bam file that includes only QNAMEs with >2 entries)... "

end="$(date +%s)"
echo ""
echo "Step 3 - Filter bam file to include only QNAMEs with >2 entries"
calculateRunTime "${start}" "${end}"  # 60 seconds (2020 MacBook Pro M1); 91 seconds, 1.52 minutes (HPC)


#  Step 4 - Filter bam file to exclude QNAMEs with >2 entries -----------------
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
displaySpinningIcon $! "Using picard to filter bam file for duplicate QNAME entries... "

samtools flagstat -@ "${parallelize}" "${bam_filter}" > "${flagstat_4}" &
displaySpinningIcon $! "Running flagstat #4 (duplicate-QNAME-filtered bam file)... "

end="$(date +%s)"
echo ""
echo "Step 4 - Filter bam file to exclude QNAMEs with >2 entries"
calculateRunTime "${start}" "${end}"  # 3584 seconds, 59.7 minutes (2020 MacBook Pro M1); 1461 seconds, 24.35 minutes (HPC)


#  Step 5 - QNAME-sort and run fixmate on the filtered bam file ---------------
start="$(date +%s)"

samtools sort -n -@ "${parallelize}" "${bam_filter}" > "${bam_sort_n_2}" &
displaySpinningIcon $! "Running samtools sort -n #2... "

samtools flagstat -@ "${parallelize}" "${bam_sort_n_2}" > "${flagstat_5}" &
displaySpinningIcon $! "Running flagstat #5 (QNAME-sorted duplicate-QNAME-filtered bam file)... "

#  Free up some space
echo "Removing ${bam_sort_n_1} and ${bam_sort_c_1} to free up space... "
rm "${bam_sort_n_1}" "${bam_sort_c_1}" && echo ""

samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_1}" &
displaySpinningIcon $! "Running samtools fixmate on the filtered, QNAME-sorted bam file... "

samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_6}" &
displaySpinningIcon $! "Running flagstat #6 (fixmate has been performed on filtered, QNAME-sorted bam file)... "

end="$(date +%s)"
echo ""
echo "Step 5 - QNAME-sort and run fixmate on the filtered bam file"
calculateRunTime "${start}" "${end}"  # 940 seconds, 15.67 minutes (2020 MacBook Pro M1); 932 seconds, 15.53 minutes (HPC)


#  Step 6 - Sort by coordinate again prior to filtering out reads... ----------
#+ ...based on pairing status and MAPQ
start="$(date +%s)"

samtools sort -@ "${parallelize}" "${bam_fixmate_1}" > "${bam_sort_c_2}" &
displaySpinningIcon $! "Sorting bam by coordinate in preparation for filtering reads for paired status and MAPQ... "

samtools flagstat -@ "${parallelize}" "${bam_sort_c_2}" > "${flagstat_7}" &
displaySpinningIcon $! "Running flagstat #7 (coordinate-sort the fixmate-status bam)... "

#  Free up some space
echo "Removing ${bam_sort_n_2} to free up space... "
rm "${bam_sort_n_2}" && echo ""

end="$(date +%s)"
echo ""
echo "Step 6 - Sort by coordinate again prior to filtering out reads"
calculateRunTime "${start}" "${end}"  # 522 seconds, 8.70 minutes (HPC)


#  Step 7 - Filter out reads based on paired status and MAPQ scores -----------
start="$(date +%s)"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sort_c_2}" -o "${bam_process}" &
displaySpinningIcon $! "Running samtools view (-f 3 -F 12 -q 30)... "

samtools flagstat -@ "${parallelize}" "${bam_process}" > "${flagstat_8}" &
displaySpinningIcon $! "Running flagstat #8 (bam has been processed for properly paired, MAPQ >= 30 reads)... "

#  Free up some space
echo "Removing ${bam_sort_c_2} to free up space... "
rm "${bam_sort_c_2}" && echo ""

end="$(date +%s)"
echo ""
echo "Step 7 - Filter out reads based on paired status and MAPQ scores"
calculateRunTime "${start}" "${end}"  # 78 seconds (2020 MacBook Pro M1); 173 seconds, 2.88 minutes (HPC)


#  Step 8 - QNAME-sort and run fixmate on bam file one last time --------------
start="$(date +%s)"

samtools sort -n -@ "${parallelize}" "${bam_process}" > "${bam_sort_n_3}" &
displaySpinningIcon $! "Running samtools sort -n #2... "

samtools fixmate -@ "${parallelize}" "${bam_sort_n_3}" "${bam_fixmate_2}" &
displaySpinningIcon $! "Running samtools fixmate #2... "

samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_9}" &
displaySpinningIcon $! "Running flagstat #9 (processed bam has been QNAME-sorted and subjected to fixmate)... "

#  Free up some space
echo "Removing ${bam_sort_n_3} to free up space... "
rm "${bam_sort_n_3}" && echo ""

end="$(date +%s)"
echo ""
echo "Step 8 - QNAME-sort and run fixmate on bam file one last time"
calculateRunTime "${start}" "${end}"  # 702 seconds, 11.7 minutes (2020 MacBook Pro M1); 791 seconds, 13.81 minutes (HPC)


#  Step 9 - Coordinate-sort the processed fixmate-status bam ------------------
start="$(date +%s)"

samtools sort -@ "${parallelize}" "${bam_fixmate_2}" > "${bam_sort_c_3}" &
displaySpinningIcon $! "Running samtools sort (by coordinate) #2... "

samtools flagstat -@ "${parallelize}" "${bam_sort_c_3}" > "${flagstat_10}" &
displaySpinningIcon $! "Running flagstat #5... "

#  Free up some space
echo "Removing ${bam_fixmate_1} and ${bam_fixmate_2} to free up space... "
rm "${bam_fixmate_1}" "${bam_fixmate_2}" && echo ""

#NOTE For now, keep these intermediate files for troubleshooting purposes
# echo "Removing ${bam_multiple}, ${bam_filter}, ${bam_process} to free up space... "
# rm "${bam_multiple}" "${bam_filter}" "${bam_process}" && echo ""

end="$(date +%s)"
echo ""
echo "Step 9 - Coordinate-sort the processed fixmate-status bam"
calculateRunTime "${start}" "${end}"  # 274 seconds, 4.57 minutes (2020 MacBook Pro M1); 452 seconds, 7.53 minutes (HPC)


#  End recording time (total) -------------------------------------------------
total_end="$(date +%s)"
echo "Calculating total run time..."
calculateRunTime "${total_start}" "${total_end}"  # 9200 seconds, 153.33 minutes, 2.56 hours (HPC)

exit 0

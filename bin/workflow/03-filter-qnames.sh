#!/bin/bash

#  filter-qnames.sh
#  KA


time_start="$(date +%s)"


#  Set working directory and source functions --------------------------------- 
#  Check for proper "$(pwd)"
if [[ "$(basename "$(pwd)")" != "2021_kga0_4dn-mouse-cross" ]]; then
    echo "Exiting: You must run this script from" \
    "\"2021_kga0_4dn-mouse-cross\", the project's top directory."
    exit 1
fi

#  Source functions into environment
# shellcheck disable=1091
. ./bin/auxiliary/functions_preprocessing.sh ||
    {
        echo "Exiting: Unable to source auxiliary information."
        echo "Are you in the correct working directory," \
        "\"2021_kga0_4dn-mouse-cross\"?"
        exit 1
    }

#  Additional functions
echo_completion_file() {
    # #TODO Description of function
    #
    # :param 1: outpath (chr)
    # :param 2: step number (int)
    # :param 3: bam infile (chr)
    echo "${1}/filter-qnames.step-${2}.$(basename "${3/.bam/.txt}")"
}


echo_completion_message() {
    # #TODO Description of function
    #
    # :param 1: step number (int)
    echo "Step ${1} completed; moving to next step..."
}


echo_exit_message() {
    # #TODO Description of function
    #
    # :param 1: step number (int)
    echo "Exiting: Step ${1}."
}


echo_flagstat() {
    # Echo a filename then print the first 25 lines of the file
    #
    # :param 1: file, e.g., a txt file (chr)
    echo ""
    echo "${1}"
    head -25 "${1}"
}


echo_loop() { for i in "${@:-*}"; do echo "${i}"; done; }


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Run pipeline to filter problematic QNAMEs from bam file."
    echo ""
    echo ""
    echo "Dependencies:"
    echo " - argparser >= 0.7.1"
    echo " - picard >= 2.27.1"
    echo " - R >= 4.0.5"
    echo " - Rsamtools >= 2.6.0"
    echo " - Rscript >= 4.0.5"
    echo " - samtools >= 1.13"
    echo " - scales >= 1.1.1"
    echo " - tidyverse >= 1.3.1"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-u use safe mode: \"TRUE\" or \"FALSE\" (logical)"
    echo "-c run on GS HPC: \"TRUE\" or \"FALSE\" (logical)"
    echo "-i bam infile, including path (chr)"
    echo "-o path for outfile(s; chr); path will be made if it does not exist"
    # echo "-r remove intermediate files: \"TRUE\" or \"FALSE\" (logical)"
    echo "-p number of cores for parallelization (int >= 1); default: 1"
    exit
}

# while getopts "h:u:c:i:o:r:p:" opt; do
while getopts "h:u:c:i:o:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        c) cluster="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        # r) remove="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${cluster}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
# [[ -z "${remove}" ]] && printUsage  #TODO Build out code for this...
[[ -z "${parallelize}" ]] && parallelize=1


#  Assign variables for completion files --------------------------------------
step_1="$(echo_completion_file "${outpath}" 1 "${infile}")"
step_2="$(echo_completion_file "${outpath}" 2 "${infile}")"
step_3="$(echo_completion_file "${outpath}" 3 "${infile}")"
step_4="$(echo_completion_file "${outpath}" 4 "${infile}")"
step_5="$(echo_completion_file "${outpath}" 5 "${infile}")"
step_6="$(echo_completion_file "${outpath}" 6 "${infile}")"
step_7="$(echo_completion_file "${outpath}" 7 "${infile}")"
step_8="$(echo_completion_file "${outpath}" 8 "${infile}")"
step_9="$(echo_completion_file "${outpath}" 9 "${infile}")"


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Check for necessary dependencies; exit if not found
check_dependency R
check_dependency samtools

#  Evaluate "${cluster}"
case "$(echo "${cluster}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-c: TRUE"
        :
        ;;
    false | f) \
        echo -e "-c: FALSE"
        check_dependency picard
        ;;
    *) \
        echo -e "Exiting: -c argument must be \"TRUE\" or \"FALSE\".\n"
        return 1
        ;;
esac

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

echo ""


#  1: Remove low-quality reads (-f 3 -F 12 -q 30) -----------------------------
if [[ ! -f "${step_1}" ]]; then
    remove_reads_low_quality_auto "${parallelize}" "${infile}" && \
    touch "${step_1}"
else
    echo_completion_message 1
    :
fi


#  2: Sort bam by coordinate, writing a tmp file that overwrites the infile ---
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" \
    "${infile/.bam/.rm.bam}" && \
    touch "${step_2}"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
    :
else
    echo_exit_message 2
    exit 1
fi


#  3: Index bam file ----------------------------------------------------------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    index_bam "${parallelize}" "${infile/.bam/.rm.bam}" && \
    touch "${step_3}"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
    :
else
    echo_exit_message 3
    exit 1
fi


#  4: Generate lists of QNAMEs to exclude -------------------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    Rscript ./bin/generate-qname-lists.R \
    --bam "${infile/.bam/.rm.bam}" \
    --bai "${infile/.bam/.rm.bam.bai}" \
    --outdir "${outpath}" \
    --chunk 100000 \
    --mated FALSE \
    --unmated TRUE \
    --ambiguous TRUE \
    --trans TRUE \
    --duplicated TRUE \
    --singleton TRUE \
    --unique TRUE \
    --tally TRUE \
    --remove TRUE && \
    touch "${step_4}"
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
    :
else
    echo_exit_message 4
    exit 1
fi


#  5: Tally number of records in outfiles -------------------------------------
unset outfiles
typeset -a outfiles
while IFS=" " read -r -d $'\0'; do
    outfiles+=( "${REPLY}" )
done < <(
    find "${outpath}" \
    -maxdepth 1 \
    -type f \
    -name "$(basename "${infile/.bam/.rm}").*.txt.gz" \
    -not -name "*tally*" \
    -print0
)
# echo_loop "${outfiles[@]}"; echo ""

if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    for i in "${outfiles[@]}"; do
        echo "$(basename "${i}"): $(zcat < "${i}" | wc -l)"
    done && \
    touch "${step_5}"
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    echo_completion_message 5
    :
else
    echo_exit_message 5
    exit 1
fi


#  6: Combine outfiles into file used for exclusion ---------------------------
if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then
    if [[ "${#outfiles[@]}" -eq 4 ]]; then
        combine_gz_qname_lists_return_unique_gzip \
        "${cluster}" "${outfiles[@]}" \
        > "${outpath}/$(basename "${infile/.bam/.rm}").to-exclude.txt.gz" && \
        touch "${step_6}"
    else
        echo "Exiting: Step 6: Number of outfiles does not equal 4."
        exit 1
    fi
elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
    echo_completion_message 6
    :
else
    echo_exit_message 6
    exit 1
fi


#  7: Tally number of records in "to-exclude" file ----------------------------
if [[ ! -f "${step_7}" && -f "${step_6}" ]]; then
    gzcat "${outpath}/$(basename "${infile/.bam/.rm}").to-exclude.txt.gz" | wc -l && \
    touch "${step_7}"
elif [[ -f "${step_7}" && -f "${step_6}" ]]; then
    echo_completion_message 7
    :
else
    echo_exit_message 7
    exit 1
fi


#  8: Exclude problematic QNAME reads from bam infile -------------------------
if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
    exclude_qname_reads_picard \
    "${infile/.bam/.rm.bam}" \
    "${outpath}/$(basename "${infile/.bam/.rm}").to-exclude.txt.gz" \
    "${outpath}/$(basename "${infile/.bam/.corrected.bam}")" \
    "${cluster}" && \
    touch "${step_8}"
elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
    echo_completion_message 8
    :
else
    echo_exit_message 8
    exit 1
fi


#  9: Run flagstat on bam in and outfiles -------------------------------------
if [[ ! -f "${step_9}" && -f "${step_8}" ]]; then
    run_flagstat 2 "${infile}"
    run_flagstat 2 "${infile/.bam/.rm.bam}"
    run_flagstat 2 "${outpath}/$(basename "${infile/.bam/.corrected.bam}")"

    echo_flagstat "${infile/.bam/.flagstat.txt}"
    echo_flagstat "${infile/.bam/.rm.flagstat.txt}"
    echo_flagstat "${outpath}/$(basename "${infile/.bam/.corrected.flagstat.txt}")"

    touch "${step_9}"
elif [[ -f "${step_9}" && -f "${step_8}" ]]; then
    echo_completion_message 9
    echo_flagstat "${infile/.bam/.flagstat.txt}"
    echo_flagstat "${infile/.bam/.rm.flagstat.txt}"
    echo_flagstat "${outpath}/$(basename "${infile/.bam/.corrected.flagstat.txt}")"
else
    echo_exit_message 9
    exit 1
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
run_time="$(echo "${time_end}" - "${time_start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "Run time: ${run_time} seconds."
echo ""

exit 0

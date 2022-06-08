#!/bin/bash

#  filter-qnames-by-assignment.sh
#  KA


time_start="$(date +%s)"


#  Source functions into environment ------------------------------------------
# shellcheck disable=1091
if [[ -f "./bin/auxiliary/functions-preprocessing-HPC.sh" ]]; then
    . ./bin/auxiliary/functions-preprocessing-HPC.sh ||
        {
            echo "Exiting: Unable to source 'functions-preprocessing-HPC.sh'."
            exit 1
        }

    . ./bin/auxiliary/functions-in-progress.sh ||
        {
            echo "Exiting: Unable to source 'functions-preprocessing-HPC.sh'."
            exit 1
        }
elif [[ -f "./functions-preprocessing-HPC.sh" ]]; then
    . ./functions-preprocessing-HPC.sh ||
        {
            echo "Exiting: Unable to source 'functions-preprocessing-HPC.sh'."
            exit 1
        }

    . ./functions-in-progress.sh ||
        {
            echo "Exiting: Unable to source 'functions-in-progress.sh'."
            exit 1
        }
fi

echo_completion_file_filter() {
    # #TODO Description of function
    #
    # :param 1: outpath (chr)
    # :param 2: step number (int)
    # :param 3: sample descriptor (chr)
    echo "${1}/$(basename "${0}" ".sh").${3}.step-${2}.txt"
}


#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Use sorted, gzipped lists of QNAMEs for alignment score-based"
    echo "assignments to generate assignment-specific bam files via the"
    echo "following steps:"
    echo "  - Step 01: Combine \"assign\" and \"unique\" files into \"combined\" files"
    echo "  - Step 02: Create a version of bam infile #1 for \"ambiguous\" assignments"
    echo "  - Step 03: Create a version of bam infile #2 for \"ambiguous\" assignments"
    echo "  - Step 04: Create a version of bam infile #1 for \"sample #1\" assignments"
    echo "  - Step 05: Create a version of bam infile #2 for \"sample #1\" assignments"
    echo "  - Step 06: Create a version of bam infile #1 for \"sample #2\" assignments"
    echo "  - Step 07: Create a version of bam infile #2 for \"sample #2\" assignments"
    echo "  - Step 08: Check that counts in bam outfiles are equal to bam infiles (optional)"
    echo "  - #TODO Step 09: Check that samtools flagstat readouts for bam infiles and outfiles are equivalent (optional)"
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - picard #TODO"
    echo "  - samtools #TODO"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-s use safe mode: \"TRUE\" or \"FALSE\" (logical)"
    echo "-l run on GS HPC: \"TRUE\" or \"FALSE\" (logical)"
    echo "-m initial memory allocation pool for JVM (chr; default \"512m\")"
    echo "-x maximum memory allocation pool for JVM (chr; default \"1g\")"
    echo "-b bam infile #1, including path (chr)"
    echo "-c bam infile #2, including path (chr)"
    echo "-a \"ambiguous\" assignment file, including path (chr)"
    echo "-i \"sample #1\" assignment file, including path (chr)"
    echo "-j \"sample #2\" assignment file, including path (chr)"
    echo "-u \"sample #1\" unique file, including path (chr)"
    echo "-v \"sample #2\" unique file, including path (chr)"
    echo "-1 string for \"sample #1\" (chr)"
    echo "-2 string for \"sample #2\" (chr)"
    echo "-o path for outfiles (chr); path will be made if it does not exist"
    echo "-p prefix for outfiles (chr)"
    echo "-n count lines: \"TRUE\" or \"FALSE\" (logical)"
    echo "#TODO -f run samtools flagstat: \"TRUE\" or \"FALSE\" (logical)"
    exit
}


while getopts "h:s:l:m:x:b:c:a:i:j:u:v:1:2:o:p:n:" opt; do
    case "${opt}" in
        h) print_usage ;;
        s) safe_mode="${OPTARG}" ;;
        l) cluster="${OPTARG}" ;;
        m) memory_min="${OPTARG}" ;;
        x) memory_max="${OPTARG}" ;;
        b) bam_1="${OPTARG}" ;;
        c) bam_2="${OPTARG}" ;;
        a) assign_ambiguous="${OPTARG}" ;;
        i) assign_sample_1="${OPTARG}" ;;
        j) assign_sample_2="${OPTARG}" ;;
        u) unique_sample_1="${OPTARG}" ;;
        v) unique_sample_2="${OPTARG}" ;;
        1) sample_1="${OPTARG}" ;;
        2) sample_2="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        p) prefix="${OPTARG}" ;;
        n) count="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && print_usage
[[ -z "${cluster}" ]] && print_usage
[[ -x "${memory_min}" ]] && memory_min="512m"
[[ -x "${memory_max}" ]] && memory_max="1g"
[[ -z "${bam_1}" ]] && print_usage
[[ -z "${bam_2}" ]] && print_usage
[[ -z "${assign_ambiguous}" ]] && print_usage
[[ -z "${assign_sample_1}" ]] && print_usage
[[ -z "${assign_sample_2}" ]] && print_usage
[[ -z "${unique_sample_1}" ]] && print_usage
[[ -z "${unique_sample_2}" ]] && print_usage
[[ -z "${sample_1}" ]] && print_usage
[[ -z "${sample_2}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${prefix}" ]] && print_usage
[[ -z "${count}" ]] && print_usage


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-u: \"Safe mode\" is TRUE." && set -Eeuxo pipefail ;;
    false | f) echo -e "-u: \"Safe mode\" is FALSE." ;;
    *) \
        echo -e "Exiting: -u \"safe mode\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Check that "${bam_1}" exists
[[ -f "${bam_1}" ]] ||
    {
        echo -e "Exiting: -b ${bam_1} does not exist.\n"
        exit 1
    }

#  Check that "${bam_2}" exists
[[ -f "${bam_2}" ]] ||
    {
        echo -e "Exiting: -c ${bam_2} does not exist.\n"
        exit 1
    }

#  Check that "${assign_ambiguous}" exists
[[ -f "${assign_ambiguous}" ]] ||
    {
        echo -e "Exiting: -a ${assign_ambiguous} does not exist.\n"
        exit 1
    }

#  Check that "${assign_sample_1}" exists
[[ -f "${assign_sample_1}" ]] ||
    {
        echo -e "Exiting: -i ${assign_sample_1} does not exist.\n"
        exit 1
    }

#  Check that "${assign_sample_2}" exists
[[ -f "${assign_sample_2}" ]] ||
    {
        echo -e "Exiting: -j ${assign_sample_2} does not exist.\n"
        exit 1
    }

#  Check that "${assign_sample_1}" exists
[[ -f "${unique_sample_1}" ]] ||
    {
        echo -e "Exiting: -u ${unique_sample_1} does not exist.\n"
        exit 1
    }

#  Check that "${assign_sample_2}" exists
[[ -f "${unique_sample_2}" ]] ||
    {
        echo -e "Exiting: -v ${unique_sample_2} does not exist.\n"
        exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo -e "-o: Directory ${outpath} does not exist; making the"
        echo -e "directory."
        mkdir -p "${outpath}"
    }

#  Evaluate "${count}"
case "$(echo "${count}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        count=1
        echo -e "-n: \"Count lines\" is TRUE."
        ;;
    false | f) \
        count=0
        echo -e "-n: \"Count lines\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -n \"count lines\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

# #  Evaluate "${flagstat}"
# case "$(echo "${flagstat}" | tr '[:upper:]' '[:lower:]')" in
#     true | t) \
#         flagstat=1
#         echo -e "-f: \"samtools flagstat\" is TRUE."
#         ;;
#     false | f) \
#         flagstat=0
#         echo -e "-f: \"samtools flagstat\" is FALSE."
#         ;;
#     *) \
#         echo -e "Exiting: -f \"samtools flagstat\" argument must be TRUE or FALSE.\n"
#         exit 1
#         ;;
# esac

echo ""


#  Assign variables needed for the pipeline -----------------------------------
combined_sample_1="${outpath}/${prefix}.${sample_1}.combined.txt.gz"
combined_sample_2="${outpath}/${prefix}.${sample_2}.combined.txt.gz"

bam_1_ambiguous="${outpath}/$(basename "${bam_1}" ".bam").ambiguous.bam"
bam_1_sample_1="${outpath}/$(basename "${bam_1}" ".bam").${sample_1}.bam"
bam_1_sample_2="${outpath}/$(basename "${bam_1}" ".bam").${sample_2}.bam"

bam_2_ambiguous="${outpath}/$(basename "${bam_2}" ".bam").ambiguous.bam"
bam_2_sample_1="${outpath}/$(basename "${bam_2}" ".bam").${sample_1}.bam"
bam_2_sample_2="${outpath}/$(basename "${bam_2}" ".bam").${sample_2}.bam"

#  Assign variables for completion files
step_1="$(echo_completion_file_filter "${outpath}" 1 "${prefix}")"
step_2="$(echo_completion_file_filter "${outpath}" 2 "${prefix}")"
step_3="$(echo_completion_file_filter "${outpath}" 3 "${prefix}")"
step_4="$(echo_completion_file_filter "${outpath}" 4 "${prefix}")"
step_5="$(echo_completion_file_filter "${outpath}" 5 "${prefix}")"
step_6="$(echo_completion_file_filter "${outpath}" 6 "${prefix}")"
step_7="$(echo_completion_file_filter "${outpath}" 7 "${prefix}")"
step_8="$(echo_completion_file_filter "${outpath}" 8 "${prefix}")"
# step_9="$(echo_completion_file_filter "${outpath}" 9 "${prefix}")"


#  01: Combine "assign" and "unique" files into "combined" files --------------
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started ${0} step 1/X: Combining \"assign\" and \"unique\"" \
    "files into \"combined\" files."

    cat "${assign_sample_1}" "${unique_sample_1}" > "${combined_sample_1}" && \
    cat "${assign_sample_2}" "${unique_sample_2}" > "${combined_sample_2}" && \
    touch "${step_1}"

    echo -e "Completed ${0} step 1/X: Combining \"assign\" and \"unique\"" \
    "files into \"combined\" files.\n"
elif [[ -f "${step_1}" ]]; then
    echo_completion_message 1
else
    echo_exit_message 1
    exit 1
fi


#  03: Create a version of bam infile #1 for "ambiguous" assignments ----------
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started ${0} step 2/X: Using ${assign_ambiguous} to filter" \
    "reads in ${bam_1}, thus creating ${bam_1_ambiguous}."

    if [[ $(LC_ALL=C gzip -l "${assign_ambiguous}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${assign_ambiguous} is empty."
        echo -e "Skipping ${0} step 2/X: Using ${assign_ambiguous} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_ambiguous}.\n"
    else
        retain_qname_reads_picard \
        "${bam_1}" \
        "${assign_ambiguous}" \
        "${bam_1_ambiguous}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_2}"

        echo -e "Completed ${0} step 2/X: Using ${assign_ambiguous} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_ambiguous}.\n"
    fi
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
else
    echo_exit_message 2
    exit 1
fi


#  03: Create a version of bam infile #2 for "ambiguous" assignments ----------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Started ${0} step 3/X: Using ${assign_ambiguous} to filter" \
    "reads in ${bam_2}, thus creating ${bam_2_ambiguous}."

    if [[ $(LC_ALL=C gzip -l "${assign_ambiguous}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${assign_ambiguous} is empty."
        echo -e "Skipping ${0} step 3/X: Using ${assign_ambiguous} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_ambiguous}.\n"
    else
        retain_qname_reads_picard \
        "${bam_2}" \
        "${assign_ambiguous}" \
        "${bam_2_ambiguous}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_3}"

        echo -e "Completed ${0} step 3/X: Using ${assign_ambiguous} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_ambiguous}.\n"
    fi
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
else
    echo_exit_message 3
    exit 1
fi


#  04: Create a version of bam infile #1 for "sample #1" assignments ----------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Started ${0} step 4/X: Using ${combined_sample_1} to filter" \
    "reads in ${bam_1}, thus creating ${bam_1_sample_1}."
    
    if [[ $(LC_ALL=C gzip -l "${combined_sample_1}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${combined_sample_1} is empty."
        echo -e "Skipping ${0} step 4/X: Using ${combined_sample_1} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_sample_1}.\n"
    else
        retain_qname_reads_picard \
        "${bam_1}" \
        "${combined_sample_1}" \
        "${bam_1_sample_1}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_4}"

        echo -e "Completed ${0} step 4/X: Using ${combined_sample_1} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_sample_1}.\n"
    fi
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
else
    echo_exit_message 4
    exit 1
fi


#  05: Create a version of bam infile #2 for "sample #1" assignments ----------
if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    echo -e "Started ${0} step 5/X: Using ${combined_sample_1} to filter" \
    "reads in ${bam_2}, thus creating ${bam_2_sample_1}."

    if [[ $(LC_ALL=C gzip -l "${combined_sample_1}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${combined_sample_1} is empty."
        echo -e "Skipping ${0} step 5/X: Using ${combined_sample_1} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_sample_1}.\n"
    else
        retain_qname_reads_picard \
        "${bam_2}" \
        "${combined_sample_1}" \
        "${bam_2_sample_1}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_5}"

        echo -e "Completed ${0} step 5/X: Using ${combined_sample_1} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_sample_1}.\n"
    fi
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    echo_completion_message 5
else
    echo_exit_message 5
    exit 1
fi


#  06: Create a version of bam infile #1 for "sample #2" assignments ----------
if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then
    echo -e "Started ${0} step 6/X: Using ${combined_sample_2} to filter" \
    "reads in ${bam_1}, thus creating ${bam_1_sample_2}."

    if [[ $(LC_ALL=C gzip -l "${combined_sample_2}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${combined_sample_2} is empty."
        echo -e "Skipping ${0} step 6/X: Using ${combined_sample_2} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_sample_2}.\n"
    else
        retain_qname_reads_picard \
        "${bam_1}" \
        "${combined_sample_2}" \
        "${bam_1_sample_2}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_6}"

        echo -e "Completed ${0} step 6/X: Using ${combined_sample_2} to filter" \
        "reads in ${bam_1}, thus creating ${bam_1_sample_2}.\n"
    fi
elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
    echo_completion_message 6
else
    echo_exit_message 6
    exit 1
fi


#  07: Create a version of bam infile #2 for "sample #2" assignments ----------
if [[ ! -f "${step_7}" && -f "${step_6}" ]]; then
    echo -e "Started ${0} step 7/X: Using ${combined_sample_2} to filter" \
    "reads in ${bam_2}, thus creating ${bam_2_sample_2}."

    if [[ $(LC_ALL=C gzip -l "${combined_sample_2}" | awk 'NR==2 {exit($2!=0)}') ]]; then
        echo -e "WARNING: ${combined_sample_2} is empty."
        echo -e "Skipping ${0} step 7/X: Using ${combined_sample_2} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_sample_2}.\n"
    else
        retain_qname_reads_picard \
        "${bam_2}" \
        "${combined_sample_2}" \
        "${bam_2_sample_2}" \
        "${memory_min}" \
        "${memory_max}" \
        "${cluster}" && \
        touch "${step_7}"

        echo -e "Completed ${0} step 7/X: Using ${combined_sample_2} to filter" \
        "reads in ${bam_2}, thus creating ${bam_2_sample_2}.\n"
    fi
elif [[ -f "${step_7}" && -f "${step_6}" ]]; then
    echo_completion_message 7
else
    echo_exit_message 7
    exit 1
fi


#  08: Check that outfile counts equal infile counts (optional) ---------------
if [[ "${count}" == 1 ]]; then
    n_bam_1="$(samtools view -c "${bam_1}")"
    n_bam_1_ambiguous="$(samtools view -c "${bam_1_ambiguous}")"
    n_bam_1_sample_1="$(samtools view -c "${bam_1_sample_1}")"
    n_bam_1_sample_2="$(samtools view -c "${bam_1_sample_2}")"
    sum_bam_1_outfiles=$(( n_bam_1_ambiguous + n_bam_1_sample_1 + n_bam_1_sample_2 ))

    n_bam_2="$(samtools view -c "${bam_2}")"
    n_bam_2_ambiguous="$(samtools view -c "${bam_2_ambiguous}")"
    n_bam_2_sample_1="$(samtools view -c "${bam_2_sample_1}")"
    n_bam_2_sample_2="$(samtools view -c "${bam_2_sample_2}")"
    sum_bam_2_outfiles=$(( n_bam_2_ambiguous + n_bam_2_sample_1 + n_bam_2_sample_2 ))
    if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
        echo -e "Started ${0} step 8/X: Check that counts in bam outfiles" \
        "are equal to bam infiles."

        echo "  - ${bam_1}: $(( n_bam_1 ))" && \
        echo "  - ${bam_1_ambiguous}: $(( n_bam_1_ambiguous ))" && \
        echo "  - ${bam_1_sample_1}: $(( n_bam_1_sample_1 ))" && \
        echo "  - ${bam_1_sample_2}: $(( n_bam_1_sample_2 ))" && \
        echo "  - Sum of ${bam_1} outfiles: $(( sum_bam_1_outfiles ))" && \
        echo "" && \
        echo "  - ${bam_2}: $(( n_bam_2 ))" && \
        echo "  - ${bam_2_ambiguous}: $(( n_bam_2_ambiguous ))" && \
        echo "  - ${bam_2_sample_1}: $(( n_bam_2_sample_1 ))" && \
        echo "  - ${bam_2_sample_2}: $(( n_bam_2_sample_2 ))" && \
        echo "  - Sum of ${bam_2} outfiles: $(( sum_bam_2_outfiles ))" && \
        echo "" && \
        touch "${step_8}"

        echo -e "Completed ${0} step 8/X: Check that counts in bam outfiles" \
        "are equal to bam infiles.\n"
    elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
        echo -e "Started ${0} step 8/X: Check that counts in bam outfiles" \
        "are equal to bam infiles."

        echo "  - ${bam_1}: $(( n_bam_1 ))" && \
        echo "  - ${bam_1_ambiguous}: $(( n_bam_1_ambiguous ))" && \
        echo "  - ${bam_1_sample_1}: $(( n_bam_1_sample_1 ))" && \
        echo "  - ${bam_1_sample_2}: $(( n_bam_1_sample_2 ))" && \
        echo "  - ${bam_1} outfiles: $(( sum_bam_1_outfiles ))" && \
        echo "" && \
        echo "  - ${bam_2}: $(( n_bam_2 ))" && \
        echo "  - ${bam_2_ambiguous}: $(( n_bam_2_ambiguous ))" && \
        echo "  - ${bam_2_sample_1}: $(( n_bam_2_sample_1 ))" && \
        echo "  - ${bam_2_sample_2}: $(( n_bam_2_sample_2 ))" && \
        echo "  - ${bam_2} outfiles: $(( sum_bam_2_outfiles ))" && \
        echo ""

        echo -e "Completed ${0} step 8/X: Check that counts in bam outfiles" \
        "are equal to bam infiles.\n"
    else
        echo_exit_message 8
        exit 1
    fi
fi


# #  09: Check that samtools flagstat readouts are equivalent (optional) --------
# if [[ "${flagstat}" == 1 ]]; then
#     #TODO Insert code here
# fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

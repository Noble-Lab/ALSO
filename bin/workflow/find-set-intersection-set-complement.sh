#!/bin/bash

#  find-set-intersection-set-complement.sh
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

echo_completion_file_sets() {
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
    echo "Filter duplicate QNAMEs from AS.txt.gz file via the following steps:"
    echo "  - Step 01: Sort AS.txt.gz files"
    echo "  - Step 02: Find set intersection between sample #1 and sample #2"
    echo "  - Step 03: Find set complement for \"sample #1\""
    echo "  - Step 04: Find set complement for \"sample #2\""
    echo "  - Step 05: Check that intersection and complement sums to set size (optional)"
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - GNU zcat >= #TODO (#TODO I think >= 1.09 but need to check...)"
    echo "  - GNU sort >= #TODO (#TODO Is BSD sort OK?)"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-u use safe mode: TRUE or FALSE [logical; default: FALSE]"
    echo "-i AS.txt.gz for \"sample #1\", including path [chr]"
    echo "-j AS.txt.gz for \"sample #2\", including path [chr]"
    echo "-1 string for \"sample #1\" [chr]"
    echo "-2 string for \"sample #2\" [chr]"
    echo "-o path for outfiles [chr]; path will be made if it does not exist"
    echo "-p prefix for outfiles [chr]"
    echo "-c count lines: TRUE or FALSE [logical]"
    exit
}


while getopts "h:u:i:j:1:2:p:o:c:" opt; do
    case "${opt}" in
        h) print_usage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile_1="${OPTARG}" ;;
        j) infile_2="${OPTARG}" ;;
        1) sample_1="${OPTARG}" ;;
        2) sample_2="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        p) prefix="${OPTARG}" ;;
        c) count="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${infile_1}" ]] && print_usage
[[ -z "${infile_2}" ]] && print_usage
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

#  Check that "${infile_1}" exists
[[ -f "${infile_1}" ]] ||
    {
        echo -e "Exiting: -i ${infile_1} does not exist.\n"
        exit 1
    }

#  Check that "${infile_2}" exists
[[ -f "${infile_2}" ]] ||
    {
        echo -e "Exiting: -j ${infile_2} does not exist.\n"
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
        echo -e "-c: \"Count lines\" is TRUE."
        ;;
    false | f) \
        count=0
        echo -e "-c: \"Count lines\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -c \"count lines\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

echo ""


#  Assign variables needed for the pipeline -----------------------------------
inter="${outpath}/${prefix}.${sample_1}-${sample_2}.intersection.txt.gz"
comp_1="${outpath}/${prefix}.${sample_1}.complement.txt.gz"
comp_2="${outpath}/${prefix}.${sample_2}.complement.txt.gz"


#  Assign variables for completion files
step_1="$(echo_completion_file_sets "${outpath}" 1 "${prefix}")"
step_2="$(echo_completion_file_sets "${outpath}" 2 "${prefix}")"
step_3="$(echo_completion_file_sets "${outpath}" 3 "${prefix}")"
step_4="$(echo_completion_file_sets "${outpath}" 4 "${prefix}")"
step_5="$(echo_completion_file_sets "${outpath}" 4 "${prefix}")"


#  01: Sort AS.txt.gz files ---------------------------------------------------
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started step 1/5: Sorting ${infile_1} and" "${infile_2}."
    
    # sort_file_AS_overwrite_infile "${infile_1}" && \
    # sort_file_AS_overwrite_infile "${infile_2}" && \
    # touch "${step_1}"
    presorted_1="${infile_1%.txt.gz}.pre-sorted.txt.gz"
    presorted_2="${infile_2%.txt.gz}.pre-sorted.txt.gz"
    sorted_1="${infile_1%.txt.gz}.sorted.txt.gz"
    sorted_2="${infile_2%.txt.gz}.sorted.txt.gz"

    cp "${infile_1}" "${presorted_1}" && \
    gzip -dk "${infile_1}" && \
    sort -k1,1 -k2n "${infile_1%.gz}" | gzip > "${sorted_1}" && \
    mv -f "${sorted_1}" "${infile_1}" && \
    cp "${infile_2}" "${presorted_2}" && \
    gzip -dk "${infile_2}" && \
    sort -k1,1 -k2n "${infile_2%.gz}" | gzip > "${sorted_2}" && \
    mv -f "${sorted_2}" "${infile_2}" && \
    touch "${step_1}"

    echo -e "Completed step 1/5: Sorting ${infile_1} and" "${infile_2}."
else
    echo_completion_message 1
fi


#  02: Find set intersection between sample #1 and sample #2 ------------------
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started step 2/5: Finding set intersection between ${infile_1}" \
    "and ${infile_2}."

    find_set_intersection "${infile_1}" "${infile_2}" "${inter}" && \
    touch "${step_2}"

    echo -e "Completed step 2/5: Finding set intersection between ${infile_1}" \
    "and ${infile_2}.\n"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
else
    echo_completion_message 2
fi


#  03: Find set complement for sample #1 --------------------------------------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Started step 3/5: Between ${infile_1} and ${infile_2}, finding" \
    "set complement for ${infile_1}."

    find_set_complement "${infile_1}" "${infile_2}" "${comp_1}" && \
    touch "${step_3}"

    echo -e "Completed step 3/5: Between ${infile_1} and ${infile_2}, finding" \
    "set complement for ${infile_1}.\n"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
else
    echo_exit_message 3
    exit 1
fi


#  04: Find set complement for sample #2 --------------------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Started step 4/5: Between ${infile_1} and ${infile_2}, finding set" \
    "complement for ${infile_2}."

    find_set_complement "${infile_2}" "${infile_1}" "${comp_2}" && \
    touch "${step_4}"

    echo -e "Completed step 4/5: Between ${infile_1} and ${infile_2}, finding set" \
    "complement for ${infile_2}.\n"
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
else
    echo_exit_message 4
    exit 1
fi


#  05: Check that intersection and complement sums to set size (optional) -----
if [[ "${count}" == 1 ]]; then
    if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
        n_inter=$(zcat "${inter}" | wc -l)
        n_comp_1=$(zcat "${comp_1}" | wc -l)
        n_comp_2=$(zcat "${comp_2}" | wc -l)
        n_1=$(zcat "${infile_1}" | wc -l)
        n_2=$(zcat "${infile_2}" | wc -l)

        echo -e "Started step 5/5: Check that the set intersection and set" \
        "complement sums to the set size."

        echo "Sum of ${inter} and ${comp_1}: $(( n_inter + n_comp_1 ))" && \
        echo "Count for ${infile_1}: $(( n_1 ))" && \
        echo "" && \
        echo "Sum of ${inter} and ${comp_2}: $(( n_inter + n_comp_2 ))" && \
        echo "Count for ${infile_2}: $(( n_2 ))" && \
        touch "${step_5}"

        echo -e "Completed step 5/5: Check that the set intersection and set" \
        "complement sums to the set size.\n"
    elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
        n_inter=$(zcat "${inter}" | wc -l)
        n_comp_1=$(zcat "${comp_1}" | wc -l)
        n_comp_2=$(zcat "${comp_2}" | wc -l)
        n_1=$(zcat "${infile_1}" | wc -l)
        n_2=$(zcat "${infile_2}" | wc -l)

        echo_completion_message 5
        echo "Sum of ${inter} and ${comp_1}: $(( n_inter + n_comp_1 ))" && \
        echo "Count for ${infile_1}: $(( n_1 ))" && \
        echo "" && \
        echo "Sum of ${inter} and ${comp_2}: $(( n_inter + n_comp_2 ))" && \
        echo "Count for ${infile_2}: $(( n_2 ))"
    else
        echo_exit_message 5
        exit 1
    fi
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

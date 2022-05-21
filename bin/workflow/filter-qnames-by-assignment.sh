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


#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Run pipeline to generate bam files for assignments."
    echo "  - Step 01: Copy files of interest to \${TMPDIR}"
    echo "  - Step 02: ..."
    echo "  - Step 03: ..."
    echo "  - Step 04: ..."
    echo "  - Step 05: ..."
    echo "  - Step 06: ..."
    echo "  - Step 07: ..."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - ..."
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-s use safe mode: \"TRUE\" or \"FALSE\" (logical)"
    echo "-c use KA conda environment: \"TRUE\" or \"FALSE\" (logical)"
    echo "-l run on GS HPC: \"TRUE\" or \"FALSE\" (logical)"
    echo "-b bam infile, including path (chr)"
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
    exit
}


while getopts "h:s:c:l:b:a:i:j:u:v:1:2:o:p:n:" opt; do
    case "${opt}" in
        h) print_usage ;;
        s) safe_mode="${OPTARG}" ;;
        c) conda="${OPTARG}" ;;
        l) cluster="${OPTARG}" ;;
        b) bam="${OPTARG}" ;;
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
[[ -z "${conda}" ]] && print_usage
[[ -z "${cluster}" ]] && print_usage
[[ -z "${bam}" ]] && print_usage
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

#  Check that "${bam}" exists
[[ -f "${bam}" ]] ||
    {
        echo -e "Exiting: -b ${bam} does not exist.\n"
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
base_bam="$(basename "${bam}")"
base_assign_ambiguous="$(basename "${assign_ambiguous}")"
base_assign_sample_1="$(basename "${assign_sample_1}")"
base_assign_sample_2="$(basename "${assign_sample_2}")"
base_unique_sample_1="$(basename "${unique_sample_1}")"
base_unique_sample_2="$(basename "${unique_sample_2}")"
base_combined_sample_1="${prefix}.${sample_1}.combined.txt.gz"
base_combined_sample_2="${prefix}.${sample_2}.combined.txt.gz"
base_bam_ambiguous="${base_bam%.bam}.ambiguous.bam"
base_bam_sample_1="${base_bam%.bam}.${sample_1}.bam"
base_bam_sample_2="${base_bam%.bam}.${sample_2}.bam"

tmp_bam="${TMPDIR}/${base_bam}"
tmp_assign_ambiguous="${TMPDIR}/${base_assign_ambiguous}"
tmp_assign_sample_1="${TMPDIR}/${base_assign_sample_1}"
tmp_assign_sample_2="${TMPDIR}/${base_assign_sample_2}"
tmp_unique_sample_1="${TMPDIR}/${base_unique_sample_1}"
tmp_unique_sample_2="${TMPDIR}/${base_unique_sample_2}"
tmp_combined_sample_1="${TMPDIR}/${base_combined_sample_1}"
tmp_combined_sample_2="${TMPDIR}/${base_combined_sample_2}"
tmp_bam_ambiguous="${TMPDIR}/${base_bam_ambiguous}"
tmp_bam_sample_1="${TMPDIR}/${base_bam_sample_1}"
tmp_bam_sample_2="${TMPDIR}/${base_bam_sample_2}"

#  Assign variables for completion files
step_1="$(echo_completion_file "${outpath}" 1 "${prefix}")"
step_2="$(echo_completion_file "${outpath}" 2 "${prefix}")"
step_3="$(echo_completion_file "${outpath}" 3 "${prefix}")"
step_4="$(echo_completion_file "${outpath}" 4 "${prefix}")"
step_5="$(echo_completion_file "${outpath}" 5 "${prefix}")"
step_6="$(echo_completion_file "${outpath}" 6 "${prefix}")"


#  01: Copy files of interest to ${TMPDIR} ------------------------------------
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started step 1/X: Copying pertinent files into ${TMPDIR}."

    cp "${bam}" \
    "${assign_ambiguous}" \
    "${assign_sample_1}" \
    "${assign_sample_2}" \
    "${unique_sample_1}" \
    "${unique_sample_2}" \
    "${TMPDIR}" && \
    touch "${step_1}"

    echo -e "Completed step 1/X: Copying pertinent files into ${TMPDIR}.\n"
else
    echo_completion_message 1
fi


#  02: Combine "assign" and "unique" files into "combined" files --------------
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started step 2/X: Combining \"assign\" and \"unique\" files"
    echo -e "into \"combined\" files."

    cat "${tmp_assign_sample_1}" "${tmp_unique_sample_1}" \
    > "${tmp_combined_sample_1}" && \
    cat "${tmp_assign_sample_2}" "${tmp_unique_sample_2}" \
    > "${tmp_combined_sample_2}" && \
    touch "${step_2}"

    echo -e "Completed step 2/X: Combining \"assign\" and \"unique\" files"
    echo -e "into \"combined\" files.\n"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
else
    echo_exit_message 2
    exit 1
fi


#  03: Create a bam file for "ambiguous" assignments --------------------------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Started step 3/X: Using ${tmp_assign_ambiguous} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_ambiguous}."

    include_qname_reads_picard \
    "${tmp_bam}" \
    "${tmp_assign_ambiguous}" \
    "${tmp_bam_ambiguous}" \
    "${cluster}" && \
    touch "${step_3}"

    echo -e "Completed step 3/X: Using ${tmp_assign_ambiguous} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_ambiguous}.\n"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
else
    echo_exit_message 3
    exit 1
fi


#  04: Create a bam file for "sample #1" assignments --------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Started step 4/X: Using ${tmp_combined_sample_1} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_sample_1}."

    include_qname_reads_picard \
    "${tmp_bam}" \
    "${tmp_combined_sample_1}" \
    "${tmp_bam_sample_1}" \
    "${cluster}" && \
    touch "${step_4}"

    echo -e "Completed step 4/X: Using ${tmp_combined_sample_1} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_sample_1}.\n"
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
else
    echo_exit_message 4
    exit 1
fi


#  05: Create a bam file for "sample #2" assignments --------------------------
if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    echo -e "Started step 5/X: Using ${tmp_combined_sample_2} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_sample_2}."

    include_qname_reads_picard \
    "${tmp_bam}" \
    "${tmp_combined_sample_2}" \
    "${tmp_bam_sample_2}" \
    "${cluster}" && \
    touch "${step_5}"

    echo -e "Completed step 5/X: Using ${tmp_combined_sample_2} to filter"
    echo -e "reads in ${tmp_bam}, thus creating ${tmp_bam_sample_2}.\n"
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    echo_completion_message 5
else
    echo_exit_message 5
    exit 1
fi


# #  Optional step 06 -----------------------------------------------------------
# if [[ "${count}" == 1 ]]; then
#     #  06: Check that intersection and complement sums to set size ------------
#     if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then
#         n_inter=$(zcat "${tmp_inter}" | wc -l)
#         n_comp_1=$(zcat "${tmp_comp_1}" | wc -l)
#         n_comp_2=$(zcat "${tmp_comp_2}" | wc -l)
#         n_1=$(zcat "${tmp_1}" | wc -l)
#         n_2=$(zcat "${tmp_2}" | wc -l)
#
#         echo -e "Started step 6/X: Check that the set intersection and set"
#         echo -e "complement sums to the set size."
#
#         echo "Sum of ${base_inter} and ${base_comp_1}: $(( n_inter + n_comp_1 ))" && \
#         echo "Count for ${tmp_1}: $(( n_1 ))" && \
#         echo "" && \
#         echo "Sum of ${base_inter} and ${base_comp_2}: $(( n_inter + n_comp_2 ))" && \
#         echo "Count for ${tmp_2}: $(( n_2 ))" && \
#         touch "${step_6}"
#
#         echo -e "Completed step 6/X: Check that the set intersection and set"
#         echo -e "complement sums to the set size.\n"
#     elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
#         n_inter=$(zcat "${tmp_inter}" | wc -l)
#         n_comp_1=$(zcat "${tmp_comp_1}" | wc -l)
#         n_comp_2=$(zcat "${tmp_comp_2}" | wc -l)
#         n_1=$(zcat "${tmp_1}" | wc -l)
#         n_2=$(zcat "${tmp_2}" | wc -l)
#
#         echo_completion_message 6
#         echo "Sum of ${base_inter} and ${base_comp_1}: $(( n_inter + n_comp_1 ))" && \
#         echo "Count for ${tmp_1}: $(( n_1 ))" && \
#         echo "" && \
#         echo "Sum of ${base_inter} and ${base_comp_2}: $(( n_inter + n_comp_2 ))" && \
#         echo "Count for ${tmp_2}: $(( n_2 ))"
#     else
#         echo_exit_message 6
#         exit 1
#     fi
# fi


#  0X: Remove temporary files, move "${TMPDIR}" outfiles to "${outpath}" ------
if [[ -f "${step_5}" && -f "${tmp_bam_sample_2}" ]]; then
    echo -e "Started step X/X: Removing temporary files, moving outfiles"
    echo -e "from ${TMPDIR} to ${outpath}."

    rm \
    "${tmp_bam}" "${tmp_assign_ambiguous}" \
    "${tmp_assign_sample_1}" "${tmp_assign_sample_2}" \
    "${tmp_unique_sample_1}" "${tmp_unique_sample_2}" && \
    mv -f "${TMPDIR}/"*.{txt.gz,bam} "${outpath}" && \
    touch "${step_6}"

    echo -e "Completed step X/X: Removing temporary files, moving outfiles"
    echo -e "from ${TMPDIR} to ${outpath}.\n"
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0


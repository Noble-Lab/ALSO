#!/bin/bash

#  remove-duplicate-qnames.sh
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
else
    echo -e "Exiting: Could not find auxiliary information."
    exit 1
fi



#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Run pipeline to filter duplicate QNAMEs from bam file."
    echo "  - Step 01: Copy files of interest to \${TMPDIR}"
    echo "  - Step 02: Sort bam by QNAME"
    echo "  - Step 03: List and tally QNAMEs in the sorted bam file"
    echo "  - Step 04: Create txt.gz outfiles for QNAME > 2"
    echo "  - Step 05: Count lines in infile, outfiles (optional)"
    echo "  - Step 06: Tally entries in infile, outfiles (optional)"
    echo "  - Step 07: Exclude problematic QNAME reads from bam infile"
    echo "  - Step 08: Sort corrected bam by QNAME (optional)"
    echo "  - Step 09: List and tally QNAMEs in the corrected bam file"
    echo "             (optional)"
    echo "  - Step 10: Create txt.gz outfiles for QNAME >, <, = 2 (optional)"
    echo "  - Step 11: Remove temporary bams, move \${TMPDIR} outfiles to"
    echo "             \${outpath}"
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - parallel >= 20200101"
    echo "  - picard >= 2.27.1"
    echo "  - samtools >= 1.13"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-u use safe mode: \"TRUE\" or \"FALSE\" (logical)"
    echo "-c run on GS HPC: \"TRUE\" or \"FALSE\" (logical)"
    echo "-m initial memory allocation pool for JVM (chr; default \"512m\")"
    echo "-x maximum memory allocation pool for JVM (chr; default \"1g\")"
    echo "-i bam infile, including path (chr)"
    echo "-o path for outfiles (chr); path will be made if it does not exist"
    echo "-n count lines: \"TRUE\" or \"FALSE\" (logical)"
    echo "-t tally entries: \"TRUE\" or \"FALSE\" (logical)"
    echo "-e evaluate corrected bam: \"TRUE\" or \"FALSE\" (logical)"
    echo "-r remove intermediate files: \"TRUE\" or \"FALSE\" (logical)"
    echo "-p number of cores for parallelization (int >= 1; default: 1)"
    exit
}


while getopts "h:u:c:m:x:i:o:n:t:e:r:p:" opt; do
    case "${opt}" in
        h) print_usage ;;
        u) safe_mode="${OPTARG}" ;;
        c) cluster="${OPTARG}" ;;
        m) memory_min="${OPTARG}" ;;
        x) memory_max="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        n) count="${OPTARG}" ;;
        t) tally="${OPTARG}" ;;
        e) evaluate="${OPTARG}" ;;
        r) remove="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && print_usage
[[ -z "${cluster}" ]] && print_usage
[[ -x "${memory_min}" ]] && memory_min="512m"
[[ -x "${memory_max}" ]] && memory_max="1g"
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${count}" ]] && print_usage
[[ -z "${tally}" ]] && print_usage
[[ -z "${evaluate}" ]] && print_usage
[[ -z "${remove}" ]] && print_usage
[[ -z "${parallelize}" ]] && parallelize=1


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

#  Check for necessary dependencies; exit if not found
check_dependency samtools

#  Evaluate "${cluster}"
case "$(echo "${cluster}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-c: \"Run on GS HPC\" is TRUE." ;;
    false | f) echo -e "-c: \"Run on GS HPC\" is FALSE." && check_dependency picard ;;
    *) \
        echo -e "Exiting: -c argument must be \"TRUE\" or \"FALSE\".\n"
        return 1
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
        echo -e "-o: Directory ${outpath} does not exist; making the directory."
        mkdir -p "${outpath}"
    }

#  Evaluate "${tally}"
case "$(echo "${tally}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        tally=1
        echo -e "-t: \"Tally entries\" is TRUE."
        ;;
    false | f) \
        tally=0
        echo -e "-t: \"Tally entries\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -t \"tally entries\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

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

#  Evaluate "${evaluate}"
case "$(echo "${evaluate}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        evaluate=1
        echo -e "-e: \"Evaluate corrected bam\" is TRUE."
        ;;
    false | f) \
        evaluate=0
        echo -e "-e: \"Evaluate corrected bam\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -e \"evaluate corrected bam\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Evaluate "${remove}"
case "$(echo "${remove}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        remove=1
        echo -e "-r: \"Remove intermediate files\" is TRUE."
        ;;
    false | f) \
        remove=0
        echo -e "-r: \"Remove intermediate files\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -r \"remove intermediate files\" must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Check "${parallelize}"
[[ ! "${parallelize}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -p \"parallelize\" argument must be an integer.\n"
        exit 1
    }

[[ ! $((parallelize)) -ge 1 ]] &&
    {
        echo -e "Exiting: -p \"parallelize\" argument must be an integer >= 1.\n"
        exit 1
    }

echo ""


#  Assign variables needed for the pipeline -----------------------------------
# infile="Disteche_sample_1.dedup.CAST.corrected.CAST.AS.txt.gz"
# outpath="$(pwd)"
# tmp="${outpath}/${infile}"

base="$(basename "${infile}")"
tmp="${TMPDIR}/${base}"
tmp_corrected="${TMPDIR}/${base%.bam}.corrected.bam"
tmp_corrected_sorted="${tmp_corrected%.bam}.sort-n.bam"
tmp_sorted="${TMPDIR}/${base%.bam}.sort-n.bam"
tmp_sorted_QNAME_full="${tmp_sorted%.bam}.QNAME.txt.gz"
tmp_sorted_QNAME_gt="${tmp_sorted%.bam}.QNAME.gt.txt.gz"
tmp_sorted_QNAME_full_tally="${tmp_sorted_QNAME_full%.txt.gz}.tally.txt.gz"
tmp_sorted_QNAME_gt_tally="${tmp_sorted_QNAME_gt%.txt.gz}.tally.txt.gz"

#  Assign variables for completion files
step_1="$(echo_completion_file "${outpath}" 1 "${infile%.gz}")"
step_2="$(echo_completion_file "${outpath}" 2 "${infile%.gz}")"
step_3="$(echo_completion_file "${outpath}" 3 "${infile%.gz}")"
step_4="$(echo_completion_file "${outpath}" 4 "${infile%.gz}")"
step_5="$(echo_completion_file "${outpath}" 5 "${infile%.gz}")"
step_6="$(echo_completion_file "${outpath}" 6 "${infile%.gz}")"
step_7="$(echo_completion_file "${outpath}" 7 "${infile%.gz}")"
step_8="$(echo_completion_file "${outpath}" 8 "${infile%.gz}")"
step_9="$(echo_completion_file "${outpath}" 9 "${infile%.gz}")"
step_10="$(echo_completion_file "${outpath}" 10 "${infile%.gz}")"
step_11="$(echo_completion_file "${outpath}" 11 "${infile%.gz}")"


#  01: Copy files of interest to ${TMPDIR} ------------------------------------
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started step 1/11: Copying ${base} into ${TMPDIR}."

    cp "${infile}" "${TMPDIR}" && \
    touch "${step_1}"

    echo -e "Completed step 1/11: Copying ${base} into ${TMPDIR}.\n"
else
    echo_completion_message 1
fi


#  02: Sort bam by QNAME ------------------------------------------------------
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started step 2/11: Sorting ${tmp} by QNAME."

    sort_bam_qname_samtools "${parallelize}" \
    "${tmp}" "${tmp_sorted}" && \
    touch "${step_2}"

    echo -e "Completed step 2/11: Sorting ${tmp} by QNAME.\n"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
else
    echo_exit_message 2
    exit 1
fi


#  03: List and tally QNAMEs in the sorted bam file ---------------------------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Started step 3/11: Listing, tallying QNAMEs in ${tmp_sorted}."

    list_tally_qnames_gzip_updated "${tmp_sorted}" && \
    touch "${step_3}"

    echo -e "Completed step 3/11: Listing, tallying QNAMEs in ${tmp_sorted}.\n"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
else
    echo_exit_message 3
    exit 1
fi


#  04: Create txt.gz outfiles for QNAME > 2 -----------------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Started step 4/11: Create txt.gz outfiles for QNAME > 2 from ${tmp_sorted%.bam}.QNAME.txt.gz."

    identify_qnames_updated "${parallelize}" \
    "gt" "${tmp_sorted_QNAME_full}" "gzip" && \
    touch "${step_4}"

    echo -e "Completed step 4/11: Create txt.gz outfiles for QNAME > 2 from ${tmp_sorted_QNAME_full}.\n"
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
else
    echo_exit_message 4
    exit 1
fi


#  Optional steps 05, 06 ------------------------------------------------------
if [[ "${count}" == 1 ]]; then
    #  05: Count lines in infile, outfiles (optional) -------------------------
    if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
        count_full="$(count_lines_gzip "${tmp_sorted_QNAME_full}")"
        count_remove="$(count_lines_gzip "${tmp_sorted_QNAME_gt}")"

        echo -e "Started step 5/11: Counting entries in files, writing out lists."

        echo "$(basename "${tmp_sorted_QNAME_full}"): ${count_full}" && \
        echo "$(basename "${tmp_sorted_QNAME_gt}"): ${count_remove}" && \
        echo "" && \
        touch "${step_5}"

        echo -e "Completed step 5/11: Counting entries in files, writing out lists.\n"
    elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
        count_full="$(count_lines_gzip "${tmp_sorted_QNAME_full}")"
        count_remove="$(count_lines_gzip "${tmp_sorted_QNAME_gt}")"

        echo_completion_message 5

        echo "$(basename "${tmp_sorted_QNAME_full}"): ${count_full}" && \
        echo "$(basename "${tmp_sorted_QNAME_gt}"): ${count_remove}" && \
        echo ""
    else
        echo_exit_message 5
        exit 1
    fi
fi

if [[ "${tally}" == 1 ]]; then
    #  06: Tally entries in infile, outfiles (optional) -----------------------
    if [[ ! -f "${step_6}" && -f "${step_4}" ]]; then
        echo -e "Started step 6/11: Tallying entries in"
        echo -e "${tmp_sorted_QNAME_full} and ${tmp_sorted_QNAME_gt}, writing out lists."

        tally_qnames_gzip \
        "${tmp_sorted_QNAME_full}" \
        "${tmp_sorted_QNAME_full_tally}" && \
        tally_qnames_gzip \
        "${tmp_sorted_QNAME_gt}" \
        "${tmp_sorted_QNAME_gt_tally}" && \
        touch "${step_6}"

        echo -e "Completed step 6/11: Tallying entries in"
        echo -e "${tmp_sorted_QNAME_full} and ${tmp_sorted_QNAME_gt}, writing out lists.\n"
    elif [[ -f "${step_6}" && -f "${step_4}" ]]; then
        echo_completion_message 6
    else
        echo_exit_message 6
        exit 1
    fi
fi


#  07: Exclude problematic QNAME reads from bam infile ------------------------
if [[ ! -f "${step_7}" && -f "${step_4}" ]]; then
    echo -e "Started step 7/11: Using ${tmp_sorted_QNAME_gt} to exclude reads from ${tmp}."

    exclude_qname_reads_picard \
    "${tmp}" \
    "${tmp_sorted_QNAME_gt}" \
    "${tmp_corrected}" \
    "${memory_min}" \
    "${memory_max}" \
    "${cluster}" && \
    touch "${step_7}"
    # :param 1: name of bam infile, including path (chr)
    # :param 2: name of txt QNAME list, including path (chr)
    # :param 3: name of bam outfile, including path (cannot be same as bam
    #           infile) (chr)
    # :param 4: initial memory allocation pool for JVM (chr)
    # :param 5: maximum memory allocation pool for JVM (chr)
    # :param 6: use the picard.jar available on the GS grid system (logical)

    echo -e "Completed step 7/11: Using ${tmp_sorted_QNAME_gt} to exclude reads from ${tmp}.\n"
elif [[ -f "${step_7}" && -f "${step_4}" ]]; then
    echo_completion_message 7
else
    echo_exit_message 7
    exit 1
fi


#  Optional steps 08-10 -------------------------------------------------------
if [[ "${evaluate}" == 1 ]]; then
    #  08: Sort corrected bam by QNAME ----------------------------------------
    if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
        echo -e "Started step 8/11: Sorting ${tmp} by QNAME."

        sort_bam_qname_samtools \
        "${parallelize}" \
        "${tmp_corrected}" \
        "${tmp_corrected_sorted}" && \
        touch "${step_8}"

        echo -e "Completed step 8/11: Sorting ${tmp} by QNAME.\n"
    elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
        echo_completion_message 8
    else
        echo_exit_message 8
        exit 1
    fi

    #  09: List and tally QNAMEs in the corrected bam file --------------------
    if [[ ! -f "${step_9}" && -f "${step_8}" ]]; then
        echo -e "Started step 9/11: Listing, tallying QNAMEs in"
        echo -e "${tmp_corrected_sorted}."

        list_tally_qnames_gzip_updated "${tmp_corrected_sorted}" && \
        touch "${step_9}"

        echo -e "Completed step 9/11: Listing, tallying QNAMEs in"
        echo -e "${tmp_corrected_sorted}.\n"
    elif [[ -f "${step_9}" && -f "${step_8}" ]]; then
        echo_completion_message 9
    else
        echo_exit_message 9
        exit 1
    fi

    #  10: Create txt.gz outfiles for QNAME >, <, = 2 -------------------------
    if [[ ! -f "${step_10}" && -f "${step_9}" ]]; then
        echo -e "Started step 10/11: Create txt.gz outfiles for QNAME > 2, QNAME < 2, QNAME = 2 from ${tmp_sorted%.bam}.QNAME.txt.gz."

        identify_qnames_updated "${parallelize}" \
        "gt" "${tmp_corrected_sorted%.bam}.QNAME.txt.gz" "gzip" && \
        identify_qnames_updated "${parallelize}" \
        "lt" "${tmp_corrected_sorted%.bam}.QNAME.txt.gz" "gzip" && \
        identify_qnames_updated "${parallelize}" \
        "eq" "${tmp_corrected_sorted%.bam}.QNAME.txt.gz" "gzip" && \
        touch "${step_10}"

        echo -e "Completed step 10/11: Create txt.gz outfiles for QNAME > 2, QNAME < 2, QNAME = 2 from ${tmp_sorted%.bam}.QNAME.txt.gz.\n"
    elif [[ -f "${step_10}" && -f "${step_9}" ]]; then
        echo_completion_message 10
    else
        echo_exit_message 10
        exit 1
    fi
fi


#  11: Remove temporary bams, move "${TMPDIR}" outfiles to "${outpath}" -------
if [[ -f "${step_7}" && -f "${tmp_corrected}" ]]; then
    echo -e "Started step 11/11: Removing temporary bams, moving outfiles from ${TMPDIR} to ${outpath}."

    rm "${tmp}" "${tmp_sorted}" "${tmp_corrected_sorted}" && \
    mv -f "${TMPDIR}/"*{.txt,.txt.gz} "${tmp_corrected}" "${outpath}" && \
    touch "${step_11}"

    echo -e "Completed step 11/11: Removing temporary bams, moving outfiles from ${TMPDIR} to ${outpath}.\n"
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

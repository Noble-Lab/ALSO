#!/bin/bash

#  filter-problematic-qnames.sh
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
    echo "Run pipeline to filter problematic QNAMEs from bam file."
    echo "  - Step 01: Copy bam of interest to \${TMPDIR}."
    echo "  - Step 02: Remove low-quality reads (-f 3 -F 12 -q 30)."
    echo "  - Step 03: Sort filtered bam by coordinate."
    echo "  - Step 04: Index filtered bam file."
    echo "  - Step 05: Identify, write lists of problematic QNAMEs."
    echo "  - Step 06: Sort filtered bam by QNAME."
    echo "  - Step 07: List and tally QNAMEs in the sorted bam file."
    echo "  - Step 08: Create txt.gz file for QNAME > 2."
    echo "  - Step 09: Create txt.gz file for QNAME < 2."
    echo "  - Step 10: Count lines in txt.gz files (optional)."
    echo "  - Step 11: Collect lists in an array."
    echo "  - Step 12: Create master list of unique QNAMEs."
    echo "  - Step 13: Using master list, filter problematic QNAMEs from bam."
    echo "  - Step 14: Run samtools flagstat on bams (optional)."
    echo "  - Step 15: Sort corrected bam."
    echo "  - Step 16: Index corrected bam."
    echo "  - Step 17: Remove unneeded intermediate files (optional)."
    echo "  - Step 18: Move outfiles from \${TMPDIR} to \${outpath}."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - argparser >= 0.7.1"
    echo "  - picard >= 2.27.1"
    echo "  - R >= 4.0.5"
    echo "  - Rsamtools >= 2.6.0"
    echo "  - Rscript >= 4.0.5"
    echo "  - samtools >= 1.13"
    echo "  - scales >= 1.1.1"
    echo "  - tidyverse >= 1.3.1"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h print this help message and exit"
    echo "-u use safe mode: \"TRUE\" or \"FALSE\" (logical; default: FALSE)"
    echo "-l run on GS HPC: \"TRUE\" or \"FALSE\" (logical; default: FALSE)"
    echo "-m initial memory allocation pool for JVM (chr; default \"512m\")"
    echo "-x maximum memory allocation pool for JVM (chr; default \"1g\")"
    echo "-i bam infile, including path (chr)"
    echo "-o path for outfiles (chr); path will be made if it does not exist"
    echo "-f run samtools flagstat on bams: \"TRUE\" or \"FALSE\" (logical)"
    echo "-c count lines: \"TRUE\" or \"FALSE\" (logical)"
    echo "-t tally entries: \"TRUE\" or \"FALSE\" (logical)"
    echo "-n run upto step (int 1-13; default: 13)"
    echo "-r remove intermediate files: \"TRUE\" or \"FALSE\" (logical)"
    echo "-p number of cores for parallelization (int >= 1; default: 1)"
    exit
}


while getopts "h:u:l:m:x:i:o:f:c:t:n:r:p:" opt; do
    case "${opt}" in
        h) print_usage ;;
        u) safe_mode="${OPTARG}" ;;
        l) cluster="${OPTARG}" ;;
        m) memory_min="${OPTARG}" ;;
        x) memory_max="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        f) flagstat="${OPTARG}" ;;
        c) count="${OPTARG}" ;;
        t) tally="${OPTARG}" ;;
        n) run_upto="${OPTARG}" ;;
        r) remove="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${cluster}" ]] && cluster=FALSE
[[ -x "${memory_min}" ]] && memory_min="512m"
[[ -x "${memory_max}" ]] && memory_max="1g"
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${flagstat}" ]] && print_usage
[[ -z "${count}" ]] && print_usage
[[ -z "${tally}" ]] && print_usage
[[ -z "${run_upto}" ]] && run_upto=13
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
check_dependency R
check_dependency samtools

#  Evaluate "${cluster}"
case "$(echo "${cluster}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-c: \"Run on GS HPC\" is TRUE." ;;
    false | f) echo -e "-c: \"Run on GS HPC\" is FALSE." && check_dependency picard ;;
    *) \
        echo -e "Exiting: -c argument must be TRUE or FALSE.\n"
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

#  Evaluate "${flagstat}"
case "$(echo "${flagstat}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        flagstat=1
        echo -e "-f: \"Run samtools flagstat\" is TRUE."
        ;;
    false | f) \
        flagstat=0
        echo -e "-f: \"Run samtools flagstat\" is FALSE."
        ;;
    *) \
        echo -e "Exiting: -f \"flagstat\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

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

#  Evaluate "${run_upto}"
[[ ! "${run_upto}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -n \"run_upto\" argument must be an integer in the range of 1-13.\n"
        exit 1
    }

case "${run_upto}" in
    1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | \
    11 | 12 | 13 | 14 | 15 | 16 | 17 | 18) \
        echo -e "-n: \"Run upto\" is set to step ${run_upto}"
        ;;
    *) \
        echo -e "Exiting: -n \"run_upto\" argument must be an integer in the range of 1-13.\n"
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


#  Assign variables needed for the pipeline -------------------------------
base="$(basename "${infile}")"
base_tmp="${TMPDIR}/${base}"
base_rm="${base%.bam}.rm.bam"
base_rm_bai="${base%.bam}.rm.bam.bai"
base_qname="${base%.bam}.rm.sort-n.bam"
out_rm="${TMPDIR}/${base_rm}"
out_rm_qname="${TMPDIR}/${base_qname}"
out_rm_bai="${TMPDIR}/${base_rm_bai}"
out_ambiguous="${base%.bam}.rm.ambiguous.txt.gz"
out_duplicated="${base%.bam}.rm.duplicated.txt.gz"
out_singleton="${base%.bam}.rm.singleton.txt.gz"
out_trans="${base%.bam}.rm.trans.txt.gz"
out_unmated="${base%.bam}.rm.unmated.txt.gz"
out_qname_full="${out_rm_qname%.bam}.QNAME.txt.gz"
out_qname_gt="${out_rm_qname%.bam}.QNAME.gt.txt.gz"
out_qname_lt="${out_rm_qname%.bam}.QNAME.lt.txt.gz"
out_exclude="${TMPDIR}/${base%.bam}.to-exclude.txt.gz"
out_corrected="${TMPDIR}/${base%.bam}.corrected.bam"
out_flagstat_infile="${TMPDIR}/${base%.bam}.flagstat.txt"
out_flagstat_rm="${TMPDIR}/${base%.bam}.rm.flagstat.txt"
out_flagstat_corrected="${TMPDIR}/${base%.bam}.corrected.flagstat.txt"

#  Assign variables for completion files
step_1="$(echo_completion_file "${outpath}" 1 "${infile}")"
step_2="$(echo_completion_file "${outpath}" 2 "${infile}")"
step_3="$(echo_completion_file "${outpath}" 3 "${infile}")"
step_4="$(echo_completion_file "${outpath}" 4 "${infile}")"
step_5="$(echo_completion_file "${outpath}" 5 "${infile}")"
step_6="$(echo_completion_file "${outpath}" 6 "${infile}")"
step_7="$(echo_completion_file "${outpath}" 7 "${infile}")"
step_8="$(echo_completion_file "${outpath}" 8 "${infile}")"
step_9="$(echo_completion_file "${outpath}" 9 "${infile}")"
step_10="$(echo_completion_file "${outpath}" 10 "${infile}")"
step_11="$(echo_completion_file "${outpath}" 11 "${infile}")"
step_12="$(echo_completion_file "${outpath}" 12 "${infile}")"
step_13="$(echo_completion_file "${outpath}" 13 "${infile}")"
step_14="$(echo_completion_file "${outpath}" 14 "${infile}")"
step_15="$(echo_completion_file "${outpath}" 15 "${infile}")"
step_16="$(echo_completion_file "${outpath}" 16 "${infile}")"
step_17="$(echo_completion_file "${outpath}" 17 "${infile}")"


#  01: Copy files of interest to ${TMPDIR} ------------------------------------
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started step 1/13: Copying ${base} into ${TMPDIR}."
    cp "${infile}" "${TMPDIR}" && \
    touch "${step_1}" && \
    echo -e "Completed step 1/13: Copying ${base} into ${TMPDIR}.\n"
else
    echo_completion_message 1
    :
fi

evaluate_run_upto "${run_upto}" 1


#  02: Remove low-quality reads (-f 3 -F 12 -q 30) ----------------------------
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started step 2/13: Running samtools view (-f 3 -F 12 -q 30) on ${base_tmp}."
    remove_reads_low_quality "${parallelize}" "${base_tmp}" "${out_rm}" && \
    touch "${step_2}" && \
    echo -e "Completed step 2/13: Running samtools view (-f 3 -F 12 -q 30) on ${base_tmp}.\n"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
    :
else
    echo_exit_message 2
    exit 1
fi

evaluate_run_upto "${run_upto}" 2


#  03: Sort bam by coordinate, writing a tmp file that overwrites the infile --
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Started step 3/13: Running samtools sort on ${out_rm}."
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" "${out_rm}" && \
    touch "${step_3}" && \
    echo -e "Completed step 3/13: Running samtools sort on ${out_rm}.\n"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
    :
else
    echo_exit_message 3
    exit 1
fi

evaluate_run_upto "${run_upto}" 3


#  04: Index bam file ---------------------------------------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Started step 4/13: Running samtools index on ${out_rm}."
    index_bam "${parallelize}" "${out_rm}" && \
    touch "${step_4}" && \
    echo -e "Completed step 4/13: Running samtools index on ${out_rm}.\n"
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
    :
else
    echo_exit_message 4
    exit 1
fi

evaluate_run_upto "${run_upto}" 4


#  05: Generate lists of QNAMEs to exclude ------------------------------------
if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    echo -e "Started step 5/13: Generating lists of QNAMEs to exclude from ${out_rm}."
    if [[ -f ./bin/generate-qname-lists.R ]]; then
        script="./bin/generate-qname-lists.R"
    elif [[ -f ./generate-qname-lists.R ]]; then
        script="./generate-qname-lists.R"
    else
        echo -e "Exiting: Could not find \"generate-qname-lists.R\"."
        exit 1
    fi

    Rscript "${script}" \
    --bam "${out_rm}" \
    --bai "${out_rm_bai}" \
    --outdir "${TMPDIR}" \
    --chunk 100000 \
    --mated FALSE \
    --unmated TRUE \
    --ambiguous TRUE \
    --trans TRUE \
    --duplicated FALSE \
    --singleton FALSE \
    --unique TRUE \
    --tally TRUE \
    --remove TRUE
    
    if [[ -f "${TMPDIR}/${out_unmated}" ]]; then
        touch "${step_5}"
    else
        echo -e "There was a problem: ${TMPDIR}/${out_unmated} was not generated."
        echo -e "Check on this."
        echo_exit_message 5
        exit 1
    fi
    echo -e "Completed step 5/13: Generating lists of QNAMEs to exclude from ${out_rm}.\n"
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    echo_completion_message 5
    :
else
    echo_exit_message 5
    exit 1
fi

evaluate_run_upto "${run_upto}" 5


#  06: Sort bam by QNAME ------------------------------------------------------
if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then
    echo -e "Started step 6/18: Sorting ${out_rm} by QNAME."

    sort_bam_qname_samtools "${parallelize}" \
    "${out_rm}" "${out_rm_qname}" && \
    touch "${step_6}"

    echo -e "Completed step 6/18: Sorting ${out_rm} by QNAME.\n"
elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
    echo_completion_message 6
else
    echo_exit_message 6
    exit 1
fi

evaluate_run_upto "${run_upto}" 6


#  07: List and tally QNAMEs in the sorted bam file ---------------------------
if [[ ! -f "${step_7}" && -f "${step_6}" ]]; then
    echo -e "Started step 7/18: Listing, tallying QNAMEs in ${out_rm_qname}."

    list_tally_qnames_gzip_updated "${out_rm_qname}" && \
    touch "${step_7}"

    echo -e "Completed step 7/18: Listing, tallying QNAMEs in ${out_rm_qname}.\n"
elif [[ -f "${step_7}" && -f "${step_6}" ]]; then
    echo_completion_message 7
else
    echo_exit_message 7
    exit 1
fi

evaluate_run_upto "${run_upto}" 7


#  08: Create txt.gz outfiles for QNAME > 2 -----------------------------------
if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
    echo -e "Started step 8/18: Create txt.gz outfiles for QNAME > 2 from ${out_qname_full}."

    identify_qnames_updated "${parallelize}" \
    "gt" "${out_qname_full}" "gzip" && \
    touch "${step_8}"

    echo -e "Completed step 8/18: Create txt.gz outfiles for QNAME > 2 from ${out_qname_full}.\n"
elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
    echo_completion_message 8
else
    echo_exit_message 8
    exit 1
fi

evaluate_run_upto "${run_upto}" 8


#  09: Create txt.gz outfiles for QNAME < 2 -----------------------------------
if [[ ! -f "${step_9}" && -f "${step_8}" ]]; then
    echo -e "Started step 9/18: Create txt.gz outfiles for QNAME < 2 from ${out_qname_full}."

    identify_qnames_updated "${parallelize}" \
    "lt" "${out_qname_full}" "gzip" && \
    touch "${step_9}"

    echo -e "Completed step 9/18: Create txt.gz outfiles for QNAME < 2 from ${out_qname_full}.\n"
elif [[ -f "${step_9}" && -f "${step_8}" ]]; then
    echo_completion_message 9
else
    echo_exit_message 9
    exit 1
fi

evaluate_run_upto "${run_upto}" 9


#  10: Count lines in outfiles with QNAME > 2, QNAME < 2 (optional) -----------
if [[ "${count}" == 1 ]]; then
    if [[ ! -f "${step_10}" && -f "${step_9}" ]]; then
        count_full="$(count_lines_gzip "${out_qname_full}")"
        count_remove_gt="$(count_lines_gzip "${out_qname_gt}")"
        count_remove_lt="$(count_lines_gzip "${out_qname_lt}")"

        echo -e "Started step 10/18: Counting entries in files, writing out lists."

        echo "$(basename "${out_qname_full}"): ${count_full}" && \
        echo "$(basename "${out_qname_gt}"): ${count_remove_gt}" && \
        echo "$(basename "${out_qname_lt}"): ${count_remove_lt}" && \
        echo "" && \
        touch "${step_10}"

        echo -e "Completed step 10/18: Counting entries in files, writing out lists.\n"
    elif [[ -f "${step_10}" && -f "${step_9}" ]]; then
        count_full="$(count_lines_gzip "${out_qname_full}")"
        count_remove_gt="$(count_lines_gzip "${out_qname_gt}")"
        count_remove_lt="$(count_lines_gzip "${out_qname_lt}")"

        echo_completion_message 10

        echo "$(basename "${out_qname_full}"): ${count_full}" && \
        echo "$(basename "${out_qname_gt}"): ${count_remove_gt}" && \
        echo "$(basename "${out_qname_lt}"): ${count_remove_lt}" && \
        echo ""
    else
        echo_exit_message 10
        exit 1
    fi
fi

evaluate_run_upto "${run_upto}" 10


#  11: Collect outfiles in an array -------------------------------------------
find_outfiles() {
    #TODO Description of function
    #
    # :param TMPDIR: path for outfiles (chr)
    # :param out_ambiguous: txt.gz list of ambiguous QNAMEs (chr)
    # :param out_duplicated: txt.gz list of duplicated QNAMEs (chr)
    # :param out_singleton: txt.gz list of singleton QNAMEs (chr)
    # :param out_trans: txt.gz list of trans QNAMEs (chr)
    # :param out_unmated: txt.gz list of unmated QNAMEs (chr)
    # :param out_qname_gt: txt.gz list of QNAMEs with >2 entries (chr)
    # :param out_qname_lt: txt.gz list of QNAMEs with <2 entries (chr)
    find "${TMPDIR}" \
    -maxdepth 1 \
    -type f \
    \( -name "${out_ambiguous}" -o \
    -name "${out_duplicated}" -o \
    -name "${out_singleton}" -o \
    -name "${out_trans}" -o \
    -name "${out_unmated}" -o \
    -name "$(basename "${out_qname_gt}")" -o \
    -name "$(basename "${out_qname_lt}")" \) \
    -print0
}


if [[ ! -f "${step_11}" && -f "${step_9}" ]]; then
    echo -e "Started step 11/18: Collecting outfiles into an array."
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo -e "There was a problem: The number of outfiles is ${#outfiles[@]}. Check on this."
        echo_exit_message 11
        exit 1
    else
        echo -e "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo -e "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo -e ""
        touch "${step_11}" && \
        echo -e "Completed step 11/18: Collecting outfiles into an array.\n"
    fi
elif [[ -f "${step_11}" && -f "${step_9}" ]]; then
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo -e "There was a problem: The number of outfiles is ${#outfiles[@]}. Check on this."
        echo_exit_message 11
        exit 1
    else
        echo -e "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo -e "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo -e ""
    fi

    echo_completion_message 11
    :
else
    echo_exit_message 11
    exit 1
fi

evaluate_run_upto "${run_upto}" 11


#  12: Combine outfiles into one file for excluding problematic QNAME reads ---
if [[ ! -f "${step_12}" && -f "${step_11}" ]]; then
    echo -e "Started step 12/18: Combining outfiles into one file, $(basename "${out_exclude}")."
    combine_gz_qname_lists_return_unique_gzip "${outfiles[@]}" \
    > "${out_exclude}" && \
    touch "${step_12}"

    records=$(count_lines_gzip "${out_exclude}" 2> /dev/null)
    if [[ $(( records )) -eq 0 ]]; then
        echo -e "There was a problem generating $(basename "${out_exclude}"), which"
        echo -e "currently contains $(( records )) lines. Check on this."
        echo_exit_message 12
    else
        echo -e "$(basename "${out_exclude}"): $(( records )) lines"
        echo -e "Completed step 12/18: Combining outfiles into one file, $(basename "${out_exclude}").\n"
    fi
elif [[ -f "${step_12}" && -f "${step_11}" ]]; then
    if [[ -f "${out_exclude}" ]]; then
        records=$(count_lines_gzip "${out_exclude}" 2> /dev/null)
        if [[ $(( records )) -eq 0 ]]; then
            echo -e "There was a problem generating $(basename "${out_exclude}"), which"
            echo -e "currently contains $(( records )) lines. Check on this."
            echo_exit_message 12
        else
            echo -e "$(basename "${out_exclude}"): $(( records )) lines"
            echo_completion_message 12
            :
        fi
    else
        echo "${out_exclude} is missing. Check on this."
        echo_exit_message 12
    fi
else
    echo_exit_message 12
    exit 1
fi

evaluate_run_upto "${run_upto}" 12


#  13: Exclude problematic QNAME reads from bam infile ------------------------
if [[ ! -f "${step_13}" && -f "${step_12}" ]]; then
    echo -e "Started step 13/18: Using ${out_exclude} to exclude reads from ${out_rm}."
    exclude_qname_reads_picard \
    "${out_rm}" \
    "${out_exclude}" \
    "${out_corrected}" \
    "${memory_min}" \
    "${memory_max}" \
    "${cluster}" && \
    touch "${step_13}" && \
    echo -e "Completed step 13/18: Using ${out_exclude} to exclude reads from ${out_rm}.\n"
elif [[ -f "${step_13}" && -f "${step_12}" ]]; then
    echo_completion_message 13
    :
else
    echo_exit_message 13
    exit 1
fi

evaluate_run_upto "${run_upto}" 13


#  14: Run flagstat on bam in and outfiles (optional) -------------------------
if [[ ! -f "${step_14}" && -f "${step_13}" ]]; then
    echo -e "Started step 14/18 (optional): Running samtools flagstat on ${infile}, ${out_rm}, and ${out_corrected}."
    if [[ "${flagstat}" == 1 ]]; then
        run_flagstat "${parallelize}" "${infile}" "${out_flagstat_infile}"
        run_flagstat "${parallelize}" "${out_rm}" "${out_flagstat_rm}"
        run_flagstat "${parallelize}" "${out_corrected}" "${out_flagstat_corrected}"

        echo_flagstat "${out_flagstat_infile}"
        echo_flagstat "${out_flagstat_rm}"
        echo_flagstat "${out_flagstat_corrected}"
        echo -e ""
    fi
    touch "${step_14}" && \
    echo -e "Completed step 14/18 (optional): Running samtools flagstat on ${infile}, ${out_rm}, and ${out_corrected}.\n"
elif [[ -f "${step_14}" && -f "${step_13}" ]]; then
    if [[ "${flagstat}" == 1 ]]; then
        echo_completion_message 14
        echo_flagstat "${out_flagstat_infile}" || echo -e "${out_flagstat_infile} not found."
        echo_flagstat "${out_flagstat_rm}" || echo -e "${out_flagstat_rm} not found."
        echo_flagstat "${out_flagstat_corrected}" || echo -e "${out_flagstat_corrected} not found."
        echo -e ""
    fi
else
    echo_exit_message 14
    exit 1
fi

evaluate_run_upto "${run_upto}" 14


#  15: Sort corrected bam -----------------------------------------------------
if [[ ! -f "${step_15}" && -f "${step_14}" ]]; then
    echo -e "Started step 15/18: Running samtools sort on ${out_corrected}."
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" "${out_corrected}" && \
    touch "${step_15}" && \
    echo -e "Completed step 15/18: Running samtools sort on ${out_corrected}.\n"
elif [[ -f "${step_15}" && -f "${step_14}" ]]; then
    echo_completion_message 15
    :
else
    echo_exit_message 15
    exit 1
fi

evaluate_run_upto "${run_upto}" 15


#  16: Index corrected bam ----------------------------------------------------
if [[ ! -f "${step_16}" && -f "${step_15}" ]]; then
    echo -e "Started step 16/18: Running samtools index on ${out_corrected}."
    index_bam "${parallelize}" "${out_corrected}" && \
    touch "${step_16}" && \
    echo -e "Completed step 16/18: Running samtools index on ${out_corrected}.\n"
elif [[ -f "${step_16}" && -f "${step_15}" ]]; then
    echo_completion_message 16
    :
else
    echo_exit_message 16
    exit 1
fi

evaluate_run_upto "${run_upto}" 16


#  17: Remove unneeded intermediate files (optional) --------------------------
rm -f "${base_tmp}"  # Remove infile copied into "${TMPDIR}"
if [[ ! -f "${step_17}" && -f "${step_16}" ]]; then
    echo -e "Started step 17/18 (optional): Removing unneeded intermediate files."
    if [[ "${remove}" == 1 ]]; then
        rm -f "${outfiles[@]}"
        "${out_rm}" "${out_rm_bai}" "${out_rm_qname}" "${out_qname_full}"
    fi
    touch "${step_17}" && \
    echo -e "Completed step 17/18 (optional): Removing unneeded intermediate files.\n"
elif [[ -f "${step_17}" && -f "${step_16}" ]]; then
    echo_completion_message 17
    :
else
    echo_exit_message 17
    exit 1
fi

evaluate_run_upto "${run_upto}" 17


#  18: Move outfiles from "${TMPDIR}" to "${outpath}" -------------------------
if [[ -f "${step_17}" && -f "${out_corrected}" ]]; then
    echo -e "Started step 18/18: Moving outfiles from ${TMPDIR} to ${outpath}."
    mv -f "${TMPDIR}/"*.{bam,bam.bai,txt,txt.gz} "${outpath}" && \
    touch "${step_13}" && \
    echo -e "Completed step 18/18: Moving outfiles from ${TMPDIR} to ${outpath}.\n"
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

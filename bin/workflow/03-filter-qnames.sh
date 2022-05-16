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
. ./bin/auxiliary/functions-preprocessing.sh ||
    {
        echo "Exiting: Unable to source auxiliary information."
        echo "Are you in the correct working directory," \
        "\"2021_kga0_4dn-mouse-cross\"?"
        exit 1
    }


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Run pipeline to filter problematic QNAMEs from bam file."
    echo "  - Step 1: Remove low quality reads."
    echo "  - Step 2: Sort resulting bam."
    echo "  - Step 3: Index bam."
    echo "  - Step 4: Identify, write lists of problematic QNAMEs."
    echo "  - Step 5: Collect lists in an array."
    echo "  - Step 6: Use array to create a master list of unique QNAMEs."
    echo "  - Step 7: Using master list, filter problematic QNAMEs from bam."
    echo "  - Step 8: Run samtools flagstat on bams (optional)."
    echo "  - Step 9: Sort corrected bam."
    echo "  - Step 10: Index corrected bam."
    echo "  - Step 11: Remove intermediate files from outdirectory (optional)."
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
    echo "-u use safe mode: \"TRUE\" or \"FALSE\" (logical)"
    echo "-c use KA conda environment: \"TRUE\" or \"FALSE\" (logical)"
    echo "-l run on GS HPC: \"TRUE\" or \"FALSE\" (logical)"
    echo "-i bam infile, including path (chr)"
    echo "-o path for outfiles (chr); path will be made if it does not exist"
    echo "-f run samtools flagstat on bams: \"TRUE\" or \"FALSE\" (logical)"
    echo "-r remove intermediate files: \"TRUE\" or \"FALSE\" (logical)"
    echo "-p number of cores for parallelization (int >= 1); default: 1"
    exit
}

while getopts "h:u:c:l:i:o:f:r:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        c) conda="${OPTARG}" ;;
        l) cluster="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        f) flagstat="${OPTARG}" ;;
        r) remove="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${conda}" ]] && print_usage
[[ -z "${cluster}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${flagstat}" ]] && printUsage
[[ -z "${remove}" ]] && printUsage
[[ -z "${parallelize}" ]] && parallelize=1


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-u: Safe mode is on." && set -Eeuxo pipefail ;;
    false | f) echo -e "-u: Safe mode is off." ;;
    *) \
        echo -e "Exiting: -u safe mode argument must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

#  Evaluate "${conda}"
case "$(echo "${conda}" | tr '[:upper:]' '[:lower:]')" in
    false | f) \
        echo -e "-c: Use KA conda environment is FALSE."
        ;;
    true | t) \
        #  Conda anvironment used by KA for writing and testing the pipeline
        echo -e "-c: Use KA conda environment is TRUE."
        conda activate pipeline-test_env  
        ;;
    *) \
        echo -e "Exiting: -c argument must be \"TRUE\" or \"FALSE\".\n"
        return 1
        ;;
esac

#  Check for necessary dependencies; exit if not found
check_dependency R
check_dependency samtools

#  Evaluate "${cluster}"
case "$(echo "${cluster}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-l: Run on GS HPC is TRUE." ;;
    false | f) echo -e "-l: Run on GS HPC is FALSE." && check_dependency picard ;;
    *) \
        echo -e "Exiting: -l argument must be \"TRUE\" or \"FALSE\".\n"
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
        echo "-o: Directory ${outpath} does not exist; making the directory."
        mkdir -p "${outpath}"
    }

#  Evaluate "${flagstat}"
case "$(echo "${flagstat}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        flagstat=1
        echo -e "-f: Run samtools flagstat is TRUE."
        ;;
    false | f) \
        flagstat=0
        echo -e "-f: Run samtools flagstat is FALSE."
        ;;
    *) \
        echo -e "Exiting: -f flagstat argument must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

#  Evaluate "${remove}"
case "$(echo "${remove}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        remove=1
        echo -e "-r: Remove intermediate files is TRUE."
        ;;
    false | f) \
        remove=0
        echo -e "-r: Remove intermediate files is FALSE."
        ;;
    *) \
        echo -e "Exiting: -r remove intermediate files must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

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


#  Assign variables needed for the pipeline -------------------------------
base="$(basename "${infile}")"
base_rm="${base%.bam}.rm.bam"
base_rm_bai="${base%.bam}.rm.bam.bai"
out_rm="${outpath}/${base_rm}"
out_rm_bai="${outpath}/${base_rm_bai}"
out_ambiguous="${base%.bam}.rm.ambiguous.txt.gz"
out_duplicated="${base%.bam}.rm.duplicated.txt.gz"
out_singleton="${base%.bam}.rm.singleton.txt.gz"
out_trans="${base%.bam}.rm.trans.txt.gz"
out_unmated="${base%.bam}.rm.unmated.txt.gz"
out_exclude="${outpath}/${base%.bam}.to-exclude.txt.gz"
out_corrected="${outpath}/${base%.bam}.corrected.bam"
out_flagstat_infile="${outpath}/${base%.bam}.flagstat.txt"
out_flagstat_rm="${outpath}/${base%.bam}.rm.flagstat.txt"
out_flagstat_corrected="${outpath}/${base%.bam}.corrected.flagstat.txt"

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


#  01: Remove low-quality reads (-f 3 -F 12 -q 30) ----------------------------
if [[ ! -f "${step_1}" ]]; then
    remove_reads_low_quality "${parallelize}" "${infile}" "${out_rm}" && \
    touch "${step_1}"
else
    echo_completion_message 1
    :
fi


#  02: Sort bam by coordinate, writing a tmp file that overwrites the infile --
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" "${out_rm}" && \
    touch "${step_2}"
elif [[ -f "${step_2}" && -f "${step_1}" ]]; then
    echo_completion_message 2
    :
else
    echo_exit_message 2
    exit 1
fi


#  03: Index bam file ---------------------------------------------------------
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    index_bam "${parallelize}" "${out_rm}" && \
    touch "${step_3}"
elif [[ -f "${step_3}" && -f "${step_2}" ]]; then
    echo_completion_message 3
    :
else
    echo_exit_message 3
    exit 1
fi


#  04: Generate lists of QNAMEs to exclude ------------------------------------
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    Rscript ./bin/generate-qname-lists.R \
    --bam "${out_rm}" \
    --bai "${out_rm_bai}" \
    --outdir "${outpath}" \
    --chunk 100000 \
    --mated FALSE \
    --unmated TRUE \
    --ambiguous TRUE \
    --trans TRUE \
    --duplicated TRUE \
    --singleton TRUE \
    --unique TRUE \
    --tally FALSE \
    --remove TRUE

    if [[ -f "${outpath}/${out_unmated}" ]]; then
        touch "${step_4}"
    else
        echo -e "There was a problem: ${outpath}/${out_unmated} was not generated."
        echo -e "Check on this."
        echo_exit_message 4
        exit 1
    fi
elif [[ -f "${step_4}" && -f "${step_3}" ]]; then
    echo_completion_message 4
    :
else
    echo_exit_message 4
    exit 1
fi


#  05: Collect outfiles in an array -----------------------------------------
find_outfiles() {
    #TODO Description of function
    #
    # :param outpath: path for outfiles (chr)
    # :param out_ambiguous: txt.gz list of ambiguous QNAMEs (chr) 
    # :param out_duplicated: txt.gz list of duplicated QNAMEs (chr) 
    # :param out_singleton: txt.gz list of singleton QNAMEs (chr) 
    # :param out_trans: txt.gz list of trans QNAMEs (chr) 
    # :param out_unmated: txt.gz list of unmated QNAMEs (chr) 
    find "${outpath}" \
    -maxdepth 1 \
    -type f \
    \( -name "${out_ambiguous}" -o \
    -name "${out_duplicated}" -o \
    -name "${out_singleton}" -o \
    -name "${out_trans}" -o \
    -name "${out_unmated}" \) \
    -print0
}


if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo_exit_message 5
        exit 1
    else
        echo "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo ""
        touch "${step_5}"
    fi
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo_exit_message 5
        exit 1
    else
        echo "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo ""
        touch "${step_5}"
    fi

    echo_completion_message 5
    :
else
    echo_exit_message 5
    exit 1
fi


#  06: Combine outfiles into one file for excluding problematic QNAME reads ---
if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then    
    combine_gz_qname_lists_return_unique_gzip "${outfiles[@]}" \
    > "${out_exclude}" && \
    touch "${step_6}"

    records=$(count_lines_gzip "${out_exclude}" 2> /dev/null)
    if [[ $(( records )) -eq 0 ]]; then
        echo -e "There was a problem generating $(basename "${out_exclude}"), which"
        echo -e "currently contains $(( records )) lines. Check on this."
        echo_exit_message 6
    else
        echo -e "$(basename "${out_exclude}"): $(( records )) lines"
    fi
elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
    if [[ -f "${out_exclude}" ]]; then
        records=$(count_lines_gzip "${out_exclude}" 2> /dev/null)
        if [[ $(( records )) -eq 0 ]]; then
            echo "There was a problem generating $(basename "${out_exclude}"), which"
            echo "currently contains $(( records )) lines. Check on this."
            echo_exit_message 6
        else
            echo_completion_message 6
            :
        fi
    else
        echo "${out_exclude} is missing. Check on this."
        echo_exit_message 6
    fi
else
    echo_exit_message 6
    exit 1
fi


#  07: Exclude problematic QNAME reads from bam infile ------------------------
if [[ ! -f "${step_7}" && -f "${step_6}" ]]; then
    exclude_qname_reads_picard \
    "${out_rm}" \
    "${out_exclude}" \
    "${out_corrected}" \
    "${cluster}" && \
    touch "${step_7}"
elif [[ -f "${step_7}" && -f "${step_6}" ]]; then
    echo_completion_message 7
    :
else
    echo_exit_message 7
    exit 1
fi


#  08: Run flagstat on bam in and outfiles (optional) -------------------------
if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
    if [[ "${flagstat}" == 1 ]]; then
        run_flagstat "${parallelize}" "${infile}" "${out_flagstat_infile}"
        run_flagstat "${parallelize}" "${out_rm}" "${out_flagstat_rm}"
        run_flagstat "${parallelize}" "${out_corrected}" "${out_flagstat_corrected}"

        echo_flagstat "${out_flagstat_infile}"
        echo_flagstat "${out_flagstat_rm}"
        echo_flagstat "${out_flagstat_corrected}"
        echo ""
    fi
    touch "${step_8}"
elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
    if [[ "${flagstat}" == 1 ]]; then
        echo_completion_message 8
        echo_flagstat "${out_flagstat_infile}" || echo "${out_flagstat_infile} not found."
        echo_flagstat "${out_flagstat_rm}" || echo "${out_flagstat_rm} not found."
        echo_flagstat "${out_flagstat_corrected}" || echo "${out_flagstat_corrected} not found."
        echo ""
    fi
else
    echo_exit_message 8
    exit 1
fi


#  09: Sort corrected bam -----------------------------------------------------
if [[ ! -f "${step_9}" && -f "${step_8}" ]]; then
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" "${out_corrected}" && \
    touch "${step_9}"
elif [[ -f "${step_9}" && -f "${step_8}" ]]; then
    echo_completion_message 9
    :
else
    echo_exit_message 9
    exit 1
fi


#  10: Index corrected bam ----------------------------------------------------
if [[ ! -f "${step_10}" && -f "${step_9}" ]]; then
    index_bam "${parallelize}" "${out_corrected}" && \
    touch "${step_10}"
elif [[ -f "${step_10}" && -f "${step_9}" ]]; then
    echo_completion_message 10
    :
else
    echo_exit_message 10
    exit 1
fi


#  11: Remove unneeded intermediate files (optional) ---------------------------
if [[ ! -f "${step_11}" && -f "${step_10}" ]]; then
    if [[ "${remove}" == 1 ]]; then
        rm -f "${outfiles[@]}" "${out_rm}" "${out_rm_bai}"
    fi
    touch "${step_11}"
elif [[ -f "${step_11}" && -f "${step_10}" ]]; then
    echo_completion_message 11
    :
else
    echo_exit_message 11
    exit 1
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo ""

exit 0

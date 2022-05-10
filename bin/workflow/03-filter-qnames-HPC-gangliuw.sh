#!/bin/bash

#  filter-qnames.sh
#  KA


time_start="$(date +%s)"


#  Source functions into environment ------------------------------------------
# shellcheck disable=1091
. ./functions-preprocessing-HPC.sh ||
    {
        echo "Exiting: Unable to source auxiliary information."
        exit 1
    }


#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Run pipeline to filter problematic QNAMEs from bam file."
    echo "  - Step 01: Copy files of interest to \${TMPDIR}."
    echo "  - Step 02: Remove low quality reads."
    echo "  - Step 03: Sort resulting bam."
    echo "  - Step 04: Index bam."
    echo "  - Step 05: Identify, write lists of problematic QNAMEs."
    echo "  - Step 06: Collect lists in an array."
    echo "  - Step 07: Create master list of unique QNAMEs."
    echo "  - Step 08: Using master list, filter problematic QNAMEs from bam."
    echo "  - Step 09: Run samtools flagstat on bams (optional)."
    echo "  - Step 10: Sort corrected bam."
    echo "  - Step 11: Index corrected bam."
    echo "  - Step 12: Remove intermediate files from outdirectory (optional)."
    echo "  - Step 13: Move outfiles from \${TMPDIR} to \${outpath}."
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
        h) print_usage ;;
        u) safe_mode="${OPTARG}" ;;
        c) conda="${OPTARG}" ;;
        l) cluster="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        f) flagstat="${OPTARG}" ;;
        r) remove="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && print_usage
[[ -z "${conda}" ]] && print_usage
[[ -z "${cluster}" ]] && print_usage
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${flagstat}" ]] && print_usage
[[ -z "${remove}" ]] && print_usage
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
        echo -e "-o: Directory ${outpath} does not exist; making the directory."
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
        echo -e "Exiting: -f flagstat argument must be \"TRUE\" or \"FALSE\".\n"
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
base_tmp="${TMPDIR}/${base}"
base_rm="${base%.bam}.rm.bam"
base_rm_bai="${base%.bam}.rm.bam.bai"
out_rm="${TMPDIR}/${base_rm}"
out_rm_bai="${TMPDIR}/${base_rm_bai}"
out_ambiguous="${base%.bam}.rm.ambiguous.txt.gz"
out_duplicated="${base%.bam}.rm.duplicated.txt.gz"
out_singleton="${base%.bam}.rm.singleton.txt.gz"
out_trans="${base%.bam}.rm.trans.txt.gz"
out_unmated="${base%.bam}.rm.unmated.txt.gz"
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


#  05: Generate lists of QNAMEs to exclude ------------------------------------
if [[ ! -f "${step_5}" && -f "${step_4}" ]]; then
    echo -e "Started step 5/13: Generating lists of QNAMEs to exclude from ${out_rm}."
    echo -e "Running generate-qname-lists.R... "
    Rscript ./generate-qname-lists.R \
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
    --remove TRUE && \
    touch "${step_5}" && \
    echo -e "Completed step 5/13: Generating lists of QNAMEs to exclude from ${out_rm}.\n"
elif [[ -f "${step_5}" && -f "${step_4}" ]]; then
    echo_completion_message 5
    :
else
    echo_exit_message 5
    exit 1
fi


#  06: Collect outfiles in an array -----------------------------------------
find_outfiles() {
    #TODO Description of function
    #
    # :param TMPDIR: path for outfiles (chr)
    # :param out_ambiguous: txt.gz list of ambiguous QNAMEs (chr) 
    # :param out_duplicated: txt.gz list of duplicated QNAMEs (chr) 
    # :param out_singleton: txt.gz list of singleton QNAMEs (chr) 
    # :param out_trans: txt.gz list of trans QNAMEs (chr) 
    # :param out_unmated: txt.gz list of unmated QNAMEs (chr) 
    find "${TMPDIR}" \
    -maxdepth 1 \
    -type f \
    \( -name "${out_ambiguous}" -o \
    -name "${out_duplicated}" -o \
    -name "${out_singleton}" -o \
    -name "${out_trans}" -o \
    -name "${out_unmated}" \) \
    -print0
}


if [[ ! -f "${step_6}" && -f "${step_5}" ]]; then
    echo -e "Started step 6/13: Collecting outfiles into an array."
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo -e "There was a problem: The number of outfiles is ${#outfiles[@]}."
        echo -e "Check on this."
        echo_exit_message 6
        exit 1
    else
        echo -e "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo -e "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo -e ""
        touch "${step_6}" && \
        echo -e "Completed step 6/13: Collecting outfiles into an array.\n"
    fi
elif [[ -f "${step_6}" && -f "${step_5}" ]]; then
    unset outfiles
    typeset -a outfiles
    while IFS=" " read -r -d $'\0'; do
        outfiles+=( "${REPLY}" )
    done < <(find_outfiles)

    if [[ ${#outfiles[@]} -eq 0 ]]; then
        echo -e "There was a problem: The number of outfiles is ${#outfiles[@]}."
        echo -e "Check on this."
        echo_exit_message 6
        exit 1
    else
        echo -e "Number of elements in \${outfiles[@]}: ${#outfiles[@]}"
        echo -e "Elements in \${outfiles[@]}:"
        echo_loop "${outfiles[@]}"; echo -e ""
        touch "${step_6}" && \
        echo -e "Completed step 6/13: Collecting outfiles into an array.\n"
    fi

    echo_completion_message 6
    :
else
    echo_exit_message 6
    exit 1
fi


#  07: Combine outfiles into one file for excluding problematic QNAME reads ---
records=$(count_lines_gzip "${out_exclude}" 2> /dev/null)
if [[ ! -f "${step_7}" && -f "${step_6}" ]]; then
    echo -e "Started step 7/13: Combining outfiles into one file, $(basename "${out_exclude}")."
    if [[ $(( records )) -eq 0 ]]; then
        echo -e "There was a problem generating $(basename "${out_exclude}"), which"
        echo -e "currently contains $(( records )) lines. Check on this."
        echo_exit_message 7
    else
        combine_gz_qname_lists_return_unique_gzip "${outfiles[@]}" \
        > "${out_exclude}" && \
        touch "${step_7}" && \
        echo -e "Completed step 7/13: Combining outfiles into one file, $(basename "${out_exclude}").\n"
    fi
elif [[ -f "${step_7}" && -f "${step_6}" ]]; then
    if [[ $(( records )) -eq 0 ]]; then
        echo -e "There was a problem generating $(basename "${out_exclude}"), which"
        echo -e "currently contains $(( records )) lines. Check on this."
        echo_exit_message 7
    else
        echo_completion_message 7
        :
    fi
else
    echo_exit_message 7
    exit 1
fi


#  08: Exclude problematic QNAME reads from bam infile ------------------------
if [[ ! -f "${step_8}" && -f "${step_7}" ]]; then
    echo -e "Started step 8/13: Using ${out_exclude} to exclude reads from ${out_rm}."
    exclude_qname_reads_picard \
    "${out_rm}" \
    "${out_exclude}" \
    "${out_corrected}" \
    "${cluster}" && \
    touch "${step_8}" && \
    echo -e "Completed step 8/13: Using ${out_exclude} to exclude reads from ${out_rm}.\n"
elif [[ -f "${step_8}" && -f "${step_7}" ]]; then
    echo_completion_message 8
    :
else
    echo_exit_message 8
    exit 1
fi


#  09: Run flagstat on bam in and outfiles (optional) -------------------------
if [[ ! -f "${step_9}" && -f "${step_8}" ]]; then
    echo -e "Started step 9/13 (optional): Running samtools flagstat on ${infile}, ${out_rm}, and ${out_corrected}."
    if [[ "${flagstat}" == 1 ]]; then
        run_flagstat "${parallelize}" "${infile}" "${out_flagstat_infile}"
        run_flagstat "${parallelize}" "${out_rm}" "${out_flagstat_rm}"
        run_flagstat "${parallelize}" "${out_corrected}" "${out_flagstat_corrected}"

        echo_flagstat "${out_flagstat_infile}"
        echo_flagstat "${out_flagstat_rm}"
        echo_flagstat "${out_flagstat_corrected}"
        echo -e ""
    fi
    touch "${step_9}" && \
    echo -e "Completed step 9/13 (optional): Running samtools flagstat on ${infile}, ${out_rm}, and ${out_corrected}.\n"
elif [[ -f "${step_9}" && -f "${step_8}" ]]; then
    if [[ "${flagstat}" == 1 ]]; then
        echo_completion_message 9
        echo_flagstat "${out_flagstat_infile}" || echo -e "${out_flagstat_infile} not found."
        echo_flagstat "${out_flagstat_rm}" || echo -e "${out_flagstat_rm} not found."
        echo_flagstat "${out_flagstat_corrected}" || echo -e "${out_flagstat_corrected} not found."
        echo -e ""
    fi
else
    echo_exit_message 9
    exit 1
fi


#  10: Sort corrected bam -----------------------------------------------------
if [[ ! -f "${step_10}" && -f "${step_9}" ]]; then
    echo -e "Started step 10/13: Running samtools sort on ${out_corrected}."
    sort_bam_coordinate_samtools_overwrite_infile \
    "${parallelize}" "${out_corrected}" && \
    touch "${step_10}" && \
    echo -e "Completed step 10/13: Running samtools sort on ${out_corrected}.\n"
elif [[ -f "${step_10}" && -f "${step_9}" ]]; then
    echo_completion_message 10
    :
else
    echo_exit_message 10
    exit 1
fi


#  11: Index corrected bam ----------------------------------------------------
if [[ ! -f "${step_11}" && -f "${step_10}" ]]; then
    echo -e "Started step 11/13: Running samtools index on ${out_corrected}."
    index_bam "${parallelize}" "${out_corrected}" && \
    touch "${step_11}" && \
    echo -e "Completed step 11/13: Running samtools index on ${out_corrected}.\n"
elif [[ -f "${step_11}" && -f "${step_10}" ]]; then
    echo_completion_message 11
    :
else
    echo_exit_message 11
    exit 1
fi


#  12: Remove unneeded intermediate files (optional) ---------------------------
rm -f "${base_tmp}"  # Remove infile copied into "${TMPDIR}"
if [[ ! -f "${step_12}" && -f "${step_11}" ]]; then
    echo -e "Started step 12/13 (optional): Removing unneeded intermediate files."
    if [[ "${remove}" == 1 ]]; then
        rm -f "${outfiles[@]}" "${out_rm}" "${out_rm_bai}"
    fi
    touch "${step_12}" && \
    echo -e "Completed step 12/13 (optional): Removing unneeded intermediate files.\n"
elif [[ -f "${step_12}" && -f "${step_11}" ]]; then
    echo_completion_message 12
    :
else
    echo_exit_message 12
    exit 1
fi


#  13: Move outfiles from "${TMPDIR}" to "${outpath}" -------------------------
if [[ -f "${step_12}" && -f "${out_corrected}" ]]; then
    echo -e "Started step 13/13: Moving outfiles from ${TMPDIR} to ${outpath}."
    mv -f "${TMPDIR}/"*.{bam,bam.bai,txt,txt.gz} "${outpath}" && \
    touch "${step_13}" && \
    echo -e "Completed step 13/13: Moving outfiles from ${TMPDIR} to ${outpath}.\n"
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

#!/bin/bash

#  mark_duplicates.sh
#  KA


#  Source functions into environment ------------------------------------------
check_dependency() {
    # Check if program is available in "${PATH}"; exit if not
    #
    # :param 1: program to be checked (chr)
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            exit 1
        }
}


check_exit() {
    # Check the exit code of a child process
    # 
    # :param 1:
    # :param 2:
    if [ "${1}" == "0" ]; then
        echo "[done] ${2} $(date)"
    else
        err "[error] ${2} returned exit code ${1}"
    fi
}


err() {
    # Print an error meassage, then exit with code 1
    # 
    # :param 1:
    echo "${1} exited unexpectedly"
    exit 1
}


print_usage() {
    # Print the help message and exit
    echo "${help}"
    exit 1
}


#  Handle arguments, assign variables -----------------------------------------
help="""
${0}:
Use picard MarkDuplicates to check duplicates in bam infile.

Dependencies:
  - samtools >= 1.9

Arguments:
  -h print this help message and exit
  -u use safe mode: \"TRUE\" or \"FALSE\" <logical; default: FALSE>
  -i bam infile, including path <chr>
  -o path for outfiles <chr>
"""

while getopts "h:u:i:o:" opt; do
    case "${opt}" in
        h) echo "${help}" ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage


#  Check variable assignments -------------------------------------------------
echo ""
echo "Running ${0}... "

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-u: \"Safe mode\" is TRUE." && set -Eeuxo pipefail ;;
    false | f) \
        echo -e "-u: \"Safe mode\" is FALSE." ;;
    *) \
        echo -e "Exiting: -u \"safe mode\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Check for necessary dependencies; exit if not found
conda activate pipeline-test_env
check_dependency picard

#  Check that "${infile}" exists
[[ -f "${infile}" ]] ||
    {
        echo -e "Exiting: -i ${infile} does not exist.\n"
        exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo -e "-o: Directory ${outpath} does not exist; making the directory.\n"
        mkdir -p "${outpath}"
    }

echo "" 


#  Assign variables needed for the pipeline -----------------------------------
base="$(basename "${infile}")"
outpath_b="${outpath}/mark_bam"
outpath_t="${outpath}/mark_txt"
outfile_m="${base%.bam}.mark"



#  Make subdirectories for mapped and unmapped outfiles -----------------------
mkdir -p "${outpath_b}"
mkdir -p "${outpath_t}"


#  Run picard MarkDuplicates --------------------------------------------------
picard MarkDuplicates \
--INPUT "${infile}" \
--OUTPUT "${outpath_b}/${outfile_m}.bam" \
--METRICS_FILE "${outpath_t}/${outfile_m}.txt"

#!/bin/bash

#  split_bams_mapped-unmapped.sh
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
Split bam infile into two separate bam outfiles: one for mapped reads, the
other for unmapped reads. Names of outfiles will be derived from the infile.

Dependencies:
  - samtools >= 1.9

Arguments:
  -h print this help message and exit
  -u use safe mode: \"TRUE\" or \"FALSE\" <logical; default: FALSE>
  -i bam infile, including path <chr>
  -o path for outfiles <chr>
  -f run samtools flagstat on bams: \"TRUE\" or \"FALSE\" <logical>
  -p number of cores for parallelization <int >= 1; default: 1>
"""

while getopts "h:u:i:o:f:p:" opt; do
    case "${opt}" in
        h) echo "${help}" ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        f) flagstat="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${flagstat}" ]] && print_usage
[[ -z "${parallelize}" ]] && parallelize=1


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
module load samtools/1.14
check_dependency samtools

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
base="$(basename "${infile}")"
outfile_m="${base%.bam}.mapped.bam"
outfile_u="${base%.bam}.unmapped.bam"
outpath_m="${outpath}/mapped"
outpath_u="${outpath}/unmapped"


#  Make subdirectories for mapped and unmapped outfiles -----------------------
mkdir -p "${outpath_m}"
mkdir -p "${outpath_u}"


#  Split bam infile into mapped bam outfile and unmapped bam outfile ----------
echo "[info] Separating bam files..."

samtools view \
-@ "${parallelize}" \
-b -F4 "${infile}" \
-o "${outpath_m}/${outfile_m}"
check_exit $? "samtools"

samtools view \
-@ "${parallelize}" \
-b -f4 "${infile}" \
-o "${outpath_u}/${outfile_u}"
check_exit $? "samtools"


#  Run flagstat on mapped and unmapped outfiles -------------------------------
if [[ $((flagstat)) -eq 1 ]]; then
    #  Make subdirectories
    flag_m="${outpath_m}/flagstat"
    flag_u="${outpath_u}/flagstat"
    mkdir -p "${flag_m}"
    mkdir -p "${flag_u}"

    #  Mapped
    samtools flagstat \
    -@ "${parallelize}" \
    "${outpath_m}/${outfile_m}" \
        > "${flag_m}/${outfile_m%.bam}.flagstat.txt" &
    check_exit $? "samtools"

    #  Unmapped
    samtools flagstat \
    -@ "${parallelize}" \
    "${outpath_u}/${outfile_u}" \
        > "${flag_u}/${outfile_u%.bam}.flagstat.txt" &
    check_exit $? "samtools"

    wait
fi

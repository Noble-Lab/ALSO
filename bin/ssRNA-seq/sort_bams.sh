#!/bin/bash

#  sort_bams.sh
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


sort_bam_coordinate_samtools_auto() {
    # Run samtools sort on a bam infile; bam outfile name is derived from the
    # bam infile name
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    samtools sort -@ "${1}" "${2}" > "${2/.bam/.sort-c.bam}"
}


#  Handle arguments, assign variables -----------------------------------------
help="""
${0}:
Run samtools sort on bam file. Sorted bam outfile will overwrite bam infile.

Arguments:
  -h print this help message and exit
  -u use safe mode: \"TRUE\" or \"FALSE\" <logical; default: FALSE>
  -i infile, including path <chr>
  -p number of cores for parallelization <int > 0>
"""

while getopts "h:u:i:p:" opt; do
    case "${opt}" in
        h) echo "${help}" ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${infile}" ]] && print_usage
[[ -z "${parallelize}" ]] && print_usage


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
check_dependency samtools

#  Check that "${infile}" exists
[[ -f "${infile}" ]] ||
    {
        echo -e "Exiting: -i ${infile} does not exist.\n"
        exit 1
    }

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


#  Run samtools sort ----------------------------------------------------------
sort_bam_coordinate_samtools_auto "${parallelize}" "${infile}" 

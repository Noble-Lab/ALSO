#!/bin/bash

#  combine_STAR_logs.sh
#  KA


#  Source functions into environment ------------------------------------------
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
Concatenate STAR *.Log.final.out files in a user-specified directory; basename
and location of outfile is specified by the user.

Also, do the same for *.Log.out and *.Log.progress.out files in a user-
specified directory; basename and location of outfile is specified by the user.

Arguments:
  -h print this help message and exit
  -u use safe mode: \"TRUE\" or \"FALSE\" <logical; default: FALSE>
  -i path for infiles <chr>
  -o basename for txt outfiles, including path <chr>
"""

while getopts "h:u:i:o:" opt; do
    case "${opt}" in
        h) echo "${help}" ;;
        u) safe_mode="${OPTARG}" ;;
        i) inpath="${OPTARG}" ;;
        o) outbase="${OPTARG}" ;;
        *) print_usage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${inpath}" ]] && print_usage
[[ -z "${outbase}" ]] && print_usage


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

#  Check that "${infile}" exists
[[ -d "${inpath}" ]] ||
    {
        echo -e "Exiting: -i ${inpath} does not exist.\n"
        exit 1
    }

#  Make "${outbase}" if it doesn't exist
[[ -d "${outbase}" ]] ||
    {
        echo -e "-o: Directory ${outbase} does not exist; making the directory.\n"
        mkdir -p "${outbase}"
    }

echo "" 


#  Assign variables needed for the pipeline -----------------------------------
outpath="$(dirname "${outbase}")"
outbase="$(basename "${outbase}")"


#  Make subdirectories for mapped and unmapped outfiles -----------------------
if [[ ! -d "${outpath}" ]]; then
    mkdir -p "${outpath}"
fi


#  DO the concatenations ------------------------------------------------------
# stackoverflow.com/questions/5917413/concatenate-multiple-files-but-include-filename-as-section-headers
tail -n +1 -- "${inpath}/"*.Log.final.out > "${outpath}/${outbase}.Log.final.out"
tail -n +1 -- "${inpath}/"*.Log.progress.out > "${outpath}/${outbase}.Log.progress.out"
tail -n +1 -- "${inpath}/"*.Log.out > "${outpath}/${outbase}.Log.out"

#  bash combine_STAR_logs.sh -i ./alignments -o ./alignments/Berletch_PRJNA256188/Berletch_PRJNA256188

#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=7:59:59
#$ -l mfree=8G
#$ -l gpgpu=FALSE
#$ -pe serial 1
#$ -q noble-short.q
#$ -cwd

#  KA

#  Description:
#+ See printUsage()
#+
#+ Usage:
#+ See printUsage()
#+
#+ Notes:
#+

printUsage() {
    echo ""
    echo "${0}"
    echo "Run bowtie2-build."
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: TRUE or FALSE (logical)>"
    echo "-e <activate conda environment: TRUE or FALSE (logical)>"
    echo "-p <number of jobs to run in parallel (int)>"
    echo "-l <list from make_list_${0} (file)>"
    exit
}

while getopts "h:u:e:p:l:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        e) environment="${OPTARG}" ;;
        p) parallel="${OPTARG}" ;;
        l) list="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${environment}" ]] && printUsage
[[ -z "${parallel}" ]] && printUsage
[[ -z "${list}" ]] && printUsage

# shellcheck disable=SC1091
. "./bin/auxiliary/auxiliary.sh" ||
    {
        echo "Exiting: Auxiliary information not found."
        echo "Are you in the correct working directory," \
        "\"2021_kga0_4dn-mouse-cross\"?"
        exit 1
    }

checkSafeMode

activateEnvironmentBowtie2

checkDependencyParallel

checkDependencyBowtie2

reportExperimentStart

#  Log experiment
parallel --header : --colsep " " -k -j 1 echo \
"bowtie2-build {path_withdraw} {path_deposit}/{prefix}" \
:::: "${list}"
echo ""

#  Run experiment
parallel --header : --colsep " " -k -j "${parallel}" \
"bowtie2-build {path_withdraw} {path_deposit}/{prefix}" \
:::: "${list}"

reportExperimentEnd
exit 0

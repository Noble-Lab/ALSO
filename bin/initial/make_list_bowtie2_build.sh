#!/bin/bash

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
    echo "Make an infile/outfile list to run bowtie2-build."
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: TRUE or FALSE (logical)>"
    echo "-w <path to directory containing *.fa infile (chr)>"
    echo "-d <path to directory to deposit outfile (chr)>"
    echo "-p <prefix for bowtie2 indices (chr)>"
    echo "-l <path and filename.extension for infile/outfile list (chr)>"
    exit
}

while getopts "h:u:w:d:p:l:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        w) path_withdraw="${OPTARG}" ;;
        d) path_deposit="${OPTARG}" ;;
        p) prefix="${OPTARG}" ;;
        l) list="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${path_withdraw}" ]] && printUsage
[[ -z "${path_deposit}" ]] && printUsage
[[ -z "${prefix}" ]] && printUsage
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

checkDependencyParallel

reportExperimentStart

where="log/kga0"
when="$(date '+%Y-%m%d')"
what="$(basename "${0}" .sh)"
where_when_what="${where}/${when}_${what}"
checkDirectoryPresenceWhereWhenWhat

checkPresenceList

#  Header
echo "path_withdraw path_deposit prefix" >> "${where_when_what}/${list}"

#  Body
parallel --header : --colsep " " -k -j 1 echo \
"{path_withdraw}/{fasta} \
{path_deposit} \
{prefix}" \
::: path_withdraw "${path_withdraw}" \
::: fasta "$(findFa)" \
::: path_deposit "${path_deposit}" \
::: prefix "${prefix}" \
>> "${where_when_what}/${list}"

reportExperimentEnd
exit 0

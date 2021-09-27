#!/bin/bash

#  KA

#  Description:
#+ Run bowtie2-build.
#+
#+ Usage:
#+
#+ Notes:
#+


# set -Eeuxo pipefail
# set -euo pipefail

#  Set working directory and source functions --------------------------------- 
#  Check for proper "$(pwd)"
if [[ "$(basename "$(pwd)")" != "2021_kga0_4dn-mouse-cross" ]]; then
    echo "Exiting: You must run this script from" \
    "\"2021_kga0_4dn-mouse-cross\", the project's top directory."
    exit 1
fi

#  Source auxiliary information
# shellcheck disable=SC1091
. "./bin/auxiliary/auxiliary.sh" ||
    {
        echo "Exiting: Auxiliary information not found."
        echo "Are you in the correct working directory," \
        "\"2021_kga0_4dn-mouse-cross\"?"
        exit 1
    }


#  Define scripts in use ------------------------------------------------------
script_driver="$(basename "${0}" .sh)"
script_main="${script_driver#*_}.sh"
script_make_list="make_list_${script_main}"


#  Set up and enable logging for the driver -----------------------------------
where="log/kga0"
when="$(date '+%Y-%m%d')"
where_when_what="${where}/${when}_${script_driver}"
export where_when_what
checkDirectoryPresenceWhereWhenWhat

stderr="${script_driver}.stderr.txt"
stdout="${script_driver}.stdout.txt"
export stderr
export stdout

checkPresenceErr
checkPresenceOut

printSaveErr
printSaveOut


#  Set up and enable logging for the grid jobs --------------------------------
job_err_out="${where}/${when}_${script_main%.sh}"
[[ -d "${job_err_out}" ]] || mkdir -p "${job_err_out}"


#  Configure variables, arrays for lists, etc. --------------------------------
list_prefix="list.${script_main%.sh}"
list_suffix="${when}.txt"

path_base="/net/noble/vol1/home/kga0/genomes"
typeset -a array_withdraw=(
    "${path_base}/ENA.GCA_001624185.129S1-SvImJ/All-Seq-FASTA"
    "${path_base}/ENA.GCA_001624445.CAST-EiJ/All-Seq-FASTA"
    "${path_base}/ENA.GCA_001624865.SPRET-EiJ/All-Seq-FASTA"
    "${path_base}/ENA.GCA_001632555.C57BL-6NJ/All-Seq-FASTA"
)
# echo "array_withdraw is ${array_withdraw[*]}"

typeset -a array_withdraw_basename
for i in "${array_withdraw[@]}"; do
    array_withdraw_basename+=( "$(basename "${i/\/All-Seq-FASTA/}")" )
done
# echo "array_withdraw_basename is ${array_withdraw_basename[*]}"

typeset -a array_deposit=(
    "${path_base}/ENA.GCA_001624185.129S1-SvImJ/bowtie2"
    "${path_base}/ENA.GCA_001624445.CAST-EiJ/bowtie2"
    "${path_base}/ENA.GCA_001624865.SPRET-EiJ/bowtie2"
    "${path_base}/ENA.GCA_001632555.C57BL-6NJ/bowtie2"
)
# echo "array_deposit is ${array_deposit[*]}"

typeset -a array_prefix=(
    "129S1-SvImJ"
    "CAST-EiJ"
    "SPRET-EiJ"
    "C57BL-6NJ"
)
# echo "array_prefix is ${array_prefix[*]}"

#  Define variables for grid jobs ---------------------------------------------
h_rt="7:59:59"
mfree="8G"
queue="noble-short.q"
pe_serial="1"
gnu_parallel_jobs="${pe_serial}"
disk_free="5G"

# safe_mode="TRUE"
safe_mode="FALSE"
conda_environment="TRUE"

#  Maximum concurrent jobs on grid
grid_jobs_max="4"


#  Create "master" lists ------------------------------------------------------

#  Create deposit path if it does not exist
parallel --header : --colsep " " -k -j 1 \
"[[ -d {deposit} ]] || mkdir -p {deposit}" \
::: deposit "${array_deposit[@]}"

#  Create multi-line lists
parallel --header : --colsep " " -k -j 1 \
"bash bin/{script_make_list} \
-u {safe_mode} \
-w {withdraw} \
-d {deposit} \
-p {prefix} \
-l {list_prefix}.{withdraw_basename}.{list_suffix}" \
::: script_make_list "${script_make_list}" \
::: safe_mode "${safe_mode}" \
::: withdraw "${array_withdraw[@]}" \
:::+ deposit "${array_deposit[@]}" \
:::+ prefix "${array_prefix[@]}" \
:::+ withdraw_basename "${array_withdraw_basename[@]}" \
::: list_prefix "${list_prefix}" \
::: list_suffix "${list_suffix}"


#  Job submission, etc. -------------------------------------------------------
findMultiLineLists | while read -r job; do
    grid_jobs_tally="$(checkJobs | grep -c "${script_main}")"

    while [[ ${grid_jobs_tally} -ge "${grid_jobs_max}" ]]; do
        sleep 5
        printf "."
        grid_jobs_tally="$(checkJobs | grep -c "${script_main}")"
    done

    echo ""
    echo "# -------------------------------------------------------"
    printf "Job:\n%s\n\n" "${job}"
    printf "Job contents:\n%s\n\n" "$(cat "${job}")"

    qsub \
    -S /bin/bash \
    -l h_rt="${h_rt}" \
    -l mfree="${mfree}" \
    -l gpgpu=FALSE \
    -l disk_free="${disk_free}" \
    -pe serial "${pe_serial}" \
    -q "${queue}" \
    -cwd \
    -o "${job_err_out}" \
    -e "${job_err_out}" \
    "bin/${script_main}" \
    -u "${safe_mode}" \
    -e "${conda_environment}" \
    -p "${gnu_parallel_jobs}" \
    -l "${job}"

    printf "Job submission time: %s\n\n" "$(date)"
    sleep 1
done

#  For debugging $(find)
# script_main="bowtie2_build.sh"
# where="log/kga0"
# when="$(date '+%Y-%m%d')"
# list_prefix="list.${script_main%.sh}"
# list_suffix="${when}.txt"

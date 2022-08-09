#!/bin/bash

#  driver_bulk-RNA-seq.sh
#  KA


#  Activate appropriate environment -------------------------------------------
# shellcheck disable=1091
. "${HOME}/.bashrc"

# shellcheck disable=1091
. activate pipeline-test_env


#  Set up appropriate variables for running ALSO jobs -------------------------
# memory="18G"
memory="67G"
email="kga0@uw.edu"
date="2022-0731"
run_up_to=4
heap_min=512
# heap_max=16384
heap_max=65536

dir_base="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross"
dir_results="${dir_base}/results"

inpath="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data"


#  Set array elements for samples, including paths ----------------------------
unset samples
typeset -a samples
while IFS=" " read -r -d $'\0'; do
    samples+=( "${REPLY}" )
done < <(find "${inpath}" -type f -name "*.primary.rm.bam" -print0 | LC_ALL=C sort -z)
# for i in "${samples[@]}"; do echo "${i}"; done  # Check

for i in "${!samples[@]}"; do
    samples[i]="$( echo "${samples[i]%.Aligned.sortedByCoord.out.primary.rm.bam}" | awk -F '.' '{ print $1"."$2 }' )"
done  # Check

unset clean
typeset -a clean
while IFS= read -r -d $'\0'; do
    clean+=( "${REPLY}" )
done < <(printf "%s\0" "${samples[@]}" | sort -uz)
# for i in "${clean[@]}"; do echo "${i}"; done  # Check


#  Submit ALSO jobs -----------------------------------------------------------
# shellcheck disable=2066
# for i in "${clean[@]}"; do
# for i in "${clean[4]}"; do
for i in "${clean[@]:5:6}"; do
    sample="${i}"
    abbrev=$(basename "${sample}" ".Aligned.sortedByCoord.out.primary.rm.bam")
    job_name="ALSO_${abbrev}"
    dir_experiment="${dir_results}/2022-0731_${job_name}"

    # [[ ! -d "${dir_experiment}" ]] && mkdir "${dir_experiment}"
    #
    # chmod 751 submit.sh

    #  Echo test for submitting the job
    echo -e "\
    qsub\n\
    -l mfree=\"${memory}\"\n\
    -m bea\n\
    -M \"${email}\"\n\
    -N \"${job_name}\"\n\
    \"${dir_base}/submit_bulk-RNA-seq.sh\"\n\
    \"${dir_experiment}\"\n\
    \"${sample}\"\n\
    ${run_up_to}\n\
    ${heap_min}\n\
    ${heap_max}\n\
    > \"${date}-ALSO_${abbrev}.submit\"
    " > "${date}-ALSO_${abbrev}.submit"
    echo "" >> "${date}-ALSO_${abbrev}.submit"

    #  Echo test for submit script
    # bash "${dir_base}/submit_bulk-RNA-seq.sh" \
    # "${dir_experiment}" \
    # "${sample}" \
    # "${run_up_to}" \
    # "${heap_min}" \
    # "${heap_max}"

    #  Submitting the job
    qsub \
    -l mfree="${memory}" \
    -m bea \
    -M "${email}" \
    -N "${job_name}" \
    "${dir_base}/submit_bulk-RNA-seq.sh" \
    "${dir_experiment}" \
    "${sample}" \
    "${run_up_to}" \
    "${heap_min}" \
    "${heap_max}" \
    >> "${date}-ALSO_${abbrev}.submit"
    echo "Submitted."
    # :param 1: name of experiment directory, including path <chr>
    # :param 2: bam infile, including path <chr>
    # :param 3: step in pipeline to run up to <int 1-4>
    # :param 4: JVM min heap value; do not include '-Xms' <int > 0>
    # :param 5: JVM min heap value; do not include '-Xmx' <int > 0 && int ≥ param 4>
done

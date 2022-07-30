#!/bin/bash

#  submit_exclude_non-primary-alignments.sh
#  KA


#  Arguments
script="exclude_non-primary-alignments.sh"
project="${1:-"sage"}"
h_rt="${2:-"7:59:59"}"
mfree="${3:-"2G"}"
pe_serial="${4:-"2"}"
queue="${5:-"sage-short.q"}"
inpath="${6:-"./alignments_None"}"
outpath="${7:-"./alignments_None/alignments_primary"}"


#  Infiles
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "${inpath}" -type f -name "*.bam" -print0)


#  Load HPC module, etc.
module load samtools/1.14


#  Submit jobs
for i in "${infiles[@]}"; do
    echo "Submitting job for ${i}"

    qsub \
    -S "/bin/bash" \
    -R "y" \
    -P "${project}" \
    -l h_rt="${h_rt}" \
    -l mfree="${mfree}" \
    -l gpgpu=FALSE \
    -pe serial "${pe_serial}" \
    -q "${queue}" \
    -cwd \
    "${script}" \
    -u FALSE \
    -i "${i}" \
    -o "${outpath}" \
    -f TRUE \
    -p "${pe_serial}"

    sleep .5
    echo ""
done

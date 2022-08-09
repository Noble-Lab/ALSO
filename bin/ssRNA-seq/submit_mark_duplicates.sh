#!/bin/bash

#  submit_mark_duplicates.sh
#  KA


#  Source functions into environment ------------------------------------------
#  Arguments
script="mark_duplicates.sh"
project="${1:-"sage"}"
h_rt="${2:-"7:59:59"}"
mfree="${3:-"8G"}"
queue="${4:-"sage-short.q"}"
inpath="${5:-"./splits/mapped"}"
outpath="${6:-"./splits/mapped"}"


#  Infiles
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "${inpath}" -type f -name "*.bam" -print0)
for i in "${infiles[@]}"; do echo "${i}"; done

#  Load HPC module, etc.
conda activate pipeline-test_env


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
    -q "${queue}" \
    -cwd \
    "${script}" \
    -u FALSE \
    -i "${i}" \
    -o "${outpath}"

    sleep .5
    echo ""
done

# bash mark_duplicates.sh \
# -i "./splits/mapped/brain_rep_2.SRR1525407.SPRET-EiJ.Aligned.out.mapped.bam" \
# -o "./splits/mapped"

# mark_duplicates.sh:
# Use picard MarkDuplicates to check duplicates in bam infile.
#
# Dependencies:
#   - samtools >= 1.9
#
# Arguments:
#   -h print this help message and exit
#   -u use safe mode: "TRUE" or "FALSE" <logical; default: FALSE>
#   -i bam infile, including path <chr>
#   -o path for outfiles <chr>

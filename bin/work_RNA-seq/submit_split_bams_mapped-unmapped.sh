#!/bin/bash

#  submit_split_bams_mapped-unmapped.sh
#  KA


#  Arguments
script="split_bams_mapped-unmapped.sh"
project="${1:-"sage"}"
h_rt="${2:-"7:59:59"}"
mfree="${3:-"2G"}"
pe_serial="${4:-"2"}"
queue="${5:-"sage-short.q"}"
inpath="${6:-"./alignments"}"
outpath="${7:-"./splits"}"


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
    # prefix="$(basename "${script}" ".sh")"
    # suffix="$(basename "${i}" ".Aligned.out.bam")"
    # job_err_out="${prefix}.${suffix}"
    # -o "${job_err_out}" \
    # -e "${job_err_out}" \

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


# split_bams_mapped-unmapped.sh:
# Split bam infile into two separate bam outfiles: one for mapped reads, the
# other for unmapped reads. Names of outfiles will be derived from the infile.
#
# Dependencies:
#   - samtools >= 1.9
#
# Arguments:
#   -h print this help message and exit
#   -u use safe mode: "TRUE" or "FALSE" <logical; default: FALSE>
#   -i bam infile, including path <chr>
#   -o path for outfiles <chr>
#   -f run samtools flagstat on bams: "TRUE" or "FALSE" <logical>
#   -p number of cores for parallelization <int >= 1; default: 1>

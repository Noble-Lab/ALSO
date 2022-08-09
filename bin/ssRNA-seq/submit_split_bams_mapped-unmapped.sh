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
inpath="${6:-"./Berletch_Fang/alignments_primary"}"
outpath="${7:-"./Berletch_Fang/alignments_primary"}"


#  Infiles
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "${inpath}" -type f -name "*.bam" -print0)


#  Load HPC module, etc.
module load samtools/1.14


#  Echo test
# for i in "${infiles[@]}"; do
#     echo "Submitting job for ${i}"
#
#     echo -e "qsub\n \
#     -S \"/bin/bash\"\n \
#     -R \"y\"\n \
#     -P \"${project}\"\n \
#     -l h_rt=\"${h_rt}\"\n \
#     -l mfree=\"${mfree}\"\n \
#     -l gpgpu=FALSE\n \
#     -pe serial \"${pe_serial}\"\n \
#     -q \"${queue}\"\n \
#     -cwd\n \
#     \"${script}\"\n \
#     -u FALSE\n \
#     -i \"${i}\"\n \
#     -o \"${outpath}\"\n \
#     -f TRUE\n \
#     -p \"${pe_serial}\"\n
#     "
#
#     echo ""
# done

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

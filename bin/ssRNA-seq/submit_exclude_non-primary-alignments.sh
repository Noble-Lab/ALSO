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
# inpath="${6:-"./alignments_None"}"
inpath="${6:-"./Berletch_Fang"}"
outpath="${7:-"./Berletch_Fang/alignments_primary"}"


#  Infiles
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "${inpath}" -type f -name "*.bam" -print0 | sort -z)
for i in "${infiles[@]}"; do echo "${i}"; done

#  Load HPC module, etc.
module load samtools/1.14


#  Echo test
# for i in "${infiles[@]}"; do
#     echo "Submitting job for ${i}"

#     echo -e "\tqsub\n \
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

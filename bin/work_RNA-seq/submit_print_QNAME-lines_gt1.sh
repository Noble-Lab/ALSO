#!/bin/bash

#  submit_print_QNAME-lines_gt1.sh
#  KA


#  Parse positional arguments and set up variables ----------------------------
script="print_QNAME-lines_gt1.sh"
project="${1:-"sage"}"
h_rt="${2:-"7:59:59"}"
mfree="${3:-"1G"}"
pe_serial="${4:-"1"}"
queue="${5:-"sage-short.q"}"
inpath="${6:-"./alignments_None/alignments_primary"}"
threshold="${7:-1}"

unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "${inpath}" -maxdepth 1 -type f -name "*.QNAME.txt.gz" -print0)
for i in "${infiles[@]}"; do echo "${i}"; done


#  Load module and submit jobs ------------------------------------------------
# shellcheck disable=1091
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble gzip/1.12


for i in "${infiles[@]}"; do
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
    "${threshold}" \
    "${i}"
done

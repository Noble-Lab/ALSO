#!/bin/bash

#  run_fastqc.sh
#  KA



#  Infiles
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "." -type f -name "*.fastq.gz" -print0)

#  Set up directory for fastqc results
mkdir -p fastqc

#  Run fastqc
for i in "${infiles[@]}"; do
    echo "${i}"
    fastqc "${i}" -o "./fastqc"
done

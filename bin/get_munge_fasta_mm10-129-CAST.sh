#!/bin/bash

#  get_munge_fasta_mm10-129-CAST.sh
#  KA

directory="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/genomes"

cd "${directory}" ||
    {
        echo "Exiting: Directory not found."
        return 1 2> /dev/null
        exit 1
    }

# sGrab() { scp -p nexus:"$@" .; }
sGrab /net/noble/vol1/home/kga0/genomes/Ensembl.129S1-SvImJ/129S1-SvImJ.2bit
sGrab /net/noble/vol1/home/kga0/genomes/Ensembl.CAST-EiJ/CAST-EiJ.2bit
sGrab /net/noble/vol1/home/kga0/genomes/UCSC.mm10/mm10.fa.gz

fasta_mm10_full="mm10.full.fa"
fasta_mm10="mm10.fa"
fasta_129="129S1-SvImJ.fa"
fasta_CAST="CAST-EiJ.fa"
gzip -cd "mm10.fa.gz" > "${fasta_mm10_full}"

twoBitInfo "${fasta_129/.fa/.2bit}" stdout
twoBitInfo "${fasta_CAST/.fa/.2bit}" stdout

twoBitToFa "${fasta_129/.fa/.2bit}" "${fasta_129}"
twoBitToFa "${fasta_CAST/.fa/.2bit}" "${fasta_CAST}"

unset chr
typeset -a chr=($(seq 1 19) "X")
for i in "${chr[@]}"; do
    samtools faidx "${fasta_mm10_full}" "chr${i}" >> "${fasta_mm10}"
done

samtools faidx "${fasta_mm10}"
samtools faidx "${fasta_129}"
samtools faidx "${fasta_CAST}"

picard CreateSequenceDictionary R="${fasta_mm10}" O="${fasta_mm10/.fa/.dict}"
picard CreateSequenceDictionary R="${fasta_129}" O="${fasta_129/.fa/.dict}"
picard CreateSequenceDictionary R="${fasta_CAST}" O="${fasta_CAST/.fa/.dict}"

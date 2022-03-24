#!/bin/bash

#  test_deduplication-issues_samtools.sh
#  KA

project="${1:-"test_deduplication-issues_samtools"}"
path_base="${2:-"/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"}"

cd "${path_base}" ||
    {
        echo -e "Exiting: Directory ${path_base} not found."
    }

[[ -d "${project}" ]] || mkdir -p "${project}"

cd "${project}" ||
    {
        echo -e "Exiting: Directory ${project} not found."
    }

bam="${3:-"test"}.full.bam"
bam_init="F121-6-CASTx129.undifferentiated.dedup.bam"
bam_source="${4:-"/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2021-1026_Bonora-et-al_129S1-SvImJ_test/get_unique_fragments"}"
[[ -f "${bam}" ]] || scp -p nexus:"${bam_source}/${bam_init}" .
[[ -f "${bam}" ]] || mv "${bam_init}" "${bam}"

parallel="${5:-"7"}"
# samtools sort -@ "${parallel}" "${bam}" -O "${bam/bam/sort.bam}"
# rm -- "${bam/bam/sort.bam}"

#  Sort by coordinate
#  htslib.org/doc/samtools-sort.html
# samtools sort -@ "${parallel}" "${bam}" -o "${bam/bam/sort.bam}"
# samtools view "${bam/bam/sort.bam}" | head -50

bam_queryname="${bam/bam/sort-queryname.bam}"
samtools sort -n -@ "${parallel}" "${bam}" -o "${bam_queryname}"
samtools view "${bam_queryname}" | head -50

#  "Fix" mates
#  htslib.org/doc/samtools-fixmate.html
bam_fixmate="${bam_queryname/bam/fixmate.bam}"
samtools fixmate -c -m -@ "${parallel}" "${bam_queryname}" "${bam_fixmate}"
samtools view "${bam_fixmate}" | head -50

#  samtools markdup needs position order
bam_coordinate="${bam_fixmate/bam/sort-coordinate.bam}"
samtools sort -@ "${parallel}" "${bam_fixmate}" -o "${bam_coordinate}"

distance="100"
mode="t"
bam_markdup="${bam_coordinate/bam/markdup}"
# $ samtools markdup -@ "${parallel}" -s -f "${bam_markdup}.txt" -d "${distance}" -c -t --mode "${mode}" "${bam_coordinate}" "${bam_markdup}.bam"
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGCGCTTATCGTGGCGGAGGTCGTACTGAC:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGGCTTGAAGAGGTAATGATCGTATAGCCT:0 for optical duplicate marking.
# [markdup] warning: cannot decipher read name ATTACTCGGCTTGAAGAGGTAATGATCGTATAGCCT:0 for optical duplicate marking.
# [markdup] warning: 10 decipher read name warnings.  New warnings will not be reported.

samtools markdup -@ "${parallel}" -s -f "${bam_markdup}.txt" -c -t --mode "${mode}" "${bam_coordinate}" "${bam_markdup}.bam"

bam_rmdup="${bam_coordinate/bam/rmdup}"
samtools markdup -@ "${parallel}" -r -s -f "${bam_rmdup}.txt" -c -t --mode "${mode}" "${bam_coordinate}" "${bam_rmdup}.bam"
# samtools view "${bam_rmdup}.bam" | head -50

samtools index -@ "${parallel}" -b "${bam_rmdup}.bam"

samtools view -h "${bam_rmdup}.bam" chrX -b > "${bam_rmdup}.chrX.bam"
samtools index -@ "${parallel}" -b "${bam_rmdup}.chrX.bam"

# samtools view -h "${bam_rmdup}.bam" chr19 -b > "${bam_rmdup}.chr19.bam"
# samtools index -@ "${parallel}" -b "${bam_rmdup}.chr19.bam"

# sGrab "${bam_source}/F121-6-CASTx129.undifferentiated.duplicate_report.txt"
# sGrab "${bam_source}/F121-6-CASTx129.undifferentiated.bed.gz"
# sGrab "${bam_source}/F121-6-CASTx129.undifferentiated.fragments.txt.gz"
# sGrab "${bam_source}/F121-6-CASTx129.undifferentiated.insert_sizes.txt"

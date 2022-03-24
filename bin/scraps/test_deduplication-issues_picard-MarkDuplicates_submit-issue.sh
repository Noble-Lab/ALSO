#!/bin/bash

#  test_deduplication-issues_picard-MarkDuplicatesWithMateCigar.sh
#  KA

project="${1:-"test_deduplication-issues_picard-MarkDuplicatesWithMateCigar"}"
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


#  Part 1
# bam_markdup="${bam/bam/MarkDuplicatesWithMateCigar}"
# picard MarkDuplicatesWithMateCigar --INPUT "${bam}" --OUTPUT "${bam_markdup}.bam" --METRICS_FILE "${bam_markdup}.txt"

bam_rmdup="${bam/bam/MarkDuplicatesWithMateCigar-remove}"
picard MarkDuplicatesWithMateCigar --INPUT "${bam}" --OUTPUT "${bam_rmdup}.bam" --METRICS_FILE "${bam_rmdup}.txt" --REMOVE_DUPLICATES true

# picard BuildBamIndex --INPUT "${bam_rmdup}.bam"

# samtools view -h "${bam_rmdup}.bam" chrX -b > "${bam_rmdup}.chrX.bam"
# picard BuildBamIndex --INPUT "${bam_rmdup}.chrX.bam"


#  Part 2
#NOTE Doesn't work b/c, for MarkDuplicatesWithMateCigar, 
# bam_queryname="${bam/bam/sort-queryname.bam}"
# picard SortSam --INPUT "${bam}" --OUTPUT "${bam_queryname}" --SORT_ORDER queryname
#
# bam_rmdup="${bam_queryname/bam/MarkDuplicatesWithMateCigar-remove}"
# picard MarkDuplicatesWithMateCigar --INPUT "${bam_queryname}" --OUTPUT "${bam_rmdup}.bam" --METRICS_FILE "${bam_rmdup}.txt" --REMOVE_DUPLICATES true
#
# bam_coordinate="${bam_rmdup}.sort-queryname.bam"
# picard SortSam --INPUT "${bam_rmdup}.bam" --OUTPUT "${bam_coordinate}" --SORT_ORDER coordinate
#
# picard BuildBamIndex --INPUT "${bam_coordinate}"
#
# samtools view -h "${bam_coordinate}" chrX -b > "${bam_coordinate/bam/chrX.bam}"
# picard BuildBamIndex --INPUT "${bam_coordinate/bam/chrX.bam}"


#  Processing for submitting as issue

#  First, index the initial bam
picard BuildBamIndex --INPUT "${bam}"

bam_chrX="picard-test.chrX.bam"
bam_chrX_rmdup="picard-test.chrX.rmdup"

picard MarkDuplicates --INPUT "${bam_chrX}" --OUTPUT "${bam_chrX_rmdup}.bam" --METRICS_FILE "${bam_chrX_rmdup}.txt" --REMOVE_DUPLICATES true
samtools view -h "${bam}" chrX -b > "${bam_chrX}"

bam_chrX="picard-test.chrX.bam"
bam_chrX_fixmate="picard-test.chrX.fixmate.bam"

picard FixMateInformation --INPUT "${bam_chrX}" --OUTPUT "${bam_chrX_fixmate}" --ADD_MATE_CIGAR true
picard MarkDuplicatesWithMateCigar --INPUT "${bam_chrX}" --OUTPUT "${bam_chrX_rmdup}.bam" --METRICS_FILE "${bam_chrX_rmdup}.txt" --REMOVE_DUPLICATES true




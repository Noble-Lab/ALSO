#!/bin/bash

#  test_deduplication-issues_picard-MarkDuplicates.sh
#  KA

project="${1:-"test_deduplication-issues_picard-MarkDuplicates"}"
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

bam_markdup="${bam/bam/MarkDuplicates}"
picard MarkDuplicates --INPUT "${bam}" --OUTPUT "${bam_markdup}.bam" --METRICS_FILE "${bam_markdup}.txt"

bam_rmdup="${bam/bam/MarkDuplicates-remove}"
picard MarkDuplicates --INPUT "${bam}" --OUTPUT "${bam_rmdup}.bam" --METRICS_FILE "${bam_rmdup}.txt" --REMOVE_DUPLICATES true

picard BuildBamIndex --INPUT "${bam_rmdup}.bam"

samtools view -h "${bam_rmdup}.bam" chrX -b > "${bam_rmdup}.chrX.bam"
picard BuildBamIndex --INPUT "${bam_rmdup}.chrX.bam"

# bam_queryname="${bam/bam/sort-queryname.bam}"
# picard SortSam --INPUT "${bam}" --OUTPUT "${bam_queryname}" --SORT_ORDER queryname

# bam_rmdup="${bam_queryname/bam/MarkDuplicates-remove}"
# picard MarkDuplicates --INPUT "${bam_queryname}" --OUTPUT "${bam_rmdup}.bam" --METRICS_FILE "${bam_rmdup}.txt" --REMOVE_DUPLICATES true

# bam_coordinate="${bam_rmdup}.sort-queryname.bam"
# picard SortSam --INPUT "${bam_rmdup}.bam" --OUTPUT "${bam_coordinate}" --SORT_ORDER coordinate

# picard BuildBamIndex --INPUT "${bam_coordinate}"

# samtools view -h "${bam_coordinate}" chrX -b > "${bam_coordinate/bam/chrX.bam}"
# picard BuildBamIndex --INPUT "${bam_coordinate/bam/chrX.bam}"

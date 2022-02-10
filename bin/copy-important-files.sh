#!/bin/bash

# crossmap-bams.sh
# KA

script="copy-important-files.sh"

directory_from="/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
directory_to="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0"

cd "${directory_from}" || 
{
    echo "Directory not found; look into this."
    return 1 2> /dev/null
    exit 1
}

bam129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam"
bamCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam"
bamNmasked="mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam"
bai129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam.bai"
baiCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam.bai"
baiNmasked="mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam.bai"
s2p129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.sam2pairwise.munged.txt"
s2pCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.sam2pairwise.munged.txt"
s2pNmasked="mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.sam2pairwise.munged.txt"
liftOver129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.bam"
liftOverCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.bam"
liftOver129sorted="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.sorted.bam"
liftOver129sortedIndex="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.sorted.bam.bai"
liftOverCASTsorted="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.sorted.bam"
liftOverCASTsortedIndex="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.liftOverMm10.sorted.bam.bai"
Rdata5="test.PE-processing.part-5.Rdata"

cp \
"${bam129}" "${bamCAST}" "${bamNmasked}" \
"${bai129}" "${baiCAST}" "${baiNmasked}" \
"${s2p129}" "${s2pCAST}" "${s2pNmasked}" \
"${liftOver129}" "${liftOver129sorted}" "${liftOver129sortedIndex}" \
"${liftOverCAST}" "${liftOverCASTsorted}" "${liftOverCASTsortedIndex}" \
"${Rdata5}" \
"${directory_to}"

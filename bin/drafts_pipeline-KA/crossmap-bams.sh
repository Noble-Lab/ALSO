#!/bin/bash

# crossmap-bams.sh
# KA

script="crossmap-bams.sh"

#  Important dependency: sam2pairwise
#+ 
#+ Can install via conda
#+ conda install -c bioconda crossmap

#  Working with the following three files for now:
#+ - 129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam
#+ - CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam

command -v CrossMap.py &> /dev/null ||
{
    echo "Exiting: CrossMap.py not found. Install CrossMap.py."
    return 1 2> /dev/null
    exit 1
}

directory="/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
cd "${directory}" || 
{
    echo "Directory not found; look into this."
    return 1 2> /dev/null
    exit 1
}

directory_liftOver="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/liftOver"

chain129toMm10="${directory_liftOver}/129S1-SvImJ-to-mm10.over.chain.munged"
chainCASTtoMm10="${directory_liftOver}/CAST-EiJ-to-mm10.over.chain.munged"

orig129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam"
origCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.extendedCIGAR.bam"

liftOver129="${orig129/.bam/.liftOverMm10}"
liftOverCAST="${origCAST/.bam/.liftOverMm10}"

CrossMap.py bam -a "${chain129toMm10}" "${orig129}" "${liftOver129}" 2> "${liftOver129}.txt"
CrossMap.py bam -a "${chainCASTtoMm10}" "${origCAST}" "${liftOverCAST}" 2> "${liftOverCAST}.txt"

echo ""
samtools view -h "${liftOver129}.bam" | head -500
echo ""
samtools view -h "${liftOverCAST}.bam" | head -500

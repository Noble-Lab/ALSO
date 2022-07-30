#$ -S /bin/bash
#$ -R y
#$ -P sage
#$ -l h_rt=48:0:0,mfree=9G
#$ -pe serial 4
#$ -q sage-long.q
#$ -M kga0@uw.edu
#$ -cwd

#  align_Berletch-data_outSAMunmapped-None.sh
#  KA

#  Description:
#+ Align single-strand RNA-seq datasets from Berletch et al. using mm11 and
#+ SPRET-EiJ fasta and gtf files
#+
#+ Usage:
#+ The script has one dependency: STAR
#+
#+ Notes:
#+

module load STAR/2.7.6a

# threads=8
# memory=72000000000
threads=4
memory=36000000000
index="${1}"
file="${2}"
align_to="${3}"
gtf="${4}"
overhang="${5}"

#  Basic parameters
# STAR \
# --runThreadN "${threads}" \
# --limitGenomeGenerateRAM "${memory}" \
# --genomeDir "${index}" \
# --readFilesIn "${file}" \
# --readFilesCommand zcat \
# --outFileNamePrefix "./alignments/$(basename "${file%.fastq.gz}").${align_to}" \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes Standard \
# --sjdbGTFfile "${gtf}" \
# --sjdbOverhang "${overhang}"

#  Parameters used by He Fang
STAR \
--limitGenomeGenerateRAM "${memory}" \
--limitBAMsortRAM 0 \
--runThreadN "${threads}" \
--genomeDir "${index}" \
--readFilesIn "${file}" \
--readFilesCommand zcat \
--outFileNamePrefix "./alignments_None/$(basename "${file%.fastq.gz}").${align_to}" \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 5000000 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--outSAMstrandField intronMotif \
--outSAMunmapped None \
--outSAMattributes NH HI NM MD AS \
--outSAMtype BAM SortedByCoordinate \
--outSAMheaderHD @HD VN:1.4 \
--sjdbScore 2 \
--sjdbGTFfile "${gtf}" \
--sjdbOverhang "${overhang}" \
--limitOutSJcollapsed 5000000 \
--quantMode GeneCounts \
--alignEndsType EndToEnd

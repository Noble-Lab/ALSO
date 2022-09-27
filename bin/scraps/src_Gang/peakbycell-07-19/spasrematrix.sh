#!/bin/bash
#$ -l h_vmem=20G
#$ -l h_rt=32:00:00
#$ -M gangliuw@uw.edu
#$ -m ae
#$ -j y
#$ -N sparse_matrix
#$ -V
#$ -cwd

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-07-19-CellByPeak
module add samtools/1.14
module add python/3.7.7
module add bedtools/2.29.2
#module add bwa/0.7.17
#module load python/3.7.7 numpy/1.21.1 bx-python/0.8.13 tbb/2020_U2 bowtie2/2.4.4 samtools/1.14 pcre2/10.39 R/4.1.2 scipy/1.7.0 pandas/1.3.1 cython/0.29.17 joblib/0.15.1 scikit-learn/0.24.2 iced/0.5.10 hic-pro/3.1.0

sample_index=1
#bam_input="bam-ALSO/Disteche_sample_"${sample_index}"_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_"${sample_index}".dedup.corrected.corrected.CAST.CAST.bam"
#samtools index bam-ALSO/Disteche_sample_1_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_1.dedup.corrected.corrected.CAST.CAST.bam 

# python bin/sciatac_pipeline/src/get_unique_fragments.py \
# --output_bam results/test/OUTPUT_BAM \
# --fragments results/test/FRAGMENTS \
# --transposition_sites_bed results/test/TRANSPOSITION_SITES_BED \
# --duplicate_read_counts results/test/DUPLICATE_READ_COUNTS \
#   bam-ALSO/Disteche_sample_1_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_1.dedup.corrected.corrected.CAST.CAST.bam 

# results/merged.cellwhitelist.txt

# mkdir results/2022-09-27-sample-specific-peak
# mkdir results/2022-09-27-sample-specific-peak/CAST
# mkdir results/2022-09-27-sample-specific-peak/ambiguous
# mkdir results/2022-09-27-sample-specific-peak/mm10

# dir="results/2022-09-27-sample-specific-peak/CAST/"

# echo "Cast Fragments Generating..."

# python bin/sciatac_pipeline/src/get_unique_fragments.py \
# --output_bam ${dir}OUTPUT_BAM \
# --fragments ${dir}FRAGMENTS \
# --transposition_sites_bed ${dir}TRANSPOSITION_SITES_BED \
# --duplicate_read_counts ${dir}DUPLICATE_READ_COUNTS \
#   bam-ALSO/Disteche_sample_${sample_index}_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_${sample_index}.dedup.corrected.corrected.CAST.CAST.bam

# echo "Cast Fragments Done."

# dir="results/2022-09-27-sample-specific-peak/mm10/"

# echo "mm10 Fragments Generating..."

# python bin/sciatac_pipeline/src/get_unique_fragments.py \
# --output_bam ${dir}OUTPUT_BAM \
# --fragments ${dir}FRAGMENTS \
# --transposition_sites_bed ${dir}TRANSPOSITION_SITES_BED \
# --duplicate_read_counts ${dir}DUPLICATE_READ_COUNTS \
#   bam-ALSO/Disteche_sample_${sample_index}_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_${sample_index}.dedup.corrected.corrected.mm10.mm10.bam

# echo "mm10 Fragments Done."

# dir="results/2022-09-27-sample-specific-peak/ambiguous/"

# echo "ambiguous Fragments Generating..."

# python bin/sciatac_pipeline/src/get_unique_fragments.py \
# --output_bam ${dir}OUTPUT_BAM \
# --fragments ${dir}FRAGMENTS \
# --transposition_sites_bed ${dir}TRANSPOSITION_SITES_BED \
# --duplicate_read_counts ${dir}DUPLICATE_READ_COUNTS \
#   bam-ALSO/Disteche_sample_${sample_index}_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_${sample_index}.dedup.corrected.corrected.mm10.ambiguous.bam

# echo "ambiguous Fragments Done."

#zcat bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz > results/2022-09-27-sample-specific-peak/Disteche_sample_1_peaks.narrowPeak.bed

## CAST

bedtools intersect -sorted -a results/2022-09-27-sample-specific-peak/Disteche_sample_1_peaks.narrowPeak.bed -b results/2022-09-27-sample-specific-peak/CAST/TRANSPOSITION_SITES_BED -wa -wb > results/2022-09-27-sample-specific-peak/CAST/transposition_sites_intersect

python bin/sciatac_pipeline/src/generate_sparse_matrix.py \
--transposition_sites_intersect results/2022-09-27-sample-specific-peak/CAST/transposition_sites_intersect \
--intervals results/2022-09-27-sample-specific-peak/Disteche_sample_1_peaks.narrowPeak.bed \
--cell_whitelist results/merged.cellwhitelist.txt \
--matrix_output results/2022-09-27-sample-specific-peak/CAST/MATRIX_OUTPUT_merged_peaks_QCed_cells_${sample_index}

echo "Cast Matrix Done."

## mm10

bedtools intersect -sorted -a results/2022-09-27-sample-specific-peak/Disteche_sample_1_peaks.narrowPeak.bed -b results/2022-09-27-sample-specific-peak/mm10/TRANSPOSITION_SITES_BED -wa -wb > results/2022-09-27-sample-specific-peak/mm10/transposition_sites_intersect

python bin/sciatac_pipeline/src/generate_sparse_matrix.py \
--transposition_sites_intersect results/2022-09-27-sample-specific-peak/mm10/transposition_sites_intersect \
--intervals bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz \
--cell_whitelist results/merged.cellwhitelist.txt \
--matrix_output results/2022-09-27-sample-specific-peak/mm10/MATRIX_OUTPUT_merged_peaks_QCed_cells_${sample_index}

echo "mm10 Matrix Done."

## ambigous
bedtools intersect -sorted -a results/2022-09-27-sample-specific-peak/Disteche_sample_1_peaks.narrowPeak.bed -b results/2022-09-27-sample-specific-peak/ambiguous/TRANSPOSITION_SITES_BED -wa -wb > results/2022-09-27-sample-specific-peak/ambiguous/transposition_sites_intersect


python bin/sciatac_pipeline/src/generate_sparse_matrix.py \
--transposition_sites_intersect results/2022-09-27-sample-specific-peak/ambiguous/transposition_sites_intersect \
--intervals bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz \
--cell_whitelist results/merged.cellwhitelist.txt \
--matrix_output results/2022-09-27-sample-specific-peak/ambiguous/MATRIX_OUTPUT_merged_peaks_QCed_cells_${sample_index}

echo "ambiguous Matrix Done."




# END

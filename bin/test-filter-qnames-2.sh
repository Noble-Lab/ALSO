#!/bin/bash

#  test-filter-qnames-2.sh
#  KA


#  Set up variables, home directory, functions --------------------------------
# shellcheck disable=1091
# shellcheck disable=SC2002

dir_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
dir_data="${dir_base}/data/files_bam"
f_mm10="${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm.bam"
f_CAST="${dir_data}/Disteche_sample_13.dedup.CAST.sort-c.rm.bam"

cd "${dir_base}" || 
    {
        echo "Exiting: cd failed. Check on this."
        exit 1
    }

. ./bin/auxiliary/functions_preprocessing.sh

#  Make test bam files --------------------------------------------------------
split_bam_chromosome 4 "${f_mm10}" "chr19"
index_bam 4 "${f_mm10/.bam/.chr19.bam}"

split_bam_chromosome 4 "${f_CAST}" "chr19"
index_bam 4 "${f_CAST/.bam/.chr19.bam}"


#  Run script -----------------------------------------------------------------
Rscript ./bin/generate-qname-lists.R \
--bam "${f_mm10}" \
--bai "${f_mm10/.bam/.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique TRUE \
--tally TRUE \
--remove TRUE
# Lines in unmated.txt.gz: 112
# Lines in ambiguous.txt.gz: 4
# Lines in duplicated.txt.gz: 92
# Lines in singleton.txt.gz: 108

Rscript ./bin/generate-qname-lists.R \
--bam "${f_CAST}" \
--bai "${f_CAST/.bam/.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique TRUE \
--tally TRUE \
--remove TRUE
# Lines in unmated.txt.gz: 114
# Lines in ambiguous.txt.gz: 2
# Lines in duplicated.txt.gz: 78
# Lines in singleton.txt.gz: 110

combine_gz_qname_lists_return_unique_gzip \
"${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm."{ambiguous,duplicated,singleton,unmated}".txt.gz" \
> "${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm.to-exclude.txt.gz"

combine_gz_qname_lists_return_unique_gzip \
"${dir_data}/Disteche_sample_13.dedup.CAST.sort-c.rm."{ambiguous,duplicated,singleton,unmated}".txt.gz" \
> "${dir_data}/Disteche_sample_13.dedup.CAST.sort-c.rm.to-exclude.txt.gz"

gzcat "${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm.to-exclude.txt.gz" | wc -l  # 136
gzcat "${dir_data}/Disteche_sample_13.dedup.CAST.sort-c.rm.to-exclude.txt.gz" | wc -l  # 134

exclude_qname_reads_picard \
"${f_mm10}" \
"${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm.to-exclude.txt.gz" \
"${f_mm10/.bam/.corrected.bam}"

exclude_qname_reads_picard \
"${f_CAST}" \
"${dir_data}/Disteche_sample_13.dedup.CAST.sort-c.rm.to-exclude.txt.gz" \
"${f_CAST/.bam/.corrected.bam}"

run_flagstat 2 "${f_mm10}"
run_flagstat 2 "${f_mm10/.bam/.corrected.bam}"

run_flagstat 2 "${f_CAST}"
run_flagstat 2 "${f_CAST/.bam/.corrected.bam}"


head -25 "${f_mm10/.bam/.flagstat.txt}"  # Before correction
# 85513712 + 0 in total (QC-passed reads + QC-failed reads)
# 85513712 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513712 + 0 mapped (100.00% : N/A)
# 85513712 + 0 primary mapped (100.00% : N/A)
# 85513712 + 0 paired in sequencing
# 42756850 + 0 read1
# 42756862 + 0 read2
# 85513712 + 0 properly paired (100.00% : N/A)
# 85513712 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

head -25 "${f_mm10/.bam/.corrected.flagstat.txt}"  # After correction
# 85513496 + 0 in total (QC-passed reads + QC-failed reads)
# 85513496 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 85513496 + 0 mapped (100.00% : N/A)
# 85513496 + 0 primary mapped (100.00% : N/A)
# 85513496 + 0 paired in sequencing
# 42756748 + 0 read1
# 42756748 + 0 read2
# 85513496 + 0 properly paired (100.00% : N/A)
# 85513496 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

head -25 "${f_CAST/.bam/.flagstat.txt}"  # Before correction
# 80638866 + 0 in total (QC-passed reads + QC-failed reads)
# 80638866 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 80638866 + 0 mapped (100.00% : N/A)
# 80638866 + 0 primary mapped (100.00% : N/A)
# 80638866 + 0 paired in sequencing
# 40319425 + 0 read1
# 40319441 + 0 read2
# 80638866 + 0 properly paired (100.00% : N/A)
# 80638866 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

head -25 "${f_CAST/.bam/.corrected.flagstat.txt}"  # After correction
# 80638660 + 0 in total (QC-passed reads + QC-failed reads)
# 80638660 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 80638660 + 0 mapped (100.00% : N/A)
# 80638660 + 0 primary mapped (100.00% : N/A)
# 80638660 + 0 paired in sequencing
# 40319330 + 0 read1
# 40319330 + 0 read2
# 80638660 + 0 properly paired (100.00% : N/A)
# 80638660 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


#  mm10
echo $(( 85513712 - 85513496 ))  # 216

#  CAST
echo $(( 80638866 - 80638660 ))  # 206


# Did --unique TRUE work properly? --------------------------------------------
for i in "${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm."{ambiguous,duplicated,singleton,unmated}".txt.gz"; do
    decompress_gzip "${i}"
done

head -10 "${dir_data}/Disteche_sample_13.dedup.mm10.sort-c.rm.duplicated.txt"
# CCGCGGAGACATGCTAACCTTATCAAGGTCCCATCCTAGT:1153827938
# CCGCGGAGACATGCTAACCTTATCAAGGTCCCATCCTAGT:1153827938
# CCGCGGAGACATGCTAACCTTATCAAGGTCCCATCCTAGT:1153827938
# CCGCGGAGACATGCTAACCTTATCAAGGTCCCATCCTAGT:1153827938
# CCTTCCTTCATTGGTTCTCAACGTACTAGGAAGGTACTAA:1517312853
# CCTTCCTTCATTGGTTCTCAACGTACTAGGAAGGTACTAA:1517312853
# CCTTCCTTCATTGGTTCTCAACGTACTAGGAAGGTACTAA:1517312853
# CCTTCCTTCATTGGTTCTCAACGTACTAGGAAGGTACTAA:1517312853
# ATAGGAGTTCGATCATGATATGGACCAAGGAAGGTACTAA:1347465485
# ATAGGAGTTCGATCATGATATGGACCAAGGAAGGTACTAA:1347465485

#  --unique TRUE did not work properly...


#  HPC work -------------------------------------------------------------------
gangliuw_CAST
cp Disteche_sample_20.dedup.bam /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0501_test-preprocessing-module
cd /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0501_test-preprocessing-module || ! echo "cd failed"
mv Disteche_sample_20.dedup.bam Disteche_sample_20.dedup.CAST.bam

gangliuw_mm10
cp Disteche_sample_20.dedup.bam /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0501_test-preprocessing-module
cd /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0501_test-preprocessing-module || ! echo "cd failed"
mv Disteche_sample_20.dedup.bam Disteche_sample_20.dedup.mm10.bam

. ./bin/auxiliary/functions_preprocessing.sh

dir_data="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0501_test-preprocessing-module"
f_CAST="${dir_data}/Disteche_sample_20.dedup.CAST.bam"
f_mm10="${dir_data}/Disteche_sample_20.dedup.mm10.bam"

mv "${f_CAST}" "${TMPDIR}" && cd "${TMPDIR}"
mv "${f_mm10}" "${TMPDIR}" && cd "${TMPDIR}"

dir_data="${TMPDIR}"
f_CAST="${TMPDIR}/Disteche_sample_20.dedup.CAST.bam"
f_mm10="${TMPDIR}/Disteche_sample_20.dedup.mm10.bam"

sort_bam_coordinate_samtools_auto 1 "${f_CAST}"
index_bam 1 "${f_CAST/.bam/.sort-c.bam}"
Rscript /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/bin/generate-qname-lists.R \
--bam "${f_CAST/.bam/.sort-c.bam}" \
--bai "${f_CAST/.bam/.sort-c.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique TRUE \
--tally TRUE \
--remove TRUE

sort_bam_coordinate_samtools_auto 1 "${f_mm10}"
index_bam 1 "${f_mm10/.bam/.sort-c.bam}"
Rscript /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/bin/generate-qname-lists.R \
--bam "${f_mm10/.bam/.sort-c.bam}" \
--bai "${f_mm10/.bam/.sort-c.bam.bai}" \
--outdir "${dir_data}" \
--chunk 100000 \
--mated FALSE \
--unmated TRUE \
--ambiguous TRUE \
--trans TRUE \
--duplicated TRUE \
--singleton TRUE \
--unique TRUE \
--tally TRUE \
--remove TRUE

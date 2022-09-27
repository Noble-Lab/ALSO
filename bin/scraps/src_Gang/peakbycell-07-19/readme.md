### Generate cell by peak matrix 

ssh grid
qlogin 
#-l 
#mfree=12G -pe serial 4 
#,cuda=1,gpgpu=TRUE

# cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO
# cd ../2022-07-19-CellByPeak
cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-07-19-CellByPeak
module add samtools/1.14
module add python/3.7.7
module add bedtools/2.29.2

######################
## set up to vol6 ####
######################

df -kh /net/noble/vol6

df -h |grep 'noble'

# cd /net/noble/vol6
mkdir /net/noble/vol6/user/gangliuw/mouse-cross/results/2022-07-19


ln -s /net/noble/vol6/user/gangliuw/mouse-cross/results/2022-07-19 results

ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/ bin

# input bam
ln -s /net/noble/vol6/user/gangliuw/mouse-cross/results/2022-06-06 bam-ALSO

# input bed
mkdir bed
ln -s /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/analysis_mouse_mm10_Cast_Nmasked/call_peaks bed/call_peaks
<!-- # from CX: the first one includes E9.5 - E13.5, and the second one includes P1
/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/reanalyze_embryo/call_peaks
/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/reanalyze_P1/call_peaks -->

head bed/call_peaks/merged_peaks.bed 
<!-- 
chr1    3094657 3095334
chr1    3119438 3120970
chr1    3121319 3122087 
-->

## try sciatac_pipeline
python bin/sciatac_pipeline/src/generate_sparse_matrix.py

<!-- usage: Script to take BAM and BED intervals and generate a sparse peak by cell matrix. Note that if the BED file of intervals has 4 or more columns, the fourth column is used as the feature and the score will be the aggregate of all intervals belonging to that group.
       [-h] --transposition_sites_intersect [TRANSPOSITION_SITES_INTERSECT]
       --intervals INTERVALS --cell_whitelist CELL_WHITELIST
       [--matrix_output MATRIX_OUTPUT]
Script to take BAM and BED intervals and generate a sparse peak by cell matrix. Note that if the BED file of intervals has 4 or more columns, the fourth column is used as the feature and the score will be the aggregate of all intervals belonging to that group.: error: the following arguments are required: --transposition_sites_intersect, --intervals, --cell_whitelist -->

samtools view -h --no-header bam-ALSO/Disteche_sample_11_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_11.dedup.corrected.corrected.CAST.CAST.bam | cut -f1 | head -100 > cell_whitelist.txt

<!-- Traceback (most recent call last):
  File "bin/sciatac_pipeline/src/generate_sparse_matrix.py", line 96, in <module>
    features, interval_columns = load_feature_whitelist(open_file(args.intervals))
  File "bin/sciatac_pipeline/src/generate_sparse_matrix.py", line 36, in load_feature_whitelist
    for line in file_object:
  File "/net/gs/vol3/software/modules-sw/python/3.7.7/Linux/CentOS7/x86_64/lib/python3.7/codecs.py", line 322, in decode
    (result, consumed) = self._buffer_decode(data, self.errors, final)
UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte -->





<!-- debugfs -w /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-07-19-CellByPeak/bam-ALSO -->

# 07.26
# run get_unique_fragments by inputing the bam from ALSO pipeline to generate fragments files. 
# Then perform generate_sparse_matrix, which take both fragments files and bed files (peak information), to generate cell by peak matrix. 

ls bin/sciatac_pipeline/src/
python -h bin/sciatac_pipeline/src/get_unique_fragments.py

<!-- usage: Script to deduplicate position sorted sciATAC BAM file.
       [-h] [--output_bam OUTPUT_BAM] [--fragments FRAGMENTS]
       [--transposition_sites_bed TRANSPOSITION_SITES_BED]
       [--duplicate_read_counts DUPLICATE_READ_COUNTS]
       [--insert_sizes INSERT_SIZES] [--flanking_distance FLANKING_DISTANCE]
       input_bam -->

mkdir results/test

samtools index bam-ALSO/Disteche_sample_1_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_1.dedup.corrected.corrected.CAST.CAST.bam 


python bin/sciatac_pipeline/src/get_unique_fragments.py \
--output_bam results/test/OUTPUT_BAM \
--fragments results/test/FRAGMENTS \
--transposition_sites_bed results/test/TRANSPOSITION_SITES_BED \
--duplicate_read_counts results/test/DUPLICATE_READ_COUNTS \
  bam-ALSO/Disteche_sample_1_ALSO/04_run_filter-qnames-by-assignment/Disteche_sample_1.dedup.corrected.corrected.CAST.CAST.bam 

qsub fragments.sh


ls -lah results/test/

head results/test/FRAGMENTS
# chr1    1048    1234    TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG        1
# chr1    1049    1245    ATAGTATTGCTCCATACCTGAGTCGCGTCGAAGTTCGCTG        1
# chr1    1063    1216    TGAATTATGAAGCTCATTCGAGTCATGGATCGTTACGTTG        1

wc -l results/test/FRAGMENTS
# 12273720

head results/test/TRANSPOSITION_SITES_BED
# chr1    1047    1049    TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG
# chr1    1048    1050    ATAGTATTGCTCCATACCTGAGTCGCGTCGAAGTTCGCTG
# chr1    1062    1064    TGAATTATGAAGCTCATTCGAGTCATGGATCGTTACGTTG

wc -l results/test/TRANSPOSITION_SITES_BED
# 24547440

head cell_whitelist.txt
# CGAATAGTTGTCAGCGCAATTAAGCTCGCCGAGACGTACG:458694603
# ATATGATGAATTGACTTCGTTAAGCTCGCCACTAGTTGCA:727747625

cut -f4 results/test/FRAGMENTS | uniq | wc -l 
# 12184613
cut -f4 results/test/FRAGMENTS | sort | uniq | wc -l
# 395044
cut -f4 results/test/FRAGMENTS | sort | uniq > results/test/cell_whitelist.txt

head results/test/cell_whitelist.txt
wc -l results/test/cell_whitelist.txt
# 395044
 
# CX:
# bedtools intersect -sorted -a {intervals} -b {transposition_sites_bed} -wa -wb
bedtools intersect -sorted -a bed/call_peaks/merged_peaks.bed -b results/test/TRANSPOSITION_SITES_BED -wa -wb > results/test/transposition_sites_intersect

head results/test/transposition_sites_intersect
<!-- chr1    1048    1234    TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG        1       chr1    1047    1049   TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG
chr1    1048    1234    TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG        1       chr1    1048    1050   ATAGTATTGCTCCATACCTGAGTCGCGTCGAAGTTCGCTG -->

wc -l results/test/transposition_sites_intersect
# 526396544

# transposition_sites_bed
# transposition_sites_intersect

python bin/sciatac_pipeline/src/generate_sparse_matrix.py \
--transposition_sites_intersect results/test/transposition_sites_intersect \
--intervals  bed/call_peaks/merged_peaks.bed \  #results/test/FRAGMENTS \
--cell_whitelist results/test/cell_whitelist.txt \
--matrix_output results/test/MATRIX_OUTPUT

<!-- Building matrix...
Traceback (most recent call last):
  File "bin/sciatac_pipeline/src/generate_sparse_matrix.py", line 99, in <module>
    feature_indices, cell_indices, values = get_feature_cell_counts(args.transposition_sites_intersect, features, cell_whitelist, interval_columns)
  File "bin/sciatac_pipeline/src/generate_sparse_matrix.py", line 61, in get_feature_cell_counts
    cell = entries[interval_columns + 3].split(":")[0]
IndexError: list index out of range -->
# solved by bedtools intersect

qsub sparsematrix.sh

ls results/test/

head results/test/MATRIX_OUTPUT.rows.txt
<!-- TCCAGCTCTCTCCATACCTGATTGAGGATACGTTACGTTG
ATAGTATTGCTCCATACCTGAGTCGCGTCGAAGTTCGCTG -->
head results/test/MATRIX_OUTPUT.columns.txt
<!-- AACCAATAACAAGTCAACGGATCGAATCTGCAATCTCGGA
AACCAATAACAAGTCAACGGATCGAATCTGCTCTGCGATC -->


python

from scipy.io import mmread

data = mmread('results/test/MATRIX_OUTPUT.mtx')

print(data)
<!-- (41277, 0)    2
  (1, 1)        2
  (9, 1)        1
  (271, 1)      2 -->


data
<!-- <395044x395044 sparse matrix of type '<class 'numpy.int64'>'
        with 31305057 stored elements in COOrdinate format> -->

data.todense()[1:5,1:5]


# CX: the input of intervals is a peak list
# also need to fix the bed file, since it used the peaks to generate the bed file.

# bedtools intersect -sorted -a {intervals} -b {transposition_sites_bed} -wa -wb
bedtools intersect -sorted -a bed/call_peaks/merged_peaks.bed -b results/test/TRANSPOSITION_SITES_BED -wa -wb > results/test/transposition_sites_intersect

<!-- (base) [gangliuw@n014 2022-07-19-CellByPeak]$ head results/test/transposition_sites_intersect
chr1    3094657 3095334 chr1    3094751 3094753 CTCGGTTGCAATCGGCTTGAGAAGAGTATTTGATTCTCGT
chr1    3094657 3095334 chr1    3094904 3094906 GCGCAGGCGAATCGGCTTGAGAATGAATAAAGCCGAAGTA -->


qsub sparsematrix.sh

head results/test/MATRIX_OUTPUT_merged_peaks.columns.txt
<!-- (base) [gangliuw@nexus1 2022-07-19-CellByPeak]$ head results/test/MATRIX_OUTPUT_merged_peaks.columns.txt
AACCAATAACAAGTCAACGGATCGAATCTGCAATCTCGGA
AACCAATAACAAGTCAACGGATCGAATCTGCTCTGCGATC
(base) [gangliuw@nexus1 2022-07-19-CellByPeak]$ head results/test/MATRIX_OUTPUT_merged_peaks.rows.txt
chr1_3094657_3095334
chr1_3119438_3120970
chr1_3121319_3122087 -->
head results/test/MATRIX_OUTPUT_merged_peaks.rows.txt
python

from scipy.io import mmread

data = mmread('results/test/MATRIX_OUTPUT_merged_peaks.mtx')

print(data)

  <!-- (247964, 2)   1
  (61579, 11)   2
  (190455, 23)  2
   -->
data

<!-- <287561x395044 sparse matrix of type '<class 'numpy.int64'>'
        with 1386610 stored elements in COOrdinate format> -->

# for mm10 08.09
mkdir results/test1/
ls results/test1/
qsub fragments.sh

bedtools intersect -sorted -a bed/call_peaks/merged_peaks.bed -b results/test1/TRANSPOSITION_SITES_BED -wa -wb > results/test1/transposition_sites_intersect

cut -f4 results/test1/FRAGMENTS | sort | uniq -c > results/test1/cell_fragment_counts.txt
cut -f4 results/test1/FRAGMENTS | sort | uniq > results/test1/cell_whitelist.txt


qsub sparsematrix.sh

ls results/test1/
wc -l results/test1/MATRIX_OUTPUT_merged_peaks.columns.txt # 506543
wc -l results/test1/MATRIX_OUTPUT_merged_peaks.rows.txt # 287561

python

from scipy.io import mmread

data = mmread('results/test1/MATRIX_OUTPUT_merged_peaks.mtx')

print(data)
data
<!-- 
<287561x506543 sparse matrix of type '<class 'numpy.int64'>'
        with 2860102 stored elements in COOrdinate format>
         -->


# ambiguous bam
mkdir results/test2/
ls results/test2/
qsub fragments.sh

bedtools intersect -sorted -a bed/call_peaks/merged_peaks.bed -b results/test2/TRANSPOSITION_SITES_BED -wa -wb > results/test2/transposition_sites_intersect

cut -f4 results/test2/FRAGMENTS | sort | uniq -c > results/test2/cell_fragment_counts.txt
cut -f4 results/test2/FRAGMENTS | sort | uniq > results/test2/cell_whitelist.txt

qsub sparsematrix.sh

python

from scipy.io import mmread

data = mmread('results/test2/MATRIX_OUTPUT_merged_peaks.mtx')

print(data)
data
<!-- <287561x597202 sparse matrix of type '<class 'numpy.int64'>'
        with 4427160 stored elements in COOrdinate format> -->
## knee plot
cut -f4 results/test/FRAGMENTS | sort | uniq -c > results/test/cell_fragment_counts.txt

# 08.23
## simple threshold based knee plot
# union the cell list
ls
qsub sparsematrix_QCed_cell.sh

# AH's whitelist
wc -l /net/gs/vol1/home/cxqiu/share/Gang/Disteche_sample_*.cell_whitelist.txt > results/AH_QC_cells.counts.txt

#2627 /net/gs/vol1/home/cxqiu/share/Gang/Disteche_sample_1.cell_whitelist.txt
ln -s /net/gs/vol1/home/cxqiu/share/Gang/ Cell_Whitelist_AH


# 09.27
# try sample specific peak
qsub sparsematrix_QCed_cell_sample-specific-peak.sh

ls bed/call_peaks/Disteche_sample_1_*
<!-- bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz
bed/call_peaks/Disteche_sample_1_peaks.xls
bed/call_peaks/Disteche_sample_1_summits.bed -->

cat bed/call_peaks/Disteche_sample_1_peaks.xls | wc -l
<!-- 45452 -->
wc -l bed/call_peaks/merged_peaks.bed 
<!-- 287561 bed/call_peaks/merged_peaks.bed -->

zcat bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz | wc -l 
<!-- 45428 -->
zcat bed/call_peaks/Disteche_sample_1_peaks.narrowPeak.gz | head -10
<!-- chr1    3670689 3671130
chr1    3671621 3672654
chr1    3671621 3672654 -->

head -30 bed/call_peaks/Disteche_sample_1_peaks.xls
<!-- # This file is generated by MACS version 2.2.7.1
# Command line: callpeak -t /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/analysis_mouse_mm10_Cast_Nmasked/get_unique_fragments/Disteche_sample_1.bed.gz -f BED -g mm --nomodel --shift -100 --extsize 200 --keep-dup all --call-summits -n Disteche_sample_1 --outdir /net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/analysis_mouse_mm10_Cast_Nmasked/call_peaks
# ARGUMENTS LIST:
# name = Disteche_sample_1
# format = BED
# ChIP-seq file = ['/net/shendure/vol10/projects/cxqiu/nobackup/Disteche_CastB6_cross/Novaseq/analysis_mouse_mm10_Cast_Nmasked/get_unique_fragments/Disteche_sample_1.bed.gz']
# control file = None
# effective genome size = 1.87e+09
# band width = 300
# model fold = [5, 50]
# qvalue cutoff = 5.00e-02
# The maximum gap between significant sites is assigned as the read length/tag size.
# The minimum length of peaks is assigned as the predicted fragment length "d".
# Larger dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 10000 bps
# Broad region calling is off
# Paired-End mode is off
# Searching for subpeak summits is on

# tag size is determined as 2 bps
# total tags in treatment: 80382286
# Sequencing ends will be shifted towards 5' by 100 bp(s)
# d = 200
chr     start   end     length  abs_summit      pileup  -log10(pvalue)  fold_enrichment       -log10(qvalue)  name
chr1    3670690 3671130 441     3670841 52      17.6875 4.10853 14.7825 Disteche_sample_1_peak_1
chr1    3671622 3672654 1033    3671799 60      22.5021 4.59337 19.4712 Disteche_sample_1_peak_2a
chr1    3671622 3672654 1033    3672263 53      18.0103 4.11585 15.0961 Disteche_sample_1_peak_2b -->



# END
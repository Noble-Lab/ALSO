#!/bin/bash

# qlogin -l mfree=4G  # Too little
# qlogin -l mfree=8G  # Too little
# qlogin -l mfree=64G  # Too little
# qlogin -l mfree=20G  # Too little
qlogin -l mfree=80G

4dn-home
pipeline-test_env
module load gzip/1.12
module load pigz/2.3
module load parallel/20200922

#  Try running the Python script again...
dir_base="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross"
dir_data="${dir_base}/data/files_bam_test"
dir_results="${dir_base}/results"
dir_experiment="${dir_results}/2022-0615_test_driver_allelic-segregation_Disteche-sample-11_fix-Picard-memory-issue"
dir_experiment_sub="${dir_experiment}/results/Disteche_sample_11"
dir_experiment_out="${dir_experiment_sub}/06_run_exclude-reads-from-bam"

cd "${dir_experiment_sub}" || echo "cd failed. Check on this."
mkdir -p "${dir_experiment_out}"

#  ( Y ) Running
#  (   ) Completed
bam_in="${dir_data}/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"  # CAST
qnames="04_run_filter-qnames-by-assignment/Disteche_sample_11.mm10-combined.CAST-combined.rv-srt.txt.gz"  # keep only ambiguous
bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.ambiguous.bam"  # CAST bam, ambiguous assignments
\time --verbose \
python exclude-reads-from-bam.py \
--bam_in "${bam_in}" \
--qnames "${qnames}" \
--bam_out "${bam_out}"

#  ( Y ) Running
#  (   ) Completed
bam_in="${dir_data}/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"  # CAST
qnames="04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.mm10-combined.rv-srt.txt.gz"  # keep only CAST
bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.CAST.bam"  # CAST bam, CAST assignments
\time --verbose \
python exclude-reads-from-bam.py \
--bam_in "${bam_in}" \
--qnames "${qnames}" \
--bam_out "${bam_out}"

#  ( Y ) Running
#  (   ) Completed
bam_in="${dir_data}/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"  # CAST
qnames="04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.CAST-combined.rv-srt.txt.gz"  # keep only mm10
bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.mm10.bam"  # CAST bam, mm10 assignments
\time --verbose \
python exclude-reads-from-bam.py \
--bam_in "${bam_in}" \
--qnames "${qnames}" \
--bam_out "${bam_out}"


#  For mfree=64G (CAST bam, ambiguous assignments)
# Command terminated by signal 7
#         Command being timed: "python exclude-reads-from-bam.py --bam_in /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/files_bam_test/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam --qnames 04_run_filter-qnames-by-assignment/Disteche_sample_11.mm10-combined.CAST-combined.rv-srt.txt.gz --bam_out /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0615_test_driver_allelic-segregation_Disteche-sample-11_fix-Picard-memory-issue/results/Disteche_sample_11/06_run_exclude-reads-from-bam/Disteche_sample_11.dedup.corrected.corrected.CAST.ambiguous.bam"
#         User time (seconds): 284.51
#         System time (seconds): 139.44
#         Percent of CPU this job got: 46%
#         Elapsed (wall clock) time (h:mm:ss or m:ss): 15:02.67
#         Average shared text size (kbytes): 0
#         Average unshared data size (kbytes): 0
#         Average stack size (kbytes): 0
#         Average total size (kbytes): 0
#         Maximum resident set size (kbytes): 67105356
#         Average resident set size (kbytes): 0
#         Major (requiring I/O) page faults: 0
#         Minor (reclaiming a frame) page faults: 15627685
#         Voluntary context switches: 69373
#         Involuntary context switches: 7844
#         Swaps: 0
#         File system inputs: 0
#         File system outputs: 0
#         Socket messages sent: 0
#         Socket messages received: 0
#         Signals delivered: 0
#         Page size (bytes): 4096
#         Exit status: 0

#  For mfree=20G (CAST bam, mm10 assignments)
# Started: Running 'ids = list(set(...))'
# Command terminated by signal 7
#         Command being timed: "python exclude-reads-from-bam.py --bam_in /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/files_bam_test/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam --qnames 04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.CAST-combined.rv-srt.txt.gz --bam_out /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0615_test_driver_allelic-segregation_Disteche-sample-11_fix-Picard-memory-issue/results/Disteche_sample_11/06_run_exclude-reads-from-bam/Disteche_sample_11.dedup.corrected.corrected.CAST.mm10.bam"
#         User time (seconds): 144.97
#         System time (seconds): 106.63
#         Percent of CPU this job got: 58%
#         Elapsed (wall clock) time (h:mm:ss or m:ss): 7:08.30
#         Average shared text size (kbytes): 0
#         Average unshared data size (kbytes): 0
#         Average stack size (kbytes): 0
#         Average total size (kbytes): 0
#         Maximum resident set size (kbytes): 20967996
#         Average resident set size (kbytes): 0
#         Major (requiring I/O) page faults: 0
#         Minor (reclaiming a frame) page faults: 8462610
#         Voluntary context switches: 24011
#         Involuntary context switches: 14192654
#         Swaps: 0
#         File system inputs: 0
#         File system outputs: 0
#         Socket messages sent: 0
#         Socket messages received: 0
#         Signals delivered: 0
#         Page size (bytes): 4096
#         Exit status: 0

#  For mfree=20G (CAST bam, CAST assignments)
# Command terminated by signal 7
#         Command being timed: "python exclude-reads-from-bam.py --bam_in /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/files_bam_test/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam --qnames 04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.mm10-combined.rv-srt.txt.gz --bam_out /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2022-0615_test_driver_allelic-segregation_Disteche-sample-11_fix-Picard-memory-issue/results/Disteche_sample_11/06_run_exclude-reads-from-bam/Disteche_sample_11.dedup.corrected.corrected.CAST.CAST.bam"
#         User time (seconds): 136.96
#         System time (seconds): 29.81
#         Percent of CPU this job got: 51%
#         Elapsed (wall clock) time (h:mm:ss or m:ss): 5:26.67
#         Average shared text size (kbytes): 0
#         Average unshared data size (kbytes): 0
#         Average stack size (kbytes): 0
#         Average total size (kbytes): 0
#         Maximum resident set size (kbytes): 20968020
#         Average resident set size (kbytes): 0
#         Major (requiring I/O) page faults: 0
#         Minor (reclaiming a frame) page faults: 4896795
#         Voluntary context switches: 23294
#         Involuntary context switches: 11687
#         Swaps: 0
#         File system inputs: 0
#         File system outputs: 8
#         Socket messages sent: 0
#         Socket messages received: 0
#         Signals delivered: 0
#         Page size (bytes): 4096
#         Exit status: 0

# bam_in="${dir_data}/mm10/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"  # mm10
# bam_in="${dir_data}/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"  # CAST

# qnames="${dir_experiment_sub}/04_run_generate-assignment-lists/Disteche_sample_11.ambiguous.CAST-combined.rv-srt.txt.gz"  # keep only mm10
# qnames="${dir_experiment_sub}/04_run_generate-assignment-lists/Disteche_sample_11.ambiguous.mm10-combined.rv-srt.txt.gz"  # keep only CAST
# qnames="${dir_experiment_sub}/04_run_generate-assignment-lists/Disteche_sample_11.mm10-combined.CAST-combined.rv-srt.txt.gz"  # keep only ambiguous

# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.mm10.bam"  # CAST bam, mm10 assignments
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.CAST.bam"  # CAST bam, CAST assignments
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.ambiguous.bam"  # CAST bam, ambiguous assignments
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.mm10.mm10.bam"  # mm10 bam, mm10 assignments
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.mm10.CAST.bam"  # mm10 bam, CAST assignments
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.mm10.ambiguous.bam"  # mm10 bam, ambiguous assignments


#  Initial work...
# #  Python-specific installations
# pip install pysam
# pip install tqdm


# python ./bin/workflow/exclude-reads-from-bam.py
# usage: exclude-reads-from-bam.py [-h] [-i BAM_IN] [-q QNAMES] [-o BAM_OUT]
#
# Filter bam infile to exclude reads based on list of QNAMEs (supplied in txt file); write to bam outfile or STDOUT.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -i BAM_IN, --bam_in BAM_IN
#                         QNAME-sorted bam infile to be filtered, including path (default: None)
#   -q QNAMES, --qnames QNAMES
#                         text infile containing QNAMEs to be excluded, including path (default: None)
#   -o BAM_OUT, --bam_out BAM_OUT
#                         bam outfile, including path (default: None)
#
# dir_base="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross"
# dir_data="${dir_base}/data/files_bam_test"
# dir_results="${dir_base}/results"
# dir_experiment="${dir_results}/2022-0615_test_driver_allelic-segregation_Disteche-sample-11_fix-Picard-memory-issue"
# dir_experiment_sub="${dir_experiment}/results/Disteche_sample_11"
# dir_experiment_out="${dir_experiment_sub}/06_run_exclude-reads-from-bam"
#
# cd "${dir_experiment_sub}" || echo "cd failed. Check on this."
# mkdir -p "${dir_experiment_out}"
#
# bam_in="${dir_data}/CAST/Disteche_sample_11.dedup.corrected.corrected.queryname.bam"
# qnames="${dir_experiment_sub}/03_run_generate-assignment-lists/Disteche_sample_11.ambiguous.sort.txt.gz"
# bam_out="${dir_experiment_out}/Disteche_sample_11.dedup.corrected.corrected.CAST.ambiguous.bam"
#
# python exclude-reads-from-bam.py \
# --bam_in "${bam_in}" \
# --qnames "${qnames}" \
# --bam_out "${bam_out}"
#
# samtools view "${bam_in}" | head
# zcat "${qnames}" | head
#
#
# #  ambiguous + combined mm10
# cat \
# ./03_run_generate-assignment-lists/Disteche_sample_11.ambiguous.sort.txt.gz \
# ./04_run_filter-qnames-by-assignment/Disteche_sample_11.mm10.combined.sort.txt.gz \
# > ./04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.mm10-combined.txt.gz
#
# #  ambiguous + combined CAST
# cat \
# ./03_run_generate-assignment-lists/Disteche_sample_11.ambiguous.sort.txt.gz \
# ./04_run_filter-qnames-by-assignment/Disteche_sample_11.CAST.combined.sort.txt.gz \
# > ./04_run_filter-qnames-by-assignment/Disteche_sample_11.ambiguous.CAST-combined.txt.gz
#
# #  combined mm10 + combined CAST
# cat \
# ./04_run_filter-qnames-by-assignment/Disteche_sample_11.mm10.combined.sort.txt.gz \
# ./04_run_filter-qnames-by-assignment/Disteche_sample_11.CAST.combined.sort.txt.gz \
# > ./04_run_filter-qnames-by-assignment/Disteche_sample_11.mm10-combined.CAST-combined.txt.gz
#
# cd ./04_run_filter-qnames-by-assignment || echo "cd failed. Check on this."
#
# gz_files=(
#     Disteche_sample_11.ambiguous.mm10-combined.txt.gz
#     Disteche_sample_11.ambiguous.CAST-combined.txt.gz
#     Disteche_sample_11.mm10-combined.CAST-combined.txt.gz
# )
# echo "${gz_files[0]}"
# echo "${gz_files[1]}"
# echo "${gz_files[2]}"
# for i in "${gz_files[@]}"; do time parsort_file_qname_auto "${i}"; done
# # real    27m43.650s
# # user    14m12.790s
# # sys     1m48.633s
# #
# # real    23m45.014s
# # user    13m15.080s
# # sys     1m40.992s
# #
# # real    18m56.667s
# # user    12m5.515s
# # sys     1m31.240s
#
# for i in "${gz_files[@]}"; do time parsort_file_qname_reverse_auto "${i}"; done
# # real    22m21.849s
# # user    14m11.641s
# # sys     1m49.080s
# #
# # real    21m21.585s
# # user    13m27.932s
# # sys     1m45.275s
# #
# # real    19m19.347s
# # user    12m7.599s
# # sys     1m33.214s

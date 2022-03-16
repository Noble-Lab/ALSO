#!/bin/bash

cd /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools
module add samtools/1.14

# repair is used to put mates next to each other; 
# it's part of the subread package (http://subread.sourceforge.net/); 
# if there's no module for subread, then you can do a local installation through sourceforge or conda

ln -s /net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/bin/subread-2.0.3-Linux-x86_64 subread
subread/bin/utilities/repair 

#  Call from repo home directory

#  Test for one chromosome, chrX
bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_chr19/" \
-c "chr19" \
-r "TRUE" \
-b "TRUE" \
-p 6

#  Test for all chromosomes
bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_all/" \
-c "all" \
-r "TRUE" \
-b "TRUE" \
-p 6

# -h <print this help message and exit>
# -u <use safe mode: TRUE or FALSE (logical)>
# -i <bam infile, including path (chr)>
# -o <path for split bam file(s) (chr)>
# -c <chromosome(s) to split out (chr); for example, "chr1" for
#     chromosome 1, "chrX" for chromosome X, "all" for all
#     chromosomes>
# -r <use Subread repair on split bam files: TRUE or FALSE (logical)>
# -b <create bed files from split bam files: TRUE or FALSE (logical)>
# -p <number of cores for parallelization (int >= 1)>

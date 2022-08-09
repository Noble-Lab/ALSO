#!/bin/bash

#  downsample-fastqs.sh
#  KA

#  The package 'BBMap' needs to be in your "${PATH}"
#+ Get it here: https://sourceforge.net/projects/bbmap/

f1="${1:-"/net/noble/vol1/home/gangliuw/proj/2022-07-HiChIP-QC/data/K27ac_Hi_ChIP_patWT_S1/K27ac_Hi_ChIP_patWT_S1_R1.fastq.gz"}"
f2="${2:-"/net/noble/vol1/home/gangliuw/proj/2022-07-HiChIP-QC/data/K27ac_Hi_ChIP_patWT_S1/K27ac_Hi_ChIP_patWT_S1_R2.fastq.gz"}"
samp="${3:-50k}"
dir_out="${4:-"/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data/files_fastq"}"

reformat.sh \
in1="${f1}" \
in2="${f2}" \
out1="${dir_out}/$(basename "${f1%.fastq.gz}.${samp}.fastq.gz")" \
out2="${dir_out}/$(basename "${f2%.fastq.gz}.${samp}.fastq.gz")" \
samplereadstarget="${samp}"

#!/bin/bash

#  get-files-from-HPC.sh
#  KA

sGrab() { scp -p nexus:"$@" .; }

dir_HPC_4DN="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross"
dir_HPC_CAST="${dir_HPC_4DN}/results/2021-1027_Disteche_CAST-EiJ_test/get_unique_fragments"
dir_HPC_mm10="${dir_HPC_4DN}/results/2021-1028_Disteche_mm10_test/get_unique_fragments"
dir_data_local="./data/files_bam_test"

cd "${dir_data_local}" || ! echo "Directory not found; check on this..."

sGrab "${dir_HPC_CAST}/Disteche_sample_1.dedup.bam"
mv "Disteche_sample_1.dedup.bam" "Disteche_sample_1.CAST-EiJ.dedup.bam"

sGrab "${dir_HPC_mm10}/Disteche_sample_1.dedup.bam"
mv "Disteche_sample_1.dedup.bam" "Disteche_sample_1.mm10.dedup.bam"

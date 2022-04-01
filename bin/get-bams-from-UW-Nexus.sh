#!/bin/bash

#  get-bams-from-UW-Nexus.sh
#  KA


#  Function -------------------------------------------------------------------
# Download a single file from a directory on the UW HPC server
# :param $1: a single file to grab, including path
sGrab() { scp -p nexus:"${1}" .; }


#  Download files -------------------------------------------------------------
dir_HPC_kga0="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results"

dir_HPC_CAST="${dir_HPC_kga0}/2021-1027_Disteche_CAST-EiJ_test/get_unique_fragments"
sample_1_CAST="Disteche_sample_1.dedup.bam"

dir_HPC_mm10="${dir_HPC_kga0}/2021-1028_Disteche_mm10_test/get_unique_fragments"
sample_1_mm10="Disteche_sample_1.dedup.bam"

dir_data_local="${1:-"./data/files_bam_test"}"

cd "${dir_data_local}" || ! echo "Directory not found; check on this..."

sGrab "${dir_HPC_CAST}/${sample_1_CAST}"
mv "${sample_1_CAST}" "${sample_1_CAST/.dedup.bam/.CAST-EiJ.dedup.bam}"

sGrab "${dir_HPC_mm10}/${sample_1_mm10}"
mv "${sample_1_mm10}" "${sample_1_CAST/.dedup.bam/.mm10.dedup.bam}"


# dir_ganliuw="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-08-samtools"
#
# dir_HPC_mm10="${dir_ganliuw}/mm10-output/sorted"
# sample_22_mm10="Disteche_sample_22.sorted.bam"
# sample_7_mm10="Disteche_sample_7.sorted.bam"
#
# dir_HPC_CAST="${dir_ganliuw}/CAST-EiJ-output/sorted"
# sample_22_CAST="Disteche_sample_22.sorted.bam"
# sample_7_CAST="Disteche_sample_7.sorted.bam"
#
# dir_data_local="./data/files_bam_test"
#
# cd "${dir_data_local}" || ! echo "Directory not found; check on this..."
#
# sGrab "${dir_HPC_mm10}/${sample_22_mm10}"
# mv "${sample_22_mm10}" "${sample_22_mm10/.sorted.bam/.mm10.sorted.bam}"
#
# sGrab "${dir_HPC_mm10}/${sample_7_mm10}"
# mv "${sample_7_mm10}" "${sample_7_mm10/.sorted.bam/.mm10.sorted.bam}"
#
# sGrab "${dir_HPC_CAST}/${sample_22_CAST}"
# mv "${sample_22_CAST}" "${sample_22_CAST/.sorted.bam/.mm10.sorted.bam}"
#
# sGrab "${dir_HPC_CAST}/${sample_7_CAST}"
# mv "${sample_7_CAST}" "${sample_7_CAST/.sorted.bam/.mm10.sorted.bam}"

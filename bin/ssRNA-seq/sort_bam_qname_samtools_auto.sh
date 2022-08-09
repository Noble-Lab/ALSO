#!/bin/bash

#  sort_bam_qname_samtools_auto.sh
#  KA


#  Load module ----------------------------------------------------------------
# shellcheck disable=1091
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble samtools/1.14


#  Source functions into environment ------------------------------------------
# shellcheck disable=1091
. ./functions_single-end-reads.sh ||
    {
        echo "Exiting: Unable to source 'functions_single-end-reads.sh'."
        exit 1
    }


#  Parse positional arguments and set up variables ----------------------------
parallelize="${1:-"4"}"
infile="${2}"


#  Run function to sort bam files by QNAME ------------------------------------
sort_bam_qname_samtools_auto "${parallelize}" "${infile}"

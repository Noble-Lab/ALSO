#!/bin/bash

#  print_QNAME-lines_gt1.sh
#  KA


#  Load module ----------------------------------------------------------------
# shellcheck disable=1091
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble gzip/1.12


#  Source functions into environment ------------------------------------------


#  Parse positional arguments and set up variables ----------------------------
threshold="${1}"
infile="${2}"


#  Run function to sort bam files by QNAME ------------------------------------
awk -v threshold="${threshold}" '$1 > threshold' <(zcat -dk "${infile}") \
    | gzip \
    > "${infile%.txt}.gt1.txt.gz"

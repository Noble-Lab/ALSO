#!/bin/bash

#  submit_driver-ALSO_brain-rep-2.sh
#  KA

strain_1="mm11"
strain_2="SPRET"

p_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
p_bam="data/files_bam"

f_bam_1="brain_rep_2.SRR1525407.mm11.Aligned.sortedByCoord.out.primary.rm.bam"
f_bam_2="brain_rep_2.SRR1525407.SPRET-EiJ.Aligned.sortedByCoord.out.primary.rm.bam"

bam_1="${p_base}/${p_bam}/${strain_1}/${f_bam_1}"
bam_2="${p_base}/${p_bam}/${strain_2}/${f_bam_2}"

prefix="brain_rep_2"

p_out="results/kga0/2022-0730_ALSO_${prefix}"

end="single"

bash ./bin/workflow/driver_ALSO.sh \
-u FALSE \
-d FALSE \
-r "${strain_1}" \
-s "${strain_2}" \
-1 "${bam_1}" \
-2 "${bam_2}" \
-p "${prefix}" \
-o "${p_base}/${p_out}" \
-e "${end}"

# -h  print this help message and exit
# -u  use safe mode: TRUE or FALSE <logical> [default: FALSE]
# -d  run pipeline in ${TMPDIR}: TRUE or FALSE <logical> [default:
#     TRUE]
# -m  initial memory allocation pool (JVM) <chr> [default: "512m"]
# -x  maximum memory allocation pool (JVM) <chr> [default: "4096m"]
# -r  string for "sample #1" <chr>
# -s  string for "sample #2" <chr>
# -1  bam infile #1, including path <chr>
# -2  bam infile #2, including path <chr>
# -p  prefix for outfiles <chr>
# -o  results directory for outfiles <chr>; path will be made if it
#     does not exist
# -e  reads from paired- or single-end sequencing run: "paired" or
#     "single" <chr> [default: "paired"]
# -b  number of records to read into memory at one time when running
#     the script for Part #1, get-AS-per-qname.R <int > 0> [default:
#     100000]
# -c  number of records to read into memory at one time when running
#     the script for Part #3, generate-assignment-lists.R <int > 0>
#     [default: 1000000]
# -t  alignment score threshold <int >= 0> [default: 0]; the absolute
#     value of the difference in alignment scores between "sample
#     #1" and "sample #2" must be greater than this value in order
#     for a sample-specific assignment to be made; if not greater than
#     this value, then the assignment will be "ambiguous"
# -a  count lines: TRUE or FALSE <logical> [default: TRUE]
# -n  step in pipeline to run up to <int 1-4> [default: 4]

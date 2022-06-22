#!/bin/bash

cd "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross" ||
    {
        echo "cd failed. Check on this."
    }

infile="./data/files_bam/Disteche_sample_13.dedup.CAST.bam"
outpath="./data"
bash ./bin/workflow/03-filter-qnames.sh \
-u FALSE \
-c FALSE \
-l FALSE \
-i "${infile}" \
-o "${outpath}" \
-f TRUE \
-r TRUE \
-p 4

infile="./data/files_bam/Disteche_sample_13.dedup.mm10.bam"
outpath="./data"
bash ./bin/workflow/03-filter-qnames.sh \
-u FALSE \
-c FALSE \
-l FALSE \
-i "${infile}" \
-o "${outpath}" \
-f TRUE \
-r TRUE \
-p 4
# -u use safe mode: "TRUE" or "FALSE" (logical)
# -c use KA conda environment: "TRUE" or "FALSE" (logical)
# -l run on GS HPC: "TRUE" or "FALSE" (logical)
# -i bam infile, including path (chr)
# -o path for outfiles (chr); path will be made if it does not exist
# -f run samtools flagstat on bams: "TRUE" or "FALSE" (logical)
# -r remove intermediate files: "TRUE" or "FALSE" (logical)
# -p number of cores for parallelization (int >= 1); default: 1


infile="./data/files_bam/Disteche_sample_13.dedup.CAST.corrected.bam"
outpath="./data"
Rscript ./bin/workflow/06-get-AS-per-qname.R \
-b "${infile}" \
-i "${infile}.bai" \
-o "${outpath}" \
-s "CAST" \
-c 100000

infile="./data/files_bam/Disteche_sample_13.dedup.mm10.corrected.bam"
outpath="./data"
Rscript ./bin/workflow/06-get-AS-per-qname.R \
-b "${infile}" \
-i "${infile}.bai" \
-o "${outpath}" \
-s "mm10" \
-c 100000
# -b, --bam     bam infile, including path <chr>
# -i, --bai     bam index, including path <chr>
# -o, --outdir  directory for saving rds outfile, including path <chr>
# -s, --strain  strain name to be appended to rds outfile columns <chr>
# -c, --chunk   number of records to read into memory at a single time
#               <even int> [default: 100000]

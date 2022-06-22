#!/bin/bash

#  test.get_unique_fragments_CX.sh
#  KA

# file="Disteche_sample_1.q10.sorted.merged.bam"
# chromosome="chr19"
# # samtools index "${file}"
# # samtools view -b "${file}" "${chromosome}" > "${file%.bam}.${chromosome}.bam"
# # samtools collate -@ "${parallelize}" "${file%.bam}.${chromosome}.bam" -o "${file%.bam}.${chromosome}.collate.bam"
# # samtools rmdup "${file%.bam}.${chromosome}.bam" "${file%.bam}.${chromosome}.rmdup.bam"

# infile="${file%.bam}.${chromosome}.bam"
# outfile="${file%.bam}.${chromosome}.CX.bam"
# fragments="${file%.bam}.${chromosome}.CX.fragments.tsv.gz"
# transposition_sites="${file%.bam}.${chromosome}.CX.transposition_sites.bed"
# duplicates="${file%.bam}.${chromosome}.CX.duplicates.txt"
# insert_sizes="${file%.bam}.${chromosome}.CX.insert_sizes.txt"

# python get_unique_fragments_CX.py \
# --output_bam "${outfile}" \
# --fragments "${fragments}" \
# --transposition_sites_bed "${transposition_sites}" \
# --duplicate_read_counts "${duplicates}" \
# --insert_sizes "${insert_sizes}" \
# "${infile}"


start="$(date +%s)"

file="Disteche_sample_1.q10.sorted.merged.bam"
infile="${file%.bam}.bam"
outfile="${file%.bam}.CX.bam"
fragments="${file%.bam}.CX.fragments.tsv.gz"
transposition_sites="${file%.bam}.CX.transposition_sites.bed"
duplicates="${file%.bam}.CX.duplicates.txt"
insert_sizes="${file%.bam}.CX.insert_sizes.txt"

python get_unique_fragments_CX.py \
--output_bam "${outfile}" \
--fragments "${fragments}" \
--transposition_sites_bed "${transposition_sites}" \
--duplicate_read_counts "${duplicates}" \
--insert_sizes "${insert_sizes}" \
"${infile}"

end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""

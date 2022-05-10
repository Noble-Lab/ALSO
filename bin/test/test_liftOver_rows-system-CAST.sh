# -------
#  Set up paths and files

# --
path_in="/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora"
path_out="/Users/kalavattam/Downloads/to-do/get_unique_fragments/Bonora"

path_liftOver="/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver"

# --
file_prefix="tmp.CAST"
file_in="${file_prefix}.bed"
file_in_rename="${file_prefix}.rename.bed"

file_out_lifted="${file_prefix}.lifted.bed"
file_out_lifted_rename="${file_prefix}.lifted.rename.bed"

file_out_unlifted="${file_prefix}.unlifted.bed"
file_out_unlifted_rename="${file_prefix}.unlifted.rename.bed"

file_liftOver="CAST-EiJ-to-mm10.over.chain.gz"

# -------
#  Convert infile chromosomes names from common to official
cat "${file_in}" \
| sed "s/chr1/CM003994.1/g" \
| sed "s/chr2/CM003995.1/g" \
| sed "s/chr3/CM003996.1/g" \
| sed "s/chr4/CM003997.1/g" \
| sed "s/chr5/CM003998.1/g" \
| sed "s/chr6/CM003999.1/g" \
| sed "s/chr7/CM004000.1/g" \
| sed "s/chr8/CM004001.1/g" \
| sed "s/chr9/CM004002.1/g" \
| sed "s/chr10/CM004003.1/g" \
| sed "s/chr11/CM004004.1/g" \
| sed "s/chr12/CM004005.1/g" \
| sed "s/chr13/CM004006.1/g" \
| sed "s/chr14/CM004007.1/g" \
| sed "s/chr15/CM004008.1/g" \
| sed "s/chr16/CM004009.1/g" \
| sed "s/chr17/CM004010.1/g" \
| sed "s/chr18/CM004011.1/g" \
| sed "s/chr19/CM004012.1/g" \
| sed "s/chrX/CM004013.1/g" \
> "${file_in_rename}"

# -------
#  Do the liftOver, saving those regions that are and are not lifted; check the
#+ run time for liftOver too
start=$(date +%s)

liftOver \
-bedPlus=3 \
-tab \
"${path_in}/${file_in_rename}" \
"${path_liftOver}/${file_liftOver}" \
"${path_out}/${file_out_lifted}" \
"${path_out}/${file_out_unlifted}"

end=$(date +%s)

run_time=$(echo "${end} - ${start}" | bc -l)
echo "Run time is ${run_time} seconds."

# -------
#  Convert the outfile chromosome names from official to common
cat "${file_out_lifted}" \
| sed "s/CM003994.1/chr1/g" \
| sed "s/CM003995.1/chr2/g" \
| sed "s/CM003996.1/chr3/g" \
| sed "s/CM003997.1/chr4/g" \
| sed "s/CM003998.1/chr5/g" \
| sed "s/CM003999.1/chr6/g" \
| sed "s/CM004000.1/chr7/g" \
| sed "s/CM004001.1/chr8/g" \
| sed "s/CM004002.1/chr9/g" \
| sed "s/CM004003.1/chr10/g" \
| sed "s/CM004004.1/chr11/g" \
| sed "s/CM004005.1/chr12/g" \
| sed "s/CM004006.1/chr13/g" \
| sed "s/CM004007.1/chr14/g" \
| sed "s/CM004008.1/chr15/g" \
| sed "s/CM004009.1/chr16/g" \
| sed "s/CM004010.1/chr17/g" \
| sed "s/CM004011.1/chr18/g" \
| sed "s/CM004012.1/chr19/g" \
| sed "s/CM004013.1/chrX/g" \
> "${file_out_lifted_rename}"

cat "${file_out_unlifted}" \
| sed "s/CM003994.1/chr1/g" \
| sed "s/CM003995.1/chr2/g" \
| sed "s/CM003996.1/chr3/g" \
| sed "s/CM003997.1/chr4/g" \
| sed "s/CM003998.1/chr5/g" \
| sed "s/CM003999.1/chr6/g" \
| sed "s/CM004000.1/chr7/g" \
| sed "s/CM004001.1/chr8/g" \
| sed "s/CM004002.1/chr9/g" \
| sed "s/CM004003.1/chr10/g" \
| sed "s/CM004004.1/chr11/g" \
| sed "s/CM004005.1/chr12/g" \
| sed "s/CM004006.1/chr13/g" \
| sed "s/CM004007.1/chr14/g" \
| sed "s/CM004008.1/chr15/g" \
| sed "s/CM004009.1/chr16/g" \
| sed "s/CM004010.1/chr17/g" \
| sed "s/CM004011.1/chr18/g" \
| sed "s/CM004012.1/chr19/g" \
| sed "s/CM004013.1/chrX/g" \
> "${file_out_unlifted_rename}"

# -------
#  Clean up
# rm "${file_in}"
rm "${file_in_rename}"
rm "${file_out_lifted}"
rm "${file_out_unlifted}"

mv "${file_out_lifted_rename}" "${file_out_lifted}"
mv "${file_out_unlifted_rename}" "${file_out_unlifted}"

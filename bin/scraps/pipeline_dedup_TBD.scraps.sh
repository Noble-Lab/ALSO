unset bam_chrX_SamSort
typeset -a bam_chrX
for i in "${bam_chrX[@]}"; do bam_chrX_SamSort+=( "${i/.chrX.bam/.chrX.sort_queryname.bam}" ); done
for i in "${bam_chrX_SamSort[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"picard SortSam INPUT={bam_initial} OUTPUT={bam_sorted_picard} SORT_ORDER=queryname" \
::: bam_initial "${bam_chrX[@]}" \
:::+ bam_sorted_picard "${bam_chrX_SamSort[@]}"

unset bam_chrX_SamSort_rmdup
for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_rmdup+=( "${i/.bam/.rmdup.bam}" ); done
for i in "${bam_chrX_SamSort_rmdup[@]}"; do echo "${i}"; done
echo ""

unset bam_chrX_SamSort_metric
for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_metric+=( "${i/.bam/.rmdup.txt}" ); done
for i in "${bam_chrX_SamSort_metric[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 2 \
"picard MarkDuplicates INPUT={bam_initial} OUTPUT={bam_rmdup} METRICS_FILE={bam_metric} REMOVE_DUPLICATES=true" \
::: bam_initial "${bam_chrX_SamSort[@]}" \
:::+ bam_rmdup "${bam_chrX_SamSort_rmdup[@]}" \
:::+ bam_metric "${bam_chrX_SamSort_metric[@]}"


unset bam_chrX_SamSort_queryname_rmdup_SamSort_coordinate
typeset -a bam_chrX
for i in "${bam_chrX_SamSort_rmdup[@]}"; do bam_chrX_SamSort_queryname_rmdup_SamSort_coordinate+=( "${i/.chrX.sort_queryname.rmdup.bam/.chrX.sort_queryname.rmdup.sort_coordinate.bam}" ); done
for i in "${bam_chrX_SamSort_queryname_rmdup_SamSort_coordinate[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"picard SortSam INPUT={bam_initial} OUTPUT={bam_sorted_picard} SORT_ORDER=coordinate" \
::: bam_initial "${bam_chrX_SamSort_rmdup[@]}" \
:::+ bam_sorted_picard "${bam_chrX_SamSort_queryname_rmdup_SamSort_coordinate[@]}"

# unset bam_chrX_SamSort_rmdup
# for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_rmdup+=( "${i/.bam/.rmdup.bam}" ); done
# for i in "${bam_chrX_SamSort_rmdup[@]}"; do echo "${i}"; done
# echo ""

# unset bam_chrX_SamSort_metric
# for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_metric+=( "${i/.bam/.rmdup.txt}" ); done
# for i in "${bam_chrX_SamSort_metric[@]}"; do echo "${i}"; done
# echo ""

# parallel --header : --colsep " " -k -j 2 \
# "picard MarkDuplicates INPUT={bam_initial} OUTPUT={bam_rmdup} METRICS_FILE={bam_metric} REMOVE_DUPLICATES=true" \
# ::: bam_initial "${bam_chrX_SamSort[@]}" \
# :::+ bam_rmdup "${bam_chrX_SamSort_rmdup[@]}" \
# :::+ bam_metric "${bam_chrX_SamSort_metric[@]}"

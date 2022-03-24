#!/bin/bash


# -------------------------------------
# https://www.biostars.org/p/415831/

# for bam in $(find . -name "*.bam" -printf "%P\n"); do
#     echo "Command will be..."
#     echo -e "cat "${bam}"\n \
#     | samtools sort -n -O BAM\n \
#     | samtools fixmate -m - -\n \
#     | samtools sort -O BAM\n \
#     | tee \"${bam/.bam/.avec_dups.bam}\"\n \
#     | samtools markdup -r - \"${bam/.bam/.sans_dups.bam}\""
#     echo ""
#
#     echo "Running command..."
#     cat "${bam}" \
#     | samtools sort -n -O BAM \
#     | samtools fixmate -m - - \
#     | samtools sort -O BAM \
#     | tee "${bam/.bam/.avec_dups.bam}" \
#     | samtools markdup -r - "${bam/.bam/.sans_dups.bam}"
#     echo "Command completed"
# done


# unset array_bam
# typeset -a array_bam=(
#     "CAST-EiJ.Disteche_sample_1.q10.sorted.merged.bam"
#     "CAST-EiJ.Disteche_sample_1.dedup.bam"
#     "129S1-SvImJ.run.0.lane1.F121-6-CASTx129.undifferentiated.bam"
#     "129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.bam"
#     "CAST-EiJ.run.0.lane1.F121-6-CASTx129.undifferentiated.bam"
#     "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.bam"
#     "mm10-CAST-129-Nmasked.run.0.lane1.F121-6-CASTx129.undifferentiated.bam"
#     "mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.bam"
#     "mm10.run.0.lane1.F121-6-CASTx129.undifferentiated.bam"
#     "mm10.F121-6-CASTx129.undifferentiated.dedup.bam"
# )

# unset array_bam
# typeset -a array_bam
# while IFS=" " read -r -d $'\n'; do
#     array_bam+=( "${REPLY}" )
# done < <(find . -name "*.bam" -printf "%P\n" | sort)
# for i in "${array_bam[@]}"; do echo "${i}"; done
# echo ""

# unset array_bam
# typeset -a array_bam=(
#     mm10.Disteche_sample_1.dedup.bam
#     mm10.Disteche_sample_1.q10.sorted.merged.bam
# )
# typeset -a array_bam=(
#     mm10-CAST-Nmasked.Disteche_sample_1.dedup.bam
#     mm10-CAST-Nmasked.Disteche_sample_1.q10.sorted.merged.bam
# )
# for i in "${array_bam[@]}"; do echo "${i}"; done
# echo ""

unset array_bam
typeset -a array_bam=(
    "129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.bam"
    "CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.bam"
    "mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.bam"
    "mm10.F121-6-CASTx129.undifferentiated.dedup.bam"
)
# typeset -a array_bam=(
#     "run.0.lane1.F121-6-CASTx129.undifferentiated.bam"
# )
# typeset -a array_bam=(
#     "sciATAC_mouseDiff.split.q10.nsorted.alt.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.ambig.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.contra.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.ref.bam"
# )
# typeset -a array_bam=(
#     "sciATAC_mouseDiff.split.q10.nsorted.alt.sort.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.ambig.sort.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.contra.sort.bam"
#     "sciATAC_mouseDiff.split.q10.nsorted.ref.sort.bam"
# )
for i in "${array_bam[@]}"; do echo "${i}"; done
echo ""


# -------------------------------------
# unset bam_avec_dups
# typeset -a bam_avec_dups
# for i in "${array_bam[@]}"; do bam_avec_dups+=( "${i/.bam/.avec_dups.bam}" ); done
# for i in "${bam_avec_dups[@]}"; do echo "${i}"; done 
# echo ""

# unset bam_sans_dups
# typeset -a bam_sans_dups
# for i in "${array_bam[@]}"; do bam_sans_dups+=( "${i/.bam/.sans_dups.bam}" ); done
# for i in "${bam_sans_dups[@]}"; do echo "${i}"; done 
# echo ""


# -------------------------------------
unset bam_MarkDuplicates
typeset -a bam_MarkDuplicates
for i in "${array_bam[@]}"; do bam_MarkDuplicates+=( "${i/.bam/.MarkDuplicates.bam}" ); done
for i in "${bam_MarkDuplicates[@]}"; do echo "${i}"; done
echo ""

unset bam_MarkDuplicates_file
typeset -a bam_MarkDuplicates_file
for i in "${bam_MarkDuplicates[@]}"; do bam_MarkDuplicates_file+=( "${i/.bam/.txt}" ); done
for i in "${bam_MarkDuplicates_file[@]}"; do echo "${i}"; done
echo ""

# #  Need to use '' (single quotes) for piping within the GNU parallel call
# parallel --header : --colsep " " -k -j 1 echo \
# 'cat {bam_initial} | samtools sort -n -O BAM | samtools fixmate -m - - | samtools sort -O BAM | tee {bam_avec_dups} | samtools markdup -r - {bam_sans_dups}' \
# ::: bam_initial "${array_bam[@]}" \
# :::+ bam_avec_dups "${bam_avec_dups[@]}" \
# :::+ bam_sans_dups "${bam_sans_dups[@]}"

# parallel --header : --colsep " " -k -j 1 \
# 'cat {bam_initial} | samtools sort -n -O BAM | samtools fixmate -m - - | samtools sort -O BAM | tee {bam_avec_dups} | samtools markdup -r - {bam_sans_dups}' \
# ::: bam_initial "${array_bam[@]}" \
# :::+ bam_avec_dups "${bam_avec_dups[@]}" \
# :::+ bam_sans_dups "${bam_sans_dups[@]}"

# parallel --header : --colsep " " -k -j 1 echo -e \
# "picard MarkDuplicates --INPUT {bam_initial} --OUTPUT {bam_MarkDuplicates} --METRICS_FILE {txt_MarkDuplicates}
# echo \"\"" \
# ::: bam_initial "${array_bam[@]}" \
# :::+ bam_MarkDuplicates "${bam_MarkDuplicates[@]}" \
# :::+ txt_MarkDuplicates "${bam_MarkDuplicates_file[@]}"

parallel --header : --colsep " " -k -j 2 \
"picard MarkDuplicates --INPUT {bam_initial} --OUTPUT {bam_MarkDuplicates} --METRICS_FILE {txt_MarkDuplicates}" \
::: bam_initial "${array_bam[@]}" \
:::+ bam_MarkDuplicates "${bam_MarkDuplicates[@]}" \
:::+ txt_MarkDuplicates "${bam_MarkDuplicates_file[@]}"


# -------------------------------------
# unset bam_MarkDuplicates
# typeset -a bam_MarkDuplicates=(
#     mm10.Disteche_sample_1.dedup.MarkDuplicates.bam
#     mm10.Disteche_sample_1.q10.sorted.merged.MarkDuplicates.bam
#     CAST-EiJ.Disteche_sample_1.dedup.MarkDuplicates.bam
#     CAST-EiJ.Disteche_sample_1.q10.sorted.merged.MarkDuplicates.bam
# )
# typeset -a bam_MarkDuplicates=(
#     mm10-CAST-Nmasked.Disteche_sample_1.dedup.MarkDuplicates.bam
#     mm10-CAST-Nmasked.Disteche_sample_1.q10.sorted.merged.MarkDuplicates.bam
# )
# for i in "${bam_MarkDuplicates[@]}"; do echo "${i}"; done
# echo ""

# unset bam_sorted
# typeset -a bam_sorted=(
#     mm10.Disteche_sample_1.dedup.MarkDuplicates.sort.bam
#     mm10.Disteche_sample_1.q10.merged.MarkDuplicates.sort.bam
#     CAST-EiJ.Disteche_sample_1.dedup.MarkDuplicates.sort.bam
#     CAST-EiJ.Disteche_sample_1.q10.merged.MarkDuplicates.sort.bam
# )
# typeset -a bam_sorted=(
#     mm10-CAST-Nmasked.Disteche_sample_1.dedup.MarkDuplicates.sort.bam
#     mm10-CAST-Nmasked.Disteche_sample_1.q10.merged.MarkDuplicates.sort.bam
# )
# for i in "${bam_sorted[@]}"; do echo "${i}"; done
# echo ""

unset bam_sorted
typeset -a bam_sorted
for i in "${bam_MarkDuplicates[@]}"; do bam_sorted+=( "${i/.bam/.sort.bam}" ); done
# for i in "${array_bam[@]}"; do bam_sorted+=( "${i/.bam/.sort.bam}" ); done
for i in "${bam_sorted[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 1 \
"samtools sort -@ 6 {bam_initial} > {bam_sorted}" \
::: bam_initial "${bam_MarkDuplicates[@]}" \
:::+ bam_sorted "${bam_sorted[@]}"
# parallel --header : --colsep " " -k -j 1 \
# "samtools sort -@ 6 {bam_initial} > {bam_sorted}" \
# ::: bam_initial "${array_bam[@]}" \
# :::+ bam_sorted "${bam_sorted[@]}"

# parallel --header : --colsep " " -k -j 4 \
# "samtools index {bam_sorted}" \
# ::: bam_sorted "${bam_MarkDuplicates[@]}"
parallel --header : --colsep " " -k -j 4 \
"samtools index {bam_sorted}" \
::: bam_sorted "${bam_sorted[@]}"
# parallel --header : --colsep " " -k -j 4 \
# "samtools index {bam_sorted}" \
# ::: bam_sorted "${bam_chr1[@]}"

# -----------------------------------------------------------------------------
unset bam_chr1
typeset -a bam_chr1
for i in "${bam_sorted[@]}"; do bam_chr1+=( "${i/.bam/.chr1.bam}" ); done
for i in "${bam_chr1[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr1 > {bam_chr1}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr1 "${bam_chr1[@]}"

parallel --header : --colsep " " -k -j 4 \
"samtools index {bam_sorted}" \
::: bam_sorted "${bam_chr1[@]}"


unset bam_chr2
typeset -a bam_chr2
for i in "${bam_sorted[@]}"; do bam_chr2+=( "${i/.bam/.chr2.bam}" ); done
for i in "${bam_chr2[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr2 > {bam_chr2}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr2 "${bam_chr2[@]}"


unset bam_chr3
typeset -a bam_chr3
for i in "${bam_sorted[@]}"; do bam_chr3+=( "${i/.bam/.chr3.bam}" ); done
for i in "${bam_chr3[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr3 > {bam_chr3}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr3 "${bam_chr3[@]}"


unset bam_chr4
typeset -a bam_chr4
for i in "${bam_sorted[@]}"; do bam_chr4+=( "${i/.bam/.chr4.bam}" ); done
for i in "${bam_chr4[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr4 > {bam_chr4}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr4 "${bam_chr4[@]}"


unset bam_chr5
typeset -a bam_chr5
for i in "${bam_sorted[@]}"; do bam_chr5+=( "${i/.bam/.chr5.bam}" ); done
for i in "${bam_chr5[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr5 > {bam_chr5}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr5 "${bam_chr5[@]}"


unset bam_chr6
typeset -a bam_chr6
for i in "${bam_sorted[@]}"; do bam_chr6+=( "${i/.bam/.chr6.bam}" ); done
for i in "${bam_chr6[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr6 > {bam_chr6}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr6 "${bam_chr6[@]}"


unset bam_chr7
typeset -a bam_chr7
for i in "${bam_sorted[@]}"; do bam_chr7+=( "${i/.bam/.chr7.bam}" ); done
for i in "${bam_chr7[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr7 > {bam_chr7}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr7 "${bam_chr7[@]}"


unset bam_chr8
typeset -a bam_chr8
for i in "${bam_sorted[@]}"; do bam_chr8+=( "${i/.bam/.chr8.bam}" ); done
for i in "${bam_chr8[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr8 > {bam_chr8}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr8 "${bam_chr8[@]}"


unset bam_chr9
typeset -a bam_chr9
for i in "${bam_sorted[@]}"; do bam_chr9+=( "${i/.bam/.chr9.bam}" ); done
for i in "${bam_chr9[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr9 > {bam_chr9}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr9 "${bam_chr9[@]}"


unset bam_chr10
typeset -a bam_chr10
for i in "${bam_sorted[@]}"; do bam_chr10+=( "${i/.bam/.chr10.bam}" ); done
for i in "${bam_chr10[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr10 > {bam_chr10}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr10 "${bam_chr10[@]}"


unset bam_chr11
typeset -a bam_chr11
for i in "${bam_sorted[@]}"; do bam_chr11+=( "${i/.bam/.chr11.bam}" ); done
for i in "${bam_chr11[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr11 > {bam_chr11}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr11 "${bam_chr11[@]}"


unset bam_chr12
typeset -a bam_chr12
for i in "${bam_sorted[@]}"; do bam_chr12+=( "${i/.bam/.chr12.bam}" ); done
for i in "${bam_chr12[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr12 > {bam_chr12}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr12 "${bam_chr12[@]}"


unset bam_chr13
typeset -a bam_chr13
for i in "${bam_sorted[@]}"; do bam_chr13+=( "${i/.bam/.chr13.bam}" ); done
for i in "${bam_chr13[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr13 > {bam_chr13}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr13 "${bam_chr13[@]}"


unset bam_chr14
typeset -a bam_chr14
for i in "${bam_sorted[@]}"; do bam_chr14+=( "${i/.bam/.chr14.bam}" ); done
for i in "${bam_chr14[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr14 > {bam_chr14}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr14 "${bam_chr14[@]}"


unset bam_chr15
typeset -a bam_chr15
for i in "${bam_sorted[@]}"; do bam_chr15+=( "${i/.bam/.chr15.bam}" ); done
for i in "${bam_chr15[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr15 > {bam_chr15}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr15 "${bam_chr15[@]}"


unset bam_chr16
typeset -a bam_chr16
for i in "${bam_sorted[@]}"; do bam_chr16+=( "${i/.bam/.chr16.bam}" ); done
for i in "${bam_chr16[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr16 > {bam_chr16}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr16 "${bam_chr16[@]}"


unset bam_chr17
typeset -a bam_chr17
for i in "${bam_sorted[@]}"; do bam_chr17+=( "${i/.bam/.chr17.bam}" ); done
for i in "${bam_chr17[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr17 > {bam_chr17}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr17 "${bam_chr17[@]}"


unset bam_chr18
typeset -a bam_chr18
for i in "${bam_sorted[@]}"; do bam_chr18+=( "${i/.bam/.chr18.bam}" ); done
for i in "${bam_chr18[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr18 > {bam_chr18}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr18 "${bam_chr18[@]}"


unset bam_chr19
typeset -a bam_chr19
for i in "${bam_sorted[@]}"; do bam_chr19+=( "${i/.bam/.chr19.bam}" ); done
for i in "${bam_chr19[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chr19 > {bam_chr19}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chr19 "${bam_chr19[@]}"


# -----------------------------------------------------------------------------
unset bam_chrX
typeset -a bam_chrX
for i in "${bam_sorted[@]}"; do bam_chrX+=( "${i/.bam/.chrX.bam}" ); done
for i in "${bam_chrX[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chrX > {bam_chrX}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chrX "${bam_chrX[@]}"

parallel --header : --colsep " " -k -j 4 \
"samtools index {bam_sorted}" \
::: bam_sorted "${bam_chrX[@]}"

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


unset bam_chrX_SamSort
typeset -a bam_chrX
for i in "${bam_chrX[@]}"; do bam_chrX_SamSort+=( "${i/.chrX.bam/.chrX.sort_coordinate.bam}" ); done
for i in "${bam_chrX_SamSort[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"picard SortSam INPUT={bam_initial} OUTPUT={bam_sorted_picard} SORT_ORDER=coordinate" \
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



unset bam_chrX_rmdup
for i in "${bam_chrX[@]}"; do bam_chrX_rmdup+=( "${i/.chrX.bam/.chrX.rmdup.bam}" ); done
for i in "${bam_chrX_rmdup[@]}"; do echo "${i}"; done
echo ""

unset bam_chrX_metric
for i in "${bam_chrX[@]}"; do bam_chrX_metric+=( "${i/.chrX.bam/.chrX.rmdup.txt}" ); done
for i in "${bam_chrX_metric[@]}"; do echo "${i}"; done
echo ""


unset bam_chrX_SamSort_markdup
for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_markdup+=( "${i/.bam/.markdup.bam}" ); done
for i in "${bam_chrX_SamSort_markdup[@]}"; do echo "${i}"; done
echo ""

unset bam_chrX_SamSort_markdup_metric
for i in "${bam_chrX_SamSort[@]}"; do bam_chrX_SamSort_markdup_metric+=( "${i/.bam/.markdup.txt}" ); done
for i in "${bam_chrX_SamSort_markdup_metric[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 2 \
"picard MarkDuplicates INPUT={bam_initial} OUTPUT={bam_markdup} METRICS_FILE={bam_metric}" \
::: bam_initial "${bam_chrX_SamSort[@]}" \
:::+ bam_markdup "${bam_chrX_SamSort_markdup[@]}" \
:::+ bam_metric "${bam_chrX_SamSort_markdup_metric[@]}"



parallel --header : --colsep " " -k -j 2 \
"picard MarkDuplicates INPUT={bam_initial} OUTPUT={bam_rmdup} METRICS_FILE={bam_metric} REMOVE_DUPLICATES=true" \
::: bam_initial "${bam_chrX[@]}" \
:::+ bam_rmdup "${bam_chrX_rmdup[@]}" \
:::+ bam_metric "${bam_chrX_metric[@]}"

parallel --header : --colsep " " -k -j 4 \
"samtools index {bam_sorted}" \
::: bam_sorted "${bam_chrX_rmdup[@]}"
# -----------------------------------------------------------------------------


unset bam_chrY
typeset -a bam_chrY
for i in "${bam_sorted[@]}"; do bam_chrY+=( "${i/.bam/.chrY.bam}" ); done
for i in "${bam_chrY[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chrY > {bam_chrY}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chrY "${bam_chrY[@]}"
# rm -- *MarkDuplicates.bam


unset bam_chrM
typeset -a bam_chrM
for i in "${bam_sorted[@]}"; do bam_chrM+=( "${i/.bam/.chrM.bam}" ); done
for i in "${bam_chrM[@]}"; do echo "${i}"; done
echo ""

parallel --header : --colsep " " -k -j 4 \
"samtools view -b {bam_sorted} chrM > {bam_chrM}" \
::: bam_sorted "${bam_sorted[@]}" \
:::+ bam_chrM "${bam_chrM[@]}"

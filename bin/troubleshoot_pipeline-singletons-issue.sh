#!/bin/bash

#  troubleshoot_pipeline-singletons-issue.sh
#  KA


#  Functions ------------------------------------------------------------------
# :param $1: a single file to grab, including path
sGrab() { scp -p nexus:"${1}" .; }


#  Check Disteche_sample_6.chr11.* --------------------------------------------
typeset path_files_to_test
path_files_to_test="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-26-split-lift-CAST/CAST-split"

unset files_to_test
typeset -a files_to_test=(
    Disteche_sample_6.chr11.bam
    Disteche_sample_6.chr11.bam.bai
    Disteche_sample_6.chr11.mpos.bed
    Disteche_sample_6.chr11.pos.bed
)

unset files_err_out
typeset -a files_err_out=(
    split_sample_6.e
    split_sample_6.o
)

typeset parallelize="${1:-"4"}"

4dn-data
mkdir -p 2022-0330_troubleshoot_pipeline-singletons-issue
cd 2022-0330_troubleshoot_pipeline-singletons-issue || ! echo "Warning: cd failed."

#  Get bam, bai, and bed files
for i in "${files_to_test[@]}"; do
    printf "Started: Grabbing %s\n" "${path_files_to_test}/${i}"
    sGrab "${path_files_to_test}/${i}"
    printf "Started: Grabbing %s\n\n" "${path_files_to_test}/${i}"
done

#  Get stderr and stdout files
for i in "${files_err_out[@]}"; do
    printf "Started: Grabbing %s\n" "${path_files_to_test}/${i}"
    sGrab "${path_files_to_test}/${i}"
    printf "Started: Grabbing %s\n\n" "${path_files_to_test}/${i}"
done

#  Get README documenting Gang's tests
sGrab "/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-03-26-split-lift-CAST/readme"
mv readme readme.md

#  Check flagstat
samtools flagstat -@ "${parallelize}" "${files_to_test[1]}" \
> "${files_to_test[1]/.bam/.samtools.flagstat.txt}"

less "${files_to_test[1]/.bam/.flagstat.txt}"
# 7240578 + 0 in total (QC-passed reads + QC-failed reads)
# 7240578 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 7240578 + 0 mapped (100.00% : N/A)
# 7240578 + 0 primary mapped (100.00% : N/A)
# 7240578 + 0 paired in sequencing
# 3620287 + 0 read1
# 3620291 + 0 read2
# 7240578 + 0 properly paired (100.00% : N/A)
# 7240578 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

#  Check stats
samtools stats -@ "${parallelize}" "${files_to_test[1]}" \
> "${files_to_test[1]/.bam/.samtools.stats.txt}"

less "${files_to_test[1]/.bam/.samtools.stats.txt}"
# ...

#  Check number of unpaired reads in split_sample_6.e
# shellcheck disable=SC2002
cat "${files_err_out[1]}" \
| grep "Unpaired reads" \
| cut -f3 -d ":" \
| paste -sd+ - \
| bc
# 201


#  Get the original Disteche_sample_6, then preprocess it ---------------------
typeset path_to_bam_input="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-23/get_unique_fragments"
typeset bam_input="Disteche_sample_6.dedup.bam"
typeset bam_processed="${bam_input/.dedup.bam/.processed.bam}"
typeset bam_processed_f1="${bam_input/.dedup.bam/.processed.f1.bam}"
typeset bam_processed_F8="${bam_input/.dedup.bam/.processed.F8.bam}"

sGrab "${path_to_bam_input}/${bam_input}"

    #  Check flagstat: initial
    samtools flagstat -@ "${parallelize}" "${bam_input}" \
    > "${bam_input/.bam/.samtools.flagstat.txt}"

    less "${bam_input/.bam/.samtools.flagstat.txt}"
    # 130725752 + 0 in total (QC-passed reads + QC-failed reads)
    # 130725752 + 0 primary
    # 0 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 0 + 0 primary duplicates
    # 130725752 + 0 mapped (100.00% : N/A)
    # 130725752 + 0 primary mapped (100.00% : N/A)
    # 130725752 + 0 paired in sequencing
    # 65362875 + 0 read1
    # 65362877 + 0 read2
    # 130725752 + 0 properly paired (100.00% : N/A)
    # 130725752 + 0 with itself and mate mapped
    # 0 + 0 singletons (0.00% : N/A)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)

    # #  Check stats: initial
    # samtools stats -@ "${parallelize}" "${bam_input}" \
    # > "${bam_input/.bam/.samtools.stats.txt}"
    #
    # less "${bam_input/.bam/.samtools.stats.txt}"
    # # ...

#  Preprocessing
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 \
"${bam_input}" \
-o "${bam_processed}"

    #  Check flagstat: preprocessed
    samtools flagstat -@ "${parallelize}" "${bam_processed}" \
    > "${bam_processed/.bam/.samtools.flagstat.txt}"

    less "${bam_processed/.bam/.samtools.flagstat.txt}"
    # 123308945 + 0 in total (QC-passed reads + QC-failed reads)
    # 123308945 + 0 primary
    # 0 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 0 + 0 primary duplicates
    # 123308945 + 0 mapped (100.00% : N/A)
    # 123308945 + 0 primary mapped (100.00% : N/A)
    # 123308945 + 0 paired in sequencing
    # 61654460 + 0 read1
    # 61654485 + 0 read2
    # 123308945 + 0 properly paired (100.00% : N/A)
    # 123308945 + 0 with itself and mate mapped
    # 0 + 0 singletons (0.00% : N/A)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)

    # #  Check stats: preprocessed
    # samtools stats -@ "${parallelize}" "${bam_processed}" \
    # > "${bam_processed/.bam/.samtools.stats.txt}"
    #
    # less "${bam_processed/.bam/.samtools.stats.txt}"
    # # ...

#  Preprocessing -f 1
samtools view -@ "${parallelize}" -h -b -f 1 \
"${bam_processed}" \
-o "${bam_processed_f1}"

    #  Check flagstat: preprocessed -f 1
    samtools flagstat -@ "${parallelize}" "${bam_processed_f1}" \
    > "${bam_processed_f1/.bam/.samtools.flagstat.txt}"

    less "${bam_processed_f1/.bam/.samtools.flagstat.txt}"
    # 123308945 + 0 in total (QC-passed reads + QC-failed reads)
    # 123308945 + 0 primary
    # 0 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 0 + 0 primary duplicates
    # 123308945 + 0 mapped (100.00% : N/A)
    # 123308945 + 0 primary mapped (100.00% : N/A)
    # 123308945 + 0 paired in sequencing
    # 61654460 + 0 read1
    # 61654485 + 0 read2
    # 123308945 + 0 properly paired (100.00% : N/A)
    # 123308945 + 0 with itself and mate mapped
    # 0 + 0 singletons (0.00% : N/A)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)

#  Preprocessing -F 8
samtools view -@ "${parallelize}" -h -b -F 8 \
"${bam_processed}" \
-o "${bam_processed_F8}"

    #  Check flagstat: preprocessed -f 1
    samtools flagstat -@ "${parallelize}" "${bam_processed_F8}" \
    > "${bam_processed_F8/.bam/.samtools.flagstat.txt}"

    less "${bam_processed_F8/.bam/.samtools.flagstat.txt}"
    # 123308945 + 0 in total (QC-passed reads + QC-failed reads)
    # 123308945 + 0 primary
    # 0 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 0 + 0 primary duplicates
    # 123308945 + 0 mapped (100.00% : N/A)
    # 123308945 + 0 primary mapped (100.00% : N/A)
    # 123308945 + 0 paired in sequencing
    # 61654460 + 0 read1
    # 61654485 + 0 read2
    # 123308945 + 0 properly paired (100.00% : N/A)
    # 123308945 + 0 with itself and mate mapped
    # 0 + 0 singletons (0.00% : N/A)
    # 0 + 0 with mate mapped to a different chr
    # 0 + 0 with mate mapped to a different chr (mapQ>=5)


#  Try running samtools fixmate to see if this updates the flagstat readout
#+ biostars.org/p/319966/#320238
samtools sort -n -@ "${parallelize}" \
"${bam_processed}" \
> "${bam_processed/.bam/.sort-name.bam}"

samtools fixmate -@ "${parallelize}" \
"${bam_processed/.bam/.sort-name.bam}" \
"${bam_processed/.bam/.sort-name.fixed.bam}"

    samtools flagstat -@ "${parallelize}" "${bam_processed/.bam/.sort-name.fixed.bam}" \
    > "${bam_processed/.bam/.sort-name.fixed.flagstat.txt}"

    less "${bam_processed/.bam/.sort-name.fixed.flagstat.txt}"
    # 123308945 + 0 in total (QC-passed reads + QC-failed reads)
    # 123308945 + 0 primary
    # 0 + 0 secondary
    # 0 + 0 supplementary
    # 0 + 0 duplicates
    # 0 + 0 primary duplicates
    # 123308945 + 0 mapped (100.00% : N/A)
    # 123308945 + 0 primary mapped (100.00% : N/A)
    # 123308744 + 0 paired in sequencing
    # 61654372 + 0 read1
    # 61654372 + 0 read2
    # 123300968 + 0 properly paired (99.99% : N/A)
    # 123308744 + 0 with itself and mate mapped
    # 0 + 0 singletons (0.00% : N/A)
    # 7332 + 0 with mate mapped to a different chr
    # 7332 + 0 with mate mapped to a different chr (mapQ>=5)

#  So far, preprocessing step should be...
#+ 
#+ 1. $ samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${bam_sciatac/.bam/.flagstat.txt}"
#+ 2. $ samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sciatac}" -o "${bam_process}"
#+ 3. $ samtools flagstat -@ "${parallelize}" "${bam_process}" > "${bam_process/.bam/.flagstat.txt}"
#+ 4. $ samtools sort -n -@ "${parallelize}" "${bam_process}" > "${bam_processed/.bam/.sort-name.bam}"
#+ 5. $ samtools fixmate -@ "${parallelize}" "${bam_process/.bam/.sort-name.bam}" "${bam_process/.bam/.sort-name.fixed.bam}"

bam_sciatac="Disteche_sample_6.dedup.bam"
bam_process="${bam_sciatac/.bam/.process.bam}"
bam_sort_name="${bam_process/.bam/.sort_name.bam}"
bam_fixmate="${bam_sort_name/.bam/.fixmate.bam}"
bam_sort_coord="${bam_fixmate/.bam/.sort_coord.bam}"
bam_inter="${bam_sort_coord/.bam/.inter.bam}"

samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "flagstat.1.sciatac.txt"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sciatac}" -o "${bam_process}"
samtools flagstat -@ "${parallelize}" "${bam_process}" > "flagstat.2.process.txt"
samtools sort -@ "${parallelize}" -n "${bam_process}" > "${bam_sort_name}"
samtools fixmate -@ "${parallelize}" "${bam_sort_name}" "${bam_fixmate}"
samtools flagstat -@ "${parallelize}" "${bam_fixmate}" > "flagstat.3.fixmate.txt"
samtools sort -@ "${parallelize}" "${bam_fixmate}" > "${bam_sort_coord}"
samtools flagstat -@ "${parallelize}" "${bam_sort_coord}" > "flagstat.4.sort_coord.txt"
samtools index -@ "${parallelize}" "${bam_sort_coord}"
samtools view -@ "${parallelize}" -h -b "${bam_sort_coord}" "chr11" -o "${bam_sort_coord/.bam/.chr11.bam}"
samtools flagstat -@ "${parallelize}" "${bam_sort_coord/.bam/.chr11.bam}" > "flagstat.5.chr11.txt"
samtools index -@ "${parallelize}" "${bam_sort_coord/.bam/.chr11.bam}"
samtools view -@ "${parallelize}" -H "${bam_sort_coord}" -o "${bam_sort_coord/.bam/.header.txt}"

rename 's/\-/\_/g' -- *


#  Extract interchromosomal mate pairs
#+ 
#+ Note that the first call to samtools view does not include -b in order to
#+ output a text stream
samtools view -h -F 14 "${bam_sort_coord}" \
| awk '$7 !~ /=/' \
| samtools view -b - \
> "${bam_inter}"

# shellcheck disable=SC2002
echo "$(( $(cat "${bam_inter}" | wc -l) / 2 )) mate pairs are interchromosomal."

# #  Remove interchromosomal mate pairs
# samtools view -@ "${parallelize}" -h "${bam_sort_coord}" \
# | awk '($3 != $7 && $7 != "=")' \
# | samtools view -b - \
# > "${bam_intra}"
#
# samtools view -@ "${parallelize}" -h "${bam_sort_coord}" \
# | awk '$7 != "=" || $1 ~ /^@/' \
# | samtools view -b - \
# "${bam_intra}"
#
# samtools flagstat -@ "${parallelize}" "${bam_sort_coord}" > "flagstat.5.intra.txt"

# samtools view -@ "${parallelize}" -f 2 -F 512 -b "${bam_sort_coord}" -o "${bam_sort_coord/.bam/.test.bam}"
#
# samtools flagstat -@ "${parallelize}" "${bam_sort_coord/.bam/.test.bam}" > "flagstat.5.test-f2-F512.txt"

#  Pipeline in progress
bam_sciatac="Disteche_sample_6.dedup.bam"
bam_preprocessed="Disteche_sample_6.dedup.preprocessed.bam"

samtools view -@ "${parallelize}" -h -b -f 3 -F 12 -q 30 "${bam_sciatac}" \
| samtools sort -n -@ "${parallelize}" - \
| samtools fixmate -@ "${parallelize}" - \
| samtools sort -@ "${parallelize}" - \
> "${bam_preprocessed}"












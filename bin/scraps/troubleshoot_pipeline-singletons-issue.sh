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

#NOTE 2/2 The above attempts to filter out interchromosomal reads are too slow
#NOTE 2/2 to be practical

#NOTE Below, attempting to follow advice fr/a biostars.org post...
#NOTE Does not work
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


#  Work: 2022-0401 ------------------------------------------------------------

#  Set up, clean up
4dn-data && \
cd 2022-0330_troubleshoot_pipeline-singletons-issue || \
! echo "Warning: cd failed."

rm -- *process*

#  Real work ------------------------------------------------------------------
#  Start recording time
start="$(date +%s)"

#  Function
displaySpinningIcon() {
    # Display "spinning icon" while a background process runs
    # 
    # :param 1: PID of the last program the shell ran in the background (int)
    # :param 2: message to be displayed next to the spinning icon (chr)
    spin="/|\\–"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


parallelize=6

bam_sciatac="Disteche_sample_6.dedup.bam"
bam_Q30="${bam_sciatac/.bam/.filter-q30.bam}"
bam_sort_n_1="${bam_Q30/.bam/.sort_name.bam}"
bam_fixmate_1="${bam_sort_n_1/.bam/.fixmate.bam}"
bam_sort_c_1="${bam_fixmate_1/.bam/.sort_coord.bam}"
bam_f3_F12="${bam_sort_c_1/.bam/.filter-f3-F12.bam}"
bam_sort_n_2="${bam_f3_F12/.bam/.sort_name.bam}"
bam_fixmate_2="${bam_sort_n_2/.bam/.fixmate.bam}"
bam_sort_c_2="${bam_fixmate_2/.bam/.sort_coord.bam}"

flagstat_0="flagstat.0.from-sciatac-pipeline.txt"
flagstat_1="flagstat.1.filter-Q30.txt"
flagstat_2="flagstat.2.fixmate.txt"
flagstat_3="flagstat.3.sort_coord.txt"
flagstat_4="flagstat.4.filter-f3-F12.txt"
flagstat_5="flagstat.5.fixmate.txt"
flagstat_6="flagstat.6.sort_coord.txt"
flagstat_7="flagstat.7.chr11.txt"

#  Step 0
echo -e "Started: Step 0: Running flagstat on data from sciatac_pipeline"
samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 0: Running flagstat on data from sciatac_pipeline\n"

#  Step 1
echo -e "Started: Step 1: Filtering for only mates with MAPQ >= 30"
samtools view -@ "${parallelize}" -h -b -q 30 "${bam_sciatac}" -o "${bam_Q30}" &
displaySpinningIcon $! "Filtering for only mates with MAPQ >= 30... "
samtools flagstat -@ "${parallelize}" "${bam_Q30}" > "${flagstat_1}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 1: Filtering for only mates with MAPQ >= 30\n"

#  Step 2
echo -e "Started: Step 2: Sorting by QNAME and performing 'fixmate' to update flags"
samtools sort -@ "${parallelize}" -n "${bam_Q30}" > "${bam_sort_n_1}" &
displaySpinningIcon $! "Sorting by QNAME... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_1}" "${bam_fixmate_1}" &
displaySpinningIcon $! "Performing 'fixmate' to update flags... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_2}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 2: Sorting by QNAME and performing 'fixmate' to update flags\n"

#  Free up some space
rm "${bam_Q30}"

#  Step 3
echo -e "Started: Step 3: Sorting by coordinate in preparation for further filtering"
samtools sort -@ "${parallelize}" "${bam_fixmate_1}" > "${bam_sort_c_1}" &
displaySpinningIcon $! "Sorting by coordinate... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_1}" > "${flagstat_3}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 3: Sorting by coordinate in preparation for further filtering\n"

#  Free up some space
rm "${bam_sort_n_1}"

#  Step 4
echo -e "Started: Step 4: Filtering for mates: '-f 3 -F 12'"
samtools view -@ "${parallelize}" -h -b -f 3 -F 12 "${bam_sort_c_1}" -o "${bam_f3_F12}" &
displaySpinningIcon $! "Filtering for mates: '-f 3 -F 12'... "
samtools flagstat -@ "${parallelize}" "${bam_f3_F12}" > "${flagstat_4}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 4: Filtering for mates: '-f 3 -F 12'\n"

#  Free up some space
rm "${bam_fixmate_1}"

#  Step 5
echo -e "Started: Step 5: Sorting by QNAME and performing 'fixmate' to update flags"
samtools sort -@ "${parallelize}" -n "${bam_f3_F12}" > "${bam_sort_n_2}" &
displaySpinningIcon $! "Sorting by QNAME... "
samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_2}" &
displaySpinningIcon $! "Performing 'fixmate' to update flags... "
samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_5}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 5: Sorting by QNAME and performing 'fixmate' to update flags\n"

#  Free up some space
rm "${bam_sort_c_1}"

#  Step 6
echo -e "Started: Step 6: Sorting by coordinate in preparation for splitting by chromosomes"
samtools sort -@ "${parallelize}" "${bam_fixmate_2}" > "${bam_sort_c_2}" &
displaySpinningIcon $! "Sorting by coordinate... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_2}" > "${flagstat_6}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 6: Sorting by coordinate in preparation for splitting by chromosomes\n"

#  Free up some space
rm "${bam_fixmate_2}"

#  Step 7
echo -e "Started: Step 7: Indexing bam, then splitting by chromosome(s)"
samtools index -@ "${parallelize}" "${bam_sort_c_2}" &
displaySpinningIcon $! "Indexing bam... "
samtools view -@ "${parallelize}" -h -b "${bam_sort_c_2}" "chr11" -o "${bam_sort_c_2/.bam/.chr11.bam}" &
displaySpinningIcon $! "Splitting bam by chromosome(s)... "
samtools flagstat -@ "${parallelize}" "${bam_sort_c_2/.bam/.chr11.bam}" > "${flagstat_7}" &
displaySpinningIcon $! "Running flagstat... "
echo -e "Completed: Step 7: Indexing bam, then splitting by chromosome(s)\n"

#  Step 8
echo -e "Started: Step 8: Indexing and creating header.txt for split chromosome bam(s)"
samtools index -@ "${parallelize}" "${bam_sort_c_2/.bam/.chr11.bam}" &
displaySpinningIcon $! "Indexing split chromosome bam(s)... "
samtools view -@ "${parallelize}" -H "${bam_sort_c_2}" -o "${bam_sort_c_2/.bam/.header.txt}" &
displaySpinningIcon $! "Creating header.txt for split chromosome bam(s)... "
echo -e "Completed: Step 8: Indexing and creating header.txt for split chromosome bam(s)\n"

#  End recording time
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""


#  Examining output -----------------------------------------------------------
# ❯ head -25 flagstat.0*  # Examine bam outfile from sciatac_pipeline
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

# ❯ head -25 flagstat.1*  # Examine outfile after filtering for MAPQ >= 30
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

# ❯ head -25 flagstat.2*  # ...after sorting by QNAME and performing 'fixmate'
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

# ❯ head -25 flagstat.3*  # ...after sorting by coordinate
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

# ❯ head -25 flagstat.4*  # ...after filtering for '-f 3 -F 12'
# 123300968 + 0 in total (QC-passed reads + QC-failed reads)
# 123300968 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 123300968 + 0 mapped (100.00% : N/A)
# 123300968 + 0 primary mapped (100.00% : N/A)
# 123300968 + 0 paired in sequencing
# 61650480 + 0 read1
# 61650488 + 0 read2
# 123300968 + 0 properly paired (100.00% : N/A)
# 123300968 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# ❯ head -25 flagstat.5*  # ...after sorting by QNAME and performing 'fixmate'
# 123300968 + 0 in total (QC-passed reads + QC-failed reads)
# 123300968 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 123300968 + 0 mapped (100.00% : N/A)
# 123300968 + 0 primary mapped (100.00% : N/A)
# 123300968 + 0 paired in sequencing
# 61650480 + 0 read1
# 61650488 + 0 read2
# 123300968 + 0 properly paired (100.00% : N/A)
# 123300968 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# ❯ head -25 flagstat.6*  # ...after sorting by coordinate
# 123300968 + 0 in total (QC-passed reads + QC-failed reads)
# 123300968 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 123300968 + 0 mapped (100.00% : N/A)
# 123300968 + 0 primary mapped (100.00% : N/A)
# 123300968 + 0 paired in sequencing
# 61650480 + 0 read1
# 61650488 + 0 read2
# 123300968 + 0 properly paired (100.00% : N/A)
# 123300968 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# ❯ head -25 flagstat.7*  # ...after splitting bam by chromosome(s)
# 7240138 + 0 in total (QC-passed reads + QC-failed reads)
# 7240138 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 7240138 + 0 mapped (100.00% : N/A)
# 7240138 + 0 primary mapped (100.00% : N/A)
# 7240138 + 0 paired in sequencing
# 3620074 + 0 read1
# 3620064 + 0 read2
# 7240138 + 0 properly paired (100.00% : N/A)
# 7240138 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


# Started: Step 0: Running flagstat on data from sciatac_pipeline
# [1] 50668
# | Running flagstat... [1]  + 50668 done       samtools flagstat -@ "${parallelize}" "${bam_sciatac}" > "${flagstat_0}"
# Completed: Step 0: Running flagstat on data from sciatac_pipeline
#
#
# Started: Step 1: Filtering for only mates with MAPQ >= 30
# [1] 50747
# \ Filtering for only mates with MAPQ >= 30... [1]  + 50747 done       samtools view -@ "${parallelize}" -h -b -q 30 "${bam_sciatac}" -o "${bam_Q30}
# [1] 51144
# / Running flagstat... [1]  + 51144 done       samtools flagstat -@ "${parallelize}" "${bam_Q30}" > "${flagstat_1}"
# Completed: Step 1: Filtering for only mates with MAPQ >= 30
#
#
# Started: Step 2: Sorting by QNAME and performing 'fixmate' to update flags
# [1] 51214
# | Sorting by QNAME... [bam_sort_core] merging from 36 files and 6 in-memory blocks...
# | Sorting by QNAME... [1]  + 51214 done       samtools sort -@ "${parallelize}" -n "${bam_Q30}" > "${bam_sort_n_1}"
# [1] 54370
# / Performing 'fixmate' to update flags... [1]  + 54370 done       samtools fixmate -@ "${parallelize}" "${bam_sort_n_1}" "${bam_fixmate_1}"
# [1] 54984
# | Running flagstat... [1]  + 54984 done       samtools flagstat -@ "${parallelize}" "${bam_fixmate_1}" > "${flagstat_2}"
# Completed: Step 2: Sorting by QNAME and performing 'fixmate' to update flags
#
#
# Started: Step 3: Sorting by coordinate in preparation for further filtering
# [1] 55146
# – Sorting by coordinate... [bam_sort_core] merging from 36 files and 6 in-memory blocks...
# \ Sorting by coordinate... [1]  + 55146 done       samtools sort -@ "${parallelize}" "${bam_fixmate_1}" > "${bam_sort_c_1}"
# [1] 57436
# – Running flagstat... [1]  + 57436 done       samtools flagstat -@ "${parallelize}" "${bam_sort_c_1}" > "${flagstat_3}"
# Completed: Step 3: Sorting by coordinate in preparation for further filtering
#
#
# Started: Step 4: Filtering for mates: '-f 3 -F 12'
# [1] 57510
# – Filtering for mates: '-f 3 -F 12'... [1]  + 57510 done       samtools view -@ "${parallelize}" -h -b -f 3 -F 12 "${bam_sort_c_1}" -o
# [1] 57896
# – Running flagstat... [1]  + 57896 done       samtools flagstat -@ "${parallelize}" "${bam_f3_F12}" > "${flagstat_4}"
# Completed: Step 4: Filtering for mates: '-f 3 -F 12'
#
#
# Started: Step 5: Sorting by QNAME and performing 'fixmate' to update flags
# [1] 57967
# – Sorting by QNAME... [bam_sort_core] merging from 36 files and 6 in-memory blocks...
# – Sorting by QNAME... [1]  + 57967 done       samtools sort -@ "${parallelize}" -n "${bam_f3_F12}" > "${bam_sort_n_2}"
# [1] 62202
# – Performing 'fixmate' to update flags... [1]  + 62202 done       samtools fixmate -@ "${parallelize}" "${bam_sort_n_2}" "${bam_fixmate_2}"
# [1] 62627
# | Running flagstat... [1]  + 62627 done       samtools flagstat -@ "${parallelize}" "${bam_fixmate_2}" > "${flagstat_5}"
# Completed: Step 5: Sorting by QNAME and performing 'fixmate' to update flags
#
#
# Started: Step 6: Sorting by coordinate in preparation for splitting by chromosomes
# [1] 62696
# \ Sorting by coordinate... [bam_sort_core] merging from 36 files and 6 in-memory blocks...
# – Sorting by coordinate... [1]  + 62696 done       samtools sort -@ "${parallelize}" "${bam_fixmate_2}" > "${bam_sort_c_2}"
# [1] 64438
# | Running flagstat... [1]  + 64438 done       samtools flagstat -@ "${parallelize}" "${bam_sort_c_2}" > "${flagstat_6}"
# Completed: Step 6: Sorting by coordinate in preparation for splitting by chromosomes
#
#
# Started: Step 7: Indexing bam, then splitting by chromosome(s)
# [1] 64509
# / Indexing bam... [1]  + 64509 done       samtools index -@ "${parallelize}" "${bam_sort_c_2}"
# [1] 64574
# – Splitting bam by chromosome(s)... [1]  + 64574 done       samtools view -@ "${parallelize}" -h -b "${bam_sort_c_2}" "chr11" -o
# [1] 64598
# / Running flagstat... [1]  + 64598 done       samtools flagstat -@ "${parallelize}" "${bam_sort_c_2/.bam/.chr11.bam}" >
# Completed: Step 7: Indexing bam, then splitting by chromosome(s)
#
#
# Started: Step 8: Indexing and creating header.txt for split chromosome bam(s)
# [1] 64603
# / Indexing split chromosome bam(s)... [1]  + 64603 done       samtools index -@ "${parallelize}" "${bam_sort_c_2/.bam/.chr11.bam}"
# [1] 64608
# | Creating header.txt for split chromosome bam(s)... [1]  + 64608 done       samtools view -@ "${parallelize}" -H "${bam_sort_c_2}" -o
# Completed: Step 8: Indexing and creating header.txt for split chromosome bam(s)
#
#
# Completed: /bin/zsh
# /bin/zsh run time: 2005 seconds.

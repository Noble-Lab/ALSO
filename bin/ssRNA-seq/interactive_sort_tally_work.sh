#!/bin/bash

#  interactive_sort_tally_work.sh
#  KA


#  Run everything below interactively...
qlogin -l mfree=2G -pe serial 8


#  Source module, functions ---------------------------------------------------
calculate_run_time() {
    # Calculate run time for chunk of code
    #
    # :param 1: start time in $(date +%s) format
    # :param 2: end time in $(date +%s) format
    # :param 3: message to be displayed when printing the run time (chr)
    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo ""
    echo "${3}"
    printf 'Run time: %dh:%dm:%ds\n' \
    $(( run_time/3600 )) $(( run_time%3600/60 )) $(( run_time%60 ))
    echo ""
}


display_spinning_icon() {
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


sort_bam_qname_samtools_auto() {
    # Run samtools sort -n on bam infile; outfile name is derived from infile
    # name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort -n on $(basename "${2}")."
}


list_tally_qnames_gzip() {
    # List and tally QNAMEs in a QNAME-sorted bam infile; function acts on a
    # bam infile to perform commands (samtools view, cut, uniq -c, sort -nr)
    # that list and tally QNAMEs; function writes the results to a txt
    # outfile, the name of which is derived from the txt infile
    #
    # :param 1: name of QNAME-sorted bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -f 1 \
        | uniq -c \
        > "${1/.bam/.QNAME.tmp.uniq.txt}"
    
    #  Do a parallel sort by number of records
    if [[ -f "${1/.bam/.QNAME.tmp.uniq.txt}" ]]; then
        # sort -nr "${1/.bam/.QNAME.tmp.uniq.txt}" \
        parsort -nr "${1/.bam/.QNAME.tmp.uniq.txt}" \
        > "${1/.bam/.QNAME.tmp.sort-nr.txt}"
    else
        echo "$(basename "${1/.bam/.QNAME.tmp.uniq.txt}") not found."
        return 1
    fi
    
    #  Trim leading whitespaces
    if [[ -f "${1/.bam/.QNAME.tmp.sort-nr.txt}" ]]; then
        cut -c 7- "${1/.bam/.QNAME.tmp.sort-nr.txt}" > "${1/.bam/.QNAME.txt}"
    else
        echo "$(basename "${1/.bam/.QNAME.tmp.sort-nr.txt}") not found."
        return 1
    fi

    #  Remove temporary intermediate file, then gzip the outfile
    if [[ -f "${1/.bam/.QNAME.txt}" ]]; then
        rm -f "${1/.bam/.QNAME.tmp.txt}" \
        "${1/.bam/.QNAME.tmp.uniq.txt}" \
        "${1/.bam/.QNAME.tmp.sort-nr.txt}"
        gzip "${1/.bam/.QNAME.txt}"
    else
        echo "$(basename "${1/.bam/.QNAME.txt}") not found."
        return 1
    fi
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally QNAMEs in $(basename "${1}")."
}


module load samtools/1.14 parallel/20200922


#  Copy files -----------------------------------------------------------------
#  Need more to do samtools sort -n with more memory for these samples:
#+     - Berletch_PRJNA256188
#+         - ovary_rep.SRR1525412-13-14.mm11.Aligned.sortedByCoord.out.primary.bam
#+         - brain_rep_1.SRR1525404-05-06.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam
#+         - spleen_rep_1.SRR1525408-09-10.mm11.Aligned.sortedByCoord.out.primary.bam
#+         - ovary_rep.SRR1525412-13-14.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam

dir_4dn_data="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data"

cd "${dir_4dn_data}/Berletch_PRJNA256188/alignments_None/alignments_primary" &&
    {
        parallel -j 4 cp {} "${TMPDIR}" ::: \
        ovary_rep.SRR1525412-13-14.mm11.Aligned.sortedByCoord.out.primary.bam \
        brain_rep_1.SRR1525404-05-06.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        spleen_rep_1.SRR1525408-09-10.mm11.Aligned.sortedByCoord.out.primary.bam \
        ovary_rep.SRR1525412-13-14.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam
    } &
display_spinning_icon $! "Copying files... "

cd "${dir_4dn_data}/Berletch_unpublished/alignments_None/alignments_primary" &&
    {
        parallel -j 8 cp {} "${TMPDIR}" ::: \
        heart_S1_R1_001.AHCG2FBGX9.mm11.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.AHCG2FBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.AHG5WVBGX9.mm11.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.AHG5WVBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.AHTTCCBGX7.mm11.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.AHTTCCBGX7.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.combined.mm11.Aligned.sortedByCoord.out.primary.bam \
        heart_S1_R1_001.combined.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHCG2FBGX9.mm11.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHCG2FBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHG5WVBGX9.mm11.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHG5WVBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHTTCCBGX7.mm11.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.AHTTCCBGX7.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.combined.mm11.Aligned.sortedByCoord.out.primary.bam \
        kidney_S2_R1_001.combined.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam
    } &
display_spinning_icon $! "Copying files... "


#  Sort bam files by QNAME ----------------------------------------------------
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "." -type f -name "*.primary.bam" -print0)

echo "Working with..."
for i in "${infiles[@]}"; do echo "${i}"; done
echo ""

for i in "${infiles[@]}"; do sort_bam_qname_samtools_auto 8 "${i}"; done &
display_spinning_icon $! "Sorting bam files by QNAME... "

# | Sorting bam files by QNAME... [bam_sort_core] merging from 32 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHG5WVBGX9.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:18m:2s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 112 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.combined.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:56m:0s
#
# \ Sorting bam files by QNAME... [bam_sort_core] merging from 24 files and 8 in-memory blocks...
# – Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.AHG5WVBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:13m:32s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 96 files and 8 in-memory blocks...
# / Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.combined.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:49m:28s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 88 files and 8 in-memory blocks...
# – Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.combined.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:46m:44s
#
# | Sorting bam files by QNAME... [bam_sort_core] merging from 16 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on ovary_rep.SRR1525412-13-14.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:9m:33s
#
# / Sorting bam files by QNAME... [bam_sort_core] merging from 32 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.AHG5WVBGX9.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:17m:18s
#
# | Sorting bam files by QNAME... [bam_sort_core] merging from 16 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.AHTTCCBGX7.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:11m:28s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 48 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHCG2FBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:24m:55s
#
# \ Sorting bam files by QNAME... [bam_sort_core] merging from 24 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on spleen_rep_1.SRR1525408-09-10.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:12m:39s
#
# | Sorting bam files by QNAME... [bam_sort_core] merging from 24 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on brain_rep_1.SRR1525404-05-06.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:13m:2s
#
# \ Sorting bam files by QNAME... [bam_sort_core] merging from 24 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on ovary_rep.SRR1525412-13-14.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:6m:49s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 104 files and 8 in-memory blocks...
# – Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.combined.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:35m:7s
#
# – Sorting bam files by QNAME... [bam_sort_core] merging from 16 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on heart_S1_R1_001.AHTTCCBGX7.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:4m:56s
#
# / Sorting bam files by QNAME... [bam_sort_core] merging from 48 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHCG2FBGX9.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:15m:54s
#
# / Sorting bam files by QNAME... [bam_sort_core] merging from 32 files and 8 in-memory blocks...
# | Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHG5WVBGX9.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:9m:54s
#
# \ Sorting bam files by QNAME... [bam_sort_core] merging from 16 files and 8 in-memory blocks...
# \ Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHTTCCBGX7.mm11.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:5m:15s
#
# | Sorting bam files by QNAME... [bam_sort_core] merging from 8 files and 8 in-memory blocks...
# – Sorting bam files by QNAME...
# Run samtools sort -n on kidney_S2_R1_001.AHTTCCBGX7.SPRET-EiJ.Aligned.sortedByCoord.out.primary.bam.
# Run time: 0h:4m:22s


#  Extract and sort QNAMEs from bam files -------------------------------------
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find "." -type f -name "*.primary.sort-n.bam" -print0)

echo "Working with..."
for i in "${infiles[@]}"; do echo "${i}"; done
echo ""

for i in "${infiles[@]}"; do list_tally_qnames_gzip "${i}"; done
display_spinning_icon $! "Listing and tallying QNAMEs in QNAME-sorted bam files... "

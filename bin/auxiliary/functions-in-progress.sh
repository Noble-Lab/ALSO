#!/bin/bash


#  Functions ------------------------------------------------------------------
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


check_file_sorted() {
    # Check if file is sorted or not
    sort -C <(zcat -df "${1}") || printf "%s" "not "; echo "sorted" #&
    # display_spinning_icon $! \
    # "Checking if file $(basename "${1}") is sorted... "
}


count_lines_gzip() {
    # Count number of records in a gzipped file
    #
    # :param 1: gzipped file, including path (chr)
    zcat "${1}" | wc -l
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


echo_loop() { for i in "${@:-*}"; do echo "${i}"; done; }


filter_duplicate_qnames_gzip(){
    # Using an infile single-column list of QNAMEs, exclude QNAMEs from an
    # AS.txt.gz file
    #
    # :param 1: infile from which to remove lines, including path (chr)
    # :param 2: infile containing entries to filter out, including path (chr)
    # :param 3: filtered outfile, including path (chr)
    start="$(date +%s)"

    grep -v -f <(zcat -df "${2}") <(zcat -df "${1}") | gzip > "${3}" & \
    display_spinning_icon $! \
    "Removing duplicate QNAME lines from $(basename "${1}")... "
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Remove duplicate QNAME lines from $(basename "${1}")."
}


#QUESTION Depending on the order of inputs here, what is in the output?
find_set_complement() {
    # Find and list the set complement between AS.txt.gz files for samples #1
    # and #2; comparing samples #1 and #2, list elements unique to sample #1;
    # function acts on and outputs only the first column
    #
    # :param 1: sample #1 AS.txt.gz file
    # :param 2: sample #2 AS.txt.gz file
    # :param 3: outputs elements unique to sample #2
    start="$(date +%s)"

    grep -vxFf <(zcat -df "${2}" | cut -f1) <(zcat -df "${1}" | cut -f1) \
    | sort \
    | gzip \
    > "${3}" &
    display_spinning_icon $! \
    "Writing out and sorting set elements unique to $(basename "${1}") in comparison to $(basename "${2}")... "
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Write out and sort set elements unique to $(basename "${1}") in comparison to $(basename "${2}")."
}


find_set_intersection() {
    # Find and list the set intersection elements between AS.txt.gz files for
    # samples #1 and #2
    #
    # :param 1: sample #1 AS.txt.gz file
    # :param 2: sample #2 AS.txt.gz file
    # :param 3: outputs set intersection between samples #1 and #2
    start="$(date +%s)"

    join <(zcat -df "${1}") <(zcat -df "${2}") \
    | tr ' ' \\t \
    | sort \
    | gzip \
    > "${3}" &
    display_spinning_icon $! \
    "Writing out and sorting set intersection elements for $(basename "${1}") and $(basename "${2}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Write out and sort set intersection elements for $(basename "${1}") and $(basename "${2}")."
}  #FIXME Output should be sorted


head_10() { zcat -d "${1}" | head -10; }


list_qnames_to_cut() {
    # Find and list QNAMEs with more than one entry
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"

    echo "Started: Find and list QNAMEs with more than one entry for $(basename "${1}")."
    zcat -df "${1}" \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | sort -nr \
    | cut -c 7- \
    | awk '$1 > 1 {print $0}' \
    | tr ' ' \\t \
    | cut -f2 \
    | gzip \
    > "${2}"
    #TODO Get Spinning icon properly working

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: Find and list QNAMEs with more than one entry for $(basename "${1}")."
}


randomly_sample_lines_from_file() {
    # Randomly samples lines from an infile
    #
    # :param 1: number of lines to samples (int)
    # :param 2: AS.txt.gz infile, including path (chr)
    # :param 3: gzipped outfile, including path (chr)
    start="$(date +%s)"

    shuf -n "${1}" <(zcat -df "${2}") | sort | gzip > "${3}" & \
    display_spinning_icon $! \
    "Randomly sampling ${1} lines from $(basename "${2}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Randomly sample ${1} lines from $(basename "${2}")."
}


tail_10() { zcat -d "${1}" | tail -10; }


tally_qnames_gzip() {
    # Tally numbers of entries per qname
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"

    echo "Started: Tally numbers of entries per QNAME in $(basename "${1}")."
    zcat -df "${1}" \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | sort -nr \
    | cut -c 7- \
    | tr ' ' \\t \
    | gzip \
    > "${2}"
    #TODO Get Spinning icon properly working
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: Tally numbers of entries per QNAME in $(basename "${1}")."
}


tally_qnames_gt_1_gzip() {
    # Tally numbers of entries per qname; retain in list only qnames with more
    # than one entry
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"

    echo "Started: Tally numbers of entries per QNAME in $(basename "${1}"), retaining only those with more than one entry in an outlist."
    zcat -df "${1}" \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | sort -nr \
    | cut -c 7- \
    | awk '$1 > 1 {print $0}' \
    | tr ' ' \\t \
    | gzip \
    > "${2}"
    #TODO Get Spinning icon properly working

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: Tally numbers of entries per QNAME in $(basename "${1}"), retaining only those with more than one entry in an outlist."
}

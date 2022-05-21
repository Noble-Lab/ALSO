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
    # :param 1: sorted sample #1 AS.txt.gz file
    # :param 2: sorted sample #2 AS.txt.gz file
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
    # :param 1: sorted sample #1 AS.txt.gz file
    # :param 2: sorted sample #2 AS.txt.gz file
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


head_20() { zcat -d "${1}" | head -20; }


identify_qnames_updated() {
    # Using txt.gz infile from list_tally_qnames, create txt.gz outfiles for
    # QNAME == 2, QNAME > 2, or QNAME < 2; outfile names are derived from
    # infile name; function runs the following commands in succession:
    #    1. filter txt infile for 'QNAME == 2', 'QNAME > 2', or 'QNAME < 2'
    #       based on user specification; write results to new txt file
    #    2. create a txt-file list of QNAMEs without tallies from the txt file
    #       output by step 1
    #    3. if 'keep', keep the txt file output by step 1; if 'gzip', gzip the
    #       txt file output by step 1; if 'delete', rm the txt file output by
    #       step 1
    #
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: identify QNAMES as follows: 'eq', QNAME = 2; 'gt', QNAME > 2;
    #           'lt', QNAME < 2; for (chr; default: 'gt')
    # :param 3: gzipped txt infile, including path (chr)
    # :param 4: 'keep', 'gzip', or 'delete' txt output from step 1 (chr;
    #           default: 'delete')
    start="$(date +%s)"

    case "$(echo "${2:-"gt"}" | tr '[:upper:]' '[:lower:]')" in
        eq | e) \
            comp="\$1 == 2"
            str="eq"
            ;;
        gt | g) \
            comp="\$1 > 2"
            str="gt"
            ;;
        lt | l) \
            comp="\$1 < 2"
            str="lt"
            ;;
        *) \
            echo "Parameter 1 is not \"eq\", \"gt\", or \"lt\""
            echo "Setting parameter 1 to \"eq\""
            comp="'\$1 == 2'"
            str="eq"
            ;;
    esac

    # #  Step 1
    # # shellcheck disable=SC2016
    # zcat -dfk "${3}" > "${3%.gz}"
    # split -n l/"${1}" "${3%.gz}" "${TMPDIR}/_pawk"$$
    # rm "${3%.gz}"
    #
    # for i in "${TMPDIR}/_pawk"$$*; do
    #     awk "${comp}" "${i}" > "${i}.out" &
    # done
    # wait
    # cat "${TMPDIR}/_pawk"$$*.out > "${3/.txt.gz/.${str}.tally.txt}"
    # rm "${TMPDIR}/_pawk"$$*

    #  Step 1
    awk "${comp}" <(zcat -dk "${3}") > "${3/.txt.gz/.${str}.tally.txt}"

    #  Step 2
    cut -c 3- "${3/.txt.gz/.${str}.tally.txt}" | gzip > "${3/.txt.gz/.${str}.txt.gz}"

    #  Step 3
    case "$(echo "${4:-"delete"}" | tr '[:upper:]' '[:lower:]')" in
        keep | k) \
            :
            ;;
        gzip | g) \
            gzip "${3/.txt.gz/.${str}.tally.txt}"
            ;;
        delete | d) \
            rm "${3/.txt.gz/.${str}.tally.txt}"
            ;;
        *) \
            echo "Parameter 4 is not \"keep\", \"gzip\", or \"delete\""
            echo "Will delete (rm) $(basename "${3/.txt.gz/.${str}.tally.txt}")"
            rm "${3/.txt.gz/.${str}.tally.txt}"
            ;;
    esac

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}" \
    "List entries with \"QNAME ${str} 2\" for $(basename "${3}")."
}


include_qname_reads_picard() {
    # Filter a bam infile to exclude reads with QNAMEs listed in a txt file;
    # write the filtered results to a bam outfile
    #
    # :param 1: name of bam infile, including path (chr)
    # :param 2: name of txt QNAME list, including path (chr)
    # :param 3: name of bam outfile, including path (cannot be same as bam
    #           infile) (chr)
    # :param 4: use the picard.jar available on the GS grid system (logical)
    start="$(date +%s)"
    dir_picard="/net/gs/vol3/software/modules-sw/picard/2.26.4/Linux/CentOS7/x86_64"

    case "$(echo "${4}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            java -jar "${dir_picard}"/picard.jar FilterSamReads \
            I="${1}" \
            O="${3}" \
            READ_LIST_FILE="${2}" \
            FILTER="includeReadList"
            ;;
        false | f) \
            picard FilterSamReads \
            I="${1}" \
            O="${3}" \
            READ_LIST_FILE="${2}" \
            FILTER="includeReadList"
            ;;
        *) \
            echo "Exiting: Parameter 4 is not \"TRUE\" or \"FALSE\"."
            return 1
            ;;
    esac

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Exclude reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
}


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


list_tally_qnames_gzip_updated() {
    # List and tally QNAMEs in a QNAME-sorted bam infile; function acts on a
    # bam infile to perform commands (samtools view, cut, uniq -c, sort -nr)
    # that list and tally QNAMEs; function writes the results to a txt
    # outfile, the name of which is derived from the txt infile
    #
    # :param 1: name of bam infile, including path (chr)
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
        rm -f "${1/.bam/.QNAME.tmp.txt}" "${1/.bam/.QNAME.tmp.uniq.txt}" "${1/.bam/.QNAME.tmp.sort-nr.txt}"
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

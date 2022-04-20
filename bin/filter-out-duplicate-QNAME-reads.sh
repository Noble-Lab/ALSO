#!/bin/bash

#  filter-out-duplicate-QNAME-reads.sh
#  KA


#  Functions ------------------------------------------------------------------
calculateRunTime() {
    # Calculate run time for processes
    # 
    # :param 1: start time in $(date +%s) format
    # :param 2: end time in $(date +%s) format
    # :param 3: message to be displayed when printing the run time
    run_time="$(echo "${2}" - "${1}" | bc -l)"
    
    echo ""
    echo "${3}"
    echo "Run time: ${run_time} seconds."
    echo ""
}


checkDependency() {
    # Check if program is available in "${PATH}"; exit if not
    # 
    # :param 1: program to be checked (chr)
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            exit 1
        }
}


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


identifyDuplicateQNAMEs() {
    # Identify duplicate QNAME reads from bam infile; the function runs the
    # following commands in succession:
    #    1. piped commands (samtools view to awk) to isolate duplicate QNAME 
    #       read entries (>2) in a txt file
    #    2. create a txt-file list of QNAMEs from the above txt file
    #    3. gzip the txt file from step 1
    #    4. (optional) if TRUE, generate a bam file comprised of only duplicate
    #       QNAME reads
    #
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: write out bam outfile comprised of reads with >2 QNAMEs: "TRUE"
    #           or "FALSE" (logical)
    # :param 4: term for sample in use (chr)
    start="$(date +%s)"
    
    #  Step 1
    samtools view "${2}" \
    | cut -f 1 \
    | parsort \
    | uniq -c \
    | parsort -nr \
    | awk '$1 > 2 {print $0}' \
    > "${2/.bam/.multiple.txt}" &
    displaySpinningIcon $! "For sample ${4}, running samtools view pipeline on ${2}... "
    
    #  Step 2
    cut -c 6- "${2/.bam/.multiple.txt}" > "${2/.bam/.multiple-QNAME.txt}"
    
    #  Step 3
    gzip "${2/.bam/.multiple.txt}"
    
    #  Step 4
    case "$(echo "${3}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            samtools view -hN "${2/.bam/.multiple-QNAME.txt}" "${2}" \
            | samtools view -b - \
            > "${2/.bam/.multiple-QNAME.bam}" &
            displaySpinningIcon $! "For sample ${4}, running samtools view -N on ${2}... "
            ;;

        false | f) \
            :
            ;;

        *) \
            echo "Parameter 3 is not \"TRUE\" or \"FALSE\", so skipping step to filter ${2}"
            :
            ;;
    esac

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${4}, ${2}, identified and listed reads based on duplicate QNAME status"
}


excludeDuplicateQNAMEs() {
    # Filter out duplicate QNAME reads from bam infile; bam and txt outfile
    # names are derived bam infile name; the function runs the following
    # commands in succession:
    #    1. (optional) coordinate-sort the bam infile
    #    2. (optional) write out flagstat report for bam infile
    #    3. filter out duplicate QNAME reads from bam infile
    #    4. (optional) write out flagstat report for filtered bam
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr); outfile name is
    #           derived from infile name
    # :param 3: single-column txt file containing QNAMEs to exclude on each
    #           line, including path (chr)
    # :param 4: coordinate-sort the bam infile prior to running picard
    #           FilterSamReads: "TRUE" or "FALSE" (logical)
    # :param 5: write out flagstat txt file prior to running picard
    #           FilterSamReads: "TRUE" or "FALSE" (logical)
    # :param 6: write out flagstat txt file after running picard
    #           FilterSamReads: "TRUE" or "FALSE" (logical)
    # :param 7: term for sample in use (chr)
    start="$(date +%s)"

    case "$(echo "${4}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            #  Step 1 (optional)
            samtools sort -@ "${1}" "${2}" > "${2/.bam/.sort-c.bam}" &
            displaySpinningIcon $! "Sorting bam by coordinate in preparation for use with picard FilterSamReads... "

            #  Step 2
            case "$(echo "${5}" | tr '[:upper:]' '[:lower:]')" in
                true | t) \
                    samtools flagstat -@ "${1}" "${2/.bam/.sort-c.bam}" > "${2/.bam/.sort-c.flagstat.txt}" &
                    displaySpinningIcon $! "Running flagstat prior to use with picard FilterSamReads... "
                    ;;
                false | f) \
                    :
                    ;;
                *) \
                    echo "Stopping: Parameter 5 is not \"TRUE\" or \"FALSE\"."
                    return 1
                    ;;
            esac

            #  Step 3
            picard FilterSamReads \
            I="${2/.bam/.sort-c.bam}" \
            O="${2/.bam/.sort-c.filter.bam}" \
            READ_LIST_FILE="${3}" \
            FILTER="excludeReadList" &
            displaySpinningIcon $! "Using picard to filter bam file for duplicate QNAME entries... "

            #  Step 4 (optional)
            case "$(echo "${6}" | tr '[:upper:]' '[:lower:]')" in
                true | t) \
                    samtools flagstat -@ "${1}" "${2/.bam/.sort-c.filter.bam}" > "${2/.bam/.sort-c.filter.flagstat.txt}" &
                    displaySpinningIcon $! "Running flagstat following use with picard FilterSamReads... "
                    ;;
                false | f) \
                    :
                    ;;
                *) \
                    echo "Stopping: Parameter 6 is not \"TRUE\" or \"FALSE\"."
                    return 1
                    ;;
            esac
            ;;

        false | f) \
            #  Step 2 (optional)
            case "$(echo "${5}" | tr '[:upper:]' '[:lower:]')" in
                true | t) \
                    samtools flagstat -@ "${1}" "${2}" > "${2/.bam/flagstat.txt}" &
                    displaySpinningIcon $! "Running flagstat prior to use with picard FilterSamReads... "
                    ;;
                false | f) \
                    :
                    ;;
                *) \
                    echo "Stopping: Parameter 5 is not \"TRUE\" or \"FALSE\"."
                    return 1
                    ;;
            esac

            #  Step 3
            picard FilterSamReads \
            I="${2}" \
            O="${2/.bam/.filter.bam}" \
            READ_LIST_FILE="${3}" \
            FILTER="excludeReadList" &
            displaySpinningIcon $! "Using picard to filter bam file for duplicate QNAME entries... "

            #  Step 4 (optional)
            case "$(echo "${6}" | tr '[:upper:]' '[:lower:]')" in
                true | t) \
                    samtools flagstat -@ "${1}" "${2/.bam/.filter.bam}" > "${2/.bam/.filter.flagstat.txt}" &
                    displaySpinningIcon $! "Running flagstat following use with picard FilterSamReads... "
                    ;;
                false | f) \
                    :
                    ;;
                *) \
                    echo "Stopping: Parameter 6 is not \"TRUE\" or \"FALSE\"."
                    return 1
                    ;;
            esac
            ;;

        *) \
            echo "Stopping: Parameter 4 is not \"TRUE\" or \"FALSE\"."
            return 1
            ;;
    esac

    end="$(date +%s)"
    echo ""
    echo "For sample ${7}, ${2}, filtering out reads based on duplicate QNAME status"
    calculateRunTime "${start}" "${end}"
}



#  Check for necessary dependencies; exit if not found ------------------------
checkDependency parallel
checkDependency picard
checkDependency samtools


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Identify and filter out duplicate QNAME reads from bam infile; output..."
    echo ""
    echo ""
    echo "Preprocessing is made up of the following steps:"
    # echo "    1. QNAME-sort bam file output by the Shendure-Lab pipeline (this"
    # echo "       is not necessary but speeds up the following step"
    # echo "       considerably)"
    echo "    2. Identify QNAMEs with >2 entries in the bam file (this step is"
    echo "       slow; there's likely room for optimization)"
    echo "    3. Create a bam file comprised of only duplicated QNAME entries"
    echo "       (this step is optional; it's not strictly necessary)"
    echo "    4. Filter bam file to exclude QNAMEs with >2 entries by"
    echo "       (a) coordinate-sorting the QNAME-sorted bam file and"
    echo "       (b) filtering the coordinate-sorted bam file with picard"
    echo "           FilterSamReads (picard FilterSamReads takes only"
    echo "           coordinate-sorted bam files as input; also, picard"
    echo "           FilterSamReads is exponentially faster than filtering"
    echo "           with grep)"
    # echo "    5. QNAME-sort the filtered bam file; then, to update flag"
    # echo "       information in the bam file, perform samtools fixmate on the"
    # echo "       bam file"
    # echo "    6. Prior to filtering out reads based on pairing status and MAPQ"
    # echo "       values, sort by coordinate again"
    # echo "    7. Filter out reads based on pairing status and MAPQ: Run" 
    # echo "       samtools view with flags -f 3 -F 12 -q 30"
    # echo "    8. Sort bam by QNAME and perform samtools fixmate to update flag"
    # echo "       information again"
    # echo "    9. Coordinate-sort and index the processed bam file, which"
    # echo "       should now be ready for the subsequent module"
    echo ""
    echo "Dependencies:"
    echo "    - parallel >= 20200101"
    echo "    - picard >= 2.26.4"
    echo "    - samtools >= 1.13"
    echo ""
    echo "Arguments:"
    echo "    -h <print this help message and exit>"
    echo "    -u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
    echo "    -i <bam infile, including path (chr)>"
    echo "    -o <directory for outfile, including path (chr)>"
    echo "    -p <number of cores for parallelization (int >= 1); default: 1>"
    echo ""
    exit
}


while getopts "h:u:i:o:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${parallelize}" ]] && parallelize=6

bam="${1:-"Disteche_sample_7.dedup.bam"}"



#  Scraps ---------------------------------------------------------------------
# printUsage() {
#     echo ""
#     echo "${0}:"
#     echo "Identify and filter out duplicate QNAME reads from bam infile; output..."
#     echo ""
#     echo ""
#     echo "Preprocessing is made up of the following steps:"
#     # echo "    1. QNAME-sort bam file output by the Shendure-Lab pipeline (this"
#     # echo "       is not necessary but speeds up the following step"
#     # echo "       considerably)"
#     echo "    2. Identify QNAMEs with >2 entries in the bam file (this step is"
#     echo "       slow; there's likely room for optimization)"
#     echo "    3. Create a bam file comprised of only duplicated QNAME entries"
#     echo "       (this step is optional; it's not strictly necessary)"
#     echo "    4. Filter bam file to exclude QNAMEs with >2 entries by"
#     echo "       (a) coordinate-sorting the QNAME-sorted bam file and"
#     echo "       (b) filtering the coordinate-sorted bam file with picard"
#     echo "           FilterSamReads (picard FilterSamReads takes only"
#     echo "           coordinate-sorted bam files as input; also, picard"
#     echo "           FilterSamReads is exponentially faster than filtering"
#     echo "           with grep)"
#     # echo "    5. QNAME-sort the filtered bam file; then, to update flag"
#     # echo "       information in the bam file, perform samtools fixmate on the"
#     # echo "       bam file"
#     # echo "    6. Prior to filtering out reads based on pairing status and MAPQ"
#     # echo "       values, sort by coordinate again"
#     # echo "    7. Filter out reads based on pairing status and MAPQ: Run" 
#     # echo "       samtools view with flags -f 3 -F 12 -q 30"
#     # echo "    8. Sort bam by QNAME and perform samtools fixmate to update flag"
#     # echo "       information again"
#     # echo "    9. Coordinate-sort and index the processed bam file, which"
#     # echo "       should now be ready for the subsequent module"
#     echo ""
#     echo "Dependencies:"
#     echo "    - parallel >= 20200101"
#     echo "    - picard >= 2.26.4"
#     echo "    - samtools >= 1.13"
#     echo ""
#     echo "Arguments:"
#     echo "    -h <print this help message and exit>"
#     echo "    -u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
#     echo "    -i <bam infile, including path (chr)>"
#     echo "    -o <directory for outfile, including path (chr)>"
#     echo "    -p <number of cores for parallelization (int >= 1); default: 1>"
#     echo ""
#     exit
# }
#
#
# while getopts "h:u:i:o:p:" opt; do
#     case "${opt}" in
#         h) printUsage ;;
#         u) safe_mode="${OPTARG}" ;;
#         i) infile="${OPTARG}" ;;
#         o) outpath="${OPTARG}" ;;
#         p) parallelize="${OPTARG}" ;;
#         *) printUsage ;;
#     esac
# done
#
# [[ -z "${safe_mode}" ]] && printUsage
# [[ -z "${infile}" ]] && printUsage
# [[ -z "${outpath}" ]] && printUsage
# [[ -z "${parallelize}" ]] && parallelize=6

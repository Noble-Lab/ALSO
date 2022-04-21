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


#TODO Split function in two?
evaluateAndFilterReadsForQnames() {
    # Identify duplicate qname reads from bam infile; name(s) of bam outfile(s)
    # is(are) derived from the name of the bam infile; the function runs the
    # following steps in succession:
    #    1. piped commands (samtools view, cut, parsort, uniq, parsort) to
    #       extract all qname read entries into a txt file
    #    2. using the txt file from step 1, create a txt-file list of all
    #       qnames equal to 2
    #    3. using the txt file from step 2, create a txt file in which all
    #       information except qnames is trimmed away
    #    4. using the txt file from step 3, generate a bam file comprised of
    #       reads with with qname entries equals 2 (not greater than 2 or less
    #       than 2)
    #    5. (optional) generate (1) a bam file comprised of reads with >2 qname
    #       entries and (2) a bam file comprised of reads with <2 qname entries
    #    6. (optional) gzip the txt files output by this function
    #
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: create bam outfiles for (1) reads with >2 qnames and (2) reads
    #           with <2 qnames: "TRUE" or "FALSE" (logical)
    # :param 4: gzip the txt files output by this function: "TRUE" or "FALSE"
    #           (logical)
    start="$(date +%s)"
    
    #  Step 1
    samtools view "${2}" \
    | cut -f 1 \
    | parsort \
    | uniq -c \
    | parsort -nr \
    > "${2/.bam/.qname.txt}" &
    displaySpinningIcon $! "Running piped commands (samtools view, cut, parsort, uniq, parsort): $(basename "${2}")"

    #  (create txt file for qname = 2)
    #  Step 2
    # shellcheck disable=SC2016
    getQnameInParallel "${1}" \
    '$1 == 2' \
    "${2/.bam/.qname.txt}" \
    "${2/.bam/.qname-eq-2.txt}"
    
    #  Step 3
    cut -c 6- "${2/.bam/.qname-eq-2.txt}" \
    > "${2/.bam/.qname-eq-2.trim.txt}"

    #  Step 4
    samtools view -hN "${2/.bam/.qname-eq-2.trim.txt}" "${2}" \
    | samtools view -b - \
    > "${2/.bam/.qname-eq-2.bam}" &
    displaySpinningIcon $! "Running samtools view -hN: $(basename "${2}"), $(basename "${2/.bam/.qname-eq-2.trim.txt}")"
    
    #  Step 5 (optional)
    case "$(echo "${3}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            #  (create txt file for qname > 2)
            #  Step 2
            # shellcheck disable=SC2016
            getQnameInParallel "${1}" \
            '$1 > 2' \
            "${2/.bam/.qname.txt}" \
            "${2/.bam/.qname-gt-2.txt}"

            #  Step 3
            cut -c 6- "${2/.bam/.qname-gt-2.txt}" \
            > "${2/.bam/.qname-gt-2.trim.txt}"

            #  Step 4
            samtools view -hN "${2/.bam/.qname-gt-2.trim.txt}" "${2}" \
            | samtools view -b - \
            > "${2/.bam/.qname-gt-2.bam}" &
            displaySpinningIcon $! "Running samtools view -hN: $(basename "${2}"), $(basename "${2/.bam/.qname-gt-2.trim.txt}")"

            #  (create txt file for qname < 2)
            #  Step 2
            # shellcheck disable=SC2016
            getQnameInParallel "${1}" \
            '$1 < 2' \
            "${2/.bam/.qname.txt}" \
            "${2/.bam/.qname-lt-2.txt}"

            #  Step 3
            cut -c 6- "${2/.bam/.qname-lt-2.txt}" \
            > "${2/.bam/.qname-lt-2.trim.txt}"

            #  Step 4
            samtools view -hN "${2/.bam/.qname-lt-2.trim.txt}" "${2}" \
            | samtools view -b - \
            > "${2/.bam/.qname-lt-2.bam}" &
            displaySpinningIcon $! "Running samtools view -hN: $(basename "${2}"), $(basename "${2/.bam/.qname-lt-2.trim.txt}")"

            #  Step 6 (optional)
            case "$(echo "${4}" | tr '[:upper:]' '[:lower:]')" in
                true | t) \
                    gzip "${2/.bam/.qname-gt-2.txt}"
                    gzip "${2/.bam/.qname-gt-2.trim.txt}"
                    gzip "${2/.bam/.qname-lt-2.txt}"
                    gzip "${2/.bam/.qname-lt-2.trim.txt}"
                    ;;
                false | f) \
                    :
                    ;;
                *) \
                    echo "Parameter 4 is not \"TRUE\" or \"FALSE\", so skipping..."
                    :
                    ;;
            esac
            ;;
        false | f) \
            :
            ;;
        *) \
            echo "Parameter 3 is not \"TRUE\" or \"FALSE\", so skipping..."
            :
            ;;
    esac

    #  Step 6 (optional)
    case "$(echo "${4}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            gzip "${2/.bam/.qname.txt}"
            gzip "${2/.bam/.qname-eq-2.txt}"
            gzip "${2/.bam/.qname-eq-2.trim.txt}"
            ;;
        false | f) \
            :
            ;;
        *) \
            echo "Parameter 4 is not \"TRUE\" or \"FALSE\", so skipping..."
            :
            ;;
    esac

    end="$(date +%s)"
    echo ""
    echo "Identified and listed reads based on duplicate qname status for $(basename "${2}")"
    calculateRunTime "${start}" "${end}"
}


evaluateGreaterEqual() {
    # Sort two values such that the greater is placed ahead of the lesser
    # 
    # :param 1: first value
    # :param 2: second value
    printf '%s\n%s\n' "${2}" "${1}" | sort -V -C
}


getQnameInParallel() {
    # Select qname entries with a awk-comparison string input by the user,
    # e.g., '$1 == 2' or '$1 > 2' by splitting txt infile into chunks,
    # processing the chunks with awk in parallel ("pawk"), then outputting a
    # txt file for qname entries < 2
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: awk evaluation for a field, e.g., '$1 == 2'
    # :param 3: txt infile, including path (chr)
    # :param 4: txt outfile, including path (chr)
    gsplit -n l/"${1}" "${3}" /tmp/_pawk$$
    # shellcheck disable=SC2231
    for i in /tmp/_pawk$$*; do
        awk "${2}" "${i}" > "${i}.out" &
    done
    wait
    cat /tmp/_pawk$$*.out > "${4}"
    rm /tmp/_pawk$$*
}


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Identify and filter out duplicate qname reads from bam infile; output..."
    echo ""
    echo ""
    echo "Preprocessing is made up of the following steps:"
    echo "    1. Evaluate qnames in the bam file, reporting those with exactly"
    echo "       2 entries, less than 2 entries (optional), and more than 2"
    echo "       entries (optional) in txt files"
    echo "    2. Using txt files generated in step 1, create bam file(s)"
    echo "       comprised of only qname = 2 read entries, qname < 2 read"
    echo "       entries (optional), qname > 2 read entries (optional)"
    echo ""
    echo "Dependencies:"
    echo "    - split (GNU coreutils) >= 8.32"
    echo "    - parallel >= 20200101"
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


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Check for necessary dependencies; exit if not found
checkDependency gsplit
checkDependency parallel
checkDependency samtools
evaluateGreaterEqual "$(parallel --version | head -1 | cut -d" " -f3)" "20200101" ||
    {
        echo -e "Exiting: GNU Parallel version must be from 2020 or later."
        exit 1
    }

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-u: Safe mode is on."
        set -Eeuxo pipefail
        ;;
    false | f) \
        echo -e "-u: Safe mode is off."
        :
        ;;
    *) \
        echo -e "Exiting: -u safe-mode argument must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

#  Check that "${infile}" exists
[[ -f "${infile}" ]] ||
    {
        echo -e "Exiting: -i ${infile} does not exist.\n"
        exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo "-o: Directory ${outpath} does not exist; making the directory."
        mkdir -p "${outpath}"
    }

#  Check "${parallelize}"
[[ ! "${parallelize}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -p parallelize argument must be an integer.\n"
        exit 1
    }

[[ ! $((parallelize)) -ge 1 ]] &&
    {
        echo -e "Exiting: -p parallelize argument must be an integer >= 1.\n"
        exit 1
    }


#  Process bam infile ---------------------------------------------------------
evaluateAndFilterReadsForQnames \
"${parallelize}" \
"${infile}" \
"TRUE" \
"TRUE"


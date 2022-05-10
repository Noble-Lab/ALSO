#!/bin/bash

#  test-preprocessing-gangliuw_2022-0416.sh
#  KA


qlogin -l mfree=2G -pe serial "${parallelize}"
pipeline-test_env


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
    # Identify duplicate QNAME reads from bam file; functions run the following
    # commands in succession:
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
    "For sample ${4}, ${2}, filtering out reads based on duplicate QNAME status"
}


# filterReadsWithPicard() {
#     # picard FilterSamReads \
#     # I="${2/.bam/.sort-c.rm-singletons.bam}" \
#     # O="${2/.bam/.sort-c.rm-singletons.filter.bam}" \
#     # READ_LIST_FILE="${2/.bam/.multiple-QNAME.txt}" \
#     # FILTER="excludeReadList" &
#     # displaySpinningIcon $! "For sample ${4}, using picard to filter ${2/.bam/.sort-c.rm-singletons.bam}... "
#
#     # picard FilterSamReads \
#     # I="${2/.bam/.sort-c.bam}" \
#     # O="${2/.bam/.sort-c.filter.bam}" \
#     # READ_LIST_FILE="${2/.bam/.multiple-QNAME.txt}" \
#     # FILTER="excludeReadList" &
#     # displaySpinningIcon $! "For sample ${4}, using picard to filter ${2/.bam/.sort-c.bam}... "
#     # #  SAM validation error: ERROR::INVALID_FLAG_FIRST_OF_PAIR:Record 728935,
#     # #+ Read name CTGGATTGGAATCGGCTTGAGAAGAGTATTCTTAAGTTCG:595833476, First of
#     # #+ pair flag should not be set for unpaired read.
#     # #+ 
#     # #+ Not sure of what to make of this...
#     #
#     # # rm "${2/.bam/.sort-c.bam}"
#     #
#     # samtools sort -n -@ "${1}" "${2/.bam/.sort-c.filter.bam}" > "${2/.bam/.sort-c.filter.sort-n.bam}" &
#     # displaySpinningIcon $! "For sample ${4}, running samtools sort -n on ${2/.bam/.sort-c.filter.bam}... "
#     #
#     # samtools fixmate -@ "${1}" "${2/.bam/.sort-c.filter.sort-n.bam}" "${2/.bam/.sort-c.filter.sort-n.fixmate.bam}" &
#     # displaySpinningIcon $! "For sample ${4}, running samtools fixmate on ${2/.bam/.sort-c.filter.sort-n.bam}... "
#     #
#     # # rm "${2/.bam/.sort-c.filter.sort-n.bam}"
#     #
#     # end="$(date +%s)"
#     # calculateRunTime \
#     # "${start}" \
#     # "${end}" \
#     # "For sample ${4}, generate bam outfiles taking into account QNAMEs with >2 entries"
# }


indexBamFile() {
    # Run samtools index on bam infile; total run time for running the function
    # is reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Term for sample in use (chr)
    start="$(date +%s)"

    samtools index -@ "${1}" "${2}" &
    displaySpinningIcon $! "Running samtools index for ${3} sample, ${2}... "

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${3}, indexing ${2}"
}


runFixmate() {
    # Run samtools sort -n and fixmate on bam infile; total run time for
    # running the function is reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Term for sample in use (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}" &
    displaySpinningIcon $! "Running samtools sort -n... "

    samtools fixmate -@ "${1}" "${2/.bam/.sort-n.bam}" "${2/.bam/.sort-n.fixmate.bam}" &
    displaySpinningIcon $! "Running samtools fixmate on the filtered, QNAME-sorted bam file... "

    rm "${2/.bam/.sort-n.bam}"

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${3}, running samtools sort -n and samtools fixmate on ${2}"
}


runFlagstat() {
    # Run samtools flagstat on a bam infile; total run time for running the
    # function is reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Name of flagstat outfile, including path (chr)
    # :param 4: Term for sample in use (chr)
    start="$(date +%s)"

    samtools flagstat -@ "${1}" "${2}" > "${3}" &
    displaySpinningIcon $! "Running samtools flagstat for ${4} sample, ${2}... "

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${4}, ${2}, generating flag statistics"
}


splitBamByChromosome() {
    # Run samtools view on a bam infile to split it by chromosome; bam outfile
    # name is derived from infile and chromosome names; total run time for
    # running the function is reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Name of chromosome to split out from bam infile (chr)  #TODO Add check
    # :param 4: Term for sample in use (chr)
    start="$(date +%s)"

    samtools view -@ "${parallelize}" -bh "${2}" "${3}" \
    > "${2/.bam/.${3}.bam}" &
    displaySpinningIcon $! "For ${4} sample, running samtools view to create ${2/.bam/.${3}.bam}... "

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${4}, ${2}, creating bam file for ${3}"
}


sortBamByCoordinateOverwriteInfile() {
    # Run samtools sort-by-coordinate on a bam infile; bam infile is
    # overwritten by sorted bam file; total run time for running the function
    # is reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Term for sample in use (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" -o "${2/.bam/.tmp.bam}" &
    displaySpinningIcon $! "Running samtools sort... "

    mv -f "${2/.bam/.tmp.bam}" "${2}" &
    displaySpinningIcon $! "${2} is being overwritten by ${2/.bam/.tmp.bam}... "

    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${3}, ${2}, running samtools sort"
}


sortBamByCoordinate() {
    # Run samtools sort-by-coordinate on a bam infile; bam outfile name is
    # derived from infile name; total run time for running the function is
    # reported too
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr); outfile name is derived from infile name
    # :param 3: Term for sample in use (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" > "${2/.bam/.sort-c.bam}" &
    displaySpinningIcon $! "For sample ${3}, running samtools sort on ${3}... "
    
    end="$(date +%s)"
    calculateRunTime \
    "${start}" \
    "${end}" \
    "For sample ${3}, ${2}, running samtools sort"
}


#  Check for necessary dependencies; exit if not found ------------------------
checkDependency parallel
checkDependency repair
checkDependency samtools


#  Variables, aliases ---------------------------------------------------------
dir_gangliuw_mm10="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-21/get_unique_fragments"
dir_gangliuw_CAST="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/results/2022-02-23/get_unique_fragments"
dir_results="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results"
dir_experiment="${dir_results}/$(date '+%Y-%m%d')_test-preprocessing-module"

bam="${1:-"Disteche_sample_7.dedup.bam"}"
parallelize="${2:-"6"}"

bam_CAST="${bam/.dedup.bam/.CAST.bam}"
bam_CAST_processed="${bam_CAST/.bam/.processed.bam}"
bam_CAST_processed_flagstat="${bam_CAST_processed/.bam/.flagstat.txt}"

bam_mm10="${bam/.dedup.bam/.mm10.bam}"
bam_mm10_processed="${bam_mm10/.bam/.processed.bam}"
bam_mm10_processed_flagstat="${bam_mm10_processed/.bam/.flagstat.txt}"


#  Set up TMPDIR, experiment directory ---------------------------------------
cp "${dir_gangliuw_mm10}/${bam}" "${TMPDIR}" && mv "${TMPDIR}/${bam}" "${TMPDIR}/${bam_CAST}"
cp "${dir_gangliuw_CAST}/${bam}" "${TMPDIR}" && mv "${TMPDIR}/${bam}" "${TMPDIR}/${bam_mm10}"

cp "${TMPDIR}/${bam_CAST}" "${TMPDIR}/bak.${bam_CAST}"
cp "${TMPDIR}/${bam_mm10}" "${TMPDIR}/bak.${bam_mm10}"

cd "${TMPDIR}" || ! echo "cd into \"\${TMPDIR}\" failed. Look into this."

[[ -d "${dir_experiment}" ]] ||
    {
        echo "Directory ${dir_experiment} does not exist; making the directory."
        mkdir -p "${dir_experiment}"
    }


#  Generate test datasets -----------------------------------------------------
#  CAST ---------------------------------------------------
sortBamByCoordinateOverwriteInfile \
"${parallelize}" \
"${bam_CAST}" \
"CAST"

indexBamFile \
"${parallelize}" \
"${bam_CAST}" \
"CAST"


#  Test dataset #1 ----------
filterReads \
"${parallelize}" \
"${bam_CAST}" \
"${bam_CAST_processed}" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed}" \
"${bam_CAST_processed_flagstat}" \
"CAST"

indexBamFile \
"${parallelize}" \
"${bam_CAST_processed}" \
"CAST"

splitBamByChromosome \
"${parallelize}" \
"${bam_CAST_processed}" \
"chr1" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.chr1.bam}" \
"${bam_CAST_processed/.bam/.chr1.flagstat.txt}" \
"CAST"

splitBamByChromosome \
"${parallelize}" \
"${bam_CAST_processed}" \
"chr16" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.chr16.bam}" \
"${bam_CAST_processed/.bam/.chr16.flagstat.txt}" \
"CAST"


#  Test dataset #2 ----------
runFixmate \
"${parallelize}" \
"${bam_CAST_processed}" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.bam}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.flagstat.txt}" \
"CAST"

sortBamByCoordinate \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.bam}" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.flagstat.txt}" \
"CAST"

indexBamFile \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"CAST"

splitBamByChromosome \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"chr1" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.chr1.bam}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.chr1.flagstat.txt}" \
"CAST"

splitBamByChromosome \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"chr16" \
"CAST"

runFlagstat \
"${parallelize}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.chr16.bam}" \
"${bam_CAST_processed/.bam/.sort-n.fixmate.sort-c.chr16.flagstat.txt}" \
"CAST"


#  mm10 ---------------------------------------------------
sortBamByCoordinateOverwriteInfile \
"${parallelize}" \
"${bam_mm10}" \
"mm10"

indexBamFile \
"${parallelize}" \
"${bam_mm10}" \
"mm10"


#  Test dataset #1 ----------
filterReads \
"${parallelize}" \
"${bam_mm10}" \
"${bam_mm10_processed}" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed}" \
"${bam_mm10_processed_flagstat}" \
"mm10"

indexBamFile \
"${parallelize}" \
"${bam_mm10_processed}" \
"mm10"

splitBamByChromosome \
"${parallelize}" \
"${bam_mm10_processed}" \
"chr1" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.chr1.bam}" \
"${bam_mm10_processed/.bam/.chr1.flagstat.txt}" \
"mm10"

splitBamByChromosome \
"${parallelize}" \
"${bam_mm10_processed}" \
"chr16" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.chr16.bam}" \
"${bam_mm10_processed/.bam/.chr16.flagstat.txt}" \
"mm10"


#  Test dataset #2 ----------
runFixmate \
"${parallelize}" \
"${bam_mm10_processed}" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.bam}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.flagstat.txt}" \
"mm10"

sortBamByCoordinate \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.bam}" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.flagstat.txt}" \
"mm10"

indexBamFile \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"mm10"

splitBamByChromosome \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"chr1" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.chr1.bam}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.chr1.flagstat.txt}" \
"mm10"

splitBamByChromosome \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.bam}" \
"chr16" \
"mm10"

runFlagstat \
"${parallelize}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.chr16.bam}" \
"${bam_mm10_processed/.bam/.sort-n.fixmate.sort-c.chr16.flagstat.txt}" \
"mm10"


#  To be determined -----------------------------------------------------------
# cp "${bam_CAST_processed_flagstat}" "${bam_mm10_processed_flagstat}" "${dir_experiment}"
# cp "${bam_CAST_processed}" "${bam_mm10_processed}" "${dir_experiment}"
# cp -- *.bam "${dir_experiment}"


#  Scraps ---------------------------------------------------------------------
# total_start="$(date +%s)"
# total_end="$(date +%s)"
# calculateRunTime \
# "${total_start}" \
# "${total_end}" \
# "For CAST and mm10 samples, calculing total run time for filtering out reads based on pairing status and MAPQ score assignment"
# #  - Previous: 78 seconds (2020 MacBook Pro M1); 173 seconds, 2.88 minutes (HPC)
# #+ - Today: 

#!/bin/bash


#TODO User-specified outfiles
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


check_dependency() {
    # Check if program is available in "${PATH}"; exit if not
    #
    # :param 1: program to be checked (chr)
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            exit 1
        }
}


combine_gz_qname_lists_return_unique_gzip() {
    # Use gzcat on txt.gz infiles to combine, sort, and select unique entries;
    # txt stream is gzipped
    #
    # :param 1: Set "TRUE" for setting with GNU coreutils installed (logical)
    # :param @: QNAME txt.gz infiles, including paths    
    # shellcheck disable=SC2002
    case "$(echo "${1}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            zcat "${@}" | sort -u | gzip
            ;;
        false | f) \
            gzcat "${@}" | sort -u | gzip
            ;;
        *) \
            echo "Exiting: Parameter 2 is not \"TRUE\" or \"FALSE\"."
            return 1
            ;;
    esac
}


combine_qname_lists() {
    # Use cat to combine txt or txt.gz files
    #
    # :param @: QNAME txt/txt.gz infiles, including paths
    cat "${@}"
}


combine_two_qname_lists() {
    # Use cat to combine two txt or txt.gz files into user-specified outfile
    #
    # :param 1: QNAME txt/txt.gz infile #1, including path (chr)
    # :param 2: QNAME txt/txt.gz infile #2, including path (chr)
    # :param 3: combined QNAME txt.gz outfile, including path (chr)
    cat "${1}" "${2}" > "${3}"
}


count_lines() {
    # Count number of records in a file
    #
    # :param 1: file, including path (chr)
    # shellcheck disable=SC2002
    cat "${1}" | wc -l
}


count_lines_bam() {
    # Count number of records in a bam file; print to console
    #
    # :param 1: bam infile, including path (chr)
    start="$(date +%s)"

    samtools view "${1}" | wc -l &
    display_spinning_icon $! \
    "Running samtools view to wc -l on $(basename "${1}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view on $(basename "${1}")."
}


count_lines_gzip() {
    # Count number of records in a file
    #
    # :param 1: Set "TRUE" for setting with GNU coreutils installed (logical)
    # :param 2: gzipped file, including path (chr)
    # shellcheck disable=SC2002
    case "$(echo "${1}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            zcat "${2}" | wc -l
            ;;
        false | f) \
            gzcat "${2}" | wc -l
            ;;
        *) \
            echo "Exiting: Parameter 2 is not \"TRUE\" or \"FALSE\"."
            return 1
            ;;
    esac
}


decompress_gzip() {
    # Decompress a gzipped infile without removing said infile
    #
    # :param 1: gzipped infile, including path (chr)
    gzip -dk "${1}"
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


exclude_qname_reads_picard() {
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
            FILTER="excludeReadList" &
            display_spinning_icon $! \
            "Running picard FilterSamReads with $(basename "${1}") filtered by $(basename "${2}")... "
            ;;
        false | f) \
            picard FilterSamReads \
            I="${1}" \
            O="${3}" \
            READ_LIST_FILE="${2}" \
            FILTER="excludeReadList" &
            display_spinning_icon $! \
            "Running picard FilterSamReads with $(basename "${1}") filtered by $(basename "${2}")... "
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


get_qname_parallel() {
    # Select QNAME entries with an awk-comparison string input by the user,
    # e.g., '$1 == 2' or '$1 > 2' by splitting txt infile into chunks,
    # processing the chunks with awk in parallel ("pawk"), then outputting a
    # txt file for QNAME entries < 2
    #
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: one of two options: 'tmp' to split files in /tmp, 'tmpdir' to
    #           split files in ${TMPDIR} (chr)
    # :param 3: awk evaluation for a field, e.g., '$1 == 2' (chr)
    # :param 4: txt infile, including path, from, e.g., list_tally_qnames
    #           (chr)
    # :param 5: txt outfile, including path (chr)
    start="$(date +%s)"

    case "$(echo "${2}" | tr '[:upper:]' '[:lower:]')" in
        tmp) \
            str="/tmp"  # e.g., for use on M1 MacBook Pro 2020
            gsplit -n l/"${1}" "${4}" "${str}/_pawk"$$
            ;;
        tmpdir) \
            str="${TMPDIR}"  # e.g., for use with GS HPC
            split -n l/"${1}" "${4}" "${str}/_pawk"$$
            ;;
        *) \
            echo "Exiting: Parameter 2 is not \"tmp\" or \"tmpdir\"."
            return 1
            ;;
    esac

    # shellcheck disable=SC2231
    for i in "${str}/_pawk"$$*; do
        awk "${3}" "${i}" > "${i}.out" &
    done
    wait
    cat "${str}/_pawk"$$*.out > "${5}"
    rm "${str}/_pawk"$$*

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}" \
    "Evaluate QNAME entries (${3}) with parallelized chunking strategy for $(basename "${4}")."
}


identify_qnames() {
    # Using txt infile from list_tally_qnames, create txt outfiles for
    # QNAME == 2; outfile names are derived from infile name; function runs the
    # following commands in succession:
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
    # :param 3: set where to store temporary split files: 'tmp', split files
    #           in /tmp; 'tmpdir', split files in ${TMPDIR}
    # :param 4: txt infile, including path (chr)
    # :param 5: 'keep', 'gzip', or 'delete' txt output from step 1 (chr;
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

    case "$(echo "${3}" | tr '[:upper:]' '[:lower:]')" in
        tmp) \
            dir_str="tmp"  # For use on local computer
            ;;
        tmpdir) \
            dir_str="tmpdir"  # For use with GS HPC
            ;;
        *) \
            echo "Exiting: Parameter 2 is not \"tmp\" or \"tmpdir\"."
            return 1
            ;;
    esac

    #  Step 1
    # shellcheck disable=SC2016
    get_qname_parallel "${1}" \
    "${dir_str}" \
    "${comp}" \
    "${4}" \
    "${4/.txt/.${str}.tally.txt}" &
    display_spinning_icon $! \
    "Listing QNAME entries in which \"QNAME ${str} 2\" for $(basename "${4}")... "

    #  Step 2
    cut -c 4- "${4/.txt/.${str}.tally.txt}" > "${4/.txt/.${str}.txt}"
    
    #  Step 3
    case "$(echo "${5:-"delete"}" | tr '[:upper:]' '[:lower:]')" in
        keep | k) \
            :
            ;;
        gzip | g) \
            gzip "${4/.txt/.${str}.tally.txt}"
            ;;
        delete | d) \
            rm "${4/.txt/.${str}.tally.txt}"
            ;;
        *) \
            echo "Parameter 5 is not \"keep\", \"gzip\", or \"delete\""
            echo "Will delete (rm) $(basename "${4/.txt/.${str}.tally.txt}")"
            rm "${4/.txt/.${str}.tally.txt}"
            ;;
    esac

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}" \
    "List entries with \"QNAME ${str} 2\" for $(basename "${4}")."
}


index_bam() {
    # Use samtools index to index a bam file; index outfile name is derived
    # from name of bam infile
    #
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools index -@ "${1}" "${2}" &
    display_spinning_icon $! \
    "Running samtools index on $(basename "${2}")... "

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}" \
    "Index $(basename "${2}")."
}


list_tally_qnames() {
    # List and tally QNAMEs in a bam infile; function acts on a bam infile to
    # perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
    # list and tally QNAMEs; function writes the results to a txt outfile, the
    # name of which is derived from the txt infile
    #
    # :param 1: name of bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools view "${1}" \
    | cut -f 1 \
    | sort \
    | uniq -c \
    | sort -nr \
    > "${1/.bam/.QNAME.tmp.txt}" &
    display_spinning_icon $! \
    "Running piped commands (samtools view, cut, sort, uniq -c, sort -nr) on $(basename "${1}")... "

    #  Trim leading whitespaces
    if [[ -f "${1/.bam/.QNAME.tmp.txt}" ]]; then
        cut -c 6- "${1/.bam/.QNAME.tmp.txt}" > "${1/.bam/.QNAME.txt}" &
        display_spinning_icon $! \
        "Trimming away leading whitespaces in $(basename "${1/.bam/.QNAME.tmp.txt}")... "
    else
        echo "$(basename "${1/.bam/.QNAME.tmp.txt}") not found."
        return 1
    fi

    #  Remove temporary intermediate file
    if [[ -f "${1/.bam/.QNAME.txt}" ]]; then
        rm "${1/.bam/.QNAME.tmp.txt}"
    else
        echo "$(basename "${1/.bam/.QNAME.txt}") not found."
        return 1
    fi
        
    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally QNAMEs in $(basename "${1}")."
}


list_tally_qnames_trans() {
    # List and tally interchromosomal QNAMEs in a bam infile; function acts on
    # a bam infile to perform piped commands (samtools view, awk, cut, sort,
    # uniq -c) that list and tally QNAMEs; function writes the results to a txt
    # outfile, the name of which is derived from the txt infile
    #
    # :param 1: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools view "${1}" \
    | awk '($3 != $7 && $7 != "=")' \
    | cut -f 1 \
    | sort \
    | uniq -c \
    > "${1/.bam/.trans-QNAME.txt}" &
    display_spinning_icon $! \
    "Running piped commands (samtools view, awk, cut, sort, uniq -c) on $(basename "${1}")... "

    # #  Trim leading whitespaces
    # if [[ -f "${1/.bam/.trans-QNAME.tmp.txt}" ]]; then
    #     cut -c 6- "${1/.bam/.trans-QNAME.tmp.txt}" > "${1/.bam/.trans-QNAME.txt}" &
    #     display_spinning_icon $! \
    #     "Trimming away leading whitespaces in $(basename "${1/.bam/.trans-QNAME.tmp.txt}")... "
    # else
    #     echo "$(basename "${1/.bam/.trans-QNAME.tmp.txt}") not found."
    #     return 1
    # fi
    #
    # #  Remove temporary intermediate file
    # if [[ -f "${1/.bam/.trans-QNAME.txt}" ]]; then
    #     rm "${1/.bam/.trans-QNAME.tmp.txt}"
    # else
    #     echo "$(basename "${1/.bam/.trans-QNAME.txt}") not found."
    #     return 1
    # fi

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "List and tally interchromosomal QNAMEs in $(basename "${1}")."
}


list_tally_qnames_trans_simple() {
    #TODO Description of function
    #
    # :param 1: 

    samtools view "${1}" \
    | awk '($3 != $7 && $7 != "=")' \
    > "${1/.bam/.trans-QNAME.txt}" &
    display_spinning_icon $! \
    "Isolating interchromosomal QNAMEs reads from $(basename "${1}")... "
}


perform_diff() {
    # Prints the lines in file #2 that do not match lines in file #1
    # 
    # :param 1: infile #1, including path (chr)
    # :param 2: infile #2, including path (chr)
    diff "${1}" "${2}"

}


perform_diff_reverse() {
    # Prints the lines in file #2 that exactly match lines in file #1; grep -F
    # searches for fixed string matches, -f uses file #1 as a list of grep
    # search patterns, and -x prints only lines matched in entirely
    #
    # :param 1: infile #1, including path (chr)
    # :param 2: infile #2, including path (chr)
    grep -Fxf "${1}" "${2}"
}


get_unique_records() {
    # Print records from file #1 that do not match records in file #2
    #
    # :param 1: infile #1, including path (chr)
    # :param 2: infile #2, including path (chr)
    grep -Fvf "${2}" "${1}"
}


remove_reads_low_quality() {
    # From a bam infile, remove reads with MAPQ < 30 and with flags for ‘read
    # unmapped’ (0x4) and ‘mate unmapped’ (0x8) while retaining reads with
    # flags for ‘mate paired’ (0x1) and ‘read mapped in proper pair’ (0x2);
    # user defines bam outfile name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: name of bam outfile, including path (chr)
    start="$(date +%s)"
    
    samtools view -@ "${1}" -h -b -f 3 -F 12 -q 30 "${2}" -o "${3}" &
    display_spinning_icon $! \
    "Running samtools view (-f 3 -F 12 -q 30) on $(basename "${2}")... "

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Filter out reads based on MAPQ scores and pairing status for $(basename "${2}")."
}


remove_reads_low_quality_auto() {
    # From a bam infile, remove reads with MAPQ < 30 and with flags for ‘read
    # unmapped’ (0x4) and ‘mate unmapped’ (0x8) while retaining reads with
    # flags for ‘mate paired’ (0x1) and ‘read mapped in proper pair’ (0x2);
    # outfile name is derived from infile name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools view -@ "${1}" -h -b -f 3 -F 12 -q 30 "${2}" -o "${2/.bam/.rm.bam}" &
    display_spinning_icon $! \
    "Running samtools view -f 3 -F 12 -q 30 on $(basename "${2}")... "

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Filter out reads based on MAPQ scores and pairing status for $(basename "${2}")."
}


repair_bam() {
    # Using Subread repair, order bam file contents such that read pairs are
    # together; user defines bam outfile name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: name of bam outfile, including path (chr)
    start="$(date +%s)"

    repair -d -T "${1}" -c -i "${2}" -o "${3}" &
    display_spinning_icon $! \
    "Running repair -d -c on $(basename "${2}")... "

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Order $(basename "${2}") such that pairs are together."
}


repair_bam_auto() {
    # Using Subread repair, order bam file contents such that read pairs are
    # together; bam outfile name is derived from the name of the bam infile
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    repair -d -T "${1}" -c -i "${2}" -o "${2/.bam/.repair.bam}" &
    display_spinning_icon $! \
    "Running repair -d -c on $(basename "${2}")... "

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Order $(basename "${2}") such that pairs are together."
}


retain_qname_reads_picard() {
    # Filter a bam infile to include reads with QNAMEs listed in a txt file;
    # write the filtered results to a bam outfile, the name and path of which
    # is user-specified
    # 
    # :param 1: name of bam infile, including path (chr)
    # :param 2: name of txt QNAME list, including path (chr)
    # :param 3: name of bam outfile, including path (chr; cannot be same as bam
    #           infile)
    # :param 4: use the picard.jar available on the GS grid system (logical)
    start="$(date +%s)"

    case "$(echo "${4}" | tr '[:upper:]' '[:lower:]')" in
        true | t) \
            java -jar "${dir_picard}"/picard.jar FilterSamReads \
            I="${1}" \
            O="${3}" \
            READ_LIST_FILE="${2}" \
            FILTER="includeReadList" &
            display_spinning_icon $! \
            "Running picard FilterSamReads with $(basename "${1}") filtered by $(basename "${2}")"
            ;;
        false | f) \
            picard FilterSamReads \
            I="${1}" \
            O="${3}" \
            READ_LIST_FILE="${2}" \
            FILTER="includeReadList" &
            display_spinning_icon $! \
            "Running picard FilterSamReads with $(basename "${1}") filtered by $(basename "${2}")"                        
            ;;
        *) \
            echo "Exiting: Parameter 4 is not \"TRUE\" or \"FALSE\"."
            return 1
            ;;
    esac

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Retain reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
}


retain_qname_reads_samtools() {
    # Filter a bam infile to include reads with QNAMEs listed in a txt file;
    # write the filtered results to a bam outfile, the name and path of which
    # is user-specified
    # 
    # :param 1: name of bam infile, including path (chr)
    # :param 2: name of txt QNAME list, including path (chr)
    # :param 3: name of bam outfile, including path (chr; cannot be same as bam
    #           infile)
    start="$(date +%s)"
    
    samtools view -hN "${2}" "${1}" \
    | samtools view -b - \
    > "${3}" &
    display_spinning_icon $! \
    "Running samtools view -hN with $(basename "${1}") filtered by $(basename "${2}")"

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Retain reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
}


run_flagstat() {
    # Run samtools flagstat on a bam infile; txt outfile name is derived from
    # the name of the bam infile
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools flagstat -@ "${1}" "${2}" > "${2/.bam/.flagstat.txt}" &
    display_spinning_icon $! \
    "Running samtools flagstat for $(basename "${2}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Generate flag statistics for $(basename "${2}")."
}


sort_bam_coordinate_samtools() {
    # Run samtools sort on a bam infile; user defines bam outfile
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Name of bam outfile, including path (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" > "${3}" &
    display_spinning_icon $! \
    "Running samtools sort (by coordinate) on $(basename "${2}")... "
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort (by coordinate) on $(basename "${2}")."
}


sort_bam_coordinate_samtools_auto() {
    # Run samtools sort on a bam infile; bam outfile name is derived from the
    # bam infile name
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" > "${2/.bam/.sort-c.bam}" &
    display_spinning_icon $! "Running samtools sort (by coordinate) on $(basename "${2}")... "
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort (by coordinate) on $(basename "${2}")."
}


sort_bam_coordinate_samtools_overwrite_infile() {
    # Run samtools sort on a bam infile; bam infile is overwritten by sorted
    # bam file
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" -o "${2/.bam/.tmp.bam}" &
    display_spinning_icon $! \
    "Running samtools sort... "

    mv -f "${2/.bam/.tmp.bam}" "${2}" &
    display_spinning_icon $! \
    "${2} is being overwritten by ${2/.bam/.tmp.bam}... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort (by coordinate) on $(basename "${2}")."
}


sort_bam_coordinate_picard() {
    # Run picard samsort SORT_ORDER="coordinate" on a bam infile; user defines
    # bam outfile
    # 
    # :param 1: Name of bam infile, including path (chr)
    # :param 2: Name of bam outfile, including path (chr)
    start="$(date +%s)"

    picard SortSam \
    I="${1}" \
    O="${2}" \
    SORT_ORDER="coordinate" &
    display_spinning_icon $! "Running picard SortSam SORT_ORDER=\"coordinate\" on $(basename "${1}")... "
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run picard SortSam SORT_ORDER=\"coordinate\" on $(basename "${1}")."
}  #UNTESTED


sort_bam_qname_samtools() {
    # Run samtools sort -n on bam infile; user defines bam outfile
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: name of bam outfile, including path (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${3}" &
    display_spinning_icon $! "Running samtools sort -n on $(basename "${2}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort -n on $(basename "${2}")."
}


sort_bam_qname_samtools_auto() {
    # Run samtools sort -n on bam infile; outfile name is derived from infile
    # name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}" &
    display_spinning_icon $! "Running samtools sort -n on $(basename "${2}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort -n on $(basename "${2}")."
}


sort_bam_qname_fixmate_auto() {
    # Run samtools sort -n and fixmate on bam infile; outfile name is derived
    # from infile name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}" &
    display_spinning_icon $! "Running samtools sort -n on $(basename "${2}")... "

    samtools fixmate -@ "${1}" "${2/.bam/.sort-n.bam}" "${2/.bam/.fixmate.bam}" &
    display_spinning_icon $! "Running samtools fixmate on $(basename "${2/.bam/.sort-n.bam}")... "

    rm "${2/.bam/.sort-n.bam}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort -n and samtools fixmate on $(basename "${2}")."
}


split_bam_chromosome() {
    # Run samtools view on a bam infile to split it by chromosome; bam outfile
    # name is derived from infile and chromosome names
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Name of chromosome to split out from bam infile (chr)
    #TODO Add check
    start="$(date +%s)"

    samtools view -@ "${1}" -bh "${2}" "${3}" \
    > "${2/.bam/.${3}.bam}" &
    display_spinning_icon $! "Running samtools view to create $(basename "${2/.bam/.${3}.bam}")... "

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view to create $(basename "${2/.bam/.${3}.bam}")."
}

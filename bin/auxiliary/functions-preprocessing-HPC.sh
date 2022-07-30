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
    # Use zcat on txt.gz infiles to combine, sort, and select unique entries;
    # txt stream is gzipped
    #
    # :param @: QNAME txt.gz infiles, including paths    
    # shellcheck disable=SC2002
    zcat "${@}" | sort -u | gzip
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

    samtools view "${1}" | wc -l

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view on $(basename "${1}")."
}


count_lines_gzip() {
    # Count number of records in a gzipped file
    #
    # :param 1: gzipped file, including path (chr)
    # shellcheck disable=SC2002
    zcat "${1}" | wc -l
}


decompress_gzip() {
    # Decompress a gzipped infile without removing the infile
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


duplicate_records_in_txt_gz() {
    # #TODO Description of function
    #
    # :param 1: gzipped text infile for which to duplicate each records
    # :param 2: gzipped text outfile in which each record is duplicated
    sed -e p <(zcat "${1}") | gzip > "${2}"
}


echo_completion_message() {
    # #TODO Description of function
    #
    # :param 1: step number (int)
    echo "Step ${1} completed; moving to next step..."
}


echo_exit_message() {
    # #TODO Description of function
    #
    # :param 1: step number (int)
    echo "Exiting: Step ${1}."
}


echo_flagstat() {
    # Echo a filename then print the first 25 lines of the file
    #
    # :param 1: file, e.g., a txt file (chr)
    echo ""
    echo "${1}"
    head -25 "${1}"
}


#TODO Write a function description
echo_loop() { for i in "${@:-*}"; do echo "${i}"; done; }


extract_n_lines_gzip_auto() {
    # Extract n number of records from a gzipped file
    #
    # :param 1: gzipped file, including path (chr)
    # :param 2: number of records to extract (int)
    # shellcheck disable=SC2002
    zcat "${1}" | head -n "${2}" | gzip > "${1/.txt.gz/.${2}.txt.gz}"
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
    get_qname_parallel \
    "${1}" \
    "${dir_str}" \
    "${comp}" \
    "${4}" \
    "${4/.txt/.${str}.tally.txt}"

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
    
    samtools index -@ "${1}" "${2}"

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
        > "${1/.bam/.QNAME.tmp.txt}"

    #  Trim leading whitespaces
    if [[ -f "${1/.bam/.QNAME.tmp.txt}" ]]; then
        cut -c 7- "${1/.bam/.QNAME.tmp.txt}" > "${1/.bam/.QNAME.txt}"
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


list_tally_qnames_gzip() {
    # List and tally QNAMEs in a bam infile; function acts on a bam infile to
    # perform piped commands (samtools view, cut, sort, uniq -c, sort -nr) that
    # list and tally QNAMEs; function writes the results to a txt outfile, the
    # name of which is derived from the txt infile; finally, the outfile is
    # gzipped
    #
    # :param 1: name of bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools view "${1}" \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | sort -nr \
        > "${1/.bam/.QNAME.tmp.txt}"

    #  Trim leading whitespaces
    if [[ -f "${1/.bam/.QNAME.tmp.txt}" ]]; then
        cut -c 7- "${1/.bam/.QNAME.tmp.txt}" > "${1/.bam/.QNAME.txt}"
    else
        echo "$(basename "${1/.bam/.QNAME.tmp.txt}") not found."
        return 1
    fi

    #  Remove temporary intermediate file, gzip outfile
    if [[ -f "${1/.bam/.QNAME.txt}" ]]; then
        rm "${1/.bam/.QNAME.tmp.txt}"
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
        > "${1/.bam/.trans-QNAME.txt}"

    # #  Trim leading whitespaces
    # if [[ -f "${1/.bam/.trans-QNAME.tmp.txt}" ]]; then
    #     cut -c 7- "${1/.bam/.trans-QNAME.tmp.txt}" > "${1/.bam/.trans-QNAME.txt}"
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
        > "${1/.bam/.trans-QNAME.txt}"
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
    
    samtools view -@ "${1}" -h -b -f 3 -F 12 -q 30 "${2}" -o "${3}"

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
    
    samtools view -@ "${1}" -h -b -f 3 -F 12 -q 30 "${2}" -o "${2/.bam/.rm.bam}"

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

    repair -d -T "${1}" -c -i "${2}" -o "${3}"

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

    repair -d -T "${1}" -c -i "${2}" -o "${2/.bam/.repair.bam}"

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Order $(basename "${2}") such that pairs are together."
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
        > "${3}"

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

    samtools flagstat -@ "${1}" "${2}" > "${3}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Generate flag statistics for $(basename "${2}")."
}


run_flagstat_auto() {
    # Run samtools flagstat on a bam infile; txt outfile name is derived from
    # the name of the bam infile
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"

    samtools flagstat -@ "${1}" "${2}" > "${2/.bam/.flagstat.txt}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Generate flag statistics for $(basename "${2}")."
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
    SORT_ORDER="coordinate"
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run picard SortSam SORT_ORDER=\"coordinate\" on $(basename "${1}")."
}  #UNTESTED


sort_bam_coordinate_samtools() {
    # Run samtools sort on a bam infile; user defines bam outfile
    # 
    # :param 1: Number of cores for parallelization (int >= 1)
    # :param 2: Name of bam infile, including path (chr)
    # :param 3: Name of bam outfile, including path (chr)
    start="$(date +%s)"

    samtools sort -@ "${1}" "${2}" > "${3}"
    
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

    samtools sort -@ "${1}" "${2}" > "${2/.bam/.sort-c.bam}"
    
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

    samtools sort -@ "${1}" "${2}" -o "${2/.bam/.tmp.bam}"

    mv -f "${2/.bam/.tmp.bam}" "${2}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools sort (by coordinate) on $(basename "${2}")."
}


sort_bam_qname_picard() {
    # Run picard samsort SORT_ORDER="coordinate" on a bam infile; user defines
    # bam outfile
    # 
    # :param 1: Name of bam infile, including path (chr)
    # :param 2: Name of bam outfile, including path (chr)
    start="$(date +%s)"

    picard SortSam \
    I="${1}" \
    O="${2}" \
    SORT_ORDER="queryname"
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run picard SortSam SORT_ORDER=\"queryname\" on $(basename "${1}")."
}


#LOOK
sort_bam_qname_picard_auto() {
    # Run picard samsort SORT_ORDER="coordinate" on a bam infile; user defines
    # bam outfile
    # 
    # :param 1: Name of bam infile, including path (chr)
    # :param 2: Name of bam outfile, including path (chr)
    start="$(date +%s)"

    picard SortSam \
    I="${1}" \
    O="${1%.bam}.queryname.bam" \
    SORT_ORDER="queryname"
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run picard SortSam SORT_ORDER=\"queryname\" on $(basename "${1}")."
}


sort_bam_qname_samtools() {
    # Run samtools sort -n on bam infile; user defines bam outfile
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: name of bam outfile, including path (chr)
    start="$(date +%s)"

    samtools sort -n -@ "${1}" "${2}" > "${3}"

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

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}"

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

    samtools sort -n -@ "${1}" "${2}" > "${2/.bam/.sort-n.bam}"

    samtools fixmate -@ "${1}" "${2/.bam/.sort-n.bam}" "${2/.bam/.fixmate.bam}"

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

    samtools view -@ "${1}" -bh "${2}" "${3}" > "${2/.bam/.${3}.bam}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view to create $(basename "${2/.bam/.${3}.bam}")."
}


sort_file_AS() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    # :param 2: txt.gz outfile, including path
    start="$(date +%s)"
    
    sort -k1,1 -k2n <(gunzip -c "${1}") \
        | gzip \
        > "${1/.txt.gz/.tmp.txt.gz}"
    
    mv -f "${1/.txt.gz/.tmp.txt.gz}" "${2}"
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view on $(basename "${1}")."
}


sort_file_AS_overwrite_infile() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    start="$(date +%s)"

    sort -k1,1 -k2n <(gunzip -c "${1}") \
        | gzip \
        > "${1/.txt.gz/.tmp.txt.gz}"

    mv -f "${1/.txt.gz/.tmp.txt.gz}" "${1}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Run samtools view on $(basename "${1}")."
}


sort_file_qname() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    # :param 2: txt.gz outfile, including path
    sort <(zcat "${1}") | gzip > "${2}"
}


sort_file_qname_auto() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    sort <(zcat "${1}") | gzip > "${1%.txt.gz}.sorted.txt.gz"
}


parsort_file_qname_auto() {
    # Sort a gzipped txt file of QNAMEs (single column)
    #
    # :param 1: txt.gz infile, including path
    # :param 2: remove infile [logical; default: FALSE]
    no_gz="${1%.gz}"
    no_gz_sort="${no_gz%.txt}.srt.txt"
    gz_sort="${no_gz_sort}.gz"

    [[ ! -f "${no_gz}" ]] && pigz -dk "${1}"
    LC_ALL=C parsort "${no_gz}" > "${no_gz_sort}"
    pigz "${no_gz_sort}"
    
    [[ -f "${gz_sort}" ]] && rm "${no_gz}"

    case "$(echo "${2:-"FALSE"}" | tr '[:upper:]' '[:lower:]')" in
        false | f) : ;;
        true | t) rm "${2}" ;;
    esac
}


sort_file_qname_reverse_auto() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    sort -r <(zcat "${1}") | gzip > "${1%.txt.gz}.reverse-sorted.txt.gz"
}


parsort_file_qname_reverse_auto() {
    # Sort a gzipped txt file of QNAMEs (single column)
    #
    # :param 1: txt.gz infile, including path
    # :param 2: remove infile [logical; default: FALSE]
    no_gz="${1%.gz}"
    no_gz_sort="${no_gz%.txt}.rv-srt.txt"
    gz_sort="${no_gz_sort}.gz"

    [[ ! -f "${no_gz}" ]] && pigz -dk "${1}"
    LC_ALL=C parsort -r "${no_gz}" > "${no_gz_sort}"
    pigz "${no_gz_sort}"
    
    [[ -f "${gz_sort}" ]] && rm "${no_gz}"

    case "$(echo "${2:-"FALSE"}" | tr '[:upper:]' '[:lower:]')" in
        false | f) : ;;
        true | t) rm "${2}" ;;
    esac
}


sort_file_duplicate_records_qname_auto() {
    #TODO Write a description
    #
    # :param 1: txt.gz infile, including path
    sort <(zcat "${1}") | sed -e p | gzip > "${1%.txt.gz}.sorted.duplicated.txt.gz"
}

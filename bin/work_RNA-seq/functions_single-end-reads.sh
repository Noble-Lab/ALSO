#!/bin/bash

#  functions_single-end-reads.sh
#  KA


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


list_qnames_to_cut() {
    # Find and list QNAMEs with more than one entry
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"
    message="Find and list QNAMEs with more than one entry for $(basename "${1}")."

    echo "Started: ${message}"
    zcat -df "${1}" \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | sort -nr \
        | cut -c 7- \
        | awk '$1 > 1 {print $0}' \
        | tr ' ' \\t \
        | cut -f 2 \
        | gzip \
        > "${2}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: ${message}."
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


remove_reads_low_quality() {
    # From a bam infile, remove reads with MAPQ < 30; user defines bam outfile
    # name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    # :param 3: name of bam outfile, including path (chr)
    start="$(date +%s)"
    
    samtools view -@ "${1}" -h -b -q 30 "${2}" -o "${3}"

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Filter out reads based on MAPQ scores and pairing status for $(basename "${2}")."
}


remove_reads_low_quality_auto() {
    # From a bam infile, remove reads with MAPQ < 30; outfile name is derived
    # from infile name
    # 
    # :param 1: number of cores for parallelization (int >= 1)
    # :param 2: name of bam infile, including path (chr)
    start="$(date +%s)"
    
    samtools view -@ "${1}" -h -b -q 30 "${2}" -o "${2/.bam/.rm.bam}"

    end="$(date +%s)"
    echo ""
    calculate_run_time "${start}" "${end}"  \
    "Filter out reads based on MAPQ scores and pairing status for $(basename "${2}")."
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


tally_qnames_gzip() {
    # Tally numbers of entries per qname
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"
    message="Tally numbers of entries per QNAME in $(basename "${1}")."

    echo "Started: ${message}"
    zcat -df "${1}" \
        | cut -f 1 \
        | sort \
        | uniq -c \
        | sort -nr \
        | cut -c 7- \
        | tr ' ' \\t \
        | gzip \
        > "${2}"
    
    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: ${message}."
}


tally_qnames_gt_1_gzip() {
    # Tally numbers of entries per qname; retain in list only qnames with more
    # than one entry
    #
    # :param 1: gzipped or uncompressed infile, including path (chr)
    # :param 2: gzipped outfile, including path (chr)
    start="$(date +%s)"
    message="""
    Started: Tally numbers of entries per QNAME in $(basename "${1}"),
    retaining only those with more than one entry in an outlist.
    """
    
    echo "Started: ${message}"
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

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Completed: ${message}"
}


#!/bin/bash

#KEEP
find_set_intersection() {
    # Find, list, and sort the set intersection elements between AS.txt.gz
    # files for samples #1 and #2
    #
    # :param 1: sorted sample #1 AS.txt.gz file
    # :param 2: sorted sample #2 AS.txt.gz file
    # :param 3: for for set intersection elements for samples #1 and #2
    start="$(date +%s)"

    echo "Writing out and sorting set intersection elements for $(basename "${1}") and $(basename "${2}")... "
    join <(zcat -df "${1}") <(zcat -df "${2}") \
        | tr ' ' \\t \
        | sort \
        | gzip \
        > "${3}"

    end="$(date +%s)"
    calculate_run_time "${start}" "${end}" \
    "Write out and sort set intersection elements for $(basename "${1}") and $(basename "${2}")."
}  #FIXME Output should be sorted


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
# for i in *.rm.bam; do echo "${i}" && list_tally_qnames_gzip "${i}"; done
# for i in *.rm.bam; do echo "${i}" && samtools view -c "${i}" && echo ""; done
# for i in *.txt.gz; do echo "${i}" && zcat -df "${i}" | wc -l && echo ""; done


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
# for i in *.txt.gz; do echo "${i}" && parsort_file_qname_auto "${i}" FALSE && echo ""; done
# for i in *.srt.txt.gz; do echo "${i}" &&  sort -c <(zcat -df "${i}") && echo ""; done
# for i in *.srt.txt.gz; do echo "${i}" && zcat -df "${i}" | wc -l && echo ""; done


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

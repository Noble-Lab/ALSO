# retain_qname_reads_picard() {
#     # Filter a bam infile to retain reads with QNAMEs listed in a txt file;
#     # write the filtered results to a bam outfile
#     #
#     # :param 1: name of bam infile, including path (chr)
#     # :param 2: name of txt QNAME list, including path (chr)
#     # :param 3: name of bam outfile, including path (cannot be same as bam
#     #           infile) (chr)
#     # :param 4: initial memory allocation pool for JVM (chr)
#     # :param 5: maximum memory allocation pool for JVM (chr)
#     start="$(date +%s)"
#
#     picard "${4}" "${5}" FilterSamReads \
#     I="${1}" \
#     O="${3}" \
#     READ_LIST_FILE="${2}" \
#     FILTER="includeReadList" \
#     SORT_ORDER="coordinate" \
#     MAX_RECORDS_IN_RAM=500000 \
#     TMP_DIR="${TMPDIR}"
#
#     end="$(date +%s)"
#     calculate_run_time "${start}" "${end}" \
#     "Exclude reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
# }

# retain_qname_reads_picard() {
#     # Filter a bam infile to retain reads with QNAMEs listed in a txt file;
#     # write the filtered results to a bam outfile
#     #
#     # :param 1: name of bam infile, including path (chr)
#     # :param 2: name of txt QNAME list, including path (chr)
#     # :param 3: name of bam outfile, including path (cannot be same as bam
#     #           infile) (chr)
#     # :param 4: initial memory allocation pool for JVM (chr)
#     # :param 5: maximum memory allocation pool for JVM (chr)
#     # :param 6: use the picard.jar available on the GS grid system (logical)
#     start="$(date +%s)"
#     dir_picard="/net/gs/vol3/software/modules-sw/picard/2.26.4/Linux/CentOS7/x86_64"
#
#     case "$(echo "${6}" | tr '[:upper:]' '[:lower:]')" in
#         true | t) \
#             java -jar -Xms"${4}" -Xmx"${5}" \
#             "${dir_picard}"/picard.jar FilterSamReads \
#             I="${1}" \
#             O="${3}" \
#             READ_LIST_FILE="${2}" \
#             FILTER="includeReadList"
#             ;;
#         false | f) \
#             picard -Xms"${4}" -Xmx"${5}" FilterSamReads \
#             I="${1}" \
#             O="${3}" \
#             READ_LIST_FILE="${2}" \
#             FILTER="includeReadList"
#             ;;
#         *) \
#             echo "Exiting: Parameter 6 is not \"TRUE\" or \"FALSE\"."
#             return 1
#             ;;
#     esac
#
#     end="$(date +%s)"
#     calculate_run_time "${start}" "${end}" \
#     "Exclude reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
# }

# exclude_qname_reads_picard() {
#     # Filter a bam infile to exclude reads with QNAMEs listed in a txt file;
#     # write the filtered results to a bam outfile
#     #
#     # :param 1: name of bam infile, including path (chr)
#     # :param 2: name of txt QNAME list, including path (chr)
#     # :param 3: name of bam outfile, including path (cannot be same as bam
#     #           infile) (chr)
#     # :param 4: initial memory allocation pool for JVM (chr)
#     # :param 5: maximum memory allocation pool for JVM (chr)
#     # :param 6: use the picard.jar available on the GS grid system (logical)
#     start="$(date +%s)"
#     dir_picard="/net/gs/vol3/software/modules-sw/picard/2.26.4/Linux/CentOS7/x86_64"
#
#     case "$(echo "${6}" | tr '[:upper:]' '[:lower:]')" in
#         true | t) \
#             java -jar -Xms"${4}" -Xmx"${5}" \
#             "${dir_picard}"/picard.jar FilterSamReads \
#             I="${1}" \
#             O="${3}" \
#             READ_LIST_FILE="${2}" \
#             FILTER="excludeReadList"
#             ;;
#         false | f) \
#             picard -Xms"${4}" -Xmx"${5}" FilterSamReads \
#             I="${1}" \
#             O="${3}" \
#             READ_LIST_FILE="${2}" \
#             FILTER="excludeReadList"
#             ;;
#         *) \
#             echo "Exiting: Parameter 6 is not \"TRUE\" or \"FALSE\"."
#             return 1
#             ;;
#     esac
#
#     end="$(date +%s)"
#     calculate_run_time "${start}" "${end}" \
#     "Exclude reads in $(basename "${1}") based on QNAMEs in $(basename "${2}")."
# }

# picard() {
#     # Call picard using the HPC installation in such a way that the call is
#     # consistent with synatx for calling a conda installation of picard
#     typeset dir_modules="/net/gs/vol3/software/modules-sw"
#     typeset dir_picard="${dir_modules}/picard/2.26.4/Linux/CentOS7/x86_64"
#     typeset jar_picard="${dir_picard}/picard.jar"
#
#     #  Capture arguments
#     arguments=("$@")
#
#     #  Stop function if insufficient numbers of arguments passed
#     if [[ "${#arguments[@]}" -eq 0 ]]; then
#         echo "Stopping: No arguments were passed"
#         return 1
#     fi
#
#     if [[ ! "${#arguments[@]}" -eq 10 ]]; then
#         echo "Stopping: ${#arguments[@]} arguments were passed; ten" \
#         "arguments should be passed"
#         return 1
#     fi
#
#     #  Stop function if passed arguments don't contain necessary substrings
#     if [[ "${arguments[0]}" != "-Xms"* ]]; then
#         echo "Stopping: Minimum heap not specified with '-Xms'"
#         return 1
#     fi
#
#     if [[ "${arguments[1]}" != "-Xmx"* ]]; then
#         echo "Stopping: Maximum heap not specified with '-Xmx'"
#         return 1
#     fi
#
#     if [[ "${arguments[2]}" != "FilterSamReads" ]]; then
#         echo "Stopping: Specified program is not 'FilterSamReads'"
#         return 1
#     fi
#
#     if [[ "${arguments[3]}" != "I="* ]]; then
#         echo "Stopping: Infile not specified with 'I='"
#         return 1
#     fi
#
#     if [[ "${arguments[4]}" != "O="* ]]; then
#         echo "Stopping: Outfile not specified with 'O='"
#         return 1
#     fi
#
#     if [[ "${arguments[5]}" != "READ_LIST_FILE="* ]]; then
#         echo "Stopping: File for filtering not specified with 'READ_LIST_FILE='"
#         return 1
#     fi
#
#     if \
#         [[ "${arguments[6]}" != "FILTER=excludeReadList" ]] && \
#         [[ "${arguments[6]}" != "FILTER=includeReadList" ]]
#     then
#         echo "Stopping: 'FILTER=' parameter must be 'FILTER=excludeReadList'" \
#         "or 'FILTER=includeReadList'"
#         return 1
#     fi
#
#     if [[ "${arguments[7]}" != "SORT_ORDER="* ]]; then
#         echo "Stopping: Argument must be specified with 'SORT_ORDER='"
#         return 1
#     fi
#
#     if [[ "${arguments[8]}" != "MAX_RECORDS_IN_RAM="* ]]; then
#         echo "Stopping: Argument must be specified with 'MAX_RECORDS_IN_RAM='"
#         return 1
#     fi
#
#     if [[ "${arguments[9]}" != "TMP_DIR="* ]]; then
#         echo "Stopping: Argument must be specified with 'TMP_DIR='"
#         return 1
#     fi
#
#     min="${arguments[0]}"
#     max="${arguments[1]}"
#     program="${arguments[2]}"
#     input="${arguments[3]}"
#     output="${arguments[4]}"
#     read_list_file="${arguments[5]}"
#     filter="${arguments[6]}"
#     sort_order="${arguments[7]}"
#     max_records_in_ram="${arguments[8]}"
#     tmp_dir="${arguments[9]}"
#
#     # echo "Running..." \
#     # "java ${min} ${max} -jar ${jar_picard}" \
#     # "${program} ${input} ${output} ${read_list_file} ${filter}"
#
#     java "${min}" "${max}" -jar "${jar_picard}" \
#     "${program}" "${input}" "${output}" "${read_list_file}" "${filter}" \
#     "${sort_order}" "${max_records_in_ram}" "${tmp_dir}"
# }

picard() {
    # Call picard using the HPC installation in such a way that the call is
    # consistent with synatx for calling a conda installation of picard
    typeset dir_modules="/net/gs/vol3/software/modules-sw"
    typeset dir_picard="${dir_modules}/picard/2.26.4/Linux/CentOS7/x86_64"
    typeset jar_picard="${dir_picard}/picard.jar"

    #  Capture arguments
    arguments=("$@")

    #  Stop function if insufficient numbers of arguments passed
    if [[ "${#arguments[@]}" -eq 0 ]]; then
        echo "Stopping: No arguments were passed"
        return 1
    fi

    if [[ ! "${#arguments[@]}" -eq 7 ]]; then
        echo "Stopping: ${#arguments[@]} arguments were passed; seven" \
        "arguments should be passed"
        return 1
    fi

    #  Stop function if passed arguments don't contain necessary substrings
    if [[ "${arguments[0]}" != "-Xms"* ]]; then
        echo "Stopping: Minimum heap not specified with '-Xms'"
        return 1
    fi

    if [[ "${arguments[1]}" != "-Xmx"* ]]; then
        echo "Stopping: Maximum heap not specified with '-Xmx'"
        return 1
    fi

    if [[ "${arguments[2]}" != "FilterSamReads" ]]; then
        echo "Stopping: Specified program is not 'FilterSamReads'"
        return 1
    fi

    if [[ "${arguments[3]}" != "I="* ]]; then
        echo "Stopping: Infile not specified with 'I='"
        return 1
    fi

    if [[ "${arguments[4]}" != "O="* ]]; then
        echo "Stopping: Outfile not specified with 'O='"
        return 1
    fi

    if [[ "${arguments[5]}" != "READ_LIST_FILE="* ]]; then
        echo "Stopping: File for filtering not specified with 'READ_LIST_FILE='"
        return 1
    fi

    if \
        [[ "${arguments[6]}" != "FILTER=excludeReadList" ]] && \
        [[ "${arguments[6]}" != "FILTER=includeReadList" ]]
    then
        echo "Stopping: 'FILTER=' parameter must be 'FILTER=excludeReadList'" \
        "or 'FILTER=includeReadList'"
        return 1
    fi

    min="${arguments[0]}"
    max="${arguments[1]}"
    program="${arguments[2]}"
    input="${arguments[3]}"
    output="${arguments[4]}"
    read_list_file="${arguments[5]}"
    filter="${arguments[6]}"

    # echo "Running..." \
    # "java ${min} ${max} -jar ${jar_picard}" \
    # "${program} ${input} ${output} ${read_list_file} ${filter}"

    java "${min}" "${max}" -jar "${jar_picard}" \
    "${program}" "${input}" "${output}" "${read_list_file}" "${filter}"
}
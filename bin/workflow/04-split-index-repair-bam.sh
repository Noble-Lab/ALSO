#!/bin/bash

#  split-index-repair-bam.sh
#  KA


start="$(date +%s)"


#  Functions ------------------------------------------------------------------
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


evaluateGreaterEqual() {
    # Sort two values such that the greater is placed ahead of the lesser
    # 
    # :param 1: first value
    # :param 2: second value
    printf '%s\n%s\n' "${2}" "${1}" | sort -V -C
}


# loopEcho() { for i in "${@:-*}"; do echo "${i}"; done; }


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Take a bam infile containing all mouse chromosomes, then output the"
    echo "following for one or all mouse chromosomes:"
    echo "  - bam file (both \"mm10 mode\" and \"strain mode\")"
    echo "  - bam index (both \"mm10 mode\" and \"strain mode\")"
    echo "  - \"POS\" bed file for RNAME, POS, POS + 49, QNAME (only \"strain"
    echo "    mode\")"
    echo "  - \"MPOS\" bed file for MRNM, MPOS, MPOS + 49, QNAME (only"
    echo "    \"strain mode\")"
    echo ""
    echo "Chromosomes in the bam infile are assumed to be in \"chrN\" format."
    echo "The bam infile is also assumed to be coordinate sorted; if not, then"
    echo "the script will detect this, copy the bam infile into a bak"
    echo "directory, and then proceed to coordinate sort the infile."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - bedtools >= 2.29.0"
    echo "  - parallel >= 20200101"
    echo "  - repair >= 2.0.1"
    echo "  - samtools >= 1.13"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
    echo "-i <bam infile, including path (chr)>"
    echo "-o <path for outfile(s; chr); path will be made if it does not exist>"
    echo "-x <prefix for outfile(s; chr)>"
    echo "-c <chromosome(s) to split out (chr); for example, \"chr1\" for"
    echo "    chromosome 1, \"chrX\" for chromosome X, \"all\" for all"
    echo "    chromosomes>"
    echo "-m <mode in which to run the script: \"M\" or \"S\" (chr); with"
    echo "    \"M\" (or \"mm10\"), singletons will be removed and Subread"
    echo "    repair will be run on split bam files, but \"POS\" and \"MPOS\""
    echo "    bed files will not be generated (since liftOver coordinate"
    echo "    conversion to mm10 will not need to be performed); with \"S\""
    echo "    (or \"strain\"), singletons will be removed, Subread repair will"
    echo "    be run, and \"POS\" and \"MPOS\" bed files will be generated (to"
    echo "    be used in subsequent liftOver coordinate conversion)>"
    # echo "-r <use Subread repair on split bam files: \"TRUE\" or \"FALSE\"" 
    # echo "    (logical); default: \"TRUE\" if \"mm10 mode\" is \"FALSE\">"
    # echo "-b <if \"-r TRUE\", create bed files from split bam files: \"TRUE\""
    # echo "    or \"FALSE\" (logical); argument \"-b\" only needed when \"-r"
    # echo "    TRUE\"; default: \"TRUE\" if \"mm10 mode\" is \"FALSE\">"
    echo "-p <number of cores for parallelization (int >= 1); default: 1>"
    exit
}

while getopts "h:u:i:o:x:c:m:r:b:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        x) prefix="${OPTARG}" ;;
        c) chromosome="${OPTARG}" ;;
        m) mode="${OPTARG}" ;;
        r) repair="${OPTARG}" ;;
        b) bed="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${prefix}" ]] && printUsage
[[ -z "${chromosome}" ]] && printUsage
[[ -z "${mode}" ]] && printUsage
[[ -z "${repair}" ]] && repair="TRUE"
[[ -z "${bed}" ]] && bed="TRUE"
[[ -z "${parallelize}" ]] && parallelize=1


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Check for necessary dependencies; exit if not found
checkDependency bedtools
checkDependency parallel
checkDependency repair
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

#  Evaluate "mode" argument, which leads to evaluating the "Subread repair"
#+ and "bed" arguments
case "$(echo "${mode}" | tr '[:upper:]' '[:lower:]')" in
    mm10 | m) \
        echo -e "-m: Will run the script in \"mm10 mode.\"\n"
        flag_bed=0
        flag_subread=1
        ;;
    strain | s) \
        echo -e "-m: Will run the script in \"strain mode.\"\n"
        flag_bed=1
        flag_subread=1
        ;;
    test | t) \
        echo -e "-m: Will run the script in undocumented \"test mode,\" which"
        echo -e "    parses input to the undocumented \"-r\" and \"-b\""
        echo -e "    arguments.\n"
        
        #  Set flag to repair bam file(s) if "${repair}" is TRUE
        case "$(echo "${repair}" | tr '[:upper:]' '[:lower:]')" in
            true | t) \
                echo -e "-r: Will use Subread repair on split bam files."
                flag_subread=1
                ;;
            false | f) \
                echo -e "-r: Will not use Subread repair on split bam files."
                flag_subread=0
                ;;
            *) \
                echo -e "Exiting: -r Subread repair argument must be \"TRUE\" or \"FALSE\".\n"
                exit 1
                ;;
        esac
        
        #  Set flag to create bed file(s) if "${bed}" is TRUE
        case "$(echo "${bed}" | tr '[:upper:]' '[:lower:]')" in
            true | t) \
                if [[ $((flag_subread)) -eq 1 ]]; then
                    echo -e "-b: Will create bed file(s) from split bam file(s).\n"
                    echo -e ""
                    flag_bed=1
                elif [[ $((flag_subread)) -eq 0 ]]; then
                    echo -e "-b: Not using Subread repair on split bam files; thus,"
                    echo -e "    cannot create bed file(s) from split bam file(s).\n"
                    echo -e ""
                    flag_bed=0
                else
                    echo -e "Exiting: There was an error processing the \"-r\" and \"-b\" arguments.\n"
                    exit 1
                fi
                ;;
            false | f) \
                if [[ $((flag_subread)) -eq 1 || $((flag_subread)) -eq 0 ]]; then
                    echo -e "-b: Will not create bed file(s) from split bam file(s).\n"
                    echo -e ""
                    flag_bed=0
                else
                    echo -e "Exiting: There was an error processing the \"-r\" and \"-b\" arguments.\n"
                    exit 1
                fi
                ;;
            *) \
                echo -e "Exiting: -b bed argument must be \"TRUE\", \"FALSE\", or undefined.\n"
                echo -e ""
                exit 1
                ;;
        esac
        ;;
    *) \
        echo -e "Exiting: -m \"mode\" argument must be \"M\" or \"S\".\n"
        exit 1
        ;;
esac

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
#  If bam infile is not sorted and/or indexed, then do so; in order to split a
#+ bam file by chromosome, the file must be sorted and indexed
if [[ $(samtools index -@ "${parallelize}" "${infile}")$? -eq 0 ]]; then
    echo -e "#  Step 0"
    echo -e "Bam file is sorted."
else
    echo -e "#  Step 0"
    echo -e "Bam file is not sorted and index is presumed absent."
    echo -e "Started: Sorting and indexing bam infile."

    echo -e "Copy infile to back-up directory prior to sorting, etc."
    mkdir -p "${infile%/*}/bak"
    cp "${infile}" "${infile%/*}/bak" &
    displaySpinningIcon $! "Copying infile to back-up directory... "

    samtools sort -@ "${parallelize}" "${infile}" -o "${infile}" &
    displaySpinningIcon $! "Sorting bam infile... "

    samtools index -@ "${parallelize}" "${infile}" &
    displaySpinningIcon $! "Indexing bam infile... "

    echo -e "Completed: Sorting and indexing bam infile.\n"
fi

#  Split bam file by chromosome, then perform additional processing
case "${chromosome}" in
    all) \
        unset all
        unset name
        unset tmp
        typeset -a all
        typeset -a name
        typeset -a tmp
        for i in $(seq 1 19) "X" "Y"; do all+=( "chr${i}" ); done
        for i in $(seq 1 19) "X" "Y"; do name+=( "${outpath}/${prefix}.chr${i}.bam" ); done
        for i in $(seq 1 19) "X" "Y"; do tmp+=( "${outpath}/${prefix}.chr${i}.bam.tmp" ); done

        #  Step 1: Split bam infile -----------------------
        echo -e "#  Step 1"
        echo -e "Started: Splitting bam infile into individual bam files, one for each chromosome."

        parallel -k -j "${parallelize}" \
        "samtools view -b {1} {2} > {3}" \
        ::: "${infile}" \
        ::: "${all[@]}" \
        :::+ "${name[@]}"

        echo -e "Completed: Splitting bam infile into individual bam files, one for each chromosome.\n"


        #  Step 2: Remove singletons from split bam file --
        echo -e "#  Step 2"
        echo -e "Started: Removing singletons from each split bam file."

        parallel -k -j "${parallelize}" \
        "samtools view -b -F 0x08 {1} > {2} && mv -f {2} {1}" \
        ::: "${name[@]}" \
        :::+ "${tmp[@]}"

        echo -e "Completed: Removing singletons from bam file for chromosome ${chromosome}.\n"


        #  Step 3: Index each split bam file --------------
        echo -e "#  Step 3"
        echo -e "Started: Indexing each split bam file."

        echo -e "Indexing each split bam file... "
        parallel -k -j "${parallelize}" \
        "samtools index {1}" \
        ::: "${name[@]}"
       
        echo -e "Completed: Indexing each split bam file.\n"


        #  Step 4: Repair each split bam file -------------
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "#  Step 4"
            echo -e "Started: Running Subread repair on each split bam file."
           
            parallel -k -j 1 \
            "repair -d -T {1} -c -i {2} -o {3} && mv -f {3} {2}" \
            ::: "${parallelize}" \
            ::: "${name[@]}" \
            :::+ "${tmp[@]}"

            echo -e "Completed: Running Subread repair on each split bam file.\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            echo -e "#  Step 4"
            echo -e "Skipping: Running Subread repair on each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread repair flag.\n"
            exit 1
        fi


        #  Step 5: Create "POS" bed files -------------------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 5"
            echo -e "Started: Creating bed files from each split bam file."

            echo -e "Creating temporary bedpe files... "
            parallel -k -j "${parallelize}" \
            "bedtools bamtobed -i {1} -bedpe > {2}" \
            ::: "${name[@]}" \
            :::+ "${name[@]/.bam/.bedpe}"

            # echo -e "Adjusting bedpe \"POS\", \"MPOS\", etc. values... "
            echo -e "Adjusting bedpe \"POS\", etc. values... "
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2 + 1, \$3 + 1, \$4, \$5 + 1, \$6 + 1, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.bedpe.tmp}"

            echo -e "Cleaning up bedpe.tmp files... "
            parallel -k -j "${parallelize}" \
            "mv -f {1} {2}" \
            ::: "${name[@]/.bam/.bedpe.tmp}" \
            :::+ "${name[@]/.bam/.bedpe}"

            echo -e "Creating \"POS\" bed files... "
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2, \$3, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.pos.bed}"

            # echo -e "Creating \"MPOS\" bed files... "
            # parallel -k -j "${parallelize}" \
            # "awk 'BEGIN{FS=OFS=\"\t\"} {print \$4, \$5, \$6, \$7}' {1} > {2}" \
            # ::: "${name[@]/.bam/.bedpe}" \
            # :::+ "${name[@]/.bam/.mpos.bed}"

            echo -e "Removing temporary bedpe files... "
            parallel -k -j "${parallelize}" \
            "rm {1}" \
            ::: "${name[@]/.bam/.bedpe}"

            echo -e "Completed: Creating bed files from each split bam file.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
            echo -e "#  Step 5"
            echo -e "Skipping: Creating bed files from each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the bed flag.\n"
            exit 1
        fi
        ;;
    chr1 | chr2 | chr3 | chr4 | chr5 | chr6 | chr7 | chr8 | chr9 | chr10 | \
    chr11 | chr12 | chr13 | chr14 | chr15 | chr16 | chr17 | chr18 | chr19 | \
    chrX | chrY) \
        #  Step 1: Split bam infile -----------------------
        echo -e "#  Step 1"
        echo -e "Started: Splitting bam infile into individual bam file for chromosome ${chromosome}."

        samtools view -@ "${parallelize}" -b "${infile}" "${chromosome}" \
        > "${outpath}/${prefix}.${chromosome}.bam" &
        displaySpinningIcon $! "Splitting bam infile... "

        echo -e "Completed: Splitting bam infile into individual bam file for chromosome ${chromosome}.\n"


        #  Step 2: Remove singletons from split bam file --
        echo -e "#  Step 2"
        echo -e "Started: Removing singletons from bam file for chromosome ${chromosome}."

        samtools view -@ "${parallelize}" -b -F 0x08 \
        "${outpath}/${prefix}.${chromosome}.bam" \
        > "${outpath}/${prefix}.${chromosome}.bam.tmp" &
        displaySpinningIcon $! "Removing singletons from split bam infile... "

        mv -f \
        "${outpath}/${prefix}.${chromosome}.bam.tmp" \
        "${outpath}/${prefix}.${chromosome}.bam"

        echo -e "Completed: Removing singletons from bam file for chromosome ${chromosome}.\n"


        #  Step 3: Index each split bam file --------------
        echo -e "#  Step 3"
        echo -e "Started: Indexing bam file for chromosome ${chromosome}."

        samtools index -@ "${parallelize}" "${outpath}/${prefix}.${chromosome}.bam" &
        displaySpinningIcon $! "Indexing split bam file... "
               
        echo -e "Completed: Indexing bam for chromosome ${chromosome}.\n"


        #  Step 4: Repair the split bam file --------------
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "#  Step 4"
            echo -e "Started: Running Subread repair on the split bam file."
            
            repair \
            -d -c -T "${parallelize}" \
            -i "${outpath}/${prefix}.${chromosome}.bam" \
            -o "${outpath}/${prefix}.${chromosome}.bam.tmp" &
            displaySpinningIcon $! "Running repair on split bam file... "

            mv -f \
            "${outpath}/${prefix}.${chromosome}.bam.tmp" \
            "${outpath}/${prefix}.${chromosome}.bam" &
            displaySpinningIcon $! "Removing unneeded intermediate file... "

            echo -e "Completed: Running Subread repair on the split bam file.\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            echo -e "#  Step 4"
            echo -e "Skipping: Running Subread repair on each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread repair flag.\n"
            exit 1
        fi


        #  Step 5: Create "POS" bed files -----------------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 5"
            echo -e "Started: Creating bed files for chromosome \"${chromosome}.\""

            bedtools bamtobed -i "${outpath}/${prefix}.${chromosome}.bam" -bedpe \
            > "${outpath}/${prefix}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Creating temporary bedpe file... "

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2 + 1, $3 + 1, $4, $5 + 1, $6 + 1, $7}' \
            "${outpath}/${prefix}.${chromosome}.bedpe" \
            > "${outpath}/${prefix}.${chromosome}.bedpe.tmp" &
            displaySpinningIcon $! "Adjusting bedpe \"POS\" etc. values... "
            # displaySpinningIcon $! "Adjusting bedpe \"POS\", \"MPOS\", etc. values... "

            mv -f \
            "${outpath}/${prefix}.${chromosome}.bedpe.tmp" \
            "${outpath}/${prefix}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Cleaning up bedpe.tmp file... "

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $7}' \
            "${outpath}/${prefix}.${chromosome}.bedpe" \
            > "${outpath}/${prefix}.${chromosome}.pos.bed" &
            displaySpinningIcon $! "Creating \"POS\" bed file... "

            # awk 'BEGIN{FS=OFS="\t"} {print $4, $5, $6, $7}' \
            # "${outpath}/${prefix}.${chromosome}.bedpe" \
            # > "${outpath}/${prefix}.${chromosome}.mpos.bed" &
            # displaySpinningIcon $! "Creating \"MPOS\" bed file... "

            rm "${outpath}/${prefix}.${chromosome}.bedpe"  &
            displaySpinningIcon $! "Removing temporary bedpe file... "

            echo -e ""
            echo -e "Completed: Creating bed files for chromosome ${chromosome}.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
            echo -e "#  Step 5"
            echo -e "Skipping: Creating bed files for each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the bed flag.\n"
            exit 1
        fi
        ;;
    *) \
        echo -e "Exiting: Script cannot handle \"${chromosome}\" as chromosome option.\n"
        exit 1
        ;;
esac


#  End recording time ---------------------------------------------------------
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "Run time: ${run_time} seconds."
echo ""

exit 0

#!/bin/bash

#  split-index-repair-bam.sh
#  KA


#  Start recording time -------------------------------------------------------
start="$(date +%s)"


#  Functions ------------------------------------------------------------------
checkDependency() {
    #  Check if program is available in "${PATH}"; exit if not
    command -v "${1}" &>/dev/null ||
        {
            echoErrOut "Exiting: ${1} not found. Install ${1}."
            # exit 1
        }
}

displaySpinningIcon() {
    #  Display "spinning icon" while a background process runs
    spin="/|\\–/|\\-"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


# loopEcho() { for i in "${@:-*}"; do echo "${i}"; done; }


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Take a bam infile containing all mouse chromosomes, then output the"
    echo "following for one or all mouse chromosomes:"
    echo " - bam file"
    echo " - bam index"
    echo " - \"pos\" bed file for RNAME, POS, POS + 49, QNAME"
    echo " - \"mpos\" bed file for MRNM, MPOS, MPOS + 49, QNAME"
    echo ""
    echo "Chromosomes in bam infile are assumed to be in \"chrN\" format."
    echo ""
    echo ""
    echo "Dependencies:"
    echo " - bedtools >= 2.30.0"
    echo " - parallel >= 20200101"
    echo " - repair >= 2.0.1"
    echo " - samtools >= 1.13"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
    echo "-i <bam infile, including path (chr)>"
    echo "-o <path for split bam file(s) (chr); path will be made if it does"
    echo "    not exist>"
    echo "-c <chromosome(s) to split out (chr); for example, \"chr1\" for"
    echo "    chromosome 1, \"chrX\" for chromosome X, \"all\" for all"
    echo "    chromosomes>"
    echo "-r <use Subread repair on split bam files: \"TRUE\" or \"FALSE\"" 
    echo "    (logical)>"
    echo "-b <if \"-r TRUE\", create bed files from split bam files: \"TRUE\""
    echo "    or \"FALSE\" (logical); argument \"-b\" only needed when \"-r"
    echo "    TRUE\">"
    echo "-p <number of cores for parallelization (int >= 1)>"
    exit
}

while getopts "h:u:i:o:c:r:b:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        c) chromosome="${OPTARG}" ;;
        r) repair="${OPTARG}" ;;
        b) bed="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${chromosome}" ]] && printUsage
[[ -z "${repair}" ]] && printUsage
[[ -z "${bed}" ]] && bed="FALSE"
[[ -z "${parallelize}" ]] && printUsage

# #  Test defaults
# safe_mode="FALSE"
# infile="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/Disteche_sample_1.dedup.bam"
# outpath="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/2022-0316_outpath_test"
# chromosome="chr19"
# repair="TRUE"
# bed="TRUE"
# parallelize=4


#  Check variable assignments -------------------------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Check for necessary dependencies; exit if not found
checkDependency bedtools
checkDependency parallel
checkDependency repair
checkDependency samtools

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
        echo -e "Exiting: -p parallelize argument must be an integer."
        exit 1
    }

[[ ! $((parallelize)) -ge 1 ]] &&
    {
        echo -e "Exiting: -p parallelize argument must be an integer >= 1."
        exit 1
    }

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
        fi
        ;;
    false | f) \
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "-b: Will not create bed file(s) from split bam file(s).\n"
            echo -e ""
            flag_bed=0
        elif [[ $((flag_subread)) -eq 0 ]]; then
            echo -e "-b: Will not create bed file(s) from split bam file(s).\n"
            echo -e ""
            flag_bed=0
        fi
        ;;
    *) \
        echo -e "Exiting: -b bed argument must be \"TRUE\", \"FALSE\", or undefined.\n"
        echo -e ""
        exit 1
        ;;
esac


#  Process bam infile ---------------------------------------------------------
#  If bam infile is not indexed, then do so
if [[ ! -f "${infile/.bam/.bam.bai}" ]]; then
    echo -e "#  Step 0"
    echo -e "Infile bam index not found."
    echo -e "Started: Indexing bam infile."

    samtools index -@ "${parallelize}" "${infile}" &
    displaySpinningIcon $! "Indexing bam infile... "

    echo -e "Completed: Indexing bam infile.\n"
else
    echo -e "#  Step 0"
    echo -e "Index for bam infile found. Moving on to next step.\n"
    :
fi

#  Run samtools view for "${chromosome}"
case "${chromosome}" in
    all) \
        unset all
        unset name
        unset tmp
        typeset -a all
        typeset -a name
        typeset -a tmp
        for i in $(seq 1 19) "X" "Y"; do all+=( "chr${i}" ); done
        for i in $(seq 1 19) "X" "Y"; do name+=( "${outpath}/$(basename "${infile%.bam}").chr${i}.bam" ); done
        for i in $(seq 1 19) "X" "Y"; do tmp+=( "${outpath}/$(basename "${infile%.bam}").chr${i}.bam.tmp" ); done

        #  Step 1: Split bam infile -----------------------
        echo -e "#  Step 1"
        echo -e "Started: Splitting bam infile into individual bam files, one for each chromosome."

        #TODO Skip this step if split bam files exist
        parallel -k -j "${parallelize}" \
        "samtools view -b {1} {2} > {3}" \
        ::: "${infile}" \
        ::: "${all[@]}" \
        :::+ "${name[@]}"

        echo -e "Completed: Splitting bam infile into individual bam files, one for each chromosome.\n"


        #  Step 2: Index each split bam file --------------
        echo -e "#  Step 2"
        echo -e "Started: Indexing each split bam file."

        #TODO Skip this step if split bam indices exist
        parallel -k -j "${parallelize}" \
        "samtools index {1}" \
        :::+ "${name[@]}"
        
        echo -e "Completed: Indexing each split bam file.\n"


        #  Step 3: Repair each split bam file -------------
        #TODO Skip this step if split bam files are already "repaired"
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "#  Step 3"
            echo -e "Started: Running Subread repair on each split bam file."
            
            parallel -k -j 1 \
            "repair -d -T {1} -c -i {2} -o {3} && mv -f {3} {2}" \
            ::: "${parallelize}" \
            ::: "${name[@]}" \
            :::+ "${tmp[@]}"

            echo -e "Completed: Running Subread repair on each split bam file.\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            echo -e "#  Step 3"
            echo -e "Skipping: Running Subread repair on each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread repair flag.\n"
            exit 1
        fi


        #  Step 4: Create "pos", "mpos" bed files ---------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 4"
            echo -e "Started: Creating bed files from each split bam file."

            echo -e "Creating temporary bedpe files... "
            parallel -k -j "${parallelize}" \
            "bamToBed -i {1} -bedpe > {2}" \
            ::: "${name[@]}" \
            :::+ "${name[@]/.bam/.bedpe}"

            echo -e "Adjusting bedpe \"pos\", \"mpos\", etc. values... "
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2 + 1, \$3 + 1, \$4, \$5 + 1, \$6 + 1, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.bedpe.tmp}"

            echo -e "Cleaning up bedpe.tmp files... "
            parallel -k -j "${parallelize}" \
            "mv -f {1} {2}" \
            ::: "${name[@]/.bam/.bedpe.tmp}" \
            :::+ "${name[@]/.bam/.bedpe}"

            echo -e "Creating \"pos\" bed files... "
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2, \$3, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.pos.bed}"

            echo -e "Creating \"mpos\" bed files... "
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$4, \$5, \$6, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.mpos.bed}"

            echo -e "Removing temporary bedpe files... "
            parallel -k -j "${parallelize}" \
            "rm {1}" \
            ::: "${name[@]/.bam/.bedpe}"

            echo -e "Completed: Creating bed files from each split bam file.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
            echo -e "#  Step 4"
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
        infile_base="$(basename "${infile}" .bam)"

        #  Step 1: Split bam infile -----------------------
        echo -e "#  Step 1"
        echo -e "Started: Splitting bam infile into individual bam file for chromosome ${chromosome}."

        samtools view -@ "${parallelize}" -b "${infile}" "${chromosome}" \
        > "${outpath}/${infile_base}.${chromosome}.bam" &
        displaySpinningIcon $! "Splitting bam infile... "

        echo -e "Completed: Splitting bam infile into individual bam file for chromosome ${chromosome}.\n"


        #  Step 2: Index the split bam file ---------------
        echo -e "#  Step 2"
        echo -e "Started: Indexing bam file for chromosome ${chromosome}."

        samtools index -@ "${parallelize}" "${outpath}/${infile_base}.${chromosome}.bam" &
        displaySpinningIcon $! "Indexing split bam file... "
                
        echo -e "Completed: Indexing bam file for chromosome ${chromosome}.\n"


        #  Step 3: Repair the split bam file --------------
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "#  Step 3"
            echo -e "Started: Running Subread repair on the split bam file."
            
            repair \
            -d -c -T "${parallelize}" \
            -i "${outpath}/${infile_base}.${chromosome}.bam" \
            -o "${outpath}/${infile_base}.${chromosome}.bam.tmp" &
            displaySpinningIcon $! "Running repair on split bam file... "

            [[ -f "${outpath}/${infile_base}.${chromosome}.bam.tmp" ]] &&
                {
                    mv -f \
                    "${outpath}/${infile_base}.${chromosome}.bam.tmp" \
                    "${outpath}/${infile_base}.${chromosome}.bam"
                }

            echo -e "Completed: Running Subread repair on the split bam file.\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            echo -e "#  Step 3"
            echo -e "Skipping: Running Subread repair on each split bam file.\n"
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread repair flag.\n"
            exit 1
        fi


        #  Step 4: Create "pos", "mpos" bed files ---------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 4"
            echo -e "Started: Creating bed files for chromosome \"${chromosome}.\""

            bamToBed -i "${outpath}/${infile_base}.${chromosome}.bam" -bedpe \
            > "${outpath}/${infile_base}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Creating temporary bedpe file... "

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2 + 1, $3 + 1, $4, $5 + 1, $6 + 1, $7}' \
            "${outpath}/${infile_base}.${chromosome}.bedpe" \
            > "${outpath}/${infile_base}.${chromosome}.bedpe.tmp" &
            displaySpinningIcon $! "Adjusting bedpe \"pos\", etc. values... "

            mv -f \
            "${outpath}/${infile_base}.${chromosome}.bedpe.tmp" \
            "${outpath}/${infile_base}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Cleaning up bedpe.tmp file... "

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $7}' \
            "${outpath}/${infile_base}.${chromosome}.bedpe" \
            > "${outpath}/${infile_base}.${chromosome}.pos.bed" &
            displaySpinningIcon $! "Creating \"pos\" bed file... "

            awk 'BEGIN{FS=OFS="\t"} {print $4, $5, $6, $7}' \
            "${outpath}/${infile_base}.${chromosome}.bedpe" \
            > "${outpath}/${infile_base}.${chromosome}.mpos.bed" &
            displaySpinningIcon $! "Creating \"mpos\" bed file... "

            rm "${outpath}/${infile_base}.${chromosome}.bedpe"  &
            displaySpinningIcon $! "Removing temporary bedpe file... "

            echo -e ""
            echo -e "Completed: Creating bed files for chromosome ${chromosome}.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
            echo -e "#  Step 4"
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
echo "${0} run time: ${run_time} seconds."
echo ""

exit 0


#  Scraps ---------------------------------------------------------------------

# [[ ${#name[@]} -eq 21 ]] &&
#     {
#         echo -e "Exiting: Split bam files already exist."
#     }

# [[ ! -f "${infile%.bam}.${chromosome}.bam" ]] &&
#     {
#         echo "..."
#     }

# 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | \
# 17 | 18 | 19 | X | Y) \
#     parallel -k -j 1 \
#     "samtools view -b {1} {2} > {3}" \
#     ::: "${infile}" \
#     ::: "chr${chromosome}" \
#     ::: "${infile%.bam}.chr${chromosome}.bam"
#     ;;

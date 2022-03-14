#!/bin/bash

#  split-index-repair-bam.sh
#  KA


#  Functions ------------------------------------------------------------------
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


loopEcho() { for i in "${@:-*}"; do echo "${i}"; done; }


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Take in a bam file containing all mouse chromosomes, output a bam"
    echo "file for one chromosome or bam files for all mouse chromosomes."
    echo "Chromosomes in bam file are assumed to be in \"chr{#,X,Y}\" format."
    echo ""
    echo "Dependencies:"
    echo " - bedtools >= 2.30.0"
    echo " - liftover >= 366"
    echo " - parallel >= 20200101"
    echo " - repair >= 2.0.1"
    echo " - samtools >= 1.13"
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: TRUE or FALSE (logical)>"
    echo "-i <bam infile, including path (chr)>"
    echo "-r <use Subread repair on bam infile: TRUE or FALSE (logical)>"
    echo "-b <create bed files from split bam files: TRUE or FALSE (logical)>"
    echo "-c <chromosome(s) to split out (chr); for example, \"chr1\" for"
    echo "    chromosome 1, \"chrX\"for chromosome X, \"all\" for all"
    echo "    chromosomes>"
    echo "-p <number of cores for parallelization (int >= 1)>"
    exit
}

while getopts "h:u:i:r:b:c:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        r) repair="${OPTARG}" ;;
        b) bed="${OPTARG}" ;;
        c) chromosome="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${repair}" ]] && printUsage
[[ -z "${bed}" ]] && printUsage
[[ -z "${chromosome}" ]] && printUsage
[[ -z "${parallelize}" ]] && printUsage

# #  Test defaults
# [[ -z "${safe_mode}" ]]   && safe_mode="FALSE"
# [[ -z "${infile}" ]]      && infile="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/Disteche_sample_1.dedup.bam"
# [[ -z "${repair}" ]]      && repair="TRUE"
# [[ -z "${bed}" ]]         && bed="TRUE"
# [[ -z "${chromosome}" ]]  && chromosome="chrX"
# [[ -z "${parallelize}" ]] && parallelize=4
#
# safe_mode="FALSE"
# infile="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/Disteche_sample_1.dedup.bam"
# repair="TRUE"
# bed="TRUE"
# chromosome="chrX"
# parallelize=4


#  Check variable assignments -------------------------------------------------
#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-s: Safe mode is on.\n"
        set -Eeuxo pipefail
        ;;
    false | f) \
        echo -e "-s: safe mode is off.\n"
        :
        ;;
    *) \
        echo -e "Exiting: -s safe-mode argument must be TRUE or FALSE.\n"
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

#  Set flag to repair bam file(s) if "${repair}" is TRUE
case "$(echo "${repair}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-r: will use Subread repair on the bam infile.\n"
        flag_subread=1
        ;;
    false | f) \
        echo -e "-r: will not use Subread repair on the bam infile.\n"
        flag_subread=0
        ;;
    *) \
        echo -e "Exiting: -r Subread repair argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Set flag to create bed file(s) if "${bed}" is TRUE
case "$(echo "${bed}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-b: will create bed file(s) from split bam file(s).\n"
        flag_bed=1
        ;;
    false | f) \
        echo -e "-b: will not create bed file(s) from split bam file(s).\n"
        flag_bed=0
        ;;
    *) \
        echo -e "Exiting: -b bed argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac


#  Process bam infile ---------------------------------------------------------
#  If bam infile is not indexed, then do so
if [[ ! -f "${infile/.bam/.bam.bai}" ]]; then
    echo -e "Infile bam index not found."
    echo -e "Started: Indexing bam infile."

    samtools index -@ "${parallelize}" "${infile}" &
    displaySpinningIcon $! "Indexing bam infile..."

    echo -e "Completed: Indexing bam infile.\n"
else
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
        for i in $(seq 1 19) "X" "Y"; do name+=( "${infile%.bam}.chr${i}.bam" ); done
        for i in $(seq 1 19) "X" "Y"; do tmp+=( "${infile%.bam}.chr${i}.bam.tmp" ); done

        #  Step 1: Split bam infile -----------------------
        echo -e "Started: Splitting bam infile into individual bam files, one"
        echo -e "         for each chromosome."

        #TODO Skip this step if split bam files exist
        parallel -k -j "${parallelize}" \
        "samtools view -b {1} {2} > {3}" \
        ::: "${infile}" \
        ::: "${all[@]}" \
        :::+ "${name[@]}"

        echo -e "Completed: Splitting bam infile into individual bam files, one"
        echo -e "           for each chromosome.\n"


        #  Step 2: Index each split bam file --------------
        echo -e "Started: Indexing each split bam file."

        #TODO Skip this step if split bam indices exist
        parallel -k -j "${parallelize}" \
        "samtools index {1}" \
        ::: "${name[@]}"
        
        echo -e "Completed: Indexing each split bam file.\n"


        #  Step 3: Repair each split bam file -------------
        
        #TODO Skip this step if split bam files are already "repaired"
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "Started: Running Subread repair on each individual bam"
            echo -e "         file."
            
            parallel -k -j 1 \
            "repair -d -T {1} -c -i {2} -o {3} && mv -f {3} {2}" \
            ::: "${parallelize}" \
            ::: "${name[@]}" \
            :::+ "${tmp[@]}"

            echo -e "Completed: Running Subread repair on each individual bam"
            echo -e "           file.\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread"
            echo -e "         repair flag.\n"
            exit 1
        fi


        #  Step 4: Create "pos", "mpos" bed files ---------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 4 ----------------------------------------"
            echo -e "Started: Creating bed files for chromosome"
            echo -e "         ${chromosome}."

            echo -e "Creating temporary bedpe files..."
            parallel -k -j "${parallelize}" \
            "bamToBed -i {1} -bedpe > {2}" \
            ::: "${name[@]}" \
            :::+ "${name[@]/.bam/.bedpe}"

            echo -e "Adjusting bedpe \"pos\", \"mpos\", etc. values..."
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2 + 1, \$3 + 1, \$4, \$5 + 1, \$6 + 1, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.bedpe.tmp}"

            echo -e "Cleaning up bedpe.tmp files..."
            parallel -k -j "${parallelize}" \
            "mv -f {1} {2}" \
            ::: "${name[@]/.bam/.bedpe.tmp}" \
            :::+ "${name[@]/.bam/.bedpe}"

            echo -e "Creating \"pos\" bed files..."
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$1, \$2, \$3, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.pos.bed}"

            echo -e "Creating \"mpos\" bed files..."
            parallel -k -j "${parallelize}" \
            "awk 'BEGIN{FS=OFS=\"\t\"} {print \$4, \$5, \$6, \$7}' {1} > {2}" \
            ::: "${name[@]/.bam/.bedpe}" \
            :::+ "${name[@]/.bam/.mpos.bed}"

            echo -e "Removing temporary bedpe files..."
            parallel -k -j "${parallelize}" \
            "rm {1}" \
            ::: "${name[@]/.bam/.bedpe}"

            echo -e "Completed: Creating bed files for chromosome"
            echo -e "           ${chromosome}.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
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
        echo -e "#  Step 1 ----------------------------------------"
        echo -e "Started: Splitting bam infile into individual bam file for"
        echo -e "         chromosome ${chromosome}."

        samtools view -@ "${parallelize}" -b "${infile}" "${chromosome}" \
        > "${infile%.bam}.${chromosome}.bam" &
        displaySpinningIcon $! "Splitting bam infile..."

        echo -e "Completed: Splitting bam infile into individual bam file for"
        echo -e "           chromosome ${chromosome}.\n"


        #  Step 2: Index the split bam file ---------------
        echo -e "#  Step 2 ----------------------------------------"
        echo -e "Started: Indexing bam file for chromosome ${chromosome}."

        samtools index -@ "${parallelize}" "${infile%.bam}.${chromosome}.bam" &
        displaySpinningIcon $! "Indexing split bam file..."
                
        echo -e "Completed: Indexing bam file for chromosome ${chromosome}.\n"


        #  Step 3: Repair the split bam file --------------
        if [[ $((flag_subread)) -eq 1 ]]; then
            echo -e "#  Step 3 ----------------------------------------"
            echo -e "Started: Running Subread repair on the split bam file."
            
            repair \
            -d -c -T "${parallelize}" \
            -i "${infile%.bam}.${chromosome}.bam" \
            -o "${infile%.bam}.${chromosome}.bam.tmp" &
            displaySpinningIcon $! "Running repair on split bam file..."

            [[ -f "${infile%.bam}.${chromosome}.bam.tmp" ]] &&
                {
                    mv -f \
                    "${infile%.bam}.${chromosome}.bam.tmp" \
                    "${infile%.bam}.${chromosome}.bam"
                }

            echo -e "Completed: Running Subread repair on the split bam file\n"
        elif [[ $((flag_subread)) -eq 0 ]]; then
            :
        else
            echo -e "Exiting: An error occurred when setting the Subread"
            echo -e "         repair flag.\n"
            exit 1
        fi


        #  Step 4: Create "pos", "mpos" bed files ---------
        if [[ $((flag_bed)) -eq 1 ]]; then
            echo -e "#  Step 4 ----------------------------------------"
            echo -e "Started: Creating bed files for chromosome"
            echo -e "         \"${chromosome}.\""

            bamToBed -i "${infile%.bam}.${chromosome}.bam" -bedpe \
            > "${infile%.bam}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Creating temporary bedpe file..."

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2 + 1, $3 + 1, $4, $5 + 1, $6 + 1, $7}' \
            "${infile%.bam}.${chromosome}.bedpe" \
            > "${infile%.bam}.${chromosome}.bedpe.tmp" &
            displaySpinningIcon $! "Adjusting bedpe \"pos\", etc. values..."

            mv -f \
            "${infile%.bam}.${chromosome}.bedpe.tmp" \
            "${infile%.bam}.${chromosome}.bedpe" &
            displaySpinningIcon $! "Cleaning up bedpe.tmp files..."

            awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $7}' \
            "${infile%.bam}.${chromosome}.bedpe" \
            > "${infile%.bam}.${chromosome}.pos.bed" &
            displaySpinningIcon $! "Creating \"pos\" bed file..."

            awk 'BEGIN{FS=OFS="\t"} {print $4, $5, $6, $7}' \
            "${infile%.bam}.${chromosome}.bedpe" \
            > "${infile%.bam}.${chromosome}.mpos.bed" &
            displaySpinningIcon $! "Creating \"mpos\" bed file..."

            rm "${infile%.bam}.${chromosome}.bedpe"  &
            displaySpinningIcon $! "Removing temporary bedpe file..."

            echo -e "Completed: Creating bed files for chromosome"
            echo -e "           ${chromosome}.\n"
        elif [[ $((flag_bed)) -eq 0 ]]; then
            :
        else
            echo -e "Exiting: An error occurred when setting the bed flag.\n"
            exit 1
        fi
        ;;
    *) \
        echo "Exiting: Script cannot handle \"${chromosome}\" as chromosome option."
        exit 1
        ;;
esac


#  Scraps ---------------------------------------------------------------------

# echo "-l <perform liftOver of bed files: TRUE or FALSE (logical); note:"
# echo "    only need to call this argument if creating bed files from"
# echo "    split bam files>"
# echo "-s <strain for performing liftOver of bed files (int); current"
# echo "    options are CAST-EiJ (\"1\") and 129S1-SvImJ (\"2\"); note:"
# echo "    only need to call this argument if creating bed files from"
# echo "    split bam files>"

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

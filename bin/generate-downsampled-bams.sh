#!/bin/bash

#  generate-downsampled-bams.sh
#  KA


time_start="$(date +%s)"

#  Set up functions -----------------------------------------------------------
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


#  Parse arguments, assign variables ------------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Downsample a paired-end bam file not containing singletons"
    echo "to a user-specified number of reads while retaining mate"
    echo "information."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - BBMap >= 38.95"
    echo "  - bedtools >= 2.29.0"
    echo "  - repair >= 2.0.1 (part of Subread)"
    echo "  - samtools >= 1.13"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h  print this help message and exit"
    echo "-u  use safe mode: TRUE or FALSE [logical]; default: FALSE"
    echo "-i  deduplicated bam infile, including path [chr]"
    echo "-o  path for downsampled bam file [chr]; path will be made if it"
    echo "    does not exist"
    echo "-d  number of paired-end reads to sample down to [even int >= 2;"
    echo "    default: 300000]"
    echo "-x  prefix for downsampled paired-end bam outfiles [chr, optional];"
    echo "    if \"-x\" is blank, then prefix is derived from the name of the"
    echo "    infile"
    echo "-s  seed number for deterministic random sampling [int >= 1;"
    echo "    default: 24]"
    echo "-p  number of cores for parallelization [int >= 1; default: 4]"
    exit
}

interactive=FALSE  # Hardcode TRUE for interactive testing; FALSE for CLI
if [[ "${interactive}" == TRUE ]]; then
    safe_mode=FALSE
    infile="./data/files_bam/Disteche_sample_1.dedup.CAST.corrected.bam"
    # infile="./data/files_bam/Disteche_sample_1.dedup.mm10.corrected.bam"
    outpath="./data/files_bam"
    downsample=3000000
    prefix=""
    seed=24
    parallelize=4
else
    while getopts "h:u:i:o:x:d:s:p:" opt; do
        case "${opt}" in
            h) print_usage ;;
            u) safe_mode="${OPTARG}" ;;
            i) infile="${OPTARG}" ;;
            o) outpath="${OPTARG}" ;;
            d) downsample="${OPTARG}" ;;
            x) prefix="${OPTARG}" ;;
            s) seed="${OPTARG}" ;;
            p) parallelize="${OPTARG}" ;;
            *) print_usage ;;
        esac
    done
fi

[[ -z "${safe_mode}" ]] && safe_mode="FALSE"
[[ -z "${infile}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${prefix}" ]] && prefix=""
[[ -z "${downsample}" ]] && downsample=300000
[[ -z "${seed}" ]] && seed=24
[[ -z "${parallelize}" ]] && parallelize=4


#  Check user input -----------------------------------------------------------
check_dependency bedtools
check_dependency reformat.sh  # part of BBMap
check_dependency repair
check_dependency samtools

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-u: \"Safe mode\" is on.\n" && set -Eeuxo pipefail ;;
    false | f) echo -e "-u: \"Safe mode\" is off.\n" ;;
    *) \
        echo -e "Exiting: -u \"safe mode\" argument must be \"TRUE\" or \"FALSE\".\n"
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
        echo "-o: Directory ${outpath} does not exist; making the directory.\n"
        mkdir -p "${outpath}"
    }

#  Check "${downsample}", $((downsample))
[[ ! "${downsample}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an integer.\n"
        exit 1
    }

[[ ! $((downsample)) -ge 2 ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an integer >= 2.\n"
        exit 1
    }

[[ ! $((downsample % 2)) -eq 0 ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an even integer.\n"
        exit 1
    }

#  Check "${seed}", $((seed))
[[ ! "${seed}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -s seed argument must be an integer.\n"
        exit 1
    }

[[ ! $((seed)) -ge 1 ]] &&
    {
        echo -e "Exiting: -s seed argument must be an integer >= 1.\n"
        exit 1
    }

#  Check "${parallelize}", $((parallelize))
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


#  Downsample a bam file while maintaining mate-pairing information -----------
echo "Started: ${0}"

#  Prior to downsampling, get the mate pairs next to each other 
echo -e "Started: \"Repairing\" $(basename "${infile}")..."
repair -d -c -T "${parallelize}" \
-i "${infile}" \
-o "${infile%.bam}.tmp.bam"
echo -e "Completed: \"Repairing\" $(basename "${infile}")...\n"

#  Split bam infile by forward and reverse strands
#+ e.g., see biostars.org/p/420892/
echo -e "Started: Splitting forward strand to separate bam..."
samtools view -@ "${parallelize}" \
-F 16 "${infile%.bam}.tmp.bam" \
-o "${infile/.bam/.only-forward.bam}"
echo -e "Completed: Splitting forward strand to separate bam...\n"

echo -e "Started: Splitting reverse strand to separate bam..."
samtools view -@ "${parallelize}" \
-f 16 "${infile%.bam}.tmp.bam" \
-o "${infile/.bam/.only-reverse.bam}"
echo -e "Completed: Splitting reverse strand to separate bam...\n"

#  Remove unneeded intermediate file
rm "${infile%.bam}.tmp.bam"

#  Move "work files" into ${outpath}, then cd into ${outpath} for further work
mv -f \
"${infile/.bam/.only-forward.bam}" \
"${infile/.bam/.only-reverse.bam}" \
"${outpath}/"

cd "${outpath}" ||
    {
        echo "Exiting: An error occurred when cd'ing into ${outpath}"
        exit 1
    }

#  Set up/adjust needed variables
downsample_split=$((downsample / 2))
infile="$(basename "${infile}")"

#  Sample forward mate
echo -e "Started: Randomly sampling ${infile/.bam/.only-forward.bam} to ${downsample_split} reads..."
reformat.sh \
in="${infile/.bam/.only-forward.bam}" \
out="${infile/.bam/.only-forward.downsample.bam}" \
sampleseed="${seed}" \
samplereadstarget="${downsample_split}"
echo -e "Completed: Randomly sampling ${infile/.bam/.only-forward.bam} to ${downsample_split} reads...\n"

#  Sample reverse mate
echo -e "Started: Randomly sampling ${infile/.bam/.only-reverse.bam} to ${downsample_split} reads..."
reformat.sh \
in="${infile/.bam/.only-reverse.bam}" \
out="${infile/.bam/.only-reverse.downsample.bam}" \
sampleseed="${seed}" \
samplereadstarget="${downsample_split}"
echo -e "Completed: Randomly sampling ${infile/.bam/.only-reverse.bam} to ${downsample_split} reads...\n"

#  Remove unneeded intermediate files
rm \
"${infile/.bam/.only-forward.bam}" \
"${infile/.bam/.only-reverse.bam}"

#  Coordinate-sort sampled forward mate
echo -e "Started: Sorting ${infile/.bam/.only-forward.downsample.bam}..."
samtools sort -@ "${parallelize}" \
"${infile/.bam/.only-forward.downsample.bam}" \
-o "${infile/.bam/.only-forward.downsample.sort.bam}"
echo -e "Completed: Sorting ${infile/.bam/.only-forward.downsample.bam}...\n"

#  Coordinate-sort sampled reverse mate
echo -e "Started: Sorting ${infile/.bam/.only-reverse.downsample.bam}..."
samtools sort -@ "${parallelize}" \
"${infile/.bam/.only-reverse.downsample.bam}" \
-o "${infile/.bam/.only-reverse.downsample.sort.bam}"
echo -e "Completed: Sorting ${infile/.bam/.only-reverse.downsample.bam}...\n"

#  Remove unneeded intermediate files
rm \
"${infile/.bam/.only-forward.downsample.bam}" \
"${infile/.bam/.only-reverse.downsample.bam}"

#  Merged the coordinate-sorted, sampled mates
echo -e "Started: Merging..."
echo -e "  - ${infile/.bam/.only-forward.downsample.sort.bam}"
echo -e "  - ${infile/.bam/.only-reverse.downsample.sort.bam}"
echo -e "...to form ${infile/.bam/.downsample-}${downsample}.bam..."
samtools merge -@ "${parallelize}" \
"${infile/.bam/.only-forward.downsample.sort.bam}" \
"${infile/.bam/.only-reverse.downsample.sort.bam}" \
-o "${infile/.bam/.downsample-}${downsample}.bam"
echo -e "Completed: Merging..."
echo -e "  - ${infile/.bam/.only-forward.downsample.sort.bam}"
echo -e "  - ${infile/.bam/.only-reverse.downsample.sort.bam}"
echo -e "...to form ${infile/.bam/.downsample-}${downsample}.bam...\n"

#  Remove unneeded intermediate files
rm \
"${infile/.bam/.only-forward.downsample.sort.bam}" \
"${infile/.bam/.only-reverse.downsample.sort.bam}"

#  Coordinate-sort the repaired bam file
samtools sort -@ "${parallelize}" \
"${infile/.bam/.downsample-}${downsample}.bam" \
> "${infile/.bam/.downsample-}${downsample}.sort.bam"

#  Remove unneeded intermediate file
rm "${infile/.bam/.downsample-}${downsample}.bam"

#  If specified by the user, rename outfiles
if [[ "${prefix}" != "" ]]; then
    if [[ "${prefix: -4}" == ".bam" ]]; then
        #  Reset "${prefix}" to version in which ".bam" is stripped off
        prefix="${prefix%.bam}"

        #  Rename the sorted bam file
        mv -f "${infile/.bam/.downsample-}${downsample}.sort.bam" "${prefix}.bam"
        echo "${infile/.bam/.downsample-}${downsample}.sort.bam renamed to ${prefix}.bam"
    else
        #  Rename the sorted bam file
        mv -f "${infile/.bam/.downsample-}${downsample}.sort.bam" "${prefix}.bam"
        echo "${infile/.bam/.downsample-}${downsample}.sort.bam renamed to ${prefix}.bam"
    fi
else
    mv -f \
    "${infile/.bam/.downsample-}${downsample}.sort.bam" \
    "${infile/.bam/.downsample-}${downsample}.bam"
fi

#  Return to starting directory
cd - || ! echo "Warning: \"cd -\" failed. Look into this..."


#  End recording time ---------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

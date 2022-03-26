#!/bin/bash

#  generate-downsampled-bam.sh
#  KA


#  Start recording time -------------------------------------------------------
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


#  Parse arguments, assign variables ------------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Take a deduplicated paired-end bam file not containing singletons"
    echo "and sample it down to a user-specified number of reads while"
    echo "retaining mate information."
    echo ""
    echo ""
    echo "Dependencies:"
    echo " - BBMap >= 38.95"
    echo " - samtools >= 1.13"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: \"TRUE\" or \"FALSE\" (logical); defualt is"
    echo "    \"FALSE\">"
    echo "-i <deduplicated bam infile, including path (chr)>"
    echo "-o <path for downsampled bam file (chr); path will be made if it"
    echo "    does not exist>"
    echo "-d <number of paired-end reads to sample down to (even int >= 2);"
    echo "    default: 300000>"
    echo "-x <prefix for final downsampled paired-end bam (chr, optional); if"
    echo "    \"-x\" is undefined, then prefix is derived from the name of the"
    echo "    infile>"
    echo "-s <seed number for deterministic random sampling (int >= 1);"
    echo "    default: 24>"
    echo "-p <number of cores for parallelization (int >= 1); default: 4>"
    exit
}

while getopts "h:u:i:o:x:d:s:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        d) downsample="${OPTARG}" ;;
        x) prefix="${OPTARG}" ;;
        s) seed="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode="FALSE"
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${prefix}" ]] && prefix=""
[[ -z "${downsample}" ]] && downsample=300000
[[ -z "${seed}" ]] && seed=24
[[ -z "${parallelize}" ]] && parallelize=4


#  Check user input -----------------------------------------------------------
checkDependency samtools
checkDependency reformat.sh

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-u: Safe mode is on.\n"
        set -Eeuxo pipefail
        ;;
    false | f) \
        echo -e "-u: Safe mode is off.\n"
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
repair -d -c -T "${parallelize}" \
-i "${infile}" \
-o "${infile%.bam}.tmp.bam" &
displaySpinningIcon $! "\"Repairing\" ${infile}... "

#  biostars.org/p/420892/
samtools view -@ "${parallelize}" \
-F 16 "${infile%.bam}.tmp.bam" \
-o "${infile/.bam/.only-forward.bam}" &
displaySpinningIcon $! "Splitting forward strand to separate bam... " 

samtools view -@ "${parallelize}" \
-f 16 "${infile%.bam}.tmp.bam" \
-o "${infile/.bam/.only-reverse.bam}" &
displaySpinningIcon $! "Splitting reverse strand to separate bam... "

#  Remove unneeded intermediate file
rm "${infile%.bam}.tmp.bam"

#  Move "work files" into ${outpath}, then cd into ${outpath} for further work
mv -f "${infile/.bam/.only-forward.bam}" "${infile/.bam/.only-reverse.bam}" "${outpath}/"

cd "${outpath}" ||
    {
        echo "Exiting: An error occurred when cd'ing into ${outpath}"
        exit 1
    }

#  Set up/adjust needed variables
downsample_split=$((downsample / 2))
infile="$(basename "${infile}")"

#  Sample forward mate
reformat.sh \
in="${infile/.bam/.only-forward.bam}" \
out="${infile/.bam/.only-forward.downsample.bam}" \
sampleseed="${seed}" \
samplereadstarget="${downsample_split}" &
displaySpinningIcon $! "Randomly sampling ${infile/.bam/.only-forward.bam} to ${downsample_split} reads... "

#  Sample reverse mate
reformat.sh \
in="${infile/.bam/.only-reverse.bam}" \
out="${infile/.bam/.only-reverse.downsample.bam}" \
sampleseed="${seed}" \
samplereadstarget="${downsample_split}" &
displaySpinningIcon $! "Randomly sampling ${infile/.bam/.only-reverse.bam} to ${downsample_split} reads... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.only-forward.bam}" "${infile/.bam/.only-reverse.bam}"

#  Coordinate-sort sampled forward mate
samtools sort -@ "${parallelize}" \
"${infile/.bam/.only-forward.downsample.bam}" \
-o "${infile/.bam/.only-forward.downsample.sort.bam}" &
displaySpinningIcon $! "Sorting ${infile/.bam/.only-forward.downsample.bam}... "

#  Coordinate-sort sampled reverse mate
samtools sort -@ "${parallelize}" \
"${infile/.bam/.only-reverse.downsample.bam}" \
-o "${infile/.bam/.only-reverse.downsample.sort.bam}" &
displaySpinningIcon $! "Sorting ${infile/.bam/.only-reverse.downsample.bam}... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.only-forward.downsample.bam}" "${infile/.bam/.only-reverse.downsample.bam}"

#  Merged the coordinate-sorted, sampled mates
samtools merge -@ "${parallelize}" \
"${infile/.bam/.only-forward.downsample.sort.bam}" \
"${infile/.bam/.only-reverse.downsample.sort.bam}" \
-o "${infile/.bam/.downsample-}${downsample}.bam" &
displaySpinningIcon $! "Merging ${infile/.bam/.only-forward.downsample.sort.bam} and ${infile/.bam/.only-reverse.downsample.sort.bam} to form ${infile/.bam/.downsample-}${downsample}.bam... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.only-forward.downsample.sort.bam}" "${infile/.bam/.only-reverse.downsample.sort.bam}"

#  Coordinate-sort the repaired bam file
samtools sort -@ "${parallelize}" \
"${infile/.bam/.downsample-}${downsample}.bam" \
> "${infile/.bam/.downsample-}${downsample}.sort.bam"

#  Remove unneeded intermediate file
rm "${infile/.bam/.downsample-}${downsample}.bam"

#  If specified by the user, rename outfile
[[ "${prefix}" != "" ]] &&
    {
        mv "${infile/.bam/.downsample-}${downsample}.sort.bam" "${prefix}"
        echo "${infile/.bam/.downsample-}${downsample}.sort.bam renamed to ${prefix}"
    }

#  Return to starting directory
cd - || ! echo "Warning: \"cd -\" failed. Look into this..."


#  End recording time ---------------------------------------------------------
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""

exit 0

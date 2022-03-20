#!/bin/bash

#  generate_downsampled-bam.sh
#  KA


#  Start recording time -------------------------------------------------------
start="$(date +%s)"


#  Functions ------------------------------------------------------------------
checkDependency() {
    #  Check if program is available in "${PATH}"; exit if not
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            # exit 1
        }
}

displaySpinningIcon() {
    #  Display "spinning icon" while a background process runs
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
    echo " - BBMap >= 38.95 (untested w/earlier versions)"
    echo " - samtools >= 1.13 (untested w/earlier versions)"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: \"TRUE\" or \"FALSE\" (logical); defualt is"
    echo "    \"FALSE\">"
    echo "-i <deduplicated bam infile, including path (chr)>"
    echo "-o <path for downsampled bam file (chr); path will be made if it"
    echo "    does not exist>"
    echo "-n <name of final downsampled paired-end bam (chr, optional); if"
    echo "    \"-n\" is undefined, then name is derived from the name of the"
    echo "    infile>"
    echo "-d <number of paired-end reads to sample down to (even int >= 2);"
    echo "    default: \"300000\">"
    echo "-s <seed number for deterministic random sampling (int >= 1);"
    echo "    default: \"24\">"
    echo "-p <number of cores for parallelization (int >= 1); defualt: \"4\">"
    exit
}

while getopts "h:u:i:o:n:d:s:p:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        n) name="${OPTARG}" ;;
        d) downsample="${OPTARG}" ;;
        s) seed="${OPTARG}" ;;
        p) parallelize="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && safe_mode="FALSE"
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${name}" ]] && name=""
[[ -z "${downsample}" ]] && downsample=300000
[[ -z "${seed}" ]] && seed=24
[[ -z "${parallelize}" ]] && parallelize=4

# #  Test defaults
# infile="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/Disteche_sample_1.q10.sorted.merged.bam"
# outpath="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_bam_test"
# downsample=750000
# name="test.not-deduplicated.${downsample}.bam"
# seed=24
# parallelize=4


#  Check user input -----------------------------------------------------------
checkDependency samtools
checkDependency reformat.sh

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

#  Check "${downsample}", $((downsample))
[[ ! "${downsample}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an integer."
        # exit 1
    }

[[ ! $((downsample)) -ge 2 ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an integer >= 2."
        # exit 1
    }

[[ ! $((downsample % 2)) -eq 0 ]] &&
    {
        echo -e "Exiting: -d downsample argument must be an even integer."
        # exit 1
    }

#  Check "${seed}", $((seed))
[[ ! "${seed}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -s seed argument must be an integer."
        exit 1
    }

[[ ! $((seed)) -ge 1 ]] &&
    {
        echo -e "Exiting: -s seed argument must be an integer >= 1."
        exit 1
    }

#  Check "${parallelize}", $((parallelize))
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


#  Downsample a bam file while maintaining mate-pairing information -----------
echo "Started: ${0}"

#  biostars.org/p/420892/
samtools view -@ "${parallelize}" -F 16 "${infile}" -o "${infile/.bam/.only-forward.bam}" &
displaySpinningIcon $! "Splitting forward strand to separate bam... " 

samtools view -@ "${parallelize}" -f 16 "${infile}" -o "${infile/.bam/.only-reverse.bam}" &
displaySpinningIcon $! "Splitting reverse strand to separate bam... "

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
-o "${infile/.bam/.only-forward.downsample.sorted.bam}" &
displaySpinningIcon $! "Sorting ${infile/.bam/.only-forward.downsample.bam}... "

#  Coordinate-sort sampled reverse mate
samtools sort -@ "${parallelize}" \
"${infile/.bam/.only-reverse.downsample.bam}" \
-o "${infile/.bam/.only-reverse.downsample.sorted.bam}" &
displaySpinningIcon $! "Sorting ${infile/.bam/.only-reverse.downsample.bam}... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.only-forward.downsample.bam}" "${infile/.bam/.only-reverse.downsample.bam}"

#  Merged the coordinate-sorted, sampled mates
samtools merge -@ "${parallelize}" \
"${infile/.bam/.only-forward.downsample.sorted.bam}" \
"${infile/.bam/.only-reverse.downsample.sorted.bam}" \
-o "${infile/.bam/.downsample-}${downsample}.bam" &
displaySpinningIcon $! "Merging ${infile/.bam/.only-forward.downsample.sorted.bam} and ${infile/.bam/.only-reverse.downsample.sorted.bam} to form ${infile/.bam/.downsample-}${downsample}.bam... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.only-forward.downsample.sorted.bam}" "${infile/.bam/.only-reverse.downsample.sorted.bam}"

#  Get the mate pairs next to each other, resulting in the outfile
repair -d -c -T "${parallelize}" \
-i "${infile/.bam/.downsample-}${downsample}.bam" \
-o "${infile/.bam/.downsample-}${downsample}.repaired.bam" &
displaySpinningIcon $! "\"Repairing\" ${infile/.bam/.downsample-}${downsample}.bam... "

#  Remove unneeded intermediate files
rm "${infile/.bam/.downsample-}${downsample}.bam"

#  If specified by the user, rename outfile
[[ "${name}" != "" ]] &&
    {
        mv "${infile/.bam/.downsample-}${downsample}.repaired.bam" "${name}"
        echo "${infile/.bam/.downsample-}${downsample}.repaired.bam renamed to ${name}"
    }

#  Return to starting directory
cd - ||
    {
        echo "Warning: \"cd -\" failed. Look into this..."
    }


#  End recording time ---------------------------------------------------------
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Completed: ${0}"
echo "${0} run time: ${run_time} seconds."
echo ""

exit 0

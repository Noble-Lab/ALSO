#!/bin/bash

#  submit_bulk-RNA-seq.sh
#  KA


# :param 1: name of experiment directory, including path <chr>
# :param 2: bam infile, including path <chr>
# :param 3: step in pipeline to run up to <int 1-4>
# :param 4: JVM min heap value; do not include '-Xms' <int > 0>
# :param 5: JVM min heap value; do not include '-Xmx' <int > 0 && int ≥ param 4>


#  Capture arguments ----------------------------------------------------------
arguments=("$@")

#  Stop function if insufficient numbers of arguments passed
if [[ ! "${#arguments[@]}" -eq 5 ]]; then
    echo "Stopping: ${#arguments[@]} arguments were passed; five arguments"\
    "should be passed:"
    echo ":param 1: name of experiment directory, including path <chr]; for" \
    "example, 2022-0617_test_driver_allelic-segregation_Disteche-sample-11_downsample-6000000"
    echo ":param 2: bam infile, including path <chr>"
    echo ":param 3: step in pipeline to run up to <int 1-4>"
    echo ":param 4: JVM min heap value; do not include '-Xms', 'm', 'g', etc. <int > 0>"
    echo ":param 5: JVM min heap value; do not include '-Xmx', 'm', 'g', etc. <int > 0 && int ≥ param 4>"
    exit 1
fi


#  Establish relevant directories ---------------------------------------------
dir_base="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross"
dir_experiment="${1}"

cd "${dir_base}" ||
    {
        echo "Exiting: Could not cd into ${dir_base}. Check on this."
        exit 1
    }


#  Activate appropriate environment -------------------------------------------
# shellcheck disable=1091
. "${HOME}/.bashrc"

# shellcheck disable=1091
. activate pipeline-test_env


#  Set up variables -----------------------------------------------------------
sample="${2}"
prefix=$(basename "${sample}")
echo "Working with sample ${prefix}"

end="single"
suffix=".Aligned.sortedByCoord.out.primary.rm.bam"
strain_1="mm11"
strain_2="SPRET"
path_file_1="${sample}.${strain_1}${suffix}"
path_file_2="${sample}.${strain_2}${suffix}"
dir_results="${dir_experiment}/results/"
run_up_to="${3}"

dir_log="${dir_base}/log/$(basename "${1}")"
mkdir -p "${dir_log}"


#  Run ALSO -------------------------------------------------------------------
#  Echo test
echo -e "\
bash \"${dir_base}/bin/workflow/driver_ALSO.sh\"\n\
-u FALSE\n\
-d FALSE\n\
-m \"-Xms${4}m\"\n\
-x \"-Xmx${5}m\"\n\
-r \"${strain_1}\"\n\
-s \"${strain_2}\"\n\
-1 \"${path_file_1}\"\n\
-2 \"${path_file_2}\"\n\
-p \"${prefix}\"\n\
-o \"${dir_results}\"\n\
-e \"${end}\" \n\
-b 100000\n\
-c 1000000\n\
-t 0\n\
-a TRUE\n\
-n \"${run_up_to}\"\n\
> \"${dir_log}/${prefix}.stdout.txt\"\n\
2> \"${dir_log}/${prefix}.stderr.txt\"
"

#  Actually running ALSO
bash "${dir_base}/bin/workflow/driver_ALSO.sh" \
-u FALSE \
-d FALSE \
-m "-Xms${4}m" \
-x "-Xmx${5}m" \
-r "${strain_1}" \
-s "${strain_2}" \
-1 "${path_file_1}" \
-2 "${path_file_2}" \
-p "${prefix}" \
-o "${dir_results}" \
-e "${end}" \
-b 100000 \
-c 1000000 \
-t 0 \
-a TRUE \
-n "${run_up_to}" \
> "${dir_log}/${prefix}.stdout.txt" \
2> "${dir_log}/${prefix}.stderr.txt"

# Arguments:
# -h  print this help message and exit
# -u  use safe mode: TRUE or FALSE <logical> [default: FALSE]
# -d  run pipeline in ${TMPDIR}: TRUE or FALSE <logical> [default:
#     TRUE]
# -m  initial memory allocation pool (JVM) <chr> [default: "-Xms512m"]
# -x  maximum memory allocation pool (JVM) <chr> [default: "-Xmx4096m"]
# -r  string for "sample #1" <chr>
# -s  string for "sample #2" <chr>
# -1  bam infile #1, including path <chr>
# -2  bam infile #2, including path <chr>
# -p  prefix for outfiles <chr>
# -o  results directory for outfiles <chr>; path will be made if it
#     does not exist
# -e  reads from paired- or single-end sequencing run: "paired" or
#     "single" <chr> [default: "paired"]
# -b  number of records to read into memory at one time when running
#     the script for Part #1, get-AS-per-qname.R <int > 0> [default:
#     100000]
# -c  number of records to read into memory at one time when running
#     the script for Part #3, generate-assignment-lists.R <int > 0>
#     [default: 1000000]
# -t  alignment score threshold <int >= 0> [default: 0]; the absolute
#     value of the difference in alignment scores between "sample
#     #1" and "sample #2" must be greater than this value in order
#     for a sample-specific assignment to be made; if not greater than
#     this value, then the assignment will be "ambiguous"
# -a  count lines: TRUE or FALSE <logical> [default: TRUE]
# -n  step in pipeline to run up to <int 1-4> [default: 4]

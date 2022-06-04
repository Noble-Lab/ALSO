#!/bin/bash

#  driver_allelic-segregation.sh
#  KA


time_start="$(date +%s)"

#  Source functions into environment ------------------------------------------
# shellcheck disable=1091
if [[ -f "./bin/auxiliary/functions-preprocessing-HPC.sh" ]]; then
    . ./bin/auxiliary/functions-preprocessing-HPC.sh ||
        {
            echo -e "Exiting: Unable to source 'functions-preprocessing-HPC.sh'.\n"
            # exit 1
        }

    . ./bin/auxiliary/functions-in-progress.sh ||
        {
            echo -e "Exiting: Unable to source 'functions-preprocessing-HPC.sh'.\n"
            # exit 1
        }
elif [[ -f "./functions-preprocessing-HPC.sh" ]]; then
    . ./functions-preprocessing-HPC.sh ||
        {
            echo -e "Exiting: Unable to source 'functions-preprocessing-HPC.sh'.\n"
            # exit 1
        }

    . ./functions-in-progress.sh ||
        {
            echo -e "Exiting: Unable to source 'functions-in-progress.sh'.\n"
            # exit 1
        }
else
    echo -e "Exiting: Couldn't find auxiliary information.\n"
    # exit 1
fi

    
echo_completion_file() {
    # :param 1: outpath (chr)
    # :param 2: string description (chr)
    # :param 3: step number (int)
    echo "${1}/${2}.allelic-segregation.step-${3}.txt"
}

# :param 1: directory in which to search
# :param 2: pattern and/or name to identify file
find_files() { find "${1}" -maxdepth 1 -type f -name "${2}"; }


#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo "Arguments:"
    echo "-h  print this help message and exit"
    echo "-u  use safe mode: TRUE or FALSE [logical; default: FALSE]"
    echo "-l  run on GS HPC: TRUE or FALSE [logical; default: FALSE]"
    echo "-m  initial memory allocation pool for JVM [chr; default: \"512m\"]"
    echo "-x  maximum memory allocation pool for JVM [chr; default: \"4096m\"]"
    echo "-r  string for sample #1 [chr]"
    echo "-s  string for sample #2 [chr]"
    echo "-1  bam infile #1, including path [chr]"
    echo "-2  bam infile #2, including path [chr]"
    echo "-p  prefix for outfiles [chr]"
    echo "-o  results directory for outfiles [chr]; path will be made if it"
    echo "    does not exist"
    echo "-c  number of records to read into memory at one time when running"
    echo "    step #1 [int > 0; default: 100000]"
    echo "-d  number of records to read into memory at one time when running"
    echo "    step #3 [int > 0; default: 1000000]"
    echo "-t  alignment score threshold [int >= 0; default: 0]; the absolute"
    echo "    value of the difference in alignment scores between samples #1"
    echo "    and #2 must be greater than this value in order for a sample-"
    echo "    specific assignment to be made; if not greater than this value,"
    echo "    then the assignment will be \"ambiguous\""
    echo "-m  count lines: TRUE or FALSE [logical; default: TRUE]"
    echo "-n  step in pipeline to run up to [int 1-4; default: 4]"
    # exit
}


interactive=TRUE  # Hardcode TRUE for interactive testing; FALSE for CLI
if [[ "${interactive}" == "TRUE" ]]; then
    safe_mode=FALSE
    cluster=FALSE
    memory_min="512m"
    memory_max="4096m"
    strain_1="mm10"
    strain_2="CAST"
    bam_1="./data/files_bam/Disteche_sample_1.dedup.mm10.corrected.downsample-3000000.bam"
    bam_2="./data/files_bam/Disteche_sample_1.dedup.CAST.corrected.downsample-3000000.bam"
    prefix="Disteche_sample_1"
    outpath="./results/kga0/2022-0603_allelic-segregation"
    chunk_step_1=100000
    chunk_step_3=1000000
    threshold=0
    count=TRUE
    run_up_to=4
else
    while getopts "h:u:l:m:x:r:s:1:2:p:o:c:d:t:a:n:" opt; do
        case "${opt}" in
            h) print_usage ;;
            u) safe_mode="${OPTARG}" ;;
            l) cluster="${OPTARG}" ;;
            m) memory_min="${OPTARG}" ;;
            x) memory_max="${OPTARG}" ;;
            r) strain_1="${OPTARG}" ;;
            s) strain_2="${OPTARG}" ;;
            1) bam_1="${OPTARG}" ;;
            2) bam_2="${OPTARG}" ;;
            p) prefix="${OPTARG}" ;;
            o) outpath="${OPTARG}" ;;
            c) chunk_step_1="${OPTARG}" ;;
            d) chunk_step_3="${OPTARG}" ;;
            t) threshold="${OPTARG}" ;;
            a) count="${OPTARG}" ;;
            n) run_up_to="${OPTARG}" ;;
            *) print_usage ;;
        esac
    done
fi

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${cluster}" ]] && cluster=FALSE
[[ -z "${memory_min}" ]] && memory_min="512m"
[[ -z "${memory_max}" ]] && memory_max="4096m"
[[ -z "${strain_1}" ]] && print_usage
[[ -z "${strain_2}" ]] && print_usage
[[ -z "${bam_1}" ]] && print_usage
[[ -z "${bam_2}" ]] && print_usage
[[ -z "${prefix}" ]] && print_usage
[[ -z "${outpath}" ]] && print_usage
[[ -z "${chunk_step_1}" ]] && chunk_step_1=100000
[[ -z "${chunk_step_3}" ]] && chunk_step_3=1000000
[[ -z "${threshold}" ]] && threshold=0
[[ -z "${count}" ]] && count=TRUE
[[ -z "${run_up_to}" ]] && run_up_to=4


#  Check variable assignments from arguments ----------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-u: \"Safe mode\" is TRUE." && set -Eeuxo pipefail ;;
    false | f) echo -e "-u: \"Safe mode\" is FALSE." ;;
    *) \
        echo -e "Exiting: -u \"safe mode\" argument must be TRUE or FALSE.\n"
        # exit 1
        ;;
esac

#  Check for necessary dependencies; exit if not found
check_dependency picard
check_dependency R
check_dependency samtools

#  Evaluate "${cluster}"
case "$(echo "${cluster}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-l: \"Run on GS HPC\" is TRUE." ;;
    false | f) echo -e "-l: \"Run on GS HPC\" is FALSE." && check_dependency picard ;;
    *) \
        echo -e "Exiting: -l \"run on GS HPC\" argument must be TRUE or FALSE.\n"
        # exit 1
        ;;
esac

#  Make sure "${strain_1}" and "${strain_2}" are not the same
[[ "${strain_1}" == "${strain_2}" ]] &&
    {
        echo -e "Exiting: -r \"strain 1\" is the same as -s \"strain 2\".\n"
        # exit 1
    }

#  Evaluate "${bam_1}"
[[ -f "${bam_1}" ]] ||
    {
        echo -e "Exiting: -1 \"bam 1\" does not exist.\n"
        # exit 1
    }

#  Evaluate "${bam_2}"
[[ -f "${bam_2}" ]] ||
    {
        echo -e "Exiting: -2 \"bam 2\" does not exist.\n"
        # exit 1
    }

#  Make sure "${bam_1}" and "${bam_2}" are not the same
[[ $(basename "${bam_1}") == $(basename "${bam_2}") ]] &&
    {
        echo -e "Exiting: -1 \"bam 1\" is apparently the same as -2 \"bam 2\".\n"
        # exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo -e "-o: Directory \"${outpath}\" does not exist; making the directory.\n"
        mkdir -p "${outpath}"
    }

#  Evaluate "${threshold}"
[[ ! "${threshold}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -n \"threshold\" argument must be an integer >= 0.\n"
        # exit 1
    }

#  Evaluate "${count}"
case "$(echo "${count}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-c: \"Count lines\" is TRUE." ;;
    false | f) echo -e "-c: \"Count lines\" is FALSE." ;;
    *) \
        echo -e "Exiting: -c \"count lines\" argument must be TRUE or FALSE.\n"
        # exit 1
        ;;
esac

#  Evaluate "${run_up_to}"
[[ ! "${run_up_to}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -n \"run_up_to\" argument must be an integer in the range of 1-4.\n"
        # exit 1
    }

case "${run_up_to}" in
    1 | 2 | 3 | 4) \
        echo -e "-n: \"Run up to\" is set to step #${run_up_to}."
        ;;
    *) \
        echo -e "Exiting: -n \"run_up_to\" argument must be an integer in the range of 1-4.\n"
        # exit 1
        ;;
esac


#  Check that pipeline scripts are accessible ---------------------------------
#+ ...if so, then assign them to variables
if [[ -f "./bin/workflow/get-AS-per-qname.R" ]]; then
    script_1="./bin/workflow/get-AS-per-qname.R"
    :
elif [[ -f "./get-AS-per-qname.R" ]]; then
    script_1="./get-AS-per-qname.R"
    :
else
    echo -e "Exiting: Couldn't find \"get-AS-per-qname.R\"."
    # exit 1
fi

if [[ -f "./bin/workflow/find-set-intersection-set-complement.sh" ]]; then
    script_2="./bin/workflow/find-set-intersection-set-complement.sh"
    :
elif [[ -f "./find-set-intersection-set-complement.sh" ]]; then
    script_2="./find-set-intersection-set-complement.sh"
    :
else
    echo -e "Exiting: Couldn't find \"find-set-intersection-set-complement.sh\"."
    # exit 1
fi

if [[ -f "./bin/workflow/generate-assignment-lists.R" ]]; then
    script_3="./bin/workflow/generate-assignment-lists.R"
    :
elif [[ -f "./generate-assignment-lists.R" ]]; then
    script_3="./generate-assignment-lists.R"
    :
else
    echo -e "Exiting: Couldn't find \"generate-assignment-lists.R\"."
    # exit 1
fi

if [[ -f "./bin/workflow/filter-qnames-by-assignment.sh" ]]; then
    script_4="./bin/workflow/filter-qnames-by-assignment.sh"
    :
elif [[ -f "./filter-qnames-by-assignment.sh" ]]; then
    script_4="./filter-qnames-by-assignment.sh"
    :
else
    echo -e "Exiting: Couldn't find \"filter-qnames-by-assignment.sh\"."
    # exit 1
fi

# echo "${script_1}"
# echo "${script_2}"
# echo "${script_3}"
# echo "${script_4}"


#  Assign variables for outdirectories ----------------------------------------
dir_experiment_1="${outpath}/01_run_get-AS-per-qname"
dir_experiment_2="${outpath}/02_run_find-set-inter-set-complement"
dir_experiment_3="${outpath}/03_run_generate-assignment-lists"
dir_experiment_4="${outpath}/04_run_filter-qnames-by-assignment"

# echo "${dir_experiment_1}"
# echo "${dir_experiment_2}"
# echo "${dir_experiment_3}"
# echo "${dir_experiment_4}"


#  Assign variables for completion files --------------------------------------
step_1="$(echo_completion_file "${outpath}" "${prefix}" 1)"
step_2="$(echo_completion_file "${outpath}" "${prefix}" 2)"
step_3="$(echo_completion_file "${outpath}" "${prefix}" 3)"
step_4="$(echo_completion_file "${outpath}" "${prefix}" 4)"

# echo "${step_1}"
# echo "${step_2}"
# echo "${step_3}"
# echo "${step_4}"


#  Create outdirectories ------------------------------------------------------
mkdir -p "${dir_experiment_1}"
mkdir -p "${dir_experiment_2}"
mkdir -p "${dir_experiment_3}"
mkdir -p "${dir_experiment_4}"


#  Step 0 ---------------------------------------------------------------------
[[ ! -f "${bam_1}.bai" ]] &&
    {
        echo -e "Started step 0/4: Running samtools index for ${bam_1}."
        
        samtools index "${bam_1}"

        echo -e "Completed step 0/4: Running samtools index for ${bam_1}."
    }

[[ ! -f "${bam_2}.bai" ]] &&
    {
        echo -e "Started step 0/4: Running samtools index for ${bam_2}."
        
        samtools index "${bam_2}"

        echo -e "Completed step 0/4: Running samtools index for ${bam_2}."
    }

evaluate_run_up_to "${run_up_to}" 0


#  Step 1 ---------------------------------------------------------------------
#   -b, --bam     bam infile, including path <chr>
#   -i, --bai     bam index, including path <chr>
#   -o, --outdir  directory for saving rds outfile, including path <chr>
#   -s, --strain  strain name to be appended to rds outfile columns <chr>
#   -c, --chunk   number of records to read into memory at a single time
#                 <even int> [default: 100000]

# script_1="./get-AS-per-qname.R"
# dir_experiment_1="${outpath}/run_get-AS-per-qname"
#TODO (   ) Provide more informative messages output by "get-AS-per-qname.R"
#TODO (   ) Clean up the message formatting
#TODO ( Y ) Fix the name of outfiles, which mention the strain twice
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started step 1/4: Running ${script_1}."
    
    #  Create key-value array for samples #1, #2 (should work with different
    #+ versions of bash, zsh, etc.)
    unset array
    declare -A array
    array["${strain_1}"]="${bam_1}"
    array["${strain_2}"]="${bam_2}"

    for i in "${strain_1}" "${strain_2}"; do
        bam="${array["${i}"]}"
        bai="${array["${i}"]}.bai"
        strain="${i}"

        if [[ -f "${bai}" && -f "${bam}" ]]; then
            echo -e \
            " \
            Running...\n \
            Rscript \"${script_1}\"\n \
            --bam \"${bam}\"\n \
            --bai \"${bai}\"\n \
            --outdir \"${dir_experiment_1}\"\n \
            --strain \"${strain}\"\n \
            --chunk \"${chunk_step_1}\"\n"

            Rscript "${script_1}" \
            --bam "${bam}" \
            --bai "${bai}" \
            --outdir "${dir_experiment_1}" \
            --strain "${strain}" \
            --chunk "${chunk_step_1}"
        else
            echo "Exiting: ${bai} and/or ${bam} not found."
            # exit 1
        fi
    done
fi

AS_1=$(find_files "${dir_experiment_1}" "*${strain_1}*AS.txt.gz")
AS_2=$(find_files "${dir_experiment_1}" "*${strain_2}*AS.txt.gz")

if [[ -f "${AS_1}" && -f "${AS_2}" ]]; then
    touch "${step_1}"
    echo -e "Completed step 1/4: Running ${script_1}.\n"
fi

evaluate_run_up_to "${run_up_to}" 1


#  Step 2 ---------------------------------------------------------------------
# -h print this help message and exit
# -u use safe mode: TRUE or FALSE (logical; default: FALSE)
# -i AS.txt.gz "infile #1", including path (chr)
# -j AS.txt.gz "infile #2", including path (chr)
# -1 string for "sample #1" (chr)
# -2 string for "sample #2" (chr)
# -o path for outfiles (chr); path will be made if it does not exist
# -p prefix for outfiles (chr)
# -c count lines: TRUE or FALSE (logical)

# script_2="./find-set-intersection-set-complement.sh"
# dir_experiment_2="${outpath}/run_find-set-inter-set-complement"
#TODO ( Y  ) Script uses display_spinning_icon(); need to remove this
if [[ ! -f "${step_2}" && -f "${step_1}" ]]; then
    echo -e "Started step 2/4: Running ${script_2}."
    
    if [[ -f "${AS_1}" ]] && [[ -f "${AS_2}" ]]; then
        echo -e \
        " \
        Running...\n \
        bash \"${script_2}\"\n \
        -u \"${safe_mode}\"\n \
        -i \"${AS_1}\"\n \
        -j \"${AS_2}\"\n \
        -1 \"${strain_1}\"\n \
        -2 \"${strain_2}\"\n \
        -o \"${dir_experiment_2}\"\n \
        -p \"${prefix}\"\n \
        -c \"${count}\"\n"

        bash "${script_2}" \
        -u "${safe_mode}" \
        -i "${AS_1}" \
        -j "${AS_2}" \
        -1 "${strain_1}" \
        -2 "${strain_2}" \
        -o "${dir_experiment_2}" \
        -p "${prefix}" \
        -c "${count}"
    else
        echo "Exiting: ${AS_1} and/or ${AS_2} not found."
        # exit 1
    fi
fi

comp_1=$(find_files "${dir_experiment_2}" "*${strain_1}*complement*")
comp_2=$(find_files "${dir_experiment_2}" "*${strain_2}*complement*")
inter=$(find_files "${dir_experiment_2}" "*${strain_1}-${strain_2}*inter*")

if [[ -f "${comp_1}" && -f "${comp_2}" && -f "${inter}" ]]; then
    touch "${step_2}"
    echo -e "Completed step 2/4: Running ${script_2}.\n"
fi

evaluate_run_up_to "${run_up_to}" 2


#  Step 3 ---------------------------------------------------------------------
#   -i, --intersection  tab-separated txt file (gzipped or not) without
#                       header, including path, containing intersecting
#                       query names (QNAME) and alignment scores (AS) for
#                       two samples (chr)
#
#                       For example, the first five rows of a given txt
#                       file appears as follows, where column 1 is the
#                       QNAME, column 2 is the AS for sample 1, and
#                       column 3 is the AS for sample 3...
#   -s, --sample_1      sample name for sample 1; corresponds to the AS
#                       in column 2 of the tab- separated inter
#                       file (chr)
#   -u, --sample_2      sample name for sample 2; corresponds to the AS
#                       in column 3 of the tab- separated inter
#                       file (chr)
#   -o, --outdir        directory for saving txt.gz outfiles, including
#                       path (chr)
#   -p, --outprefix     prefix for outfiles (chr)
#   -t, --threshold     alignment score threshold; the absolute value of
#                       the difference in alignment scores between
#                       samples 1 and 2 must be greater than this value
#                       in order for a sample-specific assignment to be
#                       made; if not, the assignment will be "ambiguous"
#                       (int > 0)
#
#                       For example...  [default: 0]
#   -c, --chunk         number of records to read into memory at a single
#                       time (int ≥ 0) [default: 1000000]
#   -r, --remove        if present, remove outfiles in the outdirectory
#                       prior to generating outfiles (logical) [default:
#                       TRUE]

# script_3="./generate-assignment-lists.R"
# dir_experiment_3="${outpath}/run_generate-assignment-lists"
#TODO (   ) Provide more informative messages output by generate-assignment-lists.R
#TODO (   ) Clean up the message formatting
if [[ ! -f "${step_3}" && -f "${step_2}" ]]; then
    echo -e "Completed step 3/4: Running ${script_3}."

    if [[ -f "${comp_1}" && -f "${comp_2}" && -f "${inter}" ]]; then
        echo -e \
        " \
        Running...\n \
        Rscript \"${script_3}\"\n \
        --intersection \"${inter}\"\n \
        --sample_1 \"${strain_1}\"\n \
        --sample_2 \"${strain_2}\"\n \
        --outdir \"${dir_experiment_3}\"\n \
        --outprefix \"${prefix}\"\n \
        --threshold \"${threshold}\"\n \
        --chunk \"${chunk_step_3}\"\n \
        --remove TRUE"

        Rscript "${script_3}" \
        --intersection "${inter}" \
        --sample_1 "${strain_1}" \
        --sample_2 "${strain_2}" \
        --outdir "${dir_experiment_3}" \
        --outprefix "${prefix}" \
        --threshold "${threshold}" \
        --chunk "${chunk_step_3}" \
        --remove TRUE
    else
        echo "Exiting: ${inter} not found."
        # exit 1
    fi
fi

assign_1=$(find_files "${dir_experiment_3}" "*${strain_1}.txt.gz")
assign_2=$(find_files "${dir_experiment_3}" "*${strain_2}.txt.gz")
ambiguous=$(find_files "${dir_experiment_3}" "*ambiguous.txt.gz")

if [[ -f "${assign_1}" && -f "${assign_2}" && -f "${ambiguous}" ]]; then
    touch "${step_3}"
    echo -e "Completed step 3/4: Running ${script_3}.\n"
fi

evaluate_run_up_to "${run_up_to}" 3


#  Step 4 ---------------------------------------------------------------------
# -h print this help message and exit
# -s use safe mode: TRUE or FALSE (logical)
# -l run on GS HPC: TRUE or FALSE (logical)
# -m initial memory allocation pool for JVM (chr; default "512m")
# -x maximum memory allocation pool for JVM (chr; default "1g")
# -b bam infile #1, including path (chr)
# -c bam infile #2, including path (chr)
# -a "ambiguous" assignment file, including path (chr)
# -i "sample #1" assignment file, including path (chr)
# -j "sample #2" assignment file, including path (chr)
# -u "sample #1" unique file, including path (chr)
# -v "sample #2" unique file, including path (chr)
# -1 string for "sample #1" (chr)
# -2 string for "sample #2" (chr)
# -o path for outfiles (chr); path will be made if it does not exist
# -p prefix for outfiles (chr)
# -n count lines: TRUE or FALSE (logical)
# #TODO -f run samtools flagstat: TRUE or FALSE (logical)

# script_4="./filter-qnames-by-assignment.sh"
# dir_experiment_4="${outpath}/run_filter-qnames-by-assignment"
#TODO (   ) Clean up the message formatting
#TODO (   ) Have picard be less verbose
if [[ ! -f "${step_4}" && -f "${step_3}" ]]; then
    echo -e "Completed step 4/4: Running ${script_3}."

    if [[ -f "${assign_1}" && -f "${assign_2}" && "${ambiguous}" ]]; then
        echo -e \
        " \
        Running...\n \
        bash \"${script_4}\"\n \
        -s \"${safe_mode}\"\n \
        -l \"${cluster}\"\n \
        -m \"${memory_min}\"\n \
        -x \"${memory_max}\"\n \
        -b \"${bam_1}\"\n \
        -c \"${bam_2}\"\n \
        -a \"${ambiguous}\"\n \
        -i \"${assign_1}\"\n \
        -j \"${assign_2}\"\n \
        -u \"${comp_1}\"\n \
        -v \"${comp_2}\"\n \
        -1 \"${strain_1}\"\n \
        -2 \"${strain_2}\"\n \
        -o \"${dir_experiment_4}\"\n \
        -p \"${prefix}\"\n \
        -n \"${count}\"\n"

        bash "${script_4}" \
        -s "${safe_mode}" \
        -l "${cluster}" \
        -m "${memory_min}" \
        -x "${memory_max}" \
        -b "${bam_1}" \
        -c "${bam_2}" \
        -a "${ambiguous}" \
        -i "${assign_1}" \
        -j "${assign_2}" \
        -u "${comp_1}" \
        -v "${comp_2}" \
        -1 "${strain_1}" \
        -2 "${strain_2}" \
        -o "${dir_experiment_4}" \
        -p "${prefix}" \
        -n "${count}"
    else
        echo "Exiting: ${assign_1}, ${assign_2}, and/or ${ambiguous} not found."
        # exit 1
    fi
fi

s1_s1=$(find_files "${dir_experiment_4}" "*${strain_1}*${strain_1}*bam")
s1_s2=$(find_files "${dir_experiment_4}" "*${strain_1}*${strain_2}*bam")
s1_ambiguous=$(find_files "${dir_experiment_4}" "*${strain_1}*ambiguous*bam")
s2_s1=$(find_files "${dir_experiment_4}" "*${strain_2}*${strain_1}*bam")
s2_s2=$(find_files "${dir_experiment_4}" "*${strain_2}*${strain_2}*bam")
s2_ambiguous=$(find_files "${dir_experiment_4}" "*${strain_2}*ambiguous*bam")

if \
    [[ -f "${s1_s1}" && -f "${s1_s2}" && -f "${s1_ambiguous}" ]] &&
    [[ -f "${s2_s1}" && -f "${s2_s2}" && -f "${s2_ambiguous}" ]]; \
then
    touch "${step_4}"
    echo -e "Completed step 4/4: Running ${script_4}.\n"
fi


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

#!/bin/bash

#  driver_allelic-segregation.sh
#  KA


time_start="$(date +%s)"

#TODO (   ) Messaging for when a step has started and completed
#TODO (   ) Replace use non-"*_updated" with use of "*_updated" functions
#  Source functions into environment ------------------------------------------
# shellcheck disable=1091
if [[ -f "./bin/auxiliary/functions-preprocessing-HPC.sh" ]]; then
    . ./bin/auxiliary/functions-preprocessing-HPC.sh ||
        {
            echo -e "Exiting: Unable to source \
            'functions-preprocessing-HPC.sh'.\n"
            exit 1
        }

    . ./bin/auxiliary/functions-in-progress.sh ||
        {
            echo -e "Exiting: Unable to source \
            'functions-preprocessing-HPC.sh'.\n"
            exit 1
        }
elif [[ -f "./functions-preprocessing-HPC.sh" ]]; then
    . ./functions-preprocessing-HPC.sh ||
        {
            echo -e "Exiting: Unable to source \
            'functions-preprocessing-HPC.sh'.\n"
            exit 1
        }

    . ./functions-in-progress.sh ||
        {
            echo -e "Exiting: Unable to source \
            'functions-in-progress.sh'.\n"
            exit 1
        }
else
    echo -e "Exiting: Couldn't find auxiliary information.\n"
    exit 1
fi

    
echo_completion_file_driver() {
    # #TODO Description of function
    #
    # :param 1: outpath (chr)
    # :param 2: string description (chr)
    # :param 3: step number (int)
    echo "${1}/$(basename "${0}" ".sh").${2}.step-${3}.txt"
}


find_files() {
    # #TODO Description of function
    #
    # :param 1: directory in which to search
    # :param 2: pattern and/or name to identify file
    find "${1}" -maxdepth 1 -type f -name "${2}"
}


#TODO (   ) For the individual scripts in this pipeline, remove transferring
#           to/from TMPDIR; then, include an option to do so in this driver
#           script
#  Handle arguments, assign variables -----------------------------------------
print_usage() {
    echo ""
    echo "${0}:"
    echo "Run the allelic-segregation pipeline."
    echo "  - Part #1  <get-AS-per-qname.R>"
    echo "             For \"sample #1\" and \"sample #2\" bam files, output"
    echo "             sorted, tab-separated, gzipped txt files of bam"
    echo "             querynames (QNAMEs) and alignment scores (AS's)."
    echo ""
    echo "  - Part #2  <find-set-intersection-set-complement.sh>"
    echo "             Find the set intersection and set complements for tab-"
    echo "             separated, gzipped txt files containing QNAMEs and"
    echo "             AS's."
    echo ""
    echo "  - Part #3  <generate-assignment-lists.R>"
    echo "             Comparing AS's, assign set-intersection QNAMEs to"
    echo "             categories \"sample #1\", \"sample #2\", or"
    echo "             \"ambiguous\", saving the categorized QNAMEs in"
    echo "             distinct sorted, gzipped txt files for each category."
    echo ""
    echo "  - Part #4  <filter-qnames-by-assignment.sh>"
    echo "             Use assignment category txt files from to generate"
    echo "             assignment-specific bam files via picard"
    echo "             FilterSamReads."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "  - argparser >= #TODO"
    echo "  - bedtools >= 2.29.0"
    echo "  - GNU echo >= #TODO"
    echo "  - GNU find >= #TODO #TODO Does BSD work also?"
    echo "  - GNU zcat >= 1.09 #TODO Check on this..."
    echo "  - picard >= 2.27.1"
    echo "  - R >= 4.0"
    echo "  - Rscript >= 4.0"
    echo "  - Rsamtools >= #TODO"
    echo "  - samtools >= 1.13"
    echo "  - Subread repair >= 2.0.1"
    echo "  - scales >= #TODO"
    echo "  - tidyverse >= #TODO"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h  print this help message and exit"
    echo "-u  use safe mode: TRUE or FALSE [logical; default: FALSE]"
    echo "-d  run pipeline in \${TMPDIR}: TRUE or FALSE [logical; default:"
    echo "    TRUE]"
    echo "-m  initial memory allocation pool for JVM [chr; default: \"512m\"]"
    echo "-x  maximum memory allocation pool for JVM [chr; default: \"4096m\"]"
    echo "-r  string for \"sample #1\" [chr]"
    echo "-s  string for \"sample #2\" [chr]"
    echo "-1  bam infile #1, including path [chr]"
    echo "-2  bam infile #2, including path [chr]"
    echo "-p  prefix for outfiles [chr]"
    echo "-o  results directory for outfiles [chr]; path will be made if it"
    echo "    does not exist"
    echo "-b  number of records to read into memory at one time when running"
    echo "    the script for Part #1, get-AS-per-qname.R [int > 0; default:"
    echo "    100000]"
    echo "-c  number of records to read into memory at one time when running"
    echo "    the script for Part #3, generate-assignment-lists.R [int > 0;"
    echo "    default: 1000000]"
    echo "-t  alignment score threshold [int >= 0; default: 0]; the absolute"
    echo "    value of the difference in alignment scores between \"sample"
    echo "    #1\" and \"sample #2\" must be greater than this value in order"
    echo "    for a sample-specific assignment to be made; if not greater than"
    echo "    this value, then the assignment will be \"ambiguous\""
    echo "-a  count lines: TRUE or FALSE [logical; default: TRUE]"
    echo "-n  step in pipeline to run up to [int 1-4; default: 4]"
    exit
}


interactive=FALSE  # Hardcode TRUE for interactive testing; FALSE for CL usage
if [[ "${interactive}" == "TRUE" ]]; then
    safe_mode=FALSE
    use_TMPDIR=TRUE
    memory_min="512m"
    memory_max="4096m"
    strain_1="mm10"
    strain_2="CAST"
    bam_1="./data/files_bam_test/${strain_1}/Disteche_sample_1.dedup.${strain_1}.corrected.bam"
    bam_2="./data/files_bam_test/${strain_2}/Disteche_sample_1.dedup.${strain_2}.corrected.bam"
    prefix="Disteche_sample_1"
    outpath="./results/kga0/2022-0604_allelic-segregation"
    chunk_step_1=100000
    chunk_step_3=1000000
    threshold=0
    count=TRUE
    run_up_to=4
else
    while getopts "h:u:d:m:x:r:s:1:2:p:o:b:c:t:a:n:" opt; do
        case "${opt}" in
            h) print_usage ;;
            u) safe_mode="${OPTARG}" ;;
            d) use_TMPDIR="${OPTARG}" ;;
            m) memory_min="${OPTARG}" ;;
            x) memory_max="${OPTARG}" ;;
            r) strain_1="${OPTARG}" ;;
            s) strain_2="${OPTARG}" ;;
            1) bam_1="${OPTARG}" ;;
            2) bam_2="${OPTARG}" ;;
            p) prefix="${OPTARG}" ;;
            o) outpath="${OPTARG}" ;;
            b) chunk_step_1="${OPTARG}" ;;
            c) chunk_step_3="${OPTARG}" ;;
            t) threshold="${OPTARG}" ;;
            a) count="${OPTARG}" ;;
            n) run_up_to="${OPTARG}" ;;
            *) print_usage ;;
        esac
    done
fi

[[ -z "${safe_mode}" ]] && safe_mode=FALSE
[[ -z "${use_TMPDIR}" ]] && use_TMPDIR=TRUE
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


#TODO (   ) Some kind of check that zcat is GNU >= 1.09 and not, say, BSD
#  Check variable assignments from arguments ----------------------------------
echo -e ""
echo -e "Running ${0}... "

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-u: \"Safe mode\" is TRUE." && set -Eeuxo pipefail ;;
    false | f) echo -e "-u: \"Safe mode\" is FALSE." ;;
    *) \
        echo -e "Exiting: -u \"safe mode\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Check for necessary dependencies; exit if not found
check_dependency R
check_dependency samtools

#  Evaluate "${use_TMPDIR}"  #REMOVE
case "$(echo "${use_TMPDIR}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-d: \"Run in \${TMPDIR}\" is TRUE." && use_TMPDIR=TRUE ;;
    false | f) echo -e "-d: \"Run in \${TMPDIR}\" is FALSE." && use_TMPDIR=FALSE ;;
    *) \
        echo -e "Exiting: -d \"run in \${TMPDIR}\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Make sure "${strain_1}" and "${strain_2}" are not the same
[[ "${strain_1}" == "${strain_2}" ]] &&
    {
        echo -e "Exiting: -r \"strain 1\" is the same as -s \"strain 2\".\n"
        exit 1
    }

#  Evaluate "${bam_1}"
[[ -f "${bam_1}" ]] ||
    {
        echo -e "Exiting: -1 \"bam 1\" does not exist.\n"
        exit 1
    }

#  Evaluate "${bam_2}"
[[ -f "${bam_2}" ]] ||
    {
        echo -e "Exiting: -2 \"bam 2\" does not exist.\n"
        exit 1
    }

#  Make sure "${bam_1}" and "${bam_2}" are not the same
[[ "${bam_1}" == "${bam_2}" ]] &&
    {
        echo -e "Exiting: -1 \"bam 1\" is apparently the same as -2 \"bam 2\".\n"
        exit 1
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
        exit 1
    }

#  Evaluate "${count}"
case "$(echo "${count}" | tr '[:upper:]' '[:lower:]')" in
    true | t) echo -e "-c: \"Count lines\" is TRUE." ;;
    false | f) echo -e "-c: \"Count lines\" is FALSE." ;;
    *) \
        echo -e "Exiting: -c \"count lines\" argument must be TRUE or FALSE.\n"
        exit 1
        ;;
esac

#  Evaluate "${run_up_to}"
[[ ! "${run_up_to}" =~ ^[0-9]+$ ]] &&
    {
        echo -e "Exiting: -n \"run_up_to\" argument must be an integer in the range of 1-4.\n"
        exit 1
    }

case "${run_up_to}" in
    1 | 2 | 3 | 4) \
        echo -e "-n: \"Run up to\" is set to step #${run_up_to}."
        ;;
    *) \
        echo -e "Exiting: -n \"run_up_to\" argument must be an integer in the range of 1-4.\n"
        exit 1
        ;;
esac

echo -e ""


#  Adjust arguments if running "-d TRUE", then report the arguments -----------
[[ ${use_TMPDIR} == TRUE ]] &&
    {
        echo -e "Started: Copying bam infiles into \${TMPDIR}"
        base_bam_1="$(basename "${bam_1}")"
        base_bam_2="$(basename "${bam_2}")"

        if [[ $(basename "${bam_1}") == $(basename "${bam_2}") ]]; then
            echo -e "WARNING: bam #1 has the same basename as bam #2; to \
            differentiate the files, will include strain information in the \
            filename when copying bams into \${TMPDIR}"

            cp_bam_1="${TMPDIR}/${base_bam_1%.bam}.${RANDOM}.${strain_1}.bam"
            cp_bam_2="${TMPDIR}/${base_bam_2%.bam}.${RANDOM}.${strain_2}.bam"

            cp "${bam_1}" "${cp_bam_1}"
            cp "${bam_2}" "${cp_bam_2}"

            bam_1="${cp_bam_1}"
            bam_2="${cp_bam_2}"
        else
            cp_bam_1="${TMPDIR}/${base_bam_1%.bam}.${RANDOM}.bam"
            cp_bam_2="${TMPDIR}/${base_bam_2%.bam}.${RANDOM}.bam"

            cp "${bam_1}" "${cp_bam_1}"
            cp "${bam_2}" "${cp_bam_2}"
            
            bam_1="${cp_bam_1}"
            bam_2="${cp_bam_2}"
        fi
        echo -e "Completed: Copying bam infiles into \${TMPDIR}\n"
    }

echo -e "Running the pipeline with the following parameters:"
echo -e "  -u ${safe_mode}"
echo -e "  -d ${use_TMPDIR}"
echo -e "  -m ${memory_min}"
echo -e "  -x ${memory_max}"
echo -e "  -r ${strain_1}"
echo -e "  -s ${strain_2}"
echo -e "  -1 ${bam_1}"
echo -e "  -2 ${bam_2}"
echo -e "  -p ${prefix}"
echo -e "  -o ${outpath}"
echo -e "  -b ${chunk_step_1}"
echo -e "  -c ${chunk_step_3}"
echo -e "  -t ${threshold}"
echo -e "  -a ${count}"
echo -e "  -n ${run_up_to}\n\n"


#  Check that pipeline scripts are accessible by the driver -------------------
#+ ...if so, then assign them to variables
if [[ -f "./bin/workflow/get-AS-per-qname.R" ]]; then
    script_1="./bin/workflow/get-AS-per-qname.R"
    :
elif [[ -f "./get-AS-per-qname.R" ]]; then
    script_1="./get-AS-per-qname.R"
    :
else
    echo -e "Exiting: Couldn't find \"get-AS-per-qname.R\"."
    exit 1
fi

if [[ -f "./bin/workflow/find-set-intersection-set-complement.sh" ]]; then
    script_2="./bin/workflow/find-set-intersection-set-complement.sh"
    :
elif [[ -f "./find-set-intersection-set-complement.sh" ]]; then
    script_2="./find-set-intersection-set-complement.sh"
    :
else
    echo -e "Exiting: Couldn't find \"find-set-intersection-set-complement.sh\"."
    exit 1
fi

if [[ -f "./bin/workflow/generate-assignment-lists.R" ]]; then
    script_3="./bin/workflow/generate-assignment-lists.R"
    :
elif [[ -f "./generate-assignment-lists.R" ]]; then
    script_3="./generate-assignment-lists.R"
    :
else
    echo -e "Exiting: Couldn't find \"generate-assignment-lists.R\"."
    exit 1
fi

if [[ -f "./bin/workflow/filter-qnames-by-assignment.sh" ]]; then
    script_4="./bin/workflow/filter-qnames-by-assignment.sh"
    :
elif [[ -f "./filter-qnames-by-assignment.sh" ]]; then
    script_4="./filter-qnames-by-assignment.sh"
    :
else
    echo -e "Exiting: Couldn't find \"filter-qnames-by-assignment.sh\"."
    exit 1
fi

# echo "${script_1}"
# echo "${script_2}"
# echo "${script_3}"
# echo "${script_4}"


#  Assign variables for outdirectories ----------------------------------------
if [[ ${use_TMPDIR} == TRUE ]]; then
    dir_experiment_1="${TMPDIR}/${prefix}/01_run_get-AS-per-qname"
    dir_experiment_2="${TMPDIR}/${prefix}/02_run_find-set-inter-set-complement"
    dir_experiment_3="${TMPDIR}/${prefix}/03_run_generate-assignment-lists"
    dir_experiment_4="${TMPDIR}/${prefix}/04_run_filter-qnames-by-assignment"
else
    dir_experiment_1="${outpath}/${prefix}/01_run_get-AS-per-qname"
    dir_experiment_2="${outpath}/${prefix}/02_run_find-set-inter-set-complement"
    dir_experiment_3="${outpath}/${prefix}/03_run_generate-assignment-lists"
    dir_experiment_4="${outpath}/${prefix}/04_run_filter-qnames-by-assignment"
fi

# echo "${dir_experiment_1}"
# echo "${dir_experiment_2}"
# echo "${dir_experiment_3}"
# echo "${dir_experiment_4}"


#  Assign variables for completion files --------------------------------------
step_1="$(echo_completion_file_driver "${outpath}" "${prefix}" 1)"
step_2="$(echo_completion_file_driver "${outpath}" "${prefix}" 2)"
step_3="$(echo_completion_file_driver "${outpath}" "${prefix}" 3)"
step_4="$(echo_completion_file_driver "${outpath}" "${prefix}" 4)"

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
#TODO (   ) Provide more informative messages from "get-AS-per-qname.R"
#TODO (   ) Clean up the message formatting
if [[ ! -f "${step_1}" ]]; then
    echo -e "Started part 1/4: Running ${script_1}."
    
    #  Create key-value array for samples #1 and #2 (this style of array should
    #+ work with different versions of bash, zsh, etc.)
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
            exit 1
        fi
    done
fi

AS_1=$(find_files "${dir_experiment_1}" "*${strain_1}*AS.txt.gz")
AS_2=$(find_files "${dir_experiment_1}" "*${strain_2}*AS.txt.gz")

if [[ -f "${AS_1}" && -f "${AS_2}" ]]; then
    touch "${step_1}"
    echo -e "Completed part 1/4: Running ${script_1}.\n"
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
    echo -e "Started part 2/4: Running ${script_2}."
    
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
        exit 1
    fi
fi

comp_1=$(find_files "${dir_experiment_2}" "*${strain_1}*complement*")
comp_2=$(find_files "${dir_experiment_2}" "*${strain_2}*complement*")
inter=$(find_files "${dir_experiment_2}" "*${strain_1}-${strain_2}*inter*")

if [[ -f "${comp_1}" && -f "${comp_2}" && -f "${inter}" ]]; then
    touch "${step_2}"
    echo -e "Completed part 2/4: Running ${script_2}.\n"
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
    echo -e "Started part 3/4: Running ${script_3}."

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
        exit 1
    fi
fi

assign_1=$(find_files "${dir_experiment_3}" "*${strain_1}.txt.gz")
assign_2=$(find_files "${dir_experiment_3}" "*${strain_2}.txt.gz")
ambiguous=$(find_files "${dir_experiment_3}" "*ambiguous.txt.gz")

if [[ -f "${assign_1}" && -f "${assign_2}" && -f "${ambiguous}" ]]; then
    touch "${step_3}"
    echo -e "Completed part 3/4: Running ${script_3}.\n"
fi

evaluate_run_up_to "${run_up_to}" 3


#  Step 4 ---------------------------------------------------------------------
# -h print this help message and exit
# -s use safe mode: TRUE or FALSE (logical)
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
    echo -e "Started part 4/4: Running ${script_4}."

    if [[ -f "${assign_1}" && -f "${assign_2}" && "${ambiguous}" ]]; then
        echo -e \
        " \
        Running...\n \
        bash \"${script_4}\"\n \
        -s \"${safe_mode}\"\n \
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
        exit 1
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
    echo -e "Completed part 4/4: Running ${script_4}.\n"
fi

[[ ${use_TMPDIR} == TRUE ]] &&
    {
        if \
            [[ -f "${step_1}" && -f "${step_2}" ]] && \
            [[ -f "${step_3}" && -f "${step_4}" ]]; \
        then
            mv -f "${TMPDIR}/${prefix}" "${outpath}" &&
            rm "${bam_1}" "${bam_2}" "${bam_1}.bai" "${bam_2}.bai" &&
            echo -e "Copied outfiles from \${TMPDIR} to ${outpath}"
        fi
    }


#  Return run time ------------------------------------------------------------
time_end="$(date +%s)"
calculate_run_time "${time_start}" "${time_end}" "Completed: ${0}"
echo -e ""

exit 0

#!/bin/bash

#  driver_allelic-segregation.sh
#  KA


#  Set directories of interest, variables -------------------------------------

#TODO (   ) Add logic to see if the user can access this from "${dir_base}"
#           and/or if the scripts are in "${dir_base}"
dir_code="./bin/workflow"
dir_data="./data/files_bam"

cd "${dir_base}" || echo "cd failed. Look into this."


#  Set up variables -----------------------------------------------------------
date=$(date '+%Y-%m%d')

ID=1
# ID="${1}"

strain_1="mm10"
# strain_1="${2}"

strain_2="CAST"
# strain_2="${3}"

memory_min="512m"
# memory_min="${4}"

memory_max="4096m"
# memory_max="${5}"

#TODO (   ) Let the user determine the base directory, which either contains
#           the scripts that are called can access them in ./bin/workflow
dir_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross"
# dir_base="${6}"

#TODO (   ) Let the user determine the outdirectory, then make new directories
#TODO       within the outdirectory
dir_results="./results/kga0"
# dir_results="${7}"


dir_experiment_1="${dir_data}"
dir_experiment_2="${dir_results}/${date}_run_get-AS-per-qname"
dir_experiment_3="${dir_results}/${date}_run_find-set-intersection-set-complement"
dir_experiment_4="${dir_results}/${date}_run_generate-assignment-lists"
dir_experiment_5="${dir_results}/${date}_run_filter-qnames-by-assignment"


#  Create outdirectories ------------------------------------------------------
mkdir -p "${dir_experiment_2}"
mkdir -p "${dir_experiment_3}"
mkdir -p "${dir_experiment_4}"


#  Step 1 ---------------------------------------------------------------------
for strain in "CAST" "mm10"; do
    samtools index \
    "${dir_experiment_1}/Disteche_sample_${ID}.dedup.${strain}.corrected.bam"
done


#  Step 2 ---------------------------------------------------------------------
#TODO (   ) Provide more informative messages output by get-AS-per-qname.R
#TODO (   ) Clean up the message formatting
#TODO (   ) Fix the name of outfiles, which mention the strain twice
for strain in "${strain_1}" "${strain_2}"; do
    if [[ -f "${dir_data}/Disteche_sample_${ID}.dedup.${strain}.corrected.bam.bai" ]]; then
        Rscript "${dir_code}/get-AS-per-qname.R" \
        --bam "${dir_data}/Disteche_sample_${ID}.dedup.${strain}.corrected.bam" \
        --bai "${dir_data}/Disteche_sample_${ID}.dedup.${strain}.corrected.bam.bai" \
        --outdir "${dir_experiment_2}" \
        --strain "${strain}" \
        --chunk 100000
    else
        echo "Exiting: ${dir_data}/Disteche_sample_${ID}.dedup.${strain}.corrected.bam.bai not found."
        exit 1
    fi
done
#   -b, --bam     bam infile, including path <chr>
#   -i, --bai     bam index, including path <chr>
#   -o, --outdir  directory for saving rds outfile, including path <chr>
#   -s, --strain  strain name to be appended to rds outfile columns <chr>
#   -c, --chunk   number of records to read into memory at a single time
#                 <even int> [default: 100000]


#  Step 3 ---------------------------------------------------------------------
#TODO (   ) Script uses display_spinning_icon(); need to remove this
if [[ -f "${dir_experiment_2}/Disteche_sample_${ID}.dedup.${strain_1}.corrected.${strain_1}.AS.txt.gz" ]]; then
    bash "${dir_code}/find-set-intersection-set-complement.sh" \
    -u FALSE \
    -i "${dir_experiment_2}/Disteche_sample_${ID}.dedup.${strain_1}.corrected.${strain_1}.AS.txt.gz" \
    -j "${dir_experiment_2}/Disteche_sample_${ID}.dedup.${strain_2}.corrected.${strain_2}.AS.txt.gz" \
    -1 "${strain_1}" \
    -2 "${strain_2}" \
    -o "${dir_experiment_3}" \
    -p "Disteche_sample_${ID}.dedup" \
    -c TRUE
else
    echo "Exiting: ${dir_experiment_2}/Disteche_sample_${ID}.dedup.${strain_1}.corrected.AS.txt.gz not found."
    # exit 1
fi
# -h print this help message and exit
# -u use safe mode: TRUE or FALSE (logical; default: FALSE)
# -i AS.txt.gz "infile #1", including path (chr)
# -j AS.txt.gz "infile #2", including path (chr)
# -1 string for "sample #1" (chr)
# -2 string for "sample #2" (chr)
# -o path for outfiles (chr); path will be made if it does not exist
# -p prefix for outfiles (chr)
# -c count lines: TRUE or FALSE (logical)


#  Step 4 ---------------------------------------------------------------------
#TODO (   ) Provide more informative messages output by generate-assignment-lists.R
#TODO (   ) Clean up the message formatting
if [[ -f "${dir_experiment_3}/Disteche_sample_${ID}.dedup.${strain_1}-${strain_2}.intersection.txt.gz" ]]; then
    Rscript "${dir_code}/generate-assignment-lists.R" \
    --intersection "${dir_experiment_3}/Disteche_sample_${ID}.dedup.${strain_1}-${strain_2}.intersection.txt.gz" \
    --sample_1 "${strain_1}" \
    --sample_2 "${strain_2}" \
    --outdir "${dir_experiment_4}" \
    --outprefix "Disteche_sample_${ID}.dedup" \
    --threshold 0 \
    --chunk 1000000 \
    --remove TRUE
else
    echo "Exiting: ${dir_experiment_3}/Disteche_sample_${ID}.${strain_1}-${strain_2}.intersection.txt.gz not found."
    # exit 1
fi
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
#                       in column 2 of the tab- separated intersection
#                       file (chr)
#   -u, --sample_2      sample name for sample 2; corresponds to the AS
#                       in column 3 of the tab- separated intersection
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


#  Step 5 ---------------------------------------------------------------------
if [[ -f "${dir_experiment_4}/Disteche_sample_${ID}.dedup.ambiguous.txt.gz" ]]; then
    bash ./bin/workflow/filter-qnames-by-assignment.sh \
    -s FALSE \
    -l FALSE \
    -m "${memory_min}" \
    -x "${memory_max}" \
    -b "${dir_data}/Disteche_sample_${ID}.dedup.${strain_1}.corrected.bam" \
    -c "${dir_data}/Disteche_sample_${ID}.dedup.${strain_2}.corrected.bam" \
    -a "${dir_experiment_4}/Disteche_sample_${ID}.dedup.ambiguous.txt.gz" \
    -i "${dir_experiment_4}/Disteche_sample_${ID}.dedup.${strain_1}.txt.gz" \
    -j "${dir_experiment_4}/Disteche_sample_${ID}.dedup.${strain_2}.txt.gz" \
    -u "${dir_experiment_3}/Disteche_sample_${ID}.dedup.${strain_1}.complement.txt.gz" \
    -v "${dir_experiment_3}/Disteche_sample_${ID}.dedup.${strain_2}.complement.txt.gz" \
    -1 "${strain_1}" \
    -2 "${strain_2}" \
    -o "${dir_experiment_5}" \
    -p "Disteche_sample_${ID}.dedup" \
    -n TRUE
else
    echo "Exiting: ./results/kga0/2022-0523_run_generate-assignment-lists/Disteche_sample_1.dedup.ambiguous.txt.gz not found."
    # exit 1
fi
# -h print this help message and exit
# -s use safe mode: "TRUE" or "FALSE" (logical)
# -l run on GS HPC: "TRUE" or "FALSE" (logical)
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
# -n count lines: "TRUE" or "FALSE" (logical)
# #TODO -f run samtools flagstat: "TRUE" or "FALSE" (logical)

#!/bin/bash

#  tally-assignments.sh
#  KA

qlogin -l mfree=2G
4dn-home

dir_log="/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/log"
dir_log_out="${dir_log}/$(date '+%Y-%m%d')_tally-assignments"

mkdir -p "${dir_log_out}"

dir_results="/net/noble/vol1/home/gangliuw/proj/2022-01-mouse-cross/src/2022-06-06-ALSO/results"  # "${1}"
cd "${dir_results}" || echo "cd failed. Check on this."

s_1="mm10"  # "${2}"
s_2="CAST"  # "${3}"
file_out="$(date '+%Y-%m%d')_tally-assignments.txt"  # "${4}"

for i in {1..22}; do
    exp="Disteche_sample_${i}_ALSO"  #  Get this into an argument
    step_2="${exp}/02_run_find-set-inter-set-complement"
    step_3="${exp}/03_run_generate-assignment-lists"

    {
        echo "Working with experiment ${exp}"
        echo "Files and directories within ${exp}"
        ls -lhaFG "${exp}" && echo ""
        ls -lhaFG "${step_2}" && echo ""
        ls -lhaFG "${step_3}" && echo ""
        echo "" && echo ""
    } \
    >> "${dir_log_out}/${file_out%.txt}.ls.txt"
done

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
"exp" "sample_1" "sample_2" "ambiguous" "total" \
"pct_sample_1" "pct_sample_2" "pct_ambiguous" \
> "${dir_log_out}/${file_out}"
# vi "${dir_log_out}/${file_out}"

for i in {1..22}; do
    exp="Disteche_sample_${i}_ALSO"  #  Get this into an argument
    step_2="${exp}/02_run_find-set-inter-set-complement"
    step_3="${exp}/03_run_generate-assignment-lists"

    echo ""
    echo "Working with experiment ${exp}"
    # echo "Files and directories within ${exp}"
    # ls -lhaFG "${exp}"

    echo "- Tallying \$comp_s_1..."
    comp_s_1=$(zcat "${step_2}/"*".${s_1}.complement.txt.gz" | wc -l)

    echo "- Tallying \$comp_s_2..."
    comp_s_2=$(zcat "${step_2}/"*".${s_2}.complement.txt.gz" | wc -l)
    
    echo "- Tallying \$assign_fr_inter_s_1..."
    assign_fr_inter_s_1=$(zcat "${step_3}/"*".${s_1}.txt.gz" | wc -l)
    
    echo "- Tallying \$assign_fr_inter_s_2..."
    assign_fr_inter_s_2=$(zcat "${step_3}/"*".${s_2}.txt.gz" | wc -l)
    
    echo "- Tallying \$assign_fr_inter_a..."
    assign_fr_inter_a=$(zcat "${step_3}/"*".ambiguous.txt.gz" | wc -l)

    echo "- Writing lines to outfile..."
    final_s_1=$(( comp_s_1 + assign_fr_inter_s_1 ))
    final_s_2=$(( comp_s_2 + assign_fr_inter_s_2 ))
    final_a=$(( assign_fr_inter_a ))
    total=$(( final_s_1 + final_s_2 + final_a ))
    pct_s_1="$(echo "scale=2; ${final_s_1}/${total}" | bc)"
    pct_s_2="$(echo "scale=2; ${final_s_2}/${total}" | bc)"
    pct_a="$(echo "scale=2; ${final_a}/${total}" | bc)"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${exp}" "${final_s_1}" "${final_s_2}" "${final_a}" "${total}" \
    "${pct_s_1}" "${pct_s_2}" "${pct_a}" \
    >> "${dir_log_out}/${file_out}"
done

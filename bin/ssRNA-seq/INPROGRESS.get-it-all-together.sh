#!/bin/bash

#  INPROGRESS.get-it-all-together.sh
#  KA

#  Dependencies
#+ - samtools/1.14
#+ - parallel/20200922
#+ - rename (util-linux package)


#  Removing the tally numbers...
unset infiles
typeset -a infiles
infiles=(
    spleen_rep_2.SNPsplit.ambiguous.srt.txt.gz
    spleen_rep_2.SNPsplit.mm10.srt.txt.gz
    spleen_rep_2.SNPsplit.SPRET.srt.txt.gz
)
for i in "${infiles[@]}"; do
    cut -d" " -f2- <(zcat -df "${i}") | gzip > "tmp.txt.gz"
    mv -f "tmp.txt.gz" "${i}"
done


#  Output elements in set 1 that are not in set 2
#+ Set 1: spleenrep2Aligned.out.sorted.primary.mapped.rm.QNAME.srt.txt.gz
#+ Set 2: spleenrep2Aligned.out.genome1.genome2.primary.mapped.rm.QNAME.srt.txt.gz

set1="spleenrep2Aligned.out.sorted.primary.mapped.rm.QNAME.srt.txt.gz"
set2="spleenrep2Aligned.out.genome1.genome2.primary.mapped.rm.QNAME.srt.txt.gz"

# comm -23 <(sort set1) <(sort set2)
time comm -23 \
<(zcat -df "${set1}" | sort) \
<(zcat -df "${set2}" | sort) \
    > SNPsplit.ambiguous.test-comm.txt
#  Unsorted "${set2}"
# real    5m39.862s
# user    7m16.116s
# sys 0m5.039s

# grep -vxF -f set2 set1
time grep -vxF -f \
<(zcat -df "${set1}" | sort) \
<(zcat -df "${set2}" | sort) \
    > SNPsplit.ambiguous.test-grep.txt
#  Unsorted "${set2}"
# Killed  # Exit 137 (>4GB RAM in use)
#
# real    4m54.968s
# user    4m45.910s
# sys 0m5.588s

parsort_file_qname_auto "${set2}"
set2="spleenrep2Aligned.out.genome1.genome2.primary.mapped.rm.QNAME.srt.srt.txt.gz"

# comm -23 <(sort set1) <(sort set2)
time comm -23 \
<(zcat -df "${set1}" | sort) \
<(zcat -df "${set2}" | sort) \
    > SNPsplit.ambiguous.test-comm-sorted.txt
#  Sorted "${set2}"
# real    5m18.670s
# user    6m51.347s
# sys 0m4.835s

time grep -vxF -f \
<(zcat -df "${set1}" | sort) \
<(zcat -df "${set2}" | sort) \
    > SNPsplit.ambiguous.test-grep-sorted.txt
#  Sorted "${set2}"
# Killed  # Exit 137 (>4GB RAM in use)
#
# real    4m52.209s
# user    4m44.408s
# sys 0m5.502s


#  Comm seems to be better, and sorting prior to set work speeds things up a little
cp SNPsplit.ambiguous.test-comm-sorted.txt spleenrep2Aligned.out.ambiguous.primary.mapped.rm.QNAME.srt.txt.gz


#  Need to find set intersections for nonredundant pairwise combinations
cp spleenrep2Aligned.out.ambiguous.primary.mapped.rm.QNAME.srt.txt.gz spleen_rep_2.SNPsplit.ambiguous.srt.txt.gz
cp spleenrep2Aligned.out.genome1.primary.mapped.rm.QNAME.srt.txt.gz spleen_rep_2.SNPsplit.mm10.srt.txt.gz
cp spleenrep2Aligned.out.genome2.primary.mapped.rm.QNAME.srt.txt.gz spleen_rep_2.SNPsplit.SPRET.srt.txt.gz
cp spleen_rep_2.SRR1525411.ambiguous.srt.txt.gz spleen_rep_2.ALSO.ambiguous.srt.txt.gz
cp spleen_rep_2.SRR1525411.mm11.combined.srt.txt.gz spleen_rep_2.ALSO.mm11.srt.txt.gz
cp spleen_rep_2.SRR1525411.SPRET.combined.srt.txt.gz spleen_rep_2.ALSO.SPRET.srt.txt.gz


#  Generating the confusion matrix by, first, getting appropriate counts
unset infiles
typeset -a infiles
while IFS=" " read -r -d $'\0'; do
    infiles+=( "${REPLY}" )
done < <(find . -maxdepth 1 -type f -name "*.srt.txt.gz" -print0 | sort -z)
for i in "${infiles[@]}"; do echo "${i}"; done

unset set_1 && typeset -a set_1
unset set_2 && typeset -a set_2
for ((i = 0; i < ${#infiles[@]}; i++)); do
    for ((j = i + 1; j < ${#infiles[@]}; j++)); do
        # echo -e "${infiles[i]}\t${infiles[j]}"
        set_1+=("${infiles[i]}")
        set_2+=("${infiles[j]}")
    done
done
for i in "${set_1[@]}"; do echo "${i}"; done
for i in "${set_2[@]}"; do echo "${i}"; done

#  Remove unneeded elements
for i in 0 1 5 12 13 14; do
    unset "set_1[${i}]"
    unset "set_2[${i}]"
done
for i in "${set_1[@]}"; do echo "${i}"; done
for i in "${set_2[@]}"; do echo "${i}"; done
echo "${#set_1[@]}"  # 9
echo "${#set_2[@]}"  # 9
echo "${!set_1[@]}"  # 2 3 4 6 7 8 9 10 11
echo "${!set_2[@]}"  # 2 3 4 6 7 8 9 10 11

#  Don't think you can directly rename bash associative array keys, so just
#+ rebuild them as new arrays
#+
#+ - stackoverflow.com/questions/12303974/copying-a-bash-array-fails
#+ - stackoverflow.com/questions/16860877/remove-an-element-from-a-bash-array
unset set_1_tmp
unset set_2_tmp
typeset -a set_1_tmp=( "${set_1[@]}" )
typeset -a set_2_tmp=( "${set_2[@]}" )
for i in "${set_1_tmp[@]}"; do echo "${i}"; done
for i in "${set_2_tmp[@]}"; do echo "${i}"; done

unset set_1 && typeset -a set_1
unset set_2 && typeset -a set_2
for i in "${!set_1_tmp[@]}"; do
    set_1+=( "${set_1_tmp[i]}" )
    set_2+=( "${set_2_tmp[i]}" )
done
for i in "${set_1[@]}"; do echo "${i}"; done
for i in "${set_2[@]}"; do echo "${i}"; done
echo "${#set_1[@]}"  # 9
echo "${#set_2[@]}"  # 9
echo "${!set_1[@]}"  # 0 1 2 3 4 5 6 7 8
echo "${!set_2[@]}"  # 0 1 2 3 4 5 6 7 8

unset set_1_tmp
unset set_2_tmp

#  stackoverflow.com/questions/12204192/using-multiple-delimiters-in-awk
unset intersection
typeset -a intersection
for i in $(seq 0 $(( ${#set_1[@]} -  1 ))); do
    # echo "${set_1[i]}"
    # echo "${set_2[i]}"

    p_1=$(echo "${set_1[i]}" | awk -F'\\/|\\.' '{ print $3 }')
    p_2=$(echo "${set_1[i]}" | awk -F'\\/|\\.' '{ print $4"-"$5 }')
    p_3=$(echo "${set_2[i]}" | awk -F'\\/|\\.' '{ print $4"-"$5 }')
    p_4=$(echo "${set_2[i]}" | awk -F'\\/|\\.' '{ print $7"."$8 }')

    echo "${p_1}.intersection.${p_2}.${p_3}.${p_4}"
    intersection+=( "${p_1}.intersection.${p_2}.${p_3}.${p_4}" )
done
for i in "${intersection[@]}"; do echo "${i}"; done
echo "${#intersection[@]}"  # 9
echo "${!intersection[@]}"  # 0 1 2 3 4 5 6 7 8


parallelize=4
parallel -k -j "${parallelize}" \
"comm -12 \
<(zcat -df {1} | sort) \
<(zcat -df {2} | sort) \
    | gzip \
    > {3}" \
::: "${set_1[@]}" \
:::+ "${set_2[@]}" \
:::+ "${intersection[@]}"

for i in *intersection*; do
    echo "${i}"
    zcat "${i}" | wc -l
    echo ""
done
# spleen_rep_2.intersection.ALSO-ambiguous.SNPsplit-ambiguous.txt.gz
# 18005061
#
# spleen_rep_2.intersection.ALSO-ambiguous.SNPsplit-mm10.txt.gz
# 9753
#
# spleen_rep_2.intersection.ALSO-ambiguous.SNPsplit-SPRET.txt.gz
# 6161
#
# spleen_rep_2.intersection.ALSO-mm11.SNPsplit-ambiguous.txt.gz
# 3020069
#
# spleen_rep_2.intersection.ALSO-mm11.SNPsplit-mm10.txt.gz
# 2790807
#
# spleen_rep_2.intersection.ALSO-mm11.SNPsplit-SPRET.txt.gz
# 64419
#
# spleen_rep_2.intersection.ALSO-SPRET.SNPsplit-ambiguous.txt.gz
# 1243674
#
# spleen_rep_2.intersection.ALSO-SPRET.SNPsplit-mm10.txt.gz
# 30433
#
# spleen_rep_2.intersection.ALSO-SPRET.SNPsplit-SPRET.txt.gz
# 2534357

unset types
typeset -a types=(
    ALSO-ambiguous
    ALSO-mm11
    ALSO-SPRET
    SNPsplit-ambiguous
    SNPsplit-mm10
    SNPsplit-SPRET
)
for i in "${types[@]}"; do echo "${i}"; done

for i in "${types[@]}"; do
    # i="${types[2]}"
    # echo "${i}"
   
    to_test=$(echo "${i}" | awk -F'-' '{ print $1 }')
    # echo "\$\"{to_test}\" is ${to_test}"
   
    if [[ "${to_test}" == "ALSO" ]]; then
        j="SNPsplit"
    elif [[ "${to_test}" == "SNPsplit" ]]; then
        j="ALSO"
    fi
    # echo "\$\"{j}\" is ${j}"

    outfile="spleen_rep_2.non-unique.${i}-in-${j}.txt.gz"
    # echo "${outfile}"
   
    # echo -- *"${i}"*
    cat -- *"${i}"* > "${outfile}"

    echo "${i}"
    zcat "${outfile}" | wc -l
    echo ""
done

#  Now, need to identify and list elements in "${infiles[@]}" that are not in
#+ associated *non-unique* files, i.e., "unique elements"; for example,
#+ elements from ALSO-related alignment not found amongst elements from
#+ SNPsplit-related alignment, and vice versa

unset non_unique
typeset -a non_unique
while IFS=" " read -r -d $'\0'; do
    non_unique+=( "${REPLY}" )
done < <(find . -maxdepth 1 -type f -name "*.non-unique.*" -print0 | sort -z)
for i in "${non_unique[@]}"; do echo "${i}"; done

for i in "${!non_unique[@]}"; do
    # i=2
    # echo "${infiles[i]}"
    # echo "${non_unique[i]}"
    # echo ""

    p_1=$(echo "${non_unique[i]}" | awk -F'\\/|\\.|-' '{ print $3 }')
    echo "${p_1}"
   
    p_2=$(echo "${non_unique[i]}" | awk -F'\\/|\\.|-' '{ print $5 }')
    echo "${p_2}"

    p_3=$(echo "${non_unique[i]}" | awk -F'\\/|\\.|-' '{ print $6"-"$7"-not-"$8"-"$9 }')
    echo "${p_3}"

    p_4=$(echo "${non_unique[i]}" | awk -F'\\/|\\.|-' '{ print $10"."$11 }')
    echo "${p_4}"

    echo "${p_1}.${p_2}.${p_3}.${p_4}"
    unique+=( "${p_1}.${p_2}.${p_3}.${p_4}" )
done
for i in "${unique[@]}"; do echo "${i}"; done


for i in "${!unique[@]}"; do
    echo "#$(( i + 1 ))"
    echo "   \${infiles[i]}  ${infiles[i]}"
    echo "\${non_unique[i]}  ${non_unique[i]}"
    echo "    \${unique[i]}  ${unique[i]}"
    echo ""
done


# comm -23 <(sort set1) <(sort set2)
parallel -k -j "${parallelize}" \
"comm -23 \
<(zcat -df {1} | sort) \
<(zcat -df {2} | sort) \
    | gzip \
    > {3}" \
::: "${infiles[@]}" \
:::+ "${non_unique[@]}" \
:::+ "${unique[@]}"

for i in "${unique[@]}"; do ls "${i}"; done
for i in "${unique[@]}"; do
    echo "${i}"
    zcat "${i}" | wc -l
    echo ""
done


#  2022-0808
#  Collecting stats on ALSO and SNPsplit assignments
cd /net/noble/vol6/user/kga0/2021_kga0_4dn-mouse-cross/data/Berletch_Fang/alignments_primary \
    || echo "cd failed."
unset bams
typeset -a bams
while IFS=" " read -r -d $'\0'; do
    bams+=( "${REPLY}" )
done < <(find . -type f -name "*.bam" -print0 | sort -z)
for i in "${bams[@]}"; do echo "${i}"; done

parallel -k -j "${parallelize}" \
"echo {1}; \
samtools view -c {1}; \
echo ''" \
::: "${bams[@]}"


cd /net/noble/vol6/user/kga0/2021_kga0_4dn-mouse-cross/data/Berletch_Fang/alignments_primary/mapped/MAPQ-gte-30 \
    || echo "cd failed."
for i in *.bam; do
    samtools flagstat "${i}" > "${i%.bam}.flagstat.txt"
done

mv -- *.flagstat.txt flagstat/


cd /net/noble/vol6/user/kga0/2021_kga0_4dn-mouse-cross/data/Berletch_Fang \
    || echo "cd failed."
unset bams
typeset -a bams
while IFS=" " read -r -d $'\0'; do
    bams+=( "${REPLY}" )
done < <(find . -maxdepth 1 -type f -name "*.bam" -print0 | sort -z)
for i in "${bams[@]}"; do echo "${i}"; done

parallel -k -j "${parallelize}" \
"echo {1}; \
samtools view -c {1}; \
echo ''" \
::: "${bams[@]}"


cd /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/data \
    || echo "cd failed."
unset bams_ALSO
typeset -a bams_ALSO
while IFS=" " read -r -d $'\0'; do
    bams_ALSO+=( "${REPLY}" )
done < <(find . -type f -name "*.bam" -not \( -path "./Bonora*" -o -path "./files*" -o -path "./*/individual/*" -prune \) -print0 | sort -z)
for i in "${bams_ALSO[@]}"; do echo "${i}"; done

parallel -k -j "${parallelize}" \
"echo {1}; \
samtools view -c {1}; \
echo ''" \
::: "${bams_ALSO[@]}"


cd /net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results \
    || echo "cd failed."
unset bams_ALSO
typeset -a bams_ALSO
while IFS=" " read -r -d $'\0'; do
    bams_ALSO+=( "${REPLY}" )
done < <(find ./*ALSO* -type f -name "*.bam" -print0 | sort -z)
for i in "${bams_ALSO[@]}"; do echo "${i}"; done

parallel -k -j "${parallelize}" \
"echo {1}; \
samtools view -c {1}; \
echo ''" \
::: "${bams_ALSO[@]}"

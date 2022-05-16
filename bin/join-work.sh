#!/bin/bash


# zcat Disteche_sample_1.dedup.CAST.corrected.CAST.AS.txt.gz \
# | head -1000 \
# | gzip \
# > test-CAST-1-1000.txt.gz

# zcat Disteche_sample_1.dedup.mm10.corrected.mm10.AS.txt.gz \
# | head -1000 \
# | gzip \
# > test-mm10-1-1000.txt.gz


#  Subset data to 100000 ------------------------------------------------------
start="$(date +%s)"
dir_in="./data/files_bam"
dir_out="./data/files_bam"
prefix="Disteche_sample_1.dedup"
suffix="AS.txt.gz"
mm10="${dir_in}/${prefix}.mm10.corrected.mm10.${suffix}"
CAST="${dir_in}/${prefix}.CAST.corrected.CAST.${suffix}"
intersection="${dir_out}/${prefix}.intersection.${suffix}"
unique_mm10="${dir_out}/${prefix}.unique-mm10.${suffix}"
unique_CAST="${dir_out}/${prefix}.unique-CAST.${suffix}"

n="100000"
mm10_100K="${mm10/.${suffix}/.${n}.${suffix}}"
CAST_100K="${CAST/.${suffix}/.${n}.${suffix}}"
intersection_100K="${intersection/.${suffix}/.sampled.${suffix}}"
unique_mm10_100K="${unique_mm10/.${suffix}/.sampled.${suffix}}"
unique_CAST_100K="${unique_CAST/.${suffix}/.sampled.${suffix}}"

randomly_sample_lines_from_file "${n}" "${mm10}" "${mm10_100K}"
randomly_sample_lines_from_file "${n}" "${CAST}" "${CAST_100K}"

check_file_sorted "${mm10_100K}"
check_file_sorted "${CAST_100K}"

head_10 "${mm10_100K}" && echo ""
head_10 "${CAST_100K}"

find_set_intersection "${mm10_100K}" "${CAST_100K}" "${intersection_100K}"
head_10 "${intersection_100K}"
# AACGTACCATCTAAGAGTTATTAGCGGAACCGTTACGTTG:1352162525 -5  0
# AACTTACGCTATGCTAACCTTACCTTAGCTCGGCGCTCAA:212010352  -3  -8
# AAGATATAGATCAGCGCAATAATATAGCGGTTATAATCAA:770876836  0   -5
# AAGGACCTTATTGACTTCGTTCTGCGTCGCCGTTACGTTG:1265379755 0   0
# AAGTCTTCTGCGTAGCCGTACCAATTCCATAAGTTCGCTG:914643066  -11 -11
# AAGTTCGTAATCAGCGCAATCCAATTCCATCGTTACGTTG:2446700132 0   0
# AATTACTAATAGCTCATTCGGCGCCTAGTTCGAGATAAGA:750992127  0   0
# AATTACTAATTGCGAATCGGACTCGGCATAGGAGTCGTCT:2470826742 0   -5
# AATTACTAATTGCGAATCGGACTCGGCATAGGAGTCGTCT:929365282  0   0
# ACCAATATGGTAACCAATCGTATCAAGGTCATCGGCATTG:1578728787 -10 0
zcat "${intersection_100K}" | wc -l  # 240
check_file_sorted "${intersection_100K}"

#  Output unique elements
find_set_complement "${mm10_100K}" "${CAST_100K}" "${unique_mm10_100K}"
check_file_sorted "${unique_mm10_100K}"

find_set_complement "${CAST_100K}" "${mm10_100K}" "${unique_CAST_100K}"
check_file_sorted "${unique_CAST_100K}"

head_10 "${unique_mm10_100K}" && echo ""
head_10 "${unique_CAST_100K}"

count_lines_gzip "${unique_mm10_100K}"  # 99760
count_lines_gzip "${unique_CAST_100K}"  # 99760
count_lines_gzip "${intersection_100K}"  # 240

echo $(( 99760 + 240 ))  # 100000


#  Subset to 100000 and 125000 ------------------------------------------------
n="100000"
m="125000"
mm10_sample="${mm10/.${suffix}/.${n}.${suffix}}"
CAST_sample="${CAST/.${suffix}/.${m}.${suffix}}"
intersection_sample="${intersection/.${suffix}/.sampled.${suffix}}"
unique_mm10_sample="${unique_mm10/.${suffix}/.sampled.${suffix}}"
unique_CAST_sample="${unique_CAST/.${suffix}/.sampled.${suffix}}"

randomly_sample_lines_from_file "${n}" "${mm10}" "${mm10_sample}"
randomly_sample_lines_from_file "${m}" "${CAST}" "${CAST_sample}"

check_file_sorted "${mm10_sample}"
check_file_sorted "${CAST_sample}"

head_10 "${mm10_sample}" && echo ""
head_10 "${CAST_sample}"

find_set_intersection "${mm10_sample}" "${CAST_sample}" "${intersection_sample}"
head_10 "${intersection_sample}"
# AACTATTGGCTAGAACCTCCAGTCGCGTCGCTCTGCGATC:2295386972 0   0
# AACTATTGGCTTGACTTCGTACGTACTAGGAAGTTCGCTG:1839539119 0   -5
# AACTGCTAAGGACTCGAGGAATTGAGGATATTATAATCAA:65283900   0   0
# AACTTACGCTATGCTAACCTTACCTTAGCTCGGCGCTCAA:370333393  -3  -3
# AACTTACGCTATGCTAACCTTACCTTAGCTCGGCGCTCAA:911638007  0   0
# AAGATATAGATCAGCGCAATAATATAGCGGTTATAATCAA:379934231  0   -4
# AATTAACCATTCATCAGCCACCGGTCCTAAAAGTTCGCTG:2474630516 0   0
# AATTACTAATCGAGTTCGCCGCGGTTGGAAAAGTAGACTA:1657963169 -3  -16
# AATTACTAATTGCGAATCGGCTGATCAGAGCGAGATAAGA:395755511  0   0
# AATTGCCGTCCGAGTTCGCCCCTTCATCCGTGATTCTCGT:172815629  0   -5
zcat "${intersection_sample}" | wc -l  # 294
check_file_sorted "${intersection_sample}"

#  Output unique elements
find_set_complement "${mm10_sample}" "${CAST_sample}" "${unique_mm10_sample}"
check_file_sorted "${unique_mm10_sample}"

find_set_complement "${CAST_sample}" "${mm10_sample}" "${unique_CAST_sample}"
check_file_sorted "${unique_CAST_sample}"

head_10 "${unique_mm10_sample}" && echo ""
head_10 "${unique_CAST_sample}"

count_lines_gzip "${unique_mm10_sample}"  # 99706
count_lines_gzip "${unique_CAST_sample}"  # 124706
count_lines_gzip "${intersection_sample}"  # 294

echo $(( 99706 + 294 ))  # 100000
echo $(( 124706 + 294 ))  # 125000
#NOTE So, these sums make sense, but the sums for the full files do not...


#  Subset to 900000 and 1000000 -----------------------------------------------
n="900000"
m="1000000"
mm10_sample="${mm10/.${suffix}/.${n}.${suffix}}"
CAST_sample="${CAST/.${suffix}/.${m}.${suffix}}"
intersection_sample="${intersection/.${suffix}/.sampled.${suffix}}"
unique_mm10_sample="${unique_mm10/.${suffix}/.sampled.${suffix}}"
unique_CAST_sample="${unique_CAST/.${suffix}/.sampled.${suffix}}"

randomly_sample_lines_from_file "${n}" "${mm10}" "${mm10_sample}"
randomly_sample_lines_from_file "${m}" "${CAST}" "${CAST_sample}"

check_file_sorted "${mm10_sample}"
check_file_sorted "${CAST_sample}"

head_10 "${mm10_sample}" && echo ""
head_10 "${CAST_sample}"

find_set_intersection "${mm10_sample}" "${CAST_sample}" "${intersection_sample}"
head_10 "${intersection_sample}"
# AACCAATAACAAGTCAACGGTATCAAGGTCCTCTGCGATC:1395717943 -5  0
# AACCAATAACAAGTCAACGGTATCAAGGTCCTCTGCGATC:1512701808 0   -5
# AACCAATAACAAGTCAACGGTATCAAGGTCCTCTGCGATC:2160241733 0   0
# AACCAATAACCGTAGCCGTAACGTACTAGGGGAGTCGTCT:2353519801 -5  0
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:114051865  0   0
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:1288156721 0   -5
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:1859554    0   0
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:188568440  0   0
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:2079772555 0   0
# AACCAATAACGAATGCCTTGGATAGCCGATGCATAGTATA:2083912767 0   0
zcat "${intersection_sample}" | wc -l  # 20648
check_file_sorted "${intersection_sample}"

#  Output unique elements
find_set_complement "${mm10_sample}" "${CAST_sample}" "${unique_mm10_sample}"
check_file_sorted "${unique_mm10_sample}"

find_set_complement "${CAST_sample}" "${mm10_sample}" "${unique_CAST_sample}"
check_file_sorted "${unique_CAST_sample}"

head_10 "${unique_mm10_sample}" && echo ""
head_10 "${unique_CAST_sample}"

count_lines_gzip "${unique_mm10_sample}"  # 879352
count_lines_gzip "${unique_CAST_sample}"  # 979352
count_lines_gzip "${intersection_sample}"  # 20648

echo $(( 879352 + 20648 ))  # 900000
echo $(( 979352 + 20648 ))  # 1000000
#NOTE So, these sums make sense, but the sums for the full files do not...


#  Run full-size data ---------------------------------------------------------
samtools view -c Disteche_sample_1.dedup.mm10.corrected.bam  # 76246086
samtools view -c Disteche_sample_1.dedup.CAST.corrected.bam  # 72044686

count_lines_gzip Disteche_sample_1.dedup.mm10.corrected.mm10.AS.txt.gz  # 38123043
count_lines_gzip Disteche_sample_1.dedup.CAST.corrected.CAST.AS.txt.gz  # 36022343

echo $(( 76246086 / 38123043 ))
echo $(( 72044686 / 36022343 ))

dir_in="./data/files_bam"
dir_out="./data/files_bam"
prefix="Disteche_sample_1.dedup"
suffix="AS.txt.gz"
mm10="${dir_in}/${prefix}.mm10.corrected.mm10.${suffix}"
CAST="${dir_in}/${prefix}.CAST.corrected.CAST.${suffix}"
intersection="${dir_out}/${prefix}.intersection.${suffix}"
unique_mm10="${dir_out}/${prefix}.unique-mm10.${suffix}"
unique_CAST="${dir_out}/${prefix}.unique-CAST.${suffix}"

count_lines_gzip "${mm10}"  # 38123043
count_lines_gzip "${CAST}"  # 36022343

comm -12 <(zcat -df "${mm10}" | cut -f1) <(zcat -df "${CAST}" | cut -f1) \
| sort \
| gzip \
> "comm-12.txt.gz"

join <(zcat -df "${mm10}" | cut -f1) <(zcat -df "${CAST}" | cut -f1) \
| sort \
| gzip \
> "join.txt.gz"

mv "comm-12.txt.gz" "join.txt.gz" "./data/files_bam"
count_lines_gzip "comm-12.txt.gz"  # 31336108
count_lines_gzip "join.txt.gz"  # 31338540
count_lines_gzip "bak.Disteche_sample_1.dedup.intersection.AS.txt.gz"  # 31338540

echo $(( 31338540 - 31336108 ))  # 2432
head_10 "comm-12.txt.gz" && echo ""
head_10 "join.txt.gz"
tail_10 "comm-12.txt.gz" && echo ""
tail_10 "join.txt.gz"

find_set_complement "join.txt.gz" "comm-12.txt.gz" "what-is-in-join-not-comm.txt.gz"
find_set_complement "comm-12.txt.gz" "join.txt.gz" "what-is-in-comm-not-join.txt.gz"

comm -12 <(zcat -df "${mm10}") <(zcat -df "${CAST}") \
| sort \
| gzip \
> "comm-12.full.txt.gz"
count_lines_gzip "comm-12.full.txt.gz"  # 16250112

join <(zcat -df "${mm10}") <(zcat -df "${CAST}") \
| sort \
| gzip \
> "join-t.full.txt.gz"
count_lines_gzip "join-t.full.txt.gz"  # 31338540

#NOTE This should give a legit readout of number of intersecting elements
sort <(zcat -df "${mm10}" | cut -f1) <(zcat -df "${CAST}" | cut -f1) \
| uniq -d \
| sort \
| gzip \
> "sort.full.txt.gz"
count_lines_gzip "sort.full.txt.gz"  # 31335301
echo $(( 31335301 + 4685948 ))  # 36021249
echo $(( 31335301 + 6786560 ))  # 38121861


#  ----------------------------------------------------------------------------
zcat -d "${unique_mm10}" | head -10 && echo ""
zcat -d "${CAST}" | head -10
# AACCAATAACAAGTCAACGGCCTTCATCCGGCTCCGCTTG:1699690380
# AACCAATAACAAGTCAACGGGAATGAATAATGATTCTCGT:300995275
# AACCAATAACAAGTCAACGGGAGCTCAGCCTTATAATCAA:1165617134
# AACCAATAACAAGTCAACGGGGCAGCAGTTAAGTTCGCTG:2352957594
# AACCAATAACAAGTCAACGGGGCAGCAGTTAAGTTCGCTG:998172925
# AACCAATAACAAGTCAACGGGGCAGCAGTTATCGGCATTG:1796828347
# AACCAATAACAAGTCAACGGGGCAGCAGTTCGGCGCTCAA:1170920050
# AACCAATAACAAGTCAACGGGGCAGCAGTTGCTCCGCTTG:158709372
# AACCAATAACAAGTCAACGGGGCAGCAGTTGCTCCGCTTG:49932068
# AACCAATAACAAGTCAACGGGGCAGCAGTTTTATAATCAA:333370883
#
# AACCAATAACAAGTCAACGGAATATAGCGGCAAGTCAGCA:1417636611 -5
# AACCAATAACAAGTCAACGGAATATAGCGGTTATAATCAA:1988125074 0
# AACCAATAACAAGTCAACGGACGCTTCGTCCGAGATAAGA:2423713332 -5
# AACCAATAACAAGTCAACGGACGCTTCGTCCTCTGCGATC:1095498111 -5
# AACCAATAACAAGTCAACGGACGTACTAGGAAGTTCGCTG:353515658  0
# AACCAATAACAAGTCAACGGACGTACTAGGAAGTTCGCTG:895780689  0
# AACCAATAACAAGTCAACGGACTCGGCATACTCTGCGATC:1138973520 -5
# AACCAATAACAAGTCAACGGAGTCATGGATATCGGCATTG:707868573  0
# AACCAATAACAAGTCAACGGAGTCGCGTCGAAGTTCGCTG:994562444  0
# AACCAATAACAAGTCAACGGAGTCGCGTCGCTCTGCGATC:662408372  -5

qname="AACCAATAACAAGTCAACGGCCTTCATCCGGCTCCGCTTG:1699690380"
grep -Fxc ${qname} <(zcat "${unique_CAST}" | cut -f1)  # 1
grep -Fxc ${qname} <(zcat "${CAST}" | cut -f1)  # 1
grep -Fxc ${qname} <(zcat "${mm10}" | cut -f1)  # 0
grep -Fxc ${qname} <(zcat "${intersection}" | cut -f1)  # 0
grep -Fxc ${qname} <(zcat "${unique_mm10}" | cut -f1)  # 0

zcat -d "${CAST}" | wc -l
zcat -d "${intersection}" | wc -l
zcat -d "${unique_CAST}" | wc -l
# 36022343
# 31338540
# 4685948
echo $(( 31338540 + 4685948 ))
# 36024488  #QUESTION Why not 36022343?
echo $(( 36024488 - 36022343 ))
# 2145  #QUESTION Where do these extra values come from?

zcat -d "${mm10}" | wc -l
zcat -d "${intersection}" | wc -l
zcat -d "${unique_mm10}" | wc -l
# 38123043
# 31338540
# 6786560
echo $(( 31338540 + 6786560 ))
# 38125100  #QUESTION Why not 38123043?
echo $(( 38125100 - 38123043 ))
# 2057  #QUESTION Where do these extra values come from?

#TODO Pick up with this tomorrow: Try grep -v with column one of the "${intersection}" data

#  Restore column 2 (AS) to lists of unique elements
find_set_intersection "${unique_mm10}" "${mm10}" "${unique_mm10/.${suffix}/.tmp.${suffix}}"
find_set_intersection "${unique_CAST}" "${CAST}" "${unique_CAST/.${suffix}/.tmp.${suffix}}"
end="$(date +%s)"

calculate_run_time "${start}" "${end}" \
"Find intersections and symmetric differences between ${mm10} and ${CAST}"

#!/bin/bash

#  KA


#  Functions ------------------------------------------------------------------
displaySpinningIcon() {
    #  Display "spinning icon" while a background process runs
    spin="/|\\–/|\\-"
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s" "${2}"
        sleep .15
    done
}


head20() { head -20 "${1}"; }


samtoolsFlagstat() {
    samtools flagstat "${1}" -O tsv > "${2}" &
    displaySpinningIcon $! "Running samtools flagstat"
}


qlogin -l mfree=1G -pe serial 4
pipeline-test_env

path_base="/net/noble/vol8/kga0"
path_project="2021_kga0_4dn-mouse-cross"
path_test="results/2021-1028_Disteche_mm10_test"


#  TEST 1
#  merge_bams -----------------------------------------------------------------
bam_merge="Disteche_sample_1.q10.sorted.merged.bam"
path_merge_bams="merge_bams"

#  queryname-sorted, fixmate-styled deduplication ---------
DOI="${path_base}/${path_project}/${path_test}/${path_merge_bams}"
cd "${DOI}" ||
    {
        echo -e "Exiting: Directory ${DOI} not found."
        exit 1
    }

bam_sort_n="test.merged.sort-n.bam"
bam_fixmate="${bam_sort_n/.bam/.fixmate.bam}"
bam_sort_c="${bam_fixmate/.bam/.sort_c.bam}"
bam_markdup_r="${bam_sort_c/.bam/.markdup_r.bam}"
bam_mapq_gte_30="${bam_markdup_r/.bam/.mapq_gte_30.bam}"
bam_mapq_sort_n="${bam_mapq_gte_30/.bam/.sort-n.bam}"
bam_mapq_fixmate="${bam_mapq_sort_n/.bam/.fixmate.bam}"
bam_mapq_proper="${bam_mapq_fixmate/.bam/.proper.bam}"

echo "${bam_merge}"
echo "${bam_sort_n}"
echo "${bam_fixmate}"
echo "${bam_sort_c}"
echo "${bam_markdup_r}"
echo "${bam_mapq_gte_30}"
echo "${bam_mapq_sort_n}"
echo "${bam_mapq_fixmate}"
echo "${bam_mapq_proper}"

# -----------------
[[ -f "${bam_merge/.bam/.flagstat.tsv}" ]] ||
    {
        samtoolsFlagstat "${bam_merge}" "${bam_merge/.bam/.flagstat.tsv}"
    }

samtools sort -n -@ 4 "${bam_merge}" -o "${bam_sort_n}"
samtoolsFlagstat "${bam_sort_n}" "${bam_sort_n/.bam/.flagstat.tsv}" &
samtools view "${bam_merge}" | head -50
samtoolsFlagstat "${bam_merge}" "${bam_merge/.bam/.flagstat.tsv}"
head20 "${bam_merge/.bam/.flagstat.tsv}"

samtools fixmate -c -m -@ 4 "${bam_sort_n}" "${bam_fixmate}" &
displaySpinningIcon $! "Running samtools fixmate -n"
samtools view "${bam_fixmate}" | head -50
samtoolsFlagstat "${bam_fixmate}" "${bam_fixmate/.bam/.flagstat.tsv}"
head20 "${bam_fixmate/.bam/.flagstat.tsv}"

samtools sort -@ 4 "${bam_fixmate}" -o "${bam_sort_c}" &
displaySpinningIcon $! "Running samtools sort (coordinate sorting)"
samtools view "${bam_sort_c}" | head -50  # Reads are no longer paired
samtoolsFlagstat "${bam_sort_c}" "${bam_sort_c/.bam/.flagstat.tsv}"
head20 "${bam_sort_c/.bam/.flagstat.tsv}"

samtools markdup -r -s -@ 4 "${bam_sort_c}" -f "${bam_markdup_r/.bam/.txt}" "${bam_markdup_r}" &
displaySpinningIcon $! "Running markdup -r -s"
samtools view "${bam_markdup_r}" | head -50  # Reads are no longer paired
samtoolsFlagstat "${bam_markdup_r}" "${bam_markdup_r/.bam/.flagstat.tsv}"
head20 "${bam_markdup_r/.bam/.flagstat.tsv}"

samtools view -h -b -q 30 -@ 4 "${bam_markdup_r}" -o "${bam_mapq_gte_30}" &
displaySpinningIcon $! "samtools view -h -b -q 30 (excluding reads with MAPQ < 30)"
samtools view "${bam_mapq_gte_30}" | head -50
samtoolsFlagstat "${bam_mapq_gte_30}" "${bam_mapq_gte_30/.bam/.flagstat.tsv}"
head20 "${bam_mapq_gte_30/.bam/.flagstat.tsv}"

samtools sort -n -@ 4 "${bam_mapq_gte_30}" -o "${bam_mapq_sort_n}" &
displaySpinningIcon $! "Running samtools fixmate"
samtools view "${bam_mapq_sort_n}" | head -50
samtoolsFlagstat "${bam_mapq_sort_n}" "${bam_mapq_sort_n/.bam/.flagstat.tsv}" &
head20 "${bam_mapq_sort_n/.bam/.flagstat.tsv}"

samtools fixmate -c -m -@ 4 "${bam_mapq_sort_n}" "${bam_mapq_fixmate}" &
displaySpinningIcon $! "Running samtools fixmate"
samtools view "${bam_mapq_fixmate}" | head -50
samtoolsFlagstat "${bam_mapq_fixmate}" "${bam_mapq_fixmate/.bam/.flagstat.tsv}" &
head20 "${bam_mapq_fixmate/.bam/.flagstat.tsv}"

samtools view -h -b -f 0x2 -@ 4 "${bam_mapq_fixmate}" > "${bam_mapq_proper}" &
displaySpinningIcon $! "Running view -h -b -f 0x2 (extracting only properly paired reads)"
samtools view "${bam_mapq_proper}" | head -50
samtoolsFlagstat "${bam_mapq_proper}" "${bam_mapq_proper/.bam/.flagstat.tsv}"
head20 "${bam_mapq_proper/.bam/.flagstat.tsv}"

mkdir -p test_1_markdup-r
mv test.merged.sort-n.* test_1_markdup

# cd test_1_markdup-r ||
#     {
#         echo -e "Exiting: Directory ${DOI} not found."
#         exit 1
#     }


#  TEST 2
#  merge_bams -----------------------------------------------------------------
bam_merge="Disteche_sample_1.q10.sorted.merged.bam"
path_merge_bams="merge_bams"

#  queryname-sorted, fixmate-styled deduplication ---------
DOI="${path_base}/${path_project}/${path_test}/${path_merge_bams}"
cd "${DOI}" ||
    {
        echo -e "Exiting: Directory ${DOI} not found."
        exit 1
    }

bam_mapq_gte_30="test.merged.mapq_gte_30.bam"
bam_mapq_sort_n="${bam_mapq_gte_30/.bam/.sort-n.bam}"
bam_mapq_fixmate="${bam_mapq_sort_n/.bam/.fixmate.bam}"
bam_mapq_proper="${bam_mapq_fixmate/.bam/.proper.bam}"
bam_mapq_sort_c="${bam_mapq_proper/.bam/.sort-c.bam}"
bam_mapq_markdup_r="${bam_mapq_sort_c/.bam/.markdup_r.bam}"
bam_mapq_repair="${bam_mapq_markdup_r/.bam/.repair.bam}"
bam_mapq_repair_no_dummy="${bam_mapq_markdup_r/.bam/.repair_no-dummy.bam}"

echo "${bam_mapq_gte_30}"
echo "${bam_mapq_sort_n}"
echo "${bam_mapq_fixmate}"
echo "${bam_mapq_proper}"
echo "${bam_mapq_sort_c}"
echo "${bam_mapq_markdup_r}"
echo "${bam_mapq_repair}"
echo "${bam_mapq_repair_no_dummy}"

# -----------------
[[ -f "${bam_merge/.bam/.flagstat.tsv}" ]] ||
    {
        samtoolsFlagstat "${bam_merge}" "${bam_merge/.bam/.flagstat.tsv}"
    }
head20 "${bam_merge/.bam/.flagstat.tsv}"

samtools view -h -b -q 30 -@ 4 "${bam_merge}" -o "${bam_mapq_gte_30}" &
displaySpinningIcon $! "samtools view -h -b -q 30 (excluding reads with MAPQ < 30)"
samtools view "${bam_mapq_gte_30}" | head -50
samtoolsFlagstat "${bam_mapq_gte_30}" "${bam_mapq_gte_30/.bam/.flagstat.tsv}"
head20 "${bam_mapq_gte_30/.bam/.flagstat.tsv}"

samtools sort -n -@ 4 "${bam_mapq_gte_30}" -o "${bam_mapq_sort_n}" &
displaySpinningIcon $! "Running samtools sort -n (queryname sorting)"
samtools view "${bam_mapq_sort_n}" | head -50
samtoolsFlagstat "${bam_mapq_sort_n}" "${bam_mapq_sort_n/.bam/.flagstat.tsv}"
head20 "${bam_mapq_sort_n/.bam/.flagstat.tsv}"

samtools fixmate -c -m -@ 4 "${bam_mapq_sort_n}" "${bam_mapq_fixmate}" &
displaySpinningIcon $! "Running samtools fixmate -n"
samtools view "${bam_mapq_fixmate}" | head -50
samtoolsFlagstat "${bam_mapq_fixmate}" "${bam_mapq_fixmate/.bam/.flagstat.tsv}"
head20 "${bam_mapq_fixmate/.bam/.flagstat.tsv}"

samtools view -h -b -f 0x2 -@ 4 "${bam_mapq_fixmate}" > "${bam_mapq_proper}" &
displaySpinningIcon $! "Running view -h -b -f 0x2 (extracting only properly paired reads)"
samtools view "${bam_mapq_proper}" | head -50
samtoolsFlagstat "${bam_mapq_proper}" "${bam_mapq_proper/.bam/.flagstat.tsv}"
head20 "${bam_mapq_proper/.bam/.flagstat.tsv}"

samtools sort -@ 4 "${bam_mapq_proper}" -o "${bam_mapq_sort_c}" &
displaySpinningIcon $! "Running samtools sort (coordinate sorting)"
samtools view "${bam_mapq_sort_c}" | head -50
samtoolsFlagstat "${bam_mapq_sort_c}" "${bam_mapq_sort_c/.bam/.flagstat.tsv}"
head20 "${bam_mapq_sort_c/.bam/.flagstat.tsv}"

samtools markdup -r -s -@ 4 "${bam_mapq_sort_c}" -f "${bam_mapq_markdup_r/.bam/.txt}" "${bam_mapq_markdup_r}" &
displaySpinningIcon $! "Running markdup -r -s"
samtools view "${bam_mapq_markdup_r}" | head -50  # Reads are no longer paired
samtoolsFlagstat "${bam_mapq_markdup_r}" "${bam_mapq_markdup_r/.bam/.flagstat.tsv}"
head20 "${bam_mapq_markdup_r/.bam/.flagstat.tsv}"

repair -c -i "${bam_mapq_markdup_r}" -o "${bam_mapq_repair}" &
displaySpinningIcon $! "Running Subread repair (getting mates together)"
samtools view "${bam_mapq_repair}" | head -50  # Reads are...
samtoolsFlagstat "${bam_mapq_repair}" "${bam_mapq_repair/.bam/.flagstat.tsv}"
head20 "${bam_mapq_repair/.bam/.flagstat.tsv}"

cat << EOF > "${bam_mapq_repair/.bam/.out.txt}"
Finished scanning the input file. Processing unpaired reads.

All finished in 2.06 minutes
Total input reads: 75324062 ; Unpaired reads: 150
EOF

repair -d -c -i "${bam_mapq_markdup_r}" -o "${bam_mapq_repair_no_dummy}" &
displaySpinningIcon $! "Running Subread repair (getting mates together; not adding dummy reads to unpaired reads)"
samtools view "${bam_mapq_repair_no_dummy}" | head -50  # Reads are...
samtoolsFlagstat "${bam_mapq_repair_no_dummy}" "${bam_mapq_repair_no_dummy/.bam/.flagstat.tsv}"
head20 "${bam_mapq_repair_no_dummy/.bam/.flagstat.tsv}"

cat << EOF > "${bam_mapq_repair_no_dummy/.bam/.out.txt}"
Finished scanning the input file. Processing unpaired reads.
All finished in 2.11 minutes
Total input reads: 75324062 ; Unpaired reads: 150
EOF

# #QUESTION How about this? #ANSWER No...
# samtools view -h -b -f 0x1 -@ 4 "${bam_mapq_markdup_r}" > "does.this.work.bam" &
# displaySpinningIcon $! "Running view -h -b -f 0x2 (extracting only properly paired reads)"
# samtools view "does.this.work.bam" | head -50
# samtoolsFlagstat "does.this.work.bam" "does.this.work.bam.flagstat.tsv"
# head20 "does.this.work.bam.flagstat.tsv"

mkdir -p test_2_markdup-r
mv test.merged.mapq_gte_30.* test_2_markdup-r

# cd test_2_markdup-r ||
#     {
#         echo -e "Exiting: Directory ${DOI} not found."
#         exit 1
#     }

path_unique="${1:-"/net/noble/vol8/kga0/2021_kga0_4dn-mouse-cross/results/2021-1028_Disteche_mm10_test/get_unique_fragments"}"
cd "${path_unique}" ||
    {
        echo -e "Exiting: Directory ${path_unique} not found."
        exit 1
    }

bam_unique="Disteche_sample_1.dedup.bam"
bam_unique_sort_n="${bam_unique/.bam/.sort-n.bam}"
bam_unique_fixmate="${bam_unique_sort_n/.bam/.fixmate.bam}"
bam_unique_markdup_r="${bam_unique_fixmate/.bam/.markdup-r.bam}"
bam_unique_sort_c="${bam_unique_markdup_r/.bam/.sort-c.bam}"
bam_unique_repair="${bam_unique_sort_c/.bam/.repair.bam}"

echo -e "${bam_unique}"
echo -e "${bam_unique_sort_n}"
echo -e "${bam_unique_markdup_r}"

samtools sort -n -@ 4 "${bam_unique}" -o "${bam_unique_sort_n}" &
displaySpinningIcon $! "Running samtools sort (coordinate sorting)"
samtools view "${bam_unique_sort_n}" | head -50
samtoolsFlagstat "${bam_unique_sort_n}" "${bam_unique_sort_n/.bam/.flagstat.tsv}"
head20 "${bam_unique_sort_n/.bam/.flagstat.tsv}"

samtools fixmate -c -m -@ 4 "${bam_unique_sort_n}" "${bam_unique_fixmate}" &
displaySpinningIcon $! "Running samtools fixmate -n"
samtools view "${bam_unique_fixmate}" | head -50
samtoolsFlagstat "${bam_unique_fixmate}" "${bam_unique_fixmate/.bam/.flagstat.tsv}"
head20 "${bam_unique_fixmate/.bam/.flagstat.tsv}"

samtools sort -@ 4 "${bam_unique_fixmate}" -o "${bam_unique_sort_c}" &
displaySpinningIcon $! "Running samtools sort (coordinate sorting)"
samtools view "${bam_unique_sort_c}" | head -50
samtoolsFlagstat "${bam_unique_sort_c}" "${bam_unique_sort_c/.bam/.flagstat.tsv}"
head20 "${bam_unique_sort_c/.bam/.flagstat.tsv}"

samtools markdup -r -s -@ 4 "${bam_unique_sort_c}" -f "${bam_unique_markdup_r/.bam/.txt}" "${bam_unique_markdup_r}" &
displaySpinningIcon $! "Running markdup -r -s"
samtools view "${bam_unique_markdup_r}" | head -50  # Reads are no longer paired
samtoolsFlagstat "${bam_unique_markdup_r}" "${bam_unique_markdup_r/.bam/.flagstat.tsv}"
head20 "${bam_unique_markdup_r/.bam/.flagstat.tsv}"

repair -c -i "${bam_unique_sort_c}" -o "${bam_unique_repair}" &
displaySpinningIcon $! "Running Subread repair (getting mates together)"
samtools view "${bam_unique_repair}" | head -50  # Reads are...
samtoolsFlagstat "${bam_unique_repair}" "${bam_unique_repair/.bam/.flagstat.tsv}"
head20 "${bam_unique_repair/.bam/.flagstat.tsv}"

mkdir -p test_4_markdup-r/test_4_markdup-r_tsv-txt
mv -- *sort-n* test_4_markdup-r
mv -- *.{tsv,txt} test_4_markdup-r/test_4_markdup-r_tsv-txt

#!/bin/bash

#  test-scraps.sh
#  KA

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


cd "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross" ||
    {
        echo -e "Exiting: Directory not found."
        # exit 1
    }

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_chr19/" \
-c "chr19" \
-r "TRUE" \
-b "TRUE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_chr18/" \
-c "chr18" \
-r "TRUE" \
-b "FALSE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_chr17/" \
-c "chr17" \
-r "FALSE" \
-b "TRUE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/kga0/Disteche_sample_1.dedup.bam" \
-o "./data/kga0/split-index-repair_Disteche_sample_1.dedup_all/" \
-c "all" \
-r "TRUE" \
-b "TRUE" \
-p 6

# -h <print this help message and exit>
# -u <use safe mode: "TRUE" or "FALSE" (logical)>
# -i <bam infile, including path (chr)>
# -o <path for split bam file(s) (chr); path will be made if it does
#     not exist>
# -c <chromosome(s) to split out (chr); for example, "chr1" for
#     chromosome 1, "chrX" for chromosome X, "all" for all
#     chromosomes>
# -r <use Subread repair on split bam files: "TRUE" or "FALSE"
#     (logical)>
# -b <if "-r TRUE", create bed files from split bam files: "TRUE"
#     or "FALSE" (logical); argument "-b" only needed when "-r
#     TRUE">
# -p <number of cores for parallelization (int >= 1)>



#ANSWER Does this work?


reformat \
in="${infile/.bam/.only-reverse.bam}" \
out="${infile/.bam/.only-reverse.sample-50.bam}" \
sampleseed=24 \
samplereadstarget=50 &
displaySpinningIcon $! "Randomly sampling reverse-strand bam... "

samtools view "${infile/.bam/.only-forward.sample-50.bam}" | head -50 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-50.bam}" | head -50 && echo ""

samtools view "${infile/.bam/.only-forward.sample-50.bam}" | tail -6 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-50.bam}" | tail -6 && echo ""
#ANSWER It works!

reformat \
in="${infile/.bam/.only-forward.bam}" \
out="${infile/.bam/.only-forward.sample-50000.bam}" \
sampleseed=24 \
samplereadstarget=50000 &
displaySpinningIcon $! "Randomly sampling forward-strand bam... "

reformat \
in="${infile/.bam/.only-reverse.bam}" \
out="${infile/.bam/.only-reverse.sample-50000.bam}" \
sampleseed=24 \
samplereadstarget=50000 &
displaySpinningIcon $! "Randomly sampling reverse-strand bam... "

samtools view "${infile/.bam/.only-forward.sample-50000.bam}" | head -50 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-50000.bam}" | head -50 && echo ""

samtools view "${infile/.bam/.only-forward.sample-50000.bam}" | tail -50 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-50000.bam}" | tail -50 && echo ""

samtools sort -@ 4 "${infile/.bam/.only-forward.sample-50.bam}" -o "${infile/.bam/.only-forward.sample-50.sorted.bam}" & displaySpinningIcon $! "Sorting... "
samtools sort -@ 4 "${infile/.bam/.only-reverse.sample-50.bam}" -o "${infile/.bam/.only-reverse.sample-50.sorted.bam}" & displaySpinningIcon $! "Sorting... "
samtools sort -@ 4 "${infile/.bam/.only-forward.sample-50000.bam}" -o "${infile/.bam/.only-forward.sample-50000.sorted.bam}" & displaySpinningIcon $! "Sorting... "
samtools sort -@ 4 "${infile/.bam/.only-reverse.sample-50000.bam}" -o "${infile/.bam/.only-reverse.sample-50000.sorted.bam}" & displaySpinningIcon $! "Sorting... "

samtools view "${infile/.bam/.only-forward.sample-50.sorted.bam}" | head -6 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-50.sorted.bam}" | head -6 && echo ""

mkdir -p files_bam_test

samtools merge -@ 4 \
"${infile/.bam/.only-forward.sample-50.sorted.bam}" \
"${infile/.bam/.only-reverse.sample-50.sorted.bam}" \
-o "${infile/.bam/.sample-100.bam}" &
displaySpinningIcon $! "Merging bam files... "

samtools merge -@ 4 \
"${infile/.bam/.only-forward.sample-50000.sorted.bam}" \
"${infile/.bam/.only-reverse.sample-50000.sorted.bam}" \
-o "${infile/.bam/.sample-100000.bam}" &
displaySpinningIcon $! "Merging bam files... "

repair -d -c -T 6 \
-i "${infile/.bam/.sample-100.bam}" \
-o "${infile/.bam/.sample-100.repaired.bam}" &
displaySpinningIcon $! "Repairing sampled bam files... "

repair -d -c -T 6 \
-i "${infile/.bam/.sample-100000.bam}" \
-o "${infile/.bam/.sample-100000.repaired.bam}" &
displaySpinningIcon $! "Repairing sampled bam files... "

samtools view "${infile/.bam/.sample-100.repaired.bam}" | head -100

samtools view "${infile/.bam/.sample-100000.bam}" | head -100
samtools view "${infile/.bam/.sample-100000.repaired.bam}" | head -100

samtools view "${infile/.bam/.sample-100000.bam}" | tail -100
samtools view "${infile/.bam/.sample-100000.repaired.bam}" | tail -100


reformat \
in="${infile/.bam/.only-forward.bam}" \
out="${infile/.bam/.only-forward.sample-150000.bam}" \
sampleseed=24 \
samplereadstarget=150000 &
displaySpinningIcon $! "Randomly sampling forward-strand bam... "

reformat \
in="${infile/.bam/.only-reverse.bam}" \
out="${infile/.bam/.only-reverse.sample-150000.bam}" \
sampleseed=24 \
samplereadstarget=150000 &
displaySpinningIcon $! "Randomly sampling reverse-strand bam... "

samtools view "${infile/.bam/.only-forward.sample-150000.bam}" | head -50 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-150000.bam}" | head -50 && echo ""

samtools view "${infile/.bam/.only-forward.sample-150000.bam}" | tail -50 && echo ""
samtools view "${infile/.bam/.only-reverse.sample-150000.bam}" | tail -50 && echo ""

samtools sort -@ 4 "${infile/.bam/.only-forward.sample-150000.bam}" -o "${infile/.bam/.only-forward.sample-150000.sorted.bam}" & displaySpinningIcon $! "Sorting... "
samtools sort -@ 4 "${infile/.bam/.only-reverse.sample-150000.bam}" -o "${infile/.bam/.only-reverse.sample-150000.sorted.bam}" & displaySpinningIcon $! "Sorting... "

samtools merge -@ 4 \
"${infile/.bam/.only-forward.sample-150000.sorted.bam}" \
"${infile/.bam/.only-reverse.sample-150000.sorted.bam}" \
-o "${infile/.bam/.sample-300000.bam}" &
displaySpinningIcon $! "Merging bam files... "

repair -d -c -T 6 \
-i "${infile/.bam/.sample-300000.bam}" \
-o "${infile/.bam/.sample-300000.repaired.bam}" &
displaySpinningIcon $! "Repairing sampled bam files... "


samtools view "${infile/.bam/.sample-300000.bam}" | head -100
samtools view "${infile/.bam/.sample-300000.repaired.bam}" | head -100

samtools view "${infile/.bam/.sample-300000.bam}" | tail -100
samtools view "${infile/.bam/.sample-300000.repaired.bam}" | tail -100

cp \
"${infile/.bam/.sample-100.bam}" \
"${infile/.bam/.sample-100000.bam}" \
"${infile/.bam/.sample-300000.bam}" \
./files_bam_test

cd "./files_bam_test" ||
    {
        echo "Exiting: Directory not found."
        # exit 1
    }

mv "${infile/.bam/.sample-100.bam}" "test.100.bam"
mv "${infile/.bam/.sample-100000.bam}" "test.100000.bam"
mv "${infile/.bam/.sample-300000.bam}" "test.300000.bam"


bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/files_bam_test/test.300000.bam" \
-o "./data/split-index-repair_test.300000_all" \
-c "all" \
-r "TRUE" \
-b "TRUE" \
-p 4

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/files_bam_test/test.300000.bam" \
-o "./data/split-index-repair_test.300000_chr1" \
-c "chr1" \
-r "TRUE" \
-b "TRUE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/files_bam_test/test.100.bam" \
-o "./data/split-index-repair_test.100_all" \
-c "all" \
-r "TRUE" \
-b "TRUE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/files_bam_test/test.100000.bam" \
-o "./data/split-index-repair_test.100000_all" \
-c "all" \
-r "TRUE" \
-b "TRUE" \
-p 6

bash ./bin/split-index-repair_bam.sh \
-u "FALSE" \
-i "./data/files_bam_test/test.100000.bam" \
-o "./data/split-index-repair_test.100000_chr8" \
-c "chr8" \
-r "TRUE" \
-b "TRUE" \
-p 6

# cd 


bash ./bin/lift_strain-to-mm10.sh \
-u "FALSE" \
-i "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/Disteche_sample_1.dedup.chr19.pos.bed" \
-o "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0" \
-s 1 \
-c "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/liftOver/GCA_001624445.1ToMm10.over.chain.gz"

bash ./bin/lift_strain-to-mm10.sh \
-u "TRUE" \
-i "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0/Disteche_sample_1.dedup.chr17.pos.bed" \
-o "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0" \
-s 1 \
-c "/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/liftOver/GCA_001624445.1ToMm10.over.chain.gz"

# -h <print this help message and exit>
# -u <use safe mode: TRUE or FALSE (logical)>
# -i <bed infile, including path (chr)>
# -o <"lifted" bed outfile, including path (chr)>
# -s <strain for performing liftOver of bed files (int = 1 | int = 2);
#     options: "1" for "CAST-EiJ", "2" for "129S1-SvImJ>"
# -c <gzipped liftOver chain file for strain, including path (chr);
#     note: for liftOver to work, the liftOver strain chain should
#     match the strain set in argument "-s">

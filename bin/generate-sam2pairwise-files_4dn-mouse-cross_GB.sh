#!/bin/bash

# generate-sam2pairwise-files_4dn-mouse-cross-tests.sh
# KA

script="generate-sam2pairwise-files_4dn-mouse-cross-tests.sh"

command -v parallel &> /dev/null ||
{
    echo "Exiting: parallel not found. Install GNU Parallel."
    return 1 2> /dev/null
    exit 1
}

reformat="/Users/kalavattam/bbmap/reformat.sh"
command -v "${reformat}" &> /dev/null ||
{
    echo "Exiting: reformat.sh not found. Install BBMap."
    return 1 2> /dev/null
    exit 1
}

command -v sam2pairwise &> /dev/null ||
{
    echo "Exiting: sam2pairwise not found. Install sam2pairwise."
    return 1 2> /dev/null
    exit 1
}

directory="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/results/kga0/2022-0207_segregated_reads.thresh_SNP_1.thresh_Q_30"
cd "${directory}" || 
{
    echo "Directory not found; look into this."
}

prefix="2022-0207_segregated_reads.thresh_SNP_1.thresh_Q_30"

#  Convert to extended CIGAR format
bamAlt="${prefix}.alt.bam"
bamAmbig="${prefix}.ambig.bam"
bamContra="${prefix}.contra.bam"
bamRef="${prefix}.ref.bam"

cigarAlt="${bamAlt/.bam/.sam2pairwise.txt}"
cigarAmbig="${bamAmbig/.bam/.sam2pairwise.txt}"
cigarContra="${bamContra/.bam/.sam2pairwise.txt}"
cigarRef="${bamRef/.bam/.sam2pairwise.txt}"

errAlt="${cigarAlt/.txt/.stderr.txt}"
errAmbig="${cigarAmbig/.txt/.stderr.txt}"
errContra="${cigarContra/.txt/.stderr.txt}"
errRef="${cigarRef/.txt/.stderr.txt}"

#  Generate CIGAR constructions with sam2pairwise
samtools view "${bamAlt}" | sam2pairwise > "${cigarAlt}" 2> "${errAlt}"
samtools view "${bamAmbig}" | sam2pairwise > "${cigarAmbig}" 2> "${errAmbig}"
samtools view "${bamContra}" | sam2pairwise > "${cigarContra}" 2> "${errContra}"
samtools view "${bamRef}" | sam2pairwise > "${cigarRef}" 2> "${errRef}"


#  Replace mismatch-spaces with mismatch-asterisks ----------------------------
replaceSpacesAsterisks() {
    awk '{gsub(" ", "*", $0); print}' "${1}" \
    > "${1/.txt/.tmp.txt}" && mv "${1/.txt/.tmp.txt}" "${1}"
}
replaceSpacesAsterisks "${cigarAlt}"
replaceSpacesAsterisks "${cigarAmbig}"
replaceSpacesAsterisks "${cigarContra}"
replaceSpacesAsterisks "${cigarRef}"


#  Munge the output of sam2pairwise: Alt --------------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigarAlt}" | awk 'NR%4==1' > "${cigarAlt/.txt/-1.txt}"
cat "${cigarAlt}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigarAlt/.txt/-2.txt}"
cat "${cigarAlt}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigarAlt/.txt/-3.txt}"
cat "${cigarAlt}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigarAlt/.txt/-4.txt}"

paste "${cigarAlt/.txt/-1.txt}" "${cigarAlt/.txt/-2.txt}" "${cigarAlt/.txt/-3.txt}" "${cigarAlt/.txt/-4.txt}" > "${cigarAlt/.txt/.munged.txt}"
rm "${cigarAlt/.txt/-1.txt}" "${cigarAlt/.txt/-2.txt}" "${cigarAlt/.txt/-3.txt}" "${cigarAlt/.txt/-4.txt}"
cat "${cigarAlt/.txt/.munged.txt}" | wc -l  # 32045
samtools view "${bamAlt}" | wc -l  # 32526

cat "${errAlt}" | awk 'NR%4==1' > "${errAlt/.txt/-1.txt}"
cat "${errAlt}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${errAlt/.txt/-2.txt}"

paste "${errAlt/.txt/-1.txt}" "${errAlt/.txt/-2.txt}" "${errAlt/.txt/-2.txt}" "${errAlt/.txt/-2.txt}" > "${errAlt/.txt/.munged.txt}"
rm "${errAlt/.txt/-1.txt}" "${errAlt/.txt/-2.txt}"
cat "${errAlt/.txt/.munged.txt}" | wc -l  # 481
echo $(( 32045 + 481 ))  # 32526

cat "${cigarAlt/.txt/.munged.txt}" "${errAlt/.txt/.munged.txt}" > "${cigarAlt/.txt/.munged.all.txt}"
cat "${cigarAlt/.txt/.munged.all.txt}" | wc -l  # 32526


#  Munge the output of sam2pairwise: Ambig ------------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigarAmbig}" | awk 'NR%4==1' > "${cigarAmbig/.txt/-1.txt}"
cat "${cigarAmbig}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigarAmbig/.txt/-2.txt}"
cat "${cigarAmbig}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigarAmbig/.txt/-3.txt}"
cat "${cigarAmbig}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigarAmbig/.txt/-4.txt}"

paste "${cigarAmbig/.txt/-1.txt}" "${cigarAmbig/.txt/-2.txt}" "${cigarAmbig/.txt/-3.txt}" "${cigarAmbig/.txt/-4.txt}" > "${cigarAmbig/.txt/.munged.txt}"
rm "${cigarAmbig/.txt/-1.txt}" "${cigarAmbig/.txt/-2.txt}" "${cigarAmbig/.txt/-3.txt}" "${cigarAmbig/.txt/-4.txt}"
cat "${cigarAmbig/.txt/.munged.txt}" | wc -l  # 170838
samtools view "${bamAmbig}" | wc -l  # 171886

cat "${errAmbig}" | awk 'NR%4==1' > "${errAmbig/.txt/-1.txt}"
cat "${errAmbig}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${errAmbig/.txt/-2.txt}"

paste "${errAmbig/.txt/-1.txt}" "${errAmbig/.txt/-2.txt}" "${errAmbig/.txt/-2.txt}" "${errAmbig/.txt/-2.txt}" > "${errAmbig/.txt/.munged.txt}"
rm "${errAmbig/.txt/-1.txt}" "${errAmbig/.txt/-2.txt}"
cat "${errAmbig/.txt/.munged.txt}" | wc -l  # 1048
echo $(( 170838 + 1048 ))  # 171886

cat "${cigarAmbig/.txt/.munged.txt}" "${errAmbig/.txt/.munged.txt}" > "${cigarAmbig/.txt/.munged.all.txt}"
cat "${cigarAmbig/.txt/.munged.all.txt}" | wc -l  # 171886


#  Munge the output of sam2pairwise: Contra -----------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigarContra}" | awk 'NR%4==1' > "${cigarContra/.txt/-1.txt}"
cat "${cigarContra}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigarContra/.txt/-2.txt}"
cat "${cigarContra}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigarContra/.txt/-3.txt}"
cat "${cigarContra}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigarContra/.txt/-4.txt}"

paste "${cigarContra/.txt/-1.txt}" "${cigarContra/.txt/-2.txt}" "${cigarContra/.txt/-3.txt}" "${cigarContra/.txt/-4.txt}" > "${cigarContra/.txt/.munged.txt}"
rm "${cigarContra/.txt/-1.txt}" "${cigarContra/.txt/-2.txt}" "${cigarContra/.txt/-3.txt}" "${cigarContra/.txt/-4.txt}"
cat "${cigarContra/.txt/.munged.txt}" | wc -l  # 27
samtools view "${bamContra}" | wc -l  # 28

cat "${errContra}" | awk 'NR%4==1' > "${errContra/.txt/-1.txt}"
cat "${errContra}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${errContra/.txt/-2.txt}"

paste "${errContra/.txt/-1.txt}" "${errContra/.txt/-2.txt}" "${errContra/.txt/-2.txt}" "${errContra/.txt/-2.txt}" > "${errContra/.txt/.munged.txt}"
rm "${errContra/.txt/-1.txt}" "${errContra/.txt/-2.txt}"
cat "${errContra/.txt/.munged.txt}" | wc -l  # 481
echo $(( 27 + 1 ))  # 28

cat "${cigarContra/.txt/.munged.txt}" "${errContra/.txt/.munged.txt}" > "${cigarContra/.txt/.munged.all.txt}"
cat "${cigarContra/.txt/.munged.all.txt}" | wc -l  # 28


#  Munge the output of sam2pairwise: Ref --------------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigarRef}" | awk 'NR%4==1' > "${cigarRef/.txt/-1.txt}"
cat "${cigarRef}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigarRef/.txt/-2.txt}"
cat "${cigarRef}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigarRef/.txt/-3.txt}"
cat "${cigarRef}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigarRef/.txt/-4.txt}"

paste "${cigarRef/.txt/-1.txt}" "${cigarRef/.txt/-2.txt}" "${cigarRef/.txt/-3.txt}" "${cigarRef/.txt/-4.txt}" > "${cigarRef/.txt/.munged.txt}"
rm "${cigarRef/.txt/-1.txt}" "${cigarRef/.txt/-2.txt}" "${cigarRef/.txt/-3.txt}" "${cigarRef/.txt/-4.txt}"
cat "${cigarRef/.txt/.munged.txt}" | wc -l  # 29811
samtools view "${bamRef}" | wc -l  # 29896

cat "${errRef}" | awk 'NR%4==1' > "${errRef/.txt/-1.txt}"
cat "${errRef}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${errRef/.txt/-2.txt}"

paste "${errRef/.txt/-1.txt}" "${errRef/.txt/-2.txt}" "${errRef/.txt/-2.txt}" "${errRef/.txt/-2.txt}" > "${errRef/.txt/.munged.txt}"
rm "${errRef/.txt/-1.txt}" "${errRef/.txt/-2.txt}"
cat "${errRef/.txt/.munged.txt}" | wc -l  # 85
echo $(( 29811 + 85 ))  # 29896

cat "${cigarRef/.txt/.munged.txt}" "${errRef/.txt/.munged.txt}" > "${cigarRef/.txt/.munged.all.txt}"
cat "${cigarRef/.txt/.munged.all.txt}" | wc -l  # 29896


#  Add headers to the pertinent files -----------------------------------------
awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigarAlt/.txt/.munged.all.txt}" \
> /tmp/out && mv /tmp/out "${cigarAlt/.txt/.munged.all.txt}"

awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigarAmbig/.txt/.munged.all.txt}" \
> /tmp/out && mv /tmp/out "${cigarAmbig/.txt/.munged.all.txt}"

awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigarContra/.txt/.munged.all.txt}" \
> /tmp/out && mv /tmp/out "${cigarContra/.txt/.munged.all.txt}"

awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigarRef/.txt/.munged.all.txt}" \
> /tmp/out && mv /tmp/out "${cigarRef/.txt/.munged.all.txt}"

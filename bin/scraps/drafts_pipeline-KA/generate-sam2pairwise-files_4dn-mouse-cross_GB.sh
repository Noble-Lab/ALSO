#!/bin/bash

#  generate-sam2pairwise-files_4dn-mouse-cross-tests.sh
#  KA

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

samfixcigar="/Users/kalavattam/jvarkit/dist/samfixcigar.jar"
[[ -f "${samfixcigar}" ]] ||
{
    echo "Exiting: samfixcigar not found. Install samfixcigar:"
    echo "http://lindenb.github.io/jvarkit/SamFixCigar.html"
    return 1 2> /dev/null
    exit 1
}

command -v sam2pairwise &> /dev/null ||
{
    echo "Exiting: sam2pairwise not found. Install sam2pairwise."
    return 1 2> /dev/null
    exit 1
}

directory_base="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/results/kga0"
directory_data="$(date '+%Y-%m%d')_segregated_reads.thresh_SNP_1.thresh_Q_30"
cd "${directory_base}/${directory_data}" || 
{
    echo "Directory not found; look into this."
}

prefix="$(date '+%Y-%m%d')_segregated_reads.thresh_SNP_1.thresh_Q_30"
chromosome="chrX"


#  Convert to extended CIGAR format -------------------------------------------
# reformatCIGARstring() {
#     java -jar "${samfixcigar}" --samoutputformat BAM --reference "${1}" "${2}"
# }

initAlt="${prefix}.${chromosome}.alt.bam"
initAmbig="${prefix}.${chromosome}.ambig.bam"
initContra="${prefix}.${chromosome}.contra.bam"
initRef="${prefix}.${chromosome}.ref.bam"

bamAlt="${initAlt/.bam/.extendedCIGAR.bam}"
bamAmbig="${initAmbig/.bam/.extendedCIGAR.bam}"
bamContra="${initContra/.bam/.extendedCIGAR.bam}"
bamRef="${initRef/.bam/.extendedCIGAR.bam}"

# mm10_fa="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/genomes/mm10.fa"
# reformatCIGARstring "${mm10_fa}" "${initAlt}"
# reformatCIGARstring "${mm10_fa}" "${initAmbig}"
# reformatCIGARstring "${mm10_fa}" "${initContra}"
# reformatCIGARstring "${mm10_fa}" "${initRef}"

reformatCIGARstring() {
    "${reformat}" in="${1}" out="${2}" sam=1.4 2> "${2/.bam/.err}"
}

reformatCIGARstring "${initAlt}" "${bamAlt}"
reformatCIGARstring "${initAmbig}" "${bamAmbig}"
reformatCIGARstring "${initContra}" "${bamContra}"
reformatCIGARstring "${initRef}" "${bamRef}"
#TODO Sometimes the above throws dydld errors, but sometimes not: why?


#  Generate CIGAR constructions with sam2pairwise -----------------------------
generateCIGARconstructions() {
    samtools view "${1}" | sam2pairwise > "${1/.bam/.sam2pairwise.txt}"
}

generateCIGARconstructions "${bamAlt}"
generateCIGARconstructions "${bamAmbig}"
generateCIGARconstructions "${bamContra}"
generateCIGARconstructions "${bamRef}"

cigarAlt="${bamAlt/.bam/.sam2pairwise.txt}"
cigarAmbig="${bamAmbig/.bam/.sam2pairwise.txt}"
cigarContra="${bamContra/.bam/.sam2pairwise.txt}"
cigarRef="${bamRef/.bam/.sam2pairwise.txt}"


#  Replace mismatch-spaces with mismatch-asterisks ----------------------------
replaceSpacesAsterisks() {
    awk '{gsub(" ", "*", $0); print}' "${1}" > \
    "${1/.txt/.tmp.txt}" && mv "${1/.txt/.tmp.txt}" "${1}"
}
replaceSpacesAsterisks "${cigarAlt}"
replaceSpacesAsterisks "${cigarAmbig}"
replaceSpacesAsterisks "${cigarContra}"
replaceSpacesAsterisks "${cigarRef}"


#  Munge the output of sam2pairwise -------------------------------------------
mungePairwise() {
    #  Munge the output of sam2pairwise

    #  Want every second, third, and fourth row to be additional tabbed entries
    #+ on the first row; thus, split each line to separate temporary files
    cat "${1}" | awk 'NR%4==1' > "${1/.txt/-1.txt}"
    cat "${1}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${1/.txt/-2.txt}"
    cat "${1}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${1/.txt/-3.txt}"
    cat "${1}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${1/.txt/-4.txt}"

    #  Column-bind the separate files
    paste "${1/.txt/-1.txt}" \
    "${1/.txt/-2.txt}" \
    "${1/.txt/-3.txt}" \
    "${1/.txt/-4.txt}" \
    > "${1/.txt/.munged.txt}"

    #  Remove the temporary files
    rm "${1/.txt/-1.txt}" \
    "${1/.txt/-2.txt}" \
    "${1/.txt/-3.txt}" \
    "${1/.txt/-4.txt}"

    #  Add headers to the pertinent files
    awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
    "${1/.txt/.munged.txt}" \
    > /tmp/out && mv /tmp/out "${1/.txt/.munged.txt}"
}

mungePairwise "${cigarAlt}"
mungePairwise "${cigarAmbig}"
mungePairwise "${cigarContra}"
mungePairwise "${cigarRef}"

cat "${cigarAlt/.txt/.munged.txt}" | wc -l  # 32581
cat "${cigarAmbig/.txt/.munged.txt}" | wc -l  # 172125
cat "${cigarContra/.txt/.munged.txt}" | wc -l  # 29
cat "${cigarRef/.txt/.munged.txt}" | wc -l  # 29951

#  Remove initial sam2pairwise output files (not useful unless they're munged)
find . -type f -name "*.sam2pairwise.txt" -print0 | xargs -0 rm

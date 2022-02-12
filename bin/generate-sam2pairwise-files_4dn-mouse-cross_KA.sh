#!/bin/bash

# generate-sam2pairwise-files_4dn-mouse-cross-tests.sh
# KA

script="generate-sam2pairwise-files_4dn-mouse-cross-tests.sh"
# echo "${script}"

#  Working with the following three files for now:
#+ - 129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam
#+ - CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam
#+ - mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam

#  Important dependency: sam2pairwise
#+ 
#+ Can install via conda thanks to the Beliveau Lab:
#+ conda install -c beliveau-lab sam2pairwise
#+
#+ More info: github.com/mlafave/sam2pairwise

#  Another important dependency: samtools
#+
#+ Can also get it from conda

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

directory="/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
cd "${directory}" ||
    {
        echo "Exiting: Directory not found."
        return 1 2> /dev/null
        exit 1
    }


#  Convert to extended CIGAR format -------------------------------------------
reformatCIGARstring() {
    "${reformat}" in="${1}" out="${2}" sam=1.4 2> "${2/.bam/.err}"
}

# init129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
# initCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"
# initNmasked="mm10-CAST-129-Nmasked.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.bam"

init129="processed.mate-paired.mm10.chrX.bam"
initCAST="processed.mate-paired.CAST.chrX.bam"
initNmasked="processed.mate-paired.129S1.chrX.bam"

# bam129="${init129/.bam/.extendedCIGAR.bam}"
# bamCAST="${initCAST/.bam/.extendedCIGAR.bam}"
# bamNmasked="${initNmasked/.bam/.extendedCIGAR.bam}"

# reformatCIGARstring "${init129}" "${bam129}"
# reformatCIGARstring "${initCAST}" "${bamCAST}"
# reformatCIGARstring "${initNmasked}" "${bamNmasked}"

# #  Check
# samtools view "${init129}" | head
# echo ""
# samtools view "${bam129}" | head

parallel samtools index ::: "${bam129}" "${bamCAST}" "${bamNmasked}"


#  Generate CIGAR constructions with sam2pairwise -----------------------------
# s2p="/Users/kalavattam/sam2pairwise/src/sam2pairwise"
# echo "s2p is ${s2p}"
generateCIGARconstructions() {
    samtools view "${1}" | sam2pairwise > "${1/.bam/.sam2pairwise.txt}"
}

bam129="processed.mate-paired.mm10.chrX.bam"
bamCAST="processed.mate-paired.CAST.chrX.bam"
bamNmasked="processed.mate-paired.129S1.chrX.bam"

generateCIGARconstructions "${bam129}"
generateCIGARconstructions "${bamCAST}"
generateCIGARconstructions "${bamNmasked}"

cigar129="${bam129/.bam/.sam2pairwise.txt}"
cigarCAST="${bamCAST/.bam/.sam2pairwise.txt}"
cigarNmasked="${bamNmasked/.bam/.sam2pairwise.txt}"


#  Replace mismatch-spaces with mismatch-asterisks ----------------------------
replaceSpacesAsterisks() {
    awk '{gsub(" ", "*", $0); print}' "${1}" \
    > "${1/.txt/.tmp.txt}" && mv "${1/.txt/.tmp.txt}" "${1}"
}

replaceSpacesAsterisks "${cigar129}"
replaceSpacesAsterisks "${cigarCAST}"
replaceSpacesAsterisks "${cigarNmasked}"


# #  Checks
# #  129
# wc129=$(cat "${cigar129}" | wc -l)
# echo "${wc129}"  # 1573032  # 937352
# echo $(( wc129 / 4 ))  # 393258  # 234338
# 
# wcBam129=$(samtools view "${bam129}" | wc -l)
# echo "${wcBam129}"  # 393258  # 234338
# 
# #  CAST
# wcCAST=$(cat "${cigarCAST}" | wc -l)
# echo "${wcCAST}"  # 1573032  # 868864
# echo $(( wcCAST / 4 ))  # 393258  # 217216
# 
# wcBamCAST=$(samtools view "${bamCAST}" | wc -l)
# echo "${wcBamCAST}"  # 393258  # 217216
# 
# #  Nmasked
# wcNmasked=$(cat "${cigarNmasked}" | wc -l)
# echo "${wcNmasked}"  # 1636848  # 937352
# echo $(( wcNmasked / 4 ))  # 409212  # 234338
# 
# wcBamNmasked=$(samtools view "${bamNmasked}" | wc -l)
# echo "${wcBamNmasked}"  # 409212  # 234338
# 
# #NOTE It all seems to check out...


#  Munge the output of sam2pairwise: 129 --------------------------------------
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

mungePairwise "${cigar129}"
mungePairwise "${cigarCAST}"
mungePairwise "${cigarNmasked}"

cat "${cigar129/.txt/.munged.txt}" | wc -l  # 393258  # 234339
cat "${cigarCAST/.txt/.munged.txt}" | wc -l  # 403192  # 217217
cat "${cigarNmasked/.txt/.munged.txt}" | wc -l  # 409212  # 234339


###############################################################################
directory="/Users/kalavattam/Dropbox/My Mac (Kriss-MacBook-Pro.local)/Downloads/to-do/get_unique_fragments/Bonora"
cd "${directory}" || 
    {
        echo "Directory not found; look into this."
    }

#  Working with the following three files for now:
#+ - 129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.liftOverMm10.sorted.bam
#+ - CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.liftOverMm10.sorted.bam

#  Important dependency: sam2pairwise
#+ 
#+ Can install via conda thanks to the Beliveau Lab:
#+ conda install -c beliveau-lab sam2pairwise
#+
#+ More info: github.com/mlafave/sam2pairwise

#  Another important dependency: samtools
#+
#+ Can also get it from conda

#TODO Regenerate the below .bam files using CIGAR-processed .bam files (above); use the script 'crossmap-bams.sh'
#TODO Need to work with processed.mate-paired.CAST.chrX.bam, processed.mate-paired.129S1.chrX.bam for liftOverMm10, not the initial, unprocessed files...

bam129="129S1-SvImJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.liftOverMm10.sorted.bam"
cigar129="${bam129/.bam/.sam2pairwise.txt}"
bamCAST="CAST-EiJ.F121-6-CASTx129.undifferentiated.dedup.MarkDuplicates.sort.chrX.rmdup.liftOverMm10.sorted.bam"
cigarCAST="${bamCAST/.bam/.sam2pairwise.txt}"

echo "${bam129}"
echo "${bamCAST}"

#  Generate CIGAR constructions with sam2pairwise
samtools view "${bam129}" | sam2pairwise > "${cigar129}"
samtools view "${bamCAST}" | sam2pairwise > "${cigarCAST}"

#  Replace mismatch-spaces with mismatch-asterisks
awk '{gsub(" ", "*", $0); print}' "${cigar129}" > "${cigar129/.txt/.tmp.txt}" && mv "${cigar129/.txt/.tmp.txt}" "${cigar129}"
awk '{gsub(" ", "*", $0); print}' "${cigarCAST}" > "${cigarCAST/.txt/.tmp.txt}" && mv "${cigarCAST/.txt/.tmp.txt}" "${cigarCAST}"

#  129
wc129=$(cat "${cigar129}" | wc -l)
echo "${wc129}"  # 1573032
echo $((${wc129} / 4))  # 393258

wcBam129=$(samtools view "${bam129}" | wc -l)
echo "${wcBam129}"  # 393258

#  CAST
wcCAST=$(cat "${cigarCAST}" | wc -l)
echo "${wcCAST}"  # 1612768
echo $((${wcCAST} / 4))  # 403192

wcBamCAST=$(samtools view "${bamCAST}" | wc -l)
echo "${wcBamCAST}"  # 403192

#  Munge the output of sam2pairwise: 129 --------------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigar129}" | awk 'NR%4==1' > "${cigar129/.txt/-1.txt}"
cat "${cigar129}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigar129/.txt/-2.txt}"
cat "${cigar129}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigar129/.txt/-3.txt}"
cat "${cigar129}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigar129/.txt/-4.txt}"

paste "${cigar129/.txt/-1.txt}" "${cigar129/.txt/-2.txt}" "${cigar129/.txt/-3.txt}" "${cigar129/.txt/-4.txt}" > "${cigar129/.txt/.munged.txt}"
rm "${cigar129/.txt/-1.txt}" "${cigar129/.txt/-2.txt}" "${cigar129/.txt/-3.txt}" "${cigar129/.txt/-4.txt}"
cat "${cigar129/.txt/.munged.txt}" | wc -l  # 393258


#  Munge the output of sam2pairwise: CAST -------------------------------------

#  Want every second, third, and fourth row to be additional tabbed entries on the first row
cat "${cigarCAST}" | awk 'NR%4==1' > "${cigarCAST/.txt/-1.txt}"
cat "${cigarCAST}" | awk 'NR > 1{ print }' | awk 'NR%4==1' > "${cigarCAST/.txt/-2.txt}"
cat "${cigarCAST}" | awk 'NR > 2{ print }' | awk 'NR%4==1' > "${cigarCAST/.txt/-3.txt}"
cat "${cigarCAST}" | awk 'NR > 3{ print }' | awk 'NR%4==1' > "${cigarCAST/.txt/-4.txt}"

paste "${cigarCAST/.txt/-1.txt}" "${cigarCAST/.txt/-2.txt}" "${cigarCAST/.txt/-3.txt}" "${cigarCAST/.txt/-4.txt}" > "${cigarCAST/.txt/.munged.txt}"
rm "${cigarCAST/.txt/-1.txt}" "${cigarCAST/.txt/-2.txt}" "${cigarCAST/.txt/-3.txt}" "${cigarCAST/.txt/-4.txt}"
cat "${cigarCAST/.txt/.munged.txt}" | wc -l  # 403192


#  Add headers to the pertinent files -----------------------------------------
awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigar129/.txt/.munged.txt}" \
> /tmp/out && mv /tmp/out "${cigar129/.txt/.munged.txt}"

awk 'BEGIN { print "qname\tflag\trname\tpos\tmapq\tcigar\tmrnm\tmpos\tisize\tread_sequence\tmatches\treference_sequence" } { print }' \
"${cigarCAST/.txt/.munged.txt}" \
> /tmp/out && mv /tmp/out "${cigarCAST/.txt/.munged.txt}"

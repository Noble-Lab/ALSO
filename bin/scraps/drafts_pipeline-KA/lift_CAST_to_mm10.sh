#!/bin/bash

#  lift_CAST_to_mm10.sh
#  KA


# -------
#  Start recording time
start="$(date +%s)"

# -------
#  Assign arguments, set up functions
mate="${1:-"pos"}"
chromosome="${2:-"chrX"}"


displaySpinningIcon() {
    #  Display \"spinning icon\" while background process runs
    spin='-\\|/'
    i=0
    while kill -0 "${1}" 2> /dev/null; do
        i=$(( (i + 1) % 4 ))
        printf "\r${spin:$i:1} %s " "${2}"
        sleep .1
    done
}


renameChrCommonCAST() {
    #  Convert CAST-EiJ chromosomes names from common to official
    cat "${1}" \
    | sed 's/chr1/CM003994.1/g' \
    | sed 's/chr2/CM003995.1/g' \
    | sed 's/chr3/CM003996.1/g' \
    | sed 's/chr4/CM003997.1/g' \
    | sed 's/chr5/CM003998.1/g' \
    | sed 's/chr6/CM003999.1/g' \
    | sed 's/chr7/CM004000.1/g' \
    | sed 's/chr8/CM004001.1/g' \
    | sed 's/chr9/CM004002.1/g' \
    | sed 's/chr10/CM004003.1/g' \
    | sed 's/chr11/CM004004.1/g' \
    | sed 's/chr12/CM004005.1/g' \
    | sed 's/chr13/CM004006.1/g' \
    | sed 's/chr14/CM004007.1/g' \
    | sed 's/chr15/CM004008.1/g' \
    | sed 's/chr16/CM004009.1/g' \
    | sed 's/chr17/CM004010.1/g' \
    | sed 's/chr18/CM004011.1/g' \
    | sed 's/chr19/CM004012.1/g' \
    | sed 's/chrX/CM004013.1/g' \
    > "${2}"
}


renameChrCASTCommon() {
    #  Convert CAST-EiJ chromosome names from official to common
    cat "${1}" \
    | sed 's/CM003994.1/chr1/g' \
    | sed 's/CM003995.1/chr2/g' \
    | sed 's/CM003996.1/chr3/g' \
    | sed 's/CM003997.1/chr4/g' \
    | sed 's/CM003998.1/chr5/g' \
    | sed 's/CM003999.1/chr6/g' \
    | sed 's/CM004000.1/chr7/g' \
    | sed 's/CM004001.1/chr8/g' \
    | sed 's/CM004002.1/chr9/g' \
    | sed 's/CM004003.1/chr10/g' \
    | sed 's/CM004004.1/chr11/g' \
    | sed 's/CM004005.1/chr12/g' \
    | sed 's/CM004006.1/chr13/g' \
    | sed 's/CM004007.1/chr14/g' \
    | sed 's/CM004008.1/chr15/g' \
    | sed 's/CM004009.1/chr16/g' \
    | sed 's/CM004010.1/chr17/g' \
    | sed 's/CM004011.1/chr18/g' \
    | sed 's/CM004012.1/chr19/g' \
    | sed 's/CM004013.1/chrX/g' \
    > "${2}"
}


# -------
#  Set up paths, files
path_in="${3:-"/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0"}"
path_out="${4:-"/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0"}"

path_liftOver="${5:-"/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver"}"

from_R=$(find . -maxdepth 1 -type f -name "tbl.CAST.${mate}.${chromosome}.bed" -printf "%f\n")
prefix_from_R=${from_R%????}

# -------
#  Run liftOver loop
for file_prefix in "${prefix_from_R[@]}"; do
    #  Declare variables
    file_in="${file_prefix}.bed"
    file_in_rename="${file_prefix}.rename.bed"

    file_out_lifted="${file_prefix}.lifted.bed"
    file_out_lifted_rename="${file_prefix}.lifted.rename.bed"

    file_out_unlifted="${file_prefix}.unlifted.bed"
    file_out_unlifted_rename="${file_prefix}.unlifted.rename.bed"

    file_liftOver="CAST-EiJ-to-mm10.over.chain.gz"

    #  Variables for liftOver
    in="${path_in}/${file_in_rename}"
    chain="${path_liftOver}/${file_liftOver}"
    out_lift="${path_out}/${file_out_lifted}"
    out_unlift="${path_out}/${file_out_unlifted}"

    # -------
    #  Convert infile chromosomes names from common to official
    renameChrCommonCAST "${file_in}" "${file_in_rename}"
    # head ${file_in_rename}

    # -------
    #  Do the liftOver in the background, saving those regions that are and are
    #+ not lifted
    liftOver -bedPlus=3 -tab "${in}" "${chain}" "${out_lift}" "${out_unlift}" # &

    #  Display a spinning icon while liftOver is taking place
    # displaySpinningIcon $! "Lifting over $(basename "${in}")"

    # -------
    #  Convert the outfile chromosome names from official to common
    renameChrCASTCommon "${file_out_lifted}" "${file_out_lifted_rename}" # &
    # displaySpinningIcon $! "Renaming $(basename "${file_out_lifted}")"

    renameChrCASTCommon "${file_out_unlifted}" "${file_out_unlifted_rename}" # &
    # displaySpinningIcon $! "Renaming $(basename "${file_out_unlifted}")"

    # -------
    #  Clean up
    # rm "${file_in}"
    rm "${file_in_rename}"
    rm "${file_out_lifted}"
    rm "${file_out_unlifted}"
    
    mv "${file_out_lifted_rename}" "${file_out_lifted}"
    mv "${file_out_unlifted_rename}" "${file_out_unlifted}"
done

# -------
#  Record end time
end="$(date +%s)"

#  Echo run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "Run time: ${run_time} seconds"
echo ""

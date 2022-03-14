#!/bin/bash

#  lift_129_to_mm10.sh
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


renameChrCommon129() {
    #  Convert 129S1-SvImJ chromosomes names from common to official
    cat "${1}" \
    | sed 's/chr1/CM003934.1/g' \
    | sed 's/chr2/CM003935.1/g' \
    | sed 's/chr3/CM003936.1/g' \
    | sed 's/chr4/CM003937.1/g' \
    | sed 's/chr5/CM003938.1/g' \
    | sed 's/chr6/CM003939.1/g' \
    | sed 's/chr7/CM003940.1/g' \
    | sed 's/chr8/CM003941.1/g' \
    | sed 's/chr9/CM003942.1/g' \
    | sed 's/chr10/CM003943.1/g' \
    | sed 's/chr11/CM003944.1/g' \
    | sed 's/chr12/CM003945.1/g' \
    | sed 's/chr13/CM003946.1/g' \
    | sed 's/chr14/CM003947.1/g' \
    | sed 's/chr15/CM003948.1/g' \
    | sed 's/chr16/CM003949.1/g' \
    | sed 's/chr17/CM003950.1/g' \
    | sed 's/chr18/CM003951.1/g' \
    | sed 's/chr19/CM003952.1/g' \
    | sed 's/chrX/CM003953.1/g' \
    > "${2}"
}


renameChr129Common() {
    #  Convert 129S1-SvImJ chromosome names from official to common
    cat "${1}" \
    | sed 's/CM003934.1/chr1/g' \
    | sed 's/CM003935.1/chr2/g' \
    | sed 's/CM003936.1/chr3/g' \
    | sed 's/CM003937.1/chr4/g' \
    | sed 's/CM003938.1/chr5/g' \
    | sed 's/CM003939.1/chr6/g' \
    | sed 's/CM003940.1/chr7/g' \
    | sed 's/CM003941.1/chr8/g' \
    | sed 's/CM003942.1/chr9/g' \
    | sed 's/CM003943.1/chr10/g' \
    | sed 's/CM003944.1/chr11/g' \
    | sed 's/CM003945.1/chr12/g' \
    | sed 's/CM003946.1/chr13/g' \
    | sed 's/CM003947.1/chr14/g' \
    | sed 's/CM003948.1/chr15/g' \
    | sed 's/CM003949.1/chr16/g' \
    | sed 's/CM003950.1/chr17/g' \
    | sed 's/CM003951.1/chr18/g' \
    | sed 's/CM003952.1/chr19/g' \
    | sed 's/CM003953.1/chrX/g' \
    > "${2}"
}


# -------
#  Set up paths, files
path_in="${3:-"/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0"}"
path_out="${4:-"/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/kga0"}"

path_liftOver="${5:-"/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver"}"

from_R=$(find . -maxdepth 1 -type f -name "tbl.129.${mate}.${chromosome}.bed" -printf "%f\n")
prefix_from_R=${from_R%????}

# -------
#  Run liftOver loop
for file_prefix in "${prefix_from_R[@]}"; do
    #  Declare variables
    file_in="${file_prefix}.bed"
    file_in_rename="${file_prefix}.bed"

    file_out_lifted="${file_prefix}.lifted.bed"
    file_out_lifted_rename="${file_prefix}.lifted.bed"

    file_out_unlifted="${file_prefix}.unlifted.bed"
    file_out_unlifted_rename="${file_prefix}.unlifted.bed"

    file_liftOver="129S1-SvImJ-to-mm10.over.chain.gz"

    #  Variables for liftOver
    in="${path_in}/${file_in_rename}"
    chain="${path_liftOver}/${file_liftOver}"
    out_lift="${path_out}/${file_out_lifted}"
    out_unlift="${path_out}/${file_out_unlifted}"

    # -------
    #  Convert infile chromosomes names from common to official
    renameChrCommon129 "${file_in}" "${file_in_rename}"
    # head ${file_in_rename}

    # -------
    #  Do the liftOver in the background, saving those regions that are and are
    #+ not lifted
    liftOver -bedPlus=3 -tab "${in}" "${chain}" "${out_lift}" "${out_unlift}" &

    #  Display a spinning icon while liftOver is taking place
    displaySpinningIcon $! "Lifting over $(basename "${in}")"

    # -------
    #  Convert the outfile chromosome names from official to common
    renameChr129Common "${file_out_lifted}" "${file_out_lifted_rename}" &
    displaySpinningIcon $! "Renaming $(basename "${file_out_lifted}")"

    renameChr129Common "${file_out_unlifted}" "${file_out_unlifted_rename}" &
    displaySpinningIcon $! "Renaming $(basename "${file_out_unlifted}")"

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

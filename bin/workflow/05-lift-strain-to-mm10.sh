#!/bin/bash

#  lift-strain-to-mm10.sh
#  KA


#  Start recording time -------------------------------------------------------
start="$(date +%s)"


#  Functions ------------------------------------------------------------------
checkDependency() {
    #  Check if program is available in "${PATH}"; exit if not
    command -v "${1}" &>/dev/null ||
        {
            echo "Exiting: ${1} not found. Install ${1}."
            exit 1
        }
}


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


renameChrCommon129() {
    #  Convert 129S1-SvImJ chromosomes names from common to official
    # shellcheck disable=SC2002
    cat "${1}" \
    | sed 's/chr1\t/CM003934.1\t/g' \
    | sed 's/chr2\t/CM003935.1\t/g' \
    | sed 's/chr3\t/CM003936.1\t/g' \
    | sed 's/chr4\t/CM003937.1\t/g' \
    | sed 's/chr5\t/CM003938.1\t/g' \
    | sed 's/chr6\t/CM003939.1\t/g' \
    | sed 's/chr7\t/CM003940.1\t/g' \
    | sed 's/chr8\t/CM003941.1\t/g' \
    | sed 's/chr9\t/CM003942.1\t/g' \
    | sed 's/chr10\t/CM003943.1\t/g' \
    | sed 's/chr11\t/CM003944.1\t/g' \
    | sed 's/chr12\t/CM003945.1\t/g' \
    | sed 's/chr13\t/CM003946.1\t/g' \
    | sed 's/chr14\t/CM003947.1\t/g' \
    | sed 's/chr15\t/CM003948.1\t/g' \
    | sed 's/chr16\t/CM003949.1\t/g' \
    | sed 's/chr17\t/CM003950.1\t/g' \
    | sed 's/chr18\t/CM003951.1\t/g' \
    | sed 's/chr19\t/CM003952.1\t/g' \
    | sed 's/chrX\t/CM003953.1\t/g' \
    > "${2}"
}


renameChr129Common() {
    #  Convert 129S1-SvImJ chromosome names from official to common
    # shellcheck disable=SC2002
    cat "${1}" \
    | sed 's/CM003934.1\t/chr1\t/g' \
    | sed 's/CM003935.1\t/chr2\t/g' \
    | sed 's/CM003936.1\t/chr3\t/g' \
    | sed 's/CM003937.1\t/chr4\t/g' \
    | sed 's/CM003938.1\t/chr5\t/g' \
    | sed 's/CM003939.1\t/chr6\t/g' \
    | sed 's/CM003940.1\t/chr7\t/g' \
    | sed 's/CM003941.1\t/chr8\t/g' \
    | sed 's/CM003942.1\t/chr9\t/g' \
    | sed 's/CM003943.1\t/chr10\t/g' \
    | sed 's/CM003944.1\t/chr11\t/g' \
    | sed 's/CM003945.1\t/chr12\t/g' \
    | sed 's/CM003946.1\t/chr13\t/g' \
    | sed 's/CM003947.1\t/chr14\t/g' \
    | sed 's/CM003948.1\t/chr15\t/g' \
    | sed 's/CM003949.1\t/chr16\t/g' \
    | sed 's/CM003950.1\t/chr17\t/g' \
    | sed 's/CM003951.1\t/chr18\t/g' \
    | sed 's/CM003952.1\t/chr19\t/g' \
    | sed 's/CM003953.1\t/chrX\t/g' \
    > "${2}"
}


renameChrCommonCAST() {
    #  Convert CAST-EiJ chromosomes names from common to official
    # shellcheck disable=SC2002
    cat "${1}" \
    | sed 's/chr1\t/CM003994.1\t/g' \
    | sed 's/chr2\t/CM003995.1\t/g' \
    | sed 's/chr3\t/CM003996.1\t/g' \
    | sed 's/chr4\t/CM003997.1\t/g' \
    | sed 's/chr5\t/CM003998.1\t/g' \
    | sed 's/chr6\t/CM003999.1\t/g' \
    | sed 's/chr7\t/CM004000.1\t/g' \
    | sed 's/chr8\t/CM004001.1\t/g' \
    | sed 's/chr9\t/CM004002.1\t/g' \
    | sed 's/chr10\t/CM004003.1\t/g' \
    | sed 's/chr11\t/CM004004.1\t/g' \
    | sed 's/chr12\t/CM004005.1\t/g' \
    | sed 's/chr13\t/CM004006.1\t/g' \
    | sed 's/chr14\t/CM004007.1\t/g' \
    | sed 's/chr15\t/CM004008.1\t/g' \
    | sed 's/chr16\t/CM004009.1\t/g' \
    | sed 's/chr17\t/CM004010.1\t/g' \
    | sed 's/chr18\t/CM004011.1\t/g' \
    | sed 's/chr19\t/CM004012.1\t/g' \
    | sed 's/chrX\t/CM004013.1\t/g' \
    > "${2}"
}


renameChrCASTCommon() {
    #  Convert CAST-EiJ chromosome names from official to common
    # shellcheck disable=SC2002
    cat "${1}" \
    | sed 's/CM003994.1\t/chr1\t/g' \
    | sed 's/CM003995.1\t/chr2\t/g' \
    | sed 's/CM003996.1\t/chr3\t/g' \
    | sed 's/CM003997.1\t/chr4\t/g' \
    | sed 's/CM003998.1\t/chr5\t/g' \
    | sed 's/CM003999.1\t/chr6\t/g' \
    | sed 's/CM004000.1\t/chr7\t/g' \
    | sed 's/CM004001.1\t/chr8\t/g' \
    | sed 's/CM004002.1\t/chr9\t/g' \
    | sed 's/CM004003.1\t/chr10\t/g' \
    | sed 's/CM004004.1\t/chr11\t/g' \
    | sed 's/CM004005.1\t/chr12\t/g' \
    | sed 's/CM004006.1\t/chr13\t/g' \
    | sed 's/CM004007.1\t/chr14\t/g' \
    | sed 's/CM004008.1\t/chr15\t/g' \
    | sed 's/CM004009.1\t/chr16\t/g' \
    | sed 's/CM004010.1\t/chr17\t/g' \
    | sed 's/CM004011.1\t/chr18\t/g' \
    | sed 's/CM004012.1\t/chr19\t/g' \
    | sed 's/CM004013.1\t/chrX\t/g' \
    > "${2}"
}


#  Handle arguments, assign variables -----------------------------------------
printUsage() {
    echo ""
    echo "${0}:"
    echo "Take a \"pos\" or \"mpos\" bed file from having run \"04-split-"
    echo "index-repair-bam.sh\" and \"lift\" its coordinates over from the"
    echo "initial alignment strain coordinates (e.g., \"CAST-EiJ\""
    echo "coordinates) to \"mm10\" coordinates."
    echo ""
    echo ""
    echo "Dependencies:"
    echo "- liftOver >= 366 (untested w/previous versions)"
    echo ""
    echo ""
    echo "Arguments:"
    echo "-h <print this help message and exit>"
    echo "-u <use safe mode: \"TRUE\" or \"FALSE\" (logical)>"
    echo "-i <bed infile, including path (chr)>"
    echo "-o <path for \"lifted\" bed outfiles (chr)>"
    echo "-s <strain for performing liftOver of bed files; currently available"
    echo "    options:"
    echo "    - \"CAST-EiJ\", \"CAST\", or \"C\" for \"CAST-EiJ\""
    echo "    - \"129S1-SvImJ\", \"129\", or \"1\" for \"129S1-SvImJ\""
    echo "    - \"CAROLI-EiJ\", \"CAROLI\", \"Ryukyu\" or \"R\" for \"CAROLI-EiJ\""
    echo "    - \"SPRET-EiJ\", \"SPRET\", or \"S\" for \"SPRET-EiJ>\""
    echo "-c <gzipped liftOver chain file for strain, including path (chr);"
    echo "    note: for liftOver to work, the liftOver strain chain should"
    echo "    match the strain set in argument \"-s\">"
    exit
}

while getopts "h:u:i:o:s:c:" opt; do
    case "${opt}" in
        h) printUsage ;;
        u) safe_mode="${OPTARG}" ;;
        i) infile="${OPTARG}" ;;
        o) outpath="${OPTARG}" ;;
        s) strain="${OPTARG}" ;;
        c) chain="${OPTARG}" ;;
        *) printUsage ;;
    esac
done

[[ -z "${safe_mode}" ]] && printUsage
[[ -z "${infile}" ]] && printUsage
[[ -z "${outpath}" ]] && printUsage
[[ -z "${strain}" ]] && printUsage
[[ -z "${chain}" ]] && printUsage

# #  Test defaults
# safe_mode="FALSE"
# infile="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/2022-0320_test_04-05/test.300000.chr19.pos.bed"
# outpath="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/2022-0320_test_04-05"
# # strain="129S1-SvImJ"
# # chain="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_chain/129S1-SvImJ-to-mm10.over.chain.gz"
# # strain="CAST-EiJ"
# # chain="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_chain/CAST-EiJ-to-mm10.over.chain.gz"
# # strain="CAROLI-EiJ"
# # chain="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_chain/CAROLI-EiJ-to-mm10.over.chain.gz"
# strain="SPRET-EiJ"
# chain="/Users/kalavattam/Dropbox/UW/projects-etc/2021_kga0_4dn-mouse-cross/data/files_chain/SPRET-EiJ-to-mm10.over.chain.gz"
#
# bash bin/lift-strain-to-mm10.sh \
# -u "${safe_mode}" \
# -i "${infile}" \
# -o "${outpath}" \
# -s "${strain}" \
# -c "${chain}"


#  Check, establish variable assignments --------------------------------------
echo -e ""
echo -e "Running ${0}..."

#  Check for necessary dependencies; exit if not found
checkDependency liftOver

#  Evaluate "${safe_mode}"
case "$(echo "${safe_mode}" | tr '[:upper:]' '[:lower:]')" in
    true | t) \
        echo -e "-u: Safe mode is on."
        set -Eeuxo pipefail
        ;;
    false | f) \
        echo -e "-u: Safe mode is off."
        :
        ;;
    *) \
        echo -e "Exiting: -u safe-mode argument must be \"TRUE\" or \"FALSE\".\n"
        exit 1
        ;;
esac

#  Check that "${infile}" exists
[[ -f "${infile}" ]] ||
    {
        echo -e "Exiting: -i ${infile} does not exist.\n"
        exit 1
    }

#  Make "${outpath}" if it doesn't exist
[[ -d "${outpath}" ]] ||
    {
        echo "-o: Directory ${outpath} does not exist; making the directory."
        mkdir -p "${outpath}"
    }

#  Evaluate "${strain}"
case "$(echo "${strain}" | tr '[:upper:]' '[:lower:]')" in
    cast-eij | cast | c) \
        flag_strain=1
        name="CAST-EiJ"
        ;;
    129s1-svimj | 129 | 1) \
        flag_strain=2
        name="129S1-SvImJ"
        ;;
    caroli-eij | caroli | ryukyu | r) \
        flag_strain=3
        name="CAROLI-EiJ"
        ;;
    spret-eij | spret | s) \
        flag_strain=4
        name="SPRET-EiJ"
        ;;
    *) \
        echo "Exiting: Script is not set up for ${strain} as -s strain"
        echo "argument; -s strain argument must be \"CAST-EiJ\","
        echo "\"129S1-SvImJ\", \"CAROLI-EiJ\", or \"SPRET-EiJ.\""
        exit 1
        ;;
esac

echo -e "-s: Strain is set to \"${name}.\""

#  Check that "${chain}" exists
[[ -f "${chain}" ]] ||
    {
        echo -e "Exiting: -c ${chain} does not exist.\n"
        exit 1
    }

echo -e "-c: Chain file is $(basename "${chain}").\n"

#  Set up tmp.bed file, which is needed for subsequent code
infile_tmp="${infile%/*}/$(basename "${infile/.bed/.${name}.tmp.bed}")"

# echo "${infile}"
# echo "${infile_tmp}"
export infile
export infile_tmp


#  Convert infile chromosomes names from common to official -------------------
if [[ $((flag_strain)) -eq 1 ]]; then
    renameChrCommonCAST "${infile}" "${infile_tmp}" &
    displaySpinningIcon $! "Converting names for $(basename "${infile}")... "
elif [[ $((flag_strain)) -eq 2 ]]; then
    renameChrCommon129 "${infile}" "${infile_tmp}" &
    displaySpinningIcon $! "Converting names for $(basename "${infile}")... "
elif [[ $((flag_strain)) -eq 3 || $((flag_strain)) -eq 4 ]]; then
    cp "${infile}" "${infile_tmp}"
else
    echo -e "Exiting: An error occurred when setting the strain flag.\n"
    exit 1
fi
# echo "\${flag_strain} is $((flag_strain))"


#  Perform liftOvers, reporting regions that are and are not lifted -----------
liftOver -bedPlus=3 -tab \
"${infile_tmp}" \
"${chain}" \
"${infile_tmp/.bed/.lifted.bed}" \
"${infile_tmp/.bed/.unlifted.bed}" &
displaySpinningIcon $! \
"Lifting over chromosome positions for $(basename "${infile_tmp}"): "


#  Convert the outfile chromosome names from official to common ---------------
#  lifted.bed file
if [[ $((flag_strain)) -eq 1 ]]; then
    renameChrCASTCommon \
    "${infile_tmp/.bed/.lifted.bed}" "${infile/.bed/.${name}.lifted.bed}" &
    displaySpinningIcon $! \
    "Converting names for $(basename "${infile/.bed/.${name}.lifted.bed}")... "
elif [[ $((flag_strain)) -eq 2 ]]; then
    renameChr129Common \
    "${infile_tmp/.bed/.lifted.bed}" "${infile/.bed/.${name}.lifted.bed}" &
    displaySpinningIcon $! \
    "Converting names for $(basename "${infile/.bed/.${name}.lifted.bed}")... "
elif [[ $((flag_strain)) -eq 3 || $((flag_strain)) -eq 4 ]]; then
    cp "${infile_tmp/.bed/.lifted.bed}" "${infile/.bed/.${name}.lifted.bed}"
else
    echo -e "Exiting: An error occurred when setting the strain flag.\n"
    exit 1
fi

#  unlifted.bed file
if [[ $((flag_strain)) -eq 1 ]]; then
    renameChrCASTCommon \
    "${infile_tmp/.bed/.unlifted.bed}" "${infile/.bed/.${name}.unlifted.bed}" &
    displaySpinningIcon $! \
    "Converting names for $(basename "${infile/.bed/.${name}.unlifted.bed}")... "
elif [[ $((flag_strain)) -eq 2 ]]; then
    renameChr129Common \
    "${infile_tmp/.bed/.unlifted.bed}" "${infile/.bed/.${name}.unlifted.bed}" &
    displaySpinningIcon $! \
    "Converting names for $(basename "${infile/.bed/.${name}.unlifted.bed}")... "
elif [[ $((flag_strain)) -eq 3 || $((flag_strain)) -eq 4 ]]; then
    cp "${infile_tmp/.bed/.unlifted.bed}" "${infile/.bed/.${name}.unlifted.bed}"
else
    echo -e "Exiting: An error occurred when setting the strain flag.\n"
    exit 1
fi

#  Remove unneeded intermediate files
rm \
"${infile_tmp}" \
"${infile_tmp/.bed/.lifted.bed}" \
"${infile_tmp/.bed/.unlifted.bed}" &
displaySpinningIcon $! "Removing temporary bed files... "


#  Munge the lifted and unlifted bed files ------------------------------------
#  Munge the lifted.bed file
# shellcheck disable=SC2002
cat "${infile/.bed/.${name}.lifted.bed}" \
| awk 'BEGIN{FS=OFS="\t"} {print $0 OFS "Liftover successful"}' \
> "${infile/.bed/.${name}.lifted.tmp.bed}" &
displaySpinningIcon $! "Munging $(basename "${infile/.bed/.${name}.lifted.tmp.bed}")... "

#  Munge the unlifted.bed file through the following steps:
#+ - Append tab to the end of each line
#+ - Move every odd row down to the beginning of the subsequent even row
#+ - Strip "#" from the beginning of each line
#+ - Move column 1 to column 5
#+ - Strip the initial space from the beginning of each line
# shellcheck disable=SC2002
cat "${infile/.bed/.${name}.unlifted.bed}" \
| sed 's/$/\t/' \
| awk 'NR%2==0{print p,$0}{p=$0}' \
| sed 's/^#//' \
| awk 'BEGIN{FS=OFS="\t"} {print $2, $3, $4, $5, $1}' \
| sed 's/^ //' \
> "${infile/.bed/.${name}.unlifted.tmp.bed}" &
displaySpinningIcon $! "Munging $(basename "${infile/.bed/.${name}.unlifted.tmp.bed}")... "


#  Concatenate lifted and unlifted bed files ----------------------------------
mv -f "${infile/.bed/.${name}.lifted.tmp.bed}" "${infile/.bed/.${name}.lifted.bed}"
mv -f "${infile/.bed/.${name}.unlifted.tmp.bed}" "${infile/.bed/.${name}.unlifted.bed}"

cat "${infile/.bed/.${name}.lifted.bed}" "${infile/.bed/.${name}.unlifted.bed}" \
> "${infile/.bed/.liftOver.${name}.bed}" &
displaySpinningIcon $! "Concatenating lifted and unlifted bed files... "

rm \
"${infile/.bed/.${name}.lifted.bed}" \
"${infile/.bed/.${name}.unlifted.bed}" &
displaySpinningIcon $! "Removing lifted and unlifted bed files... "


#  End recording time ---------------------------------------------------------
end="$(date +%s)"

#  Return run time
run_time="$(echo "${end}" - "${start}" | bc -l)"
echo ""
echo "${0} run time: ${run_time} seconds."
echo ""

exit 0

#!/bin/bash

path="/Users/kalavattam/Downloads/to-do/2021-1105-1107/liftOver"
cd "${path}" ||
    {
        echo "Failure to cd to ${path}. Exiting."
        exit 1
    }

# -----------------------------------------------------------------------------
files=(
    /net/noble/vol1/home/kga0/genomes/Ensembl.129S1-SvImJ/liftOver/129S1-SvImJ-to-mm10.over.chain.gz
    /net/noble/vol1/home/kga0/genomes/Ensembl.129S1-SvImJ/liftOver/mm10-to-129S1-SvImJ.over.chain.gz
    /net/noble/vol1/home/kga0/genomes/Ensembl.CAST-EiJ/liftOver/CAST-EiJ-to-mm10.over.chain.gz
    /net/noble/vol1/home/kga0/genomes/Ensembl.CAST-EiJ/liftOver/mm10-to-CAST-EiJ.over.chain.gz
)

for i in "${files[@]}"; do
    echo "Working with ${i}"
    sGrab "${i}"
    echo "Got ${i}"
    echo ""
done

unset files_gz
files_gz=(
    129S1-SvImJ-to-mm10.over.chain.gz
    mm10-to-129S1-SvImJ.over.chain.gz
    CAST-EiJ-to-mm10.over.chain.gz
    mm10-to-CAST-EiJ.over.chain.gz
)

unset files
files=(
    129S1-SvImJ-to-mm10.over.chain
    mm10-to-129S1-SvImJ.over.chain
    CAST-EiJ-to-mm10.over.chain
    mm10-to-CAST-EiJ.over.chain
)

parallel "gunzip < {1} > {2}" \
::: "${files_gz[@]}" \
:::+ "${files[@]}"

# -----------------------------------------------------------------------------
#  - https://stackoverflow.com/questions/5410757/how-to-delete-from-a-text-file-all-lines-that-contain-a-specific-string
#+ - https://www.unix.com/shell-programming-and-scripting/46708-delete-block-text-delimited-blank-lines-when-pattern-found.html
#+ - https://unix.stackexchange.com/questions/552188/how-to-remove-empty-lines-from-beginning-and-end-of-file
#+ - https://unix.stackexchange.com/questions/31947/how-to-add-a-newline-to-the-end-of-a-file


# -----------------
#  https://www.ebi.ac.uk/ena/browser/view/GCA_001624185?show=chromosomes
chain_129=(
    129S1-SvImJ-to-mm10.over.chain
    mm10-to-129S1-SvImJ.over.chain
)
for chain in "${chain_129[@]}"; do
    #  Cut LVX, KV, and chrM entries from the liftOver chain file
    cat "${chain}" \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /LVX/!p' \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /KV/!p' \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /CM004179/!p' \
    | sed '/./,$!d' \
    > "${chain}.tmp"

    #  Add an additional newline to the end of the file, making for a total of two newlines
    echo >> "${chain}.tmp"

    #  Convert the chromosome names from official to common
    cat "${chain}.tmp" \
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
    > "${chain}.munged"

    rm "${chain}.tmp"
done


# -----------------
#  https://www.ebi.ac.uk/ena/browser/view/GCA_001624445?show=chromosomes
chain_CAST=(
    CAST-EiJ-to-mm10.over.chain
    mm10-to-CAST-EiJ.over.chain
)
for chain in "${chain_CAST[@]}"; do
    #  Cut LVX, KV, and chrM entries from the liftOver chain file
    cat "${chain}" \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /LVX/!p' \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /KV/!p' \
    | sed -n '
    /^$/ b block
    H
    $ b block
    b
    :block
    x
    /CM004181/!p' \
    | sed '/./,$!d' \
    > "${chain}.tmp"

    #  Add an additional newline to the end of the file, making for a total of two newlines
    echo >> "${chain}.tmp"

    #  Convert the chromosome names from official to common
    cat "${chain}.tmp" \
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
    > "${chain}.munged"

    rm "${chain}.tmp"
done


#  Clean up
files=(
    129S1-SvImJ-to-mm10.over.chain
    mm10-to-129S1-SvImJ.over.chain
    CAST-EiJ-to-mm10.over.chain
    mm10-to-CAST-EiJ.over.chain
)

parallel "rm {1}" \
::: "${files[@]}"

parallel "gzip {1}" \
::: "$(ls -- *.munged)"  #TODO Find an alternative to parsing ls

# mkdir -p "gz"
# parallel "mv {1} gz/" \
# ::: "${files_gz[@]}"

# mv ./gz/* ..

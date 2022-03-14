#!/bin/bash

path="/net/noble/vol1/home/kga0/genomes"
path_GB="/net/noble/vol2/home/gbonora/refdata/mm10.Nmasked.Cast129"
# path_CX="/net/shendure/vol10/projects/cxqiu/nobackup/genome"

#  "${path}/TBD.mm10_129S1-SvImJ-SNPs-inserted"
#+ is from
#+ "${path_GB}/129S1_SvImJ_full_sequence"
#+ 
#+ "${path}/TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-N-masked"
#+ is from
#+ "${path_GB}/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_N-masked"
#+
#+ "${path}/TBD.mm10_CAST-EiJ-SNPs-inserted"
#+ is from
#+ "${path_GB}/CAST_EiJ_full_sequence"
#+
#+ "${path}/TBD.mm10_CAST-EiJ-SNPs-N-masked"
#+ is from
#+ "${path_CX}/mm10.Cast_Nmasked"
#+
#+  "${path}/TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-inserted"
#+ is from
#+ "${path_GB}/CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_full_sequence"
#+ 
#+ "${path}/TBD.mm10_129S1-SvImJ-SNPs-N-masked"
#+ is from
#+ "${path_GB}/129S1_SvImJ_N-masked"


# -----------------
unset directories
directories=(
    "${path}/TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-inserted"
    "${path}/TBD.mm10_129S1-SvImJ-SNPs-N-masked"
)
for i in "${directories[@]}"; do echo "${i}"; done

parallel "mkdir -p {}/per-chromosome" ::: "${directories[@]}"

unset directories_from_GB
typeset -A directories_from_GB=(
    ["TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-inserted"]="CAST_EiJ_129S1_SvImJ_dual_hybrid.based_on_mm10.Nmasked.Cast129_full_sequence"
    ["TBD.mm10_129S1-SvImJ-SNPs-N-masked"]="129S1_SvImJ_N-masked"
)

for i in "${!directories_from_GB[@]}"; do
    echo "                          key (i): ${i}"
    echo "value (directories_from_GB[${i}]): ${directories_from_GB[${i}]}"
    echo "Start copying."
    cp -r "${path_GB}/${directories_from_GB[${i}]}" "${path}/${i}/per-chromosome"
    echo "Copying completed."
    echo ""
done

for i in "${!directories_from_GB[@]}"; do
    echo "Begin clean up for ${path}/${i}"
    cd "${path}/${i}/per-chromosome/${directories_from_GB[${i}]}" || 
        {
            echo "Failure to cd to ${path}/${i}/per-chromosome/${directories_from_GB[${i}]}. Exiting."
            exit 1
        }
    mv -- *.fa ..
    rmdir "${path}/${i}/per-chromosome/${directories_from_GB[${i}]}"
    echo "Finished clean up for ${path}/${i}"
    echo ""
done

for i in "${directories[@]}"; do
    echo ""
    echo "Directory is ${i}/per-chromosome"
    echo ""
    cd "${i}/per-chromosome" || 
        {
            echo "Failure to cd to ${i}/per-chromosome. Exiting."
            exit 1
        }

    for j in chr{1..19}.*.fa chrX.*.fa chrY.*.fa chrMT.*.fa; do
        echo "Working with ${j}"
        echo "Old header: $(head -n 1 "${j}")"

        sed -i "1s/^>/>chr/" "${j}"
        echo "New header: $(head -n 1 "${j}")"        
        echo ""
    done
done

for i in "${!directories_from_GB[@]}"; do
    echo ""
    echo "Begin cat'ing individual .fa in ${path}/${i}/per-chromosome"
    cd "${path}/${i}/per-chromosome" || 
    {
        echo "Failure to cd to ${path}/${i}/per-chromosome. Exiting."
        exit 1
    }

    echo "Running this command: cat chr{1..19}.*.fa chrX.*.fa chrY.*.fa chrMT.*.fa > \"../${i:4}.fa\""
    cat chr{1..19}.*.fa chrX.*.fa chrY.*.fa chrMT.*.fa > "../${i:4}.fa"
    echo "Command completed."
    echo ""
done

#NOTE Make sure pipeline-test_env is sourced...
for i in "${directories[@]}"; do
    echo ""
    echo "Directory is ${i}"
    echo ""
    cd "${i}" || 
        {
            echo "Failure to cd to ${i}. Exiting."
            exit 1
        }

    fasta="$(ls -- mm10*.fa)"
    echo "Working with ${fasta}"
    echo ""
    seqkit sort -nN "${fasta}" -o "${fasta/.fa/.sorted.fa}"
done

#NOTE Make sure pipeline-test_env is sourced...
for i in "${directories[@]}"; do
    echo ""
    echo "Directory is ${i}"
    echo ""
    cd "${i}" || 
        {
            echo "Failure to cd to ${i}. Exiting."
            exit 1
        }

    fasta="$(ls -- *.sorted.fa)"
    echo "Converting ${fasta} to 2bit format."
    faToTwoBit "${fasta}" "${fasta/.fa/.2bit}"
    echo "Conversion completed."
    echo ""
done

# -----------------
directories_N_masked=(
    "${path}/TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-N-masked"
    "${path}/TBD.mm10_CAST-EiJ-SNPs-N-masked"
)

parallel "mkdir -p {}" ::: "${directories_N_masked[@]}"

for i in "${directories_N_masked[@]}"; do
    echo ""
    echo "Directory is ${i}"
    echo ""
    cd "${i}" || 
        {
            echo "Failure to cd to ${i}. Exiting."
            exit 1
        }

    fasta="$(ls -- GRCm38*.fa)"
    echo "Working with ${fasta}"
    echo ""
    seqkit sort -nN "${fasta}" -o "${fasta/.fa/.sorted.fa}"
done
# Directory is /net/noble/vol1/home/kga0/genomes/TBD.mm10_CAST-EiJ-129S1-SvImJ-SNPs-N-masked
#
# Working with GRCm38-primary-assembly.Cast129_N-masked.fa
#
# [INFO] read sequences ...
# [INFO] 22 sequences loaded
# [INFO] sorting ...
# [INFO] output ...
#
# Directory is /net/noble/vol1/home/kga0/genomes/TBD.mm10_CAST-EiJ-SNPs-N-masked
#
# Working with GRCm38-primary-assembly.Cast_N-masked.fa
#
# [INFO] read sequences ...
# [INFO] 22 sequences loaded
# [INFO] sorting ...
# [INFO] output ...


# -----------------
directories_SNPs=(
    "${path}/TBD.mm10_129S1-SvImJ-SNPs-inserted"
    "${path}/TBD.mm10_CAST-EiJ-SNPS-inserted"
)

parallel "mkdir -p {}" ::: "${directories_SNPs[@]}"

suffix_SNPs="SNPs_introduced.fa"

for i in "${directories_SNPs[@]}"; do
    echo ""
    echo "Directory is ${i}"
    echo ""
    cd "${i}" || exit

    for j in chr{1..19}.${suffix_SNPs} chrX.${suffix_SNPs}; do
        echo "Working with ${j}"
        echo "Old header: $(head -n 1 ${j})"

        sed -i "1s/^>/>chr/" "${j}"
        echo "New header: $(head -n 1 ${j})"        
        echo ""
    done
done
# Directory is /net/noble/vol1/home/kga0/genomes/TBD.mm10_129S1-SvImJ-SNPs-inserted
#
# Working with chr1.SNPs_introduced.fa
# Old header: >1
# New header: >chr1
#
# Working with chr2.SNPs_introduced.fa
# Old header: >2
# New header: >chr2
#
# Working with chr3.SNPs_introduced.fa
# Old header: >3
# New header: >chr3
#
# Working with chr4.SNPs_introduced.fa
# Old header: >4
# New header: >chr4
#
# Working with chr5.SNPs_introduced.fa
# Old header: >5
# New header: >chr5
#
# Working with chr6.SNPs_introduced.fa
# Old header: >6
# New header: >chr6
#
# Working with chr7.SNPs_introduced.fa
# Old header: >7
# New header: >chr7
#
# Working with chr8.SNPs_introduced.fa
# Old header: >8
# New header: >chr8
#
# Working with chr9.SNPs_introduced.fa
# Old header: >9
# New header: >chr9
#
# Working with chr10.SNPs_introduced.fa
# Old header: >10
# New header: >chr10
#
# Working with chr11.SNPs_introduced.fa
# Old header: >11
# New header: >chr11
#
# Working with chr12.SNPs_introduced.fa
# Old header: >12
# New header: >chr12
#
# Working with chr13.SNPs_introduced.fa
# Old header: >13
# New header: >chr13
#
# Working with chr14.SNPs_introduced.fa
# Old header: >14
# New header: >chr14
#
# Working with chr15.SNPs_introduced.fa
# Old header: >15
# New header: >chr15
#
# Working with chr16.SNPs_introduced.fa
# Old header: >16
# New header: >chr16
#
# Working with chr17.SNPs_introduced.fa
# Old header: >17
# New header: >chr17
#
# Working with chr18.SNPs_introduced.fa
# Old header: >18
# New header: >chr18
#
# Working with chr19.SNPs_introduced.fa
# Old header: >19
# New header: >chr19
#
# Working with chrX.SNPs_introduced.fa
# Old header: >X
# New header: >chrX
#
#
# Directory is /net/noble/vol1/home/kga0/genomes/TBD.mm10_CAST-EiJ-SNPS-inserted
#
# Working with chr1.SNPs_introduced.fa
# Old header: >1
# New header: >chr1
#
# Working with chr2.SNPs_introduced.fa
# Old header: >2
# New header: >chr2
#
# Working with chr3.SNPs_introduced.fa
# Old header: >3
# New header: >chr3
#
# Working with chr4.SNPs_introduced.fa
# Old header: >4
# New header: >chr4
#
# Working with chr5.SNPs_introduced.fa
# Old header: >5
# New header: >chr5
#
# Working with chr6.SNPs_introduced.fa
# Old header: >6
# New header: >chr6
#
# Working with chr7.SNPs_introduced.fa
# Old header: >7
# New header: >chr7
#
# Working with chr8.SNPs_introduced.fa
# Old header: >8
# New header: >chr8
#
# Working with chr9.SNPs_introduced.fa
# Old header: >9
# New header: >chr9
#
# Working with chr10.SNPs_introduced.fa
# Old header: >10
# New header: >chr10
#
# Working with chr11.SNPs_introduced.fa
# Old header: >11
# New header: >chr11
#
# Working with chr12.SNPs_introduced.fa
# Old header: >12
# New header: >chr12
#
# Working with chr13.SNPs_introduced.fa
# Old header: >13
# New header: >chr13
#
# Working with chr14.SNPs_introduced.fa
# Old header: >14
# New header: >chr14
#
# Working with chr15.SNPs_introduced.fa
# Old header: >15
# New header: >chr15
#
# Working with chr16.SNPs_introduced.fa
# Old header: >16
# New header: >chr16
#
# Working with chr17.SNPs_introduced.fa
# Old header: >17
# New header: >chr17
#
# Working with chr18.SNPs_introduced.fa
# Old header: >18
# New header: >chr18
#
# Working with chr19.SNPs_introduced.fa
# Old header: >19
# New header: >chr19
#
# Working with chrX.SNPs_introduced.fa
# Old header: >X
# New header: >chrX

# -----------------
cd "${path}/TBD.mm10_129S1-SvImJ-SNPs-inserted" || exit

echo "Cat'ing files..."
cat chr{1..19}.${suffix_SNPs} chrX.${suffix_SNPs} > "mm10_129S1-SvImJ-SNPs-inserted.fa"

fasta="mm10_129S1-SvImJ-SNPs-inserted.fa"
faToTwoBit "${fasta}" "${fasta/.fa/.2bit}"

twoBitInfo "${fasta/.fa/.2bit}" stdout
# chr1    195471971
# chr2    182113224
# chr3    160039680
# chr4    156508116
# chr5    151834684
# chr6    149736546
# chr7    145441459
# chr8    129401213
# chr9    124595110
# chr10   130694993
# chr11   122082543
# chr12   120129022
# chr13   120421639
# chr14   124902244
# chr15   104043685
# chr16   98207768
# chr17   94987271
# chr18   90702639
# chr19   61431566
# chrX    171031299

# -----------------
cd "${path}/TBD.mm10_CAST-EiJ-SNPS-inserted" || exit

echo "Cat'ing files..."
cat chr{1..19}.${suffix_SNPs} chrX.${suffix_SNPs} > "mm10_CAST-EiJ-SNPS-inserted.fa"

fasta="mm10_CAST-EiJ-SNPS-inserted.fa"
faToTwoBit "${fasta}" "${fasta/.fa/.2bit}"

twoBitInfo "${fasta/.fa/.2bit}" stdout
# chr1    195471971
# chr2    182113224
# chr3    160039680
# chr4    156508116
# chr5    151834684
# chr6    149736546
# chr7    145441459
# chr8    129401213
# chr9    124595110
# chr10   130694993
# chr11   122082543
# chr12   120129022
# chr13   120421639
# chr14   124902244
# chr15   104043685
# chr16   98207768
# chr17   94987271
# chr18   90702639
# chr19   61431566
# chrX    171031299

# -----------------
cd "${path}" || exit

directories=(
    "${path}/Ensembl.129S1-SvImJ"
    "${path}/Ensembl.C57BL-6NJ"
    "${path}/Ensembl.CAST-EiJ"
    "${path}/Ensembl.SPRET-EiJ"
)

for i in "${directories[@]}"; do
    echo "Directory is ${i}"
    cd "${i}" || exit
    fasta="$(ls -- *.fa)"
    echo "Fasta is ${fasta}"

    faToTwoBit "${fasta}" "${fasta/.fa/.2bit}"
    echo "Fasta is ${fasta/.fa/.2bit}"
    echo ""
done


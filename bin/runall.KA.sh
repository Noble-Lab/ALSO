#!/bin/bash

#  Remote ---------------------------------------------------------------------

#TODO Code for work with faSize
#TODO Code for renaming chromosomes and munging, etc. (chromosome_sizes.txt)
#TODO Code for setting up whitelists
#TODO Checks, e.g., don't run if file present, etc.

genomes_directory="${HOME}/genomes"


# -------------------------------------
#  Get motif database file (such as those from JASPAR) to use for motif searching in peaks
file_JASPAR="JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.pfm"
[[ -f "data/motifs/${file_JASPAR}" ]] ||
    {
        curl "http://jaspar.genereg.net/download/CORE/${file_JASPAR/.pfm/.txt}" \
        > "data/motifs/${file_JASPAR}"
    }


# -------------------------------------
#  Ensembl.129S1-SvImJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624185.1
name="Ensembl.129S1-SvImJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXH01.fasta.gz" \
> "${genomes_directory}/${name}/LVXH01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624185.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624185.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.abinitio.gtf.gz
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.chr.gtf.gz
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_129s1svimj/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done

array_files=(
    CHECKSUMS
    Mus_musculus_129s1svimj.129S1_SvImJ_v1.104.gff3.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gff3"
    mkdir -p "${genomes_directory}/${name}/gff3"
    curl "http://ftp.ensembl.org/pub/release-104/gff3/mus_musculus_129s1svimj/${i}" \
    > "${genomes_directory}/${name}/gff3/${i}"
done


# -------------------------------------
#  Ensembl.C57BL-6NJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001632555.1
name="Ensembl.C57BL-6NJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXM01.fasta.gz" \
> "${genomes_directory}/${name}/LVXM01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001632555.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001632555.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.abinitio.gtf.gz
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.chr.gtf.gz
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_c57bl6nj/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done

array_files=(
    CHECKSUMS
    Mus_musculus_c57bl6nj.C57BL_6NJ_v1.104.gff3.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gff3"
    mkdir -p "${genomes_directory}/${name}/gff3"
    curl "http://ftp.ensembl.org/pub/release-104/gff3/mus_musculus_c57bl6nj/${i}" \
    > "${genomes_directory}/${name}/gff3/${i}"
done


# -------------------------------------
#  Ensembl.CAST-EiJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624445.1
name="Ensembl.CAST-EiJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXN01.fasta.gz" \
> "${genomes_directory}/${name}/LVXN01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624445.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624445.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus_casteij.CAST_EiJ_v1.104.abinitio.gtf.gz
    Mus_musculus_casteij.CAST_EiJ_v1.104.chr.gtf.gz
    Mus_musculus_casteij.CAST_EiJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus_casteij/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done

array_files=(
    CHECKSUMS
    Mus_musculus_casteij.CAST_EiJ_v1.104.gff3.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gff3"
    mkdir -p "${genomes_directory}/${name}/gff3"
    curl "http://ftp.ensembl.org/pub/release-104/gff3/mus_musculus_casteij/${i}" \
    > "${genomes_directory}/${name}/gff3/${i}"
done

# -------------------------------------
#  Ensembl.SPRET-EiJ, https://www.ebi.ac.uk/ena/browser/view/GCA_001624865.1
name="Ensembl.SPRET-EiJ"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/lvx/LVXV01.fasta.gz" \
> "${genomes_directory}/${name}/LVXV01.fasta.gz"
curl "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_001624865.1?download=true&gzip=true" \
> "${genomes_directory}/${name}/GCA_001624865.1.fa.gz"

array_files=(
    CHECKSUMS
    Mus_spretus.SPRET_EiJ_v1.104.abinitio.gtf.gz
    Mus_spretus.SPRET_EiJ_v1.104.chr.gtf.gz
    Mus_spretus.SPRET_EiJ_v1.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_spretus/${i}" \
    > "${genomes_directory}/${name}/gtf/${i}"
done

array_files=(
    CHECKSUMS
    Mus_spretus.SPRET_EiJ_v1.104.gff3.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gff3"
    mkdir -p "${genomes_directory}/${name}/gff3"
    curl "http://ftp.ensembl.org/pub/release-104/gff3/mus_spretus/${i}" \
    > "${genomes_directory}/${name}/gff3/${i}"
done

# -------------------------------------
#  Ensembl.GRCm39
name="Ensembl.GRCm39"
mkdir -p "${genomes_directory}/${name}"
curl "http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
> "${genomes_directory}/${name}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus.GRCm39.104.abinitio.gtf.gz
    Mus_musculus.GRCm39.104.chr.gtf.gz
    Mus_musculus.GRCm39.104.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/${i}" > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
#  Ensembl.GRCm38
name="Ensembl.GRCm38"
mkdir -p "${genomes_directory}/${name}"
curl "ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz" \
> "${genomes_directory}/${name}/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"

array_files=(
    CHECKSUMS
    Mus_musculus.GRCm38.102.abinitio.gtf.gz
    Mus_musculus.GRCm38.102.chr.gtf.gz
    Mus_musculus.GRCm38.102.gtf.gz
    README
)
for i in "${array_files[@]}"; do
    # echo "mkdir -p ${genomes_directory}/${name}/gtf"
    mkdir -p "${genomes_directory}/${name}/gtf"
    curl "http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/${i}" > "${genomes_directory}/${name}/gtf/${i}"
done


# -------------------------------------
array_samples=(
    129S1-SvImJ
    C57BL-6NJ
    CAST-EiJ
    SPRET-EiJ
)

#  Munge fasta: Remove scaffolds, etc.
for i in "${array_samples[@]}"; do
    fasta_in="$(find "${genomes_directory}/Ensembl.${i}" -name "*.fa" -printf "%P\n")"
    fasta_out="original.${i}.fa"

    [[ ! -e "${fasta_out}" ]] || rm "${fasta_out}"

    cat "${genomes_directory}/Ensembl.${i}/${fasta_in}" | \
    seqkit grep -vrp "LVX|KV|CM004188|CM004181|CM004277|CM004179" | \
    seqkit rmdup > "${genomes_directory}/Ensembl.${i}/${fasta_out}"
done


#  Munge fasta: Rename headers ------------------------------------------------
typeset tab="$(printf '\t')"
typeset kv
typeset fa_in
typeset fa_out

rename_headers() {
    awk 'NR==FNR{ a[$1]=$2; next }NF==2{ $2=a[$2]; print ">" $2; next }1' \
    FS='\t' "${kv}" FS='>' "${fa_in}" > "${fa_out}"
}


# -----------------
name="Ensembl.129S1-SvImJ"
cd "${genomes_directory}/${name}" || exit 1

kv="129S1-SvImJ.kv.txt"
[[ ! -e "${kv}" ]] || rm "${kv}"

cat > "${kv}" << EOF
ENA|CM003934|CM003934.1 Mus musculus strain 129S1/SvImJ chromosome 1, whole genome shotgun sequence.${tab}chr1 1
ENA|CM003935|CM003935.1 Mus musculus strain 129S1/SvImJ chromosome 2, whole genome shotgun sequence.${tab}chr2 2
ENA|CM003936|CM003936.1 Mus musculus strain 129S1/SvImJ chromosome 3, whole genome shotgun sequence.${tab}chr3 3
ENA|CM003937|CM003937.1 Mus musculus strain 129S1/SvImJ chromosome 4, whole genome shotgun sequence.${tab}chr4 4
ENA|CM003938|CM003938.1 Mus musculus strain 129S1/SvImJ chromosome 5, whole genome shotgun sequence.${tab}chr5 5
ENA|CM003939|CM003939.1 Mus musculus strain 129S1/SvImJ chromosome 6, whole genome shotgun sequence.${tab}chr6 6
ENA|CM003940|CM003940.1 Mus musculus strain 129S1/SvImJ chromosome 7, whole genome shotgun sequence.${tab}chr7 7
ENA|CM003941|CM003941.1 Mus musculus strain 129S1/SvImJ chromosome 8, whole genome shotgun sequence.${tab}chr8 8
ENA|CM003942|CM003942.1 Mus musculus strain 129S1/SvImJ chromosome 9, whole genome shotgun sequence.${tab}chr9 9
ENA|CM003943|CM003943.1 Mus musculus strain 129S1/SvImJ chromosome 10, whole genome shotgun sequence.${tab}chr10 10
ENA|CM003944|CM003944.1 Mus musculus strain 129S1/SvImJ chromosome 11, whole genome shotgun sequence.${tab}chr11 11
ENA|CM003945|CM003945.1 Mus musculus strain 129S1/SvImJ chromosome 12, whole genome shotgun sequence.${tab}chr12 12
ENA|CM003946|CM003946.1 Mus musculus strain 129S1/SvImJ chromosome 13, whole genome shotgun sequence.${tab}chr13 13
ENA|CM003947|CM003947.1 Mus musculus strain 129S1/SvImJ chromosome 14, whole genome shotgun sequence.${tab}chr14 14
ENA|CM003948|CM003948.1 Mus musculus strain 129S1/SvImJ chromosome 15, whole genome shotgun sequence.${tab}chr15 15
ENA|CM003949|CM003949.1 Mus musculus strain 129S1/SvImJ chromosome 16, whole genome shotgun sequence.${tab}chr16 16
ENA|CM003950|CM003950.1 Mus musculus strain 129S1/SvImJ chromosome 17, whole genome shotgun sequence.${tab}chr17 17
ENA|CM003951|CM003951.1 Mus musculus strain 129S1/SvImJ chromosome 18, whole genome shotgun sequence.${tab}chr18 18
ENA|CM003952|CM003952.1 Mus musculus strain 129S1/SvImJ chromosome 19, whole genome shotgun sequence.${tab}chr19 19
ENA|CM003953|CM003953.1 Mus musculus strain 129S1/SvImJ chromosome X, whole genome shotgun sequence.${tab}chrX X
EOF

fa_in="original.129S1-SvImJ.fa"
fa_out="129S1-SvImJ.fa"

[[ ! -e "${fa_out}" ]] || rm "${fa_out}"
rename_headers


# -----------------
name="Ensembl.C57BL-6NJ"
cd "${genomes_directory}/${name}" || exit 1

kv="C57BL-6NJ.kv.txt"
[[ ! -e "${kv}" ]] || rm "${kv}"

cat > "${kv}" << EOF
ENA|CM004215|CM004215.1 Mus musculus strain C57BL/6NJ chromosome 1, whole genome shotgun sequence.${tab}chr1 1
ENA|CM004216|CM004216.1 Mus musculus strain C57BL/6NJ chromosome 2, whole genome shotgun sequence.${tab}chr2 2
ENA|CM004217|CM004217.1 Mus musculus strain C57BL/6NJ chromosome 3, whole genome shotgun sequence.${tab}chr3 3
ENA|CM004218|CM004218.1 Mus musculus strain C57BL/6NJ chromosome 4, whole genome shotgun sequence.${tab}chr4 4
ENA|CM004219|CM004219.1 Mus musculus strain C57BL/6NJ chromosome 5, whole genome shotgun sequence.${tab}chr5 5
ENA|CM004220|CM004220.1 Mus musculus strain C57BL/6NJ chromosome 6, whole genome shotgun sequence.${tab}chr6 6
ENA|CM004221|CM004221.1 Mus musculus strain C57BL/6NJ chromosome 7, whole genome shotgun sequence.${tab}chr7 7
ENA|CM004222|CM004222.1 Mus musculus strain C57BL/6NJ chromosome 8, whole genome shotgun sequence.${tab}chr8 8
ENA|CM004223|CM004223.1 Mus musculus strain C57BL/6NJ chromosome 9, whole genome shotgun sequence.${tab}chr9 9
ENA|CM004224|CM004224.1 Mus musculus strain C57BL/6NJ chromosome 10, whole genome shotgun sequence.${tab}chr10 10
ENA|CM004225|CM004225.1 Mus musculus strain C57BL/6NJ chromosome 11, whole genome shotgun sequence.${tab}chr11 11
ENA|CM004226|CM004226.1 Mus musculus strain C57BL/6NJ chromosome 12, whole genome shotgun sequence.${tab}chr12 12
ENA|CM004227|CM004227.1 Mus musculus strain C57BL/6NJ chromosome 13, whole genome shotgun sequence.${tab}chr13 13
ENA|CM004228|CM004228.1 Mus musculus strain C57BL/6NJ chromosome 14, whole genome shotgun sequence.${tab}chr14 14
ENA|CM004229|CM004229.1 Mus musculus strain C57BL/6NJ chromosome 15, whole genome shotgun sequence.${tab}chr15 15
ENA|CM004230|CM004230.1 Mus musculus strain C57BL/6NJ chromosome 16, whole genome shotgun sequence.${tab}chr16 16
ENA|CM004231|CM004231.1 Mus musculus strain C57BL/6NJ chromosome 17, whole genome shotgun sequence.${tab}chr17 17
ENA|CM004232|CM004232.1 Mus musculus strain C57BL/6NJ chromosome 18, whole genome shotgun sequence.${tab}chr18 18
ENA|CM004233|CM004233.1 Mus musculus strain C57BL/6NJ chromosome 19, whole genome shotgun sequence.${tab}chr19 19
ENA|CM004234|CM004234.1 Mus musculus strain C57BL/6NJ chromosome X, whole genome shotgun sequence.${tab}chrX X
EOF

fa_in="original.C57BL-6NJ.fa"
fa_out="C57BL-6NJ.fa"

[[ ! -e "${fa_out}" ]] || rm "${fa_out}"
rename_headers


# -----------------
name="Ensembl.CAST-EiJ"
cd "${genomes_directory}/${name}" || exit 1

kv="CAST-EiJ.kv.txt"
[[ ! -e "${kv}" ]] || rm "${kv}"

cat > "${kv}" << EOF
ENA|CM003994|CM003994.1 Mus musculus castaneus strain CAST/EiJ chromosome 1, whole genome shotgun sequence.${tab}chr1 1
ENA|CM003995|CM003995.1 Mus musculus castaneus strain CAST/EiJ chromosome 2, whole genome shotgun sequence.${tab}chr2 2
ENA|CM003996|CM003996.1 Mus musculus castaneus strain CAST/EiJ chromosome 3, whole genome shotgun sequence.${tab}chr3 3
ENA|CM003997|CM003997.1 Mus musculus castaneus strain CAST/EiJ chromosome 4, whole genome shotgun sequence.${tab}chr4 4
ENA|CM003998|CM003998.1 Mus musculus castaneus strain CAST/EiJ chromosome 5, whole genome shotgun sequence.${tab}chr5 5
ENA|CM003999|CM003999.1 Mus musculus castaneus strain CAST/EiJ chromosome 6, whole genome shotgun sequence.${tab}chr6 6
ENA|CM004000|CM004000.1 Mus musculus castaneus strain CAST/EiJ chromosome 7, whole genome shotgun sequence.${tab}chr7 7
ENA|CM004001|CM004001.1 Mus musculus castaneus strain CAST/EiJ chromosome 8, whole genome shotgun sequence.${tab}chr8 8
ENA|CM004002|CM004002.1 Mus musculus castaneus strain CAST/EiJ chromosome 9, whole genome shotgun sequence.${tab}chr9 9
ENA|CM004003|CM004003.1 Mus musculus castaneus strain CAST/EiJ chromosome 10, whole genome shotgun sequence.${tab}chr10 10
ENA|CM004004|CM004004.1 Mus musculus castaneus strain CAST/EiJ chromosome 11, whole genome shotgun sequence.${tab}chr11 11
ENA|CM004005|CM004005.1 Mus musculus castaneus strain CAST/EiJ chromosome 12, whole genome shotgun sequence.${tab}chr12 12
ENA|CM004006|CM004006.1 Mus musculus castaneus strain CAST/EiJ chromosome 13, whole genome shotgun sequence.${tab}chr13 13
ENA|CM004007|CM004007.1 Mus musculus castaneus strain CAST/EiJ chromosome 14, whole genome shotgun sequence.${tab}chr14 14
ENA|CM004008|CM004008.1 Mus musculus castaneus strain CAST/EiJ chromosome 15, whole genome shotgun sequence.${tab}chr15 15
ENA|CM004009|CM004009.1 Mus musculus castaneus strain CAST/EiJ chromosome 16, whole genome shotgun sequence.${tab}chr16 16
ENA|CM004010|CM004010.1 Mus musculus castaneus strain CAST/EiJ chromosome 17, whole genome shotgun sequence.${tab}chr17 17
ENA|CM004011|CM004011.1 Mus musculus castaneus strain CAST/EiJ chromosome 18, whole genome shotgun sequence.${tab}chr18 18
ENA|CM004012|CM004012.1 Mus musculus castaneus strain CAST/EiJ chromosome 19, whole genome shotgun sequence.${tab}chr19 19
ENA|CM004013|CM004013.1 Mus musculus castaneus strain CAST/EiJ chromosome X, whole genome shotgun sequence.${tab}chrX X
EOF

fa_in="original.CAST-EiJ.fa"
fa_out="CAST-EiJ.fa"

[[ ! -e "${fa_out}" ]] || rm "${fa_out}"
rename_headers


# -----------------
name="Ensembl.SPRET-EiJ"
cd "${genomes_directory}/${name}" || exit 1

kv="SPRET-EiJ.kv.txt"
[[ ! -e "${kv}" ]] || rm "${kv}"

cat > "${kv}" << EOF
ENA|CM004094|CM004094.1 Mus spretus strain SPRET/EiJ chromosome 1, whole genome shotgun sequence.${tab}chr1 1
ENA|CM004095|CM004095.1 Mus spretus strain SPRET/EiJ chromosome 2, whole genome shotgun sequence.${tab}chr2 2
ENA|CM004096|CM004096.1 Mus spretus strain SPRET/EiJ chromosome 3, whole genome shotgun sequence.${tab}chr3 3
ENA|CM004097|CM004097.1 Mus spretus strain SPRET/EiJ chromosome 4, whole genome shotgun sequence.${tab}chr4 4
ENA|CM004098|CM004098.1 Mus spretus strain SPRET/EiJ chromosome 5, whole genome shotgun sequence.${tab}chr5 5
ENA|CM004099|CM004099.1 Mus spretus strain SPRET/EiJ chromosome 6, whole genome shotgun sequence.${tab}chr6 6
ENA|CM004100|CM004100.1 Mus spretus strain SPRET/EiJ chromosome 7, whole genome shotgun sequence.${tab}chr7 7
ENA|CM004101|CM004101.1 Mus spretus strain SPRET/EiJ chromosome 8, whole genome shotgun sequence.${tab}chr8 8
ENA|CM004102|CM004102.1 Mus spretus strain SPRET/EiJ chromosome 9, whole genome shotgun sequence.${tab}chr9 9
ENA|CM004103|CM004103.1 Mus spretus strain SPRET/EiJ chromosome 10, whole genome shotgun sequence.${tab}chr10 10
ENA|CM004105|CM004105.1 Mus spretus strain SPRET/EiJ chromosome 11, whole genome shotgun sequence.${tab}chr11 11
ENA|CM004106|CM004106.1 Mus spretus strain SPRET/EiJ chromosome 12, whole genome shotgun sequence.${tab}chr12 12
ENA|CM004107|CM004107.1 Mus spretus strain SPRET/EiJ chromosome 13, whole genome shotgun sequence.${tab}chr13 13
ENA|CM004108|CM004108.1 Mus spretus strain SPRET/EiJ chromosome 14, whole genome shotgun sequence.${tab}chr14 14
ENA|CM004109|CM004109.1 Mus spretus strain SPRET/EiJ chromosome 15, whole genome shotgun sequence.${tab}chr15 15
ENA|CM004110|CM004110.1 Mus spretus strain SPRET/EiJ chromosome 16, whole genome shotgun sequence.${tab}chr16 16
ENA|CM004111|CM004111.1 Mus spretus strain SPRET/EiJ chromosome 17, whole genome shotgun sequence.${tab}chr17 17
ENA|CM004112|CM004112.1 Mus spretus strain SPRET/EiJ chromosome 18, whole genome shotgun sequence.${tab}chr18 18
ENA|CM004113|CM004113.1 Mus spretus strain SPRET/EiJ chromosome 19, whole genome shotgun sequence.${tab}chr19 19
ENA|CM004104|CM004104.1 Mus spretus strain SPRET/EiJ chromosome X, whole genome shotgun sequence.${tab}chrX X
EOF

fa_in="original.SPRET-EiJ.fa"
fa_out="SPRET-EiJ.fa"

[[ ! -e "${fa_out}" ]] || rm "${fa_out}"
rename_headers


# -----------------------------------------------------------------------------
cd "${genomes_directory}" || exit 1

#  Get effective genome sizes for samples
for i in "${array_samples[@]}"; do
    fasta_in="${i}.fa"
    count="$(faCount -summary "${genomes_directory}/Ensembl.${i}/${fasta_in}")"
    echo -e "${fasta_in}\n" "${count}\n" >> "data/effective-genome-sizes.txt"
done

#  Index .fa files
for i in "${array_samples[@]}"; do
    fasta_GCA="$(find "${genomes_directory}/Ensembl.${i}" -name "GCA*.fa" -printf "%P\n")"
    fasta_sample="$(find "${genomes_directory}/Ensembl.${i}" -name "${i}.fa" -printf "%P\n")"
    fasta_sample_original="$(find "${genomes_directory}/Ensembl.${i}" -name "original.${i}.fa" -printf "%P\n")"

    samtools faidx "${genomes_directory}/Ensembl.${i}/${fasta_GCA}"
    samtools faidx "${genomes_directory}/Ensembl.${i}/${fasta_sample}"
    samtools faidx "${genomes_directory}/Ensembl.${i}/${fasta_sample_original}"
done


#  Clean up, organize directories
for i in "${array_samples[@]}"; do mkdir -p "TBD.fasta_orig/fasta_orig_Ensembl.${i}"; done

for i in "${array_samples[@]}"; do
    fasta_GCA="$(find "${genomes_directory}/Ensembl.${i}" -name "GCA*.fa" -printf "%P\n")"
    fasta_GCA_fai="$(find "${genomes_directory}/Ensembl.${i}" -name "GCA*.fai" -printf "%P\n")"
    fasta_sample_original="$(find "${genomes_directory}/Ensembl.${i}" -name "original.${i}.fa" -printf "%P\n")"
    fasta_sample_original_fai="$(find "${genomes_directory}/Ensembl.${i}" -name "original.${i}.fa.fai" -printf "%P\n")"

    source_dir="${genomes_directory}/Ensembl.${i}"
    destination_dir="${genomes_directory}/TBD.fasta_orig/fasta_orig_Ensembl.${i}"

    cd "${source_dir}" || exit 1
    mv "${fasta_GCA}" "${fasta_GCA_fai}" "${destination_dir}"
    mv "${fasta_sample_original}" "${fasta_sample_original_fai}" "${destination_dir}"
    cd "${genomes_directory}" || exit 1
done


#  Prepare TSS, gene-body definitions
for i in "${array_samples[@]}"; do
    Ensembl_directory="${genomes_directory}/Ensembl.${i}"
    gtf="$(find "${Ensembl_directory}/gtf" -name "*.chr.gtf.gz" -printf "%P\n")"
    
    temp_TSS="${Ensembl_directory}/${i}.TSS.temp.bed"
    temp_gene_body="${Ensembl_directory}/${i}.temp.gene_bodies.bed"
    
    whitelist="${Ensembl_directory}/${i}.whitelist_regions.bed"

    final_TSS="${Ensembl_directory}/${i}.TSS.bed"
    final_gene_body="${Ensembl_directory}/${i}.gene_bodies.bed"

    chromosome_sizes="${Ensembl_directory}/${i}.chromosome_sizes.txt"

    final_gene_body_2k="${Ensembl_directory}/${i}.gene_bodies.plus_2kb_upstream.bed"

    [[ ! -e "${temp_TSS}" ]] || rm "${temp_TSS}"
    [[ ! -e "${temp_gene_body}" ]] || rm "${temp_gene_body}"
    [[ ! -e "${final_TSS}" ]] || rm "${final_TSS}"
    [[ ! -e "${final_gene_body}" ]] || rm "${final_gene_body}"
    [[ ! -e "${final_gene_body_2k}" ]] || rm "${final_gene_body_2k}"
    
    Rscript bin/generate_TSS_file.R \
    --infile "${Ensembl_directory}/gtf/${gtf}" \
    --outfile "${temp_TSS}"

    Rscript bin/generate_gene_body_file.R \
    --infile "${Ensembl_directory}/gtf/${gtf}" \
    --outfile "${temp_gene_body}"

    cat "${temp_TSS}" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${whitelist}" \
    | uniq > "${final_TSS}"

    cat "${temp_gene_body}" \
    | grep -i -E "protein|lincRNA|antisense|TR_V|TR_D|TR_J|TR_C|IG_V|IG_D|IG_J|IG_C|IG_LV" \
    | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$7,1,$6}' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools intersect -a stdin -b "${whitelist}" \
    | uniq > "${final_gene_body}"

    bedtools slop -i "${final_gene_body}" -s -l 2000 -r 0 -g "${chromosome_sizes}" \
    | sort -k1,1V -k2,2n -k3,3n > "${final_gene_body_2k}"
done


# -----------------------------------------------------------------------------
#  Get liftOver tracks

# -----------------
path="/net/noble/vol1/home/kga0/genomes/Ensembl.CAST-EiJ"
# cd "${path}" ||
#     {
#         echo "Failed to cd into ${path}. Exiting..."
#         exit 1
#     }

mkdir -p "${path}/liftOver"

#  liftOver: CAST-EiJ to mm10
file_liftOver="GCA_001624445.1ToMm10.over.chain.gz"
[[ -f "${path}/liftOver/${file_liftOver}" ]] ||
    {
        curl "https://hgdownload.soe.ucsc.edu/hubs/GCA/001/624/445/GCA_001624445.1/liftOver/${file_liftOver}" \
        > "${path}/liftOver/${file_liftOver}"
    }

cp "${path}/liftOver/${file_liftOver}" "${path}/liftOver/CAST-EiJ-to-mm10.over.chain.gz"

#  liftOver: mm10 to CAST-EiJ
file_liftOver="mm10ToGCA_001624445.1.over.chain.gz"
[[ -f "${path}/liftOver/${file_liftOver}" ]] ||
    {
        curl "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/${file_liftOver}" \
        > "${path}/liftOver/${file_liftOver}"
    }

cp "${path}/liftOver/${file_liftOver}" "${path}/liftOver/mm10-to-CAST-EiJ.over.chain.gz"

# -----------------
path="/net/noble/vol1/home/kga0/genomes/Ensembl.129S1-SvImJ"
# cd "${path}" ||
#     {
#         echo "Failed to cd into ${path}. Exiting..."
#         exit 1
#     }

mkdir -p "${path}/liftOver"

#  liftOver: 129S1-SvImJ to mm10
file_liftOver="GCA_001624185.1ToMm10.over.chain.gz"
[[ -f "${path}/liftOver/${file_liftOver}" ]] ||
    {
        curl "https://hgdownload.soe.ucsc.edu/hubs/GCA/001/624/185/GCA_001624185.1/liftOver/${file_liftOver}" \
        > "${path}/liftOver/${file_liftOver}"
    }

cp "${path}/liftOver/${file_liftOver}" "${path}/liftOver/129S1-SvImJ-to-mm10.over.chain.gz"

#  liftOver: mm10 to 129S1-SvImJ
file_liftOver="mm10ToGCA_001624185.1.over.chain.gz"
[[ -f "${path}/liftOver/${file_liftOver}" ]] ||
    {
        curl "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/${file_liftOver}" \
        > "${path}/liftOver/${file_liftOver}"
    }

cp "${path}/liftOver/${file_liftOver}" "${path}/liftOver/mm10-to-129S1-SvImJ.over.chain.gz"
